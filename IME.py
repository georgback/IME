#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 14:52:14 2021

@author: back1622
"""


import csv
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import itertools
import re
from scipy import stats as st
import FisherExact
import random
import scipy.stats as stat
from subprocess import call
import glob
import os
import statsmodels.stats.multitest as sts


#extraction of from respective Sequences (i.E. introns)
#gff_file = GFF3 file containing the sequences of interest

#frq_file= VCFtool summary file with the --freq option (relevant positions previously extracted with bedtools and gff file)
#threshold = how many allels diverged from the major allel are needed to count position as SNP
def get_SNP_positions_for_intron(frq_file,
                                 gff_file,
                                 output,
                                 threshhold):
    temp_lis=[]
    header=[]

    #read in frq file
    with open(frq_file,mode="rt") as frq:
        reader=csv.reader(frq,delimiter="\t")
        header=next(reader)
        for row in reader:
            #include first 4 columns
            temp_row=row[:4]
            rest_row=row[4:]

            #calculate the SNP frequency by 1- highest SNP value (the most common Variant is taken as "base")
            #split the columns at : to extract allel frequencies
            freq_list=[*map(lambda each:float(each.split(":")[1]),rest_row)]
            temp_row=[*map(int,temp_row)]
            #calculation of SNP frequency
            temp_row.append(1-max(freq_list))
            temp_lis.append(temp_row)



        frq.close()

    #create data frame from the resulting list
    df=pd.DataFrame(temp_lis)
    header=header[:4]

    header.append("SNP_percentage")
    df.columns=header
    temp_lis=0


    #incluede only entries which have at least 250 assecions measured (2 chromosoms per plant)
    df=df[df["N_CHR"]>500]


    #if threshhold is chosen, number of accession which have SNP is calculated
    if threshhold>0:
        df=df[(df["N_CHR"]*df["SNP_percentage"])>=threshhold]

    #read gff file and write results into new file
    with open(gff_file) as introns,open(output,mode="wt") as output:
        reader=csv.reader(introns,delimiter="\t")
        writer=csv.writer(output,delimiter="\t")
        writer.writerow(["Gene","Strand","Length","Pos_freq"])

        #loop through gff
        for row in reader:
            #ignore all non regular chromosoms (plastid/mitochondria)
            try:
                int(row[0])
            except:
                break

            #each line in gff resembles a intron
            #this line extracts all SNPs that are on the right chromosom and lie between start and end point of intron
            SNPs=df[(df["CHROM"]==int(row[0]))&(df["POS"]>=int(row[3]))&(df["POS"]<=int(row[4]))]
            pos=[]
            #each intron has the name of the gene, strand, length and finally the SNPs; position:frequency
            pos.append((row[8].split(";")[0]).split("=")[1])
            pos.append(row[6])
            pos.append(int(int(row[4])-int(row[3])+1))
            for index,element in SNPs.iterrows():
                pos.append(str(int(element["POS"])-int(row[3])+1)+":"+str(element["SNP_percentage"]))
            writer.writerow(pos)

        introns.close()
        output.close()


#kmer functions

#helper function; appends the given kmer in the dictionary by amout of found kmers
def count_function(k_mer,line,dictio):

    dictio[k_mer].append(line.count(k_mer))

#create reverse complement
def reverse_complement(String):
    complement={"A":"T","T":"A","G":"C","C":"G","Y":"R","R":"Y","S":"W","W":"S","*":"*","N":"N","K":"M","M":"K"}
    temp=list(String)
    temp=list(map(lambda each:complement[each],temp))
    String=''.join(temp)
    return String[::-1]

#creates two df with occurences of each kmer and a df the positions of each kmer in each gene
#autmoatically removes the splicing site (3 bp)
#seq_only if only sequences without description are given
def occurences(frq_file,fas,seq_only=False,k=6):
#crate list of each intron with SNP positions
    SNPs=[]
    bases=["A","T","G","C"]

    kmer=[''.join(p) for p in itertools.product(bases,repeat=k)]
    #kmer=np.array(kmer)
    occurrence={el:[] for el in kmer}

    after_SNPs={el:[] for el in kmer}


    SNP={el:0 for el in kmer}
    #loop through frq files and save the positions
    with open(frq_file) as frq:
        reader=csv.reader(frq,delimiter="\t")
        next(reader)
        for row in reader:
            templ=[]
            list(map(lambda each: templ.append(int(each.split(":")[0])),row[3:]))
            templ=np.array(templ,int)

            #exclude the splicing sites, if not needed disable
            #all SNPs at positions smaller 3 get to 0 or negative
            templ=templ-3
            #exclude all smaller equal 0 and all that are in the last 3 positions (-6 because the first three allready got removed)
            templ=templ[(templ>0)&(templ<=(int(row[2])-6))]
            SNPs.append(templ)
        frq.close()

    with open(fas,mode="rt") as fasta:
        row_count=0
        sum_len=0
        sum_SNP=0
        for count,line in enumerate(fasta):
            if not seq_only:
                if count%2==0:
                    continue


            #exclude splicing site, disable by comment if necessary
            line=line[3:-3]
            if len(line)<6:
                print(line)
                row_count=row_count+1
                continue
            sum_len=sum_len+len(line)
            sum_SNP=sum_SNP+len(SNPs[row_count])

            list(map(lambda each:count_function(each,line,occurrence),kmer))

            line=np.array(list(line))

            line[(SNPs[row_count]-1)]="*"
            line=''.join(line)
            list(map(lambda each:count_function(each,line,after_SNPs),kmer))
            row_count=row_count+1
        fasta.close()

    SNP_rate=sum_SNP/sum_len




    reverse_kmer=list(map(reverse_complement,kmer))
    redundant={}
    kmer_copy=kmer.copy()


    #associate the reverse complement to each
    while(len(reverse_kmer)>0):
        a=kmer_copy[0]
        b=reverse_kmer[0]
        redundant[a]=b
        try:
            kmer_copy.remove(a)
        except:
            pass
        try:
            reverse_kmer.remove(b)
        except:
            pass
        try:
            kmer_copy.remove(b)
        except:
            pass
        try:
            reverse_kmer.remove(a)
        except:
            pass

    #combine redundant kmers
    #could be done in the loop above, but its easier to understand this way





    oc_temp={}
    SNP_temp={}

    for key in redundant:
        if key== redundant[key]:
            oc_temp[key]=occurrence[key]
            SNP_temp[key]=after_SNPs[key]
            continue

        oc_temp[key]=occurrence[key]
        SNP_temp[key]=after_SNPs[key]
        for i in range(0,len(occurrence[key])):
            oc_temp[key][i]=oc_temp[key][i]+occurrence[redundant[key]][i]
            SNP_temp[key][i]=SNP_temp[key][i]+after_SNPs[redundant[key]][i]



    occurrence={}
    after_SNPs={}
    for key in oc_temp:
        occurrence[key]=sum(oc_temp[key])
        after_SNPs[key]=sum(SNP_temp[key])

    df_ml=pd.DataFrame(oc_temp)


    oc_temp={}
    SNP_temp={}


    df=pd.DataFrame(occurrence.items())
    df.columns=["kmer","occurrence"]
    df["after_SNP"]=after_SNPs.values()
    df["conservation"]= df["after_SNP"]/df["occurrence"]

    return([df,df_ml,SNP_rate])




def find_kmer_positions(k_mer,fasta_file,cutoff=3,strand=None):
    #df["Gene"]=pd.Series(dtype=str)
    reverse_kmer=reverse_complement(k_mer)
    palindrome=reverse_kmer==k_mer
    temp_lis=[]


    if strand is not None:
        with open(fasta_file) as fas:

            for count,line in enumerate(fas):
                if count%2==0:
                    continue
                if cutoff>0:
                    line=line[cutoff:-cutoff]

                if strand[int((count-1)/2)]:
                    line=reverse_complement(line)
                positions=[m.start() for m in re.finditer(k_mer,line)]
                pos_rev=[m.start() for m in re.finditer(reverse_kmer,line)]
                positions.extend(pos_rev)
                positions.sort()
                temp_lis.append([len(line),positions])


    else:

        with open(fasta_file) as fas:



            if cutoff>0:
                if palindrome:
                    for count,line in enumerate(fas):
                        if count%2==0:
                            continue
                        line=line[cutoff:-cutoff]
                        positions=[m.start() for m in re.finditer(k_mer,line)]
                        temp_lis.append([len(line),positions])
                else:

                     for count,line in enumerate(fas):
                        if count%2==0:
                            continue
                        line=line[cutoff:-cutoff]
                        positions=[m.start() for m in re.finditer(k_mer,line)]
                        pos_rev=[m.start() for m in re.finditer(reverse_kmer,line)]
                        positions.extend(pos_rev)
                        positions.sort()
                        temp_lis.append([len(line),positions])
            else:
                if palindrome:
                    for count,line in enumerate(fas):
                        if count%2==0:
                            continue

                        positions=[m.start() for m in re.finditer(k_mer,line)]
                        temp_lis.append([len(line),positions])
                else:

                     for count,line in enumerate(fas):
                        if count%2==0:
                            continue
                        positions=[m.start() for m in re.finditer(k_mer,line)]
                        pos_rev=[m.start() for m in re.finditer(reverse_kmer,line)]
                        positions.extend(pos_rev)
                        positions.sort()
                        temp_lis.append([len(line),positions])
    df=pd.DataFrame(temp_lis,columns=["length","positions"])
    fas.close()
    return df



#frq file, not summed frq file needed (see SNP_frequency.py)
def kmer_conservation(kme_pos,frq_file,cutoff=3,threshhold=100):
    temp_lis=[]
    header=[]

    #read in frq file
    with open(frq_file,mode="rt") as frq:
        reader=csv.reader(frq,delimiter="\t")
        header=next(reader)
        for row in reader:
            #include first 4 columns
            temp_row=row[:4]
            rest_row=row[4:]

            #calculate the SNP frequency by 1- highest SNP value (the most common Variant is taken as "base")
            #split the columns at : to extract allel frequencies
            freq_list=[*map(lambda each:float(each.split(":")[1]),rest_row)]
            temp_row=[*map(int,temp_row)]
            #calculation of SNP frequency
            temp_row.append(1-max(freq_list))
            temp_lis.append(temp_row)

    df=pd.DataFrame(temp_lis)
    header=header[:4]

    header.append("SNP_percentage")
    df.columns=header
    temp_lis=0


    df=df[df["N_CHR"]>500]


    #if threshhold is chosen, number of accession which have SNP is calculated
    if threshhold>0:
        df=df[(df["N_CHR"]*df["SNP_percentage"])>=threshhold]

    return df


#RNA has to be a list wiht 2gff files (for fasta1 and 2)or None
def fish_entrop(kmers,fasta1,fasta2,strand=None):
    fishers=[]
    entropy=[]
    k_mer_order=[]
    bins=np.linspace(0,1,11)
    strand1=None
    strand2=None
    if strand is not None:
        df_temp1=pd.read_csv(strand[0],sep="\t",header=None)
        df_temp2=pd.read_csv(strand[1],sep="\t",header=None)
        strand1=[x=="-" for x in df_temp1.iloc[:,6]]
        strand2=[x=="-" for x in df_temp2.iloc[:,6]]

    for kmer in kmers:

        #kmer=df_first.kmer.loc[k]
        k_mer_order.append(kmer)

        pos_first=find_kmer_positions(kmer,fasta1,3,strand1)
        pos=find_kmer_positions(kmer,fasta2,3,strand2)

        pos_first=pos_first[pos_first.positions.map(len)>0]
        pos=pos[pos.positions.map(len)>0]

        posf=[item for sub in pos_first.apply(lambda row:list(map((lambda item: item/row.length),row.positions)),axis=1) for item in sub]
        poso=[item for sub in pos.apply(lambda row:list(map((lambda item: item/row.length),row.positions)),axis=1) for item in sub]

        posf_his=np.histogram(posf,bins=bins)
        poso_his=np.histogram(poso,bins=bins)
        fishers.append( FisherExact.fisher_exact([posf_his[0],poso_his[0]],simulate_pval=True,workspace=100000))

        entrop=np.array([st.entropy(np.histogram(np.random.uniform(0,1,sum(posf_his[0])),bins=bins)[0]) for i in range(0,10000)])
        entro_f=st.entropy(posf_his[0])
        entropy.append(np.sum(entrop<=entro_f)/len(entrop))
        #print(st.ks_2samp(posf,[item for sublist in entrop2 for item in sublist]))
        #print(st.kstest(posf,"uniform"))

    df=pd.DataFrame()
    df["kmer"]=k_mer_order
    df["Fisher"]=fishers
    df["Entropy_emp_p"]=entropy
    return(df)



#gene name in fasta file required
def find_genes_containing_kmer(kmer,fastafile,cutoff=3,pos_res=[0,1],strand=None):

    #allows most nucletoid letter code
    replace={"Y":"[CT]","R":"[AG]","S":"[GC]","W":"[AT]","N":"[ATGC]","K":"[GT]","M":"[AC]"}
    pattern=re.compile("|".join(replace.keys()))

    reverse_kmer=reverse_complement(kmer)
    reverse_kmer=pattern.sub(lambda m: replace[re.escape(m.group(0))], reverse_kmer)
    kmer=pattern.sub(lambda m: replace[re.escape(m.group(0))], kmer)

    genes=[]
    genename=""

    #checking for kmer, without reversing - Strand
    if strand is None:

        with open(fastafile,mode="rt") as fas:
            for line in fas:
                if line[0]==">":
                    genename=line[1:line.find(":")]

                else:
                    if cutoff!=0:
                        line=line[cutoff:-cutoff]
                    line=line[round(pos_res[0]*len(line)):round(pos_res[1]*len(line))]
                    if (re.search(kmer,line)!=None)|(re.search(reverse_kmer,line)!=None):
                        genes.append(genename)
    else:
         with open(fastafile,mode="rt") as fas:
            count=0
            for line in fas:
                if line[0]==">":
                    genename=line[1:line.find(":")]

                else:
                    if cutoff!=0:
                        line=line[cutoff:-cutoff]
                    if strand[count]:
                        line=reverse_complement(line)
                    line=line[round(pos_res[0]*len(line)):round(pos_res[1]*len(line))]
                    if (re.search(kmer,line)!=None):
                        genes.append(genename)
                    count=count+1



    return(genes)



def find_comparable_kmers(df_with_occurances,kmer,allowed_deviation=0.1,excluded=[]):
    occ=df_with_occurances[df_with_occurances.kmer==kmer].occurrence
    occ=int(occ)
    possible_kmer=df_with_occurances[(df_with_occurances.kmer!=kmer)&(df_with_occurances.occurrence<(occ+allowed_deviation*occ))&(df_with_occurances.occurrence>(occ-allowed_deviation*occ))].kmer
    return([*filter(lambda x: x not in excluded,possible_kmer)])



#fasta must ontain information about
#should be changed from fasta to gff, but I'm too lazy
def differential_methylated_fasta(differential_meth_file,fasta_file):
    df_methylation=pd.read_csv(differential_meth_file,sep="\t")
    differential=[]
    sub_chrom=[]
    for i in df_methylation.chr.unique():
        sub_chrom.append(df_methylation[df_methylation.chr==i])
    with open(fasta_file) as fasta:
        for line in fasta:
            if ">" in line:
                pos=line.split("::")[1]
                pos=pos.split(":")
                chrom=int(pos[0])
                start=int(pos[1].split("-")[0])
                end=int(pos[1].split("-")[1])
                a=len(sub_chrom[chrom-1][(sub_chrom[chrom-1]["start"]<=start)&(sub_chrom[chrom-1]["end"]>start)])
                b=len(sub_chrom[chrom-1][(sub_chrom[chrom-1]["start"]>start)&(sub_chrom[chrom-1]["start"]<end)])
                differential.append(a+b)

            else:
                continue
    return(differential)





#requires the fasta to contain the exact positions of
def extract_relevant_methylation_with_fasta(methylation_file,fasta_file,output_file,meth_header=None):

    df_methylation=pd.read_csv(methylation_file,sep="\t",header=meth_header)
    df_methylation.columns=[0,1,2,3,4,5,6]

    meth_positions=[]
    with open(fasta_file) as fasta:
        for line in fasta:
            if ">" in line:
                pos=line.split("::")[1]
                pos=pos.split(":")
                chrom=int(pos[0])
                start=int(pos[1].split("-")[0])
                end=int(pos[1].split("-")[1])

                #get all relevant positions
                positions=(list(df_methylation[(df_methylation[0]==chrom)&(df_methylation[1]>=start)&(df_methylation[1]<=end)][1]))
                #calculate the positions in relation to the
                meth_positions.append([int(x)-start+1 for x in positions])
            else:
                continue

        #meth_positions=[[int(i) for i in x] for x in meth_positions]
        with open(output_file,"wt") as out:
            writer=csv.writer(out,delimiter="\t")
            for x in meth_positions:
                writer.writerow(x)
            out.close()
    return(pd.DataFrame(meth_positions))


#calculate how many
def chi2_cont_for_kmer(kmer,methylation_array,gene_array,fasta):
    genes=find_genes_containing_kmer(kmer,fasta)

    contained=[x in genes for x in gene_array]

    table=[[sum(methylation_array[contained]),len(methylation_array[contained])-sum(methylation_array[contained])],[sum(methylation_array),len(methylation_array)-sum(methylation_array)]]
    return(table,stat.chi2_contingency(table))



#RNA functions

def correlation_of_genes(genes,RNA_file,df_expression=None):


    #adding the possibility to give
    if df_expression is None:
        gene_list=list()

        with open(RNA_file,mode="rt") as RNA:
            reader=csv.reader(RNA,delimiter="\t")
            for row in reader:
                if row[0] in genes:
                    gene_list.append(row)
            RNA.close()

        #removed header=None

        #df_expression=pd.DataFrame(gene_list,header=None)
        df_expression=pd.DataFrame(gene_list)
    else:
        df_expression=df_expression[df_expression.iloc[:,0].isin(genes)]
    df_names=df_expression.iloc[:,0]
    corr=(((df_expression.iloc[:,1:]).astype(float)).T).corr()
    corr.columns=df_names
    corr.index=df_names
    return corr


#strand None, +/- are being ignored; For strand gff file required. Alternative is a list of booleans
def cor_for_kmer_list(kmers,fasta_file,RNA_file,df_expression=None,pos_res=[0,1],strand=None):
    cor_list=[]
    if isinstance(strand,str):
        df_temp=pd.read_csv(strand,sep="\t",header=None)
        strand=[x=="-" for x in df_temp.iloc[:,6]]

    for i in kmers:
        temp_gen=find_genes_containing_kmer(i,fasta_file,cutoff=3,pos_res=pos_res,strand=strand)

        #gen names contain the isoform, the RNA-seq file doesnt, so the isoform information is removed
        temp_gen=[i.split(".")[0] for i in temp_gen]
        if df_expression is None:

            cor_list.append(correlation_of_genes(temp_gen,RNA_file))
        else:
            cor_list.append(correlation_of_genes(temp_gen,df_expression=df_expression))
    return cor_list



def extract_go_Terms(Gene_list,GO_term_file,output_file):
    with open(GO_term_file) as GO,open(output_file,mode="wt") as out:
        reader=csv.reader(GO,delimiter="\t")
        writer=csv.writer(out,delimiter="\t")
        for line in reader:
            if line[0] in Gene_list:
                writer.writerow(line)
        GO.close()
        out.close()

def run_GO_term_analysis(gene_list_1,gene_list_2,name1,name2,output_name,
                         GO_file="",
                         output_direct="",GO_term_script=""):
    gene_list_1=np.unique(np.array(gene_list_1))
    gene_list_2=np.unique(np.array(gene_list_2))
    file1=output_direct+"/"+name1
    with open(file1, 'w') as f:
        for item in gene_list_1:
            f.write("%s\n" % item)
    file2=output_direct+"/"+name2
    with open(file2, 'w') as f:
        for item in gene_list_2:
            f.write("%s\n" % item)
    f=open(output_direct+"/"+output_name+".txt", 'a+')
    f.write(name1+" "+name2+"\n")
    call(["perl",GO_term_script,GO_file,file1,file2],stdout=f)

    f.close()



def cohensD(d1,d2,exclude_diagonal=False):


    if exclude_diagonal:
        d1=d1.copy()
        np.fill_diagonal(d1.values,-10)

        d2=d2.copy()
        np.fill_diagonal(d2.values,-10)
        d1=[*filter(lambda x:x!=-10,np.array(d1).flatten())]
        d2=[*filter(lambda x:x!=-10,np.array(d2).flatten())]


    else:

        d1=np.array(d1).flatten()
        d2=np.array(d2).flatten()
    n1, n2 = len(d1), len(d2)
    s1 = np.var(d1, ddof=1)
    s2=np.var(d2, ddof=1)
    s = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    return((np.mean(d1)-np.mean(d2))/s)



def run_dist_and_Go_random(cor_list,identifiers,genes,
                           RNA_file="",
                           pdf="",out_dir="",repetitions=1,df_expression=None):

    all_kmer_genes=[]
    all_rand_genes=[]
    bins=np.linspace(-1,1,41)
    pp=PdfPages(pdf)
    for count,i in enumerate(cor_list):

        indep_genes=[*filter(lambda x: (x not in i.columns),genes)]
        indep_genes=random.sample(indep_genes,len(i.columns))
        indep_corr=correlation_of_genes(indep_genes,RNA_file)
        d=cohensD(i,indep_corr)

        plt.figure()
        kmer=identifiers[count]
        sns.distplot(np.array(i).flatten(),bins,color="blue")
        sns.distplot(np.array(indep_corr).flatten(),bins,color="orange")

        plt.title(kmer+"\nEffec-size= "+str(d))
        plt.legend(labels=[kmer,"random"])
        if pdf!="":
            print(kmer)
            plt.savefig(pp,format="pdf")
        run_GO_term_analysis(i.columns,indep_genes,kmer,kmer+"_random",kmer+"_enrichment",output_direct=out_dir)
        all_kmer_genes.extend(i.columns)
        all_rand_genes.extend(indep_genes)
    run_GO_term_analysis(all_kmer_genes,all_rand_genes,"combined_genes_kmer","all_genes_as_comparison","all_kmer_vs_all_random_enrichment",output_direct=out_dir)
    pp.close()

#strand None, +/- wil be ignored; Otherwise gff file required
def run_dist_and_Go_sim_kmer(cor_list,identifiers,genes,df_occurence,fasta,
                           RNA_file="",
                           pdf="",out_dir="",strand=None,
                           pos_res=[0,1]):
    bins=np.linspace(-1,1,41)
    pp=PdfPages(pdf)

    all_kmer_genes=[]
    all_rand_genes=[]
    all_pos_kmers=[find_comparable_kmers(df_occurence,kmer,excluded=identifiers) for kmer in identifiers]

    if strand is not None:
        df_temp=pd.read_csv(strand,sep="\t",header=None)
        strand=[x=="-" for x in df_temp.iloc[:,6]]

    for count,i in enumerate(cor_list):
        if len(all_pos_kmers[count])==0:
            continue
        indep_kmer=random.sample(all_pos_kmers[count],1)[0]


        #fixed fasta maybe change
        indep_gen=find_genes_containing_kmer(indep_kmer,fasta,cutoff=3,pos_res=pos_res,strand=strand)
        indep_gen=[i.split(".")[0] for i in indep_gen]
        indep_corr=correlation_of_genes(indep_gen,RNA_file)


        d=cohensD(i,indep_corr)

        plt.figure()
        kmer=identifiers[count]
        sns.distplot(np.array(i).flatten(),bins,color="blue")
        sns.distplot(np.array(indep_corr).flatten(),bins,color="orange")

        plt.title(kmer+"\nEffec-size= "+str(d))
        plt.legend(labels=[kmer,indep_kmer])
        if pdf!="":
            print(kmer)
            plt.savefig(pp,format="pdf")
        run_GO_term_analysis(i.columns,indep_gen,kmer,indep_kmer+"_as_comparison",kmer+"_vs_"+indep_kmer+"_enrichment",output_direct=out_dir)
        all_kmer_genes.extend(i.columns)
        all_rand_genes.extend(indep_gen)
    run_GO_term_analysis(all_kmer_genes,all_rand_genes,"combined_genes_kmer","all_genes_from_random_kmers_as_comparison","all_kmer_vs_all_random_kmer_enrichment",output_direct=out_dir)
    pp.close()











#strand None, +/- wil be ignored; Otherwise gff file required
def run_dist_and_GO_random_kmer_all(cor_list,identifiers,genes,df_occurence,df_expression,
                           pdf="",out_dir="",
                           pos_res=[0,1],strand=None,fasta_file="",GO_component="",GO_function="",GO_process=""):
    bins=np.linspace(-1,1,41)
    #pp=PdfPages(pdf)

    all_pos_kmers=[find_comparable_kmers(df_occurence,kmer,excluded=identifiers) for kmer in identifiers]

    if strand is not None:
        df_temp=pd.read_csv(strand,sep="\t",header=None)
        strand=[x=="-" for x in df_temp.iloc[:,6]]

    for count,i in enumerate(cor_list):
        if len(all_pos_kmers[count])==0:
            continue
        kmer=identifiers[count]
        cur_dir=out_dir+"/"+kmer
        try:
            os.mkdir(cur_dir)
        except:
            pass
        pos_kmer=all_pos_kmers[count]
        print("extract all pos kmers")

        D=[]
        all_rand_genes=[]
        for indep_kmer in pos_kmer:

            #fixed fasta maybe change
            indep_gen=find_genes_containing_kmer(indep_kmer,fasta_file,cutoff=3,pos_res=pos_res,strand=strand)
            indep_gen=[i.split(".")[0] for i in indep_gen]
            indep_corr=correlation_of_genes(genes=indep_gen,df_expression=df_expression)


            d=cohensD(i,indep_corr,exclude_diagonal=True)
            D.append(d)

            #plt.figure()
            kmer=identifiers[count]
            #sns.distplot(np.array(i).flatten(),bins,color="blue")
            #sns.distplot(np.array(indep_corr).flatten(),bins,color="orange")

            #plt.title(kmer+"\nEffec-size= "+str(d))
            #plt.legend(labels=[kmer,indep_kmer])
            #if pdf!="":
             #   print(kmer)
              #  plt.savefig(pp,format="pdf")
            run_GO_term_analysis(i.columns,indep_gen,kmer,indep_kmer+"_as_comparison",kmer+"_component",output_direct=cur_dir,GO_file=GO_component)
            run_GO_term_analysis(i.columns,indep_gen,kmer,indep_kmer+"_as_comparison",kmer+"_function",output_direct=cur_dir,GO_file=GO_function)
            run_GO_term_analysis(i.columns,indep_gen,kmer,indep_kmer+"_as_comparison",kmer+"_process",output_direct=cur_dir,GO_file=GO_process)
            all_rand_genes.extend(indep_gen)
        run_GO_term_analysis(i.columns,all_rand_genes,"combined_genes_kmer","all_genes_from_random_kmers_as_comparison",kmer+"_vs_all_random_kmer_component",output_direct=cur_dir,GO_file=GO_component)
        run_GO_term_analysis(i.columns,all_rand_genes,"combined_genes_kmer","all_genes_from_random_kmers_as_comparison",kmer+"_vs_all_random_kmer_function",output_direct=cur_dir,GO_file=GO_function)
        run_GO_term_analysis(i.columns,all_rand_genes,"combined_genes_kmer","all_genes_from_random_kmers_as_comparison",kmer+"_vs_all_random_kmer_process",output_direct=cur_dir,GO_file=GO_process)
        df_d=pd.DataFrame(D)
        df_d["kmer"]=pos_kmer
        df_d.to_csv(cur_dir+"/"+kmer+"_CohensD.csv",sep="\t",index=False)
    #pp.close()


def return_W(data,conf=0.95):
        n = len(data)
        std_err = st.sem(data)
        h = std_err * st.t.ppf((1 + conf) / 2, n - 1)
        return(h)




def cohensD_go_term_random(cor_list,identifiers,genes,
                           RNA_file="",
                           out_dir="",repetitions=1,df_expression=None,GO_component="",GO_function="",GO_process=""):
    for count,i in enumerate(cor_list):
        kmer=identifiers[count]
        all_rand_genes=[]
        cur_dir=out_dir+"/"+kmer
        D=[]

        try:
            os.mkdir(cur_dir)
        except:
            pass
        for rep in range(0,repetitions):

            indep_genes=[*filter(lambda x: (x not in i.columns),genes)]
            indep_genes=random.sample(indep_genes,len(i.columns))
            indep_corr=correlation_of_genes(indep_genes,RNA_file,df_expression)
            d=cohensD(i,indep_corr)


            run_GO_term_analysis(i.columns,indep_genes,kmer,kmer+"_random"+str(rep),kmer+"_enrichment",output_direct=cur_dir)
            run_GO_term_analysis(i.columns,indep_genes,kmer,kmer+"_random"+str(rep),kmer+"_component",output_direct=cur_dir,GO_file=GO_component)
            run_GO_term_analysis(i.columns,indep_genes,kmer,kmer+"_random"+str(rep),kmer+"_function",output_direct=cur_dir,GO_file=GO_function)
            run_GO_term_analysis(i.columns,indep_genes,kmer,kmer+"_random"+str(rep),kmer+"_process",output_direct=cur_dir,GO_file=GO_process)
            all_rand_genes.extend(indep_genes)
            D.append(d)
        run_GO_term_analysis(i.columns,all_rand_genes,"combined_genes_kmer","all_genes_as_comparison","all_kmer_vs_all_random_enrichment",output_direct=cur_dir)
        run_GO_term_analysis(i.columns,all_rand_genes,"combined_genes_kmer","all_genes_from_random_as_comparison",kmer+"_vs_all_random_kmer_component",output_direct=cur_dir,GO_file="/home/mpimp-golm.mpg.de/back1622/tools/GOTermEnrichment/ATH_GO_GOSLIM_Feb18.component_slim")
        run_GO_term_analysis(i.columns,all_rand_genes,"combined_genes_kmer","all_genes_from_random_as_comparison",kmer+"_vs_all_random_kmer_function",output_direct=cur_dir,GO_file="/home/mpimp-golm.mpg.de/back1622/tools/GOTermEnrichment/ATH_GO_GOSLIM_Feb18.function_slim")
        run_GO_term_analysis(i.columns,all_rand_genes,"combined_genes_kmer","all_genes_from_random_as_comparison",kmer+"_vs_all_random_kmer_process",output_direct=cur_dir,GO_file="/home/mpimp-golm.mpg.de/back1622/tools/GOTermEnrichment/ATH_GO_GOSLIM_Feb18.process_slim")
        df_d=pd.DataFrame(D)
        df_d.to_csv(cur_dir+"/"+kmer+"_CohensD.csv",sep="\t")



def compare_to_top_gene_from_list(cor_list,identifiers,genes,gene_file_compare,
                           RNA_file="/home/mpimp-golm.mpg.de/back1622/Desktop/ExchangeMPI/AG Bioinformatics 2017/Georg/supercluster_log_base_e_quantileNorm_AGI_uniqueIDs.txt",
                           pdf="",out_dir="/home/mpimp-golm.mpg.de/back1622/tools/GOTermEnrichment/gene_files",df_expression=None,exclude_same_genes=False):

    #all_kmer_genes=[]
    #all_rand_genes=[]
    compare=pd.read_csv(gene_file_compare,sep="\t",header=None)
    compare=compare[0].map(lambda x: x.split(".")[0])
    compare=compare[compare.isin(genes)]
    compare=compare.unique()
    bins=np.linspace(-1,1,41)
    pp=PdfPages(pdf)
    for count,i in enumerate(cor_list):
        if exclude_same_genes:
            compare=[*filter(lambda x:x not in cor_list[count].index,compare)]
        indep_genes=list(compare[:len(i)])

        indep_corr=correlation_of_genes(indep_genes,RNA_file,df_expression=df_expression)
        d=cohensD(i,indep_corr,)

        plt.figure()
        kmer=identifiers[count]
        sns.distplot(np.array(i).flatten(),bins,color="blue")
        sns.distplot(np.array(indep_corr).flatten(),bins,color="orange")

        plt.title(kmer+"\nEffec-size= "+str(d))
        plt.legend(labels=[kmer,"top_imeter"])
        if pdf!="":
            print(kmer)
            plt.savefig(pp,format="pdf")
        #run_GO_term_analysis(i.columns,indep_genes,kmer,kmer+"_random",kmer+"_enrichment",output_direct=out_dir)
        #all_kmer_genes.extend(i.columns)
        #all_rand_genes.extend(indep_genes)
    #run_GO_term_analysis(all_kmer_genes,all_rand_genes,"combined_genes_kmer","all_genes_as_comparison","all_kmer_vs_all_random_enrichment",output_direct=out_dir)
    pp.close()


#more restrictive, easier to call
def compare_to_imeter_simplyfied(kmers,df_corr,gene_file_compare,fasta,
                                 pdf="",out_dir="/home/mpimp-golm.mpg.de/back1622/tools/GOTermEnrichment/gene_files",
                                 exclude_same_genes=False):


    compare=pd.read_csv(gene_file_compare,sep="\t",header=None)
    compare=compare[0].map(lambda x: x.split(".")[0])
    compare=compare[compare.isin(df_corr.index)]
    compare=compare.unique()
    bins=np.linspace(-1,1,41)
    pp=PdfPages(pdf)
    for count,kmer in enumerate(kmers):
        i=kmer_cor_sub_matrix(kmer, fasta, df_corr)
        if exclude_same_genes:
            compare=[*filter(lambda x:x not in i.index,compare)]
        indep_genes=list(compare[:len(i)])

        indep_corr=df_corr.loc[indep_genes,indep_genes]
        d=cohensD(i,indep_corr,)

        plt.figure()
        sns.distplot(np.array(i).flatten(),bins,color="blue")
        sns.distplot(np.array(indep_corr).flatten(),bins,color="orange")

        plt.title(kmer+"\nEffect-size= "+str(round(d,3)))
        plt.legend(labels=[kmer,"top_imeter"])
        plt.xlabel("Pearson's r")
        if pdf!="":
            print(kmer)
            plt.savefig(pp,format="pdf")
    pp.close()














def reverse_fasta_by_strand(fasta,gff,output):
    with open(fasta) as fas,open(gff) as gf,open(output,"wt") as out:
        reader=csv.reader(gf,delimiter="\t")
        for row in reader:
            header=next(fas)
            line=next(fas)
            if row[6]=="-":
                line=reverse_complement(line.rstrip())+"\n"
            out.writelines([header,line])
        out.close()

        fas.close()
        gf.close()



#helper function which takes in kmer,fasta file and correlation_matrix
#and returns submatrix with all genes that contain hexamer and are in the correlation matrix

def kmer_cor_sub_matrix(kmer,fasta,df_corr):
    genes=find_genes_containing_kmer(kmer,fasta)
    genes=[i.split(".")[0] for i in genes]
    genes=[*filter(lambda x: x in df_corr.index,genes)]
    return df_corr.loc[genes,genes]




def kmer_vs_comp_kmer(kmers,df_occurence,genes,df_corr,
                           pdf="",out_dir="",
                           fasta_file=""):
    bins=np.linspace(-1,1,41)
    #pp=PdfPages(pdf)

    all_pos_kmers=[find_comparable_kmers(df_occurence,kmer,excluded=kmers) for kmer in kmers]

    for count,kmer in enumerate(kmers):
        if len(all_pos_kmers[count])==0:
            continue
        cur_dir=out_dir+"/"+kmer
        try:
            os.mkdir(cur_dir)
        except:
            pass
        pos_kmer=all_pos_kmers[count]
        print("extract all pos kmers for comparison to "+kmer)

        i=kmer_cor_sub_matrix(kmer,fasta_file,df_corr)
        D=[]
        all_rand_genes=[]
        for indep_kmer in pos_kmer:

            indep_corr=kmer_cor_sub_matrix(indep_kmer,fasta_file,df_corr)


            d=cohensD(i,indep_corr,exclude_diagonal=True)
            D.append(d)

            #plt.figure()
            #sns.distplot(np.array(i).flatten(),bins,color="blue")
            #sns.distplot(np.array(indep_corr).flatten(),bins,color="orange")

            #plt.title(kmer+"\nEffec-size= "+str(d))
            #plt.legend(labels=[kmer,indep_kmer])
            #if pdf!="":
             #   print(kmer)
              #  plt.savefig(pp,format="pdf")

              #double ## for removing when goterm progrgam is installed

            run_GO_term_analysis(i.columns,indep_corr.columns,kmer,indep_kmer+"_as_comparison",kmer+"_component",output_direct=cur_dir,GO_file="/home/mpimp-golm.mpg.de/back1622/tools/GOTermEnrichment/ATH_GO_GOSLIM_Feb18.component_slim")
            run_GO_term_analysis(i.columns,indep_corr.columns,kmer,indep_kmer+"_as_comparison",kmer+"_function",output_direct=cur_dir,GO_file="/home/mpimp-golm.mpg.de/back1622/tools/GOTermEnrichment/ATH_GO_GOSLIM_Feb18.function_slim")
            run_GO_term_analysis(i.columns,indep_corr.columns,kmer,indep_kmer+"_as_comparison",kmer+"_process",output_direct=cur_dir,GO_file="/home/mpimp-golm.mpg.de/back1622/tools/GOTermEnrichment/ATH_GO_GOSLIM_Feb18.process_slim")
            all_rand_genes.extend(indep_corr.columns)
        run_GO_term_analysis(i.columns,all_rand_genes,"combined_genes_kmer","all_genes_from_random_kmers_as_comparison",kmer+"_vs_all_random_kmer_component",output_direct=cur_dir,GO_file="/home/mpimp-golm.mpg.de/back1622/tools/GOTermEnrichment/ATH_GO_GOSLIM_Feb18.component_slim")
        run_GO_term_analysis(i.columns,all_rand_genes,"combined_genes_kmer","all_genes_from_random_kmers_as_comparison",kmer+"_vs_all_random_kmer_function",output_direct=cur_dir,GO_file="/home/mpimp-golm.mpg.de/back1622/tools/GOTermEnrichment/ATH_GO_GOSLIM_Feb18.function_slim")
        run_GO_term_analysis(i.columns,all_rand_genes,"combined_genes_kmer","all_genes_from_random_kmers_as_comparison",kmer+"_vs_all_random_kmer_process",output_direct=cur_dir,GO_file="/home/mpimp-golm.mpg.de/back1622/tools/GOTermEnrichment/ATH_GO_GOSLIM_Feb18.process_slim")
        df_d=pd.DataFrame(D)
        df_d["kmer"]=pos_kmer
        df_d.to_csv(cur_dir+"/"+kmer+"_CohensD.csv",sep="\t",index=False)






#################################### Run analysis

#variables
expression_file=""
sep="\t"
Gene_ID=0
header=None
Log_norm=False
run_cohens_D=False
kmer_comp_dir=""
kmer_length=6
cohens_threshold=0.05
fasta_file=""
fasta_file_other_introns=""
pdf_folder=""
gff_file=""
gff_file_others=""
basic_frq_first=""
basic_frq_others=""

#output file
output_frq_first=""
output_frq_others=""

#reading in expression set
df_expression=pd.read_csv(expression_file,sep=sep,header=None)
df_expression.index=df_expression[Gene_ID]
df_expression=df_expression.drop(Gene_ID,1)
if Log_norm:
    df_expression=np.log10(df_expression+1)

df_expression=df_expression[df_expression.T.mean()>0]

#frq files
get_SNP_positions_for_intron(basic_frq_first, gff_file, output_frq_first, threshhold=50)
get_SNP_positions_for_intron(basic_frq_others, gff_file_others, output_frq_others, threshhold=50)

#determining all potential k_mers
df_first=occurences(output_frq_first,fasta_file)[0]
df_others=occurences(output_frq_others,fasta_file_other_introns)[0]
df_first_rel=df_first["occurrence"]/np.sum(df_first.occurrence)
df_others_rel=df_others["occurrence"]/np.sum(df_others.occurrence)


kmers=df_first.kmer

df_fis=fish_entrop(kmers,fasta_file,fasta_file_other_introns,strand=None)
df_sig=df_fis[(sts.multipletests(df_fis.Fisher,alpha=0.05,method="fdr_bh")[0])&(sts.multipletests(df_fis.Entropy_emp_p,alpha=0.05,method="fdr_bh")[0])]

df_sig.index=df_sig.kmer
df_sig_ordered=df_sig.loc[df_first[df_first.kmer.isin(df_sig.kmer)].kmer,]
df_sig_ordered.index=df_first[df_first.kmer.isin(df_sig_ordered.kmer)].index

#relative conservation for all sig kmer
df_sig_ordered["rel_conserved"]=df_first[df_first.kmer.isin(df_sig.kmer)].conservation/df_others[df_others.kmer.isin(df_sig.kmer)].conservation

#more common than expected for sig kmer
df_sig_ordered["rel_occurrence"]=(df_first_rel/df_others_rel).iloc[df_sig_ordered.index]

#remove all kmers not conserved or overrepresented
df_sig_ordered=df_sig_ordered[(df_sig_ordered.rel_conserved>1)&(df_sig_ordered.rel_occurrence>1)]

#correlation matrix
df_corr=pd.DataFrame(np.corrcoef(df_expression))
df_corr.columns=df_expression.index
df_corr.index=df_expression.index


#run cohensD analysis of kmers vs comparable kmers
if run_cohens_D:

    #this is for runs with breaks. removes kmers allready run from the list
    prerun_hexamer=[h[-kmer_length:] for h in glob.glob(kmer_comp_dir+"/*",recursive=True)]
    df_sig_ordered=df_sig_ordered[~df_sig_ordered.kmer.isin(prerun_hexamer)]

    kmer_vs_comp_kmer(df_sig_ordered.kmer,df_first,df_expression.index,df_corr, out_dir=kmer_comp_dir)
