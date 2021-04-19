#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 15:24:39 2021

@author: back1622
"""

import pandas as pd
import numpy as np
import IME as ime
import statsmodels.stats.multitest as sts
import glob
import re
import csv
import scipy.stats as st
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
ime.get_SNP_positions_for_intron(basic_frq_first, gff_file, output_frq_first, threshhold=50)
ime.get_SNP_positions_for_intron(basic_frq_others, gff_file_others, output_frq_others, threshhold=50)

#determining all potential k_mers
df_first=ime.occurences(output_frq_first,fasta_file)[0]
df_others=ime.occurences(output_frq_others,fasta_file_other_introns)[0]
df_first_rel=df_first["occurrence"]/np.sum(df_first.occurrence)
df_others_rel=df_others["occurrence"]/np.sum(df_others.occurrence)


kmers=df_first.kmer

df_fis=ime.fish_entrop(kmers,fasta_file,fasta_file_other_introns,strand=None)
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

    ime.kmer_vs_comp_kmer(df_sig_ordered.kmer,df_first,df_expression.index,df_corr, out_dir=kmer_comp_dir)


#analyzing CohensD of kmers
#without folder
cohens_D=[]
kmer_coh=[]
for i in glob.glob(kmer_comp_dir+"/*/*CohensD*",recursive=True):
    kmer_coh.append(i.split("/")[-2]+" "+i.split("/")[-3])
    cohens_D.append(pd.read_csv(i,sep="\t",index_col=1))
kmer_coh=[x.split(" ")[0] for x in kmer_coh]
df_cohens=pd.DataFrame({"kmer":kmer_coh,"mean_cohens":[float(x.mean()) for x in cohens_D]})

df_cohens.sort_values(by=['mean_cohens'], inplace=True, ascending=False)

with open(pdf_folder+"/motifs_kmer.fas","w") as f:
    for k in df_cohens[df_cohens.mean_cohens>=cohens_threshold].kmer:
        f.write(">"+k+"\n")
        f.write(k+"\n")
    f.close()

#now the kmers with mean cohens d values over threshold are chosen
#plotting

ime.compare_to_imeter_simplyfied(df_cohens[df_cohens.mean_cohens>=cohens_threshold].kmer, df_corr, "sorted_imeter_first_intron.csv", fasta_file,pdf=pdf_folder+"/kmer_vs_imeter.pdf")


ime.compare_to_imeter_simplyfied(df_cohens[df_cohens.mean_cohens>=cohens_threshold].kmer, df_corr, "sorted_imeter_first_intron.csv", fasta_file,pdf=pdf_folder+"/kmer_vs_imeter.pdf")


#comparing mean gene expression level

#chosing median expression to have a better representation
df_expression_median=df_expression.median(1)

#remove very low expressed genes
df_expression_median=df_expression_median[df_expression_median>0.05]

#machine learning data used to find overlap of introns with functional 5'UTR, which are also included in the
df_ML=pd.read_csv("ML_data_set.csv",sep="\t",index_col=0)


#exmaple for kmer
kmers=["ATCGAA"]
#WWTCG
effect_sizes=[]
for kmer in kmers:
    gen_kmer=ime.find_genes_containing_kmer(kmer,fasta_file)
    gen_kmer=[gen.split(".")[0] for gen in gen_kmer]
    gen_kmer=np.unique(gen_kmer)
    gen_kmer=pd.Series(gen_kmer)
    gen_kmer=gen_kmer[gen_kmer.isin(df_expression_median.index)]
    gen_kmer=gen_kmer[gen_kmer.isin(df_ML.index)]
    set_non_kmer=df_expression_median[(~df_expression_median.index.isin(gen_kmer))&(df_expression_median.index.isin(df_ML.index))]
    print(len(set_non_kmer))
    effect_sizes.append(ime.cohensD(df_expression_median.loc[gen_kmer],set_non_kmer))

#overlap of genes to imeter
overlapp=[]
set_len=[]
gene_file_compare="sorted_imeter_first_intron.csv"

compare=pd.read_csv(gene_file_compare,sep="\t",header=None)
compare=compare[0].map(lambda x: x.split(".")[0])
compare=compare[compare.isin(df_corr.index)]
compare=compare.unique()
for kmer in kmers:
    gen_kmer=ime.find_genes_containing_kmer(kmer,fasta_file)
    gen_kmer=[gen.split(".")[0] for gen in gen_kmer]
    gen_kmer=np.unique(gen_kmer)
    gen_kmer=pd.Series(gen_kmer)
    gen_kmer=gen_kmer[gen_kmer.isin(df_expression.index)]
    n1=len(gen_kmer)
    over=len(set(gen_kmer) & set(compare[:len(gen_kmer)]))
    overlapp.append(over)
    set_len.append(n1)


#Go-term
import os

#five identified consensmotifs and the two IMETER motifs
kmers=["GATTCG","TTTCGA","KCGAGAR","ACYCYRA","ARATCGA","TTNGATYTG","CGATT"]
out_dir=""
all_gens=df_expression.index
for kmer in kmers:
    gens=ime.find_genes_containing_kmer(kmer,fasta_file)
    gens=np.unique([*map(lambda x: x.split("."),gens)])
    gens=df_expression[df_expression.index.isin(gens)].index
    try:
        os.mkdir(out_dir+"/"+kmer)
    except:
            pass
    ime.run_GO_term_analysis(gens,all_gens,kmer,"all_genes_as_comparison",output_direct=out_dir+"/"+kmer,
                             output_name=kmer+"_component",
                             GO_file="ATH_GO_GOSLIM_Feb18.component_slim")
    ime.run_GO_term_analysis(gens,all_gens,kmer,"all_genes_as_comparison",output_direct=out_dir+"/"+kmer,
                             output_name=kmer+"_function",
                             GO_file="ATH_GO_GOSLIM_Feb18.function_slim")
    ime.run_GO_term_analysis(gens,all_gens,kmer,"all_genes_as_comparison",output_direct=out_dir+"/"+kmer,
                             output_name=kmer+"_process",
                             GO_file="ATH_GO_GOSLIM_Feb18.process_slim")




#finding natural mutatation varriants of intron motfis and compare them to the


#preperation of SNP file for comparison
#matrix was created with vcftools
df_SNP=pd.read_csv("first_intron_matrix.csv",sep="\t")



for col in df_SNP.columns[3:]:
    print(col)
    #replace homozygote with single symbol
    hom=df_SNP.loc[:,col].apply(lambda x: x.split("/")[0]==x.split("/")[1])
    df_SNP.loc[hom,col]=df_SNP.loc[hom,col].str.split("/").str[0]
    df_SNP.loc[df_SNP[col]=='.', col] = df_SNP.REF[df_SNP[col]=='.']

df_matrix=pd.DataFrame(np.zeros(np.shape(df_SNP)))

df_matrix.iloc[:,:2]=df_SNP.iloc[:,:2]
df_matrix.columns=df_SNP.columns

for col in df_SNP.columns[3:]:
    print(col)
    uneq=df_SNP[col]!=df_SNP.REF
    df_matrix.loc[uneq,col]=1

#identifying allels with mutations in IME motif :

#each motif: find all genes containing motif and the positions --> masks SNP files
#find which motifs are masked --> find SNP position --> look at matrix at that position to find accession



#modified to return dataframe with all genes with position of respective kmer (empty if none)
#fastafile orientation must be consistent (all from the same strand direction)
def find_kmer_positions_inf(k_mer,fasta_file,cutoff=3):

    replace={"Y":"[CT]","R":"[AG]","S":"[GC]","W":"[AT]","N":"[ATGC]","K":"[GT]","M":"[AC]"}
    pattern=re.compile("|".join(replace.keys()))

    #df["Gene"]=pd.Series(dtype=str)
    reverse_kmer=ime.reverse_complement(k_mer)
    palindrome=reverse_kmer==k_mer

    reverse_kmer=pattern.sub(lambda m: replace[re.escape(m.group(0))], reverse_kmer)
    k_mer=pattern.sub(lambda m: replace[re.escape(m.group(0))], k_mer)

    temp_lis=[]

    with open(fasta_file) as fas:



        if palindrome:
            gene=""
            start=0
            end=0
            for count,line in enumerate(fas):
                if count%2==0:
                    splt=line.split(":")
                    gene=splt[0][1:]
                    start=splt[-1].split("-")[0]
                    end=splt[-1].split("-")[1].strip()
                    chrom=splt[-2]
                    continue
                line=line[cutoff:-cutoff]
                positions=[m.start() for m in re.finditer(k_mer,line)]
                if len(positions)==0:
                    positions=None
                temp_lis.append([gene,chrom,start,end,len(line)+6,positions])
        else:

             for count,line in enumerate(fas):
                if count%2==0:
                    splt=line.split(":")
                    gene=splt[0][1:]
                    start=splt[-1].split("-")[0]
                    end=splt[-1].split("-")[1].strip()
                    chrom=splt[-2]
                    continue
                line=line[cutoff:-cutoff]
                positions=[m.start() for m in re.finditer(k_mer,line)]
                pos_rev=[m.start() for m in re.finditer(reverse_kmer,line)]
                positions.extend(pos_rev)
                positions.sort()
                if len(positions)==0:
                    positions=None
                temp_lis.append([gene,chrom,start,end,len(line)+6,positions])
    df=pd.DataFrame(temp_lis,columns=["gene","chr","start","end","length","positions"])
    fas.close()
    return df


#extract the genes and the respective accessions with mutation at motif
#inputs: summed_frq file containing the gene name, and the SNP with relative position
#fasta, which alos contains the position of the sequence in the title at the end with :xxx-xxx
#lastly the SNP_matrix, not as a file, cause it is large (2GB+), so reading it in evertime the function is called uses too much time
def extract_accessions_with_motif_mutation(kmer,frq_file,fastafile,SNP_matrix):
    kmer_inf=find_kmer_positions_inf(kmer,fastafile).dropna()
    kmer_inf["start"]=kmer_inf["start"].astype(int)
    kmer_inf["chr"]=kmer_inf["chr"].astype(int)
    #the positions are always the first hit of the kmer --> all positions +6 need to be checked for mutation


    #reading in frq and modifying to fit a proper dataframe, might be good idea to do that to all frq in the first place
    frq_df=[]
    with open(frq_file) as frq:
        reader=csv.reader(frq,delimiter="\t")
        next(reader)

        for row in reader:
            tmp=row[:3]
            positions=[x.split(":")[0] for x in row[3:]]
            positions=np.unique(positions)
            tmp.append(positions)
            frq_df.append(tmp)
        frq.close()
    frq_df=pd.DataFrame(frq_df)
    frq_df.columns=["gene","strand","length","positions"]

    SNPs_genes=frq_df.loc[kmer_inf.index]
    #kmer_snps=[any(pk in f for pk in [i for i in[range(1,x+1) for x in k]]) for k,f in zip(kmer_inf.positions,SNPs_genes.positions)]

    #could be done more efficient or compact, but it gets quite confusing, so "inefficient" double for loop
    res=[]
    rel_SNPS=[]
    for k,f in zip(kmer_inf.positions.values,SNPs_genes.positions.values):
        temp_truth=False
        temp_SNPS=[]
        for p in k:
            p=int(p)
            #check overlaps between kmer pos and SNPs
            #check again if not miscalculated
            overlap=[p+1<=int(y) <p+7 for y in f]

            if any(overlap):
                temp_truth=True
                temp_SNPS=np.array(list(map(int,np.array(f)[overlap])))
                break
        res.append(temp_truth)
        rel_SNPS.append(temp_SNPS)
    #all genes with a mutation in
    res=np.array(res)
    rel_SNPS=np.array(rel_SNPS)[res]
    kmer_inf=kmer_inf[res]
    accessions=[]
    for snp,cr in zip(rel_SNPS,kmer_inf.iterrows()):
        snp=np.array(snp)
        #print(cr[0])
        cr=cr[1]
        #print(np.sum(SNP_matrix.POS.isin(snp+cr.start+1)))

        #note to self : not+1 !!! GFF and fasta file count differently
        temp_slice=SNP_matrix[(SNP_matrix.CHROM==cr.chr)&(SNP_matrix.POS.isin(snp+cr.start))]
        #print(np.sum(SNP_matrix.CHROM==cr.chr))
        #print(np.sum(SNP_matrix.POS.isin(snp+cr.start)))]
        accessions.append(temp_slice.columns[3:][(temp_slice.iloc[:,3:]==1).any(axis=0)])
    return([kmer_inf.gene,accessions])


for kmer in kmers:
    genes=ime.find_genes_containing_kmer(kmer,"first_introns.fas")
    genes=[x.split(".")[0] for x in genes]
    genes=np.unique(genes)
    with open(kmer+"_containing_introns.txt","w") as f:
        for gene in genes:
            f.write(gene+"\n")
        f.close()

#1001 genome project expression file
df_expression_1001 = pd.read_csv("GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv",sep="\t",index_col="gene_id")
df_expression_1001 = df_expression_1001.astype(int)

df_matrix.rename(columns={"#CHROM":"CHROM"},inplace=True)

df_expression_1001=np.log(df_expression_1001+0.1)
df_expression_1001.columns=[*map(lambda x: x.replace("X",""),df_expression_1001.columns)]

df_expression_1001_modified=df_expression_1001.drop(df_expression_1001.columns[~np.isin(df_expression_1001.columns,df_matrix.columns)],axis=1)

df_matrix_modified=df_matrix.drop(df_matrix.columns[3:][~np.isin(df_matrix.columns[3:],df_expression_1001.columns)],axis=1)

pval=[]
effect_size=[]
all_genes=[]
for kmer in kmers:
    print("retriving information for "+kmer)
    mutated=extract_accessions_with_motif_mutation(kmer,basic_frq_first,fasta_file,df_matrix_modified)

    cohens=[]
    p_values=[]
    genes=[]
    print("analyzing p_values and effect size")
    for count,i in enumerate(mutated[0]):
        if i.split(".")[0] in df_expression_1001_modified.index:
            #choose first occurrence of gene (isoform) to be kept
            i=i.split(".")[0]
            if i in genes:
                continue

            mut_set=df_expression_1001_modified.loc[i.split(".")[0],mutated[1][count]]
            norm_set=df_expression_1001_modified.drop(mutated[1][count],axis=1).loc[i.split(".")[0]]
            p_values.append(st.ttest_ind(mut_set,norm_set)[1])
            cohens.append(cohensD(norm_set,mut_set))
            genes.append(i)
        else:
             continue
    pval.append(p_values)
    effect_size.append(cohens)
    all_genes.append(genes)

import matplotlib.pyplot as plt
import seaborn as sns



figure,ax=plt.subplots()
sns.violinplot(data=effect_size)
plt.axhline(y = 0, color = 'black', linestyle = '--',linewidth=0.5)
ax.set_xticklabels(kmers,fontsize=8)
plt.ylabel("Effect Size")

import statsmodels.stats.multitest as fdr
fdr=[fdr.multipletests(x,alpha=0.05,method="fdr_bh")[0] for x in pval]

[np.mean(x[y]) for x,y in zip(effect_size,fdr)]
