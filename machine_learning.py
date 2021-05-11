#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""



import pandas as pd
import csv
import numpy as np
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier

import shap
from sklearn.datasets.samples_generator import make_blobs




import itertools
import matplotlib.pyplot as plt



from matplotlib.backends.backend_pdf import PdfPages

import IME



def count_A_T_G_C(line):
    return [line.count("A"),line.count("T"),line.count("G"),line.count("C")]

genes = []
with open("gene_list") as f:
    genes = f.read().splitlines()


#import information from gff file
df_gff = pd.read_csv("first_intron.gff",
                   sep = "\t",header=None)
df_gff[8] = df_gff[8].apply(lambda x: x.split("=")[1])
df_gff[8] = df_gff[8].apply(lambda x: x.split(";")[0])
df_gff_reduced = df_gff[[0, 8]]
df_gff_reduced.columns=["Chrom","Gene"]
df_gff_reduced["Gene"] = df_gff_reduced["Gene"].apply(lambda x:x.split(".")[0])
df_gff_reduced = df_gff_reduced.drop_duplicates(subset = "Gene")

df_gff_reduced["length"]=df_gff.iloc[df_gff_reduced.index,:][4]-df_gff.iloc[df_gff_reduced.index,:][3]+1

df_gff_reduced = df_gff_reduced[df_gff_reduced["Chrom"].str.isnumeric()]
df_gff_reduced=df_gff_reduced.drop("Chrom",1)


#5prime UTR
df_5prime = pd.read_csv("TAIR10_5_prime.gff",sep="\t",header=None)
df_5prime[9] = df_5prime[4]-df_5prime[3]+1

mRNA_modified = pd.read_csv("mRNA.gff",sep="\t",header=None)
mRNA_modified = mRNA_modified[mRNA_modified[8].isin(df_gff[8])]
mRNA_modified.index = df_gff.index

prime = []

for gen in df_gff.loc[df_gff_reduced.index][8]:
    temp_df=df_5prime[df_5prime[8]==gen]
    if len(temp_df)==0:
        if mRNA_modified[mRNA_modified[8]==gen].iloc[0,6]=="+":
            prime.append([0,df_gff[df_gff[8]==gen].iloc[0,3]-mRNA_modified[mRNA_modified[8]==gen].iloc[0,3]])
            continue
        prime.append([0,mRNA_modified[mRNA_modified[8]==gen].iloc[0,4]-df_gff[df_gff[8]==gen].iloc[0,4]])
        continue
    #all required infromation regarding intron/utr borders
    utr_3 = temp_df.iloc[len(temp_df)-1,3]
    utr_5 = temp_df.iloc[len(temp_df)-1,4]
    intron_3 = df_gff[df_gff[8]==gen].iloc[0,3]
    intron_5 = df_gff[df_gff[8]==gen].iloc[0,4]
    if temp_df.iloc[0,6]=="+":
        if intron_5<utr_5:
            prime.append([temp_df[9].sum(),intron_5-utr_5])
        else:
            prime.append([temp_df[9].sum(),intron_3-utr_5])
    else:
        if intron_3>utr_3:
            prime.append([temp_df[9].sum(),utr_3-intron_3])
        else:
            prime.append([temp_df[9].sum(),utr_3-intron_5])

del([gen,utr_3,utr_5,intron_3,intron_5,temp_df])

df_gff_reduced["distance_CDS"] = [x[1] for x in prime]



#----------------------------------imeter score

df_imeter=pd.read_csv("imeter.csv",sep="\t",header=None)

#select imeter2 score
df_gff_reduced["imeter"]=df_imeter.loc[df_gff_reduced.index,2]


#find all UTRs
df_5prime = df_5prime[df_5prime[8].isin(df_gff.loc[df_gff_reduced.index,8])]

SNP_number = []
with open("frq_summed") as frq:
    next(frq)
    reader=csv.reader(frq,delimiter="\t")
    for row in reader:
        SNP_number.append(len(row)-3)
SNP_number = np.array(SNP_number)
#df_gff_reduced=df_gff_reduced[df_gff_reduced["Chrom"].apply(lambda x:x.isnumeric())]
df_gff_reduced["SNP_per_bp"] = SNP_number[df_gff_reduced.index]
df_gff_reduced["SNP_per_bp"] = df_gff_reduced["SNP_per_bp"]/df_gff_reduced["length"]

bases=["A","T","G","C"]
kmer=[''.join(p) for p in itertools.product(bases,repeat=2)]
reverse_kmer=list(map(IME.reverse_complement,kmer))
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



twomer_list=[]
kmer=list(redundant.keys())
ATGC=[]
with open("correct_strand_fasta") as fas:
    for count,line in enumerate(fas):
        if count%2==0:
            continue
        line=line[3:-3]
        twomer_list.append([line.count(mer)+line.count(redundant[mer]) for mer in kmer])
        ATGC.append(np.array(count_A_T_G_C(line))/len(line))

twomer_list=pd.DataFrame(twomer_list)
twomer_list.columns=kmer


df_gff_reduced=pd.concat([df_gff_reduced,twomer_list.loc[df_gff_reduced.index]],axis=1)

for mer in kmer:
    df_gff_reduced[mer]=df_gff_reduced[mer]/(df_gff_reduced["length"]-6)




#A T G C content

ATGC=pd.DataFrame(ATGC)
ATGC.columns=["A","T","G","C"]
df_gff_reduced=pd.concat([df_gff_reduced,ATGC.loc[df_gff_reduced.index]],axis=1)


#inserting distance to the TSS for all introns in the dataframe. first + and - distances are calculated. Then (in a very convoluted way),it
df_gff_reduced["distance_TSS"] = 0
plus = df_gff.loc[df_gff_reduced.index,3]-mRNA_modified.loc[df_gff_reduced.index,3]
minus = mRNA_modified.loc[df_gff_reduced.index,4]-df_gff.loc[df_gff_reduced.index,4]
sign = df_gff.loc[df_gff_reduced.index,6]=="+"

df_gff_reduced.loc[sign,"distance_TSS"] = plus.loc[sign]
df_gff_reduced.loc[~sign,"distance_TSS"] = minus.loc[~sign]

df_genes=pd.read_csv("genes.gff",
                     sep = "\t",header = None)
df_genes[8] = df_genes[8].apply(lambda x: (x.split("=")[1]).split(";")[0])
df_genes = df_genes[df_genes[8].isin(df_gff_reduced["Gene"])]
#df_gff = df_gff[df_gff[8].isin(df_gff_reduced["Gene"])]
#df_gff = df_gff.drop_duplicates(subset = 8)






#----------------------------------------------------------------------------
#methylation
methylation_C = np.array(IME.differential_methylated_fasta("C-DMR.csv","first_introns.fas"))
methylation_CG = np.array(IME.differential_methylated_fasta("CG-DMR.csv","first_introns.fas"))

df_gff_reduced["methylation_C"] = [methylation_C[x] for x in df_gff_reduced.index]
df_gff_reduced["methylation_CG"] = [methylation_CG[x] for x in df_gff_reduced.index]


all_introns = pd.read_csv("introns.gff",sep="\t",header=None)
all_introns[8] = all_introns[8].apply(lambda x: x.split("=")[1].split(";")[0])
all_introns=all_introns[8].value_counts()


#removing the fields with unknown UTR length, since those are in general lower expressed --> false signal
df_gff_reduced["5UTR_length"] = [x[0] for x in prime]
df_gff_reduced=df_gff_reduced.loc[df_gff_reduced["5UTR_length"]>0,:]


df_transpose = pd.read_csv("TAIR10_Transposable_Elements.txt",
                   sep = "\t")

df_transpose["Chrom"]=df_transpose.Transposon_Name.apply(lambda x: x[2])
df_transpose.Chrom=df_transpose.Chrom.astype(int)

df_gff=df_gff[df_gff[0].str.isnumeric()]
df_gff[0]=df_gff[0].astype(int)

def number_of_trans_elements(row,dataframe):
    return dataframe[(dataframe.Chrom==row[0])&(row[3]<dataframe["Transposon_min_Start"])&(row[4]>dataframe["Transposon_min_Start"])].shape[0]
test=df_gff.apply(lambda x:number_of_trans_elements(x,df_transpose),axis=1)

df_gff_reduced["n_transposon"]=test.loc[df_gff.index]


#------------------------------------------intron retainment
df_IR=pd.read_csv("IR_first_introns.gff",sep="\t",header=None)
df_IR[8]=df_IR[8].apply(lambda x: x.split("=")[1].split(";")[0])
#get all indices of introns which are retained
retained=df_gff[df_gff[8].isin(df_IR[8])].index

#put it into the dataframe containing all features
df_gff_reduced["IR"]=0
df_gff_reduced.loc[df_gff_reduced.index.isin(retained),"IR"]=1


#CNS frature

def extract_relevant_CNS(CNS_file="AllFootPrintsFDR0.10_scores.bed",gff_file="first_intron.gff",output_file=""):


    CNS=pd.read_csv(CNS_file,sep="\t",skiprows=1,header=None)

    #only keep entries whcih conserved in 4 species
    CNS=CNS[CNS[4]>3]

    CNS_position=[]
    with open(gff_file) as fasta:
        reader=csv.reader(fasta,delimiter="\t")
        for line in reader:
            chrom=int(line[0])
            start=int(line[3])
            end=int(line[4])
# =============================================================================
#
#             #get all relevant positions
#             positions=(list(CNS[(CNS[0]==chrom)&(CNS[1]>=start)&(CNS[1]<=end)][1]))
#             #calculate the positions in relation to the
#             #CNS_position.append([int(x)-start+1 for x in positions])
#
# =============================================================================
            a=len(CNS[CNS[0]==chrom][(CNS.loc[CNS[0]==chrom,1]<=start)&(CNS.loc[CNS[0]==chrom,2]>start)])
            b=len(CNS[CNS[0]==chrom][(CNS.loc[CNS[0]==chrom,1]>start)&(CNS.loc[CNS[0]==chrom,1]<end)])
            CNS_position.append(a+b)

            #return number of CNS
        if output_file!="":

            with open(output_file,"wt") as out:
                writer=csv.writer(out,delimiter="\t")
                for x in CNS_position:
                    writer.writerow(x)
                out.close()
    return(pd.DataFrame(CNS_position))

CNS=extract_relevant_CNS()
df_gff_reduced["CNS"]=CNS.loc[df_gff_reduced.index]

#include minimum folding energy as feature

min_f=pd.read_csv("energy.txt",sep="\t",header=None)
df_gff_reduced["min_folding_energy"]=1000

#looks convoluted, but works
#note min_f is shorter than the df_gff, since a few introns got ignored by mfold due to them beeing to long, therefore index of df_gff != index of min_f
min_f.index=df_gff[df_gff[8].isin(min_f[1])].index
#all df_gff_reduced that have an entry for min_f get the appropriate entry
df_gff_reduced.loc[min_f.loc[df_gff_reduced.index,1].index,"min_folding_energy"]=min_f.loc[df_gff_reduced.index,0]

df_gff_reduced=df_gff_reduced[df_gff_reduced.min_folding_energy!=1000]
df_gff_reduced=df_gff_reduced.dropna()


#gene expression
RNA_file=""
gene_list=[]

with open(RNA_file,mode="rt") as RNA:
    reader=csv.reader(RNA,delimiter="\t")
    for row in reader :
        if row[0] in genes:
            gene_list.append(row)
    RNA.close()


df_expression = pd.DataFrame(gene_list)
df_expression.index=df_expression[0]
df_expression = df_expression.drop(columns=0)
df_expression = df_expression.astype(float)
gene_list=None
df_expression = df_expression.median(1)
df_expression = df_expression[df_expression.index.isin(df_gff_reduced["Gene"])]
df_expression = df_expression.sort_index(0,level=df_gff_reduced["Gene"])

#removing too low expressed samples
df_expression = df_expression[df_expression>0.1]

df_gff_reduced_chip = df_gff_reduced[df_gff_reduced.Gene.isin(df_expression.index)]
#remove non intron features
df_gff_reduced_chip.index=df_gff_reduced_chip.Gene
df_gff_reduced_chip=df_gff_reduced_chip.drop(columns=["Gene","5UTR_length"])
df_gff_reduced_chip["min_folding_energy"]=df_gff_reduced_chip["min_folding_energy"]/df_gff_reduced_chip["length"]


#---------------------Binning in two classes
bin_class=pd.qcut(df_expression[df_gff_reduced_chip.index],2,labels=["low","high"])

df_test = df_gff_reduced_chip.sample(frac=0.3)
df_train = df_gff_reduced_chip[~(df_gff_reduced_chip.index.isin(df_test.index))]


#------------------------low vs high

bin_quart=pd.qcut(df_expression[df_gff_reduced_chip.index],4,labels=["low","m1","m2","high"])
low_high_micro=np.array((bin_quart=="low")|(bin_quart=="high"))
bin_lh=bin_quart[low_high_micro]
df_low_high=df_gff_reduced_chip[low_high_micro]



df_test_lh = df_low_high.sample(frac=0.3)
df_train_lh = df_low_high[~df_low_high.index.isin(list(df_test_lh.index))]


rfc=RandomForestClassifier(bootstrap=True, ccp_alpha=0.0, class_weight=None,
                          criterion='gini', max_depth=90, max_features='sqrt',
                          max_leaf_nodes=None, max_samples=None,
                          min_impurity_decrease=0.0, min_impurity_split=None,
                          min_samples_leaf=1, min_samples_split=5,
                          min_weight_fraction_leaf=0.0, n_estimators=1200,
                          n_jobs=None, oob_score=False, random_state=None,
                          verbose=0, warm_start=False)

rfc.fit(df_train_lh,bin_lh[df_train_lh.index])

explainer = shap.TreeExplainer(rfc)

shap_values = explainer.shap_values(df_train_lh)
