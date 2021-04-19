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

