#Are the Epromoters that bind ISGF3 more highly conserved than the epromoters that don't bind those TF's?
#I don’t understand the question as the majority of Epromoters binds the ISGF3.
#We could however compare the conservation 
    #Between ALL Induced genes with and without ISGF3
    #Between Induced genes and Epromoters
    #Between Epromoters with and without ISGF3

import pandas as pd 
import os
import numpy as np  

def liftover_hg19_to_hg38(bedfile):
    cmd = 'liftOver '+bedfile+' ~/mount/sacapus_remote/Data/liftOver_refseq/hg19ToHg38.over.chain.gz '+bedfile[:-3]+'hg38.bed '+bedfile[:-3]+'hg38.unlifted'
    #print(cmd)
    os.system(cmd)
'''
liftover_hg19_to_hg38('IRF9.STAT1.STAT2.merged.bed')

liftover_hg19_to_hg38('induced_epromoter.2kb.bed')
liftover_hg19_to_hg38('induced_gene_epromoter.2kb.bed')
liftover_hg19_to_hg38('induced_gene.2kb.bed')


liftover_hg19_to_hg38('induced_epromoter.250bp.bed')
liftover_hg19_to_hg38('induced_gene_epromoter.250bp.bed')
liftover_hg19_to_hg38('induced_gene.250bp.bed')
'''

def intersect_epromoter_isgf3(epromoter,peak):
    #INTERSECT INDUCED_EPROMOTER AND ISGF3
    #cmd1 = 'bedtools intersect -wa -wb -a induced_epromoter.2kb.bed -b IRF9.STAT1.STAT2.merged.bed | uniq >induced_epromoter_isgf3.intersect.bed'
    cmd1 = 'bedtools intersect -wa -wb -a '+epromoter+' -b '+peak+' | uniq >'+epromoter[:-3]+'isgf3.intersect.bed'
    os.system(cmd1)
    print(epromoter[:-3]+'isgf3.intersect.bed is generated...')

'''
intersect_epromoter_isgf3('induced_epromoter.250bp.hg38.bed','../IRF9.STAT1.STAT2.merged.hg38.bed')
intersect_epromoter_isgf3('induced_gene.250bp.hg38.bed','../IRF9.STAT1.STAT2.merged.hg38.bed')
intersect_epromoter_isgf3('induced_gene_epromoter.250bp.hg38.bed','../IRF9.STAT1.STAT2.merged.hg38.bed')

intersect_epromoter_isgf3('induced_epromoter.2kb.hg38.bed','../IRF9.STAT1.STAT2.merged.hg38.bed')
intersect_epromoter_isgf3('induced_gene.2kb.hg38.bed','../IRF9.STAT1.STAT2.merged.hg38.bed')
intersect_epromoter_isgf3('induced_gene_epromoter.2kb.hg38.bed','../IRF9.STAT1.STAT2.merged.hg38.bed')
'''

def empromoter_isgf3_binding(epromoter,epromoter_with_isgf3):
    epr = pd.read_csv(epromoter,sep="\t",header=None)
    epr_with_isgf3 = pd.read_csv(epromoter_with_isgf3,sep="\t",header=None,usecols=[0,1,2])
    epr_vs_epr_with_isgf3 = epr.merge(epr_with_isgf3,
                              indicator=True,
                              how='outer')
    epr_vs_epr_with_isgf3['TF_binding']= np.where(epr_vs_epr_with_isgf3['_merge']=='both','ISGF3','no_binding')
    epr_vs_epr_with_isgf3.drop_duplicates(subset =[0,1],keep = 'first', inplace = True)
    epr_vs_epr_with_isgf3[[0,1,2,'TF_binding']].to_csv(epromoter[:-3]+'isgf3.binding.bed',sep="\t",header=None,index=False)
    print(epromoter[:-3]+'isgf3.binding.bed is generated...')
'''
empromoter_isgf3_binding('induced_epromoter.2kb.hg38.bed','induced_epromoter.2kb.hg38.isgf3.intersect.bed')
empromoter_isgf3_binding('induced_gene.2kb.hg38.bed','induced_gene.2kb.hg38.isgf3.intersect.bed')
empromoter_isgf3_binding('induced_gene_epromoter.2kb.hg38.bed','induced_gene_epromoter.2kb.hg38.isgf3.intersect.bed')

empromoter_isgf3_binding('induced_epromoter.250bp.hg38.bed','induced_epromoter.250bp.hg38.isgf3.intersect.bed')
empromoter_isgf3_binding('induced_gene.250bp.hg38.bed','induced_gene.250bp.hg38.isgf3.intersect.bed')
empromoter_isgf3_binding('induced_gene_epromoter.250bp.hg38.bed','induced_gene_epromoter.250bp.hg38.isgf3.intersect.bed')
'''
