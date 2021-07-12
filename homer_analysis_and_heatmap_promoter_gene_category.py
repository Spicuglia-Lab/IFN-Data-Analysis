import pandas as pd
import numpy as np
from functools import reduce
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

'''
#HOMER ANALYSIS FOR BACKGROUND HG19 WHOLE GENOME
#Sequences are from same promoter/other promoter and intergenic regions. So whole genome will be background

1. Download chromosome size from ucsc hg19 (http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes)
2. Add column in middle: 0; 
chr1	0	249250621
chr2	0	243199373
3. save it as hg19.genome.bed
'''
#HOMER ANALYSIS
findMotifsGenome.pl ISGF3_peak_Other_Promoter_250bp.sorted.bed hg19 homer_isgf3_peak_other_promoter -size 200 -bg hg19.genome.bed
findMotifsGenome.pl ISGF3_peak_Same_promoter_250bp.sorted.bed hg19 homer_isgf3_peak_same_promoter -size 200 -bg hg19.genome.bed
findMotifsGenome.pl ISGF3_peak_Intergenic_250bp.sorted.bed hg19 homer_isgf3_peak_intergenic -size 200 -bg hg19.genome.bed


#heatmap for the output resluts
def top10_motifs(filename1,filename2,filename3):
    df1 = pd.read_csv(filename1,sep="\t",header=0,usecols=['Motif Name', 'Log P-value'],nrows=10)
    df1['Motif Name'] = df1['Motif Name'].str.split('/').str[0]
    df1['Log P-value'] = -df1['Log P-value']
    #df1['category'] = 'induced_epromoter'
    print(df1)
    df2 = pd.read_csv(filename2,sep="\t",header=0,usecols=['Motif Name', 'Log P-value'],nrows=10)
    df2['Motif Name'] = df2['Motif Name'].str.split('/').str[0]
    df2['Log P-value'] = -df2['Log P-value']
    #df2['category'] = 'induced_genes_epromoters'
    print(df2)
    df3 = pd.read_csv(filename3,sep="\t",header=0,usecols=['Motif Name', 'Log P-value'],nrows=10)
    df3['Motif Name'] = df3['Motif Name'].str.split('/').str[0]
    df3['Log P-value'] = -df3['Log P-value']
    #df3['category'] = 'induced_epromoters'
    print(df3)
    #df = df1.merge(df3,on='Motif Name').merge(df2,on='Motif Name')
    df = reduce(lambda left,right: pd.merge(left,right,how='outer',on='Motif Name'), [df1, df2, df3])
    df.columns = ['Motif Name','intergenic.log10Pval','same_promoter.log10Pval','other_promoter.log10Pval']
    df = df.replace(np.nan, 0)
    #print(df1.head())
    #print(df2.head())
    #print(df3.head())

    #df = df1.append([df2, df3])
    df.to_csv('ISGF3_peak_in_diff_regions_isre_binding_site_in_promoter.tsv',sep="\t",header=True,index=False)
    print(df)
top10_motifs('homer_isgf3_peak_intergenic/knownResults.txt',
             'homer_isgf3_peak_same_promoter/knownResults.txt',
             'homer_isgf3_peak_other_promoter/knownResults.txt')



#RSCRIPT TO GENERATE HEATMAP

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(ggrepel)

data <- read.csv('ISGF3_peak_in_diff_regions_isre_binding_site_in_promoter.tsv',sep='\t',header=T)
data <- data %>%
  pivot_longer(data,cols=2:4,names_to = "Category",values_to = '-log10Pvalue')
data$Category <- factor(data$Category, levels=c("same_promoter.log10Pval",
                                                "other_promoter.log10Pval",
                                                "intergenic.log10Pval"))

p <- ggplot(data = data,aes(x=Category,y=Motif.Name))+
  geom_tile(aes(fill=`-log10Pvalue`),color='white',size=0.3)+
  scale_fill_gradient(low='grey97',high = 'darkred')+
  theme(text = element_text(size = 8),axis.text.x = element_text(angle = 40,hjust=1,vjust=1))+
  xlab("")+
  ylab("")+
  theme(axis.ticks = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        #aspect.ratio = 4/3,
        #plot.title = element_blank())
        plot.title = element_text(hjust=0.5,size = 8, colour = "black"),
        plot.margin=unit(c(0.09,0.09,0.09,0.09),"cm"))

ggsave(p, file="ISGF3_peak_in_diff_regions_isre_binding_site_in_promoter.png", width = 4, height = 3, dpi = 150)
