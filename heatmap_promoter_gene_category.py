import pandas as pd
import numpy as np
from functools import reduce

def top10_motifs(filename1,filename2,filename3):
    df1 = pd.read_csv(filename1,sep="\t",header=0,usecols=['Motif Name', 'Log P-value'],nrows=10)
    df1['Motif Name'] = df1['Motif Name'].str.split('/').str[0]
    df1['Log P-value'] = -df1['Log P-value']
    #df1['category'] = 'induced_epromoter'

    df2 = pd.read_csv(filename2,sep="\t",header=0,usecols=['Motif Name', 'Log P-value'],nrows=10)
    df2['Motif Name'] = df2['Motif Name'].str.split('/').str[0]
    df2['Log P-value'] = -df2['Log P-value']
    #df2['category'] = 'induced_genes_epromoters'

    df3 = pd.read_csv(filename3,sep="\t",header=0,usecols=['Motif Name', 'Log P-value'],nrows=10)
    df3['Motif Name'] = df3['Motif Name'].str.split('/').str[0]
    df3['Log P-value'] = -df3['Log P-value']
    #df3['category'] = 'induced_epromoters'

    #df = df1.merge(df3,on='Motif Name').merge(df2,on='Motif Name')
    df = reduce(lambda left,right: pd.merge(left,right,how='outer',on='Motif Name'), [df1, df2, df3])
    df.columns = ['Motif Name','induced_genes.log10Pval','induced_genes_epromoters.log10Pval','induced_epromoters.log10Pval']
    df = df.replace(np.nan, 0)
    #print(df1.head())
    #print(df2.head())
    #print(df3.head())

    #df = df1.append([df2, df3])
    df.to_csv('IFN_three_gene_category.data.tsv',sep="\t",header=True,index=False)
    print(df)
top10_motifs('induced_genes/knownResults.txt',
             'induced_genes_epromoters/knownResults.txt',
             'induced_epromoters/knownResults.txt')
