import pandas as pd
import numpy as np
import os
import subprocess
import datetime
from itertools import chain
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

start_time = datetime.datetime.now()
print("\nThe job started at: ",start_time,"\n")

def read_expression_and_rnaseq_capstarseq(expression,rnaseq_capstarseq):
    induced_genes = pd.read_csv(expression,header=0,sep="\t")
    rnaseq_capstarseq = pd.read_csv(rnaseq_capstarseq,header=0,sep="\t",
                                    usecols=['gene','DESeq2_log2FoldChange','DESeq2_padj','IFN_vs_NS_fold_change.log2FC','IFN_vs_NS.plot_condition'])
    rnaseq_capstarseq['IFN Stimulation'] = np.where((rnaseq_capstarseq['IFN_vs_NS.plot_condition'] == 'Never Epromoter')|(rnaseq_capstarseq['IFN_vs_NS.plot_condition'] == 'IFN-'),
                                                    'Non-stimulated Epromoter','IFN-stimulated Eprmoter')
    rnaseq_capstarseq['Gene Regulation'] = np.where((rnaseq_capstarseq['DESeq2_log2FoldChange'] > 1) & (rnaseq_capstarseq['DESeq2_padj'] < 0.001),'Induced',
                                                     np.where((rnaseq_capstarseq['DESeq2_log2FoldChange'] < -1) & (rnaseq_capstarseq['DESeq2_padj'] < 0.001),'Repressed','Non-induced'))
    #print(rnaseq_capstarseq.head(4));exit(0)
    rnaseq_capstarseq = rnaseq_capstarseq[~(rnaseq_capstarseq['Gene Regulation']=='Non-induced') & ~(rnaseq_capstarseq['IFN Stimulation']== 'Non-stimulated Epromoter')]
    rnaseq_capstarseq = rnaseq_capstarseq.drop_duplicates('gene',keep=False)
    rnaseq_capstarseq.to_csv('plot_data/rnaseq_capstarseq_induced_gene_or_epromoter.tsv',sep="\t",header=True,index=False)

    induced_genes_and_epromoter = pd.merge(induced_genes,rnaseq_capstarseq,how='outer',on=['gene','DESeq2_log2FoldChange'])
    induced_genes_and_epromoter = induced_genes_and_epromoter.drop_duplicates('gene',keep='last')
    induced_genes_and_epromoter['IFN Stimulation'] = induced_genes_and_epromoter['IFN Stimulation'].replace(to_replace=np.nan, value='Non-stimulated Epromoter')
    induced_genes_and_epromoter['Gene Regulation'] = np.where(induced_genes_and_epromoter['DESeq2_log2FoldChange']>1,'Induced',
                                                              np.where(induced_genes_and_epromoter['DESeq2_log2FoldChange']<-1,'Repressed','Non-induced'))
    induced_genes_and_epromoter.to_csv('plot_data/gene_cluster/induced_genes_epromoter_expression.tsv',sep="\t",header=True,index=False)
    return(induced_genes_and_epromoter)

class Extract_induced_genes_coordinates:
    def __init__(self):
        '''Read File'''
    def induced_genes(self,coordinatefile,genelist):
        print("##### Extraction induced genes coordinates from RefSeq Database file: 'RefGene.txt' ...")
        self.coordinatefile = coordinatefile
        self.genelist = genelist

        induced_gene_data = pd.read_csv(self.genelist,sep="\t",header=0)

        refgene_col = ["bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames"]
        refgene_data = pd.read_csv(self.coordinatefile, sep="\t", header=None,names=refgene_col)

        refgene_data = refgene_data[~refgene_data['chrom'].str.contains('_')]

        refgene_data['txn_length'] = abs(refgene_data['txStart']-refgene_data['txEnd'])
        refgene_data_txn_length = refgene_data.groupby('name2', as_index=False).apply(pd.DataFrame.sort_values, 'txn_length', ascending = False).reset_index(drop=True)

        induced_refgene_data = refgene_data_txn_length[refgene_data_txn_length['name2'].isin(induced_gene_data['gene'])].copy()
        induced_refgene_data['tss'] = induced_refgene_data.groupby(['name2']).cumcount()+1
        induced_gene_tss = pd.DataFrame(induced_refgene_data, columns=['name','chrom','strand','txStart','txEnd','name2','tss'])
        #print(induced_gene_tss);exit(0)
        print("  Induced genes tss coordinates are saved in variable: 'induced_gene_tss'\n")
        induced_gene_tss.to_csv(cwd+'/plot_data/gene_cluster/induced_gene_tss',sep="\t",index=False)

        return induced_gene_tss

class Inducedgenesclustering:
    def __init__(self):
        '''Clustering on the basis of 5' coordinate of induced genes'''
    def induced_genes_clustering(self, data):
        print("##### Induced genes clustering...")
        self.data = data

        os.chdir(cwd)
        five_prime_coordinate_gene = pd.DataFrame()

        dtypes = {'chr':'str','strand':'str','5p_gene1':'int','5p_gene2':'int','3p_gene':'int','gene':'str'}
        five_prime_coordinate_gene = pd.DataFrame(columns=['chr','5p_gene1','5p_gene2','3p_gene','gene'])#gene coordinate

        for index, row in self.data.iterrows():
            if(row['strand'] == '+' and row['tss'] == 1):
                five_prime_coordinate_gene=five_prime_coordinate_gene.append(pd.Series([row['chrom'],row['txStart'],
                row['txStart'],row['txEnd'],row['name2']],index=['chr','5p_gene1','5p_gene2','3p_gene','gene']),ignore_index=True)

            elif(row['strand'] == '-' and row['tss']==1):
                five_prime_coordinate_gene=five_prime_coordinate_gene.append(pd.Series([row['chrom'],row['txEnd'],row['txEnd'],
                row['txStart'],row['name2']],index=['chr','5p_gene1','5p_gene2','3p_gene','gene']),ignore_index=True)

        sort_five_prime_coordinate_gene = five_prime_coordinate_gene.sort_values(['chr','5p_gene1','5p_gene2'])
        five_prime_coordinate = sort_five_prime_coordinate_gene[['chr','5p_gene1','5p_gene2','gene']]
        five_prime_coordinate.to_csv('plot_data/gene_cluster/induced_genes_5p_coordinates.bed',sep="\t",index=False,header=None)
        sort_five_prime_coordinate_gene.to_csv('plot_data/gene_cluster/induced_genes_5p_3p_coordinates.bed',sep='\t',index=False,header=None)

        cmd_cluster_5p_genes = "bedtools cluster -d 100000 -i plot_data/gene_cluster/induced_genes_5p_3p_coordinates.bed > plot_data/gene_cluster/induced_genes_bed_tool_clust.bed"
        subprocess.call(cmd_cluster_5p_genes,shell=True)
        #count cluster which has more than one genes
        cmd_count_cluster= "awk 'NR == FNR {CNT[$NF]++;if (!($NF in order) && CNT[$NF]>1) order[$NF]=++cnt;next} $NF in order {print $0, order[$NF]}' OFS='\t' plot_data/gene_cluster/induced_genes_bed_tool_clust.bed plot_data/gene_cluster/induced_genes_bed_tool_clust.bed > plot_data/gene_cluster/induced_genes_final_cluster.bed"
        subprocess.call(cmd_count_cluster, shell=True)

        print("  Induced Gene Cluster file is generated:\n\t"+cwd+"induced_genes_final_cluster.bed\n")
        return sort_five_prime_coordinate_gene

class Genespromotercluster:
    def __init__(self):
        '''Gene Promoter Cluster Table'''

    def genes_promoter_cluster(self,filename):
        print("#####Genes promoter cluster...")
        self.filename = filename
        gene_5p_cluster_file_name = self.filename

        columns = ['chr', 'gene_5p', 'gene_5p2', 'gene_3p', 'gene_name', 'clust1', 'clust2']
        gene_cluster_data = pd.read_csv(gene_5p_cluster_file_name, sep="\t", names=columns)

        cluster_number = pd.Series(gene_cluster_data['clust2'].unique())
        gene_promoter_cluster_table_columns = ['cluster number', 'number of genes per cluster', 'number of promoters per cluster', 'flag','genes in cluster']
        gene_promoter_cluster_table = pd.DataFrame()

        for i in cluster_number:
            number_of_genes = []
            cluster_genes = []
            for j in range(0, len(gene_cluster_data)-1):
                if(gene_cluster_data.iloc[j, 6] == i):
                    k = j + 1
                    if(gene_cluster_data.iloc[k,6] == i):
                        dist = abs(gene_cluster_data.iloc[j, 1]-gene_cluster_data.iloc[k, 1])
                        number_of_genes.append(dist)
                    cluster_genes.append(gene_cluster_data.iloc[k-1, 4])
                    if(k == len(gene_cluster_data)-1):
                        cluster_genes.append(gene_cluster_data.iloc[k, 4])
            number_of_genes_per_cluster = len(number_of_genes)+1
            number_of_promoters = sum(i > 1000 for i in number_of_genes)+1
            flag = 'U' if number_of_promoters > 1 else 'F'
            gene_promoter_cluster_table=gene_promoter_cluster_table.append(pd.Series([i,number_of_genes_per_cluster, number_of_promoters, flag, ', '.join(cluster_genes)], index=gene_promoter_cluster_table_columns), ignore_index=True)

        gene_promoter_cluster_table['cluster number'] = gene_promoter_cluster_table['cluster number'].astype(int)
        gene_promoter_cluster_table['number of genes per cluster'] = gene_promoter_cluster_table['number of genes per cluster'].astype(int)
        gene_promoter_cluster_table['number of promoters per cluster'] = gene_promoter_cluster_table['number of promoters per cluster'].astype(int)
        gene_promoter_cluster_table['flag'] = gene_promoter_cluster_table['flag'].astype(str)
        gene_promoter_cluster_table['genes in cluster'] = gene_promoter_cluster_table['genes in cluster'].astype(str)
        gene_promoter_cluster_table = gene_promoter_cluster_table.reindex(columns=gene_promoter_cluster_table_columns)

        gene_promoter_cluster_table.to_csv(cwd+'/plot_data/gene_cluster/gene_promoter_cluster.table', sep='\t',index=False)
        print("  The Gene Promoter Cluster Table is generated:\n\t",cwd+"/plot_data/gene_cluster/gene_promoter_cluster.table\n")
        return gene_promoter_cluster_table

cwd = os.getcwd()

induced_genes_epromoter = read_expression_and_rnaseq_capstarseq('data/log2FC_1_padj_0.001_k562_IFNa.induced_genes',
                                                                'data/plot_data_K562_RNAseq_Capstarseq.hg19.tsv')

extract_induced_genes=Extract_induced_genes_coordinates()
induced_genes_tss = extract_induced_genes.induced_genes('data/hg19.refGene.txt',
                                                        'plot_data/rnaseq_capstarseq_induced_gene_or_epromoter.tsv')

induced_genes_clust = Inducedgenesclustering()
induced_genes_cluster=induced_genes_clust.induced_genes_clustering(induced_genes_tss)

genes_promoter_cluster = Genespromotercluster()
genes_promoter_clust = genes_promoter_cluster.genes_promoter_cluster(cwd+'/plot_data/gene_cluster/induced_genes_final_cluster.bed')

def read_induced_gene_epromoter(filename):
    induced_gene_epromoter = pd.read_csv(filename,header=0,sep="\t")
    #print(induced_gene_epromoter.head(2));exit(0)
    return(induced_gene_epromoter)

def read_cluster_file(filename):
    cluster = pd.read_csv(filename,header=0,sep="\t",
                          usecols = ['number of genes per cluster','cluster number','genes in cluster'])
    cluster['cluster_number'] = 'clust '+cluster['cluster number'].astype(str)+': ['+cluster['genes in cluster']+']'
    cluster = pd.concat([pd.Series(row['cluster_number'],
                                   row['genes in cluster'].split(', ')) for _, row in cluster.iterrows()]).reset_index()
    cluster.columns = ['gene','cluster number']
    #print(cluster.head(2));exit(0)
    return(cluster)

def merge_rnaseq_capstarseq_cluster(df1,df2):
    df = pd.merge(df1,df2,on='gene',how='outer')
    df = df[pd.notnull(df['cluster number'])]
    df['IFN Stimulation'].fillna('Non-stimulated Epromoter', inplace=True)
    df['Gene Regulation'].fillna('Non-induced', inplace=True)
    #print(df.head(10));exit(0)
    return(df)


rnaseq_capstarseq = read_induced_gene_epromoter('plot_data/induced_genes_epromoter_expression.tsv')
cluster = read_cluster_file('plot_data/gene_cluster/gene_promoter_cluster.table')

df = merge_rnaseq_capstarseq_cluster(cluster,rnaseq_capstarseq)
df['noninduced_epromoter_induced_gene'] = np.where((df['IFN Stimulation'] == 'Non-stimulated Epromoter') & (df['Gene Regulation'] == 'Induced'),1,0)
df['induced_epromoter_induced_gene'] = np.where((df['IFN Stimulation'] == 'IFN-stimulated Eprmoter') & (df['Gene Regulation'] == 'Induced'),1,0)
df['induced_epromoter_noninduced_gene'] = np.where((df['IFN Stimulation'] == 'IFN-stimulated Eprmoter') & (df['Gene Regulation'] == 'Repressed'),1,0)
df = df[['cluster number','noninduced_epromoter_induced_gene','induced_epromoter_induced_gene','induced_epromoter_noninduced_gene']]
df = df.groupby(['cluster number']).sum().reset_index()

df = df.sort_values(by=['cluster number'], ascending=True)
df.to_csv('plot_data/heatmap_data.tsv',sep="\t",index=False,header=True)

df.set_index('cluster number', inplace=True)

df = pd.read_csv('plot_data/heatmap_data.tsv',sep="\t",header=0)
#print(df.head(4));exit(0)
df1 = df[['noninduced_epromoter_induced_gene']]
df2 = df[['induced_epromoter_induced_gene']]
df3 = df[['induced_epromoter_noninduced_gene']]

fig, (ax1,ax2,ax3) = plt.subplots(ncols=3)
fig.set_size_inches(10.0, 18.27)
fig.subplots_adjust(wspace=0.0)
sns.heatmap(df1, cmap="Blues",annot=True,fmt='g', annot_kws={'size':10},ax=ax1, cbar=False, linewidths=.5,robust=True)
ax1.set_xticklabels(ax1.get_xticklabels(),rotation=30,horizontalalignment='right')
ax1.set_ylabel('')
sns.heatmap(df2, cmap="Greens",annot=True,fmt='g', annot_kws={'size':10},ax=ax2, yticklabels=[], cbar=False, linewidths=.5,robust=True)
ax2.set_xticklabels(ax2.get_xticklabels(),rotation=30,horizontalalignment='right')
ax2.set_ylabel('')
sns.heatmap(df3, cmap="Oranges",annot=True,fmt='g', annot_kws={'size':10},ax=ax3, yticklabels=[], cbar=False, linewidths=.5,robust=True)
ax3.set_xticklabels(ax3.get_xticklabels(),rotation=30,horizontalalignment='right')
ax3.set_ylabel('')
plt.show()
exit(0)
fig.savefig("output.png")




print("The job is completed at: ",datetime.datetime.now())
print("Total time cost: ",datetime.datetime.now()-start_time)
exit(0)



















def read_rnaseq_capstarseq(filename):
    rnaseq_capstarseq1 = pd.read_csv(filename,header=0,sep="\t",
                                    usecols=['gene','DESeq2_log2FoldChange','DESeq2_padj','IFN_vs_NS_fold_change.log2FC','IFN_vs_NS.plot_condition'])
    rnaseq_capstarseq = rnaseq_capstarseq1.drop_duplicates(subset='gene', keep="first")
    #rnaseq_capstarseq.to_csv('rnaseq_capstarseq.tsv',sep="\t",header=True,index=False);exit(0)
    rnaseq_capstarseq['IFN Stimulation'] = np.where((rnaseq_capstarseq['IFN_vs_NS.plot_condition'] == 'Never Epromoter')|(rnaseq_capstarseq['IFN_vs_NS.plot_condition'] == 'IFN-'),'Non-stimulated Epromoter','IFN-stimulated Eprmoter')
    rnaseq_capstarseq['Gene Regulation'] = np.where((rnaseq_capstarseq['DESeq2_log2FoldChange'] >1) & (rnaseq_capstarseq['DESeq2_padj'] < 0.001),'Induced','Non-induced')
    induced_genes = rnaseq_capstarseq[rnaseq_capstarseq['Gene Regulation'] == 'Induced']
    induced_epromoter = rnaseq_capstarseq[rnaseq_capstarseq['IFN Stimulation'] == 'IFN-stimulated Eprmoter']
    #rnaseq_capstarseq.to_csv('ranseq_capstarseq.tsv',sep="\t",header=True,index=False);exit(0)

    #induced_genes.to_csv('induced_genes.tsv',sep="\t",header=True,index=False)
    #induced_epromoter.to_csv('induced_epromoter.tsv',sep="\t",header=True,index=False);exit(0)
    induced_genes_or_induced_epromoter = induced_genes.append(induced_epromoter,ignore_index=True)
    induced_genes_or_induced_epromoter.drop_duplicates()#data.drop_duplicates(keep=False,inplace=True)
    induced_genes_or_induced_epromoter.to_csv('induced_genes_or_induced_epromoter.tsv',sep="\t",header=True,index=False);exit(0)
    print(induced_genes_or_induced_epromoter);exit(0)

    #print(rnaseq_capstarseq.head(2))
    return(rnaseq_capstarseq)

def read_cluster_file(filename):
    cluster = pd.read_csv(filename,header=0,sep="\t",
                          usecols = ['number of genes per cluster','cluster number','genes in cluster'])
    cluster['cluster_number'] = 'clust '+cluster['cluster number'].astype(str)+': ['+cluster['genes in cluster']+']'
    #cluster['cluster_number'] = cluster['cluster number'].astype(str)+' ('+cluster['number of genes per cluster'].astype(str)+')'
    #print(cluster);exit(0)
    cluster = pd.concat([pd.Series(row['cluster_number'],
                                   row['genes in cluster'].split(', ')) for _, row in cluster.iterrows()]).reset_index()
    cluster.columns = ['gene','cluster number']
    #print(cluster.head(2));exit(0)
    return(cluster)

def merge_rnaseq_capstarseq_cluster(df1,df2):
    df = pd.merge(df1,df2,on='gene',how='outer')
    df = df[pd.notnull(df['cluster number'])]
    df['IFN Stimulation'].fillna('No Epromoter', inplace=True)
    df['Gene Regulation'].fillna('Non-induced', inplace=True)
    return(df)

#def heatmap(df):

rnaseq_capstarseq = read_rnaseq_capstarseq('data/plot_data_K562_RNAseq_Capstarseq.hg19.tsv')
cluster = read_cluster_file('data/IRF9_gene_promoter_cluster_tf_binding.table')
df = merge_rnaseq_capstarseq_cluster(cluster,rnaseq_capstarseq)
df['noninduced_epromoter_induced_gene'] = np.where((df['IFN Stimulation'] == 'No Epromoter') & (df['Gene Regulation'] == 'Induced'),1,0)
df['induced_epromoter_induced_gene'] = np.where((df['IFN Stimulation'] == 'Epromoter') & (df['Gene Regulation'] == 'Induced'),1,0)
df['induced_epromoter_noninduced_gene'] = np.where((df['IFN Stimulation'] == 'Epromoter') & (df['Gene Regulation'] == 'Non-induced'),1,0)
df = df[['cluster number','noninduced_epromoter_induced_gene','induced_epromoter_induced_gene','induced_epromoter_noninduced_gene']]

df = df.groupby(['cluster number']).sum().reset_index()
#df1 = df1.sort_values('cluster number')
df = df.sort_values(by=['cluster number'], ascending=True)
#print(df);exit(0)
#df1.to_csv('merge.tsv',sep="\t",header=True,index=False)
df.set_index('cluster number', inplace=True)


df1 = df[['noninduced_epromoter_induced_gene']]
df2 = df[['induced_epromoter_induced_gene']]
df3 = df[['induced_epromoter_noninduced_gene']]

fig, (ax1,ax2,ax3) = plt.subplots(ncols=3)
fig.subplots_adjust(wspace=0.0)
sns.heatmap(df1, cmap="Blues",annot=True,fmt='g', annot_kws={'size':10},ax=ax1, cbar=False, linewidths=.5,robust=True)
ax1.set_xticklabels(ax1.get_xticklabels(),rotation=30)
ax1.set_ylabel('')
sns.heatmap(df2, cmap="Greens",annot=True,fmt='g', annot_kws={'size':10},ax=ax2, yticklabels=[], cbar=False, linewidths=.5,robust=True)
ax2.set_xticklabels(ax2.get_xticklabels(),rotation=30)
ax2.set_ylabel('')
sns.heatmap(df3, cmap="Oranges",annot=True,fmt='g', annot_kws={'size':10},ax=ax3, yticklabels=[], cbar=False, linewidths=.5,robust=True)
ax3.set_xticklabels(ax3.get_xticklabels(),rotation=30)
ax3.set_ylabel('')
plt.show()
exit(0)

'''
df =  pd.DataFrame(np.random.rand(25,1), columns=list("A"))
df2 = pd.DataFrame(np.random.rand(25,1), columns=list("W"))
df3 = pd.DataFrame(np.random.rand(25,1), columns=list("X"))

fig, (ax,ax2,ax3) = plt.subplots(ncols=3)
fig.subplots_adjust(wspace=0.0)
sns.heatmap(df, cmap="Reds",annot=True,fmt='g', ax=ax, cbar=False)
sns.heatmap(df2, cmap="Greens",annot=True,fmt='g', ax=ax2, yticklabels=[], cbar=False)
sns.heatmap(df3, cmap="Oranges",annot=True,fmt='g', ax=ax3, yticklabels=[], cbar=False)
plt.show()
'''
