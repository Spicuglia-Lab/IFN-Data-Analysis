import pandas as pd
import os,re

isre_tf = {'MA0050.2':'IRF1',
           'MA0051.1':'IRF2',
           'MA1418.1':'IRF3',
           'MA1419.1':'IRF4',
           'MA1420.1':'IRF5',
           'MA1509.1':'IRF6',
           'MA0772.1':'IRF7',
           'MA0652.1':'IRF8',
           'MA0653.1':'IRF9',
           'MA0137.3':'STAT1',
           'MA0517.1':'STAT1::STAT2'
	      }

def download_jaspar2020_tfs():
    '''Download all TFs (hg19) from JASPAR 2020 database '''
    for item in isre_tf.keys():
        print('Reading '+item+'.tsv.gz ...')
        df = pd.read_csv('http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2020/hg19/'+item+'.tsv.gz',sep="\t",header=None,usecols=[0,1,2,3,5],names=['chr','start','end','TF','-log10Pval'])
        df = df.sort_values(['chr','start','end'], ascending = [True,True,True])
        df.to_csv([tf for id,tf in isre_tf.items() if id == item][0]+'.nofilter.bed',
                    sep="\t",header=None,index=False)

        print([tf for id,tf in isre_tf.items() if id == item][0] + ' is extracted')

download_jaspar2020_tfs()

#concatenate all ISRE TFs binding sites
cat *.nofilter.bed >isre.nofilter.bed

#intersect induced gene and induced epromoter
bedtools intersect -a induced_gene.bed -b isre.nofilter.bed -r -wa -wb >induced_gene_isre.intersect.bed
bedtools intersect -a induced_epromoter.bed -b isre.nofilter.bed -r -wa -wb >induced_epromoter_isre.intersect.bed
bedtools intersect -a induced_gene_epromoter.bed -b isre.nofilter.bed -r -wa -wb >induced_gene_epromoter_isre.intersect.bed
