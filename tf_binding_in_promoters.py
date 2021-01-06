import pandas as pd
import os
import subprocess as sp

print(
	'''
	JASPER_CODE	TF
	MA0050.2	IRF1
	MA0051.1	IRF2
	MA1418.1	IRF3
	MA1419.1	IRF4
	MA1420.1	IRF5
	MA1509.1	IRF6
	MA0772.1	IRF7
	MA0652.1	IRF8
	MA0653.1	IRF9
	MA0137.3	STAT1
	MA0517.1	STAT1::STAT2
	'''
)

def download_jaspar2020_tfs(jaspar_tf_list):
    '''Merge bed file for all the listed TFs'''
    all_tf =pd.DataFrame()
    for item in jaspar_tf_list:
        print('Reading '+item+'.tsv.gz ...')
        df = pd.read_csv('http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2020/hg19/'+item+'.tsv.gz',
                            sep="\t",header=None,usecols=[0,1,2,3,5],names=['chr','start','end','TF','-log10Pval'])
        #df.drop_duplicates(['chr','start','end'],keep= 'first', inplace = True)
        all_tf = pd.concat([all_tf, df], ignore_index=True, sort=False)
        all_tf = all_tf.sort_values(['chr', 'start','end'], ascending=[True, True,True])

    all_tf[all_tf['-log10Pval']>600].to_csv('STAT1-2_IRF1-9.filter600.bed',sep="\t",header=None,index=False)
    all_tf.to_csv('STAT1-2_IRF1-9.nofilter.bed',sep="\t",header=None,index=False)
    print('Combined bed file for the TFs: STAT1-2, IRF1-9:\n\tSTAT1-2_IRF1-9.nofilter.bed\n\tSTAT1-2_IRF1-9.filter600.bed')

#download_jaspar2020_tfs(['MA1418.1','MA0051.1','MA0050.2'])
#download_jaspar2020_tfs(['MA0050.2','MA0051.1','MA1418.1','MA1419.1','MA1420.1','MA1509.1','MA0772.1','MA0652.1','MA0653.1','MA0137.3','MA0517.1'])

def generate_bedtools_merged_file(filename):
    print('bedtools merging peaks of '+filename+'...')
    cmd = 'bedtools merge -i '+filename+' >'+filename[:-4]+'.merged.bed'
    os.system(cmd)
    print("\t"+filename[:-4]+'.merged.bed is generated...')

#generate_bedtools_merged_file('STAT1-2_IRF1-9.nofilter.bed')
#generate_bedtools_merged_file('STAT1-2_IRF1-9.filter600.bed')

def promoter_final_4line_bedfile():
    cmd1 = """awk '{print $0"\tinduced_epromoter"NR}' induced_epromoter.bed >induced_epromoter.final.bed"""
    os.system(cmd1)
    print('induced_epromoter.final.bed is generated...')

    cmd2 = """awk '{print $0"\tinduced_gene"NR}' induced_gene.bed > induced_gene.final.bed"""
    os.system(cmd2)
    print('induced_gene.final.bed is generated...')

    cmd3 = """awk '{print $0"\tinduced_gene_epromoter"NR}' induced_gene_epromoter.bed >induced_gene_epromoter.final.bed"""
    os.system(cmd3)
    print('induced_gene_epromoter.final.bed is generated...')

    cmd4 = """awk '{print $0"\tISRE"}' STAT1-2_IRF1-9.filter600.merged.bed >STAT1-2_IRF1-9.filter600.merged.final.bed"""
    os.system(cmd4)
    print('STAT1-2_IRF1-9.filter600.merged.final.bed is generated...')

#promoter_final_4line_bedfile()

def intersection_promoter_TF_binding_regions(promoter_category,bedtools_merged_tf_bed):
    #Intersection of promoter regions with all TFs
    cmd4 = """bedtools intersect -a """+promoter_category+""" -b """+bedtools_merged_tf_bed+""" -r -wa -wb >"""+promoter_category[:-10]+"""_"""+bedtools_merged_tf_bed[:-17]+'.intersect.bed'
    os.system(cmd4)
    print(promoter_category[:-10]+"""_"""+bedtools_merged_tf_bed[:-11]+'.intersect.bed is generated...')

#intersection_promoter_TF_binding_regions('induced_epromoter.final.bed','STAT1-2_IRF1-9.filter600.merged.final.bed')
#intersection_promoter_TF_binding_regions('induced_gene.final.bed','STAT1-2_IRF1-9.filter600.merged.bed')
#intersection_promoter_TF_binding_regions('induced_gene_epromoter.final.bed','STAT1-2_IRF1-9.filter600.merged.bed')

#intersection_promoter_TF_binding_regions('induced_gene_epromoter.final.bed','STAT1-2_IRF1-9.nofilter.merged.bed')
#intersection_promoter_TF_binding_regions('induced_gene_epromoter.final.bed','STAT1-2_IRF1-9.nofilter.merged.bed')
#intersection_promoter_TF_binding_regions('induced_gene_epromoter.final.bed','STAT1-2_IRF1-9.nofilter.merged.bed')

def matrix(file):
    induced_epromoter = pd.read_csv(file.split('_STAT1-2_IRF1-9')[0]+'.final.bed',header=None,sep="\t")
    intersect = pd.read_csv(file,sep="\t",header=None,usecols=[3],names=['seq'])
    total_no_of_promoters_binding_tf = len(intersect['seq'].unique())
    intersect = intersect.groupby(['seq']).size().to_frame('no of TF binding sites in '+file.split('_STAT1-2_IRF1-9')[0]).reset_index()
    #print(intersect);exit(0)
    print('no of TF binding sites in '+file.split('_STAT1-2_IRF1-9')[0])
    matrix = intersect[['no of TF binding sites in '+file.split('_STAT1-2_IRF1-9')[0]]].copy()
    matrix = matrix.groupby(['no of TF binding sites in '+file.split('_STAT1-2_IRF1-9')[0]]).size().to_frame('frequency').reset_index()
    matrix.loc[matrix.shape[0]] = [0,induced_epromoter.shape[0]-total_no_of_promoters_binding_tf]
    print(matrix)
    matrix.to_csv(file[:-4]+'.matrix',sep="\t",header=True,index=False)
    print(file[:-4]+'.matrix is generated...')

matrix('induced_epromoter_STAT1-2_IRF1-9.filter400.intersect.bed')
matrix('induced_gene_STAT1-2_IRF1-9.filter400.intersect.bed')
matrix('induced_gene_epromoter_STAT1-2_IRF1-9.filter400.intersect.bed')
