import pandas as pd

################ EPROMOTER CATEGORY ####################
def read_clust_gene_isgf3(filename):
    df = pd.read_csv(filename,sep="\t",header=None,names=['chr','tss','gene','tss_no','pChr','pStart','pEnd','tf','closest_peak','clust_no'])
    return(df)
    
#clust_gene_isgf3 = read_clust_gene_isgf3('clust_gene_tss_closest_ISGF3.bed')

def read_nonclust_gene_isgf3(filename):
    df = pd.read_csv(filename,sep="\t",header=None,names=['chr','tss','gene','tss_no','pChr','pStart','pEnd','tf','closest_peak'])
    return(df)
    
#nonclust_gene_isgf3 = read_nonclust_gene_isgf3('nonclust_gene_tss_closest_ISGF3.bed')

def read_epromoters(filename):
    epromoters = pd.read_csv(filename,sep="\t",header=None,usecols=[3],names=['gene'])
    epromoters['category'] = filename.split('.')[1].split(' IFNa ')[1].split(' (')[0]
    return(epromoters)

epromoters = pd.concat(
                        [read_epromoters('1. IFNa induced genes and epromoter (44).bed'),
                        read_epromoters('2. IFNa induced genes (498).bed'),
                        read_epromoters('3. IFNa induced epromoter (26).bed')],
                        ignore_index=True
                      )

#cat_clust_gene_isgf3 = pd.merge(clust_gene_isgf3,epromoters,on='gene',how='left')
#cat_clust_gene_isgf3 = cat_clust_gene_isgf3.drop_duplicates('gene', keep='first')
#cat_nonclust_gene_isgf3 = pd.merge(nonclust_gene_isgf3,epromoters,on='gene',how='left')
#cat_nonclust_gene_isgf3 = cat_nonclust_gene_isgf3.drop_duplicates('gene', keep='first')

#cat_nonclust_gene_isgf3.to_csv('cat_nonclust_gene_isgf3.txt',sep="\t",header=True,index=False)
#cat_clust_gene_isgf3.to_csv('cat_clust_gene_isgf3.txt',sep="\t",header=True,index=False)


################ EPROMOTER COLUMN ADDITION ####################
def read_rna_capstar_data(filename):
    #df = pd.read_excel(filename,sheet_name=0,header=0,skiprows=1,usecols=['chr','strand','start','end','gene','CapSTARR-seq'])
    df = pd.read_excel(filename,sheet_name=0,header=0,skiprows=1,usecols=['gene','CapSTARR-seq'])
    df = df[df['CapSTARR-seq'] == 'Epromoter']
    return(df)
    
rna_capstar_data = read_rna_capstar_data('Supplementary Data 1_Summary RNA-seq and CapSTARR-seq.xlsx')

############# Set of IGO outside any cluster
#		Results: Number of genes with ISGF3 in “same promoter”, “intergenic” or “other promoter”
#		File: final_table_nonclust_gene_tss_closest_ISGF3.txt

def read_nonclust_gene_closest_isgf3_peak_category(filename):
    df = pd.read_csv(filename,sep="\t",header=None, names = ['chr','tss1','tss2','gene','epromoter_category','pChr','pStart','pEnd','tf_name','distance','epromoter_location'])
    return(df)
    
nonclust_genes = read_nonclust_gene_closest_isgf3_peak_category('nonclust/nonclust_gene_closest_isgf3_peak_category.txt')

############# Set of IGO (induced gene only) inside a cluster but without Epromoter
#		Results: Number of genes with ISGF3 in “same promoter”, “intergenic” or “other promoter”
#		File: final_table_clust_gene_tss_closest_ISGF3.txt

def read_clust_gene_closest_isgf3_peak_category(filename):
    df = pd.read_csv(filename,sep="\t",header=None, names = ['chr','tss1','gene','clust_no','epromoter_category','pChr','pStart','pEnd','tf_name','distance','epromoter_location'])
    return(df)
clust_genes = read_clust_gene_closest_isgf3_peak_category('clust/clust_gene_closest_isgf3_peak_category.txt')

final_table_nonclust_gene = pd.merge(nonclust_genes,rna_capstar_data,on='gene',how='left')
final_table_nonclust_gene[['chr','tss1','gene','epromoter_category','pStart','pEnd','tf_name','distance','epromoter_location','CapSTARR-seq']].to_csv('final_table_nonclust_gene.txt',sep="\t",header=True,index=False)

final_table_clust_gene = pd.merge(clust_genes,rna_capstar_data,on='gene',how='left')
final_table_clust_gene[['chr','tss1','gene','epromoter_category','clust_no','pStart','pEnd','tf_name','distance','epromoter_location','CapSTARR-seq']].to_csv('final_table_clust_gene.txt',sep="\t",header=True,index=False)


