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

def generate_bedtools_merged_file(fliterscore=0):
    tf_bed_files = [ fname.split('.')[0] for fname in os.listdir('.') if fname.endswith('.nofilter.bed')]

    for file in tf_bed_files:
        cmd_filter = "awk '($5 >="+str(fliterscore)+")' " +file+'.nofilter.bed >'+file+'.'+str(fliterscore)+'filter.bed'
        os.system(cmd_filter)
        print(file+'.'+str(fliterscore)+'filter.bed is generated...')

#generate_bedtools_merged_file(0)#enter the filterscore -log10Pval

def generate_bedtools_merged_file(filename):
    print('bedtools merging peaks of '+filename+'...')
    cmd = 'bedtools merge -i '+filename+' >'+filename[:-4]+'.merged.bed'
    os.system(cmd)
    print("\t"+filename[:-4]+'.merged.bed is generated...')
'''
generate_bedtools_merged_file('IRF1.nofilter.bed')
generate_bedtools_merged_file('IRF5.nofilter.bed')
generate_bedtools_merged_file('IRF9.nofilter.bed')
generate_bedtools_merged_file('IRF2.nofilter.bed')
generate_bedtools_merged_file('IRF6.nofilter.bed')
generate_bedtools_merged_file('STAT1.nofilter.bed')
generate_bedtools_merged_file('IRF3.nofilter.bed')
generate_bedtools_merged_file('IRF7.nofilter.bed')
generate_bedtools_merged_file('STAT1::STAT2.nofilter.bed')
generate_bedtools_merged_file('IRF4.nofilter.bed')
generate_bedtools_merged_file('IRF8.nofilter.bed')
'''

def intersection_promoter_TF_binding_regions():
    merge_bed_files = [ fname.split('.')[0] for fname in os.listdir('.') if fname.endswith('.merged.bed')]

    for tf in merge_bed_files:
        if not os.path.exists(tf):
            os.mkdir(tf)
        cmd1 = 'bedtools intersect -a ../../induced_gene.final.bed -b '+tf+'.nofilter.merged.bed -r -wa -wb >'+tf+'/induced_gene_'+tf+'.intersect.bed'
        os.system(cmd1)
        print(tf+'/induced_gene_'+tf+'.intersect.bed is generated...')

        cmd2 = 'bedtools intersect -a ../../induced_gene_epromoter.final.bed -b '+tf+'.nofilter.merged.bed -r -wa -wb >'+tf+'/induced_gene_epromoter_'+tf+'.intersect.bed'
        os.system(cmd2)
        print(tf+'/induced_gene_epromoter_'+tf+'.intersect.bed is generated...')

        cmd3 = 'bedtools intersect -a ../../induced_epromoter.final.bed -b '+tf+'.nofilter.merged.bed -r -wa -wb >'+tf+'/induced_epromoter_'+tf+'.intersect.bed'
        os.system(cmd3)
        print(tf+'/induced_epromoter_'+tf+'.intersect.bed is generated...')

#intersection_promoter_TF_binding_regions()

def matrix():
    dir = [n for n in os.listdir('.') if os.path.isdir(os.path.join('.',n))]
    for item in dir:
        print(item)
        intersect_files = [ fname for fname in os.listdir(item) if fname.endswith('.bed')]
        df = pd.DataFrame()
        for file in intersect_files:
            print(file)
            dfname = file.split('.')[0]
            dfname = pd.read_csv(item+'/'+file,sep="\t",header=None,usecols=[3],names=['seq'])
            total_no_of_promoters_binding_tf = len(dfname['seq'].unique())
            dfname = dfname.groupby(['seq']).size().to_frame(file.split('.intersect.bed')[0]).reset_index()
            promoter_fname = dfname.columns[1].split('_'+item)[0]
            df_promoter_fname = pd.read_csv('../../'+promoter_fname+'.final.bed',sep="\t",header=None)
            dfname.loc[dfname.shape[0]] = ['No binding',df_promoter_fname.shape[0]-dfname.shape[0]]
            df = pd.concat([df,dfname],axis=1)
            print(df.head(4))
        df.to_csv(item+'/'+item+'.matrix',sep="\t",header=True,index=False)

#matrix()


def matrix2():
    dir = [n for n in os.listdir('.') if os.path.isdir(os.path.join('.',n))]
    for item in dir:
        print(item)
        intersect_files = [ fname for fname in os.listdir(item) if fname.endswith('.bed')]
        df = pd.DataFrame()
        for file in intersect_files:
            print(file)
            dfname = file.split('.')[0]#name of new dataframe as per the transcription factor
            dfname = pd.read_csv(item+'/'+file,sep="\t",header=None,usecols=[3],names=['seq'])
            total_no_of_promoters_binding_tf = len(dfname['seq'].unique())

            #dfname = dfname.groupby(['seq']).size().to_frame(file.split('.intersect.bed')[0]).reset_index()
            dfname = dfname.groupby(['seq']).size().to_frame(re.sub( r"(_[A-Z])", r" \1", file).split()[0]).reset_index()
            promoter_fname = dfname.columns[1].split('_'+item)[0]
            df_promoter_fname = pd.read_csv('../../'+promoter_fname+'.final.bed',sep="\t",header=None,usecols=[3],names=['seq'])

            no_binding = df_promoter_fname[~df_promoter_fname.seq.isin(dfname.seq)].dropna(how='all')
            no_binding[re.sub( r"(_[A-Z])", r" \1", file).split()[0]] = 0
            all_promoters = dfname.append(no_binding, ignore_index=True)
            #df1.append(df2, ignore_index = True)
            print(all_promoters.head(2))
            print(all_promoters.shape)
            print(dfname.shape)
            #df1 =
            #dfname.loc[dfname.shape[0]] = ['No binding',df_promoter_fname.shape[0]-dfname.shape[0]]
            df = pd.concat([df,all_promoters],axis=1)
            print(df.head(4))
        df.to_csv(item+'/'+item+'.all.promoter.matrix',sep="\t",header=True,index=False)


#matrix2()
