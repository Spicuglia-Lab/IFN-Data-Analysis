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

#download_jaspar2020_tfs()

def sort_bed_file(filename):
    #chr1	173640103	173640353	ANKRD45	0	+	
    df = pd.read_csv(filename,sep="\t",header=None,usecols=[0,1,2,3],names=['chr','start','end','gene'])
    df = df.sort_values(by=['chr','start','end'])
    df.to_csv(filename[:-4]+'.sorted.bed',sep="\t",header=None,index=False)
    print(filename[:-4]+'.sorted.bed')
    
#sort_bed_file('ISGF3_peak_Intergenic_250bp.bed')
#sort_bed_file('ISGF3_peak_Other_Promoter_250bp.bed')
#sort_bed_file('ISGF3_peak_Same_promoter_250bp.bed')

'''
#concatenate all ISRE TFs binding sites
cat *.nofilter.bed >isre.nofilter.bed
'''

def tfbs_merge(concat_filename):
    df = pd.read_csv(concat_filename,sep="\t",header=None)
    df = df.sort_values(by=[0,1,2])
    #no filter
    df.to_csv(concat_filename[:-4]+'.sorted.bed',sep="\t",header=None,index=False)
    print(concat_filename[:-4]+'.sorted.bed')
    #filter by score 400
    df[df[4]>=400].to_csv(concat_filename[:-12]+'400filter.sorted.bed',sep="\t",header=None,index=False)
    print(concat_filename[:-12]+'400filter.sorted.bed')
#tfbs_merge('isre.nofilter.bed')

'''
#NOFILTER: merge all TFs (IRF1-9, STAT1, STAT1::STAT2) bedtools with one nucleotide overlaps
bedtools merge -i isre.nofilter.sorted.bed -c 4 -o collapse >isre.nofilter.merge.bed

#400FILTER: merge all TFs (IRF1-9, STAT1, STAT1::STAT2) bedtools with one nucleotide overlaps
bedtools merge -i isre.400filter.sorted.bed -c 4 -o collapse >isre.400filter.merge.bed

#intersect promoters and isre motif 400 filter
bedtools intersect -a ISGF3_peak_Intergenic_250bp.sorted.bed -b isre.400filter.merge.bed -r -wa -wb >ISGF3_peak_Intergenic_250bp.isre.400filter.merge.intersect.bed
bedtools intersect -a ISGF3_peak_Other_Promoter_250bp.sorted.bed -b isre.400filter.merge.bed -r -wa -wb >ISGF3_peak_Other_Promoter_250bp.isre.400filter.merge.intersect.bed
bedtools intersect -a ISGF3_peak_Same_promoter_250bp.sorted.bed -b isre.400filter.merge.bed -r -wa -wb >ISGF3_peak_Same_promoter_250bp.isre.400filter.merge.intersect.bed

#intersect promoters and isre motif no filter
bedtools intersect -a ISGF3_peak_Intergenic_250bp.sorted.bed -b isre.nofilter.merge.bed -r -wa -wb >ISGF3_peak_Intergenic_250bp.isre.nofilter.merge.intersect.bed
bedtools intersect -a ISGF3_peak_Other_Promoter_250bp.sorted.bed -b isre.nofilter.merge.bed -r -wa -wb >ISGF3_peak_Other_Promoter_250bp.isre.nofilter.merge.intersect.bed
bedtools intersect -a ISGF3_peak_Same_promoter_250bp.sorted.bed -b isre.nofilter.merge.bed -r -wa -wb >ISGF3_peak_Same_promoter_250bp.isre.nofilter.merge.intersect.bed
'''
def count_isre_binding_freq_in_promoter(no_ISGF3_in_promoter,filename):
    df = pd.read_csv(filename,sep="\t",header=None)
    
    df = df.groupby([0, 1, 2,3]).size().reset_index(name="No_ISRE_Binding")
    freq = no_ISGF3_in_promoter-df.shape[0]
    df = df.groupby(['No_ISRE_Binding']).size().reset_index(name="Frequency")
    df = df.append({'No_ISRE_Binding': 0, 'Frequency': freq}, ignore_index = True)
    print(df)

#NO FILTER
#Intergenic 191 (#of Intergenic regions binding with ISGF3)
#count_isre_binding_freq_in_promoter(191,'isre_motif/final_results/ISGF3_peak_Intergenic_250bp.isre.nofilter.merge.intersect.bed')

#Other 134 (#of Other Promoters binding with ISGF3)
#count_isre_binding_freq_in_promoter(134,'isre_motif/final_results/ISGF3_peak_Other_Promoter_250bp.isre.nofilter.merge.intersect.bed')
#Same 68 (#of Same Promoters binding with ISGF3)
#count_isre_binding_freq_in_promoter(68,'isre_motif/final_results/ISGF3_peak_Same_promoter_250bp.isre.nofilter.merge.intersect.bed')

#400 FILTER
#Intergenic 191 (#of Intergenic regions binding with ISGF3)
#count_isre_binding_freq_in_promoter(191,'isre_motif/final_results/ISGF3_peak_Intergenic_250bp.isre.400filter.merge.intersect.bed')
#Other 134 (#of Other Promoters binding with ISGF3)
#count_isre_binding_freq_in_promoter(134,'isre_motif/final_results/ISGF3_peak_Other_Promoter_250bp.isre.400filter.merge.intersect.bed')
#Same 68 (#of Same Promoters binding with ISGF3)
count_isre_binding_freq_in_promoter(68,'isre_motif/final_results/ISGF3_peak_Same_promoter_250bp.isre.400filter.merge.intersect.bed')
