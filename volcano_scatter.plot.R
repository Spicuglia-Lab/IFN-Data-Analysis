library(ggplot2)
library(ggpubr)
library(gridExtra)
library(ggrepel)
library(dplyr)

##################### WRITING PLOT DATA ##################################################################
write_data <- function(dataframe,outputfile){
  write.table(dataframe,outputfile,sep="\t",row.names = F)
}

##################### SCATTER PLOT IFN VS NS FUNCTIONS ###################################################
data = read.csv('scatter_plot_starseq_rnaseq_plot_data.tsv',sep="\t",header=TRUE)
#data = read.csv('file.tsv',sep="\t",header=TRUE)
data = na.omit(data)
data = data[data$condition != 'Others',]
p <- ggplot(data,aes(IFN_vs_NS_fold_change.log2FC,DESeq2_log2FoldChange,col=condition))+
    geom_point(pch=20,size=0.5)+
    scale_colour_manual(values = c('orange','#00C800','#009BE6','grey60','darkcyan'))+
    theme_classic()+
    theme(axis.text=element_text(size=6),
          axis.ticks = element_blank(),
          #plot.title = element_text(hjust=0.5,size = 50, colour = "black"),
          legend.title = element_text(size=6),
          plot.margin=unit(c(0.09,2.59,0.09,0.09),"cm"),
          #panel.grid.major = element_line(size=0.1), panel.grid.minor = element_line(size = 0.1),
          #panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.text=element_text(margin = margin(r=10,unit='pt'),size=5),
          legend.position=c(1.35,0.6),
          legend.key.size = unit(8, "pt"),
          legend.spacing.x = unit(1, 'pt'))
ggsave('scatter_starseq_rnaseq.png',p, width = 3.2, height = 2.1)

##################### VOLCANO PLOT FUNCTIONS #############################################################
read_file_volcano_plot <-function(filename){
  genes <- read.table(filename, header = TRUE)[,c('gene','DESeq2_log2FoldChange','DESeq2_pvalue','DESeq2_padj')]
  colnames(genes) <- c('Gene','log2FoldChange','pValue','pAdj')
  genes$Significant <- ifelse(((genes$pAdj < 0.05) & (genes$log2FoldChange>1)) | ((genes$pAdj < 0.05) & (genes$log2FoldChange < -1)), "pAdj < 0.001", "Not Sig")
  genes <- na.omit(genes)
  
  return(genes)
}
volcano_plot<-function(filename){
  genes <- read_file_volcano_plot(filename)
  genes$labels <- ifelse(genes$Gene == 'OAS3' | genes$Gene == 'MX1',TRUE,FALSE)
  genes = genes[order(rev(genes$Significant)),]
  write_data(genes,'plot_data/volcano_plot_data.tsv')
  p<-ggplot(genes,aes(x=log2FoldChange,y=-log10(pAdj)))+
    geom_point(aes(color = Significant),size=0.001)+
    guides(fill = guide_legend(reverse = F))+
    geom_text_repel(aes(x=log2FoldChange,y=-log10(pAdj)),
                    label = ifelse(genes$labels == TRUE, as.character(genes$Genes),""), 
                    box.padding = unit(.7, "lines"),hjust= 0.30)+
    theme(legend.title=element_blank(),text = element_text(size= 13))+
    scale_color_manual(values = c('grey','darkmagenta'),name='Significance',labels = c('Not Sig','log2FC >  1 (pAdj<0.001)\nlog2FC < -1 (pAdj<0.001)'))+
    
    #scale_y_log10(breaks = trans_breaks('log10',function(x) 10^x),
    #              labels = trans_format('log10', math_format(10^.x)))+
    theme_classic()+
    theme(axis.text=element_text(size=8),
          axis.ticks = element_blank(),
          #plot.title = element_text(hjust=0.5,size = 50, colour = "black"),
          legend.title = element_text(size=6),
          plot.margin=unit(c(0.09,3.09,0.09,0.09),"cm"),
          #panel.grid.major = element_line(size=0.1), panel.grid.minor = element_line(size = 0.1),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.text=element_text(margin = margin(r=10,unit='pt'),size=5),
          legend.position=c(1.45,0.6),
          legend.key.size = unit(6, "pt"),
          legend.spacing.x = unit(1, 'pt'))
  #aspect.ratio = 4/3,
  #plot.title = element_blank())
  ggsave('final_volcano_plot.png',p, width = 3.2, height = 2.1)
  return(p)
}
volcano = volcano_plot('RNAseq_K562_DESeq2_EdgeR_pval_and_norm_count_log2.txt')



genes <- read_file_volcano_plot('data/RNAseq_K562_DESeq2_EdgeR_pval_and_norm_count_log2.txt')
genes$labels <- ifelse(genes$Gene == 'IFIT3' | genes$Gene == 'OAS1' | genes$Gene == 'OAS2' | genes$Gene == 'OAS3' | genes$Gene == 'ISG15' | genes$Gene == 'MX1' | genes$Gene == 'MX2',
                       TRUE,'')
genes = genes[order(rev(genes$Significant)),]
#write_data(genes,'plot_data/volcano_plot_data.tsv')
p<-ggplot(genes,aes(x=log2FoldChange,y=-log10(pAdj)),color = Significant,size=0.001)+
  #geom_point(aes(color = Significant),size=0.000001)+
  geom_point(aes(color = Significant),size=0.001)+
  geom_text_repel(aes(x=log2FoldChange,y=-log10(pAdj)),
                  label = ifelse(genes$labels == TRUE, as.character(genes$Gene),""))
ggsave('volcano_plot.png',p, width = 40, height = 10.1,limitsize = FALSE)