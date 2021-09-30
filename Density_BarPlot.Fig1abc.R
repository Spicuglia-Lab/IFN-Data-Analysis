library(ggplot2)
library(ggridges)
library(plyr)

################### DENSITY PLOT (FIG 1A & 1B) ###################
df <- read.csv('Muerdter_2017_density_plot_data.tsv',sep="\t",header=TRUE)
df <- df[,c('distance','category')]

#Calculate Pvalue
pval <- ks.test(no_inh$distance, inh$distance)
pval

#density plot (Fig 1a)
ggplot()+
  geom_density(data=df,aes(distance,colour=category))+
  xlim(0,1000000)+
  scale_x_log10(breaks = c(10^1,10^2,10^3,10^4,10^5,10^6))+
  #scale_x_log10()+
  scale_fill_discrete(breaks = rev(levels(df$category)))+
  theme_classic()

#density plot (Fig 1b)
ggplot()+
  geom_density(data = df, aes(distance,colour = category))+
  scale_x_log10()+
  xlim(0,5000)+
  #xlim(0,1000000)+
  theme_classic()+
  theme(legend.key.size = unit(0.4, "cm"))+
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=11, face="bold"))

########## BAR PLOT (FIG 1C) ###########################
df <- read.csv('Muerdter_2017_density_plot_data.tsv',sep="\t",header=TRUE)
df$pos <- ifelse(df$distance >1000, "Distal", "Epromoter (Proximal)")
########## Count Frequency ################
#No Inhibitor               Distal 4402
#No Inhibitor Epromoter (Proximal)  843
#With Inhibitor               Distal 8912
#With Inhibitor Epromoter (Proximal)  701

# Fig 1c

plot_df<-count(df,c('category','pos'))

plot_df_p <- plot_df[plot_df$pos == 'Epromoter (Proximal)',]
plot_df_d <- plot_df[plot_df$pos == 'Distal',]

#Plot distal regions
p_d <- 
  ggplot(plot_df_d,aes(fill=category,y=freq,x=pos))+
  geom_bar(stat='identity',width = 0.4, position = position_dodge(width=0.45))+
  theme_bw()+
  theme(text = element_text(size = 8),
        plot.margin=unit(c(0.001,0.001,0.001,0.001),"cm"),
        #panel.grid.major = element_line(size=0.01), panel.grid.minor = element_line(size = 0.01),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        #legend.position=c(0.82,0.88),
        legend.position = "none",
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"))
ggsave("distal_bar_plot.png",p_d, width = 2.9, height = 2.4)
#Plot proximal regions
p_p<- ggplot(plot_df_p,aes(fill=category,y=freq,x=pos))+
  geom_bar(stat='identity',width = 0.4, position = position_dodge(width=0.45))+
  theme_bw()+
  theme(text = element_text(size = 8),
        plot.margin=unit(c(0.01,0.01,0.01,0.01),"cm"),
        #panel.grid.major = element_line(size=0.01), panel.grid.minor = element_line(size = 0.01),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        #legend.position=c(0.82,0.88),
        legend.position = "none",
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"))
ggsave("proximal_bar_plot.png",p_p, width = 2.9, height = 2.4)
