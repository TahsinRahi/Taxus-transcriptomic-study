library(tidyverse)
library(readxl)
library(showtext)
library(svglite)
library(extrafont)


# Importing the data
df<-read_csv('/Users/tahsinrahi/Library/Mobile Documents/com~apple~CloudDocs/Transcriptomic_analysis/Processed_data/Volcano_plot/s_v_b_plus_significant.csv')


df$Regulation <- factor(df$Regulation, levels=c('UP', 'DOWN', 'NEUTRAL'))


# Creating volcano plot

df %>% 
  ggplot(aes(log2FoldChange,log_P_value,color=Regulation))+
  geom_point(size=0.5)+
  scale_color_manual(values = c("UP"="red","DOWN"="blue","NEUTRAL"="black"))+
  labs(color="Regulation:\nElicited Small\naggregates", y='-log(P value)',x='log2(Fold change)')+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  )+
  guides(color = guide_legend(override.aes = list(size = 3)))->p



ggsave(plot=p,filename='/Users/tahsinrahi/Library/Mobile Documents/com~apple~CloudDocs/Transcriptomic_analysis/Figures_and_Tables/Figure5/volcano_s_b_plus.pdf',
       device='pdf',dpi=600,width = 8, units='in')













