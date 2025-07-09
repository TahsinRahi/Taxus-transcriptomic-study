library(tidyverse)
library(readxl)
library(showtext)

path='/Users/tahsinrahi/Library/Mobile Documents/com~apple~CloudDocs/Transcriptomic_analysis/Processed_data/Ballon_plot/s_18_plus_and_s72_plus_common.csv'

df<-read_csv(path)

df$Regulation <- factor(df$Regulation, levels = c("UP", "DOWN",'NEUTRAL'))

df1<-df %>% 
  group_by(Pathway_groups,Regulation) %>% 
  summarize(gene_count=n_distinct(Transcript_ID)) %>% 
  ungroup() %>% 
  arrange(desc(gene_count)) %>% 
  na.omit()

df1 %>% 
  ggplot(aes(x = reorder(Pathway_groups, gene_count), y = gene_count)) +  # Descending order
  geom_bar(stat = "identity", position = "dodge",fill = 'mediumvioletred', width = 0.4) +
  geom_text(aes(label = gene_count), hjust = -0.3) + 
  facet_wrap(~Regulation,scales='free')+# Add text above bars
  labs(x = 'Pathway', y = 'Gene count') +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  coord_flip()

ggsave(plot=p,'/Users/tahsinrahi/Library/Mobile Documents/com~apple~CloudDocs/Transcriptomic_analysis/Figures_and_Tables/Figure4/bar_s_min_plus_only_not_in_b_min_plus.pdf',device='pdf',width=8,dpi=600)




