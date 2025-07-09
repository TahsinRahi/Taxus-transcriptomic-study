library(tidyverse)
library(showtext)


path='/Users/tahsinrahi/Library/Mobile Documents/com~apple~CloudDocs/Transcriptomic_analysis/Processed_data/Ballon_plot/s_72_plus_only_not_in_s_18_plus_pathway_enrichment.csv'

df=read_csv(path)

df$Regulation <- factor(df$Regulation, levels = c("UP", "DOWN"))

df<-df %>% 
  group_by(Regulation) %>% 
  arrange(desc(Gene),.by_group = TRUE) %>% 
  ungroup()

# Deleting primary pathway in the upregulation
df<-df %>% 
  filter(!(Regulation=='UP' & Pathway=='Primary')) %>% 
  filter(!(Regulation=='DOWN' & Pathway=='Paclitaxel'))

## Slicing top 10 pathways from both the regulation
df_sliced<-df %>% 
  group_by(Regulation) %>% 
  slice_head(n=10) %>% 
  ungroup()


up<-df %>% 
  filter(Regulation=='UP') %>% 
  filter(Pathway!="Primary") %>% 
  head(10)

down<-df %>% 
  filter(Regulation=='DOWN') %>% 
  head(10)

# Create the balloon plot
ggplot(df_sliced, aes(x = Enrichment_Score, y = reorder(Pathway, Enrichment_Score), 
               size = Gene, color = p_adjust)) +
  geom_point(alpha = 0.8) +  # Add points with transparency
  facet_wrap(~Regulation,scales='free_y')+
  scale_size(range = c(2, 5)) +  # Adjust bubble size
  scale_color_gradientn(colours = c("red", "orange", "yellow", "green")) +  # Color gradient for p-value
  theme_bw() +  # Set font & size
  theme(
    legend.position = "right") +
  scale_x_continuous(expand = expansion(mult = 0.2)) +  # Expand x-axis by 10%
  labs(
    x = "Enrichment Score",
    y = "Pathway",
    size = "Gene Count",
    color = "Adjusted p-value"
  )->p

ggsave(plot=p,"/Users/tahsinrahi/Library/Mobile Documents/com~apple~CloudDocs/Transcriptomic_analysis/Figures_and_Tables/Figure5/ballon_s_72_plus_only_not_in_s_18_plus.pdf", device='pdf',width = 8, units = "in", dpi = 600)

