library(tidyverse)
library(readxl)
library(purrr)
library(showtext)
library(ggrepel)


source('/Users/tahsinrahi/Library/Mobile Documents/com~apple~CloudDocs/Transcriptomic_R_code/custom-R-function.R')

font_add("Times New Roman", "times.ttf")  # Ensure you have this font
showtext_auto()

path=c('/Users/tahsinrahi/Library/Mobile Documents/com~apple~CloudDocs/Transcriptomic_analysis/Processed_data/b_s_18_plus.csv',
       '/Users/tahsinrahi/Library/Mobile Documents/com~apple~CloudDocs/Transcriptomic_analysis/Processed_data/b_s_72_plus.csv')
df1=read_csv(path[1])
df2=read_csv(path[2])

df<-merge_df(list(df1,df2))

pca_scores<-pca(df,c(2,3,4,5))

ggplot(pca_scores, aes(x = PC1, y = PC2, label = rownames(pca_scores))) +
  geom_point(color = "blue", size = 3) +
  geom_text_repel(size = 4, family = "Times New Roman") +  # Smart label placement
  labs(x = "PC1", y = "PC2") +
  theme_bw() +
  theme(text = element_text(family = 'Times New Roman', size = 12))
