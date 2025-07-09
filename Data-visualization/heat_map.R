library(tidyverse)
library(readxl)
library(purrr)
library(pheatmap)


source('/Users/tahsinrahi/Library/Mobile Documents/com~apple~CloudDocs/Transcriptomic_R_code/custom-R-function.R')


df<-read_excel('/Users/tahsinrahi/Library/Mobile Documents/com~apple~CloudDocs/Transcriptomic_analysis/Processed_data/Heat_map/heatmap.xlsx',
               sheet = 'filtered')

df_filtered<-df %>% 
  select(2,3,4,5)

rn<-df$Gene_names

df_filtered <- df_filtered %>%
  mutate(across(everything(), ~ replace(., is.na(.), 0)))

# Convert to matrix
mat <- as.matrix(df_filtered)

rownames(mat)<-rn
custom_colors <- colorRampPalette(c("white", "yellow", "red"))(100)


pheatmap(t(mat), cluster_rows = FALSE, cluster_cols = FALSE, color=custom_colors,
         angle_col = 45,    # Rotate x-axis labels by 45 degrees
         angle_row = 0)
