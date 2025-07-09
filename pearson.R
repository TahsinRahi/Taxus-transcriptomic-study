library(tidyverse)
library(showtext)
library(matrixStats)
library(corrplot)


path='/Users/tahsinrahi/Library/Mobile Documents/com~apple~CloudDocs/Transcriptomic_analysis/Processed_data/Ballon_plot/s_b_plus_pathway_assigned.csv'

df = read_csv(path)

geo_mean <- function(x) {
  x <- x[!is.na(x)]  # Remove NA values
  if (length(x) == 0) return(NA)  # Return NA if no valid values
  exp(mean(log(x)))  # Compute geometric mean
}

pc <- df %>% 
  group_by(Pathway_groups) %>% 
  summarise(B18 = geo_mean(as.numeric(B18_plus_Expression)),
            S18=geo_mean(as.numeric(S18_plus_Expression)),
            B72=geo_mean(as.numeric(B72_plus_Expression)),
            S72=geo_mean(as.numeric(S72_Plus_Expression))) %>% 
  ungroup()

#row names for matrix
rn<-pc$Pathway_groups

# Obtaining the matrix for pearson
mat<-pc %>% 
  select(-1) %>% 
  as.matrix()

# Setting up the row names for the matrix 
rownames(mat)<-rn


cor_matrix <- cor(t(mat), method = "pearson")

cor_df <- data.frame(RowName = rownames(cor_matrix),cor_matrix , row.names = NULL)


cor_df_filtered<-cor_df %>% 
  filter(Paclitaxel > 0.91 )



write_csv(cor_df,'/Users/tahsinrahi/Library/Mobile Documents/com~apple~CloudDocs/Transcriptomic_analysis/Processed_data/Pearson/peason_matrix.csv')


corrplot(cor_matrix, method = "color", type = "upper", tl.col = "black", tl.cex = 0.8)














