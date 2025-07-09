library(tidyverse)
library(factoextra)

rm(list=ls())
df<-read_csv('/Users/tahsinrahi/Library/Mobile Documents/com~apple~CloudDocs/Transcriptomic_analysis/Processed_data/Ballon_plot/s_b_72_plus_pathway_assigned.csv')

df<-df %>% 
  filter(!is.na(Regulation)) %>% 
  filter(!is.na(B72_plus_Expression) & !is.na(S72_Plus_Expression))



df_filterd<-df %>% 
  select(6,7) %>% 
  mutate(across(everything(), ~replace_na(.x, 0)))


df_scaled <- scale(df_filterd)

# Evaluating optimum number of cluster
mat<-as.matrix(df_scaled)
wss <- sapply(1:12, function(k) {
kmeans(mat, centers = k, nstart = 25)$tot.withinss
})

wss_df <- data.frame(k = 1:12, WSS = wss)

# Plot the Elbow Method
ggplot(wss_df, aes(x = k, y = WSS)) +
  geom_point(size = 4, color = "blue") +
  geom_line(color = "blue") +
  scale_x_continuous(breaks = 1:8) +
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  xlab("Number of Clusters (k)") +
  ylab("Total Within-Cluster Sum of Squares (WSS)")
#####################################################################################


# Running the k-means clustering
set.seed(123)  # For reproducibility
kmeans_result <- kmeans(df_scaled, centers = 7, nstart = 25)

df$cluster <- as.factor(kmeans_result$cluster)

# Compute PCA for visualization
pca_result <- prcomp(df_scaled, center = TRUE, scale. = TRUE)
pca_data <- data.frame(pca_result$x[, 1:2], cluster = df$cluster)

# Plot PCA results with clusters
ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))->p

ggsave(plot=p,'/Users/tahsinrahi/Library/Mobile Documents/com~apple~CloudDocs/Transcriptomic_analysis/Figures_and_Tables/Figure6 (k-means)/kmeans_72_hour.pdf',width=8,device='pdf',dpi=600,)

write_csv(df,'/Users/tahsinrahi/Library/Mobile Documents/com~apple~CloudDocs/Transcriptomic_analysis/Processed_data/k-means/kmeans_72_hour.csv')




