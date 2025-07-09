library(tidyverse)
library(readxl)


read_percentage<-read_excel('/Users/tahsinrahi/Library/Mobile Documents/com~apple~CloudDocs/Transcriptomic_analysis/Processed_data/Read_statistics/read_statistics.xlsx',
                            sheet = 'align_percent')
# Pie chart
ggplot(read_percentage, aes(x="", y=Percentage, fill=Type)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.position="right") +
  geom_text(aes(label = paste0(Percentage, "%")), position = position_stack(vjust = 0.5))->p

ggsave(plot=p,filename='/Users/tahsinrahi/Library/Mobile Documents/com~apple~CloudDocs/Transcriptomic_analysis/Figures_and_Tables/Figure1/pie_align_percentage.pdf',
       device='pdf',dpi=600)

