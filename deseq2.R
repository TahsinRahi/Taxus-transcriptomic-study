library(DESeq2)


# Obtainting the output from featureCounts in a table
counts <- read.table("count.txt", header=TRUE, row.names=1)

# Obtaining only the read counts by selecting columns of samples
countData <- counts[,6:9]


# Creating the data frame sampleInfo to mention samples within control or treated category
sampleInfo <- data.frame(
  row.names=colnames(normalized_count),
  condition = factor(c("control", "control", "treated", "treated"))
)

# Creating a DESeq object from normaliz
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = sampleInfo,
  design = ~ condition
)

# Finally obtaining the differential expression data set
dds <- DESeq(dds)
res <- results(dds)
print(res)



