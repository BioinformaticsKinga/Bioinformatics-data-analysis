library(DESeq2)

counts <- read.table("counts.txt", header=TRUE, row.names=1)
coldata <- data.frame(condition=c("control", "drought"))  

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)

dds <- DESeq(dds)

res <- results(dds)

write.csv(as.data.frame(res), file="zygmukin_deseq2_results.csv")
