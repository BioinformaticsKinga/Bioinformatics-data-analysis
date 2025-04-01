library(limma)
library(edgeR)

counts <- read.table("counts.txt", header=TRUE, row.names=1)
coldata <- data.frame(condition=c("control", "drought"))  

d <- DGEList(counts=counts, group=coldata$condition)
d <- calcNormFactors(d)
v <- voom(d, design=model.matrix(~ coldata$condition), plot=TRUE)

fit <- lmFit(v, design=model.matrix(~ coldata$condition))
fit <- eBayes(fit)

res <- topTable(fit, number=Inf)

write.csv(res, file="zygmukin_limma_voom_results.csv")
