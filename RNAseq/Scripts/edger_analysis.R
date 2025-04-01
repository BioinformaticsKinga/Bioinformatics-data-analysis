library(edgeR)

counts <- read.table("counts.txt", header=TRUE, row.names=1)
coldata <- data.frame(condition=c("control", "drought"))  

y <- DGEList(counts=counts, group=coldata$condition)
y <- calcNormFactors(y)

design <- model.matrix(~ coldata$condition)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

qlf <- glmQLFTest(fit)

res <- topTags(qlf, n=Inf)

write.csv(as.data.frame(res), file="zygmukin_edger_results.csv")
