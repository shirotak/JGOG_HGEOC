# DESeq from Read_counts
# https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
setwd("/Users/tshiro/Desktop/Projects/HGEOC/Github/JGOG_HGEOC/code/JGOG_C4_signature")
library(DESeq2)
library(dbplyr)
library(data.table)
library(limma)
library(edgeR)

cts=fread('/Users/tshiro/Desktop/Projects/JGOG3025_HGS/JGOG3025_data/RNAseq/gene_expression_add/JGOG_HGS_RSEM_readcounts_renamed.txt.gz')
countData <- as.matrix(cts[,-1])
rownames(countData)=cts$ens_gene_id
head(countData)

phenoData <- read.delim("../../data/JGOG_282_assigned_cluster.txt")
head(phenoData)

countData=countData[,c(phenoData$X)]

# Create a DESeqDataSet from count matrix and labels
colData = phenoData[,c('X','C4')]
rownames(colData)=colData$X
all(rownames(colData) %in% colnames(countData))
all(rownames(colData) == colnames(countData))
## Phenotype needed to be factor
colData = apply(X = colData, MARGIN = 2, FUN = as.factor)

ddsFull <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, 
                              design = ~ C4
                              )
# Run 
dds <- DESeq(ddsFull)
resultsNames(dds) 
# Contrast
cons=c("C4", "1", "0")
res1 <- results( dds, contrast = cons )
res1Ordered <- res1[order(res1$padj), ]
res1Ordered
# write out
write.table(res1Ordered, file= "./JGOG_282_C4_DESeq2.txt",
            row.names = T, col.names = NA,
            sep = "\t",quote = F)

# limma_voom
d = DGEList(countData)
colData
group = as.factor(colData[,2])
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = F)
fit <- lmFit(y, mm)
head(mm)
contr <- makeContrasts( group1 - group0, levels = colnames(coef(fit)) )
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)

write.table(top.table, file= "./JGOG_282_C4_Limma.txt",
            row.names = T, col.names = NA,
            sep = "\t",quote = F)
