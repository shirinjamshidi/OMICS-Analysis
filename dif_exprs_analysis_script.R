# dif exprs analysis _ RNA-Seq analysis 
# 12/12/2022
# author: Shirin Jamshidi

setRepositories() 

install.packages(c("GenomicRanges","DESeq2", "limma", "pheatmap", "ggplot2", "gplots", "reshape2", "plyr"))
BiocManager::install("GenomeInfoDb", force = TRUE)

library(ggplot2)
library(limma)
library(pheatmap)
library(reshape2)
library(plyr)
library(gplots)
library(DESeq2)

#setwd("C:/Users/k1507907/Downloads/diff_exprs_analysis")
#files <- list.files(".", "*.csv")

#dif <- lapply(files, read.csv)


#dif <- read.csv("C:/Users/k1507907/Downloads/diff_exprs_analysis/g2vs1_DiffRegGenesNoS103.csv")
#dif <- read.csv("C:/Users/k1507907/Downloads/diff_exprs_analysis/g3vs1_DiffRegGenesNoS103.csv")
dif <- read.csv("C:/Users/k1507907/Downloads/diff_exprs_analysis/g3vs2_DiffRegGenesNoS103.csv")

dif$padj <- p.adjust(dif$pvalue, method="BH")
dif <- dif[order(dif$padj),]

#volcano plot LogFC vs -log10 adjusted value
pdf("C:/Users/k1507907/Downloads/diff_exprs_analysis/dif_logFC-logAdjP.pdf", width=10, height=10)
ggplot(dif, aes(log2FoldChange, -log10(padj))) + geom_point() + theme_bw()
dev.off()

# volcano plot with pvalue
pdf("C:/Users/k1507907/Downloads/diff_exprs_analysis/dif_logFC-logPval.pdf", width=10, height=10)
ggplot(dif, aes(log2FoldChange, -log10(pvalue), color=log2FoldChange)) + geom_point() + theme_bw()
dev.off()

dif.up <- subset(dif, log2FoldChange > 1 & padj < 0.05)
dim(dif.up)
dif.up.genes <- unique(dif.up$GeneSymbol)
write.table(dif.up.genes, file="dif_Up.txt", quote=F, row.names=F, col.names=F)

dif.down <- subset(dif, log2FoldChange < -1 & padj < 0.05)
dim(dif.down)
dif.down.genes <- unique(dif.down$GeneSymbol) 
write.table(dif.down.genes, file="dif_Down.txt", quote=F, row.names=F, col.names=F)







