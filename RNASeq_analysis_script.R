# RNA-Seq analysis 
# 11/12/2022
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

options(stringsAsFactors=F)
setwd("~/RNA-Seq")
files <- list.files(".", "*.count")
cn <- lapply(files, read.delim, header=F, comment.char = "_")
cn <- do.call(cbind, cn)
rownames(cn) <- cn[,1]
cn <- cn[,-seq(1, ncol(cn),2)]
colnames(cn) <- sub(".count", "", files)

#head(cn)
#tail(cn)
#colSums(cn) #reads in each column

gr <- factor(c(rep("YRI",3),rep("GBR", 2),"YRI",rep("GBR",4),rep("YRI",2)))
colData <- data.frame(group=gr, type="paired-end")
#making expression matrix
cds <- DESeqDataSetFromMatrix(cn, colData, design= ~group)
#normalize
cds <- DESeq(cds)
#log2
cnt <- log2(1+counts(cds, normalized=T))

head(cnt)
boxplot(cnt)
#boxplot(1+cn)

dif <- data.frame(results(cds, c("group", "YRI", "GBR")))
dif$padj <- p.adjust(dif$pvalue, method="BH")
dif <- dif[order(dif$padj),]
#head(dif)
#class(dif)


#volcano plot LogFC vs -log10 adjusted value
ggplot(dif, aes(log2FoldChange, -log10(padj))) + geom_point() + theme_bw()

# volcano plot with pvalue
ggplot(dif, aes(log2FoldChange, -log10(pvalue), color=log2FoldChange)) + geom_point() + theme_bw()

