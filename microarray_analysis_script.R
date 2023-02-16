# microarray analysis 
# 9/12/2022
# author: Shirin Jamshidi

setRepositories() 

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("Biobase")


install.packages(c("pheatmap", "ggplot2", "gplots", "reshape2", "plyr"))

library(dplyr)
library(limma)
library(pheatmap)
library(ggplot2)
library(gplots)
library(plyr)
library(GEOquery)
library(reshape2)
library(biobase)

#######################################################################
#######################################################################
#######################################################################

setwd("D:/BulkRNA_Ales/microarray_analysis/02_GSE87806")

series <- "GSE87806"
platform <- "GPL16987" 

#### Load data
gset <- getGEO(series, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir="Data/")

class(gset)
names(gset)
length(gset)

if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]] 
# Or gset <- [[1]]
# class(gset)


gr <- c(rep("HD",10), rep("mpn",6), rep("MPN",8))
# length(gr)

#extraction of expression matrix
ex <- exprs(gset)
#ex <- exprs(gset$GSE87806_series_matrix.txt.gz)
#dim(ex)
#max(ex)
#min(ex)

#### Log2 scale, if required
#ex <- log2(ex + 1)
#exprs(gset) <- ex

##############Quality control##############################
####Box plot to see if the values are normalized
pdf("Results/boxplot.pdf", width=64)
boxplot(ex)
dev.off()

#### Normalize, if required
#ex <- normalizeQuantiles(ex)
#exprs(gset) <- ex

#### Correlation Heatmap
pdf("Results/CorHeatmap.pdf", width=15, height=15)
pheatmap(cor(ex), labels_row=gr, labels_col=gr, color=bluered(256), boarder_color=NA)
dev.off()

#### Principal Component Analysis for genes
pc <- prcomp(ex)
pdf("Results/PC.pdf")
plot(pc)

# this show those genes that have a very low expression or very high expression that is not what we want
#it means it has been sorted the genes based on their expressions
plot(pc$x[,1:2]) 
dev.off()

# we need to see the expression change
# a solution is subtracting each row (genes) from mean of expression (samples are in columns)
# we apply scale function but it works on columns so we transpose the row to columns
# genes should be in the columns to be scaled 
ex.scale <- t(scale(t(ex), scale=F)) 
pc <- prcomp(ex.scale)
pdf("Results/PC_scaled.pdf")
plot(pc)
plot(pc$x[,1:2]) 
dev.off()

names(pc) #x:genes,rotation:samples in pc
dim(pc$x)
colnames(pc$x)

#### Principal Component Analysis for samples
#class(pc$rotation)
pcr <- data.frame(pc$r[,1:3], Group=gr)
#head(pcr)
#pdf("Results/PCA_samples.pdf", width=15, height=15)
pdf("Results/PCA_samples.pdf")
ggplot(pcr, aes(PC1, PC2, color=Group)) + geom_point(size=3) + theme_bw()
dev.off()

##############Quality control is finished######################
##############Finding genes with deferential expression##############################
####Differential Expression Analysis
gr <- factor(gr)
gset$group <- gr
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gr)
#head(design)

# fit linear model
fit <- lmFit(gset, design) 
# set up contrasts of interest and recalculate model coefficients
#cts <- paste(gr, c(tail(gr, -1), head(gr, 1)), sep="-")
#cont.matrix <- makeContrasts(contrasts=cts, levels=design)
cont.matrix <- makeContrasts(HD-MPN, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant 
# B is the statistics and can be replaced by logFC
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
#head(tT)
colnames(tT)
head(ex)
head(tT$ID)

##############Finding certain genes expression##############################
####gene Expression Analysis
#information
# HERC5:
# Chromosome: 8
# Start position: 140,937,913
# End position: 141,068,711
# 
# CPN2:
# Chromosome: 1
# Start position: 140,937,913
# End position: 141,068,711
# 
# PLEKHM1:
# Chromosome: 6
# Start position: 45,369,144
# End position: 45,439,051

# Subset the rows for HERC5
HERC5_rows <- subset(tT, chromosome == 8 & start_position >= 140937913 & end_position <= 141068711)

# Subset the rows for CPN2
CPN2_rows <- subset(tT, chromosome == 1 & start_position >= 140937913 & end_position <= 141068711)

# Subset the rows for PLEKHM1
PLEKHM1_rows <- subset(tT, chromosome == 6 & start_position >= 45369144 & end_position <= 45439051)

# Bind the rows together to create a new data frame
results <- rbind(HERC5_rows, CPN2_rows, PLEKHM1_rows)

# Check if there are any rows in the results
if (nrow(results) == 0) {
  print("There are no rows in the tT dataset that match the specified chromosome and position range for HERC5, CPN2, and PLEKHM1.")
} else {
  # Write the results to a table file
  write.table(results, file = "Results/target_genes_long.txt", row.names = FALSE, sep = "\t")
  print("The relevant rows have been written to the file 'gene_expression_results.txt'.")
  
  write.table(results[, c("ID", "chromosome", "start_position", "end_position", "gene_description", "logFC", "adj.P.Val")],
              file = "Results/target_genes_short.txt", row.names = FALSE, sep = "\t")
}

###############################################################################


tT <- subset(tT, select=c("ID","chromosome","adj.P.Val","logFC"))
#head(tT)
#colnames(tT)
#write.table(tT, file=stdout(), row.names=F, sep="\t")
write.table(tT, "Results/DiffExprs_HD-MPN.txt", row.names=F, sep="\t", quote=F)

##############Deferential expression is finished##############################
#### more analysis toward Gene ontology & Pathway analysis
mpn.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
dim(mpn.up)
mpn.up.genes <- unique(mpn.up$ID)# uniqe(mpn.up$Gene.symbol)
#mpn.up.genes <- unique(as.character(strsplit2(mpn.up$ID,mpn.up$chromosome,"///"))) 
#mpn.up.genes <- sub("///.*","",mpn.up.genes)
write.table(mpn.up.genes, file="Results/MPN_Up.txt", quote=F, row.names=F, col.names=F)

mpn.down <- subset(tT, logFC < -1 & adj.P.Val < 0.05)
dim(mpn.down)
mpn.down.genes <- unique(mpn.down$ID) # uniqe(mpn.down$Gene.symbol)
#mpn.down.genes <- sub("///.*","",mpn.down.genes)
write.table(mpn.down.genes, file="Results/MPN_Down.txt", quote=F, row.names=F, col.names=F)

# Enrichr 
# kegg
# reactome
# wikipath


#x <- c("OCT4///POU5F1","SOX2///SOX2A","NANOG)
#y <- strsplit2(x,"///)
#as.character(y)



