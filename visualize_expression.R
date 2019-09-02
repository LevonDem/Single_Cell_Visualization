##############################################
## Visualize raw and normalized sc-data
## Author: Levon Demirdjian
## Email: LevonDem@gmail.com
## Last Updated: 09/01/2019 
## Description: This script performs PCA and 
## tSNE on Gene Yeo's (2017) single cell data
##############################################

library(bayesplot)
library(data.table)
library(ggplot2)
library(ggpubr)
library(plotly)
library(Rtsne)
library(scater)

###########################################################
## Read in cell identity information and gene expression ##
###########################################################

## Cell identity
cell_info  <- read.table('cell_info.txt', sep = ' ', stringsAsFactors = FALSE, header = TRUE)
cell_names <- cell_info$Run
cell_type  <- factor(cell_info$sample_type)

## Relabel cell classes
levels(cell_type)[levels(cell_type) == "fibroblast-derived"] <- "iPSC"
levels(cell_type)[levels(cell_type) == "motor"] <- "MN"
levels(cell_type)[levels(cell_type) == "neural"] <- "NPC"
cell_type <- factor(cell_type, levels = c("iPSC", "NPC", "MN"))

## Gene expression
gene_expression <- read.table('GE_matrix.txt', header = TRUE, row.names = 1, comment.char = '~')
colnames(gene_expression) <- gsub(x = colnames(gene_expression), pattern = '_Aligned.sortedByCoord.out.bam', replacement = '')


##################################
## Filter lowly expressed genes ##
##################################

## Filter genes with at least a count of 1 in at least 2 cells
filtered_genes <- apply(gene_expression,  1, function(x){length(x[x >= 1]) >= 2}) 
GE_filtered    <- gene_expression[filtered_genes, ]

##########################
## Create an SCE object ##
##########################

cdata           <- data.frame(Cell = colnames(GE_filtered), Cell_Type = cell_type)
rownames(cdata) <- colnames(GE_filtered)

sce <- SingleCellExperiment(assays = list(counts = as.matrix(GE_filtered)), colData = cdata)
sce

## Normalize sce
normalized_sce <- normalize(sce)
logcounts(normalized_sce)[1:10, 1:3]

#######################################
## PCA of raw gene expression matrix ##
#######################################

pca <- prcomp(x = t(logcounts(normalized_sce)))
df  <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], cell_type = cell_type)
pp  <- ggplot(df, aes(x = PC1, y = PC2, col = cell_type)) + geom_point() + 
  ggtitle("PCA of normalized log gene expression") + labs(col = "Cell Type") +
  scale_color_manual(values = c("black", "goldenrod4", "red4")) +
  theme(plot.title = element_text(hjust = 0.5))
ggplotly(pp)


########################################
## tSNE of raw gene expression matrix ##
########################################

tSNE <- Rtsne(t(logcounts(normalized_sce)))
df   <- data.frame(Dim1 = tSNE[[2]][,1], Dim2 = tSNE[[2]][,2], Cell_Type = cell_type)
pp   <- ggplot(df, aes(x = Dim1, y = Dim2, col = cell_type)) + geom_point() + 
  ggtitle("tSNE of raw gene expression") + labs(col = "Cell Type") +
  scale_color_manual(values = c("black", "goldenrod4", "red4")) +
  theme(plot.title = element_text(hjust = 0.5))
ggplotly(pp)

