# Single_Cell_Visualization
Visualize scRNA-seq data

---

This repository provides some handy scripts for visualizing scRNA-seq data. The script 'visualize_expression.R' creates interactive scatter plots of gene expression values after performing PCA and tSNE. The inputs to this script are:

1. cell_info.txt: A table with {number of cells} rows which contains the cell identities for each cell.
2. GE_matrix.txt: A table with {number of genes/events/isoforms} rows and {number of cells} columns of gene expression values.

These scripts are meant to be used interactively (don't source them), though they can easily be converted to be sourced if desired.
