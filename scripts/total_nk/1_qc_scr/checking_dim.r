# Load libraries
library(Seurat)
library(SeuratDisk)
library(cluster)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(patchwork)
library(cowplot)
library(RColorBrewer)
library(pheatmap)

# Load the integrated data object with clusters (assumed to use dims 1:10)
integrated_data15 <- readRDS("/home/outputs/totalNK_outputs/2_dge/dim15_20/integrated_data15.rds")

# Print available metadata columns
print(colnames(integrated_data15@meta.data))

# Check the first few rows of metadata
head(integrated_data15@meta.data)

# Check the cumulative variance explained by PCs 1:10
pca_obj <- integrated_data15@reductions$pca
var_explained <- pca_obj@stdev^2 / sum(pca_obj@stdev^2)
cumulative_var <- cumsum(var_explained)
print(cumulative_var[1:20])


# Check a silhouette score to assess cluster separation
Score > 0.5: Good separation
Subsample to 10,000 cells
set.seed(123)
subsample_cells <- sample(colnames(integrated_data15), 10000)
subsampled_data <- subset(integrated_data15, cells = subsample_cells)
# Extract UMAP and clusters
umap_coords <- subsampled_data@reductions$umap@cell.embeddings
clusters <- subsampled_data$seurat_clusters
# Compute and print silhouette score
sil <- silhouette(as.numeric(clusters), dist(umap_coords))
mean_sil <- mean(sil[, 3])
print(paste("Average silhouette score (subsample):", mean_sil))