# Load libraries
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(patchwork)
library(cowplot)
library(RColorBrewer)
library(pheatmap)

# Load the integrated data object with clusters (assumed to use dims 1:10)
integrated_data10 <- readRDS("/home/outputs/totalNK_outputs/2_dge/dim10/integrated_data10_clusters.rds")

# Validate that UMAP and clusters exist
if (!"umap" %in% names(integrated_data10@reductions)) {
  stop("Error: UMAP reduction not found in integrated_data10.")
}
if (is.null(integrated_data10$seurat_clusters)) {
  stop("Error: Cluster assignments not found in integrated_data10.")
}

# Generate a UMAP plot to verify integration by sample
dim_plot <- DimPlot(integrated_data10, group.by = "sample", pt.size = 0.5, label = FALSE) +
  ggtitle("UMAP - Integration by Sample (Dims 1:10)") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right") +
  scale_color_brewer(palette = "Set2")

# Save the DimPlot as a pdf
ggsave("/home/outputs/totalNK_outputs/2_dge/dim10/dimplot10_integration.pdf", 
       plot = dim_plot, width = 10, height = 8, dpi = 600, units = "in")
cat("UMAP plot saved.\n")

# Switch to RNA assay for gene expression analysis
DefaultAssay(integrated_data10) <- "RNA"

# Ensure layers are joined
integrated_data10 <- JoinLayers(integrated_data10, assay = "RNA")

# Get cluster identities
clusters <- unique(Idents(integrated_data10))

# Initialize a data frame to store results
gene_counts <- data.frame(
  Cluster = character(),
  Gene_Count = integer(),
  stringsAsFactors = FALSE
)

# Loop through each cluster to count genes
for (cluster in clusters) {
  # Subset cells in the current cluster
  cluster_cells <- WhichCells(integrated_data10, idents = cluster)
  cluster_subset <- subset(integrated_data10, cells = cluster_cells)
  
  # Get the RNA count matrix for this cluster
  count_matrix <- GetAssayData(cluster_subset, assay = "RNA", slot = "counts")
  
  # Count genes with non-zero expression (detected in at least one cell)
  gene_count <- sum(rowSums(count_matrix > 0) > 0)
  
  # Add to results
  gene_counts <- rbind(gene_counts, data.frame(Cluster = as.character(cluster), Gene_Count = gene_count))
}

# Print the table
print("Number of genes per UMAP cluster (dims 1:10):")
print(gene_counts)

# Save the table to a CSV file
write.csv(gene_counts, "/home/outputs/totalNK_outputs/2_dge/dim10/gene_counts_per_cluster_dims10.csv", 
          row.names = FALSE)
cat("Gene counts table saved.\n")