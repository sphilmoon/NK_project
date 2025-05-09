# ---------------------------------------------- #
# Script: extract_barcodes_by_cluster_dims25.R
# Purpose: Extract barcodes per cluster (dims=25, res=0.3) + UMAP plot
# ---------------------------------------------- #

library(Seurat)
library(dplyr)
library(ggplot2)

# Define file paths
rds_file <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/rds/integrated_data_dims25_res0.3_genecounts.rds"
output_csv <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/csv/barcodes_by_cluster_dims25_res0.3.csv"
output_umap <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/pdf/umap_clusters_dims25_res0.3_check.pdf"

# Create output directories if they don't exist
dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_umap), recursive = TRUE, showWarnings = FALSE)

# Load the Seurat object
seurat_obj <- readRDS(rds_file)
cat("✅ Loaded Seurat object from", rds_file, "\n")

# Perform clustering (dims = 1:25, resolution = 0.3) to ensure it's up to date
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:25)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)
cat("✅ Clustering complete using dims 1:25 and resolution 0.3\n")

# Extract barcodes, cluster IDs, and sample (animal ID)
barcode_cluster_df <- data.frame(
  barcode = colnames(seurat_obj),
  cluster_id = Idents(seurat_obj),
  sample = seurat_obj$sample
)

# Save to CSV
write.csv(barcode_cluster_df, output_csv, row.names = FALSE)
cat("✅ Barcodes and cluster IDs saved to", output_csv, "\n")

# # Generate and save a UMAP plot with cluster labels
# if (!"umap" %in% names(seurat_obj@reductions)) {
#   seurat_obj <- RunUMAP(seurat_obj, dims = 1:25)
#   cat("ℹ️ UMAP computed using dims 1:25\n")
# }

# umap_plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.3) +
#   ggtitle("UMAP Clustering (dims=25, res=0.3)") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_brewer(palette = "Set3")

# ggsave(output_umap, plot = umap_plot, width = 8, height = 6, dpi = 600)
# cat("✅ UMAP plot saved to", output_umap, "\n")
