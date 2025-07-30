# ---------------------------------------------- #
# Script: extract_barcodes_by_cluster_nkp46.R
# Purpose: Extract barcodes per cluster for NKp46+ and NKp46- (dims=15, res=0.25)
# ---------------------------------------------- #

library(Seurat)
library(dplyr)

# ------------------------ #
# Settings
# ------------------------ #
dims <- 15
res <- 0.25
cluster_col <- paste0("integrated_snn_res.d", dims, "_r", res)

# ------------------------ #
# Directories & Files
# ------------------------ #
output_dir <- "/home/outputs/nkp46_outs_20250730_slurm"
rds_neg_file <- file.path(output_dir, "neg", "rds", "nkp46_neg_merged.rds")
rds_pos_file <- file.path(output_dir, "pos", "rds", "nkp46_pos_merged.rds")

tsv_dir <- file.path(output_dir, "tsv", "barcode_by_cluster")
dir.create(tsv_dir, recursive = TRUE, showWarnings = FALSE)

neg_barcode_csv <- file.path(tsv_dir, "barcodes_neg_dims15_res0.25.csv")
pos_barcode_csv <- file.path(tsv_dir, "barcodes_pos_dims15_res0.25.csv")

# ------------------------ #
# Function to extract barcodes
# ------------------------ #
extract_barcodes <- function(rds_file, barcode_csv, condition_label) {
  cat("📦 Loading merged Seurat object for", condition_label, "...\n")
  seurat_obj <- readRDS(rds_file)

  # Run clustering if needed
  if (!(cluster_col %in% colnames(seurat_obj@meta.data))) {
    cat("⚠ Clustering not found → running clustering for", condition_label, "\n")
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dims)
    seurat_obj <- FindClusters(seurat_obj, resolution = res)
    seurat_obj[[cluster_col]] <- seurat_obj$seurat_clusters
  } else {
    cat("✅ Found clustering column:", cluster_col, "\n")
  }

  # Extract barcode info
  seurat_obj$cluster_id <- as.character(seurat_obj[[cluster_col]][, 1])
  barcode_df <- data.frame(
    barcode = colnames(seurat_obj),
    cluster_id = seurat_obj$cluster_id,
    sample = seurat_obj$sample,
    condition = condition_label
  )

  write.csv(barcode_df, barcode_csv, row.names = FALSE)
  cat("✅ Barcodes and cluster IDs saved to:", barcode_csv, "\n")
}

# ------------------------ #
# Run for both conditions
# ------------------------ #
extract_barcodes(rds_neg_file, neg_barcode_csv, "neg")
extract_barcodes(rds_pos_file, pos_barcode_csv, "pos")