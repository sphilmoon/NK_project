# ---------------------------------------------- #
# Script: extract_barcodes_by_cluster_dims15.R
# Purpose: Extract barcodes per cluster (dims=15, res=0.25)
# ---------------------------------------------- #

library(Seurat)
library(dplyr)

# ------------------------ #
# Paths & Settings
# ------------------------ #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
rds_output_dir <- file.path(output_dir, "rds")
tsv_output_dir <- file.path(output_dir, "tsv")
rds_file <- file.path(rds_output_dir, "merged_seurat_analysis_20250616.rds")  # ✅ Correct file
barcode_csv <- file.path(tsv_output_dir, "barcodes_by_cluster_dims15_res0.25.csv")

dims <- 15
res <- 0.25
cluster_col <- paste0("integrated_snn_res.d", dims, "_r", res)

# ------------------------ #
# Load Seurat object (SCT slot)
# ------------------------ #
cat("📦 Loading merged Seurat object...\n")
merged_rds <- readRDS(rds_file)
seurat_obj <- merged_rds[["SCT"]]  # ✅ assumes SCT normalization

# ------------------------ #
# Check for clustering metadata
# ------------------------ #
if (!(cluster_col %in% colnames(seurat_obj@meta.data))) {
  cat("⚠ Clustering not found for dims =", dims, "res =", res, "→ running clustering...\n")
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = res)
  seurat_obj[[cluster_col]] <- seurat_obj$seurat_clusters
  cat("✅ Re-clustered and stored to column:", cluster_col, "\n")
} else {
  cat("✅ Found existing clustering column:", cluster_col, "\n")
}

# ------------------------ #
# Extract barcode, cluster, sample
# ------------------------ #
seurat_obj$cluster_id <- as.character(seurat_obj[[cluster_col]][, 1])
barcode_cluster_df <- data.frame(
  barcode = colnames(seurat_obj),
  cluster_id = seurat_obj$cluster_id,
  sample = seurat_obj$animal  # assumes this column exists
)

# ------------------------ #
# Save output
# ------------------------ #
write.csv(barcode_cluster_df, barcode_csv, row.names = FALSE)
cat("✅ Barcodes and cluster IDs saved to", barcode_csv, "\n")