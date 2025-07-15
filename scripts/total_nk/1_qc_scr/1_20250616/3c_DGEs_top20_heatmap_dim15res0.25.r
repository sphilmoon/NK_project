# ----------------------------- #
# Libraries
# ----------------------------- #
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)

# ----------------------------- #
# Paths & Settings
# ----------------------------- #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
csv_path <- file.path(output_dir, "tsv", "merged", "SCT", "new_markers_merged_SCT_dims15_res0.25.csv")
rds_file <- file.path(output_dir, "rds", "merged_seurat_analysis_20250616.rds")
heatmap_pdf <- file.path(output_dir, "pdf", "3_merged_DEGs_outs", "top20_DEGs_heatmap_dim15res0.25.pdf")

dims <- 15
res <- 0.25

# ----------------------------- #
# Load Data
# ----------------------------- #
deg_df <- read.csv(csv_path, stringsAsFactors = FALSE)
seurat_obj <- readRDS(rds_file)[["SCT"]]  # Use SCT assay

# ----------------------------- #
# Select Top 20 Genes per Cluster
# ----------------------------- #
top_degs <- deg_df %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) %>%
  ungroup()

top_genes <- unique(top_degs$gene)

# ----------------------------- #
# Extract Expression Matrix
# ----------------------------- #
expr_mat <- GetAssayData(seurat_obj, slot = "data")  # log-normalized expression
expr_top <- expr_mat[top_genes, ]

# ----------------------------- #
# Optional: Scale across cells
# ----------------------------- #
expr_scaled <- t(scale(t(as.matrix(expr_top))))

# ----------------------------- #
# Generate Annotation for Columns
# ----------------------------- #
cell_annot <- seurat_obj@meta.data %>%
  dplyr::select(animal, seurat_clusters = paste0("integrated_snn_res.", res))

# ----------------------------- #
# Plot and Save Heatmap
# ----------------------------- #
pheatmap(expr_scaled,
         show_colnames = FALSE,
         show_rownames = TRUE,
         fontsize_row = 7,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = paste("Top 20 DEGs per Cluster (dims=", dims, ", res=", res, ")", sep = ""),
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         filename = heatmap_pdf,
         width = 10,
         height = 10)

cat("✅ DEG Heatmap (dims 0.5, res 0.25) saved to:", heatmap_pdf, "\n")