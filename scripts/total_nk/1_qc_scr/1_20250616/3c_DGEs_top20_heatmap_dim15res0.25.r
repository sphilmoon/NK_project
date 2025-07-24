# ----------------------------- #
# Libraries
# ----------------------------- #
library(Seurat)
library(dplyr)
library(ggplot2)

# ------------------------ #
# Paths
# ------------------------ #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
rds_file <- file.path(output_dir, "rds", "merged_seurat_analysis_20250616.rds")
marker_file <- file.path(output_dir, "tsv", "merged", "SCT", "new_markers_merged_SCT_dims15_res0.25.csv")
pdf_file <- file.path(output_dir, "pdf", "3_merged_DEGs_outs", "heatmap_top20_cluster_markers.pdf")
top20_csv_file <- file.path(output_dir, "tsv", "merged", "SCT", "top20_markers_dims15_res0.25.csv")

# ------------------------ #
# Load data
# ------------------------ #
merged_obj <- readRDS(rds_file)
markers <- read.csv(marker_file, stringsAsFactors = FALSE)

# ------------------------ #
# Set correct cluster identity
# ------------------------ #
Idents(merged_obj) <- merged_obj$`integrated_snn_res.0.25`  # Use actual metadata column

# ------------------------ #
# Filter top 20 markers per cluster
# ------------------------ #
if (all(c("cluster", "gene") %in% colnames(markers))) {
  top_markers <- markers %>%
    filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
    group_by(cluster) %>%
    slice_max(n = 20, order_by = avg_log2FC)

  # Save to CSV
  write.csv(top_markers, top20_csv_file, row.names = FALSE)
  cat("✅ Top 20 markers exported to:", top20_csv_file, "\n")
  
  # Filter gene names for heatmap
  top_genes <- unique(top_markers$gene)
  top_genes <- intersect(top_genes, rownames(merged_obj))
  if (length(top_genes) == 0) stop("❌ No matching genes found in Seurat object.")

  # ------------------------ #
  # Generate Heatmap
  # ------------------------ #
  heatmap_plot <- DoHeatmap(merged_obj, features = top_genes, group.by = "seurat_clusters") +
    ggtitle("Top 20 Markers per Cluster Heatmap")

  ggsave(
    filename = pdf_file,
    plot = heatmap_plot,
    width = 12,
    height = 16,
    dpi = 600
  )
  cat("✅ Heatmap saved to:", pdf_file, "\n")

} else {
  cat("⚠️ Required columns ('cluster', 'gene') not found in marker table.\n")
}