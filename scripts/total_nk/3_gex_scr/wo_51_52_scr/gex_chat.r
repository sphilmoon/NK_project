# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(dplyr)

# ------------------------- #
# Define Paths and Configs
# ------------------------- #

rds_file <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/rds/integrated_data_dims25_res0.3_genecounts.rds"

output_dir <- "/home/outputs/totalNK_outputs/3_gex/wo_51_52"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------- #
# Load Seurat Object
# ------------------------- #
seurat_obj <- readRDS(rds_file)
cat("✅ Loaded Seurat object from", rds_file, "\n")

# ------------------------- #
# Set Assay (Optional but Recommended)
# ------------------------- #
# Use SCT if you used SCTransform
# Use RNA if you used log-normalization
if ("SCT" %in% Assays(seurat_obj)) {
  DefaultAssay(seurat_obj) <- "SCT"
} else {
  DefaultAssay(seurat_obj) <- "RNA"
}
cat("✅ Default assay set to", DefaultAssay(seurat_obj), "\n")

# ------------------------- #
# Run Differential Expression
# ------------------------- #
cat("🔍 Running FindAllMarkers() across clusters...\n")

markers_all <- FindAllMarkers(seurat_obj,
                              only.pos = TRUE,          # only keep upregulated markers
                              min.pct = 0.25,           # gene expressed in ≥25% of cells
                              logfc.threshold = 0.25,   # log fold-change cutoff
                              test.use = "wilcox")      # default test (Wilcoxon Rank Sum)

# Save full marker results
marker_csv_all <- file.path(output_dir, "differential_expression_all_clusters.csv")
write.csv(markers_all, marker_csv_all, row.names = FALSE)
cat("✅ Full DEG results saved to", marker_csv_all, "\n")

# ------------------------- #
# Save Top 10 Markers per Cluster
# ------------------------- #
top10 <- markers_all %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

top10_csv <- file.path(output_dir, "top10_markers_per_cluster.csv")
write.csv(top10, top10_csv, row.names = FALSE)
cat("✅ Top 10 markers per cluster saved to", top10_csv, "\n")