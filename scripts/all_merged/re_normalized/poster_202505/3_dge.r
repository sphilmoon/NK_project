# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(dplyr)
library(ggplot2)

# ------------------------- #
# Define Input and Output Paths
# ------------------------- #
output_dir <- "/home/outputs/all_merged_TotalNK_nkp46/re_normalized/20250531"

rds_dir <- file.path(output_dir, "rds")
dge_dir <- file.path(output_dir, "dge")
pdf_dir <- file.path(output_dir, "pdf")

# Ensure directories exist
dir.create(dge_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# Path to the merged RDS file
input_rds_file <- file.path(rds_dir, "merged_totalNK_nkp46_dim25_res0.5.rds")

# ------------------------- #
# Load Merged Seurat Object
# ------------------------- #
cat("📦 Loading merged Seurat object...\n")
merged_obj <- readRDS(input_rds_file)
cat("✅ Loaded:", input_rds_file, "\n")

# ------------------------- #
# Validate Clustering Metadata
# ------------------------- #
if (!"seurat_clusters" %in% colnames(merged_obj@meta.data)) {
    stop("❌ 'seurat_clusters' not found in metadata. Clustering must be completed first.")
}
cat("✅ Cluster labels found:\n")
print(levels(merged_obj$seurat_clusters))


# ------------------------- #
# Run Cluster-Based DGE Analysis
# ------------------------- #
cat("📊 Running FindAllMarkers on merged object...\n")

# Set correct default assay and join layers
DefaultAssay(merged_obj) <- "RNA"

# Join RNA layers if needed (required for Seurat v5 integration)
merged_obj <- JoinLayers(merged_obj)

# Normalize RNA assay again (required before DGE if not already done)
merged_obj <- NormalizeData(merged_obj, verbose = FALSE)

# Use clustering identity
Idents(merged_obj) <- "seurat_clusters"

# Run FindAllMarkers
markers <- FindAllMarkers(
    merged_obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    assay = "RNA"
)

# Save full marker table only if non-empty
dge_path <- file.path(dge_dir, "merged_cluster_markers.tsv")

if (nrow(markers) > 0) {
    write.table(markers, dge_path, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("💾 DGE results saved to:", dge_path, "\n")
} else {
    cat("⚠️ No DGE markers were identified. Skipping save.\n")
}

# ------------------------- #
# Export Top Markers Per Cluster (if any)
# ------------------------- #
if ("cluster" %in% colnames(markers)) {
    top_markers <- markers %>%
        group_by(cluster) %>%
        top_n(n = 3, wt = avg_log2FC)

    cat("🔍 Top markers per cluster:\n")
    print(top_markers)

    # Heatmap
    top_genes <- unique(top_markers$gene)
    heatmap_plot <- DoHeatmap(merged_obj, features = top_genes, group.by = "seurat_clusters") +
        ggtitle("Top Cluster Marker Heatmap")

    ggsave(
        filename = file.path(pdf_dir, "heatmap_top_cluster_markers.pdf"),
        plot = heatmap_plot,
        width = 12, height = 10, dpi = 600
    )
    cat("✅ Heatmap saved to:", file.path(pdf_dir, "heatmap_top_cluster_markers.pdf"), "\n")
} else {
    cat("⚠️ No top markers to display or cluster column missing.\n")
}
