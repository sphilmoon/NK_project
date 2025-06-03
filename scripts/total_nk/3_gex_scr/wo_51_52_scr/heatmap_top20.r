# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(dplyr)
library(ggplot2)
# ------------------------- #
# Load Top Markers Table using base R
# ------------------------- #
top10_file <- "/home/scripts/total_nk/3_gex_scr/wo_51_52_scr/top10dge.csv"

# top10_file <- "/home/outputs/totalNK_outputs/3_gex/wo_51_52/top10_markers_per_cluster.csv"
top_markers <- read.csv(top10_file, stringsAsFactors = FALSE)
# ------------------------- #
# Extract Unique Gene Names
# ------------------------- #
# Make sure the gene column is correctly named
if (!"gene" %in% colnames(top_markers)) {
 stop("❌ The CSV must have a column named 'gene'")
}
top_genes <- unique(top_markers$gene)
# ------------------------- #
# Optional: Set cluster identity
# ------------------------- #
Idents(seurat_obj) <- "seurat_clusters"
# ------------------------- #
# Generate Heatmap
# ------------------------- #
heatmap_plot <- DoHeatmap(
 object = seurat_obj,
 features = top_genes,
 group.by = "seurat_clusters"
) +
 ggtitle("Top 10 DEG per cluster") +
 theme(plot.title = element_text(hjust = 0.5))
# ------------------------- #
# Save Heatmap to PDF
# ------------------------- #
heatmap_file <- "/home/outputs/totalNK_outputs/3_gex/wo_51_52/pdf/heatmap_top10_cluster_markers.pdf"
ggsave(
 filename = heatmap_file,
 plot = heatmap_plot,
 width = 12,
 height = 16,
 dpi = 600
)
cat("✅ Heatmap saved to:", heatmap_file, "\n")