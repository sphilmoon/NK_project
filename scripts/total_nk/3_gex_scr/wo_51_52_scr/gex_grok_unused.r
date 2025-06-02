# Script: gene_expression_analysis.R
# Purpose:
# - Load Seurat object from the specified .rds file
# - Perform differential gene expression (DGE) analysis to identify marker genes per cluster
# - Visualize top markers with a heatmap and feature plots
# - Save results to the specified output directory

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)

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

if (!file.exists(rds_file)) {
  stop("❌ RDS file not found at: ", rds_file)
}

seurat_obj <- readRDS(rds_file)
cat("✅ Loaded Seurat object from", rds_file, "\n")

# Print basic information about the object
cat("Number of cells:", ncol(seurat_obj), "\n")
cat("Number of genes:", nrow(seurat_obj), "\n")
cat("Assays available:", paste(Assays(seurat_obj), collapse = ", "), "\n")
cat("Number of clusters:", length(unique(Idents(seurat_obj))), "\n")

# ------------------------- #
# Gene Expression Analysis
# ------------------------- #

# Set the default assay for DGE analysis
# Using "RNA" assay since the .rds file likely has counts and data layers preserved
# Alternatively, you can use "SCT" if SCTransform normalization is preferred
DefaultAssay(seurat_obj) <- "RNA"
cat("✅ Default assay set to:", DefaultAssay(seurat_obj), "\n")

# Check available layers in the RNA assay
available_layers <- names(seurat_obj[["RNA"]]$layers)
cat("Available layers in RNA assay:", paste(available_layers, collapse = ", "), "\n")

# Ensure the "data" layer (normalized counts) exists for DGE analysis
if (!"data" %in% available_layers) {
  cat("⚠️ 'data' layer not found. Normalizing RNA assay...\n")
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
}

# Identify differentially expressed genes (marker genes) for each cluster
cat("📌 Performing differential gene expression analysis...\n")
dge_markers <- FindAllMarkers(
  object = seurat_obj,
  assay = "RNA",
  slot = "data",  # Use normalized data for DGE
  test.use = "wilcox",  # Wilcoxon test, commonly used for single-cell DGE
  min.pct = 0.1,  # Genes must be detected in at least 10% of cells in the cluster
  logfc.threshold = 0.25,  # Minimum log2 fold change threshold
  verbose = FALSE
)

# Filter for significant markers (adjusted p-value < 0.05)
dge_markers <- dge_markers %>% 
  filter(p_val_adj < 0.05) %>%
  arrange(cluster, desc(avg_log2FC))

# Save the DGE results to CSV
dge_file <- file.path(output_dir, "dge_markers_dims25_res0.3.csv")
write.csv(dge_markers, dge_file, row.names = FALSE)
cat("✅ DGE markers saved to", dge_file, "\n")

# ------------------------- #
# Visualization of Top Markers
# ------------------------- #

# Select top 5 marker genes per cluster (based on log2FC)
top_markers <- dge_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  ungroup()

# Extract unique genes for visualization
top_genes <- unique(top_markers$gene)
cat("Number of top marker genes for visualization:", length(top_genes), "\n")

# Heatmap of top marker genes
cat("📌 Generating heatmap of top marker genes...\n")
heatmap_plot <- DoHeatmap(
  object = seurat_obj,
  features = top_genes,
  assay = "RNA",
  slot = "scale.data",  # Use scaled data for visualization
  group.by = "seurat_clusters",
  label = TRUE
) + 
  ggtitle("Heatmap of Top Marker Genes (dims25, Res 0.3)") +
  theme(plot.title = element_text(hjust = 0.5))

# Save the heatmap
heatmap_file <- file.path(output_dir, "heatmap_top_markers_dims25_res0.3.pdf")
ggsave(heatmap_file, plot = heatmap_plot, width = 12, height = 10, dpi = 600, units = "in")
cat("✅ Heatmap saved to", heatmap_file, "\n")

# Feature plots for top 2 marker genes per cluster
cat("📌 Generating feature plots for top marker genes...\n")
top_2_markers <- dge_markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) %>%
  ungroup()

feature_plots <- lapply(unique(top_2_markers$gene), function(gene) {
  FeaturePlot(
    object = seurat_obj,
    features = gene,
    pt.size = 0.25,
    order = TRUE  # Order cells by expression level for better visualization
  ) +
    ggtitle(paste("Expression of", gene, "(dims25, Res 0.3)")) +
    theme(plot.title = element_text(hjust = 0.5))
})

# Combine feature plots into a single grid
combined_feature_plot <- cowplot::plot_grid(plotlist = feature_plots, ncol = 2)

# Save the feature plots
feature_plot_file <- file.path(output_dir, "feature_plots_top_markers_dims25_res0.3.pdf")
ggsave(feature_plot_file, plot = combined_feature_plot, width = 16, height = 8 * ceiling(length(feature_plots) / 2), dpi = 600, units = "in")
cat("✅ Feature plots saved to", feature_plot_file, "\n")

# ------------------------- #
# Summary of Marker Genes
# ------------------------- #

# Summary of number of significant markers per cluster
marker_summary <- dge_markers %>%
  group_by(cluster) %>%
  summarise(Num_Markers = n()) %>%
  ungroup()

# Save the summary
summary_file <- file.path(output_dir, "marker_summary_dims25_res0.3.csv")
write.csv(marker_summary, summary_file, row.names = FALSE)
cat("✅ Marker summary saved to", summary_file, "\n")

# Print the summary to console
cat("Summary of significant markers per cluster:\n")
print(marker_summary)

# Final completion message
cat("🎉 Gene expression analysis complete. Outputs saved in:\n", output_dir, "\n")