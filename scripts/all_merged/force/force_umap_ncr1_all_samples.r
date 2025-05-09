# Script: umap_ncr1_visualization.R
# Purpose: Visualize UMAP and NCR1 expression for Total NK, NKp46+, and NKp46-

# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# ------------------------- #
# Define Paths and Configs
# ------------------------- #
totalNK_rds_file <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/rds/integrated_data_dims25_res0.3_genecounts.rds"
nkp46_rds <- "/home/outputs/nkp46_outputs/nkp46_integrated_data.rds"
output_dir <- "/home/outputs/all_merged_TotalNK_nkp46/figures/forced"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------- #
# Custom Validation Function
# ------------------------- #
validate_counts <- function(obj, name) {
  if (is.null(obj)) {
    stop(paste("❌", name, "is NULL"))
  }
  n_cells <- ncol(obj)
  cat("✅", name, "contains", n_cells, "cells\n")
  if (n_cells == 0) warning("⚠️", name, "has no cells")
}

# ------------------------- #
# Load Seurat Objects
# ------------------------- #
totalNK <- readRDS(totalNK_rds_file)
cat("✅ Loaded Total NK object from", totalNK_rds_file, "\n")
nkp46_data <- readRDS(nkp46_rds)
cat("✅ Loaded NKp46 object from", nkp46_rds, "\n")

# Validate each object
validate_counts(totalNK, "totalNK")
validate_counts(nkp46_data, "nkp46_data")

# ------------------------- #
# Subset NKp46+ and NKp46- Populations
# ------------------------- #
nkp46_pos <- subset(nkp46_data, subset = condition == "nkp46+")
nkp46_neg <- subset(nkp46_data, subset = condition == "nkp46-")
cat("✅ Subsetted NKp46+ and NKp46- populations\n")

# Validate subsets
validate_counts(nkp46_pos, "nkp46_pos")
validate_counts(nkp46_neg, "nkp46_neg")

# ------------------------- #
# Annotate Sample Identity
# ------------------------- #
totalNK$sample_id <- "TotalNK"
nkp46_pos$sample_id <- "NKp46pos"
nkp46_neg$sample_id <- "NKp46neg"
cat("✅ Annotated sample identities\n")

# ------------------------- #
# UMAP Standardization
# ------------------------- #
# Ensure Total NK has PCA and UMAP with dims = 1:25, resolution = 0.3
if (!"pca" %in% names(totalNK@reductions) || is.null(totalNK@reductions$pca)) {
  totalNK <- RunPCA(totalNK, npcs = 25)
  cat("✅ Recalculated PCA for Total NK with 25 PCs\n")
} else {
  # Check if PCA has enough dimensions
  pca_dims <- dim(totalNK@reductions$pca)[2]
  if (pca_dims < 25) {
    stop("❌ Total NK PCA has only ", pca_dims, " dimensions, but dims = 1:25 requires 25 dimensions")
  }
  cat("✅ Total NK PCA exists with ", pca_dims, " dimensions\n")
}

if (!"umap" %in% names(totalNK@reductions) || is.null(totalNK@reductions$umap)) {
  totalNK <- FindNeighbors(totalNK, dims = 1:25)
  totalNK <- FindClusters(totalNK, resolution = 0.3)
  totalNK <- RunUMAP(totalNK, dims = 1:25)
  cat("✅ Recalculated UMAP for Total NK with dims = 1:25, resolution = 0.3\n")
} else {
  cat("✅ Total NK UMAP already exists, using existing coordinates\n")
}

# Ensure NKp46+ and NKp46- have PCA
if (!"pca" %in% names(nkp46_pos@reductions) || is.null(nkp46_pos@reductions$pca)) {
  nkp46_pos <- RunPCA(nkp46_pos, npcs = 25)
  cat("✅ Recalculated PCA for NKp46+\n")
}
if (!"pca" %in% names(nkp46_neg@reductions) || is.null(nkp46_neg@reductions$pca)) {
  nkp46_neg <- RunPCA(nkp46_neg, npcs = 25)
  cat("✅ Recalculated PCA for NKp46-\n")
}

# Map NKp46+ and NKp46- onto Total NK UMAP
nkp46_pos <- FindNeighbors(nkp46_pos, dims = 1:25)
nkp46_pos <- FindClusters(nkp46_pos, resolution = 0.3)
nkp46_pos <- RunUMAP(nkp46_pos, dims = 1:25)
anchors_pos <- FindTransferAnchors(reference = totalNK, query = nkp46_pos, dims = 1:25, reference.reduction = "pca")
nkp46_pos <- MapQuery(
  anchorset = anchors_pos,
  query = nkp46_pos,
  reference = totalNK,
  refdata = totalNK$seurat_clusters,
  reference.reduction = "umap",
  reduction.model = "umap"
)

nkp46_neg <- FindNeighbors(nkp46_neg, dims = 1:25)
nkp46_neg <- FindClusters(nkp46_neg, resolution = 0.3)
nkp46_neg <- RunUMAP(nkp46_neg, dims = 1:25)
anchors_neg <- FindTransferAnchors(reference = totalNK, query = nkp46_neg, dims = 1:25, reference.reduction = "pca")
nkp46_neg <- MapQuery(
  anchorset = anchors_neg,
  query = nkp46_neg,
  reference = totalNK,
  refdata = totalNK$seurat_clusters,
  reference.reduction = "umap",
  reduction.model = "umap"
)

cat("✅ Projected NKp46+ and NKp46- onto Total NK UMAP space\n")

# ------------------------- #
# Generate UMAP Plots
# ------------------------- #
umap_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

# UMAP plots by sample_id with same scaling
umap_totalNK <- DimPlot(totalNK, group.by = "sample_id", cols = "darkgreen", 
                        reduction = "umap", pt.size = 0.5) + 
                umap_theme + ggtitle("Total NK")

umap_nkp46pos <- DimPlot(nkp46_pos, group.by = "sample_id", cols = "blue", 
                         reduction = "umap", pt.size = 0.5) + 
                 umap_theme + ggtitle("NKp46+")

umap_nkp46neg <- DimPlot(nkp46_neg, group.by = "sample_id", cols = "red", 
                         reduction = "umap", pt.size = 0.5) + 
                 umap_theme + ggtitle("NKp46-")

# Combine UMAP plots
combined_umap_plot <- (umap_totalNK | umap_nkp46pos | umap_nkp46neg) + 
                      plot_layout(ncol = 3)

# ------------------------- #
# NCR1 Expression Visualization
# ------------------------- #
# Check if NCR1 exists
if (!"NCR1" %in% rownames(totalNK)) {
  stop("❌ NCR1 gene not found in Total NK object")
}
if (!"NCR1" %in% rownames(nkp46_pos)) {
  stop("❌ NCR1 gene not found in NKp46+ object")
}
if (!"NCR1" %in% rownames(nkp46_neg)) {
  stop("❌ NCR1 gene not found in NKp46- object")
}

# Get expression range for NCR1 across all objects
ncr1_expr_total <- GetAssayData(totalNK, assay = "RNA", slot = "data")["NCR1", ]
ncr1_expr_pos <- GetAssayData(nkp46_pos, assay = "RNA", slot = "data")["NCR1", ]
ncr1_expr_neg <- GetAssayData(nkp46_neg, assay = "RNA", slot = "data")["NCR1", ]
ncr1_range <- range(c(ncr1_expr_total, ncr1_expr_pos, ncr1_expr_neg), na.rm = TRUE)

# FeaturePlots for NCR1 with same UMAP and color scale
feature_totalNK <- FeaturePlot(totalNK, features = "NCR1", reduction = "umap", 
                               pt.size = 0.5, order = TRUE) + 
                   scale_color_gradientn(colors = c("grey", "blue"), limits = ncr1_range) + 
                   umap_theme + ggtitle("Total NK - NCR1")

feature_nkp46pos <- FeaturePlot(nkp46_pos, features = "NCR1", reduction = "umap", 
                                pt.size = 0.5, order = TRUE) + 
                    scale_color_gradientn(colors = c("grey", "blue"), limits = ncr1_range) + 
                    umap_theme + ggtitle("NKp46+ - NCR1")

feature_nkp46neg <- FeaturePlot(nkp46_neg, features = "NCR1", reduction = "umap", 
                                pt.size = 0.5, order = TRUE) + 
                    scale_color_gradientn(colors = c("grey", "blue"), limits = ncr1_range) + 
                    umap_theme + ggtitle("NKp46- - NCR1")

# Combine FeaturePlots
combined_feature_plot <- (feature_totalNK | feature_nkp46pos | feature_nkp46neg) + 
                         plot_layout(ncol = 3)

# ------------------------- #
# Extract Legend
# ------------------------- #
legend_plot <- FeaturePlot(totalNK, features = "NCR1", reduction = "umap", 
                           pt.size = 0.5) + 
               scale_color_gradientn(colors = c("grey", "blue"), limits = ncr1_range) + 
               theme(legend.position = "right")
legend <- cowplot::get_legend(legend_plot)

# Combine UMAP and FeaturePlots with legend
final_plot <- plot_grid(
  combined_umap_plot, combined_feature_plot,
  ncol = 1, rel_heights = c(1, 1)
) + 
plot_grid(NULL, legend, ncol = 2, rel_widths = c(1, 0.1))

# ------------------------- #
# Save to PDF
# ------------------------- #
output_file <- file.path(output_dir, "UMAP_NCR1_TotalNK_NKp46.pdf")
ggsave(
  filename = output_file,
  plot = final_plot,
  width = 15,
  height = 10,
  dpi = 600,
  bg = "transparent"
)
cat("✅ Final plot saved to", output_file, "\n")

# ------------------------- #
# Final Completion Message
# ------------------------- #
cat("✅ UMAP and NCR1 visualization complete. Output saved in", output_dir, "\n")