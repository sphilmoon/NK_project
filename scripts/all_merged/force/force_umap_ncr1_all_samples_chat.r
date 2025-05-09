# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# ------------------------- #
# Define Paths and Configs
# ------------------------- #
totalNK_rds_file <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/rds/integrated_data_dims25_res0.3_genecounts.rds"
nkp46_rds <- "/home/outputs/nkp46_outputs/nkp46_integrated_data.rds"
output_dir <- "/home/outputs/all_merged_TotalNK_nkp46"
figures_dir <- file.path(output_dir, "figures/forced")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------- #
# Load and Validate Seurat Objects
# ------------------------- #
validate_counts <- function(obj, label) {
  cat(sprintf("📊 %s: %d cells, %d genes\n", label, ncol(obj), nrow(obj)))
}

totalNK <- readRDS(totalNK_rds_file)
nkp46_data <- readRDS(nkp46_rds)
validate_counts(totalNK, "TotalNK")
validate_counts(nkp46_data, "NKp46 Data")

# ------------------------- #
# Subset nkp46+ and nkp46-
# ------------------------- #
nkp46_pos <- subset(nkp46_data, subset = condition == "nkp46+")
nkp46_neg <- subset(nkp46_data, subset = condition == "nkp46-")
validate_counts(nkp46_pos, "NKp46pos")
validate_counts(nkp46_neg, "NKp46neg")

# ------------------------- #
# Annotate Sample Identity
# ------------------------- #
totalNK$sample_id <- "TotalNK"
nkp46_pos$sample_id <- "NKp46pos"
nkp46_neg$sample_id <- "NKp46neg"

# ------------------------- #
# Transfer UMAP from TotalNK
# ------------------------- #
common_genes <- intersect(rownames(totalNK), rownames(nkp46_pos))
ref_features <- common_genes

# Project NKp46+ and NKp46- onto Total NK PCA space, then run UMAP using same dims
shared_features <- intersect(VariableFeatures(totalNK), VariableFeatures(nkp46_pos))

# Scale and PCA using Total NK PCA rotation
# Project NKp46+ onto TotalNK UMAP space
nkp46_pos <- ScaleData(nkp46_pos, features = shared_features, verbose = FALSE)
nkp46_pos <- RunPCA(nkp46_pos, features = shared_features, npcs = 25, reduction.name = "pca", verbose = FALSE)
nkp46_pos <- RunUMAP(nkp46_pos, reduction = "pca", dims = 1:25)
nkp46_pos <- FindNeighbors(nkp46_pos, dims = 1:25)
nkp46_pos <- FindClusters(nkp46_pos, resolution = 0.3)

# Project NKp46- onto TotalNK UMAP space
nkp46_neg <- ScaleData(nkp46_neg, features = shared_features, verbose = FALSE)
nkp46_neg <- RunPCA(nkp46_neg, features = shared_features, npcs = 25, reduction.name = "pca", verbose = FALSE)
nkp46_neg <- RunUMAP(nkp46_neg, reduction = "pca", dims = 1:25)
nkp46_neg <- FindNeighbors(nkp46_neg, dims = 1:25)
nkp46_neg <- FindClusters(nkp46_neg, resolution = 0.3)

# ------------------------- #
# UMAP Visualization by Sample
# ------------------------- #
umap_total <- DimPlot(totalNK, group.by = "seurat_clusters", label = TRUE, pt.size = 0.3) +
  ggtitle("Total NK (Res 0.3)") + theme_minimal()

umap_pos <- DimPlot(nkp46_pos, group.by = "seurat_clusters", label = TRUE, pt.size = 0.3) +
  ggtitle("NKp46+ (Projected to Total NK UMAP)") + theme_minimal()

umap_neg <- DimPlot(nkp46_neg, group.by = "seurat_clusters", label = TRUE, pt.size = 0.3) +
  ggtitle("NKp46- (Projected to Total NK UMAP)") + theme_minimal()

combined_umap <- umap_total + umap_pos + umap_neg + plot_layout(ncol = 3)

# Save combined UMAP
ggsave(
  filename = file.path(figures_dir, "forced_umap_totalNK_nkp46_res0.3.pdf"),
  plot = combined_umap,
  width = 15, height = 5, dpi = 600, bg = "transparent"
)

# ------------------------- #
# NCR1 FeaturePlot (Consistent Scale)
# ------------------------- #
gene <- "NCR1"
all_expr <- c(FetchData(totalNK, vars = gene)[, 1],
              FetchData(nkp46_pos, vars = gene)[, 1],
              FetchData(nkp46_neg, vars = gene)[, 1])
expr_range <- range(all_expr, na.rm = TRUE)

plot_gene <- function(obj, label) {
  FeaturePlot(obj, features = gene, pt.size = 0.3, order = TRUE) +
    scale_color_gradientn(colors = c("lightgrey", "blue", "red"),
                          limits = expr_range,
                          name = gene) +
    ggtitle(label) +
    theme_minimal()
}

ncr1_total <- plot_gene(totalNK, "Total NK")
ncr1_pos   <- plot_gene(nkp46_pos, "NKp46+")
ncr1_neg   <- plot_gene(nkp46_neg, "NKp46-")

combined_featureplot <- ncr1_total + ncr1_pos + ncr1_neg + plot_layout(ncol = 3)

# Save FeaturePlot
ggsave(
  filename = file.path(figures_dir, "NCR1_featureplot_totalNK_nkp46_forced.pdf"),
  plot = combined_featureplot,
  width = 15, height = 5, dpi = 600, bg = "transparent"
)

cat("✅ All UMAP and FeaturePlots generated and saved.\n")
