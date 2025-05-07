# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

# ------------------------- #
# Define Paths and Configs
# ------------------------- #
rds_file <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/rds/integrated_data_dims25_res0.3_genecounts.rds"
output_dir <- "/home/outputs/totalNK_outputs/2_umap"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------- #
# Load Seurat Object
# ------------------------- #
seurat_obj <- readRDS(rds_file)
cat("✅ Loaded Seurat object from", rds_file, "\n")

# ------------------------- #
# Set Assay
# ------------------------- #
if ("SCT" %in% Assays(seurat_obj)) {
  DefaultAssay(seurat_obj) <- "SCT"
} else {
  DefaultAssay(seurat_obj) <- "RNA"
}
cat("✅ Default assay set to", DefaultAssay(seurat_obj), "\n")

# Ensure layers are joined (Seurat v5)
# seurat_obj <- JoinLayers(seurat_obj, assay = DefaultAssay(seurat_obj))

# ------------------------- #
# Define Genes to Plot
# ------------------------- #
# explore the metadata
cat("🎯 Exploring metadata columns...\n")
print(DefaultAssay(seurat_obj))
# head(GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = "data"))
# head(GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = "counts"))
# head(seurat_obj@meta.data)

genes <- c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD40", "CD68") # CD14 is missing.

# Check if genes exist
expr_data_all <- GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = "data")
cat("Available genes (first 10):", head(rownames(expr_data_all), 10), "\n")
missing_genes <- genes[!genes %in% rownames(expr_data_all)]
if (length(missing_genes) > 0) {
  cat("⚠️ Missing genes:", paste(missing_genes, collapse = ", "), "\n")
  cat("Searching for similar names...\n")
  for (g in missing_genes) {
    similar <- grep(g, rownames(expr_data_all), value = TRUE, ignore.case = TRUE)
    if (length(similar) > 0) cat("  Possible match for", g, ":", paste(similar, collapse = ", "), "\n")
  }
  stop("❌ Genes not found in the Seurat object: ", paste(missing_genes, collapse = ", "))
}
cat("✅ All genes found in the Seurat object: ", paste(genes, collapse = ", "), "\n")

# ------------------------- #
# Check for UMAP and Animal Metadata
# ------------------------- #
# Ensure UMAP is present
if (!"umap" %in% names(seurat_obj@reductions)) {
  stop("❌ UMAP embedding not found in the Seurat object.")
}

# Check for animal metadata (assuming it's in a column named 'sample')
if (!"sample" %in% colnames(seurat_obj@meta.data)) {
  stop("❌ Metadata column 'sample' not found. Please specify the correct column for animals.")
}

# Check number of animals
n_animals <- length(unique(seurat_obj$sample))
if (n_animals != 4) {
  cat("⚠️ Expected 4 animals, but found ", n_animals, ": ", paste(unique(seurat_obj$sample), collapse = ", "), "\n")
} else {
  cat("✅ Found 4 animals: ", paste(unique(seurat_obj$sample), collapse = ", "), "\n")
}

# ------------------------- #
# Create Combined FeaturePlot
# ------------------------- #
cat("🎨 Creating combined FeaturePlot for genes across animals...\n")

# Define a custom theme for consistency
umap_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.text.x = element_text(size = 10, face = "bold"),  # Animal labels
    strip.text.y = element_text(size = 10, face = "bold"),  # Gene labels
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

# Generate FeaturePlot with split.by for animals
combined_feature_plot <- FeaturePlot(
  seurat_obj,
  features = genes,
  split.by = "sample",  # Facet by animal
  pt.size = 0.5,
  order = TRUE,  # Plot cells with higher expression on top
  ncol = 4  # Force 4 columns (one per animal)
) +
  scale_color_gradientn(
    colors = c("lightgrey", "blue", "red"),
    name = "Expression Level"
  ) +
  umap_theme

# ------------------------- #
# Save Combined Plot
# ------------------------- #
featureplot_file <- file.path(output_dir, "NK_QCmarkers_featureplot_by_animal.pdf")
ggsave(
  filename = featureplot_file,
  plot = combined_feature_plot,
  width = 4 * 4,  # 4 animals, 4 inches each
  height = 4 * 6,  # 6 genes, 4 inches each
  dpi = 600,
  bg = "transparent"
)
cat("✅ Combined FeaturePlot saved to", featureplot_file, "\n")