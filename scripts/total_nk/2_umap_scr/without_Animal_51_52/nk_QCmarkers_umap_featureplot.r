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
output_dir <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/pdf"

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
# print(DefaultAssay(seurat_obj))
head(GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = "data"))
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

# # Define a custom theme for consistency
# umap_theme <- theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 12),
#     axis.title = element_blank(),
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     strip.text.x = element_text(size = 10, face = "bold"),  # Animal labels (on top, as columns)
#     strip.text.y = element_text(size = 10, face = "bold"),  # Gene labels (on side, as rows)
#     legend.position = "right",
#     legend.title = element_text(size = 10),
#     legend.text = element_text(size = 8)
#   )

# # Generate FeaturePlot with split.by for animals
# combined_feature_plot <- FeaturePlot(
#   seurat_obj,
#   features = genes,
#   split.by = "sample",  # Facet by animal (4 columns)
#   pt.size = 0.5,
#   order = TRUE,  # Plot cells with higher expression on top
#   ncol = 4  # Force 4 columns (one per animal)
# ) +
#   scale_color_gradientn(
#     colors = c("lightgrey", "blue", "red"),
#     name = "Expression Level"
#   ) +
#   umap_theme 

# # ------------------------- #
# # Save Combined Plot
# # ------------------------- #
# featureplot_file <- file.path(output_dir, "NK_QCmarkers_featureplot_by_animal_flipped.pdf")
# ggsave(
#   filename = featureplot_file,
#   plot = combined_feature_plot,
#   width = 4 * 4,  # 4 animals (columns), 4 inches each
#   height = 4 * 6,  # 6 genes (rows), 4 inches each
#   dpi = 600,
#   bg = "transparent"
# )
# cat("✅ Combined FeaturePlot saved to", featureplot_file, "\n")




# Define a custom theme for consistency
umap_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.text.x = element_text(size = 10, face = "bold"),  # Gene labels (on top, as columns)
    strip.text.y = element_text(size = 10, face = "bold"),  # Animal labels (on side, as rows)
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

# Debug: Check CD68 expression for animal 28
cat("🔍 Checking CD68 expression for animal 28...\n")
animal_28_cells <- which(seurat_obj$sample == "28")
cd68_expr <- GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = "data")["CD68", animal_28_cells]
if (length(cd68_expr) == 0 || all(cd68_expr == 0)) {
  cat("⚠️ No CD68 expression data for animal 28 or all values are 0.\n")
} else {
  cat("✅ CD68 expression data exists for animal 28. Range:", range(cd68_expr, na.rm = TRUE), "\n")
}

# Generate individual FeaturePlots for each gene-animal combination
plot_list <- list()
animals <- unique(seurat_obj$sample)  # Should be 4 animals
# genes <- c("CD3", "CD4", "CD8a", "CD14", "CD40", "CD68")

for (animal in animals) {
  for (gene in genes) {
    # Subset the Seurat object for the current animal
    subset_obj <- subset(seurat_obj, subset = sample == animal)
    
    # Create FeaturePlot for the current gene and animal
    p <- FeaturePlot(
      subset_obj,
      features = gene,
      pt.size = 0.5,
      order = TRUE
    ) +
      scale_color_gradientn(
        colors = c("lightgrey", "blue", "red"),
        name = "Expression Level",
        limits = c(0, max(GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = "data")[genes, ], na.rm = TRUE))  # Shared scale
      ) +
      umap_theme +
      ggtitle(paste(gene, "-", animal))
    
    plot_list[[paste(animal, gene, sep = "_")]] <- p
  }
}

# Arrange plots in a 4x6 grid (4 rows for animals, 6 columns for genes)
plot_grid <- cowplot::plot_grid(
  plotlist = plot_list,
  nrow = 4,  # 4 rows (one per animal)
  ncol = 6,  # 6 columns (one per gene)
  labels = NULL  # Remove default labels; titles already include gene and animal
)

# ------------------------- #
# Save Combined Plot
# ------------------------- #
featureplot_file <- file.path(output_dir, "NK_QCmarkers_featureplot_by_animal_fixed_grid.pdf")
ggsave(
  filename = featureplot_file,
  plot = plot_grid,
  width = 4 * 6,  # 6 genes (columns), 4 inches each
  height = 4 * 4,  # 4 animals (rows), 4 inches each
  dpi = 600,
  bg = "transparent"
)
cat("✅ Combined FeaturePlot saved to", featureplot_file, "\n")