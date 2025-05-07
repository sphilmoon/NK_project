# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
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
DefaultAssay(seurat_obj) <- if ("SCT" %in% Assays(seurat_obj)) "SCT" else "RNA"
cat("✅ Default assay set to", DefaultAssay(seurat_obj), "\n")

# Debug: Check metadata and reductions
cat("🔍 Metadata columns:", paste(colnames(seurat_obj@meta.data), collapse = ", "), "\n")
cat("🔍 Reductions available:", paste(names(seurat_obj@reductions), collapse = ", "), "\n")
stopifnot("sample" %in% colnames(seurat_obj@meta.data), "umap" %in% names(seurat_obj@reductions))

# ------------------------- #
# Define Genes and Animals
# ------------------------- #
genes <- c("CD4", "CD8A", "CD40", "CD68") # "CD3D", "CD3E", "CD3G", "ENSBTAG00000015032 (cd14)" removed
animals <- c("Animal25", "Animal26", "Animal27", "Animal28")


# Check if genes exist
expr_data_all <- GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = "data")
# cat("Available genes (first 10):", head(rownames(expr_data_all), 10), "\n")
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

# Debug: Validate animals exist in the sample column
if ("sample" %in% colnames(seurat_obj@meta.data)) {
  missing_animals <- animals[!animals %in% seurat_obj$sample]
  if (length(missing_animals) > 0) {
    stop("❌ Animals not found in 'sample' column: ", paste(missing_animals, collapse = ", "))
  }
  cat("✅ All animals found in 'sample' column: ", paste(animals, collapse = ", "), "\n")
}

# Debug: Check cell counts per animal
cat("🔍 Cell counts per animal:\n")
for (animal in animals) {
  n_cells <- sum(seurat_obj$sample == animal)
  cat("  Animal", animal, ":", n_cells, "cells\n")
  if (n_cells == 0) {
    warning("⚠️ No cells found for animal ", animal, ". FeaturePlots will be empty for this animal.")
  }
}

# Validate genes
expr_data_all <- GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = "data")
missing_genes <- genes[!genes %in% rownames(expr_data_all)]
if (length(missing_genes) > 0) {
  warning("⚠️ Genes not found: ", paste(missing_genes, collapse = ", "), ". These will be skipped.")
  genes <- genes[genes %in% rownames(expr_data_all)]  # Filter out missing genes
}
cat("✅ Valid genes found: ", paste(genes, collapse = ", "), "\n")

# Check for UMAP and sample metadata
stopifnot("umap" %in% names(seurat_obj@reductions))
stopifnot("sample" %in% colnames(seurat_obj@meta.data))

# Validate animals and genes
missing_animals <- setdiff(animals, unique(seurat_obj$sample))
if (length(missing_animals) > 0) {
  stop("❌ Missing animals in metadata: ", paste(missing_animals, collapse = ", "))
}
expr_data_all <- GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = "data")
missing_genes <- genes[!genes %in% rownames(expr_data_all)]
if (length(missing_genes) > 0) {
  stop("❌ Missing genes in expression data: ", paste(missing_genes, collapse = ", "))
}

cat("✅ All animals and genes found\n")

# ------------------------- #
# Shared Expression Scale
# ------------------------- #
expr_range <- range(expr_data_all[genes, ], na.rm = TRUE)
cat("🔍 Shared expression range:", expr_range[1], "-", expr_range[2], "\n")

# ------------------------- #
# Generate FeaturePlots (Animals = Rows, Genes = Columns)
# ------------------------- #
umap_theme <- theme_minimal() +
  theme(
    plot.title = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

# Create row-wise plot list (one list per animal, each with gene-wise plots)
plot_matrix <- list()

for (animal in animals) {
  row_plots <- list()
  subset_obj <- subset(seurat_obj, subset = sample == animal)

  for (gene in genes) {
    p <- FeaturePlot(
      subset_obj,
      features = gene,
      reduction = "umap",
      pt.size = 0.5,
      order = TRUE
    ) +
    scale_color_gradientn(
      colors = c("lightgrey", "blue"),
      name = "Expression",
      limits = expr_range
    ) +
    umap_theme

    row_plots[[gene]] <- p
  }

  # Combine into one row (genes as columns)
  plot_matrix[[animal]] <- wrap_plots(row_plots, ncol = length(genes))
}

# Combine rows (animals as rows)
full_plot <- wrap_plots(plot_matrix, nrow = length(animals))

# ------------------------- #
# Add Gene Column Labels (Top)
# ------------------------- #
gene_labels <- lapply(genes, function(gene) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = gene, size = 6, fontface = "bold") +
    theme_void()
})
gene_labels_row <- wrap_plots(gene_labels, ncol = length(genes))

# ------------------------- #
# Add Animal Row Labels (Left)
# ------------------------- #
animal_labels <- lapply(animals, function(animal) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = animal, angle = 90, size = 6, fontface = "bold") +
    theme_void()
})
animal_labels_col <- wrap_plots(animal_labels, ncol = 1)

# ------------------------- #
# Extract Legend
# ------------------------- #
legend_plot <- FeaturePlot(
  subset(seurat_obj, subset = sample == animals[1]),
  features = genes[1],
  reduction = "umap",
  pt.size = 0.5
) +
scale_color_gradientn(
  colors = c("lightgrey", "blue"),
  name = "Expression",
  limits = expr_range
) +
theme(legend.position = "right")

legend <- cowplot::get_legend(legend_plot)

# ------------------------- #
# Assemble Final Plot with Labels
# ------------------------- #
grid_with_labels <- plot_grid(
  NULL, gene_labels_row,
  animal_labels_col, full_plot,
  ncol = 2, nrow = 2,
  rel_widths = c(0.05, 1),
  rel_heights = c(0.05, 1)
)

final_plot <- plot_grid(
  grid_with_labels, legend,
  ncol = 2,
  rel_widths = c(1, 0.08)
)

# ------------------------- #
# Save to PDF
# ------------------------- #
output_file <- file.path(output_dir, "FeaturePlot_animals_rows_genes_cols_NOcd3.png")
ggsave(
  filename = output_file,
  plot = final_plot,
  width = 4 * length(genes) + 2,
  height = 4 * length(animals) + 1,
  dpi = 600,
  bg = "transparent"
)
cat("✅ Final plot saved to", output_file, "\n")
