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

# Debug: Check sample metadata
cat("🔍 Available metadata columns:", paste(colnames(seurat_obj@meta.data), collapse = ", "), "\n")
if ("sample" %in% colnames(seurat_obj@meta.data)) {
  cat("🔍 Unique values in 'sample' column:", paste(unique(seurat_obj$sample), collapse = ", "), "\n")
} else {
  stop("❌ 'sample' column not found in metadata.")
}

# ------------------------- #
# Define Genes and Animals
# ------------------------- #
genes <- c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD40", "CD68") #  CD14: "ENSBTAG00000015032"
animals <- c("Animal25", "Animal26", "Animal27", "Animal28")

# Validate animals
missing_animals <- animals[!animals %in% seurat_obj$sample]
if (length(missing_animals) > 0) {
  stop("❌ Animals not found in 'sample' column: ", paste(missing_animals, collapse = ", "))
}
cat("✅ All animals found in 'sample' column.\n")

# Debug: Check cell counts
cat("🔍 Cell counts per animal:\n")
for (animal in animals) {
  n_cells <- sum(seurat_obj$sample == animal)
  cat("  ", animal, ":", n_cells, "cells\n")
  if (n_cells == 0) warning("⚠️ No cells for ", animal)
}

# Validate genes
expr_data_all <- GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = "data")
missing_genes <- genes[!genes %in% rownames(expr_data_all)]
if (length(missing_genes) > 0) {
  stop("❌ Genes not found: ", paste(missing_genes, collapse = ", "))
}
cat("✅ All genes found in assay.\n")

# UMAP check
stopifnot("umap" %in% names(seurat_obj@reductions))

# ------------------------- #
# Compute Shared Expression Scale
# ------------------------- #
cat("🔍 Computing shared expression scale...\n")
expr_range <- range(expr_data_all[genes, ], na.rm = TRUE)
cat("📈 Expression range:", expr_range[1], "to", expr_range[2], "\n")

# ------------------------- #
# Generate FeaturePlots per gene × animal
# ------------------------- #
cat("🎨 Generating FeaturePlots...\n")
plot_grid_list <- list()
umap_theme <- theme_minimal() +
  theme(
    plot.title = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

for (gene in genes) {
  row_plots <- list()
  for (animal in animals) {
    subset_obj <- subset(seurat_obj, subset = sample == animal)
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
    row_plots[[animal]] <- p
  }
  plot_grid_list[[gene]] <- wrap_plots(row_plots, ncol = length(animals))
}

# ------------------------- #
# Create Labeled Plot Grid
# ------------------------- #

# Row gene labels
row_labeled_plots <- list()
for (i in seq_along(genes)) {
  gene_label <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = genes[i], angle = 90, size = 4, fontface = "bold") +
    theme_void()

  row_with_label <- plot_grid(
    wrap_elements(gene_label), plot_grid_list[[genes[i]]],
    ncol = 2,
    rel_widths = c(0.05, 1)
  )
  row_labeled_plots[[i]] <- row_with_label
}
grid_with_row_labels <- wrap_plots(row_labeled_plots, ncol = 1)

# Top column labels
animal_labels <- lapply(animals, function(animal) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = animal, size = 4, fontface = "bold") +
    theme_void()
})
animal_labels_row <- wrap_plots(animal_labels, ncol = length(animals))
animal_labels_row <- plot_grid(NULL, animal_labels_row, ncol = 2, rel_widths = c(0.05, 1))

# Extract legend
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

# Combine all into final plot
final_plot <- plot_grid(
  animal_labels_row,
  plot_grid(grid_with_row_labels, legend, ncol = 2, rel_widths = c(1, 0.07)),
  ncol = 1,
  rel_heights = c(0.05, 1)
)

# ------------------------- #
# Save Final Combined Plot
# ------------------------- #
output_file <- file.path(output_dir, "combined_FeaturePlot_genes_rows_animals_cols.pdf")
ggsave(
  filename = output_file,
  plot = final_plot,
  width = 4 * length(animals) + 1.5,
  height = 4 * length(genes),
  dpi = 600,
  bg = "transparent"
)
cat("✅ Combined FeaturePlot saved to", output_file, "\n")
