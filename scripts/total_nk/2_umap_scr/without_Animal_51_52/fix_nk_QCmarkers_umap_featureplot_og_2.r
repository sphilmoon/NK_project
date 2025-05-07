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

# Debug: Check available reductions
cat("🔍 Available reductions:", paste(names(seurat_obj@reductions), collapse = ", "), "\n")

# Debug: Check sample metadata
cat("🔍 Available metadata columns:", paste(colnames(seurat_obj@meta.data), collapse = ", "), "\n")
if ("sample" %in% colnames(seurat_obj@meta.data)) {
  cat("🔍 Unique values in 'sample' column:", paste(unique(seurat_obj$sample), collapse = ", "), "\n")
} else {
  stop("❌ 'sample' column not found in metadata. Available columns: ", paste(colnames(seurat_obj@meta.data), collapse = ", "))
}

# ------------------------- #
# Define Genes and Animals
# ------------------------- #
genes <- c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD40", "CD68", "ENSBTAG00000015032") # Updated gene list
animals <- c("Animal25", "Animal26", "Animal27", "Animal28")

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

# ------------------------- #
# Compute Shared Expression Scale
# ------------------------- #
cat("🔍 Computing shared expression scale...\n")
expr_range <- range(expr_data_all[genes, ], na.rm = TRUE)
cat("Expression range for all genes:", expr_range[1], "to", expr_range[2], "\n")

# ------------------------- #
# Generate Plots by (animal × gene)
# ------------------------- #
cat("🎨 Generating FeaturePlots for genes across animals...\n")
plot_list <- list()

# Define a minimal theme for plots
umap_theme <- theme_minimal() +
  theme(
    plot.title = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

# Generate all individual FeaturePlots
for (animal in animals) {
  row_plots <- list()
  
  for (gene in genes) {
    subset_obj <- subset(seurat_obj, subset = sample == animal)
    
    # Generate FeaturePlot with shared expression scale
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
      limits = expr_range  # Shared scale
    ) +
    umap_theme
    
    row_plots[[gene]] <- p
  }
  
  # Combine row (one row per animal, columns are genes)
  row_patch <- wrap_plots(row_plots, ncol = length(genes))
  plot_list[[animal]] <- row_patch
}

# Combine all rows (animals = rows, genes = columns)
full_grid <- wrap_plots(plot_list, nrow = length(animals))

# Create label plots for animal IDs on the left
label_plots <- lapply(animals, function(a) {
  ggplot() +
    theme_void() +
    labs(title = a) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5, vjust = 0.5),
      plot.margin = margin(0, 0, 0, 10)  # Add margin on the right to align with UMAP plots
    )
})

# Combine into a single list with labels and UMAPs
all_plots <- c(label_plots, plot_list)
n_plots <- length(all_plots)  # 4 labels + 4 * length(genes) UMAPs

# Arrange in a 4x9 grid (4 rows, 9 columns: 1 for labels, 8 for genes)
final_grid <- wrap_plots(all_plots, nrow = 4, ncol = length(genes) + 1)

# Add column labels (genes at the top)
final_grid <- final_grid +
  plot_annotation(
    theme = theme(
      plot.margin = margin(10, 10, 10, 10)
    ),
    # Add column labels (genes at the top, skipping the first column which is labels)
    tag_levels = list(c("", genes)),
    tag_prefix = "",
    tag_suffix = "",
    tag_sep = ""
  ) &
  theme(
    plot.tag = element_text(size = 10, face = "bold", hjust = 0.5),
    plot.tag.position = "top"
  )

# Add legend from a representative plot
legend_plot <- FeaturePlot(
  subset(seurat_obj, subset = sample == "Animal25"),
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

# Extract legend
legend <- cowplot::get_legend(legend_plot)

# Combine grid and legend
final_plot <- final_grid + legend + plot_layout(widths = c(length(genes) + 1, 0.5))

# ------------------------- #
# Save Final Combined Plot
# ------------------------- #
output_file <- file.path(output_dir, "combined_FeaturePlot_animals_rows_genes_cols.pdf")
ggsave(
  filename = output_file,
  plot = final_plot,
  width = 4 * length(genes),  # 8 genes (columns)
  height = 4 * length(animals) + 1,  # 4 animals (rows) + 1 for labels
  dpi = 600,
  bg = "transparent"
)
cat("✅ Combined FeaturePlot saved to", output_file, "\n")