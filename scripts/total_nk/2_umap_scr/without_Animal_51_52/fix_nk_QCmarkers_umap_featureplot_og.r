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

# # Debug: Join layers if necessary (Seurat v5)
# if (length(Layers(seurat_obj, assay = DefaultAssay(seurat_obj))) > 1) {
#   cat("🔍 Joining layers for assay", DefaultAssay(seurat_obj), "\n")
#   seurat_obj <- JoinLayers(seurat_obj, assay = DefaultAssay(seurat_obj))
# }

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
genes <- c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD40", "CD68")
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
  stop("❌ Genes not found: ", paste(missing_genes, collapse = ", "))
}
cat("✅ All genes found: ", paste(genes, collapse = ", "), "\n")

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
# Generate Plots by (gene × animal)
# ------------------------- #
cat("🎨 Generating FeaturePlots for genes across animals...\n")
plot_grid_list <- list()

# Define a minimal theme for plots
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
    
    row_plots[[animal]] <- p
  }
  
  # Combine row (one row per gene, columns are animals)
  row_patch <- wrap_plots(row_plots, ncol = length(animals))
  plot_grid_list[[gene]] <- row_patch
}

# Combine all rows (genes = rows, animals = columns)
full_grid <- wrap_plots(plot_grid_list, nrow = length(genes))

# Add legend from a representative plot
legend_plot <- FeaturePlot(
  subset(seurat_obj, subset = sample == "Animal25"),
  features = genes[1],
  reduction = "umap",
  pt.size = 0.5
) +
scale_color_gradientn(
  colors = c("lightgrey", "blue", "red"),
  name = "Expression",
  limits = expr_range
) +
theme(legend.position = "right")

# Extract legend
legend <- cowplot::get_legend(legend_plot)

# Combine grid with labels and legend
final_plot <- full_grid +
  plot_layout(guides = "collect") +
  plot_annotation(
    theme = theme(
      plot.margin = margin(10, 10, 10, 10)
    ),
    # Add column labels (animals at the top)
    tag_levels = list(animals),
    tag_prefix = "",
    tag_suffix = "",
    tag_sep = ""
  ) &
  theme(
    legend.position = "right",
    plot.tag = element_text(size = 10, face = "bold", hjust = 0.5),
    plot.tag.position = "top"
  )

# Add row labels (genes on the left) manually by adjusting the plot layout
# Since patchwork doesn't directly support row labels, we'll use a workaround with plot titles
for (i in seq_along(genes)) {
  final_plot[[i]] <- final_plot[[i]] + plot_annotation(
    title = genes[i],
    theme = theme(
      plot.title = element_text(size = 10, face = "bold", hjust = -0.1, vjust = 0.5)
    )
  )
}

# ------------------------- #
# Save Final Combined Plot
# ------------------------- #
output_file <- file.path(output_dir, "combined_FeaturePlot_genes_rows_animals_cols.pdf")
ggsave(
  filename = output_file,
  plot = final_plot,
  width = 4 * length(animals),  # 4 animals (columns)
  height = 4 * length(genes),   # 7 genes (rows)
  dpi = 600,
  bg = "transparent"
)
cat("✅ Combined FeaturePlot saved to", output_file, "\n")