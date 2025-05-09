# Script: umap_ncr1_by_animal_sample.R
# Purpose: Visualize UMAP and NCR1 expression by animal and sample type

# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(ggplot2)
library(cowplot)

# ------------------------- #
# Define Paths and Configs
# ------------------------- #
output_dir <- "/home/outputs/all_merged_TotalNK_nkp46"

# Create output subdirectories
rds_dir <- file.path(output_dir, "rds")
figures_dir <- file.path(output_dir, "figures")
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------- #
# Load Merged Seurat Object
# ------------------------- #
merged_obj <- readRDS(file.path(rds_dir, "all_totalNK_nkp46_dim25_res0.5_clustered.rds"))
cat("✅ Loaded merged Seurat object\n")

# ------------------------- #
# Annotate Metadata
# ------------------------- #
# Debug: Inspect cell barcodes and current metadata
cat("🔍 First 10 cell barcodes:", head(colnames(merged_obj), 10), "\n")
cat("🔍 Available metadata columns:", paste(colnames(merged_obj@meta.data), collapse = ", "), "\n")

# Attempt to assign animal_id from cell barcodes
merged_obj$animal_id <- sapply(strsplit(colnames(merged_obj), "_"), function(x) {
  if (length(x) > 1) x[2] else NA
})
cat("🔍 Assigned animal_id values (first 10):", head(merged_obj$animal_id, 10), "\n")

# Fallback to 'sample' column if animal_id is invalid
if (all(is.na(merged_obj$animal_id)) || !any(merged_obj$animal_id %in% c("Animal25", "Animal26", "Animal27", "Animal28"))) {
  if ("sample" %in% colnames(merged_obj@meta.data)) {
    cat("⚠️ animal_id extraction failed. Attempting to use 'sample' column...\n")
    merged_obj$animal_id <- merged_obj$sample
    cat("🔍 Updated animal_id from 'sample' column (first 10):", head(merged_obj$animal_id, 10), "\n")
  } else {
    stop("❌ No valid method to assign animal_id. Check metadata or cell barcode format.")
  }
}

if (!"sample_id" %in% colnames(merged_obj@meta.data)) {
  stop("❌ 'sample_id' not found in metadata. Ensure it is set before merging.")
}

# ------------------------- #
# Define Plotting Order
# ------------------------- #
animals <- c("Animal25", "Animal26", "Animal27", "Animal28")
samples <- c("NKp46neg", "NKp46pos", "TotalNK")  # Order: NKp46-, NKp46+, Total NK

# Validate and debug animal and sample presence
cat("🔍 Unique animal_id values:", unique(merged_obj$animal_id), "\n")
cat("🔍 Unique sample_id values:", unique(merged_obj$sample_id), "\n")

missing_animals <- animals[!animals %in% unique(merged_obj$animal_id)]
if (length(missing_animals) > 0) {
  warning("⚠️ Missing animals in metadata: ", paste(missing_animals, collapse = ", "), 
          ". Using available animals: ", paste(intersect(animals, unique(merged_obj$animal_id)), collapse = ", "))
  animals <- intersect(animals, unique(merged_obj$animal_id))
  if (length(animals) == 0) stop("❌ No valid animals found in metadata.")
}
missing_samples <- samples[!samples %in% unique(merged_obj$sample_id)]
if (length(missing_samples) > 0) {
  stop("❌ Missing samples in metadata: ", paste(missing_samples, collapse = ", "))
}

# Define a minimal theme for plots
umap_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

# ------------------------- #
# Objective 1: UMAP by Animal × Sample
# ------------------------- #
umap_plot_list <- list()

for (animal in animals) {
  row_plots <- list()
  for (sample in samples) {
    subset_cells <- WhichCells(merged_obj, expression = animal_id == animal & sample_id == sample)
    cat("🔍 Cells for", animal, "and", sample, ":", length(subset_cells), "\n")
    if (length(subset_cells) == 0) {
      warning("⚠️ No cells found for ", animal, " and ", sample, ". Creating empty plot.")
      p <- ggplot() + 
           geom_blank() + 
           ggtitle(paste(animal, sample)) +
           umap_theme +
           coord_fixed(ratio = 1)  # Maintain aspect ratio for empty plot
    } else {
      # Define colors for each sample type
      sample_color <- switch(sample,
                             "NKp46neg" = "red",
                             "NKp46pos" = "blue",
                             "TotalNK" = "darkgreen")
      p <- DimPlot(merged_obj, cells = subset_cells, group.by = "sample_id", 
                   cols = sample_color, reduction = "umap", pt.size = 0.5) +
           ggtitle(paste(animal, sample)) +
           umap_theme
    }
    row_plots[[sample]] <- p
  }
  # Combine row (one row per animal, columns are samples)
  row_patch <- wrap_plots(row_plots, ncol = length(samples))
  umap_plot_list[[animal]] <- row_patch
}

# Combine all rows (animals = rows, samples = columns)
umap_full_grid <- wrap_plots(umap_plot_list, nrow = length(animals))

# Save UMAP plot
umap_output_file <- file.path(figures_dir, "UMAP_by_animal_and_sample.pdf")
ggsave(
  filename = umap_output_file,
  plot = umap_full_grid,
  width = 12,  # 3 samples × 4 inches
  height = 16, # 4 animals × 4 inches
  dpi = 600,
  bg = "transparent"
)
cat("✅ Saved UMAP plot by animal and sample to", umap_output_file, "\n")

# ------------------------- #
# Objective 2: NCR1 Expression on UMAP
# ------------------------- #
# Validate NCR1 presence
if (!"NCR1" %in% rownames(merged_obj)) {
  stop("❌ NCR1 gene not found in the Seurat object")
}

# Compute global expression range for NCR1
expr_vals <- FetchData(merged_obj, vars = "NCR1")[, 1]
expr_range <- range(expr_vals, na.rm = TRUE)
cat("🔍 NCR1 expression range:", expr_range[1], "to", expr_range[2], "\n")

# Generate FeaturePlots
feature_plot_list <- list()

for (animal in animals) {
  row_plots <- list()
  for (sample in samples) {
    subset_cells <- WhichCells(merged_obj, expression = animal_id == animal & sample_id == sample)
    cat("🔍 Cells for", animal, "and", sample, ":", length(subset_cells), "\n")
    if (length(subset_cells) == 0) {
      warning("⚠️ No cells found for ", animal, " and ", sample, ". Creating empty plot.")
      p <- ggplot() + 
           geom_blank() + 
           ggtitle(paste(animal, sample)) +
           umap_theme +
           coord_fixed(ratio = 1)  # Maintain aspect ratio for empty plot
    } else {
      p <- FeaturePlot(merged_obj, features = "NCR1", cells = subset_cells, 
                       reduction = "umap", pt.size = 0.5, order = TRUE) +
           scale_color_gradientn(
             colors = c("lightgrey", "blue"),
             limits = expr_range,
             name = "NCR1 Expression"
           ) +
           ggtitle(paste(animal, sample)) +
           umap_theme
    }
    row_plots[[sample]] <- p
  }
  # Combine row (one row per animal, columns are samples)
  row_patch <- wrap_plots(row_plots, ncol = length(samples))
  feature_plot_list[[animal]] <- row_patch
}

# Combine all rows (animals = rows, samples = columns)
feature_full_grid <- wrap_plots(feature_plot_list, nrow = length(animals))

# Extract a single legend for NCR1 expression
legend_plot <- FeaturePlot(merged_obj, features = "NCR1", 
                           reduction = "umap", pt.size = 0.5, order = TRUE) +
               scale_color_gradientn(
                 colors = c("lightgrey", "blue"),
                 limits = expr_range,
                 name = "NCR1 Expression"
               ) +
               theme(legend.position = "right")
legend <- cowplot::get_legend(legend_plot)

# Combine grid with legend
final_feature_plot <- plot_grid(
  feature_full_grid, legend,
  ncol = 2,
  rel_widths = c(1, 0.1)
)

# Save NCR1 FeaturePlot
ncr1_output_file <- file.path(figures_dir, "NCR1_expression_by_animal_and_sample.pdf")
ggsave(
  filename = ncr1_output_file,
  plot = final_feature_plot,
  width = 13,  # Extra width for legend
  height = 16,
  dpi = 600,
  bg = "transparent"
)
cat("✅ Saved NCR1 FeaturePlot by animal and sample to", ncr1_output_file, "\n")

# ------------------------- #
# Final Completion Message
# ------------------------- #
cat("✅ Visualization complete. Outputs saved in", figures_dir, "\n")