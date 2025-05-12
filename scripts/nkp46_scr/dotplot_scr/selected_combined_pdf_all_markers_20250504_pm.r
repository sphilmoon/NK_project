# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(cowplot)

# ------------------------- #
# Define Input and Output Paths
# ------------------------- #
data_dir <- "/home/rawdata/nkp46_pos_neg/"
dge_output_dir <- "/home/outputs/nkp46_outputs"
dimplot_output_dir <- "/home/outputs/nkp46_outputs/dimplot_poster_figures"

dir.create(dge_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dimplot_output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------- #
# Load Seurat Object
# ------------------------- #
rds_path <- file.path(dge_output_dir, "nkp46_integrated_data.rds")
integrated_data <- readRDS(rds_path)

if (is.null(integrated_data)) stop("❌ Failed to load Seurat object from ", rds_path)
cat("✅ Loaded Seurat object from", rds_path, "\n")

# ------------------------- #
# Preprocessing
# ------------------------- #
DefaultAssay(integrated_data) <- "RNA"

# Join layers if not already joined
if ("counts" %in% names(integrated_data@assays$RNA@layers)) {
  integrated_data <- JoinLayers(integrated_data, assay = "RNA")
  cat("✅ Layers joined in RNA assay\n")
}

# Remove Animal52 from metadata and object
if ("animal" %in% colnames(integrated_data@meta.data)) {
  integrated_data <- subset(integrated_data, subset = animal != "Animal52")
  cat("✅ Removed Animal52 from dataset\n")
}

# ------------------------- #
# Define Marker Groups
# ------------------------- #
high_markers <- c("CD69", "CD4", "CXCR4", "SELL", "CX3CR1", "CD7")
low_markers <- c("FCER1G", "KIR2DL1", "KLRD1", "CCL5", "ZEB2", "IL2RB", "NCR1", "PRF1", "GZMA")
similar_marker <- c("NCR3", "JUN", "JUNB", "CD52", "IL2RG", "TPT1", "S100A4", "VIM", "CD2")

marker_groups <- list(
  high_markers = high_markers,
  low_markers = low_markers,
  similar_marker = similar_marker
)

# ------------------------- #
# Generate Combined Dot Plots
# ------------------------- #
cat("🌟 Generating combined dot plots...\n")

for (group_name in names(marker_groups)) {
  cat(sprintf("📌 Processing marker group: %s\n", group_name))
  markers <- marker_groups[[group_name]]
  combined_summary_data <- data.frame()

  for (marker in markers) {
    cat(sprintf("  📥 Processing marker: %s\n", marker))

    all_genes <- rownames(integrated_data)
    if (!(marker %in% all_genes)) {
      cat(sprintf("  ⚠️ Marker %s not found in RNA assay. Skipping...\n", marker))
      next
    }

    expr_vec <- FetchData(integrated_data, vars = marker)[, 1]

    if (!all(c("animal", "condition") %in% colnames(integrated_data@meta.data))) {
      cat("  ⚠️ Metadata 'animal' or 'condition' missing. Using cluster ID instead.\n")
      plot_data <- data.frame(
        group = as.character(Idents(integrated_data)),
        expression = expr_vec
      )

      summary_data <- plot_data %>%
        group_by(group) %>%
        summarise(
          avg_expression = mean(expression, na.rm = TRUE),
          percent_expressed = mean(expression > 0, na.rm = TRUE) * 100,
          .groups = "drop"
        ) %>%
        mutate(marker = marker)
    } else {
      plot_data <- data.frame(
        animal = integrated_data$animal,
        condition = integrated_data$condition,
        expression = expr_vec
      )

      summary_data <- plot_data %>%
        group_by(animal, condition) %>%
        summarise(
          avg_expression = mean(expression, na.rm = TRUE),
          percent_expressed = mean(expression > 0, na.rm = TRUE) * 100,
          .groups = "drop"
        ) %>%
        mutate(marker = marker)
    }

    combined_summary_data <- bind_rows(combined_summary_data, summary_data)
  }

  # ------------------------- #
  # Create and Save Dot Plot
  # ------------------------- #
  if (nrow(combined_summary_data) == 0) {
    cat(sprintf("  ⚠️ No valid data found for group %s. Skipping plot.\n", group_name))
    next
  }

  n_markers <- length(unique(combined_summary_data$marker))
  n_rows <- ceiling(n_markers / 3)
  plot_height <- max(6, n_rows * 2)

  # Custom color scale: 0 = white, >0 = gradient
  combined_summary_data$expression_bin <- ifelse(combined_summary_data$avg_expression == 0, "zero", "nonzero")

  if ("animal" %in% colnames(combined_summary_data)) {
    combined_plot <- ggplot(combined_summary_data, aes(x = animal, y = condition, size = percent_expressed)) +
      geom_point(aes(color = avg_expression)) +
      facet_wrap(~ marker, ncol = 3, scales = "free") +
      scale_size_continuous(range = c(2, 8), breaks = c(20, 40, 60), name = "Percent Expressed") +
      scale_color_gradientn(colors = c("white", "lightgray", "red"), values = scales::rescale(c(0, 0.001, max(combined_summary_data$avg_expression, na.rm = TRUE))),
                            name = "Average Expression") +
      ggtitle(paste("Expression of", group_name, "Markers in NKp46+ and NKp46-")) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "right"
      ) +
      xlab("Animal") +
      ylab("Condition")
  } else {
    combined_plot <- ggplot(combined_summary_data, aes(x = group, y = 1, size = percent_expressed)) +
      geom_point(aes(color = avg_expression)) +
      facet_wrap(~ marker, ncol = 3, scales = "free") +
      scale_size_continuous(range = c(2, 8), breaks = c(20, 40, 60), name = "Percent Expressed") +
      scale_color_gradientn(colors = c("white", "lightgray", "red"), values = scales::rescale(c(0, 0.001, max(combined_summary_data$avg_expression, na.rm = TRUE))),
                            name = "Average Expression") +
      ggtitle(paste("Expression of", group_name, "Markers Across Clusters")) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "right"
      ) +
      xlab("Cluster") +
      ylab("")
  }

  output_file <- file.path(dimplot_output_dir, paste0("dotplot_", group_name, "_nkp46_markers_filtered.pdf"))
  ggsave(output_file, plot = combined_plot, width = 12, height = plot_height, dpi = 600)
  cat(sprintf("📊 Dot plot saved to: %s\n", output_file))
}

cat("🎉 All dot plots generated successfully.\n")