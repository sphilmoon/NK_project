# Load libraries
library(Seurat)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(cowplot)

# Define directories
data_dir <- "/home/rawdata/nkp46_pos_neg/"
dge_output_dir <- "/home/outputs/nkp46_outputs"
dimplot_output_dir <- "/home/outputs/nkp46_outputs/dimplot_poster_figures"

# Create output directory if it doesn't exist
dir.create(dge_output_dir, recursive = TRUE, showWarnings = FALSE)

# Load your integrated Seurat object
integrated_data <- readRDS(file.path(dge_output_dir, "nkp46_integrated_data.rds"))

# Step 7: Prepare for Feature Plots
DefaultAssay(integrated_data) <- "RNA"
integrated_data <- JoinLayers(integrated_data, assay = "RNA")

# Define the list of markers to plot
high_markers <- c(
    "CD69", "CD4", "CXCR4", "SELL", "CX3R1", "CD7"
)
  
low_markers <- c(
    "FCER1G", "KIR", "KLR", "CCL5", "ZEB2", "IL2RB", "NCR1", "PRF1", "GZMA"
)

similar_marker <- c(
    "NCR3", "JUNA", "JUNB", "CD52", "IL2RG", "TPT1", "S100A4", "VIM", "CD2"
)

# Define marker groups
marker_groups <- list(
    "high_markers" = high_markers,
    "low_markers" = low_markers,
    "similar_marker" = similar_marker
)

# Step 8: Generate Custom Dot Plots
cat("🌟 Generating custom dot plots for NKp46+ and NKp46- across all animals...\n")

# Loop over each marker group
for (group_name in names(marker_groups)) {
  markers <- marker_groups[[group_name]]
  cat(sprintf("📌 Processing marker group: %s\n", group_name))
  
  # Data frame to store summary data for all markers in the group (for combined plot)
  combined_summary_data <- data.frame()
  
  # Generate individual plot for each marker in the group
  for (marker in markers) {
    # Extract expression data
    expr_data <- GetAssayData(integrated_data, assay = "RNA", slot = "data")[marker, , drop = FALSE]
    
    # Check if the marker exists in the data
    if (nrow(expr_data) == 0) {
      cat(sprintf("⚠️ Marker %s not found in the data, skipping...\n", marker))
      next
    }
    
    # Create a data frame with metadata
    plot_data <- data.frame(
      animal = integrated_data$animal,
      condition = integrated_data$condition,
      expression = as.vector(expr_data[marker, ])
    )
    
    # Calculate average expression and percent expressed per group
    summary_data <- plot_data %>%
      group_by(animal, condition) %>%
      summarise(
        avg_expression = mean(expression, na.rm = TRUE),
        percent_expressed = mean(expression > 0, na.rm = TRUE) * 100,
        .groups = "drop"
      ) %>%
      mutate(marker = marker)  # Add marker name for combined plot
    
    # Add to combined summary data
    combined_summary_data <- rbind(combined_summary_data, summary_data)
    
    # Create the individual dot plot (swapped axes: x = animal, y = condition)
    dot_plot <- ggplot(summary_data, aes(x = animal, y = condition, size = percent_expressed, color = avg_expression)) +
      geom_point() +
      scale_size_continuous(range = c(2, 8), breaks = c(20, 40, 60), name = "Percent Expressed") +
      scale_color_gradient(low = "lightgrey", high = "red", name = "Average Expression") +
      ggtitle(paste("Expression of", marker, "in NKp46+ and NKp46- Across Animals")) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, margin = margin(t = 10)),
        axis.title.y = element_text(size = 12, margin = margin(r = 10)),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
      ) +
      xlab("Animal") +
      ylab("Condition")
    
    # Save the individual plot
    output_file <- file.path(dimplot_output_dir, paste0("dotplot_", group_name, "_", marker, "_nkp46_all_animals.pdf"))
    ggsave(filename = output_file, plot = dot_plot, width = 8, height = 6, dpi = 600)
    cat(sprintf("📊 Individual dot plot for %s in group %s saved to %s \n", marker, group_name, output_file))
  }
  
  # Generate combined plot for all markers in the group (if there are markers to plot)
  if (nrow(combined_summary_data) > 0) {
    # Create the combined dot plot (swapped axes: x = animal, y = condition, facet by marker)
    combined_plot <- ggplot(combined_summary_data, aes(x = animal, y = condition, size = percent_expressed, color = avg_expression)) +
      geom_point() +
      facet_wrap(~ marker, ncol = 3, scales = "free") +  # Facet by marker, 3 columns
      scale_size_continuous(range = c(2, 8), breaks = c(20, 40, 60), name = "Percent Expressed") +
      scale_color_gradient(low = "lightgrey", high = "red", name = "Average Expression") +
      ggtitle(paste("Expression of", group_name, "Markers in NKp46+ and NKp46- Across Animals")) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10, margin = margin(t = 10)),
        axis.title.y = element_text(size = 10, margin = margin(r = 10)),
        strip.text = element_text(size = 10, face = "bold"),  # Marker labels in facets
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.spacing = unit(1, "lines")  # Space between facets
      ) +
      xlab("") +
      ylab("")
    
    # Adjust plot dimensions based on the number of markers
    n_markers <- length(unique(combined_summary_data$marker))
    n_rows <- ceiling(n_markers / 3)  # 3 columns, calculate rows needed
    plot_height <- max(6, n_rows * 2)  # Adjust height based on number of rows
    
    # Save the combined plot
    output_file <- file.path(dimplot_output_dir, paste0("dotplot_", group_name, "_all_markers_nkp46_all_animals.pdf"))
    ggsave(filename = output_file, plot = combined_plot, width = 12, height = plot_height, dpi = 600)
    cat(sprintf("📊 Combined dot plot for group %s saved to %s \n", group_name, output_file))
  } else {
    cat(sprintf("⚠️ No valid markers found for group %s, skipping combined plot...\n", group_name))
  }
}