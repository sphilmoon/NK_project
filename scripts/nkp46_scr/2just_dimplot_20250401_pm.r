# Load libraries
library(Seurat)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(cowplot)

# Define directories
data_dir <- "/home/rawdata/nkp46_pos_neg/"
dge_output_dir <- "/home/outputs/nkp46_outputs"

# Create output directory if it doesn't exist
dir.create(dge_output_dir, recursive = TRUE, showWarnings = FALSE)

# Load your integrated Seurat object
integrated_data <- readRDS(file.path(dge_output_dir, "nkp46_integrated_data.rds"))

# Step 7: Prepare for Dot Plots
DefaultAssay(integrated_data) <- "RNA"
integrated_data <- JoinLayers(integrated_data, assay = "RNA")

# Define animals and markers
animals <- paste0("Animal", c("25", "26", "27", "28", "52"))
markers <- c("NCR1", "CD3E", "KLRB1", "CD2")  # Marker genes to plot

# Check and define NKp46 condition if not already present
if (!"condition" %in% colnames(integrated_data@meta.data)) {
  ncr1_exp <- GetAssayData(integrated_data, assay = "RNA", slot = "data")["NCR1", ]
  integrated_data$condition <- ifelse(ncr1_exp > median(ncr1_exp, na.rm = TRUE), "nkp46+", "nkp46-")
  cat("🔧 Defined NKp46+/- conditions based on NCR1 expression median.\n")
}

# Ensure only desired animals are included
integrated_data <- subset(integrated_data, subset = animal %in% animals)
cat("🌟 Subset data to include only Animals 25, 26, 27, 28, 52.\n")

# Step 8: Generate Custom Dot Plots (One Marker per Figure)
cat("🌟 Generating custom dot plots for NKp46+ and NKp46- across all animals...\n")
for (marker in markers) {
  # Extract expression data
  expr_data <- GetAssayData(integrated_data, assay = "RNA", slot = "data")[marker, , drop = FALSE]
  
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
    )
  
  # Create the dot plot using ggplot2
  dot_plot <- ggplot(summary_data, aes(x = animal, y = condition, size = percent_expressed, color = avg_expression)) +
    geom_point() +
    scale_size_continuous(range = c(2, 8), breaks = c(20, 40, 60)) +
    scale_color_gradient(low = "lightgrey", high = "red") +
    ggtitle(paste("Expression of", marker, "in NKp46+ and NKp46- Across Animals")) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12, margin = margin(t = 10)),
      axis.title.y = element_text(size = 12, margin = margin(r = 10)),
      legend.position = "right"
    ) +
    xlab("Animal") +
    ylab("Condition")
  
  # Save the plot
  output_file <- file.path(dge_output_dir, paste0("dotplot_", marker, "_nkp46_all_animals.pdf"))
  ggsave(filename = output_file, plot = dot_plot, width = 8, height = 10, dpi = 600)
  cat(sprintf("📊 Dot plot for %s saved to %s \n", marker, output_file))
}