# Load libraries
library(Seurat)
library(hdf5r) 
library(dplyr)
library(ggplot2)
library(cowplot) 

# Define directories
data_dir <- "/home/rawdata/nkp46_pos_neg/" 
dge_output_dir <- "/home/outputs/nkp46_outputs" 

# Load the integrated Seurat object
integrated_data <- readRDS(file.path(dge_output_dir, "nkp46_integrated_data.rds"))

# Set the default assay
DefaultAssay(integrated_data) <- "RNA"

# Define animals and markers
animals <- paste0("Animal", c("25", "26", "27", "28", "52"))
markers <- c("NCR1", "CD3E", "KLRB1", "CD2")

# Define NKp46+ and NKp46- conditions if not already present
if (!"condition" %in% colnames(integrated_data@meta.data)) {
  ncr1_exp <- GetAssayData(integrated_data, assay = "RNA", slot = "data")["NCR1", ]
  integrated_data$condition <- ifelse(ncr1_exp > median(ncr1_exp, na.rm = TRUE), "nkp46+", "nkp46-")
}

# Ensure only the desired animals are included
integrated_data <- subset(integrated_data, subset = animal %in% animals)

# Set the factors for plotting
integrated_data$animal <- factor(integrated_data$animal, levels = animals)
integrated_data$condition <- factor(integrated_data$condition, levels = c("nkp46+", "nkp46-"))

# Generate Dot Plots (One Marker per Plot)
for (marker in markers) {
  dot_plot <- DotPlot(
    object = integrated_data,
    features = marker,
    group.by = "condition",    # x-axis: Conditions (NKp46+, NKp46-)
    split.by = "animal",       # y-axis: Animals (Animal25, Animal26, etc.)
    dot.scale = 8,
    cols = c("lightblue", "blue", "darkblue", "purple", "black") # Color gradient for expression
  ) +
    ggtitle(paste("Expression of", marker, "in NKp46+ and NKp46- Across Animals")) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    ) +
    xlab("Condition") +
    ylab("Animal")
  
  # Save the plot
  output_file <- file.path(dge_output_dir, paste0("dotplot_", marker, "_nkp46_conditions_fixed.pdf"))
  ggsave(filename = output_file, plot = dot_plot, width = 10, height = 8, dpi = 300)
  cat(sprintf("✅ Dot plot for %s saved to %s \n", marker, output_file))
}