# Load libraries
library(Seurat)
library(ggplot2)

# Load the integrated Seurat object
integrated_data <- readRDS("/home/outputs/nkp46_outputs/nkp46_integrated_data.rds")

# Set the default assay
DefaultAssay(integrated_data) <- "RNA"

# Define animals and markers
animals <- c("Animal25", "Animal26", "Animal27", "Animal28", "Animal52")
markers <- c("NCR1", "CD3E", "KLRB1", "CD2")

# Define NKp46+ and NKp46- conditions if not already present
if (!"condition" %in% colnames(integrated_data@meta.data)) {
  ncr1_exp <- FetchData(integrated_data, vars = "NCR1")
  integrated_data$condition <- ifelse(ncr1_exp > median(ncr1_exp), "nkp46+", "nkp46-")
}

# Set the factors for animals and conditions
integrated_data$animal <- factor(integrated_data$animal, levels = animals)
integrated_data$condition <- factor(integrated_data$condition, levels = c("nkp46+", "nkp46-"))

# Generate Dot Plots (One Marker per Plot)
for (marker in markers) {
  dot_plot <- DotPlot(
    object = integrated_data,
    features = marker,
    group.by = "condition",     # x-axis: Conditions (NKp46+, NKp46-)
    split.by = "animal"         # y-axis: Animals (Animal25, Animal26, etc.)
  ) +
    ggtitle(paste("Expression of", marker, "in NKp46+ and NKp46- Across Animals")) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14)
    ) +
    scale_color_gradient(low = "lightblue", high = "red") +
    xlab("Condition") +
    ylab("Animal")
  
  # Save the plot
  output_file <- paste0("/home/outputs/nkp46_outputs/dotplot_", marker, "_fixed.pdf")
  ggsave(filename = output_file, plot = dot_plot, width = 8, height = 6)
  cat(sprintf("✅ Dot plot for %s saved to %s\n", marker, output_file))
}