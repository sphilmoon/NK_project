# Load libraries
library(Seurat)
library(hdf5r) 
library(dplyr)
library(ggplot2)
library(cowplot) 

# Define directories
data_dir <- "/home/rawdata/nkp46_pos_neg/" 
dge_output_dir <- "/home/outputs/nkp46_outputs" 

# Load your integrated Seurat object
integrated_data <- readRDS(file.path(dge_output_dir, "nkp46_integrated_data.rds"))

# Step 7: Prepare for Dot Plots
DefaultAssay(integrated_data) <- "RNA"
integrated_data <- JoinLayers(integrated_data, assay = "RNA")

# Define animals and markers
animals <- paste0("Animal", c("25", "26", "27", "28", "52"))
markers <- c("NCR1", "CD3E", "KLRB1", "CD2")

# Check and define NKp46 condition if not already present
if (!"condition" %in% colnames(integrated_data@meta.data)) {
  # Define NKp46+ vs NKp46- based on NCR1 expression median
  ncr1_exp <- GetAssayData(integrated_data, assay = "RNA", slot = "data")["NCR1", ]
  integrated_data$condition <- ifelse(ncr1_exp > median(ncr1_exp, na.rm = TRUE), "nkp46+", "nkp46-")
  cat("🔧 Defined NKp46+/- conditions based on NCR1 expression median.\n")
}

# Ensure only desired animals are included
integrated_data <- subset(integrated_data, subset = animal %in% animals)
cat("🌟 Subset data to include only Animals 25, 26, 27, 28, 52.\n")

# Create a combined group for plotting (animal_condition)
# Switch animal and condition to have condition on the x-axis
all_combinations <- expand.grid(condition = c("nkp46+", "nkp46-"), animal = animals)
all_levels <- apply(all_combinations, 1, paste, collapse = "_")

# Generate the animal_condition factor with all combinations as levels
integrated_data$animal_condition <- factor(
  paste(integrated_data$condition, integrated_data$animal, sep = "_"),
  levels = all_levels
)

# Step 8: Generate Dot Plots (One Marker per Figure)
cat("🌟 Generating dot plots for NKp46+ and NKp46- across all animals...\n")
for (marker in markers) {
  dot_plot <- DotPlot(integrated_data, 
                      features = marker, 
                      group.by = "animal_condition", 
                      dot.scale = 8, 
                      cols = c("lightblue", "darkblue")) +
              ggtitle(paste("Expression of", marker, "in NKp46+ and NKp46- Across Animals")) +
              theme(plot.title = element_text(hjust = 0.5),
                    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                    axis.title.x = element_blank()) +
              coord_flip()  # Flip for better readability (condition on x-axis)
  
  # Save the plot
  output_file <- file.path(dge_output_dir, paste0("dotplot_", marker, "_nkp46_all_animals.pdf"))
  ggsave(filename = output_file, plot = dot_plot, width = 10, height = 8, dpi = 600)
  cat(sprintf("📊 Dot plot for %s saved to %s \n", marker, output_file))
}