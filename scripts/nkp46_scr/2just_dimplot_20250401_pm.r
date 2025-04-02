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
# Format: condition_animal (for proper x-axis and y-axis arrangement)
integrated_data$animal_condition <- factor(
  paste(integrated_data$condition, integrated_data$animal, sep = "_"),
  levels = c(
    "nkp46+_Animal25", "nkp46-_Animal25",
    "nkp46+_Animal26", "nkp46-_Animal26",
    "nkp46+_Animal27", "nkp46-_Animal27",
    "nkp46+_Animal28", "nkp46-_Animal28",
    "nkp46+_Animal52", "nkp46-_Animal52"
  )
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
              coord_flip()  # Keeps y-axis as animals and x-axis as conditions
  
  # Save the plot
  output_file <- file.path(dge_output_dir, paste0("dotplot_", marker, "_nkp46_all_animals.pdf"))
  ggsave(filename = output_file, plot = dot_plot, width = 10, height = 8, dpi = 600)
  cat(sprintf("📊 Dot plot for %s saved to %s \n", marker, output_file))
}