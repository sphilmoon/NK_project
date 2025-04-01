# Load libraries
library(Seurat)
library(hdf5r) 
library(dplyr)
library(ggplot2)
library(cowplot) 

# Define directories
data_dir <- "/home/rawdata/nkp46_pos_neg/" 
dge_output_dir <- "/home/outputs/nkp46_outputs" 

# Create output directory if it doesn’t exist
dir.create(dge_output_dir, recursive = TRUE, showWarnings = FALSE)

# List of animals and conditions
animals <- c("25", "26", "27", "28", "52")
conditions <- c("nkp46+", "nkp46-")

# Initialize a list to store Seurat objects
seurat_list <- list()

# Load each .h5 file and create a Seurat object
for (condition in conditions) {
  for (animal in animals) {
    # Construct file name
    file_name <- paste0(data_dir, condition, "_animal", animal, ".h5")
    
    # Check if file exists
    if (!file.exists(file_name)) {
      warning(paste("File not found:", file_name))
      next
    }
    
    # Load the 10X .h5 file
    data <- Read10X_h5(file_name)
    
    # Create a Seurat object
    seurat_obj <- CreateSeuratObject(counts = data, 
                                     project = paste(condition, "animal", animal, sep = "_"),
                                     min.cells = 3, min.features = 200)
    
    # Add metadata for condition and animal
    seurat_obj$condition <- condition
    seurat_obj$animal <- paste0("Animal", animal)  # Changed to 'animal' from 'sample' later
    
    # Store in the list
    seurat_list[[paste(condition, animal, sep = "_")]] <- seurat_obj
  }
}
# Check the list
print(names(seurat_list))

# Preprocess each Seurat object
seurat_list <- lapply(seurat_list, function(seurat_obj) {
  # Calculate mitochondrial gene percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj <- subset(seurat_obj, 
                       subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  return(seurat_obj)
})
# Check the number of cells after filtering
lapply(seurat_list, function(x) ncol(x))

# Step 3: Scale Data, Run PCA, and Generate Elbow Plots
# Select features for integration
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 2000)
# Initialize a list to store elbow plots
combined_elbow_plots <- list()

# Process each Seurat object: Scale, PCA, and Elbow Plot
for (i in 1:length(seurat_list)) {
  # Get the sample name
  sample_name <- names(seurat_list)[i]
  
  # Scale data and run PCA
  seurat_list[[i]] <- ScaleData(seurat_list[[i]], features = features)
  seurat_list[[i]] <- RunPCA(seurat_list[[i]], features = features)
  
  # Generate and customize the elbow plot
  elbow_plot <- ElbowPlot(seurat_list[[i]], ndims = 50) +
    ggtitle(paste("Elbow Plot for PCA -", sample_name)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Store the plot
  combined_elbow_plots[[sample_name]] <- elbow_plot
  cat("✅ Elbow plot generated for", sample_name, "\n")
}

# Stack the plots in a single PDF
stacked_elbow_plot <- cowplot::plot_grid(plotlist = combined_elbow_plots, ncol = 1)
output_file <- file.path(dge_output_dir, "combined_elbow_plots.pdf")  # Fixed path
ggsave(output_file, plot = stacked_elbow_plot, 
       width = 10, height = 6 * length(seurat_list), dpi = 600, units = "in", 
       limitsize = FALSE)
cat("✅ Combined elbow plots saved to", output_file, "\n")

# Step 4: Find Integration Anchors
dims_to_use <- 1:22  # Adjusted based on elbow plot
anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features, 
                                  reduction = "cca", dims = dims_to_use)

# Step 5: Integrate the Data
integrated_data <- IntegrateData(anchorset = anchors, dims = dims_to_use)

# Step 6: Dimensionality Reduction and Clustering
DefaultAssay(integrated_data) <- "integrated"
integrated_data <- ScaleData(integrated_data)
integrated_data <- RunPCA(integrated_data, npcs = 30)
integrated_data <- RunUMAP(integrated_data, dims = dims_to_use)
integrated_data <- FindNeighbors(integrated_data, dims = dims_to_use)
integrated_data <- FindClusters(integrated_data, resolution = 0.5)

# Save the integrated data object
saveRDS(integrated_data, file.path(dge_output_dir, "nkp46_integrated_data.rds"))

# Step 7: Prepare for Dot Plots
DefaultAssay(integrated_data) <- "RNA"
integrated_data <- JoinLayers(integrated_data, assay = "RNA")

# Define animals and markers
animals <- paste0("Animal", c("25", "26", "27", "28", "52"))
markers <- c("NCR1", "CD3E", "KLRB1", "CD2")  # Marker genes to plot

# Check and define NKp46 condition if not already present
if (!"condition" %in% colnames(integrated_data@meta.data)) {
  # Define NKp46+ vs NKp46- based on NCR1 expression threshold
  ncr1_exp <- GetAssayData(integrated_data, assay = "RNA", slot = "data")["NCR1", ]
  integrated_data$condition <- ifelse(ncr1_exp > median(ncr1_exp, na.rm = TRUE), "nkp46+", "nkp46-")
  cat("🔧 Defined NKp46+/- conditions based on NCR1 expression median.\n")  # Replaced log_print
}

# Ensure only desired animals are included
integrated_data <- subset(integrated_data, subset = animal %in% animals)  # Changed 'sample' to 'animal'
cat("🌟 Subset data to include only Animals 25, 26, 27, 28, 52.\n")

# Create a combined group for plotting (animal_condition)
integrated_data$animal_condition <- factor(paste(integrated_data$animal, integrated_data$condition, sep = "_"),
                                           levels = paste(rep(animals, each = 2), c("nkp46+", "nkp46-"), sep = "_"))

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
                    axis.title.x = element_blank()) +  # Remove x-axis label for clarity
              coord_flip()  # Flip for better readability with many groups
  
  # Save the plot
  output_file <- file.path(dge_output_dir, paste0("dotplot_", marker, "_nkp46_all_animals.pdf"))
  ggsave(filename = output_file, plot = dot_plot, width = 10, height = 8, dpi = 600)
  cat(sprintf("📊 Dot plot for %s saved to %s 🎉\n", marker, output_file))
}