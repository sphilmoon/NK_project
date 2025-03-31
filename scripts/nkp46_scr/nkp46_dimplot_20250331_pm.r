# Load libraries
library(Seurat)
library(hdf5r)
library(dplyr)
library(ggplot2)

# Define the directory containing your .h5 files
data_dir <- "/home/rawdata/nkp46_pos_neg/"  # Replace with your actual path

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
    seurat_obj$animal <- paste0("Animal", animal)
    
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
output_file <- "combined_elbow_plots.pdf"
ggsave(output_file, plot = stacked_elbow_plot, limitsize = FALSE, dpi = 600, units = "in")
cat("✅ Combined elbow plots saved to", output_file, "\n")