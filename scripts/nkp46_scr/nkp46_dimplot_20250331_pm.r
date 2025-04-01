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
ggsave(output_file, plot = stacked_elbow_plot, 
       width = 10, height = 6 * length(seurat_list), dpi = 600, units = "in", 
       limitsize = FALSE)
cat("✅ Combined elbow plots saved to", output_file, "\n")



# Step 4: Find Integration Anchors
anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features, 
                                  reduction = "cca", dims = 1:20)
# Step 5: Integrate the Data
integrated_data <- IntegrateData(anchorset = anchors, dims = 1:20)

# Step 6: Dimensionality Reduction and Clustering
DefaultAssay(integrated_data) <- "integrated"
integrated_data <- ScaleData(integrated_data)
integrated_data <- RunPCA(integrated_data, npcs = 30)
integrated_data <- RunUMAP(integrated_data, dims = 1:20)
integrated_data <- FindNeighbors(integrated_data, dims = 1:20)
integrated_data <- FindClusters(integrated_data, resolution = 0.5)

# Save the integrated data object
saveRDS(integrated_data, "nkp46_integrated_data.rds")

# Step 7: Prepare for Dot Plots
DefaultAssay(integrated_data) <- "RNA"
integrated_data <- JoinLayers(integrated_data, assay = "RNA")

# Step 8: Generate Dot Plots for NKp46+ and NKp46- Only
markers <- c("NCR1", "CD3E", "KLRB1", "CD2")  # Marker genes to plot

# Loop through each animal to create a dot plot
for (animal in paste0("Animal", c("25", "26", "27", "28", "52"))) {
  # Subset data for the current animal
  animal_data <- subset(integrated_data, subset = animal == animal)
  
  # Define the plotting group as condition (NKp46+ and NKp46-)
  animal_data$plot_group <- factor(animal_data$condition, levels = c("nkp46+", "nkp46-"))
  
  # Generate the dot plot
  dot_plot <- DotPlot(animal_data, features = markers, 
                      group.by = "plot_group", 
                      dot.scale = 8, 
                      cols = c("lightblue", "darkblue")) +
    ggtitle(paste(animal, "Marker Expression in nkp46+, nkp46-")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save the plot
  ggsave(filename = paste0("dotplot_", animal, "_markers.pdf"), 
         plot = dot_plot, width = 8, height = 6, dpi = 600)
}

# Optional: Combined plot for all animals
dot_plot_all <- DotPlot(integrated_data, features = markers, 
                        group.by = "condition", split.by = "animal",
                        dot.scale = 8, cols = c("lightblue", "darkblue")) +
  ggtitle("Marker Expression Across Animals (nkp46+, nkp46-)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("dotplot_all_animals.pdf", width = 12, height = 8, dpi = 600)