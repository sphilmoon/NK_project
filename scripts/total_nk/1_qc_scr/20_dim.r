# Load libraries
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(patchwork)
library(cowplot)
library(RColorBrewer)
library(pheatmap)
library(cluster)
 
# Define file paths and sample names
file_paths <- list(
  "/home/rawdata/total_nk/animal25_totalnk_filtered_feature_bc_matrix.h5",
  "/home/rawdata/total_nk/animal27_totalnk_filtered_feature_bc_matrix.h5",
  "/home/rawdata/total_nk/animal51_totalnk_filtered_feature_bc_matrix.h5",
  "/home/rawdata/total_nk/animal26_totalnk_filtered_feature_bc_matrix.h5",
  "/home/rawdata/total_nk/animal28_totalnk_filtered_feature_bc_matrix.h5",
  "/home/rawdata/total_nk/animal52_totalnk_filtered_feature_bc_matrix.h5"
)
sample_names <- c("Animal25", "Animal27", "Animal51", "Animal26", "Animal28", "Animal52")
 
# 1. Load individual H5 files into Seurat objects
seurat_objects <- mapply(function(file, name) {
  data <- Read10X_h5(file)
  seurat_obj <- CreateSeuratObject(counts = data, project = name, min.cells = 3, min.features = 200)
  seurat_obj$sample <- name  # Add sample metadata
  return(seurat_obj)
}, file_paths, sample_names, SIMPLIFY = FALSE)
 
# Assign meaningful names to the list
names(seurat_objects) <- sample_names
 
# Print sample information with details
for (name in names(seurat_objects)) {
  cat("Sample:", name, "\n")
  print(seurat_objects[[name]])
  cat("Number of cells:", ncol(seurat_objects[[name]]), "\n")
  cat("Number of genes:", nrow(seurat_objects[[name]]), "\n\n")
}
 

# 2. quality control and pre-processing

# Verify the loaded objects (optional)
cat("Number of samples loaded:", length(seurat_objects), "\n")
lapply(names(seurat_objects), function(name) {
  cat("Sample:", name, " - Cells:", ncol(seurat_objects[[name]]), " Genes:", nrow(seurat_objects[[name]]), "\n")
})

# to see the gene names
head(rownames(seurat_objects[[1]]), 20)

# Create a list to store all combined plots
combined_plots <- list()

# Perform QC filtering and generate violin plots
for (i in 1:length(seurat_objects)) {
  # Define animal name
  animal_name <- names(seurat_objects)[i]
  
  # Apply QC filtering based on animal name
  if (animal_name == "Animal27") {
    seurat_objects[[i]] <- subset(seurat_objects[[i]], 
                                  subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & 
                                           nCount_RNA > 250 & nCount_RNA < 7500)
  } else if (animal_name == "Animal28") {
    seurat_objects[[i]] <- subset(seurat_objects[[i]], 
                                  subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & 
                                           nCount_RNA > 250 & nCount_RNA < 10000)
  } else if (animal_name == "Animal51") {
    seurat_objects[[i]] <- subset(seurat_objects[[i]], 
                                  subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & 
                                           nCount_RNA > 250 & nCount_RNA < 5000)
  } else if (animal_name == "Animal52") {
    seurat_objects[[i]] <- subset(seurat_objects[[i]], 
                                  subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & 
                                           nCount_RNA > 250 & nCount_RNA < 7500)
  } else {
    # Default threshold for other animals (e.g., Animal25, Animal26)
    seurat_objects[[i]] <- subset(seurat_objects[[i]], 
                                  subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & 
                                           nCount_RNA > 500 & nCount_RNA < 10000)
  }
  
  # Create separate violin plots for each feature, removing the legend
  p1 <- VlnPlot(seurat_objects[[i]], features = "nFeature_RNA", ncol = 1) +
        ggtitle("Gene counts QC") +
        theme(plot.title = element_text(hjust = 0.5),  # Center the title
              legend.position = "none")               # Remove legend
  
  p2 <- VlnPlot(seurat_objects[[i]], features = "nCount_RNA", ncol = 1) +
        ggtitle("UMI counts QC") +
        theme(plot.title = element_text(hjust = 0.5),  # Center the title
              legend.position = "none")               # Remove legend
  
  # Combine plots for this sample using patchwork
  combined_plot <- p1 + p2 + plot_layout(ncol = 2)
  
  # Add an overall title for the combined plot
  combined_plot <- combined_plot + plot_annotation(title = paste("QC for", animal_name))
  
  # Store the combined plot in the list
  combined_plots[[animal_name]] <- combined_plot
  
  # Print status message for this sample
  cat("✅ Plot generated for", animal_name, "\n")
}

# Stack all plots vertically using cowplot
stacked_plot <- cowplot::plot_grid(plotlist = combined_plots, ncol = 1)

# Normalize with SCTransform
seurat_objects <- lapply(seurat_objects, SCTransform, verbose = FALSE)
 
# Final completion message
cat("QC and figure generation process is finished.\n")


# 3. Batch correction and integration

# Step 1: Run PCA on the first sample to determine dims
# Create a list to store all elbow plots
combined_elbow_plots <- list()

for (i in 1:length(seurat_objects)) {
  # Define animal name
  animal_name <- names(seurat_objects)[i]
  
  # Run PCA on the current sample
  seurat_objects[[i]] <- RunPCA(seurat_objects[[i]], verbose = FALSE)
  
  # Generate and customize the elbow plot
  elbow_plot <- ElbowPlot(seurat_objects[[i]], ndims = 50) +
                ggtitle(paste("Elbow Plot for PCA -", animal_name)) +  # Dynamic title
                theme(plot.title = element_text(hjust = 0.5))  # Center the title
  
  # Store the plot in the list
  combined_elbow_plots[[animal_name]] <- elbow_plot
  
  # Print status message
  cat("✅ Elbow plot generated for", animal_name, "\n")
}

# Step 2: Determine the number of dimensions for integration
dims_to_use <- 1:20  # Replace with your chosen dims

# Step 3: Select integration features
features <- SelectIntegrationFeatures(object.list = seurat_objects, nfeatures = 3000)

# Step 4: Find integration anchors
anchors <- FindIntegrationAnchors(object.list = seurat_objects, anchor.features = features, dims = dims_to_use)

# Step 5: Integrate data
integrated_data20 <- IntegrateData(anchorset = anchors, dims = dims_to_use)

# 4. Clustering and visualizing the integration
DefaultAssay(integrated_data20) <- "integrated"
integrated_data20 <- ScaleData(integrated_data20)
integrated_data20 <- RunPCA(integrated_data20)
integrated_data20 <- RunUMAP(integrated_data20, dims = dims_to_use)


# 5. Clustering to Identify Cell Groups
# Set default assay now that all six samples are integrated with batch correction
DefaultAssay(integrated_data20) <- "integrated"

# Now identifying clusters with the integrated data and then identify DEGs for each cluster
integrated_data20 <- FindNeighbors(integrated_data20, dims = dims_to_use)
integrated_data20 <- FindClusters(integrated_data20, resolution = 0.5)

# Save the integrated data object with multiple clusters identified across six samples
saveRDS(integrated_data20, "/home/outputs/totalNK_outputs/2_dge/dim15_20/integrated_data20_clusters.rds")



# Check the cumulative variance explained by PCs 1:10
pca_obj <- integrated_data20@reductions$pca
var_explained <- pca_obj@stdev^2 / sum(pca_obj@stdev^2)
cumulative_var <- cumsum(var_explained)
print(cumulative_var[1:20])

# Check a silhouette score to assess cluster separation
# Score > 0.5: Good separation
# Subsample to 10,000 cells
set.seed(123)
subsample_cells <- sample(colnames(integrated_data20), 10000)
subsampled_data <- subset(integrated_data20, cells = subsample_cells)
# Extract UMAP and clusters
umap_coords <- subsampled_data@reductions$umap@cell.embeddings
clusters <- subsampled_data$seurat_clusters
# Compute and print silhouette score
sil <- silhouette(as.numeric(clusters), dist(umap_coords))
mean_sil <- mean(sil[, 3])
print(paste("Average silhouette score (subsample):", mean_sil))

