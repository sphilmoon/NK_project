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

# output_path <- "home/outputs/totalNK_outputs/2_dge/dim30"
 
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
 
# Optional: Save the list for later use
# saveRDS(seurat_objects, file.path(output_path, "seurat_objects_list.rds"))
 
# 2. quality control and pre-processing
# Load the saved RDS file
# seurat_objects <- readRDS("seurat_objects_list.rds")
 
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

# Save the stacked plot as a single PNG
# output_file <- "qc_violin_stacked.png"
# ggsave(output_file, plot = stacked_plot, width = 10, height = 6 * length(seurat_objects), dpi = 600, units = "in")

# Print final message
# cat("✅ All QC plots combined into", output_file, "\n")

# Normalize with SCTransform
seurat_objects <- lapply(seurat_objects, SCTransform)
 
# Save the QC'd and normalized objects
# saveRDS(seurat_objects, "seurat_objects_qc_sct.rds")
 
# Final completion message
cat("QC and figure generation process is finished.\n")


# 3. Batch correction and integration
# Load the saved RDS file
# seurat_objects <- readRDS("seurat_objects_qc_sct.rds")

# Step 1: Run PCA on the first sample to determine dims
# Create a list to store all elbow plots
combined_elbow_plots <- list()

for (i in 1:length(seurat_objects)) {
  # Define animal name
  animal_name <- names(seurat_objects)[i]
  
  # Run PCA on the current sample
  seurat_objects[[i]] <- RunPCA(seurat_objects[[i]])
  
  # Generate and customize the elbow plot
  elbow_plot <- ElbowPlot(seurat_objects[[i]], ndims = 50) +
                ggtitle(paste("Elbow Plot for PCA -", animal_name)) +  # Dynamic title
                theme(plot.title = element_text(hjust = 0.5))  # Center the title
  
  # Store the plot in the list
  combined_elbow_plots[[animal_name]] <- elbow_plot
  
  # Print status message
  cat("✅ Elbow plot generated for", animal_name, "\n")
}

stacked_elbow_plot <- cowplot::plot_grid(plotlist = combined_elbow_plots, ncol = 1)

output_file <- "combined_elbow_plots.png"
ggsave(output_file, plot = stacked_elbow_plot, 
       width = 10, height = 6 * length(seurat_objects), dpi = 600, units = "in")

# Print final message
cat("✅ All elbow plots combined into", output_file, "\n")
# Save the plot
ggsave("elbow_plot.png", plot = elbow_plot, width = 10, height = 6, dpi = 600, units = "in")

# Optional: Display the plot
# print(elbow_plot)

# Step 2: Determine the number of dimensions for integration
dims_to_use <- 1:30  # Replace with your chosen dims

# Step 3: Select integration features
features <- SelectIntegrationFeatures(object.list = seurat_objects, nfeatures = 3000)

# Step 4: Find integration anchors
anchors <- FindIntegrationAnchors(object.list = seurat_objects, anchor.features = features, dims = dims_to_use)

# Step 5: Integrate data
integrated_data30 <- IntegrateData(anchorset = anchors, dims = dims_to_use)

# 4. Clustering and visualizing the integration
DefaultAssay(integrated_data30) <- "integrated"
integrated_data30 <- ScaleData(integrated_data30)
integrated_data30 <- RunPCA(integrated_data30)
integrated_data30 <- RunUMAP(integrated_data30, dims = dims_to_use)

# Save the integrated data object
saveRDS(integrated_data30, file.path("integrated_data30.rds"))

# Print final message
cat("✅ File saved for 4. Clustering and visualizing the integration.")

# Generate and customize the DimPlot
dim_plot <- DimPlot(integrated_data30, group.by = "sample", label = TRUE, pt.size = 0.5) +
            ggtitle("UMAP - Integration by Sample") +
            theme(plot.title = element_text(hjust = 0.5),
                  legend.position = "right") +  # Adjust legend position (or "none" to remove)
            scale_color_brewer(palette = "Set2")  # Optional: Change color palette

# Save the DimPlot as a PNG
ggsave("30_dimplot_integration.png", plot = dim_plot, width = 10, height = 8, dpi = 600, units = "in")


# 5. Clustering to Identify Cell Groups
# Load your integrated Seurat object
# integrated_data30 <- readRDS("integrated_data30.rds")

# Set default assay now that all six samples are integrated with batch correction
DefaultAssay(integrated_data30) <- "integrated"

# Now identifying clusters with the integrated data and then identify DEGs for each cluster
integrated_data30 <- FindNeighbors(integrated_data30, dims = dims_to_use)
integrated_data30 <- FindClusters(integrated_data30, resolution = 0.5)

# Visualize the clusters
cluster_plot <- DimPlot(integrated_data30, label = TRUE, group.by = "seurat_clusters")
ggsave("30_cluster_plot.png", plot = cluster_plot, width = 10, height = 8, dpi = 600, units = "in")

# Validate the integrated clusters with well-known Markers
feature_plot <- FeaturePlot(integrated_data30, features = c("NCR1", "CD3E"), ncol = 2)
ggsave("30_feature_plot_markers.png", plot = feature_plot, width = 12, height = 6, dpi = 600, units = "in")

# Save the integrated data object with multiple clusters identified across six samples
saveRDS(integrated_data30, "integrated_data30_clusters.rds")



# 6. DGE Analysis for finding commonly expressed genes per sample (not cluster) using PCA dim 30. 
