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

# Define output directory and create it if it doesn't exist
output_dir <- "/home/outputs/totalNK_outputs/2_dge/wo_51_52"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define file paths and sample names (excluding Animal51 and Animal52)
file_paths <- list(
  "/home/rawdata/total_nk/animal25_totalnk_filtered_feature_bc_matrix.h5",
  "/home/rawdata/total_nk/animal26_totalnk_filtered_feature_bc_matrix.h5",
  "/home/rawdata/total_nk/animal27_totalnk_filtered_feature_bc_matrix.h5",
  "/home/rawdata/total_nk/animal28_totalnk_filtered_feature_bc_matrix.h5"
)
sample_names <- c("Animal25", "Animal26", "Animal27", "Animal28")

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

# Save the list for later use
saveRDS(seurat_objects, file.path(output_dir, "seurat_objects_list_wo_51_52.rds"))
cat("✅ Seurat objects list saved to", file.path(output_dir, "seurat_objects_list_wo_51_52.rds"), "\n")

# 2. Quality Control and Preprocessing
# Verify the loaded objects
cat("Number of samples loaded:", length(seurat_objects), "\n")
lapply(names(seurat_objects), function(name) {
  cat("Sample:", name, " - Cells:", ncol(seurat_objects[[name]]), " Genes:", nrow(seurat_objects[[name]]), "\n")
})

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
  } else {
    # Default threshold for Animal25 and Animal26
    seurat_objects[[i]] <- subset(seurat_objects[[i]], 
                                  subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & 
                                           nCount_RNA > 500 & nCount_RNA < 10000)
  }
  
  # Create separate violin plots for each feature, removing the legend
  p1 <- VlnPlot(seurat_objects[[i]], features = "nFeature_RNA", ncol = 1) +
        ggtitle("Gene counts QC") +
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = "none")
  
  p2 <- VlnPlot(seurat_objects[[i]], features = "nCount_RNA", ncol = 1) +
        ggtitle("UMI counts QC") +
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = "none")
  
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

# Save the stacked plot as a single pdf file
output_file <- file.path(output_dir, "qc_violin_stacked_wo_51_52.pdf")
ggsave(output_file, plot = stacked_plot, width = 10, height = 6 * length(seurat_objects), dpi = 600, units = "in")
cat("✅ All QC plots combined into", output_file, "\n")

# Normalize with SCTransform
seurat_objects <- lapply(seurat_objects, SCTransform, verbose = FALSE)

# Save the QC'd and normalized objects
saveRDS(seurat_objects, file.path(output_dir, "seurat_objects_qc_sct_wo_51_52.rds"))
cat("✅ SCTransform normalized objects saved to", file.path(output_dir, "seurat_objects_qc_sct_wo_51_52.rds"), "\n")

# 3. Batch Correction and Integration
# Select integration features
features <- SelectIntegrationFeatures(object.list = seurat_objects, nfeatures = 3000)

# Define the dimensions to test
dims_list <- c(10, 15, 25)

# Create a named list to store integrated objects
integrated_data_list <- list()

# Loop over each dimensionality
for (dim in dims_list) {
  dim_name <- paste0("dims", dim)
  cat("📌 Processing integration for", dim_name, "\n")
  
  # Set dims_to_use for this iteration
  dims_to_use <- 1:dim
  
  # Find integration anchors
  anchors <- FindIntegrationAnchors(object.list = seurat_objects, anchor.features = features, dims = dims_to_use)
  
  # Integrate data
  integrated_data <- IntegrateData(anchorset = anchors, dims = dims_to_use)
  
  # Store in the list
  integrated_data_list[[dim_name]] <- integrated_data
  
  # Save the integrated object
#   saveRDS(integrated_data, file.path(output_dir, paste0("integrated_data_wo_51_52_", dim_name, ".rds")))
#   cat("✅ Integrated data for", dim_name, "saved to", file.path(output_dir, paste0("integrated_data_wo_51_52_", dim_name, ".rds")), "\n")
}

# 4. UMAP and Visualization
for (dim_name in names(integrated_data_list)) {
  cat("📌 Processing UMAP and clustering for", dim_name, "\n")
  
  # Extract the integrated object
  integrated_data <- integrated_data_list[[dim_name]]
  
  # Extract the number of dimensions from dim_name (e.g., "dims10" -> 10)
  dim_value <- as.numeric(gsub("dims", "", dim_name))
  dims_to_use <- 1:dim_value
  
  # Set default assay
  DefaultAssay(integrated_data) <- "integrated"
  
  # Scale data, run PCA, and UMAP
  integrated_data <- ScaleData(integrated_data)
  integrated_data <- RunPCA(integrated_data, verbose = FALSE)
  integrated_data <- RunUMAP(integrated_data, dims = dims_to_use)
  
  # Find neighbors and clusters
  integrated_data <- FindNeighbors(integrated_data, dims = dims_to_use)
  integrated_data <- FindClusters(integrated_data, resolution = 0.5)
  
  # Save the UMAP result
  saveRDS(integrated_data, file.path(output_dir, paste0("integrated_data_wo_51_52_", dim_name, "_umap.rds")))
  cat("✅ UMAP result for", dim_name, "saved to", file.path(output_dir, paste0("integrated_data_wo_51_52_", dim_name, "_umap.rds")), "\n")
  
  # Generate and customize the UMAP plot
  umap_plot <- DimPlot(integrated_data,
                       group.by = "sample",
                       label = TRUE,
                       label.size = 5,
                       repel = TRUE,
                       pt.size = 0.3) +
               ggtitle(paste("UMAP - Integration by Sample (", dim_name, ")", sep = "")) +
               theme(plot.title = element_text(hjust = 0.5),
                     legend.position = "right") +
               scale_color_brewer(palette = "Set2")
  
  # Save the UMAP plot as a PDF
  output_file <- file.path(output_dir, paste0("umap_", dim_name, ".pdf"))
  ggsave(output_file,
         plot = umap_plot,
         width = 10,
         height = 8,
         dpi = 600,
         units = "in")
  cat("✅ UMAP plot for", dim_name, "saved to", output_file, "\n")
  
  # 5. Integrity Checks
  # Extract PCA object and compute variance explained
  pca_obj <- integrated_data@reductions$pca
  var_explained <- pca_obj@stdev^2 / sum(pca_obj@stdev^2)
  cumulative_var <- cumsum(var_explained)
  cat("Cumulative variance explained for", dim_name, "(first 30 PCs):\n")
  print(cumulative_var[1:30])
  
  # Check silhouette score to assess cluster separation
  # Subsample to 10,000 cells
  set.seed(123)
  total_cells <- ncol(integrated_data)
  subsample_size <- min(10000, total_cells)  # Use all cells if less than 10,000
  subsample_cells <- sample(colnames(integrated_data), subsample_size)
  subsampled_data <- subset(integrated_data, cells = subsample_cells)
  # Extract UMAP and clusters
  umap_coords <- subsampled_data@reductions$umap@cell.embeddings
  clusters <- subsampled_data$seurat_clusters
  # Compute and print silhouette score
  sil <- silhouette(as.numeric(clusters), dist(umap_coords))
  mean_sil <- mean(sil[, 3])
  cat("Average silhouette score for", dim_name, "(subsample):", mean_sil, "\n")
  if (mean_sil > 0.5) {
    cat("✅ Good cluster separation (score > 0.5)\n")
  } else {
    cat("⚠️ Poor cluster separation (score <= 0.5)\n")
  }
}

# Final completion message
cat("✅ Analysis complete. All outputs saved in", output_dir, "\n")