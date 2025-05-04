# Script: check_pca_resolution_integrity.R
# Purpose: Check PCA and clustering integrity for different dims_to_use and resolutions
#          using seurat_objects_qc_sct_wo_51_52.rds.

# Load required libraries
library(Seurat)
library(cluster)  # For silhouette score
library(dplyr)

# Define output directory
output_dir <- "/home/outputs/totalNK_outputs/2_dge/wo_51_52/"  # Adjust to your path
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the dimensionalities to test
dims_list <- c(15, 20, 25)

# Define the resolutions to test
resolution_list <- c(0.1, 0.2, 0.3, 0.4, 0.5)

# Load the SCTransform-normalized Seurat objects
seurat_objects_file <- file.path(output_dir, "seurat_objects_qc_sct_wo_51_52.rds")
if (!file.exists(seurat_objects_file)) {
  stop("Seurat objects file not found at ", seurat_objects_file)
}
seurat_objects <- readRDS(seurat_objects_file)
cat("✅ Loaded seurat_objects from", seurat_objects_file, "\n")

# Verify the loaded objects
cat("Number of samples loaded:", length(seurat_objects), "\n")
lapply(names(seurat_objects), function(name) {
  cat("Sample:", name, " - Cells:", ncol(seurat_objects[[name]]), " Genes:", nrow(seurat_objects[[name]]), "\n")
})

# Select integration features
features <- SelectIntegrationFeatures(object.list = seurat_objects, nfeatures = 3000)

# Create a list to store integrated objects
integrated_data_list <- list()

# Perform integration for each dimensionality
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
  saveRDS(integrated_data, file.path(output_dir, paste0("integrated_data_wo_51_52_", dim_name, ".rds")))
  cat("✅ Integrated data for", dim_name, "saved to", file.path(output_dir, paste0("integrated_data_wo_51_52_", dim_name, ".rds")), "\n")
}

# Create a data frame to store integrity check results
integrity_results <- data.frame()

# Perform PCA, UMAP, clustering, and integrity checks
for (dim_name in names(integrated_data_list)) {
  cat("📌 Processing PCA and clustering for", dim_name, "\n")
  
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
  
  # Extract PCA object and compute variance explained
  pca_obj <- integrated_data@reductions$pca
  var_explained <- pca_obj@stdev^2 / sum(pca_obj@stdev^2)
  cumulative_var <- cumsum(var_explained)
  cat("Cumulative variance explained for", dim_name, "(first 30 PCs):\n")
  print(cumulative_var[1:30])
  
  # Store PCA variance explained for the dims used
  var_explained_dims <- sum(var_explained[1:dim_value])
  cat("Variance explained by", dim_name, "PCs:", var_explained_dims, "\n")
  
  # Test different resolutions
  for (res in resolution_list) {
    cat("📌 Testing resolution", res, "for", dim_name, "\n")
    
    # Find neighbors and clusters
    integrated_data <- FindNeighbors(integrated_data, dims = dims_to_use)
    integrated_data <- FindClusters(integrated_data, resolution = res)
    
    # Check silhouette score to assess cluster separation
    set.seed(123)
    total_cells <- ncol(integrated_data)
    subsample_size <- min(15000, total_cells)  # Use all cells if less than 15,000
    subsample_cells <- sample(colnames(integrated_data), subsample_size)
    subsampled_data <- subset(integrated_data, cells = subsample_cells)
    
    # Extract UMAP and clusters
    umap_coords <- subsampled_data@reductions$umap@cell.embeddings
    clusters <- subsampled_data$seurat_clusters
    
    # Compute and print silhouette score
    sil <- silhouette(as.numeric(clusters), dist(umap_coords))
    mean_sil <- mean(sil[, 3])
    cat("Average silhouette score for", dim_name, "at resolution", res, "(subsample):", mean_sil, "\n")
    if (mean_sil > 0.5) {
      cat("✅ Good cluster separation (score > 0.5)\n")
    } else {
      cat("⚠️ Poor cluster separation (score <= 0.5)\n")
    }
    
    # Store results
    integrity_results <- rbind(integrity_results, data.frame(
      Dims = dim_name,
      Resolution = res,
      VarianceExplained = var_explained_dims,
      SilhouetteScore = mean_sil,
      NumberOfClusters = length(unique(clusters))
    ))
  }
  
  # Save the UMAP result (optional, for future use)
  saveRDS(integrated_data, file.path(output_dir, paste0("integrated_data_wo_51_52_", dim_name, "_umap.rds")))
  cat("✅ UMAP result for", dim_name, "saved to", file.path(output_dir, paste0("integrated_data_wo_51_52_", dim_name, "_umap.rds")), "\n")
}

# Save integrity check results to CSV
write.csv(integrity_results, file.path(output_dir, "new_integrity_check_results.csv"), row.names = FALSE)
cat("✅ Integrity check results saved to", file.path(output_dir, "new_integrity_check_results.csv"), "\n")

# Final completion message
cat("✅ Analysis complete. All outputs saved in", output_dir, "\n")