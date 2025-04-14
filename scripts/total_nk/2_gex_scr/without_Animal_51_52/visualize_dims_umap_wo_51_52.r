# Script: visualize_final_umap.R
# Purpose: Visualize UMAP plots for batch correction and clustering assessment
#          for dims20 (resolutions 0.4 and 0.3) and dims25 (resolutions 0.3 and 0.2).

# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)

# Define output directory
output_dir <- "/home/outputs/totalNK_outputs/2_dge/wo_51_52"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the dimensionalities to process
dims_list <- c(20, 25)
dim_names <- paste0("dims", dims_list)

# Define the resolutions for each dimensionality
resolution_configs <- list(
  "dims20" = c(0.4, 0.3),  # Resolutions for dims20
  "dims25" = c(0.3, 0.2)   # Resolutions for dims25
)

# Create a list to store integrated data objects
integrated_data_list <- list()

# Load the integrated Seurat objects for each dimensionality
for (dim_name in dim_names) {
  rds_file <- file.path(output_dir, paste0("integrated_data_wo_51_52_", dim_name, "_umap.rds"))
  if (file.exists(rds_file)) {
    integrated_data_list[[dim_name]] <- readRDS(rds_file)
    cat("✅ Loaded integrated data for", dim_name, "from", rds_file, "\n")
  } else {
    cat("⚠️ RDS file for", dim_name, "not found at", rds_file, "\n")
    next
  }
}

# Check if any data was loaded
if (length(integrated_data_list) == 0) {
  stop("No integrated data objects were loaded. Please check the RDS files.")
}

# Visualize UMAP for each dimensionality and resolution
for (dim_name in names(integrated_data_list)) {
  # Extract the integrated object
  integrated_data <- integrated_data_list[[dim_name]]
  
  # Extract the number of dimensions from dim_name (e.g., "dims20" -> 20)
  dim_value <- as.numeric(gsub("dims", "", dim_name))
  dims_to_use <- 1:dim_value
  
  # Get the resolutions for this dimensionality
  res_list <- resolution_configs[[dim_name]]
  
  # Process each resolution
  for (res in res_list) {
    cat("📌 Processing visualization for", dim_name, "at resolution", res, "\n")
    
    # Set the resolution for clustering
    integrated_data <- FindNeighbors(integrated_data, dims = dims_to_use)
    integrated_data <- FindClusters(integrated_data, resolution = res)
    
    # Generate UMAP plot by sample (to assess batch correction)
    umap_by_sample <- DimPlot(integrated_data,
                              group.by = "sample",
                              label = TRUE,
                              label.size = 5,
                              repel = TRUE,
                              pt.size = 0.3) +
                      ggtitle(paste("UMAP - Batch Correction (", dim_name, ", Res ", res, ")", sep = "")) +
                      theme(plot.title = element_text(hjust = 0.5),
                            legend.position = "right") +
                      scale_color_brewer(palette = "Set2")
    
    # Generate UMAP plot by cluster
    umap_by_cluster <- DimPlot(integrated_data,
                               group.by = "seurat_clusters",
                               label = TRUE,
                               label.size = 5,
                               repel = TRUE,
                               pt.size = 0.3) +
                       ggtitle(paste("UMAP - Clustering (", dim_name, ", Res ", res, ")", sep = "")) +
                       theme(plot.title = element_text(hjust = 0.5),
                             legend.position = "right") +
                       scale_color_brewer(palette = "Set3")
    
    # Combine the plots
    combined_umap_plot <- umap_by_sample + umap_by_cluster + plot_layout(ncol = 2)
    
    # Save the combined plot as a PDF with 600 DPI
    output_file <- file.path(output_dir, paste0("final_umap_", dim_name, "_res", res, ".pdf"))
    ggsave(output_file,
           plot = combined_umap_plot,
           width = 20,  # Wider to accommodate two plots
           height = 8,
           dpi = 600,
           units = "in")
    cat("✅ Combined UMAP plot for", dim_name, "at resolution", res, "saved to", output_file, "\n")
  }
}

# Final completion message
cat("✅ Visualization complete. All outputs saved in", output_dir, "\n")