# Script: visualize_final_umap_with_genes.R
# Purpose: Visualize UMAP plots for batch correction and clustering assessment,
#          calculate the average number of genes expressed per cluster,
#          and save results to CSV for dims10 (resolution 0.5), dims20 (resolutions 0.4 and 0.3),
#          and dims25 (resolutions 0.3 and 0.2).

# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Define output directory
output_dir <- "/home/outputs/totalNK_outputs/2_dge/wo_51_52"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the dimensionalities to process
dims_list <- c(10, 20, 25)
dim_names <- paste0("dims", dims_list)

# Define the resolutions for each dimensionality
resolution_configs <- list(
  "dims10" = 0.5,           # Resolution for dims10
  "dims20" = c(0.4, 0.3),   # Resolutions for dims20
  "dims25" = c(0.3, 0.2)    # Resolutions for dims25
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

# Initialize a data frame to store cluster information across all combinations
all_cluster_info <- data.frame()

# Visualize UMAP and calculate genes per cluster for each dimensionality and resolution
for (dim_name in names(integrated_data_list)) {
  # Extract the integrated object
  integrated_data <- integrated_data_list[[dim_name]]
  
  # Extract the number of dimensions from dim_name (e.g., "dims20" -> 20)
  dim_value <- as.numeric(gsub("dims", "", dim_name))
  dims_to_use <- 1:dim_value
  
  # Get the resolutions for this dimensionality
  res_list <- resolution_configs[[dim_name]]
  
  # Ensure res_list is a vector (for dims10, it's a single value)
  if (length(res_list) == 1) {
    res_list <- as.vector(res_list)
  }
  
  # Process each resolution
  for (res in res_list) {
    cat("📌 Processing visualization and gene calculation for", dim_name, "at resolution", res, "\n")
    
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
    
    # Calculate the number of cells per cluster
    cluster_counts <- as.data.frame(table(integrated_data$seurat_clusters))
    colnames(cluster_counts) <- c("Cluster", "CellCount")
    
    # Calculate the average number of expressed genes per cluster
    DefaultAssay(integrated_data) <- "RNA"
    
    # Check available layers in the RNA assay
    cat("Available layers in RNA assay for", dim_name, ":\n")
    available_layers <- names(integrated_data[["RNA"]]$layers)
    print(available_layers)
    
    # Use the 'counts' layer if available, otherwise fall back to 'data'
    if (!"counts" %in% available_layers) {
      cat("⚠️ 'counts' layer not found. Checking for 'data' layer...\n")
      if ("data" %in% available_layers) {
        cat("✅ 'data' layer found. Using 'data' layer instead.\n")
        count_data <- GetAssayData(integrated_data, layer = "data")
      } else {
        stop("No 'counts' or 'data' layer found in RNA assay. Cannot proceed.")
      }
    } else {
      cat("✅ 'counts' layer found. Using 'counts' layer.\n")
      count_data <- GetAssayData(integrated_data, layer = "counts")
    }
    
    gene_counts <- sapply(levels(integrated_data$seurat_clusters), function(cluster) {
      cells_in_cluster <- WhichCells(integrated_data, idents = cluster)
      cluster_counts <- count_data[, cells_in_cluster, drop = FALSE]
      # Number of genes with non-zero counts per cell
      genes_per_cell <- colSums(cluster_counts > 0)
      # Average across cells in the cluster
      mean(genes_per_cell)
    })
    
    # Create a data frame with cluster info
    cluster_info <- data.frame(
      Dims = dim_name,
      Resolution = res,
      Cluster = levels(integrated_data$seurat_clusters),
      CellCount = cluster_counts$CellCount,
      AvgGenesPerCell = gene_counts
    )
    
    # Append to the overall results
    all_cluster_info <- rbind(all_cluster_info, cluster_info)
  }
}

# Save cluster info to CSV
write.csv(all_cluster_info, file.path(output_dir, "cluster_info_all_combinations.csv"), row.names = FALSE)
cat("✅ Cluster info saved to", file.path(output_dir, "cluster_info_all_combinations.csv"), "\n")

# Final completion message
cat("✅ Visualization and analysis complete. All outputs saved in", output_dir, "\n")