# Script: visualize_umap_clusters.R
# Purpose: Visualize UMAP plots by sample and clusters for each dims_to_use, 
#          label clusters with cell counts, and save cluster info as CSV.

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(ggrepel)

# Define output directory (replace with your actual path)
output_dir <- "/home/outputs/totalNK_outputs/2_dge/wo_51_52"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the dimensionalities to process
dims_list <- c(10, 15, 25)
dim_names <- paste0("dims", dims_list)

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

# Visualize UMAP and clusters for each dimensionality
for (dim_name in names(integrated_data_list)) {
  cat("📌 Processing visualization for", dim_name, "\n")
  
  # Extract the integrated object
  integrated_data <- integrated_data_list[[dim_name]]
  
  # Generate the UMAP plot (by sample)
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
  
  # Calculate the number of cells per cluster
  cluster_counts <- as.data.frame(table(integrated_data$seurat_clusters))
  colnames(cluster_counts) <- c("Cluster", "CellCount")
  
  # Calculate the average number of expressed genes per cluster
  DefaultAssay(integrated_data) <- "RNA"
  print(names(integrated_data[["RNA"]]$layers))

  DefaultAssay(integrated_data) <- "RNA"
  str(integrated_data[["RNA"]])
  
  count_data <- GetAssayData(integrated_data, layer = "counts")
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
    Cluster = levels(integrated_data$seurat_clusters),
    CellCount = cluster_counts$CellCount,
    AvgGenesPerCell = gene_counts
  )
  
  # Save cluster info to CSV
  write.csv(cluster_info, file.path(output_dir, paste0("cluster_info_wo_51_52", dim_name, ".csv")), row.names = FALSE)
  cat("✅ Cluster info for", dim_name, "saved to", file.path(output_dir, paste0("cluster_info_wo_51_52", dim_name, ".csv")), "\n")
  
  # Get UMAP coordinates and compute cluster centroids
  umap_data <- as.data.frame(integrated_data@reductions$umap@cell.embeddings)
  umap_data$Cluster <- integrated_data$seurat_clusters
  centroids <- umap_data %>%
    group_by(Cluster) %>%
    summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2)) %>%
    ungroup()
  
  # Merge centroids with cluster counts for labeling
  centroids <- merge(centroids, cluster_counts, by = "Cluster")
  centroids$Label <- paste0("Cluster ", centroids$Cluster, " (", centroids$CellCount, " cells)")
  
  # Generate the cluster plot with custom labels
  cluster_plot <- DimPlot(integrated_data,
                          group.by = "seurat_clusters",
                          label = FALSE,  # Disable default labels
                          pt.size = 0.3) +
                  ggtitle(paste("UMAP - Clustering (", dim_name, ")", sep = "")) +
                  theme(plot.title = element_text(hjust = 0.5),
                        legend.position = "right") +
                  scale_color_brewer(palette = "Set3") +
                  ggrepel::geom_text_repel(data = centroids,
                                            aes(x = UMAP_1, y = UMAP_2, label = Label),
                                            size = 5,
                                            box.padding = 0.5,
                                            point.padding = 0.5,
                                            max.overlaps = Inf)
  
  # Combine UMAP and cluster plots
  combined_umap_plot <- umap_plot + cluster_plot + plot_layout(ncol = 2)
  
  # Save the combined plot as a PDF
  output_file <- file.path(output_dir, paste0("new_combined_umap_cluster_wo_51_52", dim_name, ".pdf"))
  ggsave(output_file,
         plot = combined_umap_plot,
         width = 20,  # Wider to accommodate two plots
         height = 8,
         dpi = 600,
         units = "in")
  cat("✅ Combined UMAP and cluster plot for", dim_name, "saved to", output_file, "\n")
}

# Final completion message
cat("✅ Visualization complete. All outputs saved in", output_dir, "\n")