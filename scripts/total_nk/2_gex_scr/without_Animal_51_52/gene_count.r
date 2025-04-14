# Script: visualize_final_umap_with_genes.R
# Purpose:
# - Visualize UMAP plots for batch correction and clustering assessment
# - Count number of cells expressing each gene per cluster
# - Save gene-cluster expression matrices and cluster cell counts to CSV

# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# ------------------------- #
# Function: Gene Expression Count per Cluster
# ------------------------- #
compute_gene_cluster_matrix <- function(seurat_obj, assay_priority = "SCT", threshold = 0, output_path = NULL) {
  assay_to_use <- assay_priority[assay_priority %in% Assays(seurat_obj)][1]
  if (is.null(assay_to_use)) stop("❌ None of the specified assays are available in the Seurat object.")

  DefaultAssay(seurat_obj) <- assay_to_use

  available_layers <- Layers(seurat_obj[[assay_to_use]])
  layer_to_use <- NULL

  # Choose the best available layer
  for (candidate in c("counts", "data", "logcounts", "corrected")) {
    if (candidate %in% available_layers) {
      expr_data <- GetAssayData(seurat_obj, layer = candidate, assay = assay_to_use)
      if (sum(expr_data) > 0) {
        layer_to_use <- candidate
        break
      }
    }
  }

  if (is.null(layer_to_use)) {
    stop(paste("❌ No usable expression layer found in assay:", assay_to_use))
  } else {
    cat("✅ Using layer:", layer_to_use, "from assay:", assay_to_use, "\n")
  }

  # Expression matrix now available
  clusters <- Idents(seurat_obj)
  cluster_levels <- levels(clusters)

  binary_expr <- expr_data > threshold

  gene_cluster_mat <- matrix(0, nrow = nrow(binary_expr), ncol = length(cluster_levels))
  rownames(gene_cluster_mat) <- rownames(binary_expr)
  colnames(gene_cluster_mat) <- cluster_levels

  for (clust in cluster_levels) {
    cells_in_clust <- WhichCells(seurat_obj, idents = clust)
    if (!all(cells_in_clust %in% colnames(binary_expr))) {
      stop("❌ Some cluster cells are missing in expression matrix. Possible mismatch in cell names.")
    }
    gene_cluster_mat[, clust] <- Matrix::rowSums(binary_expr[, cells_in_clust, drop = FALSE])
  }

  gene_cluster_df <- as.data.frame(gene_cluster_mat)
  gene_cluster_df$Gene <- rownames(gene_cluster_df)
  gene_cluster_df <- gene_cluster_df[, c("Gene", cluster_levels)]

  if (!is.null(output_path)) {
    write.csv(gene_cluster_df, output_path, row.names = FALSE)
    cat("✅ Gene-by-cluster expression count saved to", output_path, "\n")
  }

  return(gene_cluster_df)
}

# ------------------------- #
# Define Paths and Configs
# ------------------------- #

output_dir <- "/home/outputs/totalNK_outputs/2_dge/wo_51_52"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

dims_list <- c(20, 25)
dim_names <- paste0("dims", dims_list)

resolution_configs <- list(
  "dims20" = 0.4,
  "dims25" = 0.3
)

# Load Seurat objects
integrated_data_list <- list()
for (dim_name in dim_names) {
  rds_file <- file.path(output_dir, paste0("integrated_data_wo_51_52_", dim_name, "_umap.rds"))
  if (file.exists(rds_file)) {
    integrated_data_list[[dim_name]] <- readRDS(rds_file)
    cat("✅ Loaded integrated data for", dim_name, "from", rds_file, "\n")
  } else {
    cat("⚠️ RDS file for", dim_name, "not found at", rds_file, "\n")
  }
}
if (length(integrated_data_list) == 0) {
  stop("❌ No integrated data objects were loaded.")
}

# ------------------------- #
# Main Loop: Process each combination
# ------------------------- #

for (dim_name in names(integrated_data_list)) {
  integrated_data <- integrated_data_list[[dim_name]]
  dim_value <- as.numeric(gsub("dims", "", dim_name))
  dims_to_use <- 1:dim_value
  res_list <- resolution_configs[[dim_name]]
  if (length(res_list) == 1) res_list <- as.vector(res_list)

  for (res in res_list) {
    cat("📌 Processing", dim_name, "at resolution", res, "\n")

    # Clustering
    integrated_data <- FindNeighbors(integrated_data, dims = dims_to_use)
    integrated_data <- FindClusters(integrated_data, resolution = res)

    # checking assays and layers
    Assays(integrated_data)
    # Layers(integrated_data[["RNA"]])  # Or replace "RNA" with "SCT" if using SCTransform
    Layers(integrated_data[["SCT"]])

    # UMAP Plots
    umap_by_sample <- DimPlot(integrated_data,
                              group.by = "sample",
                              label = TRUE,
                              label.size = 5,
                              repel = TRUE,
                              pt.size = 0.3) +
                      ggtitle(paste("UMAP - Batch Correction (", dim_name, ", Res ", res, ")", sep = "")) +
                      theme(plot.title = element_text(hjust = 0.5),
                            legend.position = "right") + scale_color_brewer(palette = "Set2")

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
                       scale_color_viridis_d(option = "viridis", direction = -1)

    combined_umap_plot <- umap_by_sample + umap_by_cluster + plot_layout(ncol = 2)

    # Save UMAP PDF
    umap_file <- file.path(output_dir, paste0("final_umap_", dim_name, "_res", res, ".pdf"))
    ggsave(umap_file, plot = combined_umap_plot,
           width = 20, height = 8, dpi = 600, units = "in")
    cat("✅ UMAP plot saved to", umap_file, "\n")

    # Gene-by-cluster expression matrix
    # gene_cluster_csv <- file.path(output_dir, paste0("gene_counts_", dim_name, "_res", res, ".csv"))
    # gene_cluster_counts <- compute_gene_cluster_matrix(
    #   seurat_obj = integrated_data,
    #   assay_priority = "SCT",
    #   threshold = 0,
    #   output_path = gene_cluster_csv
    # )

    # # Cluster cell counts
    # cluster_counts <- as.data.frame(table(Idents(integrated_data)))
    # colnames(cluster_counts) <- c("Cluster", "CellCount")
    # cluster_csv <- file.path(output_dir, paste0("cluster_cell_counts_", dim_name, "_res", res, ".csv"))
    # write.csv(cluster_counts, cluster_csv, row.names = FALSE)
    # cat("✅ Cluster cell counts saved to", cluster_csv, "\n")

    # cat("✅ Completed processing for", dim_name, "at resolution", res, "\n\n")
  }
}

cat("🎉 All visualizations and analyses complete. Outputs saved in:\n", output_dir, "\n")