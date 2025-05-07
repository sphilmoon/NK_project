# Script: extract_cluster_barcodes.R
# Purpose: Extract cell barcodes assigned to clusters for dims25 and resolution 0.3

# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(dplyr)

# ------------------------- #
# Define Output Directory
# ------------------------- #
output_dir <- "/home/outputs/totalNK_outputs/2_dge/wo_51_52/csv"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------- #
# Define Dimensionalities and Resolutions
# ------------------------- #
rds_file <- file.path(output_dir, paste0("integrated_data_dims25_res0.3_genecounts.rds")



# Check if any data was loaded
if (length(integrated_data_list) == 0) {
  stop("⚠️ No integrated data objects were loaded. Please check the RDS files.")
}

# ------------------------- #
# Extract Barcodes for dims25 and resolution 0.3
# ------------------------- #
dim_name <- "dims25"  # Target dimensionality
if (dim_name %in% names(integrated_data_list)) {
  cat("📌 Extracting barcodes for", dim_name, "with resolution", resolutions[[dim_name]], "\n")
  
  # Extract the integrated object
  integrated_data <- integrated_data_list[[dim_name]]
  
  # Extract the number of dimensions from dim_name (e.g., "dims25" -> 25)
  dim_value <- as.numeric(gsub("dims", "", dim_name))
  dims_to_use <- 1:dim_value
  
  # Ensure clustering is performed (already done in your script, but re-run for consistency)
  integrated_data <- FindNeighbors(integrated_data, dims = dims_to_use)
  integrated_data <- FindClusters(integrated_data, resolution = resolutions[[dim_name]])
  
  # Extract cell barcodes and cluster assignments
  barcodes <- rownames(integrated_data@meta.data)
  cluster_ids <- integrated_data@meta.data[["seurat_clusters"]]
  
  # Create a data frame with barcodes and cluster IDs
  barcode_cluster_df <- data.frame(
    barcode = barcodes,
    cluster_id = as.character(cluster_ids),
    stringsAsFactors = FALSE
  )
  
  # ------------------------- #
  # Save to CSV
  # ------------------------- #
  output_file <- file.path(output_dir, paste0("cluster_barcodes_", dim_name, "_res", resolutions[[dim_name]], ".csv"))
  write.csv(barcode_cluster_df, file = output_file, row.names = FALSE, quote = FALSE)
  cat("✅ Cluster barcodes saved to", output_file, "\n")
} else {
  stop("❌ Integrated data for", dim_name, "not found in loaded objects.")
}

# ------------------------- #
# Final Completion Message
# ------------------------- #
cat("✅ Barcode extraction complete. Output saved in", output_dir, "\n")