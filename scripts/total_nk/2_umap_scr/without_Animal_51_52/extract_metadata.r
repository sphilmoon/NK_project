# Script: extract_metadata_with_genes.R
# Purpose: Extract specified metadata columns and gene list from a Seurat object and save as CSV

# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)

# ------------------------- #
# Define Paths and Configs
# ------------------------- #
rds_file <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/rds/integrated_data_dims25_res0.3_genecounts.rds"
output_dir <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/csv"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------- #
# Load Seurat Object
# ------------------------- #
seurat_obj <- readRDS(rds_file)
cat("✅ Loaded Seurat object from", rds_file, "\n")

# ------------------------- #
# Set Assay
# ------------------------- #
DefaultAssay(seurat_obj) <- if ("SCT" %in% Assays(seurat_obj)) "SCT" else "RNA"
cat("✅ Default assay set to", DefaultAssay(seurat_obj), "\n")

# ------------------------- #
# Extract Specified Metadata Columns
# ------------------------- #
metadata_columns <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "sample", 
                      "nCount_SCT", "nFeature_SCT", "integrated_snn_res.0.1", 
                      "seurat_clusters", "integrated_snn_res.0.2", "integrated_snn_res.0.3", 
                      "integrated_snn_res.0.4", "integrated_snn_res.0.5")

# Check if all columns exist
missing_columns <- metadata_columns[!metadata_columns %in% colnames(seurat_obj@meta.data)]
if (length(missing_columns) > 0) {
  warning("⚠️ Missing metadata columns: ", paste(missing_columns, collapse = ", "), 
          ". These will be excluded.")
  metadata_columns <- metadata_columns[metadata_columns %in% colnames(seurat_obj@meta.data)]
}

# Extract the metadata
metadata_df <- seurat_obj@meta.data[, metadata_columns, drop = FALSE]

# ------------------------- #
# Extract Gene Names from Assay Data
# ------------------------- #
expr_data_all <- GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = "data")
gene_names <- rownames(expr_data_all)
cat("🔍 Total number of genes found:", length(gene_names), "\n")
# cat("🔍 First 10 genes (for reference):", paste(head(gene_names, 10), collapse = ", "), "\n")

# Create a single comma-separated string of all gene names
gene_names_string <- paste(gene_names, collapse = ",")

# Add a new column to the metadata data frame with the gene list
metadata_df$all_genes <- gene_names_string

# ------------------------- #
# Save to CSV
# ------------------------- #
output_file <- file.path(output_dir, "metadata_with_genes_dims25_res0.3.csv")
write.csv(metadata_df, file = output_file, row.names = TRUE, quote = FALSE)
cat("✅ Metadata with gene list saved to", output_file, "\n")

# ------------------------- #
# Final Completion Message
# ------------------------- #
cat("✅ Metadata and gene extraction complete. Output saved in", output_dir, "\n")
