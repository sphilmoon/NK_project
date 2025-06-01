# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(ggplot2)
library(cowplot)

# ------------------------- #
# Define Paths and Configs
# ------------------------- #
output_dir <- "/home/outputs/all_merged_TotalNK_nkp46"

# Create output subdirectories
rds_dir <- file.path(output_dir, "rds")
figures_dir <- file.path(output_dir, "figures")
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------- #
# Load Merged Seurat Object
# ------------------------- #
merged_obj <- readRDS(file.path(rds_dir, "all_totalNK_nkp46_dim25_res0.5_clustered.rds"))
cat("✅ Loaded merged Seurat object\n")

# view all metadata column names
print(colnames(merged_obj@meta.data))

# check unique values in sample_id
print(unique(merged_obj$sample_id))

# view frequency of each sample_id
print(table(merged_obj$sample_id))

# inspect the first few rows of metadata
print(head(merged_obj@meta.data[, c("sample_id", "sample", "animal_id")]))

cat("🔍 Inspecting metadata columns: sample, sample_id, animal, animal_id\n\n")

# Unique values
cat("🔍 Unique values in 'sample':\n", paste(unique(merged_obj$sample), collapse = ", "), "\n")
cat("🔍 Unique values in 'sample_id':\n", paste(unique(merged_obj$sample_id), collapse = ", "), "\n")
cat("🔍 Unique values in 'animal':\n", paste(unique(merged_obj$animal), collapse = ", "), "\n")
cat("🔍 Unique values in 'animal_id':\n", paste(unique(merged_obj$animal_id), collapse = ", "), "\n\n")
