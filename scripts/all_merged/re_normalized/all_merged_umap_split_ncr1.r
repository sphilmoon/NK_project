# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(ggplot2)
library(cowplot)

# ------------------------- #
# Define Paths and Configs
# ------------------------- #
totalNK_rds <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/rds/integrated_data_dims25_res0.3_genecounts.rds"
nkp46_rds <- "/home/outputs/nkp46_outputs/nkp46_integrated_data.rds"
output_dir <- "/home/outputs/all_merged_TotalNK_nkp46"

# Create output subdirectories
rds_dir <- file.path(output_dir, "rds")
figures_dir <- file.path(output_dir, "figures")
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------- #
# Load Merged Seurat Object
# ------------------------- #
merged_obj <- readRDS(file.path(rds_dir, "all_totalNK_nkp46.rds"))
cat("✅ Loaded merged Seurat object\n")

# ------------------------- #
# Run PCA, Clustering, and UMAP
# ------------------------- #
merged_obj <- RunPCA(merged_obj, npcs = 25)
cat("✅ Ran PCA with 25 dimensions\n")

merged_obj <- FindNeighbors(merged_obj, dims = 1:25)
merged_obj <- FindClusters(merged_obj, resolution = 0.5)
merged_obj <- RunUMAP(merged_obj, dims = 1:25)
cat("✅ Completed clustering and UMAP with resolution 0.5\n")

# save the merged object
saveRDS(merged_obj, file.path(rds_dir, "all_totalNK_nkp46_dim25_res0.5_clustered.rds"))
cat("✅ Saved merged Seurat object\n")

# read the merged object
merged_obj <- readRDS(file.path(rds_dir, "all_totalNK_nkp46_dim25_res0.5_clustered.rds"))
cat("✅ Loaded merged Seurat object\n")
