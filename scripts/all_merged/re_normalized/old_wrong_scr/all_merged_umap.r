# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

# ------------------------- #
# Define Paths and Configs
# ------------------------- #
totalNK_rds <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/rds/integrated_data_dims25_res0.3_genecounts.rds"
nkp46_rds <- "/home/outputs/nkp46_outputs/nkp46_integrated_data.rds"
output_dir <- "/home/outputs/all_merged_TotalNK_nkp46"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------- #
# Load Seurat Objects
# ------------------------- #
totalNK <- readRDS(totalNK_rds)
cat("✅ Loaded Total NK Seurat object from", totalNK_rds, "\n")

nkp46_data <- readRDS(nkp46_rds)
cat("✅ Loaded NKp46 Seurat object from", nkp46_rds, "\n")

# ------------------------- #
# Set Consistent Assay
# ------------------------- #
# Ensure both objects use the RNA assay
DefaultAssay(totalNK) <- "RNA"
DefaultAssay(nkp46_data) <- "RNA"
cat("✅ Set default assay to RNA for both objects\n")

# Join layers if using Seurat v5
totalNK <- JoinLayers(totalNK, assay = "RNA")
nkp46_data <- JoinLayers(nkp46_data, assay = "RNA")
cat("✅ Joined layers for both objects\n")

# ------------------------- #
# Validate Counts Slot Before Subsetting
# ------------------------- #
validate_counts <- function(obj, name) {
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  if (is.null(counts)) {
    stop(sprintf("❌ No counts slot found in %s", name))
  }
  
  # Check for non-numeric values using sparse matrix properties
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(counts, "dgCMatrix")
  }
  
  # Extract values and check for non-numeric (use sparse operations to avoid dense conversion)
  values <- summary(counts)$x  # Extract non-zero values
  if (!all(is.numeric(values))) {
    # Convert to dense temporarily for detailed inspection (small subset for logging)
    dense_counts <- as.matrix(counts[1:min(5, nrow(counts)), 1:min(5, ncol(counts))])
    cat(sprintf("⚠️ First 5x5 subset of %s counts slot:\n", name))
    print(dense_counts)
    stop(sprintf("❌ Non-numeric values found in counts slot of %s", name))
  }
  cat(sprintf("✅ %s counts slot contains only numeric values\n", name))
}

# Validate each object
validate_counts(totalNK, "totalNK")
validate_counts(nkp46_data, "nkp46_data")

# ------------------------- #
# Subset nkp46+ and nkp46- populations
# ------------------------- #
nkp46_pos <- subset(nkp46_data, subset = condition == "nkp46+")
nkp46_neg <- subset(nkp46_data, subset = condition == "nkp46-")
cat("✅ Subsetted NKp46+ and NKp46- populations\n")

# Validate subsets
validate_counts(nkp46_pos, "nkp46_pos")
validate_counts(nkp46_neg, "nkp46_neg")

# ------------------------- #
# Annotate Sample Identity
# ------------------------- #
totalNK$sample_id <- "TotalNK"
nkp46_pos$sample_id <- "NKp46pos"
nkp46_neg$sample_id <- "NKp46neg"
cat("✅ Annotated sample identities\n")

# ------------------------- #
# Merge Seurat Objects
# ------------------------- #
merged_obj <- merge(totalNK, y = c(nkp46_pos, nkp46_neg), 
                    add.cell.ids = c("TotalNK", "NKp46pos", "NKp46neg"), 
                    project = "NK_Combined")
cat("✅ Merged Seurat objects\n")

# Ensure RNA assay for merged object
DefaultAssay(merged_obj) <- "RNA"
merged_obj <- JoinLayers(merged_obj, assay = "RNA")
cat("✅ Set default assay to RNA and joined layers for merged object\n")

# Validate merged counts
validate_counts(merged_obj, "merged_obj")

# ------------------------- #
# Normalize and Scale Data
# ------------------------- #
merged_obj <- NormalizeData(merged_obj)
cat("✅ Normalized data\n")

# Calculate gene variances manually to debug
cat("🔍 Checking gene variances...\n")
data <- GetAssayData(merged_obj, assay = "RNA", layer = "data")
gene_vars <- apply(data, 1, var, na.rm = TRUE)

# Check for zero variance genes
if (all(gene_vars == 0)) {
  stop("❌ All genes have zero variance after normalization")
}
cat("✅ Found", sum(gene_vars > 0), "genes with non-zero variance\n")

# Filter out zero-variance genes
non_zero_var_genes <- names(gene_vars[gene_vars > 0])
VariableFeatures(merged_obj) <- non_zero_var_genes
cat("✅ Set", length(non_zero_var_genes), "variable features with non-zero variance\n")

merged_obj <- FindVariableFeatures(merged_obj, selection.method = "vst", nfeatures = 2000)
cat("✅ Found variable features\n")

merged_obj <- ScaleData(merged_obj)
cat("✅ Scaled data\n")

# ------------------------- #
# Run PCA with Consistent Dims
# ------------------------- #
merged_obj <- RunPCA(merged_obj, npcs = 25)
cat("✅ Completed PCA with 25 dimensions\n")

# ------------------------- #
# Cluster and UMAP with Consistent Resolution
# ------------------------- #
merged_obj <- FindNeighbors(merged_obj, dims = 1:25)
merged_obj <- FindClusters(merged_obj, resolution = 0.3)
merged_obj <- RunUMAP(merged_obj, dims = 1:25)
cat("✅ Completed clustering and UMAP with resolution 0.3\n")

# ------------------------- #
# Create Combined UMAP Visualization
# ------------------------- #
cat("🎨 Creating combined UMAP visualization...\n")

umap_plot <- DimPlot(merged_obj, group.by = "sample_id", pt.size = 0.5) +
  theme_minimal() +
  ggtitle("Combined UMAP of Total NK, NKp46+, and NKp46- (dims25, res 0.3)") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ------------------------- #
# Save Plot
# ------------------------- #
output_file <- file.path(output_dir, "figures/combined_UMAP_TotalNK_NKp46_d25_res0.3.png")
ggsave(filename = output_file, plot = umap_plot, width = 10, height = 8, dpi = 600, bg = "transparent")
cat("✅ Combined UMAP saved to", output_file, "\n")
cat("🎉 UMAP generation complete. Outputs saved in:\n", output_dir, "\n")




# # ------------------------- #
# # ElbowPlot for PCA Dimension Selection
# # ------------------------- #
# elbow_plot <- ElbowPlot(combined_integrated, ndims = 30) +
#   ggtitle("Elbow Plot: PCA") +
#   theme_minimal()

# ggsave(output_dir, "/figures/ElbowPlot_PCA_dims.pdf",
#        plot = elbow_plot, width = 6, height = 4, dpi = 300)
# cat("📉 ElbowPlot saved\n")

# # ------------------------- #
# # Clustering at Multiple Resolutions
# # ------------------------- #
# resolutions <- c(0.3, 0.4, 0.5)

# # Find Neighbors once
# combined_integrated <- FindNeighbors(combined_integrated, dims = 1:25)

# # Perform clustering for each resolution
# for (res in resolutions) {
#   combined_integrated <- FindClusters(combined_integrated, resolution = res)
# }
# cat("✅ Clustering completed for resolutions:", paste(resolutions, collapse = ", "), "\n")

# # ------------------------- #
# # UMAP Plots for Each Resolution
# # ------------------------- #
# library(patchwork)

# umap_list <- list()
# for (res in resolutions) {
#   colname <- paste0("integrated_snn_res.", res)
#   p <- DimPlot(combined_integrated, group.by = colname, label = TRUE, pt.size = 0.5) +
#     ggtitle(paste("UMAP: Clustering at res =", res)) +
#     theme_minimal() +
#     theme(panel.background = element_rect(fill = "transparent", color = NA),
#           plot.background = element_rect(fill = "transparent", color = NA),
#           legend.background = element_rect(fill = "transparent"),
#           legend.box.background = element_rect(fill = "transparent"))
#   umap_list[[as.character(res)]] <- p
# }

# umap_res_plot <- wrap_plots(umap_list, ncol = 1)

# ggsave(output_dir, "/figures/UMAP_res_0.3_0.4_0.5.pdf",
#        plot = umap_res_plot, width = 8, height = 12, dpi = 600, bg = "transparent")
# cat("✅ UMAPs for multiple resolutions saved\n")

# # ------------------------- #
# # Save Updated Seurat Object
# # ------------------------- #
# saveRDS(combined_integrated, output_dir, "/rds/integrated_NK_all_clustered.rds")
# cat("💾 Final clustered object saved with clustering metadata\n")