# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(Matrix)

# ------------------------- #
# Define Paths and Configs
# ------------------------- #
totalNK_rds <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/rds/integrated_data_dims25_res0.3_genecounts.rds"
nkp46_rds <- "/home/outputs/nkp46_outputs/nkp46_integrated_data.rds"
output_dir <- "/home/outputs/all_merged_TotalNK_nkp46"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create subdirectories for RDS and figures
rds_dir <- file.path(output_dir, "rds")
figures_dir <- file.path(output_dir, "figures")
dir.create(rds_dir, recursive = TRUE)
dir.create(figures_dir, recursive = TRUE)

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
  
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(counts, "dgCMatrix")
  }
  
  values <- summary(counts)$x
  if (!all(is.numeric(values))) {
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
nkp46_pos$sample_id <- "NKp46+"
nkp46_neg$sample_id <- "NKp46-"
cat("✅ Annotated sample identities\n")

# ------------------------- #
# Merge Seurat Objects
# ------------------------- #
merged_obj <- merge(totalNK, y = c(nkp46_pos, nkp46_neg), 
                    add.cell.ids = c("TotalNK", "NKp46+", "NKp46-"), 
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
# Save RDS File Before PCA
# ------------------------- #
rds_file <- file.path(rds_dir, "all_totalNK_nkp46.rds")
saveRDS(merged_obj, file = rds_file)
cat("✅ Saved merged Seurat object to", rds_file, "\n")

# Validate saved RDS by reloading
cat("🔍 Validating saved RDS file...\n")
reloaded_obj <- readRDS(rds_file)
if (all(dim(reloaded_obj) == dim(merged_obj)) && 
    identical(DefaultAssay(reloaded_obj), DefaultAssay(merged_obj))) {
  cat("✅ RDS file validated: Dimensions and assay match\n")
} else {
  stop("❌ RDS validation failed: Dimensions or assay mismatch")
}

# ------------------------- #
# Run PCA with Consistent Dims
# ------------------------- #
merged_obj <- RunPCA(merged_obj, npcs = 25)
cat("✅ Ran PCA with 25 dimensions\n")

# ------------------------- #
# Cluster and UMAP with Consistent Resolution
# ------------------------- #
merged_obj <- FindNeighbors(merged_obj, dims = 1:25)
merged_obj <- FindClusters(merged_obj, resolution = 0.5)
merged_obj <- RunUMAP(merged_obj, dims = 1:25)
cat("✅ Ran clustering and UMAP with resolution 0.5\n")

# ------------------------- #
# Define Common Theme for UMAP Plots
# ------------------------- #
umap_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

# ------------------------- #
# Create Combined UMAP Visualization
# ------------------------- #
cat("🎨 Creating combined UMAP visualization...\n")

combined_umap <- DimPlot(merged_obj, group.by = "sample_id", pt.size = 0.5) +
  ggtitle("Combined UMAP of Total NK, NKp46+, and NKp46- (dims25, res 0.5)") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  umap_theme

# Save combined UMAP
combined_output_file <- file.path(figures_dir, "combined_UMAP_TotalNK_NKp46_d25_res0.5.pdf")
ggsave(filename = combined_output_file, plot = combined_umap, width = 10, height = 8, dpi = 600, bg = "transparent")
cat("✅ Combined UMAP saved to", combined_output_file, "\n")

# ------------------------- #
# Create Separate UMAPs for Each Sample Identity
# ------------------------- #
cat("🎨 Creating separate UMAP visualizations for each sample identity...\n")

# Split UMAP by sample_id
split_umap <- DimPlot(merged_obj, group.by = "sample_id", split.by = "sample_id", pt.size = 0.5, ncol = 3) +
  ggtitle("UMAP Split by Sample Identity (dims25, res 0.5)") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  umap_theme

# Save split UMAPs
split_output_file <- file.path(figures_dir, "split_UMAP_TotalNK_NKp46_d25_res0.5.pdf")
ggsave(filename = split_output_file, plot = split_umap, width = 15, height = 5, dpi = 600, bg = "transparent")
cat("✅ Split UMAPs saved to", split_output_file, "\n")

cat("🎉 UMAP generation complete. Outputs saved in:\n", output_dir, "\n")