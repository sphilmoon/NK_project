# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(ggplot2)
library(cowplot)

# ------------------------- #
# Define Input and Output Paths
# ------------------------- #
nkp_neg_rds <- "/home/outputs/nkp46_outputs/chat/neg/rds/nkp46_neg.rds"
nkp_pos_rds <- "/home/outputs/nkp46_outputs/chat/pos/rds/nkp46_pos.rds"
totalNK_rds <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/rds/integrated_data_dims25_res0.3_genecounts.rds"

output_dir <- "/home/outputs/all_merged_TotalNK_nkp46/chat"
all_merged_rds_dir <- file.path(output_dir, "rds")
pdf_dir <- file.path(output_dir, "pdf")

# Create output directories if they don't exist
dir.create(all_merged_rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------- #
# Load Seurat Objects
# ------------------------- #
nkp_neg <- readRDS(nkp_neg_rds)
nkp_pos <- readRDS(nkp_pos_rds)
totalNK <- readRDS(totalNK_rds)

# ------------------------- #
# Define UMAP Plotting Function
# ------------------------- #
plot_umap <- function(seurat_obj, condition_label, pdf_dir) {
    # Ensure 'sample' metadata exists
    if (!"sample" %in% colnames(seurat_obj@meta.data)) {
        stop(paste0("The 'sample' metadata column is missing in ", condition_label, " object."))
    }

    # UMAP split by sample
    umap_split <- DimPlot(seurat_obj, group.by = "sample", pt.size = 0.3, split.by = "sample") +
        ggtitle(paste0("UMAP Split by Sample - ", condition_label)) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(file.path(pdf_dir, paste0(condition_label, "_umap_split.pdf")), plot = umap_split, width = 10, height = 8, dpi = 600)

    # UMAP combined clusters
    umap_combined <- DimPlot(seurat_obj, group.by = "seurat_clusters", pt.size = 0.3) +
        ggtitle(paste0("UMAP (Clusters) - ", condition_label)) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(file.path(pdf_dir, paste0(condition_label, "_umap_combined.pdf")), plot = umap_combined, width = 10, height = 8, dpi = 600)
}

# ------------------------- #
# Generate UMAPs for Each Condition
# ------------------------- #
plot_umap(nkp_neg, "nkp46_neg", pdf_dir)
plot_umap(nkp_pos, "nkp46_pos", pdf_dir)
plot_umap(totalNK, "totalNK", pdf_dir)

# ------------------------- #
# Merge All Three Conditions
# ------------------------- #
# Add condition labels
nkp_neg$condition <- "nkp46_neg"
nkp_pos$condition <- "nkp46_pos"
totalNK$condition <- "totalNK"

# Merge Seurat objects
merged_obj <- merge(nkp_neg, y = list(nkp_pos, totalNK), add.cell.ids = c("nkp46_neg", "nkp46_pos", "totalNK"))

# ------------------------- #
# Re-normalize, Scale, and Cluster
# ------------------------- #
# Normalize and scale data
merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj, selection.method = "vst", nfeatures = 2000)
merged_obj <- ScaleData(merged_obj)

# Perform PCA
merged_obj <- RunPCA(merged_obj, npcs = 50)

# Generate and save ElbowPlot
elbow_plot <- ElbowPlot(merged_obj, ndims = 50)
ggsave(file.path(pdf_dir, "elbow_plot.pdf"), plot = elbow_plot, width = 8, height = 6, dpi = 600)

# Based on ElbowPlot, select number of PCs (e.g., 25)
dims_to_use <- 1:25

# Run UMAP and clustering
merged_obj <- RunUMAP(merged_obj, dims = dims_to_use)
merged_obj <- FindNeighbors(merged_obj, dims = dims_to_use)
merged_obj <- FindClusters(merged_obj, resolution = 0.5)

# Save merged Seurat object
saveRDS(merged_obj, file = file.path(all_merged_rds_dir, "merged_totalNK_nkp46.rds"))

# ------------------------- #
# Generate UMAP for Merged Object
# ------------------------- #
# UMAP by condition
umap_condition <- DimPlot(merged_obj, group.by = "condition", pt.size = 0.3) +
    ggtitle("UMAP by Condition") +
    theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(pdf_dir, "merged_umap_by_condition.pdf"), plot = umap_condition, width = 10, height = 8, dpi = 600)

# UMAP by clusters
umap_clusters <- DimPlot(merged_obj, group.by = "seurat_clusters", pt.size = 0.3) +
    ggtitle("UMAP by Clusters") +
    theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(pdf_dir, "merged_umap_by_clusters.pdf"), plot = umap_clusters, width = 10, height = 8, dpi = 600)
