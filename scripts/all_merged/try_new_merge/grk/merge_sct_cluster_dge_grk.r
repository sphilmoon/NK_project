# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)

# ------------------------- #
# Define Directory Structure
# ------------------------- #
output_dir <- "/home/outputs/new_nkp46_outputs_20250531"
rds_output_dir <- file.path(output_dir, "rds")
pdf_output_dir <- file.path(output_dir, "pdf")

dir.create(rds_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------- #
# Define Input Directories
# ------------------------- #
totalNK_dir <- "/home/rawdata/totalNK"
neg_nkp46_dir <- "/home/rawdata/neg_nkp46"
pos_nkp46_dir <- "/home/rawdata/pos_nkp46"
sample_names <- c("animal25", "animal26", "animal27", "animal28")

# ------------------------- #
# Build File Paths (Validated)
# ------------------------- #
build_h5_paths <- function(base_dir, condition_prefix) {
    fnames <- paste0(tolower(sample_names), "_", condition_prefix, "_filtered_feature_bc_matrix.h5")
    full_paths <- file.path(base_dir, fnames)

    missing <- !file.exists(full_paths)
    if (any(missing)) {
        cat("❌ Missing files:\n")
        print(full_paths[missing])
        stop("Some .h5 files are missing. Please check the paths or filenames.")
    }

    return(full_paths)
}

# Build paths and validate existence
totalNK_file_paths <- build_h5_paths(totalNK_dir, "totalnk")
neg_file_paths <- build_h5_paths(neg_nkp46_dir, "ncr1neg")
pos_file_paths <- build_h5_paths(pos_nkp46_dir, "ncr1pos")

# ------------------------- #
# Load and Process Individual Samples
# ------------------------- #
cat("📦 Loading and processing individual samples...\n")

# Function to create Seurat objects
create_seurat_object <- function(h5_path, sample_id, condition) {
    seurat_obj <- Read10X_h5(h5_path)
    seurat_obj <- CreateSeuratObject(counts = seurat_obj, min.cells = 3, min.features = 200)
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 &
        nCount_RNA > 200 & nCount_RNA < 10000)
    seurat_obj$sample_id <- sample_id
    seurat_obj$condition <- condition
    seurat_obj$animal_id <- gsub("animal(\\d+).*", "\\1", sample_id)
    return(seurat_obj)
}

# Load all samples into a list
sample_list <- c(
    mapply(create_seurat_object, totalNK_file_paths, paste0(sample_names, "_totalnk"), "totalNK", SIMPLIFY = FALSE),
    mapply(create_seurat_object, neg_file_paths, paste0(sample_names, "_ncr1neg"), "nkp46_neg", SIMPLIFY = FALSE),
    mapply(create_seurat_object, pos_file_paths, paste0(sample_names, "_ncr1pos"), "nkp46_pos", SIMPLIFY = FALSE)
)

cat("✅ Loaded and processed", length(sample_list), "samples\n")

# ------------------------- #
# Merge and Integrate Data
# ------------------------- #
cat("🔗 Merging and integrating data...\n")

# Initial merge (optional, for reference)
merged_raw <- merge(sample_list[[1]], y = sample_list[-1], add.cell.ids = sapply(sample_list, function(x) x$sample_id[1]))

# Split by animal for integration
merged_list <- SplitObject(merged_raw, split.by = "animal_id")

# Apply SCTransform for normalization and batch correction
merged_list <- lapply(merged_list, SCTransform, vars.to.regress = "percent.mt", verbose = FALSE)

# Select integration features
features <- SelectIntegrationFeatures(object.list = merged_list, nfeatures = 3000)
merged_list <- PrepSCTIntegration(object.list = merged_list, anchor.features = features)

# Find integration anchors and integrate
anchors <- FindIntegrationAnchors(object.list = merged_list, normalization.method = "SCT", anchor.features = features)
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
DefaultAssay(integrated) <- "integrated"

cat("✅ Integration complete: total cells =", ncol(integrated), "\n")

# ------------------------- #
# Dimensionality Reduction and Clustering
# ------------------------- #
cat("🔬 Performing dimensionality reduction and clustering...\n")

# Run PCA
integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)
cat("📊 Saving ElbowPlot...\n")
elbow_plot <- ElbowPlot(integrated, ndims = 50)
ggsave(file.path(pdf_output_dir, "elbow_plot_20250531.pdf"), elbow_plot, width = 8, height = 6, dpi = 600)

# Run UMAP and clustering
integrated <- RunUMAP(integrated, dims = 1:25, reduction = "pca")
integrated <- FindNeighbors(integrated, dims = 1:25)
integrated <- FindClusters(integrated, resolution = seq(0.2, 1.0, by = 0.2)) # Test multiple resolutions

# Select optimal resolution (e.g., based on silhouette score or visual inspection)
# For simplicity, use resolution 0.5 (adjust based on downstream analysis)
integrated$seurat_clusters <- integrated$SCT_snn_res.0.5
cat("✅ Clustering complete with resolution 0.5\n")

cat("💾 Saving integrated Seurat object...\n")
saveRDS(integrated, file.path(rds_output_dir, "integrated_NK_dims25_res0.3_grk_20250531.rds"))
cat("✅ Integrated object saved\n")



# ------------------------- #
# Visualize Clusters and Conditions
# ------------------------- #
cat("🎨 Generating visualizations...\n")

# UMAP by clusters
umap_clusters <- DimPlot(
    integrated,
    group.by = "seurat_clusters",
    label = TRUE,
    label.size = 5,
    pt.size = 0.3,
    repel = TRUE
) +
    ggtitle("UMAP with Numbered Clusters (Merged NK Cells)") +
    theme(plot.title = element_text(hjust = 0.5, size = 14))
ggsave(file.path(pdf_output_dir, "UMAP_clusters_20250531.pdf"), umap_clusters, width = 10, height = 8, dpi = 600)

# UMAP by condition with display labels
condition_display <- c("nkp46_neg" = "NKp46-", "nkp46_pos" = "NKp46+", "totalNK" = "Total NK")
umap_condition <- DimPlot(
    integrated,
    group.by = "condition",
    label = FALSE,
    pt.size = 0.3
) +
    ggtitle("UMAP by Condition") +
    scale_color_manual(values = c("NKp46-" = "#1b9e77", "NKp46+" = "#d95f02", "Total NK" = "#7570b3")) +
    theme(plot.title = element_text(hjust = 0.5, size = 14))
ggsave(file.path(pdf_output_dir, "UMAP_condition_20250531.pdf"), umap_condition, width = 10, height = 8, dpi = 600)

# UMAP grid by animal and condition
animals <- c("animal25", "animal26", "animal27", "animal28")
plot_list <- list()

for (animal in animals) {
    for (cond in names(condition_display)) {
        cells <- WhichCells(integrated, expression = animal_id == animal & condition == cond)
        cat(sprintf("🔍 Found %d cells for %s x %s\n", length(cells), animal, cond))

        if (length(cells) == 0) {
            cat("⚠️ No cells found for", animal, cond, "- skipping\n")
            p <- ggplot() +
                theme_void()
        } else {
            p <- DimPlot(
                integrated,
                cells = cells,
                group.by = "condition",
                label = FALSE,
                pt.size = 0.3
            ) +
                ggtitle(NULL) +
                theme_void() +
                theme(legend.position = "none") +
                scale_color_manual(values = condition_display, labels = condition_display)
        }
        plot_list[[paste(animal, cond, sep = "_")]] <- p
    }
}

# Add column titles (condition display names)
column_titles <- lapply(condition_display, function(label) {
    ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = label, size = 6, fontface = "bold") +
        theme_void()
})
top_row <- wrap_plots(column_titles, ncol = length(condition_display))

# Add row titles (animal IDs)
row_titles <- lapply(animals, function(animal) {
    ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = animal, angle = 90, size = 6, fontface = "bold") +
        theme_void()
})
left_col <- wrap_plots(row_titles, ncol = 1)

# Assemble UMAP grid
umap_matrix <- lapply(animals, function(animal) {
    plots_row <- lapply(names(condition_display), function(cond) {
        plot_list[[paste(animal, cond, sep = "_")]]
    })
    wrap_plots(plots_row, ncol = length(condition_display))
})
umap_full <- wrap_plots(umap_matrix, ncol = 1)

# Add labels: top (column), left (row)
umap_labeled <- plot_grid(
    plot_grid(NULL, top_row, ncol = 2, rel_widths = c(0.12, 1)),
    plot_grid(left_col, umap_full, ncol = 2, rel_widths = c(0.12, 1)),
    nrow = 2, rel_heights = c(0.08, 1)
)

ggsave(file.path(pdf_output_dir, "UMAP_animal_condition_grid_20250531.pdf"),
    umap_labeled,
    width = 12, height = 16, dpi = 600
)
cat("✅ Visualizations complete\n")

# ------------------------- #
# Identify NK Subpopulations (Marker Genes)
# ------------------------- #
cat("🔍 Identifying NK subpopulations...\n")
Idents(integrated) <- "seurat_clusters"
cluster_markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster_markers, file.path(output_dir, "cluster_markers_20250531.csv"), row.names = FALSE)
cat("✅ Marker genes identified and saved\n")

# ------------------------- #
# Differential Gene Expression (DGE) Analysis
# ------------------------- #
cat("📊 Performing DGE analysis...\n")

# DGE by cluster (across all conditions)
cluster_dge <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster_dge, file.path(output_dir, "cluster_dge_20250531.csv"), row.names = FALSE)

# Optional: DGE by condition within clusters (example for cluster 0)
Idents(integrated) <- "condition"
condition_dge <- FindMarkers(integrated, ident.1 = "nkp46_pos", ident.2 = "nkp46_neg", subset = seurat_clusters == "0")
write.csv(condition_dge, file.path(output_dir, "condition_dge_cluster0_20250531.csv"), row.names = FALSE)

cat("✅ DGE analysis complete\n")

# ------------------------- #
# Save Integrated Object
# ------------------------- #
cat("💾 Saving integrated Seurat object...\n")
saveRDS(integrated, file.path(rds_output_dir, "integrated_NK_subpopulations_20250531.rds"))
cat("✅ Integrated object saved\n")
