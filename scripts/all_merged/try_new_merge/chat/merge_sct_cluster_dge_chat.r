# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(ggplot2)
library(cowplot)
library(hdf5r)
library(patchwork)
library(stringr)

# ------------------------- #
# Define Directory Structure
# ------------------------- #
output_dir <- "/home/outputs/new_nkp46_outputs_20250531/chat"
rds_output_dir <- file.path(output_dir, "rds")
pdf_output_dir <- file.path(output_dir, "pdf")
tsv_output_dir <- file.path(output_dir, "tsv")

dir.create(rds_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tsv_output_dir, recursive = TRUE, showWarnings = FALSE)

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

totalNK_file_paths <- build_h5_paths(totalNK_dir, "totalnk")
neg_file_paths <- build_h5_paths(neg_nkp46_dir, "ncr1neg")
pos_file_paths <- build_h5_paths(pos_nkp46_dir, "ncr1pos")

# ------------------------- #
# Load and Create Seurat Objects
# ------------------------- #
load_seurat_objects <- function(file_paths, sample_names, condition_label) {
    objs <- mapply(function(fp, sample) {
        mat <- Read10X_h5(fp)
        obj <- CreateSeuratObject(counts = mat, project = sample, min.features = 200)
        obj$sample <- str_to_title(sample)
        obj$sample_id <- obj$sample
        obj$condition <- condition_label
        obj$animal_id <- obj$sample
        return(obj)
    }, file_paths, sample_names, SIMPLIFY = FALSE)
    return(objs)
}

neg_objs <- load_seurat_objects(neg_file_paths, sample_names, "nkp46_neg")
pos_objs <- load_seurat_objects(pos_file_paths, sample_names, "nkp46_pos")
totalNK_objs <- load_seurat_objects(totalNK_file_paths, sample_names, "totalNK")

# ------------------------- #
# QC Filtering and Violin Plot
# ------------------------- #
qc_filter <- function(obj_list, condition, output_dir) {
    qc_plots <- list()
    for (i in seq_along(obj_list)) {
        obj <- obj_list[[i]]
        sample <- unique(obj$sample)
        obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 &
            nCount_RNA > 200 & nCount_RNA < 10000)
        p1 <- VlnPlot(obj, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
        p2 <- VlnPlot(obj, features = "nCount_RNA", pt.size = 0.1) + NoLegend()
        combined <- p1 + p2 + plot_annotation(title = paste("QC -", sample))
        qc_plots[[sample]] <- combined
        obj_list[[i]] <- obj
    }
    g <- wrap_plots(qc_plots, ncol = 1)
    ggsave(file.path(output_dir, paste0(condition, "_qc_violin.pdf")), g, width = 10, height = 6 * length(obj_list), dpi = 600)
    return(obj_list)
}

neg_objs <- qc_filter(neg_objs, "nkp46_neg", pdf_output_dir)
pos_objs <- qc_filter(pos_objs, "nkp46_pos", pdf_output_dir)
totalNK_objs <- qc_filter(totalNK_objs, "totalNK", pdf_output_dir)

# ------------------------- #
# SCTransform and Integration
# ------------------------- #
process_and_integrate <- function(obj_list, condition, dims = 1:25) {
    cat("📦 Running SCTransform + Integration for", condition, "\n")
    obj_list <- lapply(obj_list, SCTransform, verbose = FALSE)
    features <- SelectIntegrationFeatures(obj_list, nfeatures = 3000)
    anchors <- FindIntegrationAnchors(obj_list, anchor.features = features, dims = dims)
    integrated <- IntegrateData(anchors, dims = dims)
    DefaultAssay(integrated) <- "integrated"
    integrated <- ScaleData(integrated)
    integrated <- RunPCA(integrated, dims = dims)
    integrated <- RunUMAP(integrated, dims = dims)
    integrated <- FindNeighbors(integrated, dims = dims)
    integrated <- FindClusters(integrated, resolution = 0.3)
    return(integrated)
}

nkp46_neg_integrated <- process_and_integrate(neg_objs, "nkp46_neg")
nkp46_pos_integrated <- process_and_integrate(pos_objs, "nkp46_pos")
totalNK_integrated <- process_and_integrate(totalNK_objs, "totalNK")

# ------------------------- #
# Save .rds Objects
# ------------------------- #
saveRDS(nkp46_neg_integrated, file.path(rds_output_dir, "nkp46_neg.rds"))
saveRDS(nkp46_pos_integrated, file.path(rds_output_dir, "nkp46_pos.rds"))
saveRDS(totalNK_integrated, file.path(rds_output_dir, "totalNK.rds"))

# ------------------------- #
# Generate UMAP Plots by Sample
# ------------------------- #
label_map <- c("nkp46_neg" = "NKp46-", "nkp46_pos" = "NKp46+", "totalNK" = "Total NK")

plot_umap_split <- function(seurat_obj, condition_name) {
    p <- DimPlot(seurat_obj, group.by = "sample", label = TRUE, pt.size = 0.3) +
        ggtitle(paste("UMAP Split -", label_map[[condition_name]])) +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_brewer(palette = "Pastel2")
    ggsave(file.path(pdf_output_dir, paste0(condition_name, "_umap_by_sample.pdf")),
        p,
        width = 10, height = 8, dpi = 600
    )
}

plot_umap_split(nkp46_neg_integrated, "nkp46_neg")
plot_umap_split(nkp46_pos_integrated, "nkp46_pos")
plot_umap_split(totalNK_integrated, "totalNK")

# ------------------------- #
# Generate Combined UMAP by Cluster
# ------------------------- #
plot_umap_combined <- function(seurat_obj, condition_name) {
    p <- DimPlot(seurat_obj, group.by = "seurat_clusters", label = TRUE, pt.size = 0.3) +
        ggtitle(paste("UMAP Clusters -", label_map[[condition_name]])) +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_brewer(palette = "Pastel1")
    ggsave(file.path(pdf_output_dir, paste0(condition_name, "_umap_combined.pdf")),
        p,
        width = 10, height = 8, dpi = 600
    )
}

plot_umap_combined(nkp46_neg_integrated, "nkp46_neg")
plot_umap_combined(nkp46_pos_integrated, "nkp46_pos")
plot_umap_combined(totalNK_integrated, "totalNK")


# ------------------------- #
# Load Integrated Seurat Objects
# ------------------------- #
nkp46_neg <- readRDS(file.path(rds_output_dir, "nkp46_neg.rds"))
nkp46_pos <- readRDS(file.path(rds_output_dir, "nkp46_pos.rds"))
totalNK <- readRDS(file.path(rds_output_dir, "totalNK.rds"))

# ------------------------- #
# Rename Cells (avoid collisions)
# ------------------------- #
nkp46_neg <- RenameCells(nkp46_neg, add.cell.id = "nkp46_neg")
nkp46_pos <- RenameCells(nkp46_pos, add.cell.id = "nkp46_pos")
totalNK <- RenameCells(totalNK, add.cell.id = "totalNK")

# ------------------------- #
# Merge Integrated Objects
# ------------------------- #
cat("🔗 Merging all integrated objects...\n")
merged <- merge(nkp46_neg, y = list(nkp46_pos, totalNK), project = "NK_Merged", merge.data = TRUE)
DefaultAssay(merged) <- "RNA" # Set default assay to RNA for DGE
cat("✅ Merged object cells:", ncol(merged), "\n")

saveRDS(merged, file.path(rds_output_dir, "merged_NK_dim25_res0.3_chat.rds"))
cat("💾 Saved merged object to:", file.path(rds_output_dir, "merged_NK_dim25_res0.3_chat.rds"), "\n")







# ------------------------- #
# Set Idents and Run DGE
# ------------------------- #
Idents(merged) <- "seurat_clusters"
cat("📊 Performing cluster-based differential gene expression...\n")

# DGE: Each cluster vs. all others
cluster_markers <- FindAllMarkers(
    merged,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    assay = "RNA" # Use raw RNA expression for DGE
)

# ------------------------- #
# Save DGE Output
# ------------------------- #
write.table(cluster_markers, dge_output_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("💾 Saved DGE results to:", dge_output_path, "\n")

# Preview Top Markers
top_markers <- cluster_markers %>%
    group_by(cluster) %>%
    top_n(n = 3, wt = avg_log2FC)
print(top_markers)
