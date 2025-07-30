# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(ggplot2)
library(cowplot)
library(hdf5r)
library(dplyr)
library(cluster)

# ------------------------- #
# Define Directory Structure
# ------------------------- #
output_dir <- "/home/outputs/nkp46_outs_20250730_slurm"

pos_output_dir <- file.path(output_dir, "pos")
neg_output_dir <- file.path(output_dir, "neg")

pos_rds_output_dir <- file.path(pos_output_dir, "rds")
neg_rds_output_dir <- file.path(neg_output_dir, "rds")

pos_pdf_output_dir <- file.path(pos_output_dir, "pdf")
neg_pdf_output_dir <- file.path(neg_output_dir, "pdf")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pos_rds_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(neg_rds_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pos_pdf_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(neg_pdf_output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------- #
# Define Input Directories
# ------------------------- #
nkp46_neg_dir <- "/home/rawdata/neg_nkp46"
nkp46_pos_dir <- "/home/rawdata/pos_nkp46"
sample_names <- c("animal25", "animal26", "animal27", "animal28")

build_h5_paths <- function(base_dir, condition_prefix) {
  fnames <- paste0(tolower(sample_names), "_", condition_prefix, "_filtered_feature_bc_matrix.h5")
  full_paths <- file.path(base_dir, fnames)
  missing <- !file.exists(full_paths)
  if (any(missing)) {
    cat("\u274c Missing files:\n")
    print(full_paths[missing])
    stop("Some .h5 files are missing. Please check the paths or filenames.")
  }
  return(full_paths)
}

neg_file_paths <- build_h5_paths(nkp46_neg_dir, "ncr1neg")
pos_file_paths <- build_h5_paths(nkp46_pos_dir, "ncr1pos")

dims_to_use <- 1:15
resolutions <- c(0.25)
dims_list <- c(15)

qc_process_condition <- function(h5_files, sample_names, condition_label, rds_dir, pdf_dir) {
  seurat_objects <- mapply(function(file, name) {
    data <- Read10X_h5(file)
    obj <- CreateSeuratObject(counts = data, project = name, min.cells = 3, min.features = 200)
    obj$sample <- name
    obj$condition <- condition_label
    obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & nCount_RNA > 200 & nCount_RNA < 10000)
    return(obj)
  }, h5_files, sample_names, SIMPLIFY = FALSE)

  names(seurat_objects) <- sample_names

  plots <- lapply(seurat_objects, function(obj) {
    p1 <- VlnPlot(obj, features = "nFeature_RNA") + NoLegend()
    p2 <- VlnPlot(obj, features = "nCount_RNA") + NoLegend()
    plot_grid(p1, p2, ncol = 2)
  })

  qc_plot <- plot_grid(plotlist = plots, ncol = 1)
  ggsave(file.path(pdf_dir, paste0("nkp46_", condition_label, "_qc_violin.pdf")), qc_plot, width = 10, height = 6 * length(plots))

  seurat_objects <- lapply(seurat_objects, SCTransform, verbose = FALSE)
  saveRDS(seurat_objects, file.path(rds_dir, paste0("nkp46_", condition_label, "_objects.rds")))
  return(seurat_objects)
}

neg_objs <- qc_process_condition(neg_file_paths, sample_names, "neg", neg_rds_output_dir, neg_pdf_output_dir)
pos_objs <- qc_process_condition(pos_file_paths, sample_names, "pos", pos_rds_output_dir, pos_pdf_output_dir)

merge_conditions <- function(seurat_objects, condition_label, rds_output_dir) {
  cat("\U0001F517 Merging all animals for condition:", condition_label, "...\n")

  features <- SelectIntegrationFeatures(object.list = seurat_objects, nfeatures = 3000)
  seurat_objects <- PrepSCTIntegration(object.list = seurat_objects, anchor.features = features)
  anchors <- FindIntegrationAnchors(object.list = seurat_objects, normalization.method = "SCT", anchor.features = features)
  merged <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

  DefaultAssay(merged) <- "integrated"
  merged <- ScaleData(merged)
  merged <- RunPCA(merged, npcs = max(dims_list))

  for (dim in dims_list) {
    dims_to_use <- 1:dim
    merged <- RunUMAP(merged, dims = dims_to_use, reduction.name = paste0("umap_dims", dim))
    for (res in resolutions) {
      merged <- FindNeighbors(merged, dims = dims_to_use)
      merged <- FindClusters(merged, resolution = res)
    }
  }

  saveRDS(merged, file.path(rds_output_dir, paste0("nkp46_", condition_label, "_merged.rds")))
  cat("\u2705 Saved merged Seurat object for", condition_label, "to", rds_output_dir, "\n")
}

merge_conditions(neg_objs, "neg", neg_rds_output_dir)
merge_conditions(pos_objs, "pos", pos_rds_output_dir)


# ----------------------------- #
# UMAP Plot for Merged SCT Output (neg and pos)
# ----------------------------- #

plot_merged_umap <- function(rds_file, pdf_dir, condition_label) {
  merged <- readRDS(rds_file)
  merged$sample <- as.factor(merged$sample)

  p1 <- DimPlot(merged, group.by = "seurat_clusters", label = TRUE, pt.size = 0.3) +
    ggtitle(paste("UMAP: Merged Clusters (SCT) -", condition_label)) +
    theme_minimal()

  p2 <- DimPlot(merged, group.by = "sample", pt.size = 0.3) +
    ggtitle(paste("UMAP: By Sample (SCT) -", condition_label)) +
    theme_minimal()

  ggsave(file.path(pdf_dir, paste0("nkp46_", condition_label, "_umap_clusters.pdf")),
         plot = p1, width = 9, height = 7, dpi = 600)

  ggsave(file.path(pdf_dir, paste0("nkp46_", condition_label, "_umap_by_sample.pdf")),
         plot = p2, width = 9, height = 7, dpi = 600)

  cat("✅ UMAPs saved for condition:", condition_label, "\n")
}

# Run for both conditions
plot_merged_umap(
  rds_file = file.path(neg_rds_output_dir, "nkp46_neg_merged.rds"),
  pdf_dir = neg_pdf_output_dir,
  condition_label = "neg"
)

plot_merged_umap(
  rds_file = file.path(pos_rds_output_dir, "nkp46_pos_merged.rds"),
  pdf_dir = pos_pdf_output_dir,
  condition_label = "pos"
)
