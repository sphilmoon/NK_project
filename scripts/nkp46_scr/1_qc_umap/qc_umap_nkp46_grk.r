# Script: process_nkp46_conditions.R
# Purpose: Process NKp46-negative and NKp46-positive conditions with Seurat, mirroring Total NK naming

# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)

# ------------------------- #
# Define Directory Structure
# ------------------------- #
output_dir <- "/home/outputs/nkp46_outputs/grok"

pos_output_dir <- file.path(output_dir, "pos")
neg_output_dir <- file.path(output_dir, "neg")

pos_rds_output_dir <- file.path(pos_output_dir, "rds")
neg_rds_output_dir <- file.path(neg_output_dir, "rds")

pos_pdf_output_dir <- file.path(pos_output_dir, "pdf")
neg_pdf_output_dir <- file.path(neg_output_dir, "pdf")

# Create directories
dir.create(pos_rds_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(neg_rds_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pos_pdf_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(neg_pdf_output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------- #
# Define Paths and Configs
# ------------------------- #
nkp46_neg_h5 <- "/home/rawdata/neg_nkp46/"
nkp46_pos_h5 <- "/home/rawdata/pos_nkp46/"
sample_names <- c("animal25", "animal26", "animal27", "animal28")

# Define and validate .h5 file paths
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
neg_file_paths <- build_h5_paths(nkp46_neg_h5, "ncr1neg")
pos_file_paths <- build_h5_paths(nkp46_pos_h5, "ncr1pos")

# ------------------------- #
# 1. Data Loading & Sample Assignment
# ------------------------- #

# Function to load H5 files and create Seurat objects
load_h5_to_seurat <- function(file_paths, condition, sample_names) {
  seurat_objects <- mapply(function(file, name) {
    data <- Read10X_h5(file)
    seurat_obj <- CreateSeuratObject(
      counts = data, project = paste(condition, name, sep = "_"),
      min.cells = 3, min.features = 200
    )
    seurat_obj$sample <- name # Tag with sample identity (mirrors Total NK)
    seurat_obj$condition <- condition # Additional field for condition
    return(seurat_obj)
  }, file_paths, sample_names, SIMPLIFY = FALSE)
  names(seurat_objects) <- sample_names
  return(seurat_objects)
}

# Load NKp46-negative data
nkp46_neg_objects <- load_h5_to_seurat(neg_file_paths, "NKp46neg", sample_names)
cat("✅ Loaded NKp46-negative Seurat objects\n")

# Load NKp46-positive data
nkp46_pos_objects <- load_h5_to_seurat(pos_file_paths, "NKp46pos", sample_names)
cat("✅ Loaded NKp46-positive Seurat objects\n")

# Print sample information
for (name in names(nkp46_neg_objects)) {
  cat("NKp46neg - Sample:", name, "\n")
  cat("Number of cells:", ncol(nkp46_neg_objects[[name]]), "\n")
  cat("Number of genes:", nrow(nkp46_neg_objects[[name]]), "\n\n")
}
for (name in names(nkp46_pos_objects)) {
  cat("NKp46pos - Sample:", name, "\n")
  cat("Number of cells:", ncol(nkp46_pos_objects[[name]]), "\n")
  cat("Number of genes:", nrow(nkp46_pos_objects[[name]]), "\n\n")
}

# Save initial objects (mirroring Total NK naming)
saveRDS(nkp46_neg_objects, file.path(neg_rds_output_dir, "seurat_objects_list_wo_51_52.rds"))
saveRDS(nkp46_pos_objects, file.path(pos_rds_output_dir, "seurat_objects_list_wo_51_52.rds"))
cat("✅ Initial Seurat objects saved to", neg_rds_output_dir, "and", pos_rds_output_dir, "\n")

# ------------------------- #
# 2. Quality Control & Filtering
# ------------------------- #

# Function to perform QC and generate violin plots
perform_qc <- function(seurat_list, condition, pdf_output_dir) {
  combined_plots <- list()
  for (i in seq_along(seurat_list)) {
    animal_name <- names(seurat_list)[i]
    seurat_list[[i]] <- subset(seurat_list[[i]],
      subset = nFeature_RNA > 200 & nFeature_RNA < 3000 &
        nCount_RNA > 200 & nCount_RNA < 10000
    )
    cat("✅ Filtered", condition, animal_name, "\n")

    # Violin plots (mirroring Total NK titles)
    p1 <- VlnPlot(seurat_list[[i]], features = "nFeature_RNA", ncol = 1) +
      ggtitle("Gene counts QC") +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
    p2 <- VlnPlot(seurat_list[[i]], features = "nCount_RNA", ncol = 1) +
      ggtitle("UMI counts QC") +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
    combined_plot <- p1 + p2 + plot_layout(ncol = 2) +
      plot_annotation(title = paste("QC for", animal_name)) # Match Total NK title
    combined_plots[[animal_name]] <- combined_plot
  }
  stacked_plot <- plot_grid(plotlist = combined_plots, ncol = 1)
  output_file <- file.path(pdf_output_dir, "qc_violin_stacked_wo_51_52.pdf") # Match Total NK naming
  ggsave(output_file, plot = stacked_plot, width = 10, height = 6 * length(seurat_list), dpi = 600)
  cat("✅ Saved QC violin plots for", condition, "to", output_file, "\n")
  return(seurat_list)
}

nkp46_neg_objects <- perform_qc(nkp46_neg_objects, "NKp46neg", neg_pdf_output_dir)
nkp46_pos_objects <- perform_qc(nkp46_pos_objects, "NKp46pos", pos_pdf_output_dir)

# Save QC'd objects (mirroring Total NK naming)
saveRDS(nkp46_neg_objects, file.path(neg_rds_output_dir, "seurat_objects_qc_sct_wo_51_52.rds"))
saveRDS(nkp46_pos_objects, file.path(pos_rds_output_dir, "seurat_objects_qc_sct_wo_51_52.rds"))
cat("✅ QC'd Seurat objects saved to", neg_rds_output_dir, "and", pos_rds_output_dir, "\n")

# ------------------------- #
# 3. Integration per Condition
# ------------------------- #

# Function to integrate and process a condition
process_condition <- function(seurat_list, condition, rds_output_dir) {
  # Normalize with SCTransform
  seurat_list <- lapply(seurat_list, SCTransform, verbose = FALSE)
  saveRDS(seurat_list, file.path(rds_output_dir, "seurat_objects_qc_sct_wo_51_52.rds")) # Already saved above, but keeping for consistency
  cat("✅ SCTransformed", condition, "objects saved to", file.path(rds_output_dir, "seurat_objects_qc_sct_wo_51_52.rds"), "\n")

  # Select integration features
  features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)

  # Find integration anchors and integrate data
  anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features, dims = 1:25)
  integrated_data <- IntegrateData(anchorset = anchors, dims = 1:25)
  saveRDS(integrated_data, file.path(rds_output_dir, "integrated_data_wo_51_52_dims25.rds")) # Match Total NK naming
  cat("✅ Integrated", condition, "data saved to", file.path(rds_output_dir, "integrated_data_wo_51_52_dims25.rds"), "\n")

  # Set default assay and scale data
  DefaultAssay(integrated_data) <- "integrated"
  integrated_data <- ScaleData(integrated_data)

  # Run PCA
  integrated_data <- RunPCA(integrated_data, verbose = FALSE, npcs = 25)
  saveRDS(integrated_data, file.path(rds_output_dir, "integrated_data_wo_51_52_dims25_pca.rds"))
  cat("✅ PCA for", condition, "saved to", file.path(rds_output_dir, "integrated_data_wo_51_52_dims25_pca.rds"), "\n")

  # Find neighbors and clusters
  integrated_data <- FindNeighbors(integrated_data, dims = 1:25)
  integrated_data <- FindClusters(integrated_data, resolution = 0.3) # As requested, differs from Total NK (0.5)
  saveRDS(integrated_data, file.path(rds_output_dir, "integrated_data_wo_51_52_dims25_clusters.rds"))
  cat("✅ Clusters for", condition, "saved to", file.path(rds_output_dir, "integrated_data_wo_51_52_dims25_clusters.rds"), "\n")

  # Run UMAP
  integrated_data <- RunUMAP(integrated_data, dims = 1:25)
  saveRDS(integrated_data, file.path(rds_output_dir, "integrated_data_wo_51_52_dims25_umap.rds")) # Match Total NK naming
  cat("✅ UMAP for", condition, "saved to", file.path(rds_output_dir, "integrated_data_wo_51_52_dims25_umap.rds"), "\n")

  return(integrated_data)
}

nkp46_neg_integrated <- process_condition(nkp46_neg_objects, "nkp46_neg", neg_rds_output_dir)
nkp46_pos_integrated <- process_condition(nkp46_pos_objects, "nkp46_pos", pos_rds_output_dir)

# ------------------------- #
# 5. Visualization
# ------------------------- #

# Function to generate UMAP plots
generate_umap_plots <- function(seurat_obj, condition, pdf_output_dir) {
  # Split UMAP by sample (mirroring Total NK style)
  umap_by_sample <- DimPlot(seurat_obj,
    group.by = "sample", split.by = "sample", ncol = 2,
    label = TRUE, label.size = 5, repel = TRUE, pt.size = 0.3
  ) +
    ggtitle(paste("UMAP - Integration by Sample (dims25)")) + # Match Total NK title
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_brewer(palette = "Set2")
  ggsave(file.path(pdf_output_dir, "umap_dims25_split_by_sample.pdf"), # Adjusted naming
    plot = umap_by_sample, width = 12, height = 8, dpi = 600
  )

  # Combined UMAP (mirroring Total NK style)
  umap_combined <- DimPlot(seurat_obj,
    group.by = "sample", label = TRUE, label.size = 5,
    repel = TRUE, pt.size = 0.3
  ) +
    ggtitle(paste("UMAP - Integration by Sample (dims25)")) + # Match Total NK title
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right") +
    scale_color_brewer(palette = "Set2")
  ggsave(file.path(pdf_output_dir, "umap_dims25.pdf"), # Match Total NK naming
    plot = umap_combined, width = 10, height = 8, dpi = 600
  )

  cat("✅ Saved UMAP plots for", condition, "to", pdf_output_dir, "\n")
}

generate_umap_plots(nkp46_neg_integrated, "nkp46_neg", neg_pdf_output_dir)
generate_umap_plots(nkp46_pos_integrated, "nkp46_pos", pos_pdf_output_dir)

# ------------------------- #
# Final Message
# ------------------------- #
cat("🎉 Processing complete. Outputs saved in:", output_dir, "\n")
