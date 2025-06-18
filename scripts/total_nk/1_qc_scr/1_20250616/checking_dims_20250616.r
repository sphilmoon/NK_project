# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(stringr)


# ------------------------- #
# Define Directory Structure
# ------------------------- #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
rds_output_dir <- file.path(output_dir, "rds")
pdf_output_dir <- file.path(output_dir, "pdf")
tsv_output_dir <- file.path(output_dir, "tsv")
dir.create(rds_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tsv_output_dir, recursive = TRUE, showWarnings = FALSE)


# # ------------------------- #
# # Define Input Directories
# # ------------------------- #
# totalNK_h5_dir <- "/home/rawdata/totalNK"
# animal_names <- c("animal25", "animal26", "animal27", "animal28")


# # ------------------------- #
# # Build File Paths
# # ------------------------- #
# build_h5_paths <- function(base_dir, condition_prefix) {
#  fnames <- paste0(tolower(animal_names), "_", condition_prefix, "_filtered_feature_bc_matrix.h5")
#  full_paths <- file.path(base_dir, fnames)
#  missing <- !file.exists(full_paths)
#  if (any(missing)) {
#    cat("❌ Missing files:\n")
#    print(full_paths[missing])
#    stop("Some .h5 files are missing. Please check the paths or filenames.")
#  }
#  names(full_paths) <- animal_names
#  return(full_paths)
# }
# totalNK_file_paths <- build_h5_paths(totalNK_h5_dir, "totalnk")


# # ------------------------- #
# # Load and QC Each Animal
# # ------------------------- #
# seurat_objects <- list()
# combined_plots <- list()
# perform_qc <- function(seurat_obj, sample_name) {
#  qc_params <- list(
#    animal27 = list(nFeature_min = 200, nFeature_max = 3000, nCount_min = 250, nCount_max = 7500),
#    animal28 = list(nFeature_min = 200, nFeature_max = 3000, nCount_min = 250, nCount_max = 10000),
#    default  = list(nFeature_min = 200, nFeature_max = 3000, nCount_min = 500, nCount_max = 10000)
#  )
#  params <- if (sample_name %in% names(qc_params)) qc_params[[sample_name]] else qc_params$default
#  seurat_obj <- subset(
#    seurat_obj,
#    subset = nFeature_RNA > params$nFeature_min & nFeature_RNA < params$nFeature_max &
#             nCount_RNA > params$nCount_min & nCount_RNA < params$nCount_max
#  )
#  return(seurat_obj)
# }

# for (i in seq_along(totalNK_file_paths)) {
#  sample <- names(totalNK_file_paths)[i]
#  filepath <- totalNK_file_paths[[i]]
#  cat("📥 Reading:", sample, "\n")
#  seurat_obj <- Read10X_h5(filepath) %>%
#    CreateSeuratObject(project = sample, min.cells = 3, min.features = 200)
#  seurat_obj$animal <- sample
#  # Add mitochondrial percentage
#  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
#  # QC filter
#  seurat_obj <- perform_qc(seurat_obj, sample)
#  # Normalize
#  seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
#  # Save object
#  seurat_objects[[sample]] <- seurat_obj
#  # Generate QC plots
#  p1 <- VlnPlot(seurat_obj, features = "nFeature_RNA", pt.size = 0) + ggtitle("Gene counts")
#  p2 <- VlnPlot(seurat_obj, features = "nCount_RNA", pt.size = 0) + ggtitle("UMI counts")
#  combined_plot <- p1 + p2 + plot_layout(ncol = 2) + plot_annotation(title = paste("QC for", sample))
#  combined_plots[[sample]] <- combined_plot
#  cat("✅ QC and normalization complete for", sample, "\n")
# }
# # Save all QC plots
# pdf(file.path(pdf_output_dir, "qc_totalNK_animals_20250616.pdf"), width = 10, height = 6 * length(combined_plots))
# for (plot in combined_plots) print(plot)
# dev.off()
# cat("✅ QC plots saved.\n")
# # Save individual objects
# saveRDS(seurat_objects, file.path(rds_output_dir, "1_qc_sctransform_20250616.rds"))
# cat("✅ QC'ed Seurat objects saved.\n")


# # ------------------------- #
# # JackStraw Analysis: Print Significant PCs Before Plotting
# # ------------------------- #
# jackstraw_plots <- list()

# for (animal in names(seurat_objects)) {
#   cat("📊 Running JackStraw (RNA assay) for", animal, "\n")

#   obj <- seurat_objects[[animal]]

#   # Use raw RNA assay
#   DefaultAssay(obj) <- "RNA"
#   obj <- NormalizeData(obj)
#   obj <- FindVariableFeatures(obj)
#   obj <- ScaleData(obj)
#   obj <- RunPCA(obj)

#   # Determine valid dimensions
#   available_pcs <- ncol(obj@reductions$pca@cell.embeddings)
#   use_dims <- 1:min(30, available_pcs)

#   # Run JackStraw
#   obj <- JackStraw(obj, num.replicate = 100)
#   obj <- ScoreJackStraw(obj, dims = use_dims)

#   # Print significant PCs
#   jack_df <- obj@reductions$pca@jackstraw$overall.p.values
#   colnames(jack_df) <- c("PC", "p_value")
#   sig_pcs <- jack_df %>% dplyr::filter(p_value < 0.05)

#   cat("📌 Significant PCs for", animal, "(p < 0.05):\n")
#   print(sig_pcs)

#   # Store plot
#   p <- JackStrawPlot(obj, dims = use_dims) +
#     ggtitle(paste("JackStraw Plot -", animal)) +
#     theme(plot.title = element_text(hjust = 0.5))
#   jackstraw_plots[[animal]] <- p

#   # Save back
#   seurat_objects[[animal]] <- obj
# }

# # Save plots to one file
# pdf(file.path(pdf_output_dir, "jackstraw_totalNK_animals_20250616.pdf"), width = 10, height = 6 * length(jackstraw_plots))
# for (plot in jackstraw_plots) print(plot)
# dev.off()
# cat("✅ JackStraw plots saved.\n")


# # ------------------------- #
# # Elbow Plot for Each Animal
# # ------------------------- #
# cat("📊 Generating elbow plots...\n")
# elbow_plots <- list()

# for (animal in names(seurat_objects)) {
#   seurat_obj <- seurat_objects[[animal]]
#   DefaultAssay(seurat_obj) <- "SCT"
#   seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
#   elbow_plots[[animal]] <- ElbowPlot(seurat_obj, ndims = 30) + 
#     ggtitle(paste("Elbow Plot -", animal)) +
#     theme(plot.title = element_text(hjust = 0.5, size = 12))
# }

# # Combine and save elbow plots
# pdf(file.path(pdf_output_dir, "elbow_plots_totalNK_each_animals_20250616.pdf"), width = 10, height = 8, onefile = TRUE)
# print(wrap_plots(elbow_plots, ncol = 2))
# dev.off()
# cat("✅ Elbow plots saved to PDF.\n")

# # ------------------------- #
# # UMAP & Clustering Per Animal (combined PDF per dim+res)
# # ------------------------- #
# dims_list <- c(10, 15, 20, 25, 30)
# resolutions <- c(0.25, 0.5, 0.75)

# for (dim in dims_list) {
#   dims_to_use <- 1:dim
#   for (res in resolutions) {
#     plot_list <- list()

#     for (animal in names(seurat_objects)) {
#       cat("📌 Processing:", animal, "dims:", dim, "res:", res, "\n")
#       seurat_obj <- seurat_objects[[animal]]
#       DefaultAssay(seurat_obj) <- "SCT"
#       seurat_obj <- ScaleData(seurat_obj)
#       seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
#       seurat_obj <- RunUMAP(seurat_obj, dims = dims_to_use, reduction.name = paste0("umap_dims", dim))
#       seurat_obj <- FindNeighbors(seurat_obj, dims = dims_to_use)
#       seurat_obj <- FindClusters(seurat_obj, resolution = res)

#       # Store updated object back
#       seurat_objects[[animal]] <- seurat_obj

#       # Generate UMAP plot
#       p <- DimPlot(
#         seurat_obj,
#         reduction = paste0("umap_dims", dim),
#         group.by = "seurat_clusters",
#         label = TRUE,
#         pt.size = 0.3,
#         repel = TRUE
#       ) +
#         ggtitle(paste("UMAP:", animal, "dims =", dim, "res =", res)) +
#         theme(plot.title = element_text(hjust = 0.5)) +
#         scale_color_viridis_d(option = "turbo")

#         plot_list[[animal]] <- p
#     }

#     # Combine plots for current dims+res
#     combined_plot <- wrap_plots(plot_list, ncol = 2)
#     pdf_name <- file.path(pdf_output_dir, paste0("combined_umap_dims", dim, "_res", res, "_20250616.pdf"))
#     ggsave(
#       filename = pdf_name,
#       plot = combined_plot,
#       width = 12,
#       height = 10,
#       dpi = 600
#     )
#     cat("✅ Saved:", pdf_name, "\n")
#   }
# }


# # ------------------------- #
# # Track Cluster Cell Counts
# # ------------------------- #
# cluster_summary <- data.frame()

# for (dim in dims_list) {
#   dims_to_use <- 1:dim
#   for (res in resolutions) {
#     for (animal in names(seurat_objects)) {
#       cat("📊 Counting cells for", animal, "dims =", dim, "res =", res, "\n")

#       obj <- seurat_objects[[animal]]
#       DefaultAssay(obj) <- "SCT"
#       obj <- ScaleData(obj)
#       obj <- RunPCA(obj, verbose = FALSE)
#       obj <- RunUMAP(obj, dims = dims_to_use, reduction.name = paste0("umap_dims", dim))
#       obj <- FindNeighbors(obj, dims = dims_to_use)
#       obj <- FindClusters(obj, resolution = res)

#       seurat_objects[[animal]] <- obj  # update stored object

#       cluster_table <- table(Idents(obj))
#       df <- data.frame(
#         animal = animal,
#         dims = dim,
#         resolution = res,
#         cluster = as.integer(names(cluster_table)),
#         cell_count = as.vector(cluster_table)
#       )

#       cluster_summary <- rbind(cluster_summary, df)
#     }
#   }
# }

# # Save as CSV
# csv_path <- file.path(tsv_output_dir, "cluster_cell_counts_summary.csv")
# write.csv(cluster_summary, csv_path, row.names = FALSE)
# cat("✅ Cluster cell count summary saved to:", csv_path, "\n")


# # ------------------------- #
# # Merge and Integrate All Animals
# # ------------------------- #
# cat("🔗 Integrating all animals...\n")
# features <- SelectIntegrationFeatures(object.list = seurat_objects, nfeatures = 3000)
# seurat_objects <- PrepSCTIntegration(object.list = seurat_objects, anchor.features = features)
# anchors <- FindIntegrationAnchors(object.list = seurat_objects, normalization.method = "SCT", anchor.features = features)
# merged_seurat_obj <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
# # Save merged object
# merged_rds <- file.path(rds_output_dir, "merged_totalNK_20250616.rds")
# saveRDS(merged_seurat_obj, merged_rds)
# cat("✅ Merged Seurat object saved to", merged_rds, "\n")


# ------------------------- #
# Elbow Plot for Merged Seurat Object
# ------------------------- #

# Load merged object if not in memory
merged_rds <- file.path(rds_output_dir, "merged_totalNK_20250616.rds")
merged_seurat_obj <- readRDS(merged_rds)

# Ensure PCA has been run (optional if already present)
if (!"pca" %in% names(merged_seurat_obj@reductions)) {
  DefaultAssay(merged_seurat_obj) <- "integrated"
  merged_seurat_obj <- ScaleData(merged_seurat_obj, verbose = FALSE)
  merged_seurat_obj <- RunPCA(merged_seurat_obj, npcs = 30, verbose = FALSE)
}

# Create Elbow Plot
elbow_plot <- ElbowPlot(merged_seurat_obj, ndims = 30) +
  ggtitle("Elbow Plot - Merged Total NK") +
  theme(plot.title = element_text(hjust = 0.5))

# Save to PDF
ggsave(
  filename = file.path(pdf_output_dir, "elbow_plot_merged_totalNK_20250616.pdf"),
  plot = elbow_plot,
  width = 8,
  height = 6,
  dpi = 600
)
cat("✅ Elbow plot saved to:", file.path(pdf_output_dir, "elbow_plot_merged_totalNK_20250616.pdf"), "\n")


# # ------------------------- #
# # UMAP & Clustering on Merged Data (Combined PDF per dims)
# # ------------------------- #
# DefaultAssay(merged_seurat_obj) <- "integrated"
# merged_seurat_obj <- ScaleData(merged_seurat_obj, verbose = FALSE)
# merged_seurat_obj <- RunPCA(merged_seurat_obj, npcs = max(dims_list))  # ensure enough PCs

# for (dim in dims_list) {
#   dims_to_use <- 1:dim
#   merged_seurat_obj <- RunUMAP(merged_seurat_obj, dims = dims_to_use, reduction.name = paste0("umap_dims", dim))

#   umap_plots <- list()
#   for (res in resolutions) {
#     merged_seurat_obj <- FindNeighbors(merged_seurat_obj, dims = dims_to_use)
#     merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = res)

#     p <- DimPlot(
#       merged_seurat_obj,
#       reduction = paste0("umap_dims", dim),
#       group.by = "seurat_clusters",
#       label = TRUE,
#       pt.size = 0.3,
#       repel = TRUE
#     ) +
#       ggtitle(paste("Merged UMAP - dims =", dim, "res =", res)) +
#       theme(plot.title = element_text(hjust = 0.5)) +
#       scale_color_viridis_d(option = "turbo")

#     umap_plots[[paste0("res_", res)]] <- p
#   }

#   # Combine all resolution plots for this dim into one PDF
#   combined_pdf <- wrap_plots(umap_plots, ncol = 2)
#   pdf_path <- file.path(pdf_output_dir, paste0("merged_umap_dims", dim, "_multiRes_20250616.pdf"))

#   ggsave(
#     filename = pdf_path,
#     plot = combined_pdf,
#     width = 12,
#     height = 4 * ceiling(length(umap_plots) / 2),
#     dpi = 600
#   )
#   cat("✅ Merged multi-resolution UMAP saved:", pdf_path, "\n")
# }