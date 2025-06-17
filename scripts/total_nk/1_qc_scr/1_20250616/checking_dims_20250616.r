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


# ------------------------- #
# Define Input Directories
# ------------------------- #
totalNK_h5_dir <- "/home/rawdata/totalNK"
animal_names <- c("animal25", "animal26", "animal27", "animal28")


# ------------------------- #
# Build File Paths
# ------------------------- #
build_h5_paths <- function(base_dir, condition_prefix) {
 fnames <- paste0(tolower(animal_names), "_", condition_prefix, "_filtered_feature_bc_matrix.h5")
 full_paths <- file.path(base_dir, fnames)
 missing <- !file.exists(full_paths)
 if (any(missing)) {
   cat("❌ Missing files:\n")
   print(full_paths[missing])
   stop("Some .h5 files are missing. Please check the paths or filenames.")
 }
 names(full_paths) <- animal_names
 return(full_paths)
}
totalNK_file_paths <- build_h5_paths(totalNK_h5_dir, "totalnk")


# ------------------------- #
# Load and QC Each Animal
# ------------------------- #
seurat_objects <- list()
combined_plots <- list()
perform_qc <- function(seurat_obj, sample_name) {
 qc_params <- list(
   animal27 = list(nFeature_min = 200, nFeature_max = 3000, nCount_min = 250, nCount_max = 7500),
   animal28 = list(nFeature_min = 200, nFeature_max = 3000, nCount_min = 250, nCount_max = 10000),
   default  = list(nFeature_min = 200, nFeature_max = 3000, nCount_min = 500, nCount_max = 10000)
 )
 params <- if (sample_name %in% names(qc_params)) qc_params[[sample_name]] else qc_params$default
 seurat_obj <- subset(
   seurat_obj,
   subset = nFeature_RNA > params$nFeature_min & nFeature_RNA < params$nFeature_max &
            nCount_RNA > params$nCount_min & nCount_RNA < params$nCount_max
 )
 return(seurat_obj)
}
for (i in seq_along(totalNK_file_paths)) {
 sample <- names(totalNK_file_paths)[i]
 filepath <- totalNK_file_paths[[i]]
 cat("📥 Reading:", sample, "\n")
 seurat_obj <- Read10X_h5(filepath) %>%
   CreateSeuratObject(project = sample, min.cells = 3, min.features = 200)
 seurat_obj$animal <- sample
 # Add mitochondrial percentage
 seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
 # QC filter
 seurat_obj <- perform_qc(seurat_obj, sample)
 # Normalize
 seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
 # Save object
 seurat_objects[[sample]] <- seurat_obj
 # Generate QC plots
 p1 <- VlnPlot(seurat_obj, features = "nFeature_RNA", pt.size = 0) + ggtitle("Gene counts")
 p2 <- VlnPlot(seurat_obj, features = "nCount_RNA", pt.size = 0) + ggtitle("UMI counts")
 combined_plot <- p1 + p2 + plot_layout(ncol = 2) + plot_annotation(title = paste("QC for", sample))
 combined_plots[[sample]] <- combined_plot
 cat("✅ QC and normalization complete for", sample, "\n")
}
# Save all QC plots
pdf(file.path(pdf_output_dir, "qc_totalNK_animals_20250616.pdf"), width = 10, height = 6 * length(combined_plots))
for (plot in combined_plots) print(plot)
dev.off()
cat("✅ QC plots saved.\n")
# Save individual objects
saveRDS(seurat_objects, file.path(rds_output_dir, "1_qc_sctransform_20250616.rds"))
cat("✅ QC'ed Seurat objects saved.\n")


# ------------------------- #
# Elbow Plot for Each Animal
# ------------------------- #
cat("📊 Generating elbow plots...\n")
elbow_plots <- list()

for (animal in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[animal]]
  DefaultAssay(seurat_obj) <- "SCT"
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  elbow_plots[[animal]] <- ElbowPlot(seurat_obj, ndims = 30) + 
    ggtitle(paste("Elbow Plot -", animal)) +
    theme(plot.title = element_text(hjust = 0.5, size = 12))
}

# Combine and save elbow plots
pdf(file.path(pdf_output_dir, "elbow_plots_totalNK_each_animals_20250616.pdf"), width = 10, height = 8, onefile = TRUE)
print(wrap_plots(elbow_plots, ncol = 2))
dev.off()
cat("✅ Elbow plots saved to PDF.\n")


# ------------------------- #
# UMAP & Clustering Per Animal
# ------------------------- #
dims_list <- c(15, 20, 25, 30)
resolutions <- c(0.25, 0.5, 0.75)
for (animal in names(seurat_objects)) {
 cat("📌 Running UMAP + clustering for", animal, "\n")
 seurat_obj <- seurat_objects[[animal]]
 DefaultAssay(seurat_obj) <- "SCT"
 seurat_obj <- ScaleData(seurat_obj)
 seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
 for (dim in dims_list) {
   dims_to_use <- 1:dim
   seurat_obj <- RunUMAP(seurat_obj, dims = dims_to_use, reduction.name = paste0("umap_dims", dim))
   for (res in resolutions) {
     seurat_obj <- FindNeighbors(seurat_obj, dims = dims_to_use)
     seurat_obj <- FindClusters(seurat_obj, resolution = res)
     umap_plot <- DimPlot(
       seurat_obj,
       reduction = paste0("umap_dims", dim),
       group.by = "seurat_clusters",
       label = TRUE,
       pt.size = 0.3,
       repel = TRUE
     ) +
       ggtitle(paste("UMAP for", animal, "(dims =", dim, ", res =", res, ")")) +
       theme(plot.title = element_text(hjust = 0.5)) +
       scale_color_brewer(palette = "Set2")
     ggsave(
       filename = file.path(pdf_output_dir, paste0("umap_", animal, "_dims", dim, "_res", res, "_20250616.pdf")),
       plot = umap_plot,
       width = 10,
       height = 8,
       dpi = 600
     )
     cat("✅ UMAP saved for", animal, "- dims:", dim, "res:", res, "\n")
   }
 }
}


# ------------------------- #
# Merge and Integrate All Animals
# ------------------------- #
cat("🔗 Integrating all animals...\n")
features <- SelectIntegrationFeatures(object.list = seurat_objects, nfeatures = 3000)
seurat_objects <- PrepSCTIntegration(object.list = seurat_objects, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seurat_objects, normalization.method = "SCT", anchor.features = features)
merged_seurat_obj <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
# Save merged object
merged_rds <- file.path(rds_output_dir, "merged_totalNK_20250616.rds")
saveRDS(merged_seurat_obj, merged_rds)
cat("✅ Merged Seurat object saved to", merged_rds, "\n")


# ------------------------- #
# UMAP & Clustering on Merged Data
# ------------------------- #
DefaultAssay(merged_seurat_obj) <- "integrated"
merged_seurat_obj <- ScaleData(merged_seurat_obj, verbose = FALSE)
merged_seurat_obj <- RunPCA(merged_seurat_obj, npcs = 30)
for (dim in dims_list) {
 dims_to_use <- 1:dim
 merged_seurat_obj <- RunUMAP(merged_seurat_obj, dims = dims_to_use, reduction.name = paste0("umap_dims", dim))
 for (res in resolutions) {
   merged_seurat_obj <- FindNeighbors(merged_seurat_obj, dims = dims_to_use)
   merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = res)
   umap_plot <- DimPlot(
     merged_seurat_obj,
     reduction = paste0("umap_dims", dim),
     group.by = "seurat_clusters",
     label = TRUE,
     pt.size = 0.3,
     repel = TRUE
   ) +
     ggtitle(paste("Merged UMAP (dims =", dim, ", res =", res, ")")) +
     theme(plot.title = element_text(hjust = 0.5)) +
     scale_color_brewer(palette = "Set2")
   ggsave(
     filename = file.path(pdf_output_dir, paste0("merged_umap_dims", dim, "_res", res, "_20250616.pdf")),
     plot = umap_plot,
     width = 10,
     height = 8,
     dpi = 600
   )
   cat("✅ Merged UMAP saved for dims =", dim, "res =", res, "\n")
 }
}