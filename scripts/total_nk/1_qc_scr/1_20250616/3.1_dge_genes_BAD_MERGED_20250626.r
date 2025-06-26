# --------------------------- #
# Load Libraries
# --------------------------- #
library(Seurat)
library(dplyr)
# --------------------------- #
# Paths
# --------------------------- #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
rds_output_dir <- file.path(output_dir, "rds")
tsv_output_dir <- file.path(output_dir, "tsv")
# Create output directories
dir.create(rds_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tsv_output_dir, showWarnings = FALSE, recursive = TRUE)
for (type in c("individual", "merged")) {
 for (method in c("SCT", "LogNormalize")) {
   dir.create(file.path(tsv_output_dir, type, method), recursive = TRUE, showWarnings = FALSE)
 }
}
# --------------------------- #
# Load Seurat Objects
# --------------------------- #
cat("Loading Seurat objects...\n")
indie_rds <- readRDS(file.path(rds_output_dir, "individual_seurat_analysis_20250616.rds"))
merged_rds <- readRDS(file.path(rds_output_dir, "merged_seurat_analysis_20250616.rds"))
# --------------------------- #
# Settings
# --------------------------- #
dims_list <- c(10, 15, 20, 25, 30)
resolutions <- c(0.25, 0.5, 0.75)
set.seed(2025)
# --------------------------- #
# DGE Function
# --------------------------- #
perform_dge <- function(seurat_obj, sample_id, method, dims, res, output_type) {
 cat("Running DGE:", sample_id, "|", method, "| dims =", dims, "| res =", res, "\n")
 cluster_col <- paste0(ifelse(method == "SCT", "SCT_snn_res.", "RNA_snn_res."), res)
 assay <- ifelse(method == "SCT", "SCT", "RNA")
 if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
   message("Cluster column not found:", cluster_col)
   return(NULL)
 }
 tryCatch({
   DefaultAssay(seurat_obj) <- assay
   if (assay == "RNA" && "layers" %in% names(seurat_obj@assays$RNA)) {
     seurat_obj <- JoinLayers(seurat_obj)
   }
   Idents(seurat_obj) <- seurat_obj[[cluster_col]][, 1]
   markers <- FindAllMarkers(
     seurat_obj,
     logfc.threshold = 0.25,
     min.pct = 0.1,
     only.pos = TRUE
   )
   if (!is.null(markers) && "cluster" %in% colnames(markers)) {
     markers$animal <- sample_id
     markers$method <- method
     markers$dims <- dims
     markers$resolution <- res
     out_path <- file.path(tsv_output_dir, output_type, method)
     file_out <- file.path(out_path, paste0("markers_", sample_id, "_", method, "_dims", dims, "_res", res, ".csv"))
     write.csv(markers, file_out, row.names = FALSE)
     cat("Saved:", file_out, "\n")
     return(markers)
   } else {
     message("No DEGs found for:", sample_id, "|", method, "| dims =", dims, "| res =", res)
     return(NULL)
   }
 }, error = function(e) {
   message("DGE failed for:", sample_id, "|", method, "| dims =", dims, "| res =", res)
   return(NULL)
 })
}
# --------------------------- #
# Run DGE for Individual Samples
# --------------------------- #
cat("Running DGE for individual samples...\n")
individual_dge_results <- list()
for (animal in names(indie_rds)) {
 for (method in c("SCT", "LogNormalize")) {
   for (dims in dims_list) {
     for (res in resolutions) {
       result <- perform_dge(indie_rds[[animal]], animal, method, dims, res, "individual")
       if (!is.null(result)) {
         individual_dge_results[[paste(animal, method, dims, res, sep = "_")]] <- result
       }
     }
   }
 }
}
saveRDS(individual_dge_results, file = file.path(rds_output_dir, "individual_dge_results_20250616.rds"))
cat("Saved RDS: individual_dge_results_20250616.rds\n")
# --------------------------- #
# Run DGE for Merged Sample
# --------------------------- #
cat("Running DGE for merged sample...\n")
merged_dge_results <- list()
for (method in c("SCT", "LogNormalize")) {
 merged_obj <- merged_rds[[method]]
 for (dims in dims_list) {
   for (res in resolutions) {
     result <- perform_dge(merged_obj, "merged", method, dims, res, "merged")
     if (!is.null(result)) {
       merged_dge_results[[paste("merged", method, dims, res, sep = "_")]] <- result
     }
   }
 }
}
saveRDS(merged_dge_results, file = file.path(rds_output_dir, "merged_dge_results_20250616.rds"))
cat("Saved RDS: merged_dge_results_20250616.rds\n")