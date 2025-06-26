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
dir.create(tsv_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rds_output_dir, showWarnings = FALSE, recursive = TRUE)
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
# merged_rds <- readRDS(file.path(rds_output_dir, "merged_seurat_analysis_20250616.rds"))
# --------------------------- #
# Load Cluster Summaries
# --------------------------- #
cat("Loading cluster summaries...\n")
indie_summary <- read.csv(file.path(tsv_output_dir, "indie_cluster_summary_sct_vs_log.csv"), stringsAsFactors = FALSE)
# merged_summary <- read.csv(file.path(tsv_output_dir, "merged_cluster_summary_sct_vs_log.csv"), stringsAsFactors = FALSE)

# Check if cluster columns exist
# colnames(merged_rds[["SCT"]]@meta.data)


# --------------------------- #
# Settings
# --------------------------- #
dims_list <- c(10, 15, 20, 25, 30)
resolutions <- c(0.25, 0.5, 0.75)
set.seed(2025)
dge_summary <- data.frame()
# --------------------------- #
# DGE Function
# --------------------------- #
perform_dge <- function(seurat_obj, sample_id, method, dims, res, output_type) {
 cat("Running DGE:", sample_id, "|", method, "| dims =", dims, "| res =", res, "\n")


 # Use correct clustering metadata column
  if (sample_id == "merged") {
  cluster_col <- paste0("integrated_snn_res.", res)
  } else {
  cluster_col <- paste0(ifelse(method == "SCT", "SCT_snn_res.", "RNA_snn_res."), res)
  }

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
   Idents(seurat_obj) <- factor(seurat_obj[[cluster_col]][, 1])
   markers <- FindAllMarkers(
     seurat_obj,
     logfc.threshold = 0.25,
     min.pct = 0.1,
     only.pos = TRUE
   )
   if (!is.null(markers) && "cluster" %in% colnames(markers)) {
     # Add metadata for tracking
     markers$animal <- sample_id
     markers$method <- method
     markers$dims <- dims
     markers$resolution <- res
     # Count DEGs per cluster
     deg_counts <- markers %>%
       count(cluster, name = "n_DEGs") %>%
       mutate(
         animal = sample_id,
         method = method,
         dims = dims,
         resolution = res
       )
     # Save markers to structured folder
     out_dir <- file.path(tsv_output_dir, output_type, method)
     file_out <- file.path(out_dir, paste0("new_markers_", sample_id, "_", method, "_dims", dims, "_res", res, ".csv"))
     write.csv(markers, file_out, row.names = FALSE)
     cat("Saved:", file_out, "\n")
     return(deg_counts)
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
dge_summary_individual <- data.frame()
for (animal in names(indie_rds)) {
 for (method in c("SCT", "LogNormalize")) {
   for (dims in dims_list) {
     for (res in resolutions) {
       result <- perform_dge(indie_rds[[animal]], animal, method, dims, res, "individual")
       if (!is.null(result)) {
         dge_summary_individual <- rbind(dge_summary_individual, result)
       }
     }
   }
 }
}
# --------------------------- #
# Join with Individual Summary
# --------------------------- #
cat("Merging with individual summary...\n")
indie_summary$cluster <- as.character(indie_summary$cluster)
dge_summary_individual$cluster <- as.character(dge_summary_individual$cluster)
indie_augmented <- indie_summary %>%
 left_join(dge_summary_individual, by = c("animal", "method", "dims", "resolution", "cluster"))
write.csv(indie_augmented, file.path(tsv_output_dir, "indie_cluster_summary_sct_vs_log_dge.csv"), row.names = FALSE)
saveRDS(indie_augmented, file.path(rds_output_dir, "indie_cluster_summary_sct_vs_log_dge.rds"))
cat("Saved: indie_cluster_summary_sct_vs_log_dge\n")



# # --------------------------- #
# # Run DGE for Merged Object
# # --------------------------- #
# cat("Running DGE for merged sample...\n")
# dge_summary_merged <- data.frame()
# for (method in c("SCT")) {
#  merged_obj <- merged_rds[[method]]
#  for (dims in dims_list) {
#    for (res in resolutions) {
#      result <- perform_dge(merged_obj, "merged", method, dims, res, "merged")
#      if (!is.null(result)) {
#        dge_summary_merged <- rbind(dge_summary_merged, result)
#      }
#    }
#  }
# }
# # --------------------------- #
# # Join with Merged Summary
# # --------------------------- #
# cat("Merging with merged summary...\n")
# merged_summary$cluster <- as.character(merged_summary$cluster)
# dge_summary_merged$cluster <- as.character(dge_summary_merged$cluster)
# merged_augmented <- merged_summary %>%
#  left_join(dge_summary_merged, by = c("animal", "method", "dims", "resolution", "cluster"))
# write.csv(merged_augmented, file.path(tsv_output_dir, "merged_cluster_summary_sct_vs_log_dge.csv"), row.names = FALSE)
# saveRDS(merged_augmented, file.path(rds_output_dir, "merged_cluster_summary_sct_vs_log_dge.rds"))
# cat("Saved: merged_cluster_summary_sct_vs_log_dge\n")