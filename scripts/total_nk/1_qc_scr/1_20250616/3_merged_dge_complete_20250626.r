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
dir.create(file.path(tsv_output_dir, "merged", "SCT"), recursive = TRUE, showWarnings = FALSE)

# --------------------------- #
# Load Merged Seurat Object
# --------------------------- #
cat("📦 Loading merged Seurat object...\n")
merged_rds <- readRDS(file.path(rds_output_dir, "merged_seurat_analysis_20250616.rds"))

# --------------------------- #
# Load Cluster Summary (Optional)
# --------------------------- #
cat("📄 Loading cluster summary for merged...\n")
summary_file <- file.path(tsv_output_dir, "merged_cluster_summary_sct_vs_log.csv")
merged_summary <- if (file.exists(summary_file)) {
  read.csv(summary_file, stringsAsFactors = FALSE)
} else {
  NULL
}

# --------------------------- #
# Settings
# --------------------------- #
dims_list <- c(10, 15, 20, 25, 30)
resolutions <- c(0.25, 0.5, 0.75)
set.seed(2025)
dge_summary_merged <- data.frame()

# --------------------------- #
# DGE Function
# --------------------------- #
perform_dge <- function(seurat_obj, cluster_col, sample_id, method, dims, res, output_type) {
  cat("🔬 Running DGE:", sample_id, "|", method, "| dims =", dims, "| res =", res, "\n")
  
  assay <- ifelse(method == "SCT", "SCT", "RNA")
  if (!(cluster_col %in% colnames(seurat_obj@meta.data))) {
    message("❌ Cluster column not found:", cluster_col)
    return(NULL)
  }
  
  tryCatch({
    DefaultAssay(seurat_obj) <- assay
    if (assay == "RNA" && "layers" %in% names(seurat_obj@assays$RNA)) {
      seurat_obj <- JoinLayers(seurat_obj)
    }
    if (assay == "SCT") {
      seurat_obj <- PrepSCTFindMarkers(seurat_obj)
    }

    Idents(seurat_obj) <- seurat_obj[[cluster_col]][, 1]

    markers <- FindAllMarkers(
      seurat_obj,
      logfc.threshold = 0.25,
      min.pct = 0.1,
      only.pos = TRUE
    )

    if (!is.null(markers) && "cluster" %in% colnames(markers)) {
      # Add metadata
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

      # Save markers
      out_dir <- file.path(tsv_output_dir, output_type, method)
      file_out <- file.path(out_dir, paste0("new_markers_", sample_id, "_", method, "_dims", dims, "_res", res, ".csv"))
      write.csv(markers, file_out, row.names = FALSE)
      cat("✅ Saved:", file_out, "\n")
      return(deg_counts)
    } else {
      message("⚠️ No DEGs found for:", sample_id, "|", method, "| dims =", dims, "| res =", res)
      return(NULL)
    }
  }, error = function(e) {
    message("❌ DGE failed for:", sample_id, "|", method, "| dims =", dims, "| res =", res)
    return(NULL)
  })
}

# --------------------------- #
# Run DGE for Merged Object (SCT only)
# --------------------------- #
cat("🚀 Running DGE for merged sample using SCT...\n")

method <- "SCT"
merged_obj <- merged_rds[[method]]

for (dims in dims_list) {
  for (res in resolutions) {
    # Generate unique clustering column name
    cluster_col <- paste0("integrated_snn_res.d", dims, "_r", res)

    # Recompute neighbors and clustering
    merged_obj <- FindNeighbors(merged_obj, dims = 1:dims, verbose = FALSE)
    merged_obj <- FindClusters(merged_obj, resolution = res, verbose = FALSE)

    # Save to metadata
    merged_obj[[cluster_col]] <- merged_obj$seurat_clusters

    # Run DEG
    result <- perform_dge(merged_obj, cluster_col, "merged", method, dims, res, "merged")
    if (!is.null(result)) {
      dge_summary_merged <- rbind(dge_summary_merged, result)
    }
  }
}

# --------------------------- #
# Join DGE with Summary Metadata
# --------------------------- #
if (!is.null(merged_summary)) {
  cat("🧾 Joining DGE results with merged summary...\n")
  merged_summary$cluster <- as.character(merged_summary$cluster)
  dge_summary_merged$cluster <- as.character(dge_summary_merged$cluster)
  
  merged_augmented <- merged_summary %>%
    left_join(dge_summary_merged, by = c("animal", "method", "dims", "resolution", "cluster"))

  write.csv(merged_augmented, file.path(tsv_output_dir, "merged_cluster_summary_sct_vs_log_dge.csv"), row.names = FALSE)
  saveRDS(merged_augmented, file.path(rds_output_dir, "merged_cluster_summary_sct_vs_log_dge.rds"))
  cat("💾 Saved merged DGE summary with metadata.\n")
} else {
  write.csv(dge_summary_merged, file.path(tsv_output_dir, "merged_dge_counts_summary.csv"), row.names = FALSE)
  cat("💾 Saved merged DGE summary without metadata.\n")
}