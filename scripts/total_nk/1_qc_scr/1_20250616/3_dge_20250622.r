# --------------------------- #
# Setup
# --------------------------- #
library(Seurat)
library(dplyr)

# Paths
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
rds_output_dir <- file.path(output_dir, "rds")
tsv_output_dir <- file.path(output_dir, "tsv")
dir.create(tsv_output_dir, showWarnings = FALSE, recursive = TRUE)

# Load Seurat Objects
cat("📦 Loading Seurat objects...\n")
indie_rds <- readRDS(file.path(rds_output_dir, "individual_seurat_analysis_20250616.rds"))
merged_rds <- readRDS(file.path(rds_output_dir, "merged_seurat_analysis_20250616.rds"))

# Load existing cluster summaries
cat("📄 Loading existing cluster summaries...\n")
indie_summary <- read.csv(file.path(tsv_output_dir, "indie_cluster_summary_sct_vs_log.csv"))
merged_summary <- read.csv(file.path(tsv_output_dir, "merged_cluster_summary_sct_vs_log.csv"))

# Settings
dims_list <- c(10, 15, 20, 25, 30)
resolutions <- c(0.25, 0.5, 0.75)

# Set seed
set.seed(2025) # 🔒 For reproducibility

# --------------------------- #
# DGE Function
# --------------------------- #
dge_summary <- data.frame()

perform_dge <- function(seurat_obj, sample_id, method, dims, res) {
    cat("🧪 Running DGE:", sample_id, "|", method, "| dims =", dims, "| res =", res, "\n")
  reduction_name <- paste0("umap_", method, "_dims", dims)
  assay <- ifelse(method == "SCT", "SCT", "RNA")
  cluster_col <- paste0(ifelse(method == "SCT", "SCT_snn_res.", "RNA_snn_res."), res)

  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    message("⚠️ ", cluster_col, " not found in ", sample_id)
    return(NULL)
  }

  Idents(seurat_obj) <- seurat_obj[[cluster_col]][, 1]
  DefaultAssay(seurat_obj) <- assay

  # 🔧 Fix for JoinLayers warning (only needed for RNA assay)
  if (assay == "RNA") {
    seurat_obj <- JoinLayers(seurat_obj)
  }

  markers <- tryCatch({
    FindAllMarkers(seurat_obj, logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE)
  }, error = function(e) {
    message("❌ DGE failed for:", sample_id, "(", method, ")", "- dims:", dims, " res:", res)
    return(NULL)
  })

  # 🔐 Check that markers has valid structure
  if (!is.null(markers) && "cluster" %in% colnames(markers)) {
    n_genes <- markers %>% group_by(cluster) %>% tally(name = "n_DEGs")
    n_clusters <- length(unique(Idents(seurat_obj)))
    n_genes$animal <- sample_id
    n_genes$method <- method
    n_genes$dims <- dims
    n_genes$resolution <- res
    n_genes$n_clusters <- n_clusters
    return(n_genes)
  } else {
    message("⚠️ No valid DEGs found for:", sample_id, "|", method, "| dims =", dims, "| res =", res)
    return(NULL)
  }
}

# --------------------------- #
# Run DGE for Individual Samples
# --------------------------- #
cat("🔬 Running DGE for individual samples...\n")
for (animal in names(indie_rds)) {
  for (method in c("SCT", "LogNormalize")) {
    for (dims in dims_list) {
      for (res in resolutions) {
        dge_result <- perform_dge(indie_rds[[animal]], animal, method, dims, res)
        if (!is.null(dge_result)) {
          dge_summary <- rbind(dge_summary, dge_result)
        }
      }
    }
  }
}

# Append DEG counts to summary tables
dge_indie <- dge_summary %>% filter(animal != "merged")

# Join on shared keys
indie_summary_augmented <- indie_summary %>%
  left_join(dge_indie, by = c("animal", "method", "dims", "resolution", "cluster"))


# ✅ Now safe to save to CSV and RDS
# --------------------------- #
write.csv(indie_summary_augmented, file.path(tsv_output_dir, "indie_cluster_summary_sct_vs_log_dge.csv"), row.names = FALSE)
saveRDS(indie_summary_augmented, file.path(tsv_output_dir, "indie_cluster_summary_sct_vs_log_dge.rds"))
cat("💾 Individual samples DGE - RDS files saved for downstream analysis.\n")


# --------------------------- #
# Run DGE for Merged Object
# --------------------------- #
cat("🔗 Running DGE for merged samples...\n")
for (method in c("SCT", "LogNormalize")) {
  merged_obj <- merged_rds[[method]]
  for (dims in dims_list) {
    for (res in resolutions) {
      dge_result <- perform_dge(merged_obj, "merged", method, dims, res)
      if (!is.null(dge_result)) {
        dge_summary <- rbind(dge_summary, dge_result)
      }
    }
  }
}

# --------------------------- #
# Append DEG counts to summary tables
# --------------------------- #
cat("🧾 Merging DEG results with cluster summaries...\n")

# Separate into individual vs merged
dge_merged <- dge_summary %>% filter(animal == "merged")

# Join on shared keys
merged_summary_augmented <- merged_summary %>%
  left_join(dge_merged, by = c("animal", "method", "dims", "resolution", "cluster"))

# --------------------------- #
# ✅ Now safe to save to CSV and RDS
# --------------------------- #
write.csv(merged_summary_augmented, file.path(tsv_output_dir, "merged_cluster_summary_sct_vs_log_dge.csv"), row.names = FALSE)
saveRDS(merged_summary_augmented, file.path(tsv_output_dir, "merged_cluster_summary_sct_vs_log_dge.rds"))
cat("💾 Merged samples DGE - RDS files saved for downstream analysis.\n")