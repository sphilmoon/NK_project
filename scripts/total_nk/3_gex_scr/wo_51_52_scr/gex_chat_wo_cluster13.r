# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(dplyr)
# ------------------------- #
# Define Paths and Configs
# ------------------------- #
rds_file <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/rds/integrated_data_dims25_res0.3_genecounts.rds"
output_dir <- "/home/outputs/totalNK_outputs/3_gex/wo_51_52"

pdf_dir <- file.path(output_dir, "pdf")
dge_dir <- file.path(output_dir, "dge")
rds_dir <- file.path(output_dir, "rds")
# ------------------------- #

if (!dir.exists(dge_dir)) {
  cat("📂 Creating directory for DGE results:", dge_dir, "\n")
  # Create the directory if it does not exist
 dir.create(dge_dir, recursive = TRUE)
}
if (!dir.exists(rds_dir)) {
  cat("📂 Creating directory for RDS files:", rds_dir, "\n")
  # Create the directory if it does not exist
 dir.create(rds_dir, recursive = TRUE)
}

# ------------------------- #
# Load Seurat Object
# ------------------------- #
seurat_obj <- readRDS(rds_file)
cat("✅ Loaded Seurat object from", rds_file, "\n")
# ------------------------- #
# Set Assay (Optional but Recommended)
# ------------------------- #
if ("SCT" %in% Assays(seurat_obj)) {
 DefaultAssay(seurat_obj) <- "SCT"
} else {
 DefaultAssay(seurat_obj) <- "RNA"
}
cat("✅ Default assay set to", DefaultAssay(seurat_obj), "\n")

# ------------------------- #
# Check if 'seurat_clusters' exists
colnames(seurat_obj@meta.data)
head(seurat_obj@meta.data$seurat_clusters)
# table(seurat_obj$orig.ident)
# table(seurat_obj$sample_id)

grep("res", colnames(seurat_obj@meta.data), value = TRUE)
seurat_obj$seurat_clusters <- seurat_obj$integrated_snn_res.0.3

# ------------------------- #
# Set Identity to Clustering Resolution
# ------------------------- #

# Ensure clustering is stored as factor
seurat_obj$seurat_clusters <- as.factor(seurat_obj$integrated_snn_res.0.3)
Idents(seurat_obj) <- "seurat_clusters"
cat("✅ Idents set to 'seurat_clusters'.\n")

# ------------------------- #
# Display Cluster Information
# ------------------------- #
clusters <- unique(Idents(seurat_obj))
cat("📊 Total clusters identified:", length(clusters), "\n")
cat("📋 Cluster IDs:", paste(sort(as.character(clusters)), collapse = ", "), "\n")

# ------------------------- #
# Exclude Cluster "13"
# ------------------------- #
# Subset to exclude cluster "13"
if ("13" %in% clusters) {
 seurat_obj <- subset(seurat_obj, idents = setdiff(clusters, "13"))
 cat("❌ Cluster '13' excluded from Seurat object.\n")
}




# ------------------------- #
# Prepare SCT Model for DGE
# ------------------------- #
cat("🔧 Preparing SCT model offsets for differential expression...\n")
seurat_obj <- PrepSCTFindMarkers(seurat_obj, split.by = "orig.ident")

# ------------------------- #
# Run Differential Expression
# ------------------------- #
cat("🔍 Running FindAllMarkers() across remaining clusters...\n")
markers_all <- FindAllMarkers(
 seurat_obj,
 only.pos = TRUE,
 min.pct = 0.25,
 logfc.threshold = 0.25,
 test.use = "wilcox"
)
# Check and save
if (nrow(markers_all) == 0) {
 cat("⚠️ No DEGs were identified after filtering. Skipping file export.\n")
} else {
 saveRDS(markers_all, file.path(rds_dir, "markers_all_wo_cluster13.rds"))
 write.csv(markers_all, file.path(dge_dir, "differential_expression_wo_cluster13.csv"), row.names = FALSE)
 cat("✅ DEG results saved.\n")

 # ------------------------- #
 # Save Top 20 Markers per Cluster
 # ------------------------- #
 top20 <- markers_all %>%
   group_by(cluster) %>%
   slice_max(n = 20, order_by = avg_log2FC)
 write.csv(top20, file.path(dge_dir, "top20_markers_per_cluster_wo_cluster13.csv"), row.names = FALSE)
 cat("✅ Top 20 markers per cluster saved.\n")
}


