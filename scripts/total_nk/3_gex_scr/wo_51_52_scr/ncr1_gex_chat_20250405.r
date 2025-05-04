# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(dplyr)
library(ggplot2)

# ------------------------- #
# Define Paths and Configs
# ------------------------- #
rds_file <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/rds/integrated_data_dims25_res0.3_genecounts.rds"
output_dir <- "/home/outputs/totalNK_outputs/3_gex/wo_51_52"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------- #
# Load Seurat Object
# ------------------------- #
seurat_obj <- readRDS(rds_file)
cat("✅ Loaded Seurat object from", rds_file, "\n")

# ------------------------- #
# Set Assay
# ------------------------- #
if ("SCT" %in% Assays(seurat_obj)) {
  DefaultAssay(seurat_obj) <- "SCT"
} else {
  DefaultAssay(seurat_obj) <- "RNA"
}
cat("✅ Default assay set to", DefaultAssay(seurat_obj), "\n")

# ------------------------- #
# Check Gene Exists
# ------------------------- #
gene <- "NCR1"
if (!gene %in% rownames(seurat_obj)) {
  stop("❌ Gene ", gene, " not found in the Seurat object.")
}

# ------------------------- #
# Plot FeaturePlot (Split by Cluster)
# ------------------------- #
cat("🎨 Creating FeaturePlot split by cluster for", gene, "\n")

plot <- FeaturePlot(
  object = seurat_obj,
  features = gene,
  split.by = "seurat_clusters",
  cols = c("white", "lightgray", "red"),
  pt.size = 0.5
) +
  theme(
    strip.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )

# ------------------------- #
# Save PDF
# ------------------------- #
pdf_path <- file.path(output_dir, paste0(gene, "_DimPlotStyle_splitByCluster.pdf"))
ggsave(pdf_path, plot, width = 12, height = 8, dpi = 600, bg = "transparent")
cat("✅ Split-by-cluster FeaturePlot saved to", pdf_path, "\n")