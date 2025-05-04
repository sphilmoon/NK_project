# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

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

# Ensure layers are joined (Seurat v5)
# seurat_obj <- JoinLayers(seurat_obj, assay = DefaultAssay(seurat_obj))

# ------------------------- #
# Define Gene to Plot
# ------------------------- #
gene <- "NCR1"

# Check if gene exists
expr_data_all <- GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = "data")
if (!gene %in% rownames(expr_data_all)) {
  stop("❌ Gene ", gene, " not found in the Seurat object.")
}
cat("✅ Gene", gene, "found in the Seurat object\n")

# ------------------------- #
# Extract Expression Data Per Cluster
# ------------------------- #
cat("🔍 Extracting", gene, "expression per cluster...\n")

# Fetch expression data and cluster identities
expression_data <- FetchData(seurat_obj, vars = c("seurat_clusters", gene))
colnames(expression_data) <- c("cluster", "expression")

# ------------------------- #
# Create DimPlot-Style Plot
# ------------------------- #
cat("🎨 Creating DimPlot-style plot for", gene, "expression per cluster...\n")

plot <- ggplot(expression_data, aes(x = cluster, y = expression)) +
  geom_violin(trim = FALSE, fill = "lightgray", alpha = 0.5) +  # Violin plot for distribution
  geom_jitter(size = 0.5, alpha = 0.6, width = 0.2, color = "black") +  # Jitter points for individual cells
  scale_y_continuous(name = paste(gene, "Expression Level")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle(paste("NCR1 Expression Per Cluster (dims25, res 0.3)"))

# ------------------------- #
# Save Plot
# ------------------------- #
output_file <- file.path(output_dir, "NCR1_DimPlotStyle_expressionPerCluster_d25_res0.3.pdf")
ggsave(filename = output_file, plot = plot, width = 10, height = 6, dpi = 600, bg = "transparent")
cat("✅ DimPlot-style plot saved to", output_file, "\n")

cat("🎉 Plot generation complete. Outputs saved in:\n", output_dir, "\n")