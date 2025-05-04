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
# Create Seurat DotPlot
# ------------------------- #
cat("🎨 Creating DotPlot for", gene, "expression across clusters...\n")

dot_plot <- DotPlot(
  seurat_obj,
  features = gene,
  cols = c("white", "red"),  # Color gradient: white (low/zero) to red (high)
  dot.scale = 8  # Adjust dot size scaling
) +
  theme_minimal() +
  scale_y_discrete(labels = gene) +  # Ensure y-axis shows the gene name
  ggtitle("NCR1 Expression Per Cluster (dims25, res 0.3)") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ------------------------- #
# Save Plot
# ------------------------- #
output_file <- file.path(output_dir, "DotPlot_NCR1_perCluster_d25_res0.3.pdf")
ggsave(filename = output_file, plot = dot_plot, width = 8, height = 4, dpi = 600, bg = "transparent")
cat("✅ DotPlot saved to", output_file, "\n")

cat("🎉 Plot generation complete. Outputs saved in:\n", output_dir, "\n")