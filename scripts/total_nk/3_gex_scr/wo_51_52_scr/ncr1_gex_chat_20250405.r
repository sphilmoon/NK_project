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
# Define Gene and Check
# ------------------------- #
gene <- "NCR1"
expr_data_all <- GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), slot = "data")
if (!gene %in% rownames(expr_data_all)) {
  stop("❌ Gene ", gene, " not found in the Seurat object.")
}
cat("✅ Gene", gene, "found in the Seurat object\n")

# ------------------------- #
# Extract Data: Clusters + Gene Expression
# ------------------------- #
expression_data <- FetchData(seurat_obj, vars = c("seurat_clusters", gene))
colnames(expression_data) <- c("cluster", "expression")
expression_data$cluster <- factor(expression_data$cluster, levels = sort(unique(expression_data$cluster)))

# ------------------------- #
# Create Plot: Gene Expression by Cluster
# ------------------------- #
plot <- ggplot(expression_data, aes(x = cluster, y = expression)) +
  geom_violin(fill = "lightgray", color = NA, alpha = 0.6, trim = FALSE) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6, color = "black") +
  labs(
    title = paste0("Expression of ", gene, " Across Clusters"),
    x = "Seurat Cluster",
    y = paste(gene, "Expression Level")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

# ------------------------- #
# Save Plot
# ------------------------- #
pdf_path <- file.path(output_dir, paste0(gene, "_Expression_by_Cluster_violin.pdf"))
ggsave(pdf_path, plot, width = 8, height = 6, dpi = 600)
cat("✅ Expression plot saved to", pdf_path, "\n")