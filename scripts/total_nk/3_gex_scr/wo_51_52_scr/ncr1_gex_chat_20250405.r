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
# Prepare for DE (only for SCT)
# ------------------------- #
if (DefaultAssay(seurat_obj) == "SCT") {
  seurat_obj <- PrepSCTFindMarkers(seurat_obj)
  cat("✅ Ran PrepSCTFindMarkers() for SCT-based DE analysis.\n")
}

# ------------------------- #
# Differential Expression Analysis
# ------------------------- #
cat("🔍 Running FindAllMarkers() across clusters...\n")

markers_all <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# Save full marker results
marker_csv_all <- file.path(output_dir, "differential_expression_all_clusters.csv")
write.csv(markers_all, marker_csv_all, row.names = FALSE)
cat("✅ Full DEG results saved to", marker_csv_all, "\n")

# Save top 20 markers per cluster — if DE results exist
if (nrow(markers_all) > 0 && "cluster" %in% colnames(markers_all)) {
  top20 <- markers_all %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 20)
  
  top20_csv <- file.path(output_dir, "top20_markers_per_cluster.csv")
  write.csv(top20, top20_csv, row.names = FALSE)
  cat("✅ Top 20 markers per cluster saved to", top20_csv, "\n")
} else {
  cat("⚠️ No DE genes identified. Skipping top 20 marker export.\n")
}

# ------------------------- #
# Extract UMAP + Expression
# ------------------------- #
gene <- "NCR1"
if (!gene %in% rownames(seurat_obj)) {
  stop("❌ Gene ", gene, " not found in the Seurat object.")
}

# Extract UMAP coordinates
umap_coords <- Embeddings(seurat_obj, reduction = "umap") %>%
  as.data.frame() %>%
  `colnames<-`(c("UMAP_1", "UMAP_2")) %>%
  mutate(cell = rownames(.))

# Extract gene expression
expr_vec <- FetchData(seurat_obj, vars = gene)[, 1]
umap_coords$expression <- expr_vec

# ------------------------- #
# Plot using ggplot2
# ------------------------- #
cat("🎨 Creating custom UMAP FeaturePlot for", gene, "\n")

plot <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = expression)) +
  geom_point(size = 0.5) +
  scale_color_gradientn(
    colors = c("white", "lightgray", "red"),
    values = scales::rescale(c(0, 0.01, max(umap_coords$expression))),
    name = gene
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  ggtitle(paste0(gene, " Expression on UMAP"))

# ------------------------- #
# Save PDF
# ------------------------- #
pdf_path <- file.path(output_dir, paste0(gene, "_UMAP_d25_res0.3_custom.pdf"))
ggsave(pdf_path, plot, width = 6, height = 5, dpi = 600, bg = "transparent")
cat("✅ Custom FeaturePlot saved to", pdf_path, "\n")