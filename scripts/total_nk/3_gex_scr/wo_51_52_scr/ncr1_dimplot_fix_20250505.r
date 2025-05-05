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
# Extract DotPlot Data
# ------------------------- #
cat("🔍 Extracting DotPlot data...\n")
dp_data <- DotPlot(seurat_obj, features = gene)$data

# ------------------------- #
# Custom ggplot with manual color mapping
# ------------------------- #
cat("🎨 Creating custom DotPlot with white for zero expression...\n")
custom_dotplot <- ggplot(dp_data, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp)) +  # Use avg.exp instead of avg.exp.scaled
  scale_size(range = c(2, 6), name = "Percent Expressed") +
  scale_color_gradientn(colors = c("white", "lightgrey", "red"), 
                        values = scales::rescale(c(0, 0.01, max(dp_data$avg.exp, na.rm = TRUE))),
                        name = "Avg. Expression (Normalized)") +
  theme_minimal() +
  ggtitle(paste("DotPlot of", gene, "Expression Across Clusters")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

# ------------------------- #
# Save Plot
# ------------------------- #
output_file <- file.path(output_dir, paste0("dotplot_custom_", gene, "_clusters_fix.pdf"))
ggsave(filename = output_file, plot = custom_dotplot, width = 8, height = 4, dpi = 600, bg = "transparent")
cat("✅ Custom DotPlot saved to", output_file, "\n")





# ------------------------- #
# UMAP: Expression of NCR1 with Clusters
# ------------------------- #
cat("🎨 Creating FeaturePlot for", gene, "on UMAP with cluster labels...\n")

# Ensure UMAP and clustering are present
if (!"umap" %in% names(seurat_obj@reductions)) {
  stop("❌ UMAP embedding not found in the Seurat object.")
}
if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  stop("❌ Clustering not found. Run FindClusters() first.")
}

# Create FeaturePlot
ncr1_umap <- FeaturePlot(seurat_obj, features = gene, reduction = "umap", pt.size = 0.5) +
  ggtitle(paste("NCR1 Expression on UMAP")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# Extract UMAP embeddings and add cluster info
umap_df <- as.data.frame(Embeddings(seurat_obj, reduction = "umap"))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")  # Explicitly rename columns to ensure consistency
umap_df$cluster <- seurat_obj$seurat_clusters

# Compute cluster centers
cluster_centers <- umap_df %>%
  group_by(cluster) %>%
  summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2), .groups = "drop")

# Add cluster labels to the FeaturePlot
ncr1_umap_labeled <- ncr1_umap +
  geom_text(data = cluster_centers, aes(x = UMAP_1, y = UMAP_2, label = cluster),
            color = "black", size = 4, fontface = "bold")

# ------------------------- #
# Save Plot
# ------------------------- #
featureplot_file <- file.path(output_dir, paste0("featureplot_", gene, "_with_clusters.pdf"))
ggsave(filename = featureplot_file, plot = ncr1_umap_labeled, width = 8, height = 6, dpi = 600, bg = "transparent")
cat("✅ FeaturePlot saved to", featureplot_file, "\n")