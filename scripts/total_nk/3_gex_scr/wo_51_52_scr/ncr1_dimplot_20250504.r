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
dp_data <- dp_data %>%
  mutate(color = ifelse(avg.exp.scaled == 0, "white", avg.exp.scaled))

# ------------------------- #
# Custom ggplot with manual color mapping
# ------------------------- #
cat("🎨 Creating custom DotPlot with white for zero expression...\n")
custom_dotplot <- ggplot(dp_data, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = color)) +
  scale_size(range = c(2, 6), name = "Percent Expressed") +
  scale_color_gradientn(colors = c("lightgrey", "red"), 
                        limits = c(min(dp_data$color[dp_data$color != "white"]), max(dp_data$color[dp_data$color != "white"])),
                        na.value = "white", name = "Avg. Expression") +
  coord_flip() +
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
output_file <- file.path(output_dir, paste0("dotplot_custom_", gene, "_clusters.pdf"))
ggsave(filename = output_file, plot = custom_dotplot, width = 8, height = 4, dpi = 600)
cat("✅ Custom DotPlot saved to", output_file, "\n")