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
output_dir <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/pdf"

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
# Define Genes to Plot
# ------------------------- #
genes <- c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD40", "CD68") # CD14 is missing.

# Check if genes exist
expr_data_all <- GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = "data")
missing_genes <- genes[!genes %in% rownames(expr_data_all)]
if (length(missing_genes) > 0) {
  cat("⚠️ Missing genes:", paste(missing_genes, collapse = ", "), "\n")
  stop("❌ Genes not found in the Seurat object: ", paste(missing_genes, collapse = ", "))
}
cat("✅ All genes found in the Seurat object: ", paste(genes, collapse = ", "), "\n")

# ------------------------- #
# Check for UMAP and Sample Metadata
# ------------------------- #
if (!"umap" %in% names(seurat_obj@reductions)) {
  stop("❌ UMAP embedding not found in the Seurat object.")
}

if (!"sample" %in% colnames(seurat_obj@meta.data)) {
  stop("❌ Metadata column 'sample' not found. Please specify the correct column for samples.")
}

samples <- unique(seurat_obj$sample)
n_samples <- length(samples)
cat("✅ Found", n_samples, "samples: ", paste(samples, collapse = ", "), "\n")

# ------------------------- #
# Create Individual FeaturePlots
# ------------------------- #
plot_list <- list()

for (sample in samples) {
  for (gene in genes) {
    subset_obj <- subset(seurat_obj, subset = sample == !!sample)
    expr_values <- FetchData(subset_obj, vars = gene)
    if (all(expr_values[[gene]] == 0)) {
      p <- ggplot() + 
        theme_void() + 
        ggtitle(paste(gene, "in", sample, "(No Expression)"))
    } else {
      p <- FeaturePlot(
        subset_obj,
        features = gene,
        pt.size = 0.5,
        order = TRUE
      ) +
      scale_color_gradientn(
        colors = c("lightgrey", "blue", "red"),
        name = "Expression Level"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none"
      ) +
      ggtitle(paste(gene, "in", sample))
    }
    plot_list[[paste(sample, gene, sep = "_")]] <- p
  }
}

# ------------------------- #
# Arrange and Save Combined Plot
# ------------------------- #
# Arrange plots: rows = samples, columns = genes
ordered_plots <- list()
for (sample in samples) {
  for (gene in genes) {
    plot_name <- paste(sample, gene, sep = "_")
    ordered_plots <- c(ordered_plots, list(plot_list[[plot_name]]))
  }
}

combined_plot <- plot_grid(plotlist = ordered_plots, nrow = n_samples, ncol = length(genes))

# Save the combined plot
featureplot_file <- file.path(output_dir, "NK_QCmarkers_featureplot_by_sample_FIXED.pdf")
ggsave(
  filename = featureplot_file,
  plot = combined_plot,
  width = 4 * length(genes),  # Width per gene
  height = 4 * n_samples,     # Height per sample
  dpi = 600,
  bg = "transparent"
)
cat("✅ Combined FeaturePlot saved to", featureplot_file, "\n")
