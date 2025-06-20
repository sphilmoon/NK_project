# ------------------------- #
# Load Libraries
# ------------------------- #
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)


# ------------------------- #
# Define Directory Structure
# ------------------------- #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
rds_output_dir <- file.path(output_dir, "rds")
pdf_output_dir <- file.path(output_dir, "pdf")
tsv_output_dir <- file.path(output_dir, "tsv")

# ------------------------- #
# Load Saved Seurat Analysis Objects
# ------------------------- #
individual_rds <- file.path(rds_output_dir, "individual_seurat_analysis_20250616.rds")
merged_rds <- file.path(rds_output_dir, "merged_seurat_analysis_20250616.rds")

# new individual and merged Seurat objects
# These objects should contain all the analyses performed in the previous steps, including clustering, UMAP, 
# and dimensionality reduction.
# Ensure these objects are saved with the latest analyses before running this script.
indie_obj <- readRDS(individual_rds)
cat("✅ Loaded individual Seurat objects.\n")
merged_obj <- readRDS(merged_rds)
cat("✅ Loaded merged Seurat object.\n")


# ------------------------- #
# Plotting Feature Plots for Individual and Merged Samples (SCT and LogNormalize)
# ------------------------- #
dims_list <- c(10, 15, 20, 25, 30)
resolutions <- c(0.25, 0.5, 0.75)
features <- c("NCR1", "KLRB1", "PRF1")

# ------------------------- #
# Individual Feature Plots
# ------------------------- #
for (method in c("SCT", "LogNormalize")) {
  for (animal in names(indie_obj)) {
    obj <- indie_obj[[animal]]

    for (dim in dims_list) {
      for (res in resolutions) {
        reduction_name <- paste0("umap_", method, "_dims", dim)

        if (!reduction_name %in% names(obj@reductions)) next

        plots <- FeaturePlot(obj,
                             features = features,
                             reduction = reduction_name,
                             pt.size = 0.3,
                             ncol = 3) +
                 plot_annotation(title = paste(animal, "-", method, "Dims:", dim, "Res:", res))

        file_name <- paste0("featureplots_", animal, "_", method, "_dims", dim, "_res", res, ".pdf")
        ggsave(file.path(pdf_output_dir, file_name), plots, width = 12, height = 6, dpi = 600)
        cat("✅ Saved:", file_name, "\n")
      }
    }
  }
}

# ------------------------- #
# Merged Feature Plots
# ------------------------- #
for (method in c("SCT", "LogNormalize")) {
  obj <- if (method == "SCT") merged_obj$SCT else merged_obj$LogNormalize
  assay_name <- ifelse(method == "SCT", "integrated", "RNA")
  DefaultAssay(obj) <- assay_name

  for (dim in dims_list) {
    for (res in resolutions) {
      reduction_name <- paste0("umap_", method, "_dims", dim)

      if (!reduction_name %in% names(obj@reductions)) next

      plots <- FeaturePlot(obj,
                           features = features,
                           reduction = reduction_name,
                           pt.size = 0.3,
                           ncol = 3) +
               plot_annotation(title = paste("Merged -", method, "Dims:", dim, "Res:", res))

      file_name <- paste0("featureplots_merged_", method, "_dims", dim, "_res", res, ".pdf")
      ggsave(file.path(pdf_output_dir, file_name), plots, width = 12, height = 6, dpi = 600)
      cat("✅ Saved:", file_name, "\n")
    }
  }
}




# ------------------------- #
# Run JackStraw on Individual Samples (RNA assay)
# ------------------------- #
jackstraw_plots <- list()

for (sample in names(indie_obj)) {
  cat("📊 JackStraw for:", sample, "\n")
  obj <- indie_obj[[sample]]
  DefaultAssay(obj) <- "RNA"

  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj <- JackStraw(obj, num.replicate = 100)
  obj <- ScoreJackStraw(obj, dims = 1:30)

  p <- JackStrawPlot(obj, dims = 1:30) + ggtitle(paste("JackStraw -", sample))
  jackstraw_plots[[sample]] <- p

  indie_obj[[sample]] <- obj  # update object
}

# Save to PDF
pdf(file.path(pdf_output_dir, "jackstraw_individual_samples_20250616.pdf"), width = 10, height = 6 * ceiling(length(jackstraw_plots) / 2))
print(wrap_plots(jackstraw_plots, ncol = 2))
dev.off()
cat("✅ JackStraw plots saved for individual samples.\n")

# ------------------------- #
# Run JackStraw on Merged Sample (RNA assay)
# ------------------------- #
cat("📊 JackStraw for: merged object\n")
DefaultAssay(merged_obj) <- "RNA"

merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj)
merged_obj <- JackStraw(merged_obj, num.replicate = 100)
merged_obj <- ScoreJackStraw(merged_obj, dims = 1:30)

jack_plot <- JackStrawPlot(merged_obj, dims = 1:30) + ggtitle("JackStraw - Merged")

ggsave(
  filename = file.path(pdf_output_dir, "jackstraw_merged_totalNK_20250616.pdf"),
  plot = jack_plot,
  width = 10,
  height = 6,
  dpi = 600
)
cat("✅ JackStraw plot saved for merged object.\n")



# # ------------------------- #
# # Optional: Save updated objects
# # ------------------------- #
# saveRDS(indie_obj, file.path(rds_output_dir, "individual_seurat_analysis_20250616.rds"))
# saveRDS(merged_obj, file.path(rds_output_dir, "merged_seurat_analysis_20250616.rds"))
# cat("✅ Updated RDS files with JackStraw results saved.\n")