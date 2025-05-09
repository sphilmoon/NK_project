# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(ggplot2)
library(cowplot)

# ------------------------- #
# Define Paths and Configs
# ------------------------- #
# totalNK_rds <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/rds/integrated_data_dims25_res0.3_genecounts.rds"
# nkp46_rds <- "/home/outputs/nkp46_outputs/nkp46_integrated_data.rds"
output_dir <- "/home/outputs/all_merged_TotalNK_nkp46"

# Create output subdirectories
rds_dir <- file.path(output_dir, "rds")
figures_dir <- file.path(output_dir, "figures")
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# # ------------------------- #
# # Load Merged Seurat Object
# # ------------------------- #
# merged_obj <- readRDS(file.path(rds_dir, "all_totalNK_nkp46.rds"))
# cat("✅ Loaded merged Seurat object\n")

# # ------------------------- #
# # Run PCA, Clustering, and UMAP
# # ------------------------- #
# merged_obj <- RunPCA(merged_obj, npcs = 25)
# cat("✅ Ran PCA with 25 dimensions\n")

# merged_obj <- FindNeighbors(merged_obj, dims = 1:25)
# merged_obj <- FindClusters(merged_obj, resolution = 0.5)
# merged_obj <- RunUMAP(merged_obj, dims = 1:25)
# cat("✅ Completed clustering and UMAP with resolution 0.5\n")

# # save the merged object
# saveRDS(merged_obj, file.path(rds_dir, "all_totalNK_nkp46_dim25_res0.5_clustered.rds"))
# cat("✅ Saved merged Seurat object\n")

# read the merged object
merged_obj <- readRDS(file.path(rds_dir, "all_totalNK_nkp46_dim25_res0.5_clustered.rds"))
cat("✅ Loaded merged Seurat object\n")

# Validate each object
validate_counts(totalNK, "totalNK")
validate_counts(nkp46_data, "nkp46_data")

# ------------------------- #
# Subset nkp46+ and nkp46- populations
# ------------------------- #
nkp46_pos <- subset(nkp46_data, subset = condition == "nkp46+")
nkp46_neg <- subset(nkp46_data, subset = condition == "nkp46-")
cat("✅ Subsetted NKp46+ and NKp46- populations\n")

# Validate subsets
validate_counts(nkp46_pos, "nkp46_pos")
validate_counts(nkp46_neg, "nkp46_neg")

# ------------------------- #
# Annotate Sample Identity
# ------------------------- #
totalNK$sample_id <- "TotalNK"
nkp46_pos$sample_id <- "NKp46pos"
nkp46_neg$sample_id <- "NKp46neg"
cat("✅ Annotated sample identities\n")

# ------------------------- #
# Merge Seurat Objects
# ------------------------- #
merged_obj <- merge(totalNK, y = c(nkp46_pos, nkp46_neg), 
                    add.cell.ids = c("TotalNK", "NKp46pos", "NKp46neg"), 
                    project = "NK_Combined")
cat("✅ Merged Seurat objects\n")

# ------------------------- #
# Annotate Metadata
# ------------------------- #
merged_obj$animal_id <- merged_obj$sample
cat("🔍 Assigned animal_id values (first 10):", head(merged_obj$animal_id, 10), "\n")
table(merged_obj$animal_id, merged_obj$sample_id)


if (!"sample_id" %in% colnames(merged_obj@meta.data)) {
  stop("❌ 'sample_id' not found in metadata. Ensure it is set before merging.")
}

# ------------------------- #
# Define Plotting Order
# ------------------------- #
animals <- c("Animal25", "Animal26", "Animal27", "Animal28")
samples <- c("NKp46-", "NKp46+", "TotalNK")

# ------------------------- #
# Objective 1: DimPlot by Animal × Sample
# ------------------------- #
plot_list <- list()
for (animal in animals) {
  for (sample in samples) {
    subset_cells <- WhichCells(merged_obj, idents = NULL)[
  merged_obj$animal_id == animal & merged_obj$sample_id == sample
]
    if (length(subset_cells) == 0) {
      cat("⚠️ No cells found for", animal, sample, "- skipping plot \n")
      next
    }

    p <- DimPlot(merged_obj, cells = subset_cells, group.by = "sample_id") +
      ggtitle(paste(animal, sample)) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 10)
      )
    plot_list[[paste(animal, sample, sep = "_")]] <- p
  }
}

# Combine into grid
umap_grid <- plot_grid(plotlist = plot_list, ncol = length(samples))
ggsave(
  filename = file.path(figures_dir, "UMAP_by_animal_and_sample.pdf"),
  plot = umap_grid, width = 12, height = 16, dpi = 600
)
cat("✅ Saved UMAP plot by animal and sample\n")

# ------------------------- #
# Objective 2: NCR1 Expression on UMAP
# ------------------------- #
# Compute global expression range for NCR1
expr_vals <- FetchData(merged_obj, vars = "NCR1")[, 1]
expr_range <- range(expr_vals, na.rm = TRUE)

feature_plot_list <- list()
for (animal in animals) {
  for (sample in samples) {
    subset_cells <- WhichCells(merged_obj, idents = NULL)[
  merged_obj$animal_id == animal & merged_obj$sample_id == sample
]
    if (length(subset_cells) == 0) {
      cat("⚠️ No cells found for", animal, sample, "- skipping plot \n")
      next
    }

    p <- FeaturePlot(merged_obj, features = "NCR1", cells = subset_cells, order = TRUE) +
      ggtitle(paste(animal, sample)) +
      scale_color_gradientn(
        colors = c("lightgrey", "blue"),
        limits = expr_range,
        name = "NCR1"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 10)
      )
    feature_plot_list[[paste(animal, sample, sep = "_")]] <- p
  }
}

# Combine into grid
ncr1_grid <- plot_grid(plotlist = feature_plot_list, ncol = length(samples))
ggsave(
  filename = file.path(figures_dir, "NCR1_expression_by_animal_and_sample.pdf"),
  plot = ncr1_grid, width = 12, height = 16, dpi = 600
)
cat("✅ Saved NCR1 FeaturePlot by animal and sample\n")
