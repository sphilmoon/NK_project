# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(ggplot2)
library(cowplot)

# ------------------------- #
# Define Paths and Configs
# ------------------------- #
totalNK_rds <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/rds/integrated_data_dims25_res0.3_genecounts.rds"
nkp46_rds <- "/home/outputs/nkp46_outputs/nkp46_integrated_data.rds"
output_dir <- "/home/outputs/all_merged_TotalNK_nkp46"

rds_dir <- file.path(output_dir, "rds")
figures_dir <- file.path(output_dir, "figures")
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------- #
# Load Seurat Objects
# ------------------------- #
totalNK <- readRDS(totalNK_rds)
nkp46_data <- readRDS(nkp46_rds)


# ------------------------- #
# Subset nkp46+ and nkp46- populations
# ------------------------- #
nkp46_pos <- subset(nkp46_data, subset = condition == "nkp46+")
nkp46_neg <- subset(nkp46_data, subset = condition == "nkp46-")
cat("✅ Subsetted NKp46+ and NKp46- populations\n")


# ------------------------- #
# Annotate Sample Identity
# ------------------------- #
totalNK$sample_id <- "TotalNK"
nkp46_pos$sample_id <- "NKp46pos"
nkp46_neg$sample_id <- "NKp46neg"
cat("✅ Annotated sample identities\n")

# Copy 'sample' column into a new column for tracking
totalNK$sample_id <- totalNK$sample_id
nkp46_pos$sample_id <- nkp46_pos$sample_id
nkp46_neg$sample_id <- nkp46_neg$sample_id
cat("✅ Copied 'sample' column to 'sample_id'\n")


cat("✅ Unique sample values in totalNK:\n")
print(unique(totalNK$sample_id))

cat("✅ Unique sample values in nkp46_pos:\n")
print(unique(nkp46_pos$sample_id))

cat("✅ Unique sample values in nkp46_neg:\n")
print(unique(nkp46_neg$sample_id))


# ------------------------- #
# Merge Seurat Objects
# ------------------------- #
merged_obj <- merge(totalNK, y = c(nkp46_pos, nkp46_neg), 
                    add.cell.ids = c("TotalNK", "NKp46pos", "NKp46neg"), 
                    project = "NK_Combined")
cat("✅ Merged Seurat objects\n")
DefaultAssay(merged_obj) <- "RNA"



# ------------------------- #
# Run PCA, Clustering, and UMAP
# ------------------------- #

merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj, npcs = 25)
merged_obj <- FindNeighbors(merged_obj, dims = 1:25)
merged_obj <- FindClusters(merged_obj, resolution = 0.5)
merged_obj <- RunUMAP(merged_obj, dims = 1:25)
cat("✅ Completed preprocessing and dimensional reduction\n")

# Save merged object
saveRDS(merged_obj, file.path(rds_dir, "all_totalNK_nkp46_dim25_res0.5_clustered.rds"))
cat("✅ Saved merged object\n")


# # read the merged object
# merged_obj <- readRDS(file.path(rds_dir, "all_totalNK_nkp46_dim25_res0.5_clustered.rds"))
# cat("✅ Loaded merged Seurat object\n")

# ------------------------- #
# Define Plotting Variables
# ------------------------- #
animals <- c("Animal25", "Animal26", "Animal27", "Animal28")
samples <- c("NKp46neg", "NKp46pos", "TotalNK")

# Validate and debug animal and sample presence
cat("🔍 Unique sample_id values:", unique(merged_obj$sample_id), "\n")
cat("📊 Cross-tabulation of sample_id vs sample_id:\n")
print(table(merged_obj$sample_id, merged_obj$sample_id))


# ------------------------- #
# Objective 1: UMAP by Animal and Sample
# ------------------------- #
plot_list <- list()
for (animal in animals) {
  for (sample in samples) {
    cells <- WhichCells(merged_obj, expression = sample_id == animal & sample_id == sample)
    if (length(cells) == 0) {
      cat("⚠️ No cells found for", animal, sample, "- skipping plot\n")
      next
    }
    p <- DimPlot(merged_obj, cells = cells, group.by = "sample_id") +
      ggtitle(paste(animal, sample)) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()
      )
    plot_list[[paste(animal, sample, sep = "_")]] <- p
  }
}
umap_grid <- plot_grid(plotlist = plot_list, ncol = length(samples))
ggsave(file.path(figures_dir, "UMAP_by_animal_and_sample.pdf"),
       plot = umap_grid, width = 12, height = 16, dpi = 600)
cat("✅ Saved UMAP plot grid\n")

# ------------------------- #
# Objective 2: NCR1 Expression FeaturePlot
# ------------------------- #
expr_vals <- FetchData(merged_obj, vars = "NCR1")[, 1]
expr_range <- range(expr_vals, na.rm = TRUE)

feature_plot_list <- list()
for (animal in animals) {
  for (sample in samples) {
    cells <- WhichCells(merged_obj, expression = sample_id == animal & sample_id == sample)
    if (length(cells) == 0) {
      cat("⚠️ No cells found for", animal, sample, "- skipping FeaturePlot\n")
      next
    }
    p <- FeaturePlot(merged_obj, features = "NCR1", cells = cells, order = TRUE) +
      ggtitle(paste(animal, sample)) +
      scale_color_gradientn(
        colors = c("lightgrey", "blue"),
        limits = expr_range,
        name = "NCR1"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()
      )
    feature_plot_list[[paste(animal, sample, sep = "_")]] <- p
  }
}
ncr1_grid <- plot_grid(plotlist = feature_plot_list, ncol = length(samples))
ggsave(file.path(figures_dir, "NCR1_expression_by_animal_and_sample.pdf"),
       plot = ncr1_grid, width = 12, height = 16, dpi = 600)
cat("✅ Saved NCR1 FeaturePlot grid\n")
