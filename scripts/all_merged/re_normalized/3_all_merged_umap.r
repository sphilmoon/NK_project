# Script: umap_ncr1_by_animal_sample.R
# Purpose: Generate UMAP and NCR1 FeaturePlots by animal and sample type

# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)

# ------------------------- #
# Define Paths and Configs
# ------------------------- #
totalNK_rds <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/rds/integrated_data_dims25_res0.3_genecounts.rds"
nkp46_rds <- "/home/outputs/nkp46_outputs/nkp46_integrated_data.rds"
output_dir <- "/home/outputs/all_merged_TotalNK_nkp46"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

rds_dir <- file.path(output_dir, "rds")
figures_dir <- file.path(output_dir, "figures")
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------- #
# Load Seurat Objects
# ------------------------- #
totalNK <- readRDS(totalNK_rds)
cat("✅ Loaded Total NK Seurat object from", totalNK_rds, "\n")

nkp46_data <- readRDS(nkp46_rds)
cat("✅ Loaded NKp46 Seurat object from", nkp46_rds, "\n")

# ------------------------- #
# Set Consistent Assay
# ------------------------- #
DefaultAssay(totalNK) <- "RNA"
DefaultAssay(nkp46_data) <- "RNA"
cat("✅ Set default assay to RNA for both objects\n")

totalNK <- JoinLayers(totalNK, assay = "RNA")
nkp46_data <- JoinLayers(nkp46_data, assay = "RNA")
cat("✅ Joined layers for both objects\n")

# ------------------------- #
# Validate Counts Slot
# ------------------------- #
validate_counts <- function(obj, name) {
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  if (is.null(counts)) stop(sprintf("❌ No counts slot found in %s", name))
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
  values <- summary(counts)$x
  if (!all(is.numeric(values))) {
    dense_counts <- as.matrix(counts[1:min(5, nrow(counts)), 1:min(5, ncol(counts))])
    cat(sprintf("⚠️ First 5x5 subset of %s counts slot:\n", name))
    print(dense_counts)
    stop(sprintf("❌ Non-numeric values found in counts slot of %s", name))
  }
  cat(sprintf("✅ %s counts slot contains only numeric values\n", name))
}

validate_counts(totalNK, "totalNK")
validate_counts(nkp46_data, "nkp46_data")

# ------------------------- #
# Subset NKp46+ and NKp46- Populations
# ------------------------- #
nkp46_pos <- subset(nkp46_data, subset = condition == "nkp46+")
nkp46_neg <- subset(nkp46_data, subset = condition == "nkp46-")
cat("✅ Subsetted NKp46+ and NKp46- populations\n")

validate_counts(nkp46_pos, "nkp46_pos")
validate_counts(nkp46_neg, "nkp46_neg")

# ------------------------- #
# Annotate Sample Identity
# ------------------------- #
totalNK$sample_id <- "TotalNK"
nkp46_pos$sample_id <- "NKp46pos"
nkp46_neg$sample_id <- "NKp46neg"
cat("✅ Annotated sample_id\n")

# Assign animal_id from sample (assuming sample contains animal info)
totalNK$animal_id <- totalNK$sample
nkp46_pos$animal_id <- nkp46_pos$sample
nkp46_neg$animal_id <- nkp46_neg$sample
cat("✅ Assigned animal_id from sample\n")

# Clean animal_id to remove NA and unexpected values (e.g., Animal52)
merged_obj$animal_id <- ifelse(merged_obj$animal_id %in% c("Animal25", "Animal26", "Animal27", "Animal28"), 
                              merged_obj$animal_id, NA)

# ------------------------- #
# Merge Seurat Objects
# ------------------------- #
merged_obj <- merge(totalNK, y = c(nkp46_pos, nkp46_neg), 
                    add.cell.ids = c("TotalNK", "NKp46pos", "NKp46neg"), 
                    project = "NK_Combined")
cat("✅ Merged Seurat objects\n")

DefaultAssay(merged_obj) <- "RNA"
merged_obj <- JoinLayers(merged_obj, assay = "RNA")
cat("✅ Set default assay to RNA and joined layers for merged object\n")

validate_counts(merged_obj, "merged_obj")

# ------------------------- #
# Normalize and Scale Data
# ------------------------- #
merged_obj <- NormalizeData(merged_obj)
cat("✅ Normalized data\n")

cat("🔍 Checking gene variances...\n")
data <- GetAssayData(merged_obj, assay = "RNA", layer = "data")
gene_vars <- apply(data, 1, var, na.rm = TRUE)
if (all(gene_vars == 0)) stop("❌ All genes have zero variance after normalization")
cat("✅ Found", sum(gene_vars > 0), "genes with non-zero variance\n")

non_zero_var_genes <- names(gene_vars[gene_vars > 0])
VariableFeatures(merged_obj) <- non_zero_var_genes
cat("✅ Set", length(non_zero_var_genes), "variable features with non-zero variance\n")

merged_obj <- FindVariableFeatures(merged_obj, selection.method = "vst", nfeatures = 2000)
cat("✅ Found variable features\n")

merged_obj <- ScaleData(merged_obj)
cat("✅ Scaled data\n")

# ------------------------- #
# Run PCA, Clustering, and UMAP
# ------------------------- #
merged_obj <- RunPCA(merged_obj, npcs = 25)
cat("✅ Completed PCA with 25 dimensions\n")

merged_obj <- FindNeighbors(merged_obj, dims = 1:25)
merged_obj <- FindClusters(merged_obj, resolution = 0.3)
merged_obj <- RunUMAP(merged_obj, dims = 1:25)
cat("✅ Completed clustering and UMAP with resolution 0.3\n")

# ------------------------- #
# Objective 1: UMAP by Animal and Sample Type
# ------------------------- #
animals <- c("Animal25", "Animal26", "Animal27", "Animal28")
samples <- c("NKp46neg", "NKp46pos", "TotalNK")

# Validate animals and samples
cat("🔍 Unique animal_id values:", unique(merged_obj$animal_id), "\n")
cat("🔍 Unique sample_id values:", unique(merged_obj$sample_id), "\n")
missing_animals <- animals[!animals %in% unique(merged_obj$animal_id)]
if (length(missing_animals) > 0) {
  warning("⚠️ Missing animals: ", paste(missing_animals, collapse = ", "))
  animals <- intersect(animals, unique(merged_obj$animal_id))
}
missing_samples <- samples[!samples %in% unique(merged_obj$sample_id)]
if (length(missing_samples) > 0) {
  stop("❌ Missing samples: ", paste(missing_samples, collapse = ", "))
}

umap_theme <- theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

umap_plot_list <- list()
for (animal in animals) {
  row_plots <- list()
  for (sample in samples) {
    subset_cells <- WhichCells(merged_obj, expression = animal_id == animal & sample_id == sample)
    cat("🔍 Cells for", animal, "and", sample, ":", length(subset_cells), "\n")
    if (length(subset_cells) == 0) {
      warning("⚠️ No cells found for ", animal, " and ", sample, ". Creating empty plot.")
      p <- ggplot() + geom_blank() + ggtitle(paste(animal, sample)) + umap_theme + coord_fixed(ratio = 1)
    } else {
      sample_color <- switch(sample, "NKp46neg" = "red", "NKp46pos" = "blue", "TotalNK" = "darkgreen")
      p <- DimPlot(merged_obj, cells = subset_cells, group.by = "sample_id", cols = sample_color,
                   reduction = "umap", pt.size = 0.5) + ggtitle(paste(animal, sample)) + umap_theme
    }
    row_plots[[sample]] <- p
  }
  row_patch <- wrap_plots(row_plots, ncol = length(samples))
  umap_plot_list[[animal]] <- row_patch
}

umap_full_grid <- wrap_plots(umap_plot_list, nrow = length(animals))
ggsave(file.path(figures_dir, "UMAP_by_animal_and_sample.pdf"), plot = umap_full_grid,
       width = 12, height = 16, dpi = 600, bg = "transparent")
cat("✅ Saved UMAP plot to", file.path(figures_dir, "UMAP_by_animal_and_sample.pdf"), "\n")

# ------------------------- #
# Objective 2: NCR1 Expression FeaturePlot
# ------------------------- #
if (!"NCR1" %in% rownames(merged_obj)) stop("❌ NCR1 gene not found")
expr_vals <- FetchData(merged_obj, vars = "NCR1")[, 1]
expr_range <- range(expr_vals, na.rm = TRUE)
cat("🔍 NCR1 expression range:", expr_range[1], "to", expr_range[2], "\n")

feature_plot_list <- list()
for (animal in animals) {
  row_plots <- list()
  for (sample in samples) {
    subset_cells <- WhichCells(merged_obj, expression = animal_id == animal & sample_id == sample)
    cat("🔍 Cells for", animal, "and", sample, ":", length(subset_cells), "\n")
    if (length(subset_cells) == 0) {
      warning("⚠️ No cells found for ", animal, " and ", sample, ". Creating empty plot.")
      p <- ggplot() + geom_blank() + ggtitle(paste(animal, sample)) + umap_theme + coord_fixed(ratio = 1)
    } else {
      p <- FeaturePlot(merged_obj, features = "NCR1", cells = subset_cells, reduction = "umap",
                       pt.size = 0.5, order = TRUE) +
           scale_color_gradientn(colors = c("lightgrey", "blue"), limits = expr_range,
                                 name = "NCR1 Expression") +
           ggtitle(paste(animal, sample)) + umap_theme
    }
    row_plots[[sample]] <- p
  }
  row_patch <- wrap_plots(row_plots, ncol = length(samples))
  feature_plot_list[[animal]] <- row_patch
}

feature_full_grid <- wrap_plots(feature_plot_list, nrow = length(animals))

legend_plot <- FeaturePlot(merged_obj, features = "NCR1", reduction = "umap", pt.size = 0.5, order = TRUE) +
               scale_color_gradientn(colors = c("lightgrey", "blue"), limits = expr_range,
                                     name = "NCR1 Expression") +
               theme(legend.position = "right")
legend <- cowplot::get_legend(legend_plot)

final_feature_plot <- plot_grid(feature_full_grid, legend, ncol = 2, rel_widths = c(1, 0.1))
ggsave(file.path(figures_dir, "NCR1_expression_by_animal_and_sample.pdf"), plot = final_feature_plot,
       width = 13, height = 16, dpi = 600, bg = "transparent")
cat("✅ Saved NCR1 FeaturePlot to", file.path(figures_dir, "NCR1_expression_by_animal_and_sample.pdf"), "\n")

# ------------------------- #
# Final Message
# ------------------------- #
cat("🎉 Visualization complete. Outputs saved in", figures_dir, "\n")