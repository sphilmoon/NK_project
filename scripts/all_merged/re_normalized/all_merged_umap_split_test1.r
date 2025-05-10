# Script: umap_per_animal.R
# Purpose: Generate UMAP plots per animal for TotalNK, NKp46+, and NKp46- conditions

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
# Annotate Sample Identity and Animal ID
# ------------------------- #
totalNK$sample_id <- "TotalNK"
nkp46_pos$sample_id <- "NKp46pos"
nkp46_neg$sample_id <- "NKp46neg"
cat("✅ Annotated sample_id\n")

# Assign animal_id from sample column (assuming sample contains Animal25–28)
if ("sample" %in% colnames(totalNK@meta.data)) {
    totalNK$animal_id <- totalNK$sample
    cat("✅ Assigned animal_id from 'sample' for totalNK\n")
} else {
    stop("❌ 'sample' column not found in totalNK metadata")
}

if ("animal" %in% colnames(nkp46_data@meta.data)) {
    nkp46_pos$animal_id <- nkp46_pos$animal
    nkp46_neg$animal_id <- nkp46_neg$animal
    cat("✅ Assigned animal_id from 'animal' for nkp46_pos and nkp46_neg\n")
} else if ("sample" %in% colnames(nkp46_data@meta.data)) {
    nkp46_pos$animal_id <- nkp46_pos$sample
    nkp46_neg$animal_id <- nkp46_neg$sample
    cat("✅ Assigned animal_id from 'sample' for nkp46_pos and nkp46_neg\n")
} else {
    stop("❌ No 'animal' or 'sample' column found in nkp46_data metadata")
}

# Clean animal_id to include only Animal25–28
animals <- c("Animal25", "Animal26", "Animal27", "Animal28")
totalNK$animal_id <- ifelse(totalNK$animal_id %in% animals, totalNK$animal_id, NA)
nkp46_pos$animal_id <- ifelse(nkp46_pos$animal_id %in% animals, nkp46_pos$animal_id, NA)
nkp46_neg$animal_id <- ifelse(nkp46_neg$animal_id %in% animals, nkp46_neg$animal_id, NA)

# ------------------------- #
# Process Each Animal Separately
# ------------------------- #
umap_plots <- list()

for (animal in animals) {
    # Subset data for the current animal across all sample types
    totalNK_subset <- subset(totalNK, subset = animal_id == animal)
    nkp46_pos_subset <- subset(nkp46_pos, subset = animal_id == animal)
    nkp46_neg_subset <- subset(nkp46_neg, subset = animal_id == animal)

    # Combine subsets for this animal
    animal_obj <- merge(totalNK_subset,
        y = c(nkp46_pos_subset, nkp46_neg_subset),
        add.cell.ids = c("TotalNK", "NKp46pos", "NKp46neg")
    )

    # Normalize, Scale, and Run Dimensionality Reduction
    animal_obj <- NormalizeData(animal_obj)
    animal_obj <- FindVariableFeatures(animal_obj, selection.method = "vst", nfeatures = 2000)
    animal_obj <- ScaleData(animal_obj)
    animal_obj <- RunPCA(animal_obj, npcs = 25)
    animal_obj <- FindNeighbors(animal_obj, dims = 1:25)
    animal_obj <- FindClusters(animal_obj, resolution = 0.3)
    animal_obj <- RunUMAP(animal_obj, dims = 1:25)
    cat("✅ Processed UMAP for", animal, "\n")

    # Create UMAP plot for this animal
    umap_plot <- DimPlot(animal_obj, group.by = "sample_id", pt.size = 0.5) +
        theme_minimal() +
        ggtitle(paste("UMAP for", animal)) +
        theme(
            plot.title = element_text(hjust = 0.5, size = 14),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 12),
            legend.position = "right",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 8),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        )
    umap_plots[[animal]] <- umap_plot
}

# ------------------------- #
# Arrange and Save Plots
# ------------------------- #
umap_grid <- plot_grid(plotlist = umap_plots, ncol = 2)
ggsave(file.path(figures_dir, "UMAP_per_animal.pdf"),
    plot = umap_grid,
    width = 12, height = 12, dpi = 600, bg = "transparent"
)
cat("✅ Saved UMAP per animal to", file.path(figures_dir, "UMAP_per_animal.pdf"), "\n")

# ------------------------- #
# Final Message
# ------------------------- #
cat("🎉 UMAP generation complete. Outputs saved in:", figures_dir, "\n")
