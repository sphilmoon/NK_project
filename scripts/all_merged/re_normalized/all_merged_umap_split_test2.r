# Script: umap_per_animal_condition.R
# Purpose: Generate separate UMAP plots per animal for each condition and save RDS files

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
figures_dir <- file.path(output_dir, "re_normalized/figures")
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

# Assign animal_id from sample column
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
# Process Each Animal and Condition Separately
# ------------------------- #
umap_plots <- list()
objects_to_save <- list()

for (animal in animals) {
    for (condition in c("TotalNK", "NKp46pos", "NKp46neg")) {
        # Determine the source object based on condition
        if (condition == "TotalNK") {
            obj <- subset(totalNK, subset = animal_id == animal)
        } else if (condition == "NKp46pos") {
            obj <- subset(nkp46_pos, subset = animal_id == animal)
        } else if (condition == "NKp46neg") {
            obj <- subset(nkp46_neg, subset = animal_id == animal)
        }

        if (ncol(obj) == 0) {
            cat("⚠️ No cells found for", animal, "and", condition, "- skipping\n")
            next
        }

        # Normalize, Scale, and Run Dimensionality Reduction
        obj <- NormalizeData(obj)
        obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
        obj <- ScaleData(obj)
        obj <- RunPCA(obj, npcs = 25)
        obj <- FindNeighbors(obj, dims = 1:25)
        obj <- FindClusters(obj, resolution = 0.3)
        obj <- RunUMAP(obj, dims = 1:25)
        cat("✅ Processed UMAP for", animal, "and", condition, "\n")

        # Create UMAP plot
        umap_plot <- DimPlot(obj, group.by = "sample_id", pt.size = 0.5) +
            theme_minimal() +
            ggtitle(paste(animal, condition, sep = " - ")) +
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
        umap_plots[[paste(animal, condition, sep = "_")]] <- umap_plot

        # Save RDS file
        rds_filename <- file.path(rds_dir, paste0(animal, "_", condition, ".rds"))
        saveRDS(obj, file = rds_filename)
        cat("✅ Saved RDS for", animal, "and", condition, "to", rds_filename, "\n")
        objects_to_save[[paste(animal, condition)]] <- rds_filename
    }
}

# ------------------------- #
# Arrange and Save Plots
# ------------------------- #
umap_grid <- wrap_plots(umap_plots, nrow = length(animals), ncol = length(c("TotalNK", "NKp46pos", "NKp46neg")))
ggsave(file.path(figures_dir, "UMAP_per_animal_condition.pdf"),
    plot = umap_grid,
    width = 18, height = 16, dpi = 600, bg = "transparent"
)
cat("✅ Saved UMAP per animal and condition to", file.path(figures_dir, "UMAP_per_animal_condition.pdf"), "\n")

# ------------------------- #
# Final Message
# ------------------------- #
cat("🎉 UMAP generation complete. RDS files saved in:", rds_dir, "\n")
cat("Outputs and RDS files available in:", output_dir, "\n")
