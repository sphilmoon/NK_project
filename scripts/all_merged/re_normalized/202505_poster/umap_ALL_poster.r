# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(ggplot2)
library(cowplot)
library(stringr)

# ------------------------- #
# Define Input and Output Paths
# ------------------------- #
nkp_neg_rds <- "/home/outputs/nkp46_outputs/chat/neg/rds/nkp46_neg.rds"
nkp_pos_rds <- "/home/outputs/nkp46_outputs/chat/pos/rds/nkp46_pos.rds"
totalNK_rds <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/rds/integrated_data_dims25_res0.3_genecounts.rds"

output_dir <- "/home/outputs/all_merged_TotalNK_nkp46/chat"
all_merged_rds_dir <- file.path(output_dir, "rds")
pdf_dir <- file.path(output_dir, "pdf")

dir.create(all_merged_rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------- #
# Load Seurat Objects
# ------------------------- #
cat("📦 Loading Seurat objects...\n")
nkp_neg <- readRDS(nkp_neg_rds)
cat("✅ Loaded nkp46_neg\n")
nkp_pos <- readRDS(nkp_pos_rds)
cat("✅ Loaded nkp46_pos\n")
totalNK <- readRDS(totalNK_rds)
cat("✅ Loaded totalNK\n")

# ------------------------- #
# Standardize Metadata
# ------------------------- #
cat("🔧 Standardizing metadata (sample, sample_id, condition, animal_id)...\n")

for (obj_name in c("nkp_neg", "nkp_pos", "totalNK")) {
    obj <- get(obj_name)

    # Convert to title case
    obj$sample <- str_to_title(obj$sample)

    # Copy fields
    obj$sample_id <- obj$sample
    obj$animal_id <- obj$sample

    # Set condition
    obj$condition <- switch(obj_name,
        "nkp_neg" = "nkp46_neg",
        "nkp_pos" = "nkp46_pos",
        "totalNK" = "totalNK"
    )

    assign(obj_name, obj)

    # Verify
    cat("🔍 Checking 'sample' field in ", obj_name, "\n")
    print(colnames(obj@meta.data))
    cat("✅ 'sample' values:\n")
    print(head(obj$sample))
}

# # diagnostic checks before the merge
# cat("🔍 Previewing cell names from each object...\n")
# cat("nkp_neg cell names example:\n")
# print(head(colnames(nkp_neg)))
# cat("nkp_pos cell names example:\n")
# print(head(colnames(nkp_pos)))
# cat("totalNK cell names example:\n")
# print(head(colnames(totalNK)))

# # Check for overlaps
# overlap1 <- intersect(colnames(nkp_neg), colnames(nkp_pos))
# overlap2 <- intersect(colnames(nkp_neg), colnames(totalNK))
# overlap3 <- intersect(colnames(nkp_pos), colnames(totalNK))

# cat("🔎 Number of overlapping cell names:\n")
# cat("nkp_neg vs nkp_pos:", length(overlap1), "\n")
# cat("nkp_neg vs totalNK:", length(overlap2), "\n")
# cat("nkp_pos vs totalNK:", length(overlap3), "\n")

# ------------------------- #
# Merge Seurat Objects
# ------------------------- #

# ✅ Rename cells to avoid name collisions
nkp_neg <- RenameCells(nkp_neg, add.cell.id = "nkp46neg")
nkp_pos <- RenameCells(nkp_pos, add.cell.id = "nkp46pos")
totalNK <- RenameCells(totalNK, add.cell.id = "totalNK")

# removing SCTransform normalization separately.
nkp_neg[["SCT"]] <- NULL
nkp_pos[["SCT"]] <- NULL
totalNK[["SCT"]] <- NULL

# 🔗 Merge all three objects safely
cat("🔗 Merging Seurat objects...\n")
merged_obj <- merge(nkp_neg, y = list(nkp_pos, totalNK))


DefaultAssay(merged_obj) <- "RNA"
cat("✅ Merge complete: total cells =", ncol(merged_obj), "\n")

# ------------------------- #
# Normalize, PCA, and Clustering
# ------------------------- #
cat("⚙️ Normalizing and scaling merged object...\n")
merged_obj <- SCTransform(merged_obj, verbose = TRUE)

# merged_obj <- NormalizeData(merged_obj)
# merged_obj <- FindVariableFeatures(merged_obj)
# merged_obj <- ScaleData(merged_obj)

cat("🔬 Running PCA...\n")
merged_obj <- RunPCA(merged_obj, npcs = 50)
cat("✅ PCA complete\n")

cat("📊 Saving ElbowPlot for dimensionality selection...\n")
elbow_plot <- ElbowPlot(merged_obj, ndims = 50)
ggsave(file.path(pdf_dir, "elbow_plot.pdf"), elbow_plot, width = 8, height = 6, dpi = 600)

dims_to_use <- 1:25

cat("🧭 Running UMAP, neighbors, and clustering (dims 1:25, res 0.5)...\n")
merged_obj <- RunUMAP(merged_obj, dims = dims_to_use)
merged_obj <- FindNeighbors(merged_obj, dims = dims_to_use)
merged_obj <- FindClusters(merged_obj, resolution = 0.5)
cat("✅ UMAP and clustering complete\n")

saveRDS(merged_obj, file = file.path(all_merged_rds_dir, "merged_totalNK_nkp46_dim25_res0.5.rds"))
cat("💾 Merged object saved\n")

# ------------------------- #
# UMAP by Condition
# ------------------------- #
cat("🎨 Generating combined UMAP by condition...\n")
umap_by_condition <- DimPlot(merged_obj, group.by = "condition", pt.size = 0.3) +
    ggtitle("UMAP by Condition") +
    theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(pdf_dir, "UMAP_by_condition_combined.pdf"),
    umap_by_condition,
    width = 10, height = 8, dpi = 600
)

# ------------------------- #
# UMAP Grid: Animal x Condition
# ------------------------- #
cat("🎨 Generating UMAP grid by animal x condition...\n")
conditions <- c("nkp46_neg", "nkp46_pos", "totalNK")
animals <- c("Animal25", "Animal26", "Animal27", "Animal28")
plot_list <- list()

for (animal in animals) {
    for (cond in conditions) {
        cells <- WhichCells(merged_obj, expression = animal_id == animal & condition == cond)
        cat(sprintf("🔍 Found %d cells for %s x %s\n", length(cells), animal, cond))

        if (length(cells) == 0) {
            cat("⚠️ No cells found for", animal, cond, "- skipping\n")
            next
        }

        p <- DimPlot(merged_obj, cells = cells, group.by = "condition") +
            ggtitle(paste(animal, cond)) +
            theme_minimal() +
            theme(
                legend.position = "none",
                plot.title = element_text(size = 10, hjust = 0.5),
                axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank()
            )
        plot_list[[paste(animal, cond, sep = "_")]] <- p
    }
}

cat("🧱 Assembling grid layout...\n")
umap_grid <- plot_grid(plotlist = plot_list, ncol = length(conditions))
ggsave(file.path(pdf_dir, "UMAP_animal_condition_grid.pdf"),
    umap_grid,
    width = 4 * length(conditions), height = 4 * length(animals), dpi = 600
)

cat("✅ All UMAP visualizations complete\n")
