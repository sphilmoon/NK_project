# # ------------------------- #
# # Load Required Libraries
# # ------------------------- #
# library(Seurat)
# library(ggplot2)
# library(cowplot)
# library(stringr)

# # ------------------------- #
# # Define Input and Output Paths
# # ------------------------- #
# nkp_neg_rds <- "/home/outputs/nkp46_outputs/chat/neg/rds/nkp46_neg.rds"
# nkp_pos_rds <- "/home/outputs/nkp46_outputs/chat/pos/rds/nkp46_pos.rds"
# totalNK_rds <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/rds/integrated_data_dims25_res0.3_genecounts.rds"

# output_dir <- "/home/outputs/all_merged_TotalNK_nkp46/chat"
# all_merged_rds_dir <- file.path(output_dir, "rds")
# pdf_dir <- file.path(output_dir, "pdf")

# dir.create(all_merged_rds_dir, recursive = TRUE, showWarnings = FALSE)
# dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

# # ------------------------- #
# # Load Seurat Objects
# # ------------------------- #
# cat("📦 Loading Seurat objects...\n")
# nkp_neg <- readRDS(nkp_neg_rds)
# cat("✅ Loaded nkp46_neg\n")
# nkp_pos <- readRDS(nkp_pos_rds)
# cat("✅ Loaded nkp46_pos\n")
# totalNK <- readRDS(totalNK_rds)
# cat("✅ Loaded totalNK\n")

# # ------------------------- #
# # Standardize Metadata
# # ------------------------- #
# cat("🔧 Standardizing metadata (sample, sample_id, condition, animal_id)...\n")

# for (obj_name in c("nkp_neg", "nkp_pos", "totalNK")) {
#     obj <- get(obj_name)

#     # Convert to title case
#     obj$sample <- str_to_title(obj$sample)

#     # Copy fields
#     obj$sample_id <- obj$sample
#     obj$animal_id <- obj$sample

#     # Set condition
#     obj$condition <- switch(obj_name,
#         "nkp_neg" = "nkp46_neg",
#         "nkp_pos" = "nkp46_pos",
#         "totalNK" = "totalNK"
#     )

#     assign(obj_name, obj)

#     # Verify
#     cat("🔍 Checking 'sample' field in ", obj_name, "\n")
#     print(colnames(obj@meta.data))
#     cat("✅ 'sample' values:\n")
#     print(head(obj$sample))
# }

# # # diagnostic checks before the merge
# # cat("🔍 Previewing cell names from each object...\n")
# # cat("nkp_neg cell names example:\n")
# # print(head(colnames(nkp_neg)))
# # cat("nkp_pos cell names example:\n")
# # print(head(colnames(nkp_pos)))
# # cat("totalNK cell names example:\n")
# # print(head(colnames(totalNK)))

# # # Check for overlaps
# # overlap1 <- intersect(colnames(nkp_neg), colnames(nkp_pos))
# # overlap2 <- intersect(colnames(nkp_neg), colnames(totalNK))
# # overlap3 <- intersect(colnames(nkp_pos), colnames(totalNK))

# # cat("🔎 Number of overlapping cell names:\n")
# # cat("nkp_neg vs nkp_pos:", length(overlap1), "\n")
# # cat("nkp_neg vs totalNK:", length(overlap2), "\n")
# # cat("nkp_pos vs totalNK:", length(overlap3), "\n")

# # ------------------------- #
# # Merge Seurat Objects
# # ------------------------- #

# # ✅ Rename cells to avoid name collisions
# nkp_neg <- RenameCells(nkp_neg, add.cell.id = "nkp46neg")
# nkp_pos <- RenameCells(nkp_pos, add.cell.id = "nkp46pos")
# totalNK <- RenameCells(totalNK, add.cell.id = "totalNK")

# # removing SCTransform normalization separately, to enable merging with a new normalization.
# nkp_neg[["SCT"]] <- NULL
# nkp_pos[["SCT"]] <- NULL
# totalNK[["SCT"]] <- NULL

# # 🔗 Merge all three objects safely
# cat("🔗 Merging Seurat objects...\n")
# merged_obj <- merge(nkp_neg, y = list(nkp_pos, totalNK), merge.data = FALSE)

# DefaultAssay(merged_obj) <- "RNA"
# cat("✅ Merge complete: total cells =", ncol(merged_obj), "\n")

# # ------------------------- #
# # Normalize, PCA, and Clustering
# # ------------------------- #
# cat("⚙️ Running SCTransform normalization and scaling merged object...\n")
# merged_obj <- SCTransform(merged_obj, verbose = TRUE)

# cat("🔬 Running PCA...\n")
# merged_obj <- RunPCA(merged_obj, npcs = 50)
# cat("✅ PCA complete\n")

# cat("📊 Saving ElbowPlot for dimensionality selection...\n")
# elbow_plot <- ElbowPlot(merged_obj, ndims = 50)
# ggsave(file.path(pdf_dir, "elbow_plot.pdf"), elbow_plot, width = 8, height = 6, dpi = 600)

# dims_to_use <- 1:25

# cat("🧭 Running UMAP, neighbors, and clustering (dims 1:25, res 0.5)...\n")
# merged_obj <- RunUMAP(merged_obj, dims = dims_to_use)
# merged_obj <- FindNeighbors(merged_obj, dims = dims_to_use)
# merged_obj <- FindClusters(merged_obj, resolution = 0.5)
# cat("✅ UMAP and clustering complete\n")

# saveRDS(merged_obj, file = file.path(all_merged_rds_dir, "merged_totalNK_nkp46_dim25_res0.5.rds"))
# cat("💾 Merged object saved\n")

# # open the rds files
# merged_obj <- file.path(all_merged_rds_dir, "merged_totalNK_nkp46_dim25_res0.5.rds")

# # ------------------------- #
# # UMAP by Condition
# # ------------------------- #
# cat("🎨 Generating combined UMAP by condition...\n")
# umap_by_condition <- DimPlot(merged_obj, group.by = "condition", pt.size = 0.3) +
#     ggtitle("UMAP by Condition") +
#     theme(plot.title = element_text(hjust = 0.5))

# ggsave(file.path(pdf_dir, "UMAP_by_condition_combined.pdf"),
#     umap_by_condition,
#     width = 10, height = 8, dpi = 600
# )

# # ------------------------- #
# # UMAP Grid: Animal x Condition
# # ------------------------- #
# cat("🎨 Generating UMAP grid by animal x condition...\n")
# conditions <- c("nkp46_neg", "nkp46_pos", "totalNK")
# animals <- c("Animal25", "Animal26", "Animal27", "Animal28")
# plot_list <- list()

# for (animal in animals) {
#     for (cond in conditions) {
#         cells <- WhichCells(merged_obj, expression = animal_id == animal & condition == cond)
#         cat(sprintf("🔍 Found %d cells for %s x %s\n", length(cells), animal, cond))

#         if (length(cells) == 0) {
#             cat("⚠️ No cells found for", animal, cond, "- skipping\n")
#             next
#         }

#         p <- DimPlot(merged_obj, cells = cells, group.by = "condition") +
#             ggtitle(paste(animal, cond)) +
#             theme_minimal() +
#             theme(
#                 legend.position = "none",
#                 plot.title = element_text(size = 10, hjust = 0.5),
#                 axis.text = element_blank(),
#                 axis.title = element_blank(),
#                 axis.ticks = element_blank()
#             )
#         plot_list[[paste(animal, cond, sep = "_")]] <- p
#     }
# }

# cat("🧱 Assembling grid layout...\n")
# umap_grid <- plot_grid(plotlist = plot_list, ncol = length(conditions))
# ggsave(file.path(pdf_dir, "UMAP_animal_condition_grid.pdf"),
#     umap_grid,
#     width = 4 * length(conditions), height = 4 * length(animals), dpi = 600
# )

# cat("✅ All UMAP visualizations complete\n")




# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(ggplot2)
library(cowplot)
library(stringr)
library(patchwork)

# ------------------------- #
# Define Input and Output Paths
# ------------------------- #
output_dir <- "/home/outputs/all_merged_TotalNK_nkp46/chat"
all_merged_rds_dir <- file.path(output_dir, "rds")
pdf_dir <- file.path(output_dir, "pdf")
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

# Load the merged object
merged_obj <- readRDS(file.path(all_merged_rds_dir, "merged_totalNK_nkp46_dim25_res0.5.rds"))

# ------------------------- #
# Define Labels and Colors
# ------------------------- #
conditions <- c("nkp46_neg", "nkp46_pos", "totalNK")
condition_labels <- c("NKp46-", "NKp46+", "Total NK")
names(condition_labels) <- conditions

animals <- c("Animal25", "Animal26", "Animal27", "Animal28")

# Pastel tone color palette
condition_colors <- c("nkp46_neg" = "#e79088", "nkp46_pos" = "#91c58a", "totalNK" = "#79cadc")

# ------------------------- #
# UMAP by Condition
# ------------------------- #
cat("🎨 Generating combined UMAP by condition...\n")
umap_by_condition <- DimPlot(merged_obj,
    group.by = "condition", pt.size = 0.3,
    cols = condition_colors, label = FALSE
) +
    scale_color_manual(values = condition_colors, labels = condition_labels, name = "Condition") +
    ggtitle("UMAP by Condition") +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right")

ggsave(file.path(pdf_dir, "UMAP_by_condition_combined.png"),
    umap_by_condition,
    width = 10, height = 8, dpi = 600
)
cat("✅ UMAP by Condition grid saved with pastel color coding\n")


# ------------------------- #
# UMAP Grid: Animal × Condition
# ------------------------- #
plot_list <- list()

for (animal in animals) {
    for (cond in conditions) {
        cells <- WhichCells(merged_obj, expression = animal_id == animal & condition == cond)
        cat(sprintf("🔍 Found %d cells for %s x %s\n", length(cells), animal, cond))

        if (length(cells) == 0) {
            cat("⚠️ No cells found for", animal, cond, "- skipping\n")
            p <- ggplot() +
                theme_void()
        } else {
            p <- DimPlot(merged_obj, cells = cells, group.by = "condition", cols = condition_colors) +
                ggtitle(NULL) +
                theme_void() +
                theme(legend.position = "none")
        }

        plot_list[[paste(animal, cond, sep = "_")]] <- p
    }
}

# ------------------------- #
# Add Column Labels (Condition Names)
# ------------------------- #
column_titles <- lapply(condition_labels, function(label) {
    ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = label, size = 6, fontface = "bold") +
        theme_void()
})
top_row <- wrap_plots(column_titles, ncol = length(conditions))

# ------------------------- #
# Add Row Labels (Animal IDs)
# ------------------------- #
row_titles <- lapply(animals, function(animal) {
    ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = animal, angle = 90, size = 6, fontface = "bold") +
        theme_void()
})
left_col <- wrap_plots(row_titles, ncol = 1)

# ------------------------- #
# Assemble UMAP Grid
# ------------------------- #
# Arrange plots in a matrix: rows = animals, cols = conditions
umap_matrix <- lapply(animals, function(animal) {
    plots_row <- lapply(conditions, function(cond) {
        plot_list[[paste(animal, cond, sep = "_")]]
    })
    wrap_plots(plots_row, ncol = length(conditions))
})
umap_full <- wrap_plots(umap_matrix, ncol = 1)

# Add labels: top (column), left (row)
umap_labeled <- plot_grid(
    plot_grid(NULL, top_row, ncol = 2, rel_widths = c(0.12, 1)),
    plot_grid(left_col, umap_full, ncol = 2, rel_widths = c(0.12, 1)),
    nrow = 2, rel_heights = c(0.08, 1)
)

# ------------------------- #
# Export
# ------------------------- #
ggsave(file.path(pdf_dir, "UMAP_grid_Animal_by_Condition_pastel.png"),
    umap_labeled,
    width = 12, height = 16, dpi = 600
)

cat("✅ UMAP Animal × Condition grid saved with pastel color coding\n")
