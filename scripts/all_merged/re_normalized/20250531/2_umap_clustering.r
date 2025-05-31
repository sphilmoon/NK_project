# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(ggplot2)
library(cowplot)
library(stringr)
library(cluster) # For silhouette score
library(dplyr) # For data frame manipulation

# ------------------------- #
# Define Input and Output Paths
# ------------------------- #
output_dir <- "/home/outputs/all_merged_TotalNK_nkp46/20250531"
all_merged_rds_dir <- file.path(output_dir, "rds")
pdf_dir <- file.path(output_dir, "pdf")

# Ensure directories exist
dir.create(all_merged_rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------- #
# Load Saved Seurat Object
# ------------------------- #
cat("📦 Loading saved Seurat object...\n")
# Correctly load the RDS file
# Note: The file was saved in /home/outputs/all_merged_TotalNK_nkp46/chat/rds/,
# but we're saving plots in /home/outputs/all_merged_TotalNK_nkp46/20250531/,
# so we need to adjust the path to point to the correct directory
input_rds_dir <- "/home/outputs/all_merged_TotalNK_nkp46/chat/rds"
merged_obj <- readRDS(file.path(input_rds_dir, "merged_totalNK_nkp46_dim25_res0.5.rds"))
cat("✅ Loaded merged object\n")

# Verify clustering metadata
cat("🔍 Checking clustering metadata...\n")
if (!"seurat_clusters" %in% colnames(merged_obj@meta.data)) {
    stop("❌ Cluster labels not found in the Seurat object. Please ensure clustering was performed.")
}
cat("✅ Found cluster labels:", levels(merged_obj$seurat_clusters), "\n")

# ------------------------- #
# UMAP by Condition with Cluster Numbers
# ------------------------- #
cat("🎨 Generating combined UMAP by condition with cluster labels...\n")
umap_by_condition <- DimPlot(
    merged_obj,
    group.by = "condition",
    label = TRUE, # Add cluster numbers
    label.size = 5, # Size of the cluster labels
    pt.size = 0.3,
    repel = TRUE # Prevent label overlap
) +
    ggtitle("UMAP by Condition (Colored by Condition, Labeled by Cluster)") +
    theme(plot.title = element_text(hjust = 0.5, size = 14))

ggsave(
    file.path(pdf_dir, "UMAP_by_condition_combined_with_clusters_20250531.pdf"),
    umap_by_condition,
    width = 10, height = 8, dpi = 600
)
cat("✅ UMAP by Condition with cluster labels saved\n")

# ------------------------- #
# UMAP Grid: Animal x Condition with Cluster Numbers
# ------------------------- #
cat("🎨 Generating UMAP grid by animal x condition with cluster labels...\n")
conditions <- c("nkp46_neg", "nkp46_pos", "totalNK")
animals <- c("Animal25", "Animal26", "Animal27", "Animal28")
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
            p <- DimPlot(
                merged_obj,
                cells = cells,
                group.by = "condition",
                label = TRUE, # Add cluster numbers
                label.size = 4, # Smaller label size for grid
                pt.size = 0.3,
                repel = TRUE # Prevent label overlap
            ) +
                ggtitle(NULL) +
                theme_void() +
                theme(legend.position = "none")
        }
        plot_list[[paste(animal, cond, sep = "_")]] <- p
    }
}

# Add column titles (condition names)
condition_labels <- c("NKp46-", "NKp46+", "Total NK")
names(condition_labels) <- conditions
column_titles <- lapply(condition_labels, function(label) {
    ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = label, size = 6, fontface = "bold") +
        theme_void()
})
top_row <- wrap_plots(column_titles, ncol = length(conditions))

# Add row titles (animal IDs)
row_titles <- lapply(animals, function(animal) {
    ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = animal, angle = 90, size = 6, fontface = "bold") +
        theme_void()
})
left_col <- wrap_plots(row_titles, ncol = 1)

# Assemble UMAP grid
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

# Export UMAP grid
ggsave(
    file.path(pdf_dir, "UMAP_animal_condition_grid_with_clusters_20250531.pdf"),
    umap_labeled,
    width = 12, height = 16, dpi = 600
)

cat("✅ UMAP grid by Animal x Condition with cluster labels saved\n")
cat("✅ All UMAP visualizations complete\n")
