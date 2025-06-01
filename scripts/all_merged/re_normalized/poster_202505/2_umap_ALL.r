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
    scale_color_manual(values = condition_colors, labels = condition_labels) +
    ggtitle(NULL) +
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
        cells <- WhichCells(merged_obj, expression = sample == animal & condition == cond)
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

# ------------------------- #
# UMAP Grid: Animal × Condition for NCR1 Gene Expression
# ------------------------- #
plot_list_ncr1 <- list()

for (animal in animals) {
    for (cond in conditions) {
        cells <- WhichCells(merged_obj, expression = sample == animal & condition == cond)
        cat(sprintf("🔍 Found %d cells for %s x %s (NCR1 plot)\n", length(cells), animal, cond))

        if (length(cells) == 0) {
            cat("⚠️ No cells found for", animal, cond, "- skipping\n")
            p <- ggplot() +
                theme_void()
        } else {
            p <- FeaturePlot(merged_obj, features = "NCR1", cells = cells, pt.size = 0.3) +
                ggtitle(NULL) +
                theme_void() +
                theme(legend.position = "none")
        }

        plot_list_ncr1[[paste(animal, cond, sep = "_")]] <- p
    }
}

# ------------------------- #
# UMAP by Condition for NCR1 Gene Expression
# ------------------------- #
cat("🎨 Generating UMAP by condition for NCR1 expression...\n")
ncr1_plots <- list()

for (cond in conditions) {
    cells <- WhichCells(merged_obj, expression = condition == cond)
    cat(sprintf("🔍 Found %d cells for %s (NCR1 by condition plot)\n", length(cells), cond))

    if (length(cells) == 0) {
        cat("⚠️ No cells found for", cond, "- skipping\n")
        p <- ggplot() +
            theme_void()
    } else {
        p <- FeaturePlot(merged_obj, features = "NCR1", cells = cells, pt.size = 0.3) +
            ggtitle(NULL) +
            theme(
                plot.title = element_text(hjust = 0.5, size = 10),
                axis.text = element_blank(),
                axis.title = element_blank(),
                legend.position = "none" # Remove individual legends
            )
    }
    ncr1_plots[[cond]] <- p
}

# ------------------------- #
# Extract Legend for NCR1 Expression
# ------------------------- #
# Determine expression range for NCR1 across all cells using the SCT assay
ncr1_data <- GetAssayData(merged_obj, assay = "SCT", layer = "data")
if ("NCR1" %in% rownames(ncr1_data)) {
    expr_range <- range(ncr1_data["NCR1", ], na.rm = TRUE)
} else {
    stop("NCR1 not found in the SCT assay data layer. Please check the gene name or assay.")
}

# Create a plot to extract the legend
legend_plot <- FeaturePlot(
    merged_obj,
    features = "NCR1",
    pt.size = 0.3
) +
    scale_color_gradientn(
        colors = c("lightgrey", "blue"),
        name = "NCR1",
        limits = expr_range
    ) +
    theme(
        legend.position = "right",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1.5, "cm"),
        legend.key.width = unit(0.5, "cm")
    )

legend <- cowplot::get_legend(legend_plot)

# ------------------------- #
# Assemble UMAP by Condition for NCR1 with a Single Legend
# ------------------------- #
# Arrange plots with condition labels on top
top_row_ncr1 <- wrap_plots(column_titles, ncol = length(conditions))
ncr1_combined <- wrap_plots(ncr1_plots, ncol = length(conditions))
ncr1_labeled <- plot_grid(
    top_row_ncr1,
    ncr1_combined,
    nrow = 2,
    rel_heights = c(0.1, 1)
)

# Add the single legend to the right
ncr1_final <- plot_grid(
    ncr1_labeled,
    legend,
    ncol = 2,
    rel_widths = c(1, 0.15)
)

# Export NCR1 UMAP by Condition with a Single Legend
ggsave(file.path(pdf_dir, "UMAP_by_Condition_NCR1_expression.png"),
    ncr1_final,
    width = 16, height = 5, dpi = 600
)

cat("✅ UMAP by Condition for NCR1 expression saved with a single legend\n")

# ------------------------- #
# Assemble UMAP Grid for NCR1 Expression with Legend
# ------------------------- #
# Arrange plots in a matrix: rows = animals, cols = conditions
umap_matrix_ncr1 <- lapply(animals, function(animal) {
    plots_row <- lapply(conditions, function(cond) {
        plot_list_ncr1[[paste(animal, cond, sep = "_")]]
    })
    wrap_plots(plots_row, ncol = length(conditions))
})
umap_full_ncr1 <- wrap_plots(umap_matrix_ncr1, ncol = 1)

# Add labels: top (column), left (row)
umap_labeled_ncr1 <- plot_grid(
    plot_grid(NULL, top_row, ncol = 2, rel_widths = c(0.12, 1)),
    plot_grid(left_col, umap_full_ncr1, ncol = 2, rel_widths = c(0.12, 1)),
    nrow = 2, rel_heights = c(0.08, 1)
)

# Add the legend to the right of the grid
umap_final_ncr1 <- plot_grid(
    umap_labeled_ncr1,
    legend,
    ncol = 2,
    rel_widths = c(1, 0.15)
)

# ------------------------- #
# Export NCR1 UMAP Grid with Legend
# ------------------------- #
ggsave(file.path(pdf_dir, "UMAP_grid_Animal_by_Condition_NCR1_expression.png"),
    umap_final_ncr1,
    width = 14, height = 16, dpi = 600
)

cat("✅ UMAP Animal × Condition grid for NCR1 expression saved with legend on the right\n")
