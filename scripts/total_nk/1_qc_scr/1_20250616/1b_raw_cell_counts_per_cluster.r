library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# --------------------------- #
# Paths
# --------------------------- #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
pdf_output_dir <- file.path(output_dir, "pdf")
dir.create(pdf_output_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------- #
# Load data
# --------------------------- #
base_dir <- file.path(output_dir, "tsv", "indie_cluster_summary_sct_vs_log_dge.csv")
if (!file.exists(base_dir)) stop("File not found at: ", base_dir)
df <- read_csv(base_dir, show_col_types = FALSE)

# --------------------------- #
# Filter and prepare data
# --------------------------- #
cell_summary <- df %>%
  mutate(
    cluster = paste0("C", cluster),  # Format cluster as C0, C1, etc.
    dims = factor(dims),
    resolution = factor(resolution),
    animal = factor(animal),
    method = factor(method)
  )

# Order cluster labels numerically: C0, C1, ..., C20
unique_clusters <- sort(as.numeric(str_extract(unique(cell_summary$cluster), "\\d+")))
ordered_cluster_levels <- paste0("C", unique_clusters)
cell_summary$cluster <- factor(cell_summary$cluster, levels = ordered_cluster_levels)

# --------------------------- #
# Plot & save for each method
# --------------------------- #
methods <- c("SCT", "LogNormalize")
for (method in methods) {
  subset_df <- cell_summary %>% filter(method == method)  # Removed unnecessary !! operator

  p <- ggplot(subset_df, aes(x = cluster, y = cell_count, group = animal, color = animal)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.8, alpha = 0.7) +
    facet_grid(resolution ~ dims, scales = "free_x", space = "free_x") +
    labs(title = paste("Cell Count per Cluster (", method, ")", sep = ""),
         x = "Cluster", y = "Cell Count", color = "Animal") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 10))

  pdf_file <- file.path(pdf_output_dir, paste0("3_cell_count_per_cluster_", method, ".pdf"))
  ggsave(pdf_file, plot = p, width = 16, height = 6, units = "in")

  cat("✅ Cell count plot saved to:", pdf_file, "\n")
}