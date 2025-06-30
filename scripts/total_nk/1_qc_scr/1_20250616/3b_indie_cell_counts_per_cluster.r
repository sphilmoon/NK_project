library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# --------------------------- #
# Paths
# --------------------------- #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
pdf_output_dir <- file.path(output_dir, "pdf", "3_dge_heatmap")
dir.create(pdf_output_dir, recursive = TRUE, showWarnings = FALSE)

# Load cluster summary file
summary_file <- file.path(output_dir, "tsv", "indie_cluster_summary_sct_vs_log_dge.csv")
cluster_data <- read_csv(summary_file, show_col_types = FALSE)

# --------------------------- #
# Preprocessing
# --------------------------- #
# Make sure cluster is factor with ordered levels
cluster_data <- cluster_data %>%
  mutate(
    cluster = paste0("C", cluster),
    dims = factor(dims),
    res = factor(res),
    animal = factor(animal, levels = c("animal25", "animal26", "animal27", "animal28")),
    method = factor(method, levels = c("SCT", "LogNormalize"))
  )

# Order cluster levels numerically (C0, C1, ..., C20)
unique_clusters <- sort(as.numeric(str_extract(unique(cluster_data$cluster), "\\d+")))
ordered_cluster_levels <- paste0("C", unique_clusters)
cluster_data$cluster <- factor(cluster_data$cluster, levels = ordered_cluster_levels)

# --------------------------- #
# Plot & save for each method
# --------------------------- #
for (method in c("SCT", "LogNormalize")) {
  df_plot <- cluster_data %>% filter(method == !!method)

  p <- ggplot(df_plot, aes(x = cluster, y = cell_count, group = animal, color = animal)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.8, alpha = 0.7) +
    facet_grid(res ~ dims, scales = "free_x", space = "free_x") +
    labs(title = paste("Cell Count per Cluster (", method, ")", sep = ""),
         x = "Cluster", y = "Cell Count", color = "Animal") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 10))

  pdf_file <- file.path(pdf_output_dir, paste0("5_cell_count_per_cluster_", method, ".pdf"))
  ggsave(pdf_file, plot = p, width = 16, height = 6, units = "in")

  cat("✅ Cell count plot saved to:", pdf_file, "\n")
}
