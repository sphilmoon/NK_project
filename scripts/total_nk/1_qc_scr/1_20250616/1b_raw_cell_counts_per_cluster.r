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

# --------------------------- #
# Load CSV
# --------------------------- #
csv_file <- file.path(output_dir, "tsv", "indie_cluster_summary_sct_vs_log_dge.csv")
df <- read_csv(csv_file, show_col_types = FALSE)

# --------------------------- #
# Clean & Prepare
# --------------------------- #
df_clean <- df %>%
  mutate(
    cluster_id = paste0(animal, "_C", cluster),         # Unique x-axis value
    dims = factor(dims),
    resolution = factor(resolution),
    method = factor(method),
    animal = factor(animal),
    cluster = factor(paste0("C", cluster))
  )

# Reorder cluster_id factor levels numerically by cluster
df_clean$cluster_num <- as.numeric(str_extract(df_clean$cluster, "\\d+"))
df_clean <- df_clean %>%
  arrange(animal, cluster_num) %>%
  mutate(cluster_id = factor(cluster_id, levels = unique(cluster_id)))  # Keep ordering

# --------------------------- #
# Plot and Save
# --------------------------- #
methods <- c("SCT", "LogNormalize")

for (m in methods) {
  plot_df <- df_clean %>% filter(method == m)

  p <- ggplot(plot_df, aes(x = cluster_id, y = cell_count, color = animal, group = animal)) +
    geom_point(size = 2.2) +
    facet_grid(resolution ~ dims, scales = "free_x", space = "free_x") +
    labs(
      title = paste("Cell Count per Cluster (", m, ")"),
      x = "Cluster (Animal-wise)", y = "Cell Count", color = "Animal"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      strip.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )

  pdf_file <- file.path(pdf_output_dir, paste0("3_cell_count_per_cluster_Fixed_", m, ".pdf"))
  ggsave(pdf_file, p, width = 18, height = 6, units = "in")

  cat("✅ Saved:", pdf_file, "\n")
}
