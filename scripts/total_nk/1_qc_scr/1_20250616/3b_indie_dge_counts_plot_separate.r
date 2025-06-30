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

base_dir <- file.path(output_dir, "tsv", "individual")
methods <- c("SCT", "LogNormalize")
animals <- c("animal25", "animal26", "animal27", "animal28")
dims_list <- c(10, 15, 20, 25, 30)
res_list <- c(0.25, 0.5, 0.75)

# --------------------------- #
# Initialize storage
# --------------------------- #
deg_summary <- data.frame()

# --------------------------- #
# Loop through CSVs to collect DEG counts
# --------------------------- #
for (method in methods) {
  method_dir <- file.path(base_dir, method)
  for (animal in animals) {
    for (dims in dims_list) {
      for (res in res_list) {
        file <- file.path(method_dir, paste0("markers_", animal, "_", method, "_dims", dims, "_res", res, ".csv"))
        if (!file.exists(file)) next

        df <- read_csv(file, show_col_types = FALSE)

        # Filter DEGs based on significance (adjust p-value < 0.05 and log2FC > 0.25)
        df_filtered <- df %>%
          filter(p_val_adj < 0.05, avg_log2FC > 0.25)

        # Add unique cluster ID per dim/res
        deg_counts <- df_filtered %>%
          group_by(cluster) %>%
          summarise(n_DEGs = n(), .groups = "drop") %>%
          mutate(
            dims = dims,
            res = res,
            method = method,
            animal = animal,
            cluster_label = paste0("D", dims, "_R", res, "_C", cluster)
          )

        deg_summary <- bind_rows(deg_summary, deg_counts)
      }
    }
  }
}

# --------------------------- #
# Prepare for plotting
# --------------------------- #
deg_summary <- deg_summary %>%
  mutate(
    dims = factor(dims, levels = dims_list),
    res = factor(res, levels = res_list),
    animal = factor(animal, levels = animals),
    cluster_label = factor(cluster_label, levels = unique(cluster_label))
  )

# --------------------------- #
# Plot & save for each method
# --------------------------- #
for (method in methods) {
  subset_df <- deg_summary %>% filter(method == !!method)

  p <- ggplot(subset_df, aes(x = cluster_label, y = n_DEGs, group = animal, color = animal)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.8, alpha = 0.7) +
    facet_wrap(~ interaction(dims, res), scales = "free_x", ncol = 1) +
    labs(
      title = paste("Number of Significant DEGs per Cluster (", method, ")", sep = ""),
      x = "Cluster (Dims × Res)", y = "# DEGs (adj p < 0.05, logFC > 0.25)", color = "Animal"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 10)
    )

  pdf_file <- file.path(pdf_output_dir, paste0("4_DEG_counts_per_cluster_", method, "_FIXED.pdf"))
  ggsave(pdf_file, plot = p, width = 16, height = 10, units = "in")

  cat("✅ Fixed DEG plot saved to:", pdf_file, "\n")
}
