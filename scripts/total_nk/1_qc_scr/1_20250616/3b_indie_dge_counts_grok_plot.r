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
        if (!file.exists(file)) {
          warning("File not found: ", file)
          next
        }

        df <- read_csv(file, show_col_types = FALSE)

        # Filter DEGs based on significance
        df_filtered <- df %>%
          filter(p_val_adj < 0.05, avg_log2FC > 0.25)

        deg_counts <- df_filtered %>%
          group_by(cluster) %>%
          summarise(n_DEGs = n(), .groups = "drop") %>%
          mutate(animal = animal,
                 dims = dims,
                 res = res,
                 method = method)

        deg_summary <- bind_rows(deg_summary, deg_counts)
      }
    }
  }
}

# --------------------------- #
# Format metadata & order clusters
# --------------------------- #
deg_summary <- deg_summary %>%
  mutate(
    cluster = paste0("C", cluster),
    dims = factor(dims, levels = dims_list),
    res = factor(res, levels = res_list),
    animal = factor(animal, levels = animals),
    method = factor(method, levels = methods)
  ) %>%
  # Validate that all expected dim-res combinations are present
  complete(dims, res, method, animal, cluster, fill = list(n_DEGs = 0))

# Order cluster labels numerically: C0, C1, ..., C20
unique_clusters <- sort(as.numeric(str_extract(unique(deg_summary$cluster), "\\d+")))
ordered_cluster_levels <- paste0("C", unique_clusters)
deg_summary$cluster <- factor(deg_summary$cluster, levels = ordered_cluster_levels)

# Diagnostic check
expected_combinations <- length(dims_list) * length(res_list) * length(methods) * length(animals) * length(unique_clusters)
actual_combinations <- nrow(distinct(deg_summary, dims, res, method, animal, cluster))
cat("Expected dim-res-method-animal-cluster combinations:", expected_combinations, "\n")
cat("Actual combinations in deg_summary:", actual_combinations, "\n")
if (actual_combinations < expected_combinations) {
  warning("Missing combinations detected. Check input CSV files.")
}

# --------------------------- #
# Plot & save for each method
# --------------------------- #
for (method in methods) {
  subset_df <- deg_summary %>% filter(method == method)

  # Create plot with data filtered per facet
  p <- ggplot(subset_df, aes(x = cluster, y = n_DEGs, color = animal)) +
    geom_point(size = 2) +
    # Filter data for each facet dynamically
    facet_grid(res ~ dims, scales = "free_x", space = "free_x") +
    labs(title = paste("Number of Significant DEGs per Cluster (", method, ")", sep = ""),
         x = "Cluster", y = "# DEGs (adj p < 0.05, logFC > 0.25)", color = "Animal") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 10))

  pdf_file <- file.path(pdf_output_dir, paste0("4_DEG_counts_per_cluster_grk_fixed_", method, ".pdf"))
  ggsave(pdf_file, plot = p, width = 16, height = 6, units = "in")

  cat("✅ DEG plot saved to:", pdf_file, "\n")
}