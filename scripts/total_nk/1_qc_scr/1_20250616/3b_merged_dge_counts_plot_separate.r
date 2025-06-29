library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# --------------------------- #
# Paths
# --------------------------- #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
merged_dir <- file.path(output_dir, "tsv", "merged", "SCT")
pdf_output_dir <- file.path(output_dir, "pdf", "3_dge_heatmap")
dir.create(pdf_output_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------- #
# Initialize storage
# --------------------------- #
deg_summary <- data.frame()

# --------------------------- #
# Define parameters
# --------------------------- #
dims_list <- c(10, 15, 20, 25, 30)
res_list <- c(0.25, 0.5, 0.75)
method <- "SCT"  # Fixed method for merged data

# --------------------------- #
# Loop through CSVs to collect DEG counts
# --------------------------- #
for (dims in dims_list) {
  for (res in res_list) {
    file <- file.path(merged_dir, paste0("new_markers_merged_SCT_dims", dims, "_res", res, ".csv"))
    if (!file.exists(file)) {
      warning("File not found: ", file)
      next
    }

    df <- read_csv(file, show_col_types = FALSE)

    # Filter DEGs based on significance (adjust p-value < 0.05 and log2FC > 0.25)
    df_filtered <- df %>%
      filter(p_val_adj < 0.05, avg_log2FC > 0.25)

    deg_counts <- df_filtered %>%
      group_by(cluster) %>%
      summarise(n_DEGs = n(), .groups = "drop") %>%
      mutate(dims = dims,
             resolution = res,
             method = method,
             animal = "merged")

    deg_summary <- bind_rows(deg_summary, deg_counts)
  }
}

# --------------------------- #
# Format metadata & order clusters
# --------------------------- #
deg_summary <- deg_summary %>%
  mutate(
    cluster = paste0("C", cluster),
    dims = factor(dims),
    resolution = factor(resolution),
    method = factor(method),
    animal = factor(animal)
  )

# Order cluster labels numerically: C0, C1, ..., C20
unique_clusters <- sort(as.numeric(str_extract(unique(deg_summary$cluster), "\\d+")))
ordered_cluster_levels <- paste0("C", unique_clusters)
deg_summary$cluster <- factor(deg_summary$cluster, levels = ordered_cluster_levels)

# --------------------------- #
# Plot & save
# --------------------------- #
subset_df <- deg_summary  # Single method, no need for loop

p <- ggplot(subset_df, aes(x = cluster, y = n_DEGs, group = dims, color = dims)) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  facet_grid(resolution ~ dims, scales = "free_x", space = "free_x") +
  labs(title = "Number of Significant DEGs per Cluster (Merged SCT)",
       x = "Cluster", y = "# DEGs (adj p < 0.05, logFC > 0.25)", color = "Dims") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 10))

pdf_file <- file.path(pdf_output_dir, "4_DEG_counts_per_cluster_merged_SCT.pdf")
ggsave(pdf_file, plot = p, width = 16, height = 6, units = "in")

cat("✅ DEG plot saved to:", pdf_file, "\n")