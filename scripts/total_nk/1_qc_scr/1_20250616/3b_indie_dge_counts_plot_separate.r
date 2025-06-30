# library(ggplot2)
# library(dplyr)
# library(readr)
# library(stringr)

# # --------------------------- #
# # Paths
# # --------------------------- #
# output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
# pdf_output_dir <- file.path(output_dir, "pdf", "3_dge_heatmap")
# dir.create(pdf_output_dir, recursive = TRUE, showWarnings = FALSE)

# base_dir <- file.path(output_dir, "tsv", "individual")
# methods <- c("SCT", "LogNormalize")
# animals <- c("animal25", "animal26", "animal27", "animal28")
# dims_list <- c(10, 15, 20, 25, 30)
# res_list <- c(0.25, 0.5, 0.75)

# # --------------------------- #
# # Initialize storage
# # --------------------------- #
# deg_summary <- data.frame()

# # --------------------------- #
# # Loop through CSVs to collect DEG counts
# # --------------------------- #
# for (method in methods) {
#   method_dir <- file.path(base_dir, method)
#   for (animal in animals) {
#     for (dims in dims_list) {
#       for (res in res_list) {
#         file <- file.path(method_dir, paste0("markers_", animal, "_", method, "_dims", dims, "_res", res, ".csv"))
#         if (!file.exists(file)) next

#         df <- read_csv(file, show_col_types = FALSE)

#         # Filter DEGs based on significance
#         df_filtered <- df %>%
#           filter(p_val_adj < 0.05, avg_log2FC > 0.25)

#         # Count DEGs per cluster
#         deg_counts <- df_filtered %>%
#           group_by(cluster) %>%
#           summarise(n_DEGs = n(), .groups = "drop") %>%
#           mutate(
#             animal = animal,
#             dims = dims,
#             resolution = res,
#             method = method,
#             cluster = paste0("C", cluster)
#           )

#         deg_summary <- bind_rows(deg_summary, deg_counts)
#       }
#     }
#   }
# }

# # --------------------------- #
# # Factor formatting
# # --------------------------- #
# deg_summary <- deg_summary %>%
#   mutate(
#     cluster = factor(cluster, levels = paste0("C", 0:30)),  # support up to C30
#     dims = factor(dims, levels = dims_list),
#     resolution = factor(resolution, levels = as.character(res_list)),
#     animal = factor(animal, levels = animals),
#     method = factor(method, levels = methods)
#   )

# # --------------------------- #
# # Plot & save for each method
# # --------------------------- #
# for (method in methods) {
#   subset_df <- deg_summary %>% filter(method == method)

#   p <- ggplot(subset_df, aes(x = cluster, y = n_DEGs, color = animal, group = animal)) +
#     geom_point(position = position_dodge(width = 0.6), size = 2) +
#     geom_line(position = position_dodge(width = 0.6), linewidth = 0.8, alpha = 0.7) +
#     facet_grid(resolution ~ dims, scales = "free_x", space = "free_x") +
#     labs(
#       title = paste("Number of Significant DEGs per Cluster (", method, ")", sep = ""),
#       x = "Cluster", y = "# DEGs (adj p < 0.05, logFC > 0.25)", color = "Animal"
#     ) +
#     theme_minimal() +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
#       strip.text = element_text(size = 10),
#       panel.grid.minor = element_blank()
#     )

#   pdf_file <- file.path(pdf_output_dir, paste0("4_DEG_counts_per_cluster_", method, "_fixed.pdf"))
#   ggsave(pdf_file, plot = p, width = 16, height = 6, units = "in")

#   cat("✅ Fixed DEG scatter plot saved to:", pdf_file, "\n")
# }


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

        # Filter DEGs based on significance
        df_filtered <- df %>%
          filter(p_val_adj < 0.05, avg_log2FC > 0.25)

        # Count DEGs per cluster
        deg_counts <- df_filtered %>%
          group_by(cluster) %>%
          summarise(n_DEGs = n(), .groups = "drop") %>%
          mutate(
            animal = animal,
            dims = dims,
            resolution = res,
            method = method,
            cluster = paste0("C", cluster)
          )

        deg_summary <- bind_rows(deg_summary, deg_counts)
      }
    }
  }
}

# --------------------------- #
# Plot (correct cluster scoping)
# --------------------------- #
deg_summary <- deg_summary %>%
  mutate(
    dims = factor(dims, levels = dims_list),
    resolution = factor(resolution, levels = as.character(res_list)),
    cluster = factor(cluster, levels = paste0("C", 0:30)),
    animal = factor(animal, levels = animals),
    method = factor(method, levels = methods)
  )

for (method in methods) {
  subset_df <- deg_summary %>% filter(method == method)

  # Plot correctly scoped clusters
  p <- ggplot(subset_df, aes(x = cluster, y = n_DEGs, color = animal, group = animal)) +
    geom_point(position = position_dodge(width = 0.6), size = 2) +
    geom_line(position = position_dodge(width = 0.6), linewidth = 0.8, alpha = 0.7) +
    facet_grid(rows = vars(resolution), cols = vars(dims), scales = "free_x", space = "free_x") +
    labs(
      title = paste("Number of Significant DEGs per Cluster (", method, ")", sep = ""),
      x = "Cluster (within dims × res)", y = "# DEGs (adj p < 0.05, logFC > 0.25)", color = "Animal"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 10)
    )

  pdf_file <- file.path(pdf_output_dir, paste0("4_DEG_counts_per_cluster_", method, "_final_fixed.pdf"))
  ggsave(pdf_file, plot = p, width = 16, height = 8, units = "in")

  cat("✅ Final fixed DEG scatter plot saved to:", pdf_file, "\n")
}
