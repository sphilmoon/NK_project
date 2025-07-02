library(Seurat)
library(dplyr)
library(ggplot2)

# --------------------------- #
# Settings & Paths
# --------------------------- #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
rds_file <- file.path(output_dir, "rds", "merged_seurat_analysis_20250616.rds")
pdf_dir <- file.path(output_dir, "pdf", "3_merged_DEGs_outs")
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
pdf_out <- file.path(pdf_dir, "cluster_fraction_scatterplot_merged_fixed.pdf")

dims_list <- c(10, 15, 20, 25, 30)
res_list <- c(0.25, 0.5, 0.75)

# --------------------------- #
# Load Merged Seurat Object
# --------------------------- #
merged_rds <- readRDS(rds_file)
merged_obj <- merged_rds[["SCT"]]  # or [["LogNormalize"]] if needed

# --------------------------- #
# Initialize collector
# --------------------------- #
all_frac <- data.frame()

# --------------------------- #
# Loop over dims × resolution
# --------------------------- #
for (dims in dims_list) {
  for (res in res_list) {
    cat("📊 Processing dims =", dims, "| res =", res, "\n")

    # Re-run dimensional reduction and clustering
    merged_obj <- FindNeighbors(merged_obj, dims = 1:dims, verbose = FALSE)
    merged_obj <- FindClusters(merged_obj, resolution = res, verbose = FALSE)

    # Save unique clustering label
    cluster_col <- paste0("integrated_snn_res.d", dims, "_r", res)
    merged_obj[[cluster_col]] <- merged_obj$seurat_clusters

    # Assign identities for downstream table building
    merged_obj$cluster_tmp <- as.character(merged_obj[[cluster_col]][, 1])
    Idents(merged_obj) <- "cluster_tmp"

    # Cell count per cluster/sample → compute fraction
    df <- table(Cluster = merged_obj$cluster_tmp, Animal = merged_obj$animal) %>%
      as.data.frame() %>%
      group_by(Cluster) %>%
      mutate(
        Fraction = Freq / sum(Freq),
        dims = dims,
        resolution = res
      ) %>%
      ungroup()

    # Format cluster as C0, C1, ...
    df <- df %>%
      mutate(
        Cluster_numeric = as.numeric(as.character(Cluster)),
        cluster = paste0("C", Cluster_numeric),
        cluster = factor(cluster, levels = paste0("C", sort(unique(Cluster_numeric))))
      )

    # Append
    all_frac <- bind_rows(all_frac, df)
  }
}

# --------------------------- #
# Plot: stacked scatter plot
# --------------------------- #
p <- ggplot(all_frac, aes(x = cluster, y = Fraction, color = Animal, group = Animal)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.2) +
  geom_line(position = position_dodge(width = 0.6), linewidth = 0.6, alpha = 0.7) +
  facet_grid(resolution ~ dims, scales = "free_x", space = "free_x") +
  labs(
    title = "Fraction of Cells per Cluster by Sample (Merged SCT)",
    x = "Cluster", y = "Fraction", color = "Animal"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

# --------------------------- #
# Save output
# --------------------------- #
ggsave(pdf_out, plot = p, width = 16, height = 6)
cat("✅ Saved scatter plot to:", pdf_out, "\n")