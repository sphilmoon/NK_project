library(Seurat)
library(dplyr)
library(ggplot2)

# --------------------------- #
# Settings & Paths
# --------------------------- #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
rds_file <- file.path(output_dir, "rds", "merged_seurat_analysis_20250616.rds")

# Define and create PDF directory
pdf_dir <- file.path(output_dir, "pdf", "3_merged_DEGs_outs")
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
# Define output PDF file path
pdf_out <- file.path(pdf_dir, "cluster_fraction_scatterplot_merged_TEST.pdf")

# Define and create TSV directory
tsv_dir <- file.path(output_dir, "tsv", "merged", "3_merged_DEGs_outs")
dir.create(tsv_dir, recursive = TRUE, showWarnings = FALSE)
# Define output CSV files
csv_summary <- file.path(tsv_dir, "cluster_fraction_summary.csv")
csv_flags <- file.path(tsv_dir, "low_fraction_flags.csv")

# Define dimensions and resolutions to test
dims_list <- c(10, 15, 20, 25, 30)
res_list <- c(0.25, 0.5, 0.75)

# --------------------------- #
# Load Merged Seurat Object
# --------------------------- #
merged_rds <- readRDS(rds_file)
merged_obj <- merged_rds[["SCT"]]

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

    # Re-run clustering
    merged_obj <- FindNeighbors(merged_obj, dims = 1:dims, verbose = FALSE)
    merged_obj <- FindClusters(merged_obj, resolution = res, verbose = FALSE)

    # Create unique label
    cluster_col <- paste0("integrated_snn_res.d", dims, "_r", res)
    merged_obj[[cluster_col]] <- merged_obj$seurat_clusters
    merged_obj$cluster_tmp <- as.character(merged_obj[[cluster_col]][, 1])
    Idents(merged_obj) <- "cluster_tmp"

    # Cell count per cluster per animal
    df <- table(Cluster = merged_obj$cluster_tmp, Animal = merged_obj$animal) %>%
      as.data.frame() %>%
      group_by(Cluster) %>%
      mutate(
        Count = Freq,
        Fraction = Count / sum(Count),
        dims = dims,
        resolution = res
      ) %>%
      ungroup()

    # Format cluster labels
    df <- df %>%
      mutate(
        Cluster_numeric = as.numeric(as.character(Cluster)),
        cluster = paste0("C", Cluster_numeric),
        cluster = factor(cluster, levels = paste0("C", sort(unique(Cluster_numeric))))
      )

    all_frac <- bind_rows(all_frac, df)
  }
}

# --------------------------- #
# Save summary & flags
# --------------------------- #
cluster_summary <- all_frac %>%
  select(Cluster, Animal, Count, Fraction, dims, resolution)

write.csv(cluster_summary, csv_summary, row.names = FALSE)
cat("✅ Saved cluster fraction summary to:", csv_summary, "\n")

low_fraction_flags <- cluster_summary %>% filter(Fraction < 0.05)
write.csv(low_fraction_flags, csv_flags, row.names = FALSE)
cat("⚠️ Saved low fraction warnings to:", csv_flags, "\n")

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

# Save plot
ggsave(pdf_out, plot = p, width = 16, height = 6)
cat("✅ Saved scatter plot to:", pdf_out, "\n")