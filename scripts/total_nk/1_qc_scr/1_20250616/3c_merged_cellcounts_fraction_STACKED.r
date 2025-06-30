library(Seurat)
library(dplyr)
library(ggplot2)

# --------------------------- #
# Settings & Paths
# --------------------------- #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
rds_file <- file.path(output_dir, "rds", "merged_seurat_analysis_20250616.rds")
pdf_dir <- file.path(output_dir, "pdf", "4_merged_cluster_fraction_by_sample")
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

pdf_out <- file.path(pdf_dir, "cluster_fraction_scatterplot_merged_FIXED.pdf")

dims_list <- c(10, 15, 20, 25, 30)
res_list <- c(0.25, 0.5, 0.75)

# --------------------------- #
# Load Merged Seurat Object
# --------------------------- #
merged_rds <- readRDS(rds_file)
merged_obj <- merged_rds[["SCT"]]  # or [["LogNormalize"]]

# --------------------------- #
# Initialize data collector
# --------------------------- #
all_frac <- data.frame()

# --------------------------- #
# Loop over dims × resolution
# --------------------------- #
for (dims in dims_list) {
  for (res in res_list) {
    cluster_col <- paste0("integrated_snn_res.", res)

    if (!(cluster_col %in% colnames(merged_obj@meta.data))) {
      cat("❌ Skipping dims =", dims, "res =", res, " (no clustering column)\n")
      next
    }

    # Add cluster info
    merged_obj$cluster_tmp <- as.character(merged_obj[[cluster_col]][, 1])
    Idents(merged_obj) <- "cluster_tmp"

    # Table of cluster vs sample
    df <- table(Cluster = merged_obj$cluster_tmp, Animal = merged_obj$animal) %>%
      as.data.frame() %>%
      group_by(Cluster) %>%
      mutate(
        Fraction = Freq / sum(Freq),
        dims = dims,
        resolution = res,
        Cluster_numeric = as.numeric(as.character(Cluster))
      ) %>%
      ungroup() %>%
      mutate(
        cluster = paste0("C", Cluster_numeric),
        cluster_combined = paste0("D", dims, "_R", res, "_C", Cluster_numeric)
      )

    # Order factor levels per combination
    df$cluster_combined <- factor(df$cluster_combined,
                                  levels = unique(df$cluster_combined))

    all_frac <- bind_rows(all_frac, df)
  }
}

# --------------------------- #
# Plot stacked scatter plot
# --------------------------- #
p <- ggplot(all_frac, aes(x = cluster_combined, y = Fraction, color = Animal, group = Animal)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.2) +
  geom_line(position = position_dodge(width = 0.6), linewidth = 0.6, alpha = 0.7) +
  facet_wrap(~ resolution + dims, scales = "free_x", ncol = 1) +
  labs(
    title = "Fraction of Cells per Cluster by Sample (Merged SCT)",
    x = "Cluster (dims × res)", y = "Fraction", color = "Animal"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

# Save to PDF
ggsave(pdf_out, plot = p, width = 16, height = 12)

cat("✅ Fixed scatter plot saved to:", pdf_out, "\n")
