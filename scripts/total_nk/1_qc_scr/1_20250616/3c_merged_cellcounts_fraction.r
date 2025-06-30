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

dims_list <- c(10, 15, 20, 25, 30)
res_list <- c(0.25, 0.5, 0.75)

# --------------------------- #
# Load Merged Seurat Object
# --------------------------- #
merged_rds <- readRDS(rds_file)
merged_obj <- merged_rds[["SCT"]]  # or [["LogNormalize"]] if needed

# --------------------------- #
# Loop Over Dims × Resolution
# --------------------------- #
for (dims in dims_list) {
  for (res in res_list) {
    cat("📊 Plotting for dims =", dims, "| res =", res, "...\n")

    # Metadata column for this clustering
    cluster_col <- paste0("integrated_snn_res.", res)

    # Skip if the clustering column doesn't exist
    if (!(cluster_col %in% colnames(merged_obj@meta.data))) {
      cat("❌ Missing clustering column:", cluster_col, "\n")
      next
    }

    # Assign clustering and identities
    merged_obj$cluster_tmp <- as.character(merged_obj[[cluster_col]][, 1])
    Idents(merged_obj) <- "cluster_tmp"

    # Compute cell fraction per cluster × sample
    df_frac <- table(Cluster = merged_obj$cluster_tmp, Sample = merged_obj$animal) %>%
      as.data.frame() %>%
      group_by(Cluster) %>%
      mutate(Fraction = Freq / sum(Freq))

    # Format Cluster as C0, C1, C2...
    df_frac <- df_frac %>%
      mutate(Cluster_numeric = as.numeric(as.character(Cluster))) %>%
      arrange(Cluster_numeric) %>%
      mutate(Cluster_label = paste0("C", Cluster_numeric))

    # Ensure cluster label is a factor with ordered levels
    df_frac$Cluster_label <- factor(df_frac$Cluster_label, levels = unique(df_frac$Cluster_label))

    # Plot
    p <- ggplot(df_frac, aes(x = Cluster_label, y = Fraction, fill = Sample)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(
        title = paste("Cluster Composition by Sample\nDims =", dims, "| Resolution =", res),
        x = "Cluster", y = "Fraction of Cells", fill = "Animal"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10)
      )

    # Save PDF
    file_name <- paste0("cluster_fraction_dims", dims, "_res", res, ".pdf")
    ggsave(file.path(pdf_dir, file_name), p, width = 10, height = 5)

    cat("✅ Saved:", file_name, "\n")
  }
}
