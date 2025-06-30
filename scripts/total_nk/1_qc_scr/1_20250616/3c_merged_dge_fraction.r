library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

# --------------------------- #
# Paths and settings
# --------------------------- #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
rds_file <- file.path(output_dir, "rds", "merged_seurat_analysis_20250616.rds")

pdf_dir <- file.path(output_dir, "pdf", "4_cluster_fraction_by_sample")
tsv_dir <- file.path(output_dir, "tsv", "merged")
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tsv_dir, recursive = TRUE, showWarnings = FALSE)

pdf_file <- file.path(pdf_dir, "combined_cluster_fractions.pdf")
csv_summary <- file.path(tsv_dir, "cluster_fraction_summary.csv")
csv_flags <- file.path(tsv_dir, "low_fraction_flags.csv")

dims_list <- c(10, 15, 20, 25, 30)
res_list <- c(0.25, 0.5, 0.75)

# --------------------------- #
# Load Seurat Object
# --------------------------- #
merged_rds <- readRDS(rds_file)
merged_obj <- merged_rds[["SCT"]]  # or [["LogNormalize"]] if needed

# --------------------------- #
# Init collectors
# --------------------------- #
all_plots <- list()
all_summary <- data.frame()
low_fraction_flags <- data.frame()

# --------------------------- #
# Loop over dims × res
# --------------------------- #
for (dims in dims_list) {
  for (res in res_list) {
    cat("📊 Processing dims =", dims, "res =", res, "\n")

    cluster_col <- paste0("integrated_snn_res.", res)
    if (!(cluster_col %in% colnames(merged_obj@meta.data))) {
      cat("❌ Skipped — missing column:", cluster_col, "\n")
      next
    }

    # Assign temp identity
    merged_obj$cluster_tmp <- as.character(merged_obj[[cluster_col]][, 1])
    Idents(merged_obj) <- "cluster_tmp"

    # Compute fraction table
    tbl <- as.data.frame(table(Cluster = merged_obj$cluster_tmp,
                               Animal = merged_obj$animal)) %>%
      group_by(Cluster) %>%
      mutate(Fraction = Freq / sum(Freq),
             dims = dims,
             resolution = res)

    # Sort clusters numerically
    tbl <- tbl %>%
      mutate(Cluster = factor(Cluster, levels = sort(as.numeric(as.character(unique(Cluster))))))

    # Add to collectors
    all_summary <- bind_rows(all_summary, tbl)
    low_fraction_flags <- bind_rows(low_fraction_flags, filter(tbl, Fraction < 0.05))

    # Plot
    p <- ggplot(tbl, aes(x = Cluster, y = Fraction, fill = Animal)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(title = paste("Cluster Composition | dims =", dims, "| res =", res),
           x = "Cluster", y = "Fraction of Cells", fill = "Animal") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    all_plots[[paste0("dims", dims, "_res", res)]] <- p
  }
}

# --------------------------- #
# Save multi-page PDF
# --------------------------- #
cat("💾 Saving combined PDF:", pdf_file, "\n")
pdf(pdf_file, width = 11, height = 6)
for (name in names(all_plots)) {
  print(all_plots[[name]])
}
dev.off()

# --------------------------- #
# Save CSV tables
# --------------------------- #
write.csv(all_summary, csv_summary, row.names = FALSE)
cat("✅ Saved cell fraction summary to:", csv_summary, "\n")

if (nrow(low_fraction_flags) > 0) {
  write.csv(low_fraction_flags, csv_flags, row.names = FALSE)
  cat("⚠️  Saved low-fraction cluster flags to:", csv_flags, "\n")
} else {
  cat("✅ No low-fraction clusters (<5%) found.\n")
}
