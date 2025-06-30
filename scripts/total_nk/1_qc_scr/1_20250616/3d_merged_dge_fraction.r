library(Seurat)
library(dplyr)
library(ggplot2)

# --------------------------- #
# Settings & Paths
# --------------------------- #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
rds_file <- file.path(output_dir, "rds", "merged_cluster_summary_sct_vs_log_dge.rds")
pdf_dir <- file.path(output_dir, "pdf", "4_merged_cluster_fraction_by_sample")
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)



