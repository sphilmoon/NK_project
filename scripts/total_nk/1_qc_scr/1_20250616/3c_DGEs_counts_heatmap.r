# Load required libraries
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)

# ------------------------ #
# Settings
# ------------------------ #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
csv_path <- file.path(output_dir, "tsv", "merged", "3_merged_DEGs_outs", "cluster_DEG_counts_per_cluster_per_sample.csv")
heatmap_pdf <- file.path(output_dir, "pdf", "3_merged_DEGs_outs", "DEG_counts_heatmap_total.pdf")
summary_csv <- file.path(output_dir, "tsv", "merged", "3_merged_DEGs_outs", "DEG_summary_per_combination.csv")

# ------------------------ #
# Load data
# ------------------------ #
df <- read.csv(csv_path, stringsAsFactors = FALSE)

# ------------------------ #
# Aggregate DEG counts
# ------------------------ #
# Total DEGs per dims × resolution combination
summary_df <- df %>%
  group_by(dims, resolution) %>%
  summarise(total_DEGs = sum(n_DEGs), .groups = "drop") %>%
  arrange(dims, resolution)

# Save summary CSV
write.csv(summary_df, summary_csv, row.names = FALSE)

# ------------------------ #
# Create DEG matrix for heatmap
# ------------------------ #
heatmap_df <- summary_df %>%
  pivot_wider(names_from = dims, values_from = total_DEGs)

rownames(heatmap_df) <- heatmap_df$resolution
heatmap_mat <- as.matrix(heatmap_df[ , -which(names(heatmap_df) == "resolution")])

# ------------------------ #
# Plot Heatmap
# ------------------------ #
pheatmap(
  heatmap_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.0f",
  main = "Total # of DEGs per Cluster (All Samples)",
  color = colorRampPalette(c("white", "#3182bd"))(100),
  filename = heatmap_pdf,
  width = 7, height = 5
)

cat("✅ Heatmap saved to:", heatmap_pdf, "\n")
cat("✅ Summary table saved to:", summary_csv, "\n")