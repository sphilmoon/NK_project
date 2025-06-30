library(Seurat)
library(dplyr)
library(ggplot2)

# --------------------------- #
# Load Seurat merged object
# --------------------------- #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
merged_rds <- readRDS(file.path(output_dir, "rds", "merged_seurat_analysis_20250616.rds"))

# Use SCT slot (or "LogNormalize" if you're checking that instead)
merged_obj <- merged_rds[["SCT"]]  # or merged_rds[["LogNormalize"]]

# --------------------------- #
# Set cluster identity
# --------------------------- #
# Set to a specific clustering (e.g., resolution 0.5)
res <- 0.5
cluster_col <- paste0("integrated_snn_res.", res)

if (!(cluster_col %in% colnames(merged_obj@meta.data))) stop("Cluster column not found.")

merged_obj$cluster <- merged_obj[[cluster_col]][, 1]
Idents(merged_obj) <- "cluster"

# --------------------------- #
# Calculate cell fraction per cluster per animal
# --------------------------- #
df_frac <- table(Cluster = merged_obj$cluster, Sample = merged_obj$animal) %>%
  as.data.frame() %>%
  group_by(Cluster) %>%
  mutate(Fraction = Freq / sum(Freq))

# --------------------------- #
# Plot: Stacked Bar Chart
# --------------------------- #
p <- ggplot(df_frac, aes(x = as.factor(Cluster), y = Fraction, fill = Sample)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = paste("Cluster Composition by Sample (Resolution", res, ")"),
       x = "Cluster", y = "Fraction of Cells", fill = "Animal") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save to PDF
pdf_file <- file.path(output_dir, "pdf", paste0("merged_cluster_fraction_by_animal_res", res, ".pdf"))
ggsave(pdf_file, p, width = 10, height = 5)
cat("✅ Cluster composition plot saved to:", pdf_file, "\n")
