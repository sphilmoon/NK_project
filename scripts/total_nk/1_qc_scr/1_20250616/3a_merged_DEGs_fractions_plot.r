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
pdf_out <- file.path(pdf_dir, "cluster_DEG_counts_scatterplot_merged.pdf")

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
all_deg_counts <- data.frame()

# --------------------------- #
# Loop over dims × resolution
# --------------------------- #
for (dims in dims_list) {
  for (res in res_list) {
    cat("🔬 Running DGE: dims =", dims, "res =", res, "\n")

    # Recompute clustering
    merged_obj <- FindNeighbors(merged_obj, dims = 1:dims, verbose = FALSE)
    merged_obj <- FindClusters(merged_obj, resolution = res, verbose = FALSE)

    # Store unique cluster identities
    cluster_col <- paste0("integrated_snn_res.d", dims, "_r", res)
    merged_obj[[cluster_col]] <- merged_obj$seurat_clusters
    merged_obj$cluster_tmp <- as.character(merged_obj[[cluster_col]][, 1])
    Idents(merged_obj) <- "cluster_tmp"

    # Prep SCT assay
    merged_obj <- PrepSCTFindMarkers(merged_obj)

    # Run FindAllMarkers
    markers <- tryCatch({
      FindAllMarkers(merged_obj, logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE)
    }, error = function(e) {
      message("⚠️ Failed DGE at dims =", dims, "res =", res)
      return(NULL)
    })

    if (!is.null(markers) && "cluster" %in% colnames(markers)) {
      deg_counts <- markers %>%
        filter(p_val_adj < 0.05, avg_log2FC > 0.25) %>%
        group_by(cluster) %>%
        summarise(n_DEGs = n(), .groups = "drop") %>%
        mutate(
          dims = dims,
          resolution = res
        )

      # Format cluster label
      deg_counts <- deg_counts %>%
        mutate(
          cluster = paste0("C", cluster),
          cluster = factor(cluster, levels = paste0("C", sort(unique(as.numeric(cluster)))))
        )

      all_deg_counts <- bind_rows(all_deg_counts, deg_counts)
    }
  }
}

# --------------------------- #
# Plot scatterplot of DEG counts
# --------------------------- #
p <- ggplot(all_deg_counts, aes(x = cluster, y = n_DEGs, group = 1)) +
  geom_point(color = "#2C7BB6", size = 2) +
  geom_line(color = "#2C7BB6", linewidth = 0.7, alpha = 0.7) +
  facet_grid(resolution ~ dims, scales = "free_x", space = "free_x") +
  labs(
    title = "# DEGs per Cluster (adj p < 0.05, logFC > 0.25)",
    x = "Cluster", y = "# DEGs"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

# --------------------------- #
# Save to PDF
# --------------------------- #
ggsave(pdf_out, plot = p, width = 16, height = 6)
cat("✅ Saved DEG scatter plot to:", pdf_out, "\n")