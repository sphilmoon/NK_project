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
pdf_out <- file.path(pdf_dir, "cluster_DEG_counts_per_cluster_per_sample.pdf")

dims_list <- c(10, 15, 20, 25, 30)
res_list <- c(0.25, 0.5, 0.75)

# --------------------------- #
# Load Merged Seurat Object
# --------------------------- #
merged_rds <- readRDS(rds_file)
merged_obj <- merged_rds[["SCT"]]
animals <- unique(merged_obj$animal)

# --------------------------- #
# Initialize storage
# --------------------------- #
all_deg_counts <- data.frame()

# --------------------------- #
# Loop over dims × resolution × animal
# --------------------------- #
for (dims in dims_list) {
  for (res in res_list) {
    cat("📊🔬 Processing DGE analysis at dims =", dims, "| res =", res, "\n")

    # Cluster entire object
    merged_obj <- FindNeighbors(merged_obj, dims = 1:dims, verbose = FALSE)
    merged_obj <- FindClusters(merged_obj, resolution = res, verbose = FALSE)

    cluster_col <- paste0("integrated_snn_res.d", dims, "_r", res)
    merged_obj[[cluster_col]] <- merged_obj$seurat_clusters

    # Assign identities globally
    merged_obj$cluster_tmp <- as.character(merged_obj[[cluster_col]][, 1])

    for (animal in animals) {
      obj_subset <- subset(merged_obj, subset = animal == !!animal)
      if (ncol(obj_subset) < 100) next  # Skip small samples

      obj_subset$cluster_tmp <- as.character(obj_subset[[cluster_col]][, 1])
      Idents(obj_subset) <- "cluster_tmp"

      # Prep SCT
      obj_subset <- PrepSCTFindMarkers(obj_subset)

      markers <- tryCatch({
        FindAllMarkers(obj_subset, logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE)
      }, error = function(e) {
        message("❌ DGE analysis failed for animal:", animal)
        return(NULL)
      })

      if (!is.null(markers)) {
        deg_counts <- markers %>%
          filter(p_val_adj < 0.05, avg_log2FC > 0.25) %>%
          group_by(cluster) %>%
          summarise(n_DEGs = n(), .groups = "drop") %>%
          mutate(
            dims = dims,
            resolution = res,
            animal = animal
          ) %>%
          mutate(
            cluster = paste0("C", cluster),
            cluster = factor(cluster, levels = paste0("C", sort(unique(as.numeric(cluster)))))
          )

        all_deg_counts <- bind_rows(all_deg_counts, deg_counts)
      }
    }
  }
}

# --------------------------- #
# Plot
# --------------------------- #
p <- ggplot(all_deg_counts, aes(x = cluster, y = n_DEGs, color = animal, group = animal)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.2) +
  geom_line(position = position_dodge(width = 0.6), linewidth = 0.6, alpha = 0.7) +
  facet_grid(resolution ~ dims, scales = "free_x", space = "free_x") +
  labs(
    title = "# DEGs per Cluster per Sample (adj p < 0.05, logFC > 0.25)",
    x = "Cluster", y = "# DEGs", color = "Animal"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

ggsave(pdf_out, plot = p, width = 16, height = 6)
cat("✅ Saved: ", pdf_out, "\n")
