library(Seurat)
library(dplyr)
library(ggplot2)

# --------------------------- #
# Paths & Settings
# --------------------------- #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
rds_file <- file.path(output_dir, "rds", "merged_seurat_analysis_20250616.rds")
pdf_dir <- file.path(output_dir, "pdf", "3_merged_DEGs_outs")
csv_dir <- file.path(output_dir, "tsv", "merged", "3_merged_DEGs_outs")
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

pdf_out <- file.path(pdf_dir, "cluster_DEG_counts_per_cluster_per_sample_v2.pdf")
csv_out <- file.path(csv_dir, "cluster_DEG_counts_per_cluster_per_sample.csv")

dims_list <- c(10, 15, 20, 25, 30)
res_list <- c(0.25, 0.5, 0.75)

# --------------------------- #
# Load merged object
# --------------------------- #
merged_rds <- readRDS(rds_file)
merged_obj <- merged_rds[["SCT"]]

# --------------------------- #
# DEG collection
# --------------------------- #
deg_summary <- data.frame()

for (dims in dims_list) {
  for (res in res_list) {
    cat("🔬 Running DGE for dims =", dims, "res =", res, "\n")
    
    merged_obj <- FindNeighbors(merged_obj, dims = 1:dims, verbose = FALSE)
    merged_obj <- FindClusters(merged_obj, resolution = res, verbose = FALSE)
    
    cluster_col <- paste0("integrated_snn_res.d", dims, "_r", res)
    merged_obj[[cluster_col]] <- merged_obj$seurat_clusters
    merged_obj$cluster_tmp <- as.character(merged_obj[[cluster_col]][, 1])
    Idents(merged_obj) <- "cluster_tmp"
    
    merged_obj <- PrepSCTFindMarkers(merged_obj)
    
    try({
      markers <- FindAllMarkers(
        merged_obj,
        logfc.threshold = 0.25,
        min.pct = 0.1,
        only.pos = TRUE
      )
      
      if (!is.null(markers) && "cluster" %in% colnames(markers)) {
        # Count DEGs per cluster × sample
        markers <- markers %>% 
          filter(p_val_adj < 0.05, avg_log2FC > 0.25)
        
        cluster_animals <- merged_obj@meta.data %>%
          dplyr::select(cluster_tmp, animal) %>%
          dplyr::mutate(cluster_tmp = as.character(cluster_tmp))
        
        markers <- markers %>%
          left_join(cluster_animals, by = c("cluster" = "cluster_tmp")) %>%
          group_by(cluster, animal) %>%
          summarise(n_DEGs = n(), .groups = "drop") %>%
          mutate(dims = dims, resolution = res)
        
        deg_summary <- bind_rows(deg_summary, markers)
      }
    })
  }
}

# --------------------------- #
# Format metadata
# --------------------------- #
deg_summary <- deg_summary %>%
  mutate(
    cluster = paste0("C", cluster),
    cluster = factor(cluster, levels = paste0("C", sort(unique(as.numeric(gsub("C", "", cluster)))))),
    dims = factor(dims),
    resolution = factor(resolution)
  )

# --------------------------- #
# Save CSV
# --------------------------- #
write.csv(deg_summary, csv_out, row.names = FALSE)
cat("✅ Saved summary CSV to:", csv_out, "\n")

# --------------------------- #
# Plot
# --------------------------- #
p <- ggplot(deg_summary, aes(x = cluster, y = n_DEGs, color = animal, group = animal)) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
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
cat("✅ Saved scatter plot to:", pdf_out, "\n")