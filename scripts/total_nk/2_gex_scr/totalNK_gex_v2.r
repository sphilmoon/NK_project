# Load required libraries
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(patchwork)
library(cowplot)
library(RColorBrewer)
library(pheatmap)
install.packages("logr")
library(logr)

# Set up logging
log_file <- "dge_analysis_log.log"
if (file.exists(log_file)) file.remove(log_file)
log_open(log_file)
log_print("🚀 Starting differential gene expression analysis...", console = TRUE)

# Define configurable parameters
base_dir <- "/home/rawdata/total_nk"
sample_info <- list(
  Animal25 = "animal25_totalnk_filtered_feature_bc_matrix.h5",
  Animal26 = "animal26_totalnk_filtered_feature_bc_matrix.h5",
  Animal27 = "animal27_totalnk_filtered_feature_bc_matrix.h5",
  Animal28 = "animal28_totalnk_filtered_feature_bc_matrix.h5", 
  Animal51 = "animal51_totalnk_filtered_feature_bc_matrix.h5",
  Animal52 = "animal52_totalnk_filtered_feature_bc_matrix.h5"
)
min_cells <- 3
min_features <- 200
dims_to_use <- 1:30
n_features_integration <- 3000
cluster_resolution <- 0.5

# Define output directories
qc_output_dir <- "/home/outputs/totalNK_outputs/1_qc"
dge_output_dir <- "/home/outputs/totalNK_outputs/2_dge"

# Function to load H5 files into Seurat objects with error handling
load_seurat_objects <- function(base_dir, sample_info) {
  seurat_objects <- list()
  for (sample in names(sample_info)) {
    file_path <- file.path(base_dir, sample_info[[sample]])
    if (!file.exists(file_path)) {
      log_print(sprintf("❌ Error: File %s not found.", file_path), console = TRUE)
      stop("Missing input file.")
    }
    log_print(sprintf("📥 Loading sample: %s", sample), console = TRUE)
    data <- tryCatch(
      Read10X_h5(file_path),
      error = function(e) {
        log_print(sprintf("❌ Error loading %s: %s", file_path, e$message), console = TRUE)
        NULL
      }
    )
    if (is.null(data)) next
    seurat_obj <- CreateSeuratObject(
      counts = data,
      project = sample,
      min.cells = min_cells,
      min.features = min_features
    )
    seurat_obj$sample <- sample
    seurat_objects[[sample]] <- seurat_obj
    log_print(sprintf("✅ Sample %s: %d cells, %d genes", sample, ncol(seurat_obj), nrow(seurat_obj)), console = TRUE)
  }
  if (length(seurat_objects) == 0) stop("No valid Seurat objects created.")
  return(seurat_objects)
}

# Load data
seurat_objects <- load_seurat_objects(base_dir, sample_info)
saveRDS(seurat_objects, file.path(qc_output_dir, "seurat_objects_list_dim30.rds"))
log_print("📂 Seurat objects saved to 1_qc/seurat_objects_list_dim30.rds 🎉", console = TRUE)

# Function for QC filtering with sample-specific thresholds
perform_qc <- function(seurat_obj, sample_name) {
  qc_params <- list(
    Animal27 = list(nFeature_min = 200, nFeature_max = 3000, nCount_min = 250, nCount_max = 7500),
    Animal28 = list(nFeature_min = 200, nFeature_max = 3000, nCount_min = 250, nCount_max = 10000),
    default = list(nFeature_min = 200, nFeature_max = 3000, nCount_min = 500, nCount_max = 10000)
  )
  params <- if (sample_name %in% names(qc_params)) qc_params[[sample_name]] else qc_params$default
  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA > params$nFeature_min & nFeature_RNA < params$nFeature_max &
             nCount_RNA > params$nCount_min & nCount_RNA < params$nCount_max
  )
  return(seurat_obj)
}

# QC and visualization
combined_plots <- list()
for (sample in names(seurat_objects)) {
  log_print(sprintf("🛠️ Performing QC for %s", sample), console = TRUE)
  seurat_objects[[sample]] <- perform_qc(seurat_objects[[sample]], sample)
  
  p1 <- VlnPlot(seurat_objects[[sample]], features = "nFeature_RNA") +
        ggtitle("Gene counts QC") +
        theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  p2 <- VlnPlot(seurat_objects[[sample]], features = "nCount_RNA") +
        ggtitle("UMI counts QC") +
        theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  combined_plots[[sample]] <- p1 + p2 + plot_annotation(title = sprintf("QC for %s", sample))
}

stacked_plot <- plot_grid(plotlist = combined_plots, ncol = 1)
ggsave(file.path(qc_output_dir, "qc_violin_stacked_dim30.png"), stacked_plot, width = 10, height = 6 * length(seurat_objects), dpi = 600)
log_print("📊 QC plots saved to 1_qc/qc_violin_stacked_dim30.png 🎉", console = TRUE)

# Normalize with SCTransform
seurat_objects <- lapply(seurat_objects, function(x) {
  SCTransform(x)
})
saveRDS(seurat_objects, file.path(qc_output_dir, "seurat_objects_qc_sct_dim30.rds"))
log_print("📂 Normalized Seurat objects saved to 1_qc/seurat_objects_qc_sct_dim30.rds 🎉", console = TRUE)

# Batch correction and integration
# for (sample in names(seurat_objects)) {
#   seurat_objects[[sample]] <- RunPCA(seurat_objects[[sample]])
#   plot <- ElbowPlot(seurat_objects[[sample]], ndims = 50) +
#           ggtitle(sprintf("Elbow Plot for PCA - %s", sample)) +
#           theme(plot.title = element_text(hjust = 0.5))
#   combined_plots[[sample]] <- plot
# }
# stacked_elbow_plot <- plot_grid(plotlist = combined_plots, ncol = 1)
# ggsave(file.path(qc_output_dir, "combined_elbow_plots_dim30.png"), stacked_elbow_plot, width = 10, height = 6 * length(seurat_objects), dpi = 600)
# log_print("📊 Elbow plots saved to 1_qc/combined_elbow_plots_dim30.png 🎉", console = TRUE)

features <- SelectIntegrationFeatures(object.list = seurat_objects, nfeatures = n_features_integration)
anchors <- FindIntegrationAnchors(object.list = seurat_objects, anchor.features = features, dims = dims_to_use)
integrated_data <- IntegrateData(anchorset = anchors, dims = dims_to_use)

# Clustering and visualization
DefaultAssay(integrated_data) <- "integrated"
integrated_data <- ScaleData(integrated_data) %>%
                  RunPCA() %>%
                  RunUMAP(dims = dims_to_use) %>%
                  FindNeighbors(dims = dims_to_use) %>%
                  FindClusters(resolution = cluster_resolution)
saveRDS(integrated_data, file.path(dge_output_dir, "integrated_data_clusters_dim30.rds"))
log_print("📂 Integrated data with clusters saved to 2_dge/integrated_data_clusters_dim30.rds 🎉", console = TRUE)

dim_plot <- DimPlot(integrated_data, group.by = "sample", pt.size = 0.5) +
            ggtitle("UMAP - Integration by All 6 animals") +
            theme(plot.title = element_text(hjust = 0.5)) +
            scale_color_brewer(palette = "Set2")
ggsave(file.path(dge_output_dir, "dimplot_integration_dim30.png"), dim_plot, width = 10, height = 8, dpi = 600)
log_print("📊 UMAP plot saved to 2_dge/dimplot_integration_dim30.png 🎉", console = TRUE)

cluster_plot <- DimPlot(integrated_data, label = TRUE, group.by = "seurat_clusters")
ggsave(file.path(dge_output_dir, "cluster_plot_dim30.png"), cluster_plot, width = 10, height = 8, dpi = 600)
log_print("📊 Cluster plot saved to 2_dge/cluster_plot_dim30.png 🎉", console = TRUE)

# # Differential Expression Analysis by Cluster
# DefaultAssay(integrated_data) <- "RNA"
# integrated_data <- JoinLayers(integrated_data, assay = "RNA") %>%
#                   NormalizeData(assay = "RNA") %>%
#                   ScaleData(assay = "RNA")
# degs_all_clusters <- FindAllMarkers(
#   integrated_data,
#   assay = "RNA",
#   min.pct = 0.25,
#   logfc.threshold = 0.25,
#   test.use = "wilcox"
# )
# degs_filtered <- degs_all_clusters %>%
#                 filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)

# shared_degs_list <- lapply(unique(degs_filtered$cluster), function(cl) {
#   degs_cluster <- degs_filtered %>% filter(cluster == cl)
#   genes <- unique(degs_cluster$gene)
#   cluster_data <- subset(integrated_data, idents = cl)
#   avg_exp <- AverageExpression(cluster_data, assays = "RNA", features = genes, group.by = "sample")$RNA
#   shared_genes <- sapply(genes, function(gene) {
#     gene_logfc <- degs_cluster %>% filter(gene == !!gene) %>% pull(avg_log2FC)
#     direction <- ifelse(gene_logfc > 0, "up", "down")
#     exp_values <- avg_exp[gene, ]
#     consistent <- if (direction == "up") all(exp_values > 0) else all(exp_values < 0)
#     return(if (consistent) gene else NULL)
#   })
#   degs_cluster %>% filter(gene %in% shared_genes)
# })
# shared_degs <- do.call(rbind, shared_degs_list)
# write.csv(shared_degs, file.path(dge_output_dir, "shared_degs_across_samples_dim30.csv"), row.names = FALSE)
# saveRDS(shared_degs, file.path(dge_output_dir, "shared_degs_gex_analysis_dim30.rds"))
# log_print("🔬 Shared DEGs analysis completed and saved to 2_dge 🎉", console = TRUE)

# # Visualization of shared DEGs
# top_degs <- shared_degs %>%
#             group_by(cluster) %>%
#             top_n(5, wt = abs(avg_log2FC)) %>%
#             ungroup() %>%
#             pull(gene) %>%
#             unique()
# dot_plot <- DotPlot(integrated_data, features = top_degs, group.by = "seurat_clusters", assay = "RNA") +
#             ggtitle("Dot Plot of Shared DEGs") +
#             theme(plot.title = element_text(hjust = 0.5))
# ggsave(file.path(dge_output_dir, "dotplot_shared_degs_dim30.png"), dot_plot, width = 14, height = 8, dpi = 600)
# log_print("📊 Dot plot saved to 2_dge/dotplot_shared_degs_dim30.png 🎉", console = TRUE)




# re-run from here in R console.

install.packages("logr")
library(logr)

integrated_data <- readRDS("/home/outputs/totalNK_outputs/2_dge/integrated_data_clusters_dim30.rds")

# Differential Expression Analysis by Sample
Idents(integrated_data) <- "sample"
sample_degs <- lapply(unique(integrated_data$sample), function(s) {
  log_print(sprintf("🔍 Computing DEGs for sample: %s", s), console = TRUE)
  degs <- FindMarkers(
    integrated_data,
    ident.1 = s,
    assay = "RNA",
    min.pct = 0.25,
    logfc.threshold = 0.25,
    test.use = "wilcox"
  )
  degs$gene <- rownames(degs)
  degs$sample <- s
  degs
})
all_sample_degs <- do.call(rbind, sample_degs)
write.csv(all_sample_degs, file.path(dge_output_dir, "all_sample_degs.csv"), row.names = FALSE)
log_print("📂 Saved all sample DEGs to 2_dge/all_sample_degs.csv 🎉", console = TRUE)

global_degs_df <- all_sample_degs %>%
                  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5) %>%
                  group_by(gene) %>%
                  summarise(mean_log2FC = mean(abs(avg_log2FC)), min_p_val_adj = min(p_val_adj))
write.csv(global_degs_df, file.path(dge_output_dir, "global_degs_summary.csv"), row.names = FALSE)
log_print("📂 Saved global DEGs summary to 2_dge/global_degs_summary.csv 🎉", console = TRUE)

global_degs_top20 <- global_degs_df %>%
                     arrange(min_p_val_adj, desc(mean_log2FC)) %>%
                     slice_head(n = 20) %>%
                     pull(gene)
avg_exp_top20 <- AverageExpression(integrated_data, assays = "RNA", features = global_degs_top20, group.by = "sample")$RNA
write.csv(avg_exp_top20, file.path(dge_output_dir, "avg_exp_top20_degs.csv"), row.names = TRUE)
log_print("📂 Saved average expression of top 20 DEGs to 2_dge/avg_exp_top20_degs.csv 🎉", console = TRUE)

pheatmap(avg_exp_top20, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("#00a2fa", "#53008e"))(50),
         filename = file.path(dge_output_dir, "heatmap_global_degs_top20_dim30.png"), width = 14, height = 8, dpi = 600)
log_print("📊 Heatmap saved to 2_dge/heatmap_global_degs_top20_dim30.png 🎉", console = TRUE)

log_print("🏁 Analysis complete! 🚀", console = TRUE)
log_close()