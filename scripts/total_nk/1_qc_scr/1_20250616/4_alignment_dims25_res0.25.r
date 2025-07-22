# ----------------------------- #
# Libraries
# ----------------------------- #
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)

# ----------------------------- #
# Paths & Parameters
# ----------------------------- #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
rds_file <- file.path(output_dir, "rds", "merged_seurat_analysis_20250616.rds")
marker_file <- file.path(output_dir, "tsv", "merged", "SCT", "new_markers_merged_SCT_dims15_res0.25.csv")

dims <- 15
res <- 0.25
res_label <- paste0("res", res)
pdf_dir <- file.path(output_dir, "pdf", "4_cluster_align_filtering", paste0("dims", dims, "_", res_label))
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

heatmap_pdf <- file.path(pdf_dir, "heatmap_top20_cluster_markers_filtering.pdf")
corr_pdf <- file.path(pdf_dir, "cluster_expression_correlation_filtering.pdf")
umap_pdf <- file.path(pdf_dir, "umap_cluster_identity_filtering.pdf")

# ----------------------------- #
# Load data
# ----------------------------- #
cat("🔄 Loading merged Seurat object...\n")
merged_rds <- readRDS(rds_file)
merged_obj <- merged_rds[["SCT"]]  # access the SCT object

# ----------------------------- #
# Recompute clustering at correct dims & res
# ----------------------------- #
cat("♻️ Recomputing clustering (dims =", dims, ", res =", res, ")...\n")
DefaultAssay(merged_obj) <- "integrated"

if (!"pca" %in% Reductions(merged_obj)) {
  merged_obj <- RunPCA(merged_obj, verbose = FALSE)
}
merged_obj <- FindNeighbors(merged_obj, dims = 1:dims)
merged_obj <- FindClusters(merged_obj, resolution = res)

# Store into metadata column
cluster_col <- paste0("integrated_snn_res.", res)
merged_obj[[cluster_col]] <- merged_obj$seurat_clusters

# Set identities
Idents(merged_obj) <- merged_obj[[cluster_col]]
merged_obj$seurat_clusters <- merged_obj[[cluster_col]]

# ----------------------------- #
# Load markers
# ----------------------------- #
markers <- read.csv(marker_file, stringsAsFactors = FALSE)

# ----------------------------- #
# Filter top 20 markers per cluster
# ----------------------------- #
if (all(c("cluster", "gene") %in% colnames(markers))) {
    cat("📊 Filtering top 20 markers per cluster...\n")
    top_markers <- markers %>%
        filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>% # Filter markers based on adjusted p-value and log2FC
        group_by(cluster) %>%
        slice_max(n = 20, order_by = avg_log2FC)

    top_genes <- unique(top_markers$gene)
    top_genes <- intersect(top_genes, rownames(merged_obj))
    if (length(top_genes) == 0) stop("❌ No matching genes found in Seurat object.")

    # ----------------------------- #
    # 1️⃣ Heatmap of Top Markers
    # ----------------------------- #
    cat("🧬 Creating heatmap of top markers...\n")
    heatmap_plot <- DoHeatmap(merged_obj, features = top_genes, group.by = "seurat_clusters") +
        ggtitle(paste("Top 20 Markers (dims =", dims, ", res =", res, ")"))

    ggsave(
        filename = heatmap_pdf,
        plot = heatmap_plot,
        width = 12,
        height = 16,
        dpi = 600
    )
    cat("✅ Heatmap saved to:", heatmap_pdf, "\n")

  # ----------------------------- #
    # 2️⃣ Correlation of Average Expression
    # ----------------------------- #
    cat("🔗 Calculating average expression per cluster...\n")
    avg_expr_all <- AverageExpression(
    merged_obj,
    group.by = "seurat_clusters",
    assays = "SCT",
    slot = "data"
    )

    if (!"SCT" %in% names(avg_expr_all)) {
    stop("❌ Assay 'SCT' not found in AverageExpression result.")
    }

    avg_expr <- avg_expr_all[["SCT"]]
    if (!is.matrix(avg_expr)) {
    avg_expr <- as.matrix(avg_expr)
    }

    # Filter expression matrix to only include top marker genes
    top_genes_filtered <- intersect(rownames(avg_expr), top_genes)
    if (length(top_genes_filtered) < 5) {
      stop("❌ Too few marker genes found in average expression matrix.")
    }
    cat("🔍 Filtering average expression matrix to top marker genes (p_val_adj < 0.05, avg_log2FC > 0.25)...\n")

    avg_expr <- avg_expr[top_genes_filtered, , drop = FALSE]

    # Now compute correlation
    cor_mat <- cor(avg_expr)

    # Plot
    pdf(corr_pdf, width = 7, height = 6)
    pheatmap(
    cor_mat,
    main = paste("Cluster Expression Correlation\n(dims =", dims, ", res =", res, ")"),
    display_numbers = TRUE,
    number_format = "%.2f",
    cluster_rows = TRUE,
    cluster_cols = TRUE
    )
    dev.off()
    cat("✅ Correlation heatmap saved to:", corr_pdf, "\n")

    # ----------------------------- #
    # 3️⃣ UMAP Plot of Clusters
    # ----------------------------- #
    cat("🎨 Generating UMAP...\n")
    p_umap <- DimPlot(merged_obj, group.by = "seurat_clusters", label = TRUE, label.size = 5) +
        ggtitle(paste("UMAP (dims =", dims, ", res =", res, ")")) +
        theme_minimal()

    ggsave(
        filename = umap_pdf,
        plot = p_umap,
        width = 7,
        height = 6,
        dpi = 600
    )
    cat("✅ UMAP saved to:", umap_pdf, "\n")

} else {
  cat("⚠️ Required columns ('cluster', 'gene') not found in marker file.\n")
}