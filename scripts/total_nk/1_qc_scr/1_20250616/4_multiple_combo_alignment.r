# --------------------------------------------- #
# Automatic Cluster Alignment & Label Transfer
# Multiple dims × res combinations
# --------------------------------------------- #

library(Seurat)
library(dplyr)
library(pheatmap)
library(igraph)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# ------------------------- #
# Settings & Paths
# ------------------------- #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
pdf_base_dir <- file.path(output_dir, "pdf", "4_cluster_align_filtering")
dir.create(pdf_base_dir, recursive = TRUE, showWarnings = FALSE)

rds_file <- file.path(output_dir, "rds", "merged_seurat_analysis_20250616.rds")
merged_rds <- readRDS(rds_file)
merged_obj_full <- merged_rds[["SCT"]]

dims_list <- c(10, 15, 20, 25, 30)
res_list <- c(0.25, 0.5, 0.75)

# ------------------------- #
# Loop Over dims × res
# ------------------------- #
for (dims in dims_list) {
  for (res in res_list) {
    
    cat("\n==============================\n")
    cat("🚀 Running alignment for dims =", dims, "| res =", res, "\n")
    cat("==============================\n")
    
    # Set output folder
    combo_dir <- file.path(pdf_base_dir, paste0("dims", dims, "_res", res))
    dir.create(combo_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Reload object for clean state
    merged_obj <- merged_obj_full
    
    # Run clustering
    merged_obj <- FindNeighbors(merged_obj, dims = 1:dims, verbose = FALSE)
    merged_obj <- FindClusters(merged_obj, resolution = res, verbose = FALSE)
    cluster_col <- paste0("integrated_snn_res.d", dims, "_r", res)
    merged_obj[[cluster_col]] <- merged_obj$seurat_clusters
    
    # Get sample info
    animals <- unique(merged_obj$animal)
    cat("🐾 Found samples:", paste(animals, collapse = ", "), "\n")
    
    # Split & prep SCT
    split_objs <- SplitObject(merged_obj, split.by = "animal")
    features <- SelectIntegrationFeatures(split_objs, nfeatures = 3000)
    split_objs <- PrepSCTIntegration(split_objs, anchor.features = features)
    
    # ------------------------- #
    # Marker Discovery (with filtering)
    # ------------------------- #
    marker_list <- list()
    for (a in animals) {
      cat("🔬 Finding markers for", a, "\n")
      Idents(split_objs[[a]]) <- split_objs[[a]][[cluster_col]][, 1]
      marker_list[[a]] <- FindAllMarkers(
        split_objs[[a]],
        only.pos = TRUE,
        logfc.threshold = 0.25,
        min.pct = 0.1
      ) %>% 
        filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25)  # ← Added filtering here
    }
    
    # ------------------------- #
    # Jaccard Similarity Matrices
    # ------------------------- #
    compute_jaccard <- function(df1, df2) {
      cl1 <- split(df1$gene, df1$cluster)
      cl2 <- split(df2$gene, df2$cluster)
      mat <- outer(names(cl1), names(cl2), Vectorize(function(i, j) {
        length(intersect(cl1[[i]], cl2[[j]])) / length(union(cl1[[i]], cl2[[j]]))
      }))
      rownames(mat) <- paste0("C", names(cl1))
      colnames(mat) <- paste0("C", names(cl2))
      return(mat)
    }
    
    alignment_matrices <- list()
    for (i in 1:(length(animals) - 1)) {
      for (j in (i + 1):length(animals)) {
        a1 <- animals[i]; a2 <- animals[j]
        cat("📊 Aligning:", a1, "vs", a2, "\n")
        mat <- compute_jaccard(marker_list[[a1]], marker_list[[a2]])
        pair_name <- paste0(a1, "_vs_", a2)
        alignment_matrices[[pair_name]] <- mat
        pheatmap(mat,
                 cluster_rows = FALSE, cluster_cols = FALSE,
                 display_numbers = TRUE, number_format = "%.2f",
                 main = pair_name,
                 filename = file.path(combo_dir, paste0("alignment_", pair_name, ".pdf")))
      }
    }
    
    # ------------------------- #
    # Label Transfer Example
    # ------------------------- #
    cat("🔁 Label transfer:", animals[1], "→", animals[2], "\n")
    ref <- animals[1]; query <- animals[2]
    
    # Reprocess for anchors
    for (i in c(ref, query)) {
      split_objs[[i]] <- SCTransform(split_objs[[i]], verbose = FALSE)
      split_objs[[i]] <- FindVariableFeatures(split_objs[[i]])
      split_objs[[i]] <- RunPCA(split_objs[[i]], verbose = FALSE)
    }
    
    anchors <- FindTransferAnchors(
      reference = split_objs[[ref]],
      query = split_objs[[query]],
      normalization.method = "SCT",
      dims = 1:dims,
      recompute.residuals = FALSE  # important for SCT integration
    )
    
    transfer <- TransferData(
      anchorset = anchors,
      refdata = split_objs[[ref]]$seurat_clusters,
      dims = 1:dims
    )
    
    split_objs[[query]]$predicted_cluster <- transfer$predicted.id
    
    # Plot UMAPs
    p1 <- DimPlot(split_objs[[query]], group.by = cluster_col, label = TRUE) +
      ggtitle(paste("True Clusters:", query))
    p2 <- DimPlot(split_objs[[query]], group.by = "predicted_cluster", label = TRUE) +
      ggtitle(paste("Predicted from", ref))
    
    pdf(file.path(combo_dir, paste0("label_transfer_", ref, "_to_", query, ".pdf")),
        width = 10, height = 5)
    print(p1 + p2)
    dev.off()
    
    # ------------------------- #
    # Global Network Plot
    # ------------------------- #
    cat("🌐 Generating network plot...\n")
    edges <- do.call(rbind, lapply(names(alignment_matrices), function(pair) {
      mat <- alignment_matrices[[pair]]
      pair_split <- unlist(strsplit(pair, "_vs_"))
      a1 <- pair_split[1]
      a2 <- pair_split[2]
      edge_df <- which(mat > 0.1, arr.ind = TRUE)
      if (nrow(edge_df) == 0) return(NULL)
      data.frame(
        from = paste0(a1, "_", rownames(mat)[edge_df[, 1]]),
        to   = paste0(a2, "_", colnames(mat)[edge_df[, 2]]),
        weight = mat[edge_df],
        pair = pair
      )
    }))
    
    if (nrow(edges) > 0) {
      g <- graph_from_data_frame(edges, directed = FALSE)
      V(g)$sample <- sapply(strsplit(V(g)$name, "_"), `[`, 1)
      animal_colors <- brewer.pal(length(animals), "Set2")
      names(animal_colors) <- animals
      V(g)$color <- animal_colors[V(g)$sample]
      
      pdf(file.path(combo_dir, "cluster_alignment_network_fixed.pdf"), width = 10, height = 7)
      plot(g,
           layout = layout_with_fr(g),
           vertex.size = 8,
           vertex.label.cex = 0.8,
           vertex.label.color = "black",
           vertex.color = V(g)$color,
           edge.width = E(g)$weight * 5,
           edge.color = "gray60",
           main = paste("Cluster Alignment Network\nDims =", dims, ", Res =", res))
      dev.off()
    }
    
    cat("✅ Finished:", combo_dir, "\n")
  }
}