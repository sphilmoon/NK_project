# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(stringr)
library(cluster) # silhouette()
library(dplyr)

# ------------------------- #
# Define Directory Structure
# ------------------------- #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
rds_output_dir <- file.path(output_dir, "rds")
pdf_output_dir <- file.path(output_dir, "pdf")
tsv_output_dir <- file.path(output_dir, "tsv")
dir.create(rds_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tsv_output_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------- #
# Define Input Directories
# ------------------------- #
totalNK_h5_dir <- "/home/rawdata/totalNK"
animal_names <- c("animal25", "animal26", "animal27", "animal28")


# ------------------------- #
# Build File Paths
# ------------------------- #
build_h5_paths <- function(base_dir, condition_prefix) {
 fnames <- paste0(tolower(animal_names), "_", condition_prefix, "_filtered_feature_bc_matrix.h5")
 full_paths <- file.path(base_dir, fnames)
 missing <- !file.exists(full_paths)
 if (any(missing)) {
   cat("❌ Missing files:\n")
   print(full_paths[missing])
   stop("Some .h5 files are missing. Please check the paths or filenames.")
 }
 names(full_paths) <- animal_names
 return(full_paths)
}
totalNK_file_paths <- build_h5_paths(totalNK_h5_dir, "totalnk")


# ------------------------- #
# Load and QC Each Animal
# ------------------------- #
seurat_objects <- list()
combined_plots <- list()
perform_qc <- function(seurat_obj, sample_name) {
 qc_params <- list(
   animal27 = list(nFeature_min = 200, nFeature_max = 3000, nCount_min = 250, nCount_max = 7500),
   animal28 = list(nFeature_min = 200, nFeature_max = 3000, nCount_min = 250, nCount_max = 10000),
   default  = list(nFeature_min = 200, nFeature_max = 3000, nCount_min = 500, nCount_max = 10000)
 )
 params <- if (sample_name %in% names(qc_params)) qc_params[[sample_name]] else qc_params$default
 seurat_obj <- subset(
   seurat_obj,
   subset = nFeature_RNA > params$nFeature_min & nFeature_RNA < params$nFeature_max &
            nCount_RNA > params$nCount_min & nCount_RNA < params$nCount_max
 )
 return(seurat_obj)
}

for (i in seq_along(totalNK_file_paths)) {
 sample <- names(totalNK_file_paths)[i]
 filepath <- totalNK_file_paths[[i]]
 cat("📥 Reading:", sample, "\n")
 seurat_obj <- Read10X_h5(filepath) %>%
   CreateSeuratObject(project = sample, min.cells = 3, min.features = 200)
 seurat_obj$animal <- sample
 # Add mitochondrial percentage
 seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
 # QC filter
 seurat_obj <- perform_qc(seurat_obj, sample)
 # Normalize
 seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
 # Save object
 seurat_objects[[sample]] <- seurat_obj
 # Generate QC plots
 p1 <- VlnPlot(seurat_obj, features = "nFeature_RNA", pt.size = 0) + ggtitle("Gene counts")
 p2 <- VlnPlot(seurat_obj, features = "nCount_RNA", pt.size = 0) + ggtitle("UMI counts")
 combined_plot <- p1 + p2 + plot_layout(ncol = 2) + plot_annotation(title = paste("QC for", sample))
 combined_plots[[sample]] <- combined_plot
 cat("✅ QC and normalization complete for", sample, "\n")
}

# # Save all QC plots
# pdf(file.path(pdf_output_dir, "qc_totalNK_animals_20250616.pdf"), width = 10, height = 6 * length(combined_plots))
# for (plot in combined_plots) print(plot)
# dev.off()
# cat("✅ QC plots saved.\n")

# Save SCT-normalized baseline
seurat_objects_sct <- seurat_objects
# Save individual objects
saveRDS(seurat_objects_sct, file.path(rds_output_dir, "1_qc_sctransform_20250616.rds"))
cat("✅ QC'ed Seurat objects saved.\n")

# # ------------------------- #
# # JackStraw Analysis: Print Significant PCs Before Plotting
# # ------------------------- #
# jackstraw_plots <- list()

# for (animal in names(seurat_objects)) {
#   cat("📊 Running JackStraw (RNA assay) for", animal, "\n")

#   obj <- seurat_objects[[animal]]

#   # Use raw RNA assay
#   DefaultAssay(obj) <- "RNA"
#   obj <- NormalizeData(obj)
#   obj <- FindVariableFeatures(obj)
#   obj <- ScaleData(obj)
#   obj <- RunPCA(obj)

#   # Determine valid dimensions
#   available_pcs <- ncol(obj@reductions$pca@cell.embeddings)
#   use_dims <- 1:min(30, available_pcs)

#   # Run JackStraw
#   obj <- JackStraw(obj, num.replicate = 100)
#   obj <- ScoreJackStraw(obj, dims = use_dims)

#   # Print significant PCs
#   jack_df <- obj@reductions$pca@jackstraw$overall.p.values
#   colnames(jack_df) <- c("PC", "p_value")
#   sig_pcs <- jack_df %>% dplyr::filter(p_value < 0.05)

#   cat("📌 Significant PCs for", animal, "(p < 0.05):\n")
#   print(sig_pcs)

#   # Store plot
#   p <- JackStrawPlot(obj, dims = use_dims) +
#     ggtitle(paste("JackStraw Plot -", animal)) +
#     theme(plot.title = element_text(hjust = 0.5))
#   jackstraw_plots[[animal]] <- p

#   # Save back
#   seurat_objects[[animal]] <- obj
# }

# # Save plots to one file
# pdf(file.path(pdf_output_dir, "jackstraw_totalNK_animals_20250616.pdf"), width = 10, height = 6 * length(jackstraw_plots))
# for (plot in jackstraw_plots) print(plot)
# dev.off()
# cat("✅ JackStraw plots saved.\n")


# # ------------------------- #
# # Elbow Plot for Each Animal using SCT (SCT assay)
# # ------------------------- #
# cat("📊 Generating elbow plots...\n")
# elbow_plots <- list()

# for (animal in names(seurat_objects)) {
#   seurat_obj <- seurat_objects[[animal]]
#   DefaultAssay(seurat_obj) <- "SCT"
#   seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
#   elbow_plots[[animal]] <- ElbowPlot(seurat_obj, ndims = 30) + 
#     ggtitle(paste("Elbow Plot -", animal)) +
#     theme(plot.title = element_text(hjust = 0.5, size = 12))
# }

# # Combine and save elbow plots
# pdf(file.path(pdf_output_dir, "elbowplots_SCT_totalNK_each_animals_20250616.pdf"), width = 10, height = 8, onefile = TRUE)
# print(wrap_plots(elbow_plots, ncol = 2))
# dev.off()
# cat("✅ Elbow plots saved to PDF.\n")


# # ------------------------- #
# # Elbow Plot for Merged Seurat Object
# # ------------------------- #

# # Ensure PCA has been run (optional if already present)
# if (!"pca" %in% names(merged_obj@reductions)) {
#   DefaultAssay(merged_obj) <- "integrated"
#   merged_obj <- ScaleData(merged_obj, verbose = FALSE)
#   merged_obj <- RunPCA(merged_obj, npcs = 30, verbose = FALSE)
# }

# # Create Elbow Plot
# elbow_plot <- ElbowPlot(merged_obj, ndims = 30) +
#   ggtitle("Elbow Plot - Merged Total NK") +
#   theme(plot.title = element_text(hjust = 0.5))

# # Save to PDF
# ggsave(
#   filename = file.path(pdf_output_dir, "elbow_plot_merged_totalNK_20250616.pdf"),
#   plot = elbow_plot,
#   width = 8,
#   height = 6,
#   dpi = 600
# )
# cat("✅ Elbow plot saved to:", file.path(pdf_output_dir, "elbow_plot_merged_totalNK_20250616.pdf"), "\n")



# # ------------------------- #
# # Elbow Plot for Each Animal using LogNormalized Data (RNA assay)
# # ------------------------- #

# cat("📊 Generating LogNormalize-based Elbow Plots...\n")

# # Create list to store elbow plots
# lognorm_elbow_plots <- list()

# # Loop through each animal/sample
# for (animal in names(seurat_objects)) {
#   cat("🔄 Processing:", animal, "\n")

#   obj <- seurat_objects[[animal]]

#   # Set to RNA assay and log-normalize
#   DefaultAssay(obj) <- "RNA"
#   obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
#   obj <- FindVariableFeatures(obj)
#   obj <- ScaleData(obj)
#   obj <- RunPCA(obj)

#   # Generate elbow plot
#   elbow_plot <- ElbowPlot(obj, ndims = 30) +
#     ggtitle(paste("Elbow Plot (LogNormalized) -", animal)) +
#     theme(plot.title = element_text(hjust = 0.5, size = 12))

#   lognorm_elbow_plots[[animal]] <- elbow_plot

#   # Optionally store PCA-modified object
#   seurat_objects[[animal]] <- obj
# }

# # Save all elbow plots in one combined PDF
# elbow_pdf_path <- file.path(pdf_output_dir, "elbow_log_totalNK_animals_20250616.pdf")
# pdf(elbow_pdf_path, width = 10, height = 6 * ceiling(length(lognorm_elbow_plots) / 2))
# print(wrap_plots(lognorm_elbow_plots, ncol = 2))
# dev.off()

# cat("✅ LogNormalize-based elbow plots saved to:", elbow_pdf_path, "\n")


# # ------------------------- #
# # Elbow Plot for Merged Object (LogNormalize)
# # ------------------------- #

# cat("📦 Loading merged object for LogNormalize-based elbow plot...\n")

# # Set to RNA assay and normalize
# DefaultAssay(merged_obj) <- "RNA"
# merged_obj <- NormalizeData(merged_obj, normalization.method = "LogNormalize", scale.factor = 10000)
# merged_obj <- FindVariableFeatures(merged_obj)
# merged_obj <- ScaleData(merged_obj)
# merged_obj <- RunPCA(merged_obj)

# # Generate elbow plot
# elbow_plot_merged <- ElbowPlot(merged_obj, ndims = 30) +
#   ggtitle("Elbow Plot (LogNormalized) - Merged") +
#   theme(plot.title = element_text(hjust = 0.5, size = 12))

# # Save as PDF
# pdf_merged_path <- file.path(pdf_output_dir, "elbow_lognormalized_merged_totalNK_20250616.pdf")
# ggsave(
#   filename = pdf_merged_path,
#   plot = elbow_plot_merged,
#   width = 10,
#   height = 6,
#   dpi = 600
# )

# cat("✅ Merged LogNormalize-based elbow plot saved to:", pdf_merged_path, "\n")


# ------------------------- #
# Individual Animal Analysis: Testing Various Dimensions and Resolutions for both SCT and LogNormalize methods
# ------------------------- #

dims_list <- c(10, 15, 20, 25, 30)
resolutions <- c(0.25, 0.5, 0.75)
small_cluster_threshold <- 100

cluster_summary <- data.frame()

for (method in c("SCT", "LogNormalize")) {
  for (dim in dims_list) {
    dims_to_use <- 1:dim

    for (res in resolutions) {
      plot_list <- list()

      for (animal in names(seurat_objects)) {
        cat("📌 Processing:", animal, "method:", method, "dims:", dim, "res:", res, "\n")
        obj <- seurat_objects[[animal]]

        # Assay setup
        if (method == "SCT") {
          DefaultAssay(obj) <- "SCT"
          obj <- ScaleData(obj, verbose = FALSE)
        } else {
          DefaultAssay(obj) <- "RNA"
          obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
          obj <- FindVariableFeatures(obj)
          obj <- ScaleData(obj)
        }

        # Run PCA
        obj <- RunPCA(obj, npcs = max(dims_list), verbose = FALSE)

        # Calculate cumulative variance explained
        var_explained <- obj@reductions$pca@stdev^2 / sum(obj@reductions$pca@stdev^2)
        cumulative_variance <- sum(var_explained[dims_to_use])

        # Continue analysis
        obj <- RunUMAP(obj, dims = dims_to_use, reduction.name = paste0("umap_", method, "_dims", dim))
        obj <- FindNeighbors(obj, dims = dims_to_use)
        obj <- FindClusters(obj, resolution = res)

        seurat_objects[[animal]] <- obj  # update object

        # UMAP plot
        p <- DimPlot(
          obj,
          reduction = paste0("umap_", method, "_dims", dim),
          group.by = "seurat_clusters",
          label = TRUE,
          pt.size = 0.3,
          repel = TRUE
        ) +
          ggtitle(paste0("UMAP: ", animal, " | ", method, " | dims=", dim, " res=", res)) +
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_color_viridis_d(option = "turbo")
        plot_list[[animal]] <- p

        # Cell count
        clust_tab <- table(Idents(obj))
        count_df <- data.frame(
          method = method,
          animal = animal,
          dims = dim,
          resolution = res,
          cluster = as.integer(names(clust_tab)),
          cell_count = as.vector(clust_tab)
        )

        # Silhouette
        pca_mat <- Embeddings(obj, "pca")[, dims_to_use]
        cluster_ids <- Idents(obj)

        sil_df <- if (length(unique(cluster_ids)) > 1) {
          sil <- silhouette(as.numeric(cluster_ids), dist(pca_mat))
          as.data.frame(sil) %>%
            mutate(cluster = as.integer(cluster_ids)) %>%
            group_by(cluster) %>%
            summarise(silhouette_score = mean(sil_width), .groups = "drop")
        } else {
          data.frame(cluster = as.integer(names(clust_tab)), silhouette_score = NA_real_)
        }

        count_df <- left_join(count_df, sil_df, by = "cluster")
        count_df$is_small_cluster <- count_df$cell_count < small_cluster_threshold

        # Average silhouette (entire animal/method/dim/res group)
        avg_sil <- tryCatch({
          if (length(unique(cluster_ids)) > 1) {
            sil <- silhouette(as.numeric(cluster_ids), dist(pca_mat))
            mean(sil[, "sil_width"], na.rm = TRUE)
          } else {
            NA_real_
          }
        }, error = function(e) NA_real_)

        count_df$average_silhouette_score <- avg_sil
        count_df$cumulative_variance <- cumulative_variance

        # Append
        cluster_summary <- rbind(cluster_summary, count_df)
      }

      # Save combined UMAP PDF
      pdf_file <- file.path(
        pdf_output_dir,
        paste0("individual_umap_", method, "_dims", dim, "_res", res, "_20250616.pdf")
      )
      ggsave(pdf_file, wrap_plots(plot_list, ncol = 2), width = 12, height = 10, dpi = 600)
      cat("✅ Saved:", pdf_file, "\n")
    }
  }
}

# Save final individual-level Seurat objects after all analyses
saveRDS(seurat_objects, file.path(rds_output_dir, "individual_seurat_analysis_20250616.rds"))
cat("✅ Saved individual Seurat objects with SCT/LogNorm clustering to RDS.\n")

# ------------------------- #
# Save Final Summary CSV
# ------------------------- #
csv_path <- file.path(tsv_output_dir, "indie_cluster_summary_sct_vs_log.csv")
write.csv(cluster_summary, csv_path, row.names = FALSE)
cat("✅ Final cluster summary saved to:", csv_path, "\n")




# ------------------------- #
# Merged Seurat Object - SCT Integration
# ------------------------- #

# Define dims/resolution for downstream UMAP/clustering
dims_list <- c(10, 15, 20, 25, 30)
resolutions <- c(0.25, 0.3, 0.5, 0.75)
small_cluster_threshold <- 100  # Customize as needed

cat("🔗 Integrating all animals using SCT...\n")
sct_features <- SelectIntegrationFeatures(object.list = seurat_objects_sct, nfeatures = 3000)
sct_prepped <- PrepSCTIntegration(object.list = seurat_objects_sct, anchor.features = sct_features)
sct_anchors <- FindIntegrationAnchors(object.list = sct_prepped, normalization.method = "SCT", anchor.features = sct_features)
merged_sct <- IntegrateData(anchorset = sct_anchors, normalization.method = "SCT")

DefaultAssay(merged_sct) <- "integrated"
merged_sct <- ScaleData(merged_sct, verbose = FALSE)
merged_sct <- RunPCA(merged_sct, npcs = max(dims_list))

saveRDS(merged_sct, file.path(rds_output_dir, "merged_totalNK_SCT_20250616.rds"))
cat("✅ Merged SCT-integrated Seurat object saved.\n")


# ------------------------- #
# Merged Seurat Object - LogNormalize (non-integrated)
# ------------------------- #
cat("📦 Merging all animals using LogNormalize (no integration)...\n")
lognorm_seurat_list <- lapply(seurat_objects, function(obj) {
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj
})

merged_lognorm <- merge(x = lognorm_seurat_list[[1]], y = lognorm_seurat_list[-1], add.cell.ids = names(lognorm_seurat_list))
merged_lognorm <- ScaleData(merged_lognorm)
merged_lognorm <- RunPCA(merged_lognorm, npcs = max(dims_list))

saveRDS(merged_lognorm, file.path(rds_output_dir, "merged_totalNK_LogNorm_20250616.rds"))
cat("✅ Merged LogNormalized Seurat object saved.\n")



# ------------------------- #
# Cluster Summary: Merged Seurat Object
# ------------------------- #

merged_summary <- data.frame()

for (method in c("SCT", "LogNormalize")) {
  cat("🔁 Summarizing merged for", method, "\n")
  merged_obj <- if (method == "SCT") merged_sct else merged_lognorm
  DefaultAssay(merged_obj) <- ifelse(method == "SCT", "integrated", "RNA")

  merged_obj <- ScaleData(merged_obj, verbose = FALSE)
  merged_obj <- RunPCA(merged_obj, npcs = max(dims_list), verbose = FALSE)

  var_explained <- merged_obj@reductions$pca@stdev^2
  var_explained <- var_explained / sum(var_explained)

  for (dim in dims_list) {
    dims_to_use <- 1:dim
    cumvar <- sum(var_explained[dims_to_use])

    merged_obj <- RunUMAP(merged_obj, dims = dims_to_use,
                          reduction.name = paste0("umap_", method, "_dims", dim))

    for (res in resolutions) {
      merged_obj <- FindNeighbors(merged_obj, dims = dims_to_use)
      merged_obj <- FindClusters(merged_obj, resolution = res)

      cluster_ids <- Idents(merged_obj)
      clust_tab <- table(cluster_ids)
      cl_ids <- as.integer(names(clust_tab))
      n_cells <- as.integer(clust_tab)
      avg_sil <- NA_real_

      # Silhouette calculation
      sil_df <- if (length(unique(cluster_ids)) > 1) {
        pca_mat <- Embeddings(merged_obj, "pca")[, dims_to_use]
        sil <- silhouette(as.numeric(cluster_ids), dist(pca_mat))
        avg_sil <- mean(sil[, "sil_width"], na.rm = TRUE)
        aggregate(sil[, "sil_width"], by = list(cluster = as.integer(cluster_ids)), FUN = mean) %>%
          rename(silhouette_score = x)
      } else {
        data.frame(cluster = cl_ids, silhouette_score = NA_real_)
      }

      # Main stats
      df <- data.frame(
        method = method,
        animal = "merged",
        dims = dim,
        resolution = res,
        cluster = cl_ids,
        cell_count = n_cells,
        is_small_cluster = n_cells < small_cluster_threshold,
        average_silhouette_score = avg_sil,
        cumulative_variance = cumvar
      )

      df <- left_join(df, sil_df, by = "cluster")
      merged_summary <- rbind(merged_summary, df)
    }
  }
}

# Save summary
summary_csv <- file.path(tsv_output_dir, "merged_cluster_summary_sct_vs_log.csv")
write.csv(merged_summary, summary_csv, row.names = FALSE)
cat("✅ Merged summary saved to:", summary_csv, "\n")


# Save both merged SCT and LogNormalized objects
merged_analysis_list <- list(
  SCT = merged_obj_sct,
  LogNormalize = merged_obj_log
)
saveRDS(merged_analysis_list, file.path(rds_output_dir, "merged_seurat_analysis_20250616.rds"))
cat("✅ Saved merged SCT and LogNormalized Seurat objects to RDS.\n")



# # ------------------------- #
# # UMAP & Clustering for Both Merged Datasets
# # ------------------------- #
# for (method in c("SCT", "LogNormalize")) {
#   cat("🚀 Running UMAP + clustering for merged:", method, "\n")
#   merged_obj <- if (method == "SCT") merged_sct else merged_lognorm
#   DefaultAssay(merged_obj) <- ifelse(method == "SCT", "integrated", "RNA")

#   merged_obj <- ScaleData(merged_obj, verbose = FALSE)
#   merged_obj <- RunPCA(merged_obj, npcs = max(dims_list))

#   for (dim in dims_list) {
#     dims_to_use <- 1:dim
#     merged_obj <- RunUMAP(merged_obj, dims = dims_to_use, reduction.name = paste0("umap_", method, "_dims", dim))

#     umap_plots <- list()

#     for (res in resolutions) {
#       merged_obj <- FindNeighbors(merged_obj, dims = dims_to_use)
#       merged_obj <- FindClusters(merged_obj, resolution = res)

#       p <- DimPlot(
#         merged_obj,
#         reduction = paste0("umap_", method, "_dims", dim),
#         group.by = "seurat_clusters",
#         label = TRUE,
#         pt.size = 0.3,
#         repel = TRUE
#       ) +
#         ggtitle(paste("Merged UMAP -", method, "- dims =", dim, "res =", res)) +
#         theme(plot.title = element_text(hjust = 0.5)) +
#         scale_color_viridis_d(option = "turbo")

#       umap_plots[[paste0("res_", res)]] <- p
#     }

#     # Save PDF with all resolutions for this dims
#     combined_plot <- wrap_plots(umap_plots, ncol = 2)
#     pdf_name <- file.path(pdf_output_dir, paste0("merged_umap_", method, "_dims", dim, "_multiRes_20250616.pdf"))
#     ggsave(
#       filename = pdf_name,
#       plot = combined_plot,
#       width = 12,
#       height = 4 * ceiling(length(resolutions) / 2),
#       dpi = 600
#     )
#     cat("✅ Saved:", pdf_name, "\n")
#   }
# }

