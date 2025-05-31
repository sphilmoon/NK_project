# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(ggplot2)
library(cowplot)
library(stringr)
library(cluster) # For silhouette score
library(dplyr) # For data frame manipulation

# ------------------------- #
# Define Input and Output Paths
# ------------------------- #
nkp_neg_rds <- "/home/outputs/nkp46_outputs/chat/neg/rds/nkp46_neg.rds"
nkp_pos_rds <- "/home/outputs/nkp46_outputs/chat/pos/rds/nkp46_pos.rds"
totalNK_rds <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/rds/integrated_data_dims25_res0.3_genecounts.rds"

output_dir <- "/home/outputs/all_merged_TotalNK_nkp46/20250531"
all_merged_rds_dir <- file.path(output_dir, "rds")
pdf_dir <- file.path(output_dir, "pdf")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(all_merged_rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

# # ------------------------- #
# # Load Seurat Objects
# # ------------------------- #
cat("📦 Loading Seurat objects...\n")
nkp_neg <- readRDS(nkp_neg_rds)
cat("✅ Loaded nkp46_neg\n")
nkp_pos <- readRDS(nkp_pos_rds)
cat("✅ Loaded nkp46_pos\n")
totalNK <- readRDS(totalNK_rds)
cat("✅ Loaded totalNK\n")

# ------------------------- #
# Standardize Metadata
# ------------------------- #
cat("🔧 Standardizing metadata (sample, sample_id, condition, animal_id)...\n")

for (obj_name in c("nkp_neg", "nkp_pos", "totalNK")) {
    obj <- get(obj_name)

    # Convert to title case
    obj$sample <- str_to_title(obj$sample)

    # Copy fields
    obj$sample_id <- obj$sample
    obj$animal_id <- obj$sample

    # Set condition
    obj$condition <- switch(obj_name,
        "nkp_neg" = "nkp46_neg",
        "nkp_pos" = "nkp46_pos",
        "totalNK" = "totalNK"
    )

    assign(obj_name, obj)

    # Verify
    cat("🔍 Checking 'sample' field in ", obj_name, "\n")
    print(colnames(obj@meta.data))
    cat("✅ 'sample' values:\n")
    print(head(obj$sample))
}

# ------------------------- #
# Merge Seurat Objects
# ------------------------- #
# ✅ Rename cells to avoid name collisions
nkp_neg <- RenameCells(nkp_neg, add.cell.id = "nkp46neg")
nkp_pos <- RenameCells(nkp_pos, add.cell.id = "nkp46pos")
totalNK <- RenameCells(totalNK, add.cell.id = "totalNK")

# Remove SCTransform normalization separately to enable merging with a new normalization
nkp_neg[["SCT"]] <- NULL
nkp_pos[["SCT"]] <- NULL
totalNK[["SCT"]] <- NULL

# 🔗 Merge all three objects safely
cat("🔗 Merging Seurat objects...\n")
merged_obj <- merge(nkp_neg, y = list(nkp_pos, totalNK), merge.data = FALSE)

DefaultAssay(merged_obj) <- "RNA"
cat("✅ Merge complete: total cells =", ncol(merged_obj), "\n")

# ------------------------- #
# Normalize, PCA, and Save RDS
# ------------------------- #
cat("⚙️ Running SCTransform normalization and scaling merged object...\n")
merged_obj <- SCTransform(merged_obj, verbose = TRUE)

cat("🔬 Running PCA...\n")
merged_obj <- RunPCA(merged_obj, npcs = 50)
cat("✅ PCA complete\n")

# Save RDS file after SCTransform and PCA
cat("💾 Saving merged object after SCTransform and PCA...\n")
saveRDS(merged_obj, file = file.path(all_merged_rds_dir, "merged_totalNK_nkp46_post_sct_pca_20250531.rds"))
cat("✅ Saved merged object as 'merged_totalNK_nkp46_post_sct_pca.rds'\n")

# ------------------------- #
# ElbowPlot for Dimensionality Selection
# ------------------------- #
cat("📊 Saving ElbowPlot for dimensionality selection...\n")
elbow_plot <- ElbowPlot(merged_obj, ndims = 50)
ggsave(file.path(pdf_dir, "elbow_plot.pdf"), elbow_plot, width = 8, height = 6, dpi = 600)


cat("📦 Loading saved Seurat object...\n")
merged_obj <- readRDS(file.path(all_merged_rds_dir, "merged_totalNK_nkp46_post_sct_pca_20250531.rds"))
cat("✅ Loaded merged object\n")




# # ------------------------- #
# # Test Different Dims and Resolutions with Silhouette Scores
# # ------------------------- #
# # Define the dimensionalities and resolutions to test
# dims_list <- c(20, 25, 30, 33)
# resolution_list <- c(0.25, 0.5, 0.75)

# # Initialize a data frame to store silhouette scores
# silhouette_results <- data.frame(
#     Dims = integer(),
#     Resolution = numeric(),
#     Avg_Silhouette_Score = numeric(),
#     stringsAsFactors = FALSE
# )

# cat("🔍 Testing different dimensions and resolutions...\n")
# for (dims in dims_list) {
#     dims_to_use <- 1:dims

#     # Run UMAP with the current dims
#     cat(sprintf("🧭 Running UMAP with dims 1:%d...\n", dims))
#     merged_obj <- RunUMAP(merged_obj, dims = dims_to_use, reduction = "pca")

#     # Run FindNeighbors with the current dims
#     cat(sprintf("🔗 Running FindNeighbors with dims 1:%d...\n", dims))
#     merged_obj <- FindNeighbors(merged_obj, dims = dims_to_use, reduction = "pca")

#     for (res in resolution_list) {
#         # Run FindClusters with the current resolution
#         cat(sprintf("🗂️ Running FindClusters with resolution %.2f...\n", res))
#         merged_obj <- FindClusters(merged_obj, resolution = res)

#         # Get cluster labels
#         cluster_labels <- as.numeric(merged_obj$seurat_clusters)

#         # Extract PCA embeddings for silhouette score calculation
#         pca_embeddings <- Embeddings(merged_obj, reduction = "pca")[, dims_to_use]

#         # Subsample cells for silhouette score (e.g., 10,000 cells)
#         set.seed(123) # For reproducibility
#         sample_cells <- sample(1:ncol(pca_embeddings), min(10000, ncol(pca_embeddings)))
#         pca_embeddings_sub <- pca_embeddings[, sample_cells]
#         cluster_labels_sub <- cluster_labels[sample_cells]

#         # Compute silhouette scores
#         cat(sprintf("📈 Computing silhouette score for dims 1:%d, res %.2f (subsample)... \n", dims, res))
#         sil_scores <- silhouette(cluster_labels, dist(pca_embeddings))

#         # Calculate average silhouette score
#         avg_sil_score <- mean(sil_scores[, "sil_width"])

#         # Store the results
#         silhouette_results <- rbind(
#             silhouette_results,
#             data.frame(
#                 Dims = dims,
#                 Resolution = res,
#                 Avg_Silhouette_Score = avg_sil_score
#             )
#         )

#         # Save UMAP plot for this combination
#         umap_plot <- DimPlot(merged_obj, group.by = "seurat_clusters", pt.size = 0.3) +
#             ggtitle(sprintf("UMAP: Dims 1:%d, Res %.2f", dims, res)) +
#             theme(plot.title = element_text(hjust = 0.5))
#         ggsave(
#             file.path(pdf_dir, sprintf("UMAP_dims%d_res%.2f.pdf", dims, res)),
#             umap_plot,
#             width = 10, height = 8, dpi = 600
#         )

#         cat(sprintf("✅ UMAP plot saved for dims 1:%d, res %.2f\n", dims, res))
#     }
# }

# # ------------------------- #
# # Save Silhouette Scores to CSV
# # ------------------------- #
# cat("💾 Saving silhouette scores to CSV...\n")
# write.csv(
#     silhouette_results,
#     file = file.path(output_dir, "silhouette_scores.csv"),
#     row.names = FALSE
# )
# cat("✅ Silhouette scores saved as 'silhouette_scores.csv'\n")





# cat("From now on, it is automated umap clustering... \n")
# cat("Need to start from 'merged_totalNK_nkp46_post_sct_pca_20250531.rds' for umap re-clustering... \n")

# # ------------------------- #
# # Final Clustering with Selected Parameters
# # ------------------------- #
# # Use the best dims and resolution based on silhouette score
# best_params <- silhouette_results %>%
#     arrange(desc(Avg_Silhouette_Score)) %>%
#     slice(1)
# best_dims <- best_params$Dims
# best_res <- best_params$Resolution

# cat(sprintf(
#     "🌟 Best parameters: Dims 1:%d, Resolution %.2f (Silhouette Score: %.3f)\n",
#     best_dims, best_res, best_params$Avg_Silhouette_Score
# ))

# dims_to_use <- 1:best_dims
# merged_obj <- RunUMAP(merged_obj, dims = dims_to_use)
# merged_obj <- FindNeighbors(merged_obj, dims = dims_to_use)
# merged_obj <- FindClusters(merged_obj, resolution = best_res)

# # Save the final merged object with the best parameters
# saveRDS(merged_obj, file = file.path(all_merged_rds_dir, sprintf("merged_totalNK_nkp46_dim%d_res%.2f.rds", best_dims, best_res)))
# cat(sprintf("💾 Final merged object saved as 'merged_totalNK_nkp46_dim%d_res%.2f.rds'\n", best_dims, best_res))

# # ------------------------- #
# # UMAP by Condition
# # ------------------------- #
# cat("🎨 Generating combined UMAP by condition...\n")
# umap_by_condition <- DimPlot(merged_obj, group.by = "condition", pt.size = 0.3) +
#     ggtitle("UMAP by Condition") +
#     theme(plot.title = element_text(hjust = 0.5))

# ggsave(file.path(pdf_dir, "UMAP_by_condition_combined.pdf"),
#     umap_by_condition,
#     width = 10, height = 8, dpi = 600
# )

# # ------------------------- #
# # UMAP Grid: Animal x Condition
# # ------------------------- #
# cat("🎨 Generating UMAP grid by animal x condition...\n")
# conditions <- c("nkp46_neg", "nkp46_pos", "totalNK")
# animals <- c("Animal25", "Animal26", "Animal27", "Animal28")
# plot_list <- list()

# for (animal in animals) {
#     for (cond in conditions) {
#         cells <- WhichCells(merged_obj, expression = animal_id == animal & condition == cond)
#         cat(sprintf("🔍 Found %d cells for %s x %s\n", length(cells), animal, cond))

#         if (length(cells) == 0) {
#             cat("⚠️ No cells found for", animal, cond, "- skipping\n")
#             p <- ggplot() +
#                 theme_void()
#         } else {
#             p <- DimPlot(merged_obj, cells = cells, group.by = "condition", pt.size = 0.3) +
#                 ggtitle(NULL) +
#                 theme_void() +
#                 theme(legend.position = "none")
#         }
#         plot_list[[paste(animal, cond, sep = "_")]] <- p
#     }
# }

# # Add column titles (condition names)
# column_titles <- lapply(conditions, function(label) {
#     ggplot() +
#         annotate("text", x = 0.5, y = 0.5, label = label, size = 6, fontface = "bold") +
#         theme_void()
# })
# top_row <- wrap_plots(column_titles, ncol = length(conditions))

# # Add row titles (animal IDs)
# row_titles <- lapply(animals, function(animal) {
#     ggplot() +
#         annotate("text", x = 0.5, y = 0.5, label = animal, angle = 90, size = 6, fontface = "bold") +
#         theme_void()
# })
# left_col <- wrap_plots(row_titles, ncol = 1)

# # Assemble UMAP grid
# umap_matrix <- lapply(animals, function(animal) {
#     plots_row <- lapply(conditions, function(cond) {
#         plot_list[[paste(animal, cond, sep = "_")]]
#     })
#     wrap_plots(plots_row, ncol = length(conditions))
# })
# umap_full <- wrap_plots(umap_matrix, ncol = 1)

# # Add labels: top (column), left (row)
# umap_labeled <- plot_grid(
#     plot_grid(NULL, top_row, ncol = 2, rel_widths = c(0.12, 1)),
#     plot_grid(left_col, umap_full, ncol = 2, rel_widths = c(0.12, 1)),
#     nrow = 2, rel_heights = c(0.08, 1)
# )

# # Export UMAP grid
# ggsave(file.path(pdf_dir, "UMAP_animal_condition_grid.pdf"),
#     umap_labeled,
#     width = 12, height = 16, dpi = 600
# )

# cat("✅ All UMAP visualizations complete\n")
