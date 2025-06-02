library(dplyr)
# Ensure identity is set to your cluster resolution
Idents(seurat_obj) <- "seurat_clusters"
# Extract metadata into a data frame
meta <- seurat_obj@meta.data
# Summarize counts by cluster, animal, and condition
cell_counts <- meta %>%
 group_by(seurat_clusters, sample) %>%
 summarise(n_cells = n(), .groups = "drop")
# Print to console
print(cell_counts)
# Cross-tabulation: cluster × animal × condition
table(Idents(seurat_obj), seurat_obj$sample)