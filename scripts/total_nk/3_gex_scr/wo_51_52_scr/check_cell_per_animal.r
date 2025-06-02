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


library(dplyr)
# Assume clusters and animal IDs are in metadata
cluster_counts <- seurat_obj@meta.data %>%
 group_by(seurat_clusters, sample) %>%
 summarise(n = n(), .groups = "drop") %>%
 group_by(sample) %>%
 mutate(proportion = n / sum(n))
print(cluster_counts)

# 1. Perform Chi-squared test for independence of clusters distribution across animals (samples)
chisq.test(table(seurat_obj$seurat_clusters, seurat_obj$sample))

# data:  table(seurat_obj$seurat_clusters, seurat_obj$sample)
# X-squared = 10501, df = 42, p-value < 2.2e-16


# 2. Calculate average expression per animal (pseudo-bulk)
avg_expr <- AggregateExpression(seurat_obj, group.by = "sample")$SCT
cor_matrix <- cor(avg_expr)
print(cor_matrix)


# 3. Check variance of cell-type proportions across animals
prop_table <- table(seurat_obj$sample, seurat_obj$seurat_clusters)
prop_df <- prop.table(prop_table, margin = 1)  # normalize per animal
apply(prop_df, 2, sd)  # SD per cluster

#           0           1           2           3           4           5 
# 0.120321541 0.063248814 0.081289135 0.058945828 0.042528238 0.054919528 
#           6           7           8           9          10          11 
# 0.034764634 0.022417331 0.026684528 0.034422420 0.019506254 0.018779338 
#          12          13          14 
# 0.019958164 0.013978006 0.003946984 