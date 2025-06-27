library(dplyr)
library(tidyr)
library(ggplot2)
# library(pheatmap)

# ------------------------ #
# Settings
# ------------------------ #
SCT_marker_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs/tsv/individual/SCT"
animal25 <- "animal25"
animal26 <- "animal26"
dims <- 15
res <- 0.5

# ------------------------ #
# Load Marker Files
# ------------------------ #
file1 <- file.path(SCT_marker_dir, paste0("markers_", animal25, "_SCT_dims", dims, "_res", res, ".csv"))
file2 <- file.path(SCT_marker_dir, paste0("markers_", animal26, "_SCT_dims", dims, "_res", res, ".csv"))

markers1 <- read.csv(file1, stringsAsFactors = FALSE)
markers2 <- read.csv(file2, stringsAsFactors = FALSE)

# ------------------------ #
# Split into cluster-wise DEG lists
# ------------------------ #
deg_list1 <- split(markers1$gene, markers1$cluster)
deg_list2 <- split(markers2$gene, markers2$cluster)

# ------------------------ #
# Compute Jaccard Index
# ------------------------ #
jaccard_matrix <- matrix(0, nrow = length(deg_list1), ncol = length(deg_list2),
                         dimnames = list(paste0(animal25, "_C", names(deg_list1)),
                                         paste0(animal26, "_C", names(deg_list2))))

for (i in names(deg_list1)) {
  for (j in names(deg_list2)) {
    set1 <- deg_list1[[i]]
    set2 <- deg_list2[[j]]
    intersect_size <- length(intersect(set1, set2))
    union_size <- length(union(set1, set2))
    jaccard_matrix[paste0(animal25, "_C", i), paste0(animal26, "_C", j)] <- ifelse(union_size > 0, intersect_size / union_size, 0)
  }
}