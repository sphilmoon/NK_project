install.packages("dplyr")
library(dplyr)
library(pheatmap)
library(stringr)
library(tidyr)
library(purrr)

# ------------------------ #
# Settings
# ------------------------ #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
pdf_output_dir <- file.path(output_dir, "pdf")
heatmap_output_dir <- file.path(pdf_output_dir, "3_dge_heatmap")
deg_dir_base <- file.path(output_dir, "tsv", "individual")


animals <- c("animal25", "animal26", "animal27", "animal28")
dims_list <- c(10, 15, 20, 25, 30)
resolutions <- c(0.25, 0.5, 0.75)
methods <- c("SCT", "LogNormalize")

# Create PDF output directories
for (method in methods) {
  dir.create(file.path(heatmap_output_dir, method), recursive = TRUE, showWarnings = FALSE)
}

# ------------------------ #
# Jaccard Calculation Function
# ------------------------ #
compute_and_save_jaccard <- function(animal1, animal2, dims, res, method) {
  marker_dir <- file.path(deg_dir_base, method)
  file1 <- file.path(marker_dir, paste0("markers_", animal1, "_", method, "_dims", dims, "_res", res, ".csv"))
  file2 <- file.path(marker_dir, paste0("markers_", animal2, "_", method, "_dims", dims, "_res", res, ".csv"))

  if (!file.exists(file1) || !file.exists(file2)) return(NULL)

  markers1 <- read.csv(file1, stringsAsFactors = FALSE)
  markers2 <- read.csv(file2, stringsAsFactors = FALSE)

  deg_list1 <- split(markers1$gene, markers1$cluster)
  deg_list2 <- split(markers2$gene, markers2$cluster)

  jaccard_matrix <- matrix(0, nrow = length(deg_list1), ncol = length(deg_list2),
                           dimnames = list(paste0(animal1, "_C", names(deg_list1)),
                                           paste0(animal2, "_C", names(deg_list2))))

  for (i in names(deg_list1)) {
    for (j in names(deg_list2)) {
      set1 <- deg_list1[[i]]
      set2 <- deg_list2[[j]]
      intersect_size <- length(intersect(set1, set2))
      union_size <- length(union(set1, set2))
      jaccard_matrix[paste0(animal1, "_C", i), paste0(animal2, "_C", j)] <-
        ifelse(union_size > 0, intersect_size / union_size, 0)
    }
  }

  # Save PDF
  pdf_file <- file.path(heatmap_output_dir, method,
                        paste0("jaccard_heatmap_", animal1, "_vs_", animal2,
                               "_dims", dims, "_res", res, ".pdf"))

  pdf(file = pdf_file, width = 8, height = 8)
  pheatmap(jaccard_matrix,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           main = paste("Jaccard Index Heatmap\n", animal1, "vs", animal2,
                        "| dims =", dims, "res =", res),
           color = colorRampPalette(c("white", "#6BAED6", "#08519C"))(100),
           display_numbers = TRUE,
           number_format = "%.2f")
  dev.off()
  cat("Saved:", pdf_file, "\n")
}

# ------------------------ #
# Iterate over combinations
# ------------------------ #
animal_pairs <- combn(animals, 2, simplify = FALSE)

for (method in methods) {
  for (dims in dims_list) {
    for (res in resolutions) {
      for (pair in animal_pairs) {
        compute_and_save_jaccard(pair[1], pair[2], dims, res, method)
      }
    }
  }
}

# Part 1. Generate Jaccard Matrices from DEG CSVs

animal_pairs <- combn(animals, 2, simplify = FALSE)

# ------------------------- #
# Initialize collector
# ------------------------- #
top_matches_all <- data.frame()

# ------------------------- #
# Jaccard from DEG files
# ------------------------- #
for (method in methods) {
  deg_dir <- file.path(deg_dir_base, method)

  for (dims in dims_list) {
    for (res in resolutions) {
      for (pair in animal_pairs) {
        a1 <- pair[1]; a2 <- pair[2]

        f1 <- file.path(deg_dir, paste0("markers_", a1, "_", method, "_dims", dims, "_res", res, ".csv"))
        f2 <- file.path(deg_dir, paste0("markers_", a2, "_", method, "_dims", dims, "_res", res, ".csv"))

        if (!file.exists(f1) || !file.exists(f2)) next

        deg1 <- read.csv(f1, stringsAsFactors = FALSE)
        deg2 <- read.csv(f2, stringsAsFactors = FALSE)

        clusters1 <- split(deg1$gene, deg1$cluster)
        clusters2 <- split(deg2$gene, deg2$cluster)

        for (c1 in names(clusters1)) {
          g1 <- unique(clusters1[[c1]])
          for (c2 in names(clusters2)) {
            g2 <- unique(clusters2[[c2]])

            jaccard <- length(intersect(g1, g2)) / length(union(g1, g2))
            if (jaccard > 0.3) {
              top_matches_all <- rbind(top_matches_all, data.frame(
                method = method,
                dims = dims,
                res = res,
                animal1 = a1,
                animal2 = a2,
                cluster1 = paste0(a1, "_C", c1),
                cluster2 = paste0(a2, "_C", c2),
                jaccard = jaccard
              ))
            }
          }
        }
      }
    }
  }
}

# Save output
write.csv(top_matches_all,
          file.path(deg_dir_base, "1_summary_top_cluster_matches_from_csv.csv"),
          row.names = FALSE)

# ------------------------- #
# Part 2. Count strong overlaps per combo
# ------------------------- #
summary_counts <- top_matches_all %>%
  group_by(method, dims, res) %>%
  summarise(n_strong_matches = n(), .groups = "drop")

write.csv(summary_counts,
          file.path(deg_dir_base, "2_summary_overlap_counts_per_combo.csv"),
          row.names = FALSE)

# ------------------------- #
# Part 3.  Identify stable clusters across comparisons
# ------------------------- #
stable_clusters <- top_matches_all %>%
  pivot_longer(cols = c("cluster1", "cluster2"), names_to = "source", values_to = "cluster") %>%
  count(cluster, sort = TRUE) %>%
  filter(n >= 3)

write.csv(stable_clusters,
          file.path(deg_dir_base, "3_summary_stable_clusters.csv"),
          row.names = FALSE)

cat("✅ Finished processing. All summaries saved to:\n", deg_dir_base, "\n")