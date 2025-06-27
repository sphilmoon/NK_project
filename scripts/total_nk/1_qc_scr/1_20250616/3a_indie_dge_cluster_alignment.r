install.packages("dplyr")
library(dplyr)
library(pheatmap)

# ------------------------ #
# Settings
# ------------------------ #
output_dir <- "/home/outputs/totalNK_outputs/1_qc/1_20250616_outs"
pdf_output_dir <- file.path(output_dir, "pdf")
heatmap_output_dir <- file.path(pdf_output_dir, "3_dge_heatmap")
tsv_dir_base <- file.path(output_dir, "tsv", "individual")

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
  marker_dir <- file.path(tsv_dir_base, method)
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