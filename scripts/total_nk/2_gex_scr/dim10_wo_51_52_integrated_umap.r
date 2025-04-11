# Load libraries
library(Seurat)
library(SeuratDisk)
library(cluster)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(patchwork)
library(cowplot)
library(RColorBrewer)
library(pheatmap)

# 7. Visualize Integrated UMAP Excluding Animal51 and Animal52
# Load the integrated data object with clusters
integrated_data <- readRDS("/home/outputs/totalNK_outputs/2_dge/dim10/integrated_data10_clusters.rds")
output_dir <- "/home/outputs/totalNK_outputs/2_dge/dim10/"

# Subset the integrated data to exclude Animal51 and Animal52
samples_to_keep <- c("Animal25", "Animal26", "Animal27", "Animal28")
integrated_data_subset <- subset(integrated_data, subset = sample %in% samples_to_keep)

# Verify the subset
cat("Samples in subset:\n")
print(table(integrated_data_subset$sample))

# Generate and customize the UMAP plot for the subset
umap_plot <- DimPlot(integrated_data_subset,
                     group.by = "sample",
                     label = TRUE,  # Enable sample labels on the plot
                     label.size = 5,  # Adjust label size for readability
                     repel = TRUE,  # Repel labels to avoid overlap
                     pt.size = 0.3) +
             ggtitle("UMAP - Integration by Sample (Excluding Animal51 and Animal52)") +
             theme(plot.title = element_text(hjust = 0.5),
                   legend.position = "right") +
             scale_color_brewer(palette = "Set2")  # Consistent color palette

# Save the UMAP plot as a PDF in the specified output directory
output_file <- file.path(output_dir, "umap_integration_no_animal51_52.pdf")
ggsave(output_file,
       plot = umap_plot,
       width = 10,
       height = 8,
       dpi = 600,
       units = "in")

# Print completion message
cat("✅ UMAP plot excluding Animal51 and Animal52 saved as", output_file, "\n")