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

# (Previous script sections remain unchanged up to the last line)

# 7. Visualize Integrated UMAP Excluding Animal51 and Animal52
# Load the integrated data object with clusters
integrated_data <- readRDS("/home/outputs/totalNK_outputs/2_dge/dim10/integrated_data10_clusters.rds")

# Subset the integrated data to exclude Animal51 and Animal52
samples_to_keep <- c("Animal25", "Animal26", "Animal27", "Animal28")
integrated_data_subset <- subset(integrated_data, subset = sample %in% samples_to_keep)

# Verify the subset
cat("Samples in subset:\n")
print(table(integrated_data_subset$sample))

# Generate and customize the UMAP plot for the subset
umap_plot <- DimPlot(integrated_data_subset, 
                     group.by = "sample", 
                     label = FALSE, 
                     pt.size = 0.5) +
             ggtitle("UMAP - Integration by Sample (Excluding Animal51 and Animal52)") +
             theme(plot.title = element_text(hjust = 0.5),
                   legend.position = "right") +
             scale_color_brewer(palette = "Set2")  # Consistent color palette

# Save the UMAP plot as a PNG
ggsave("umap_integration_no_animal51_52.pdf", 
       plot = umap_plot, 
       width = 10, 
       height = 8, 
       dpi = 600, 
       units = "in")

# Print completion message
cat("✅ UMAP plot excluding Animal51 and Animal52 saved as umap_integration_no_animal51_52.pdf\n")

