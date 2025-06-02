# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(Matrix)

# ------------------------- #
# Define Paths and Configs
# ------------------------- #
totalNK_rds <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/rds/integrated_data_dims25_res0.3_genecounts.rds"
nkp46_rds <- "/home/outputs/nkp46_outputs/nkp46_integrated_data.rds"
output_dir <- "/home/outputs/all_merged_TotalNK_nkp46"

rds_dir <- file.path(output_dir, "old/rds")
pdf_dir <- file.path(output_dir, "old/old_pdf")

# ------------------------- #
# Load Merged Seurat Object
# ------------------------- #
merged_obj <- readRDS(file.path(rds_dir, "all_totalNK_nkp46.rds"))
cat("✅ Loaded merged Seurat object\n")


how can i verify the column names of the merged object?
# Check column names of the merged object
cat("Column names in merged object:\n")
print(colnames(merged_obj))

# Check sample_id and animal_id in the merged object
sample_ids <- unique(merged_obj$sample_id)
animal_ids <- unique(merged_obj$animal)
cat("Sample IDs in merged object:", paste(sample_ids, collapse = ", "), "\n")
cat("Animal IDs in merged object:", paste(animal_ids, collapse = ", "), "\n")   






# ------------------------- #
# Define Custom UMAP Theme with Black Axes
# ------------------------- #
umap_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    legend.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )


# ------------------------- #
# Create Combined UMAP Visualization
# ------------------------- #
cat("🎨 Creating combined UMAP visualization...\n")

combined_umap <- DimPlot(merged_obj, group.by = "sample_id", pt.size = 0.5) +
  ggtitle("Combined UMAP of Total NK, NKp46+, and NKp46- (dims25, res 0.5)") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  umap_theme

# Save combined UMAP
combined_output_file <- file.path(figures_dir, "merged_UMAP_TotalNK_NKp46_d25_res0.5.pdf")
ggsave(filename = combined_output_file, plot = combined_umap, width = 10, height = 8, dpi = 600, bg = "transparent")
cat("✅ Combined UMAP saved to", combined_output_file, "\n")

# ------------------------- #
# Create Separate UMAPs for Each Sample Identity
# ------------------------- #
cat("🎨 Creating separate UMAP visualizations for each sample identity...\n")

# Split UMAP by sample_id
split_umap <- DimPlot(merged_obj, group.by = "sample_id", split.by = "sample_id", pt.size = 0.5, ncol = 3) +
  ggtitle("UMAP Split by Sample Identity (dims25, res 0.5)") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  umap_theme

# Save split UMAPs
split_output_file <- file.path(figures_dir, "split_UMAP_TotalNK_NKp46_d25_res0.5.pdf")
ggsave(filename = split_output_file, plot = split_umap, width = 15, height = 5, dpi = 600, bg = "transparent")
cat("✅ Split UMAPs saved to", split_output_file, "\n")

cat("🎉 UMAP generation complete. Outputs saved in:\n", output_dir, "\n")












here is my metadata information: 
> sample_ids <- unique(merged_obj$sample_id)
ids <- unique(merged_obj$animal)> animal_ids <- unique(merged_obj$animal)
> cat("Sample IDs in merged object:", paste(sample_ids, collapse = ", "), "\n")
Sample IDs in merged object: TotalNK, NKp46+, NKp46- 
> cat("Animal IDs in merged object:", paste(animal_ids, collapse = ", "), "\n")  
Animal IDs in merged object: NA, Animal25, Animal26, Animal27, Animal28, Animal52 

I want to see each animal displayed in a different color in the UMAP plot from each sample_id. 
you can generate three different umap plots, one for each animal, and then combine them into a single plot.






library(ggplot2)
library(patchwork)
# Ensure animal is a factor with NA removed
merged_obj$animal <- as.factor(merged_obj$animal)
# Get list of sample IDs
sample_ids <- unique(merged_obj$sample_id)
# Generate a named list of plots for each sample_id
umap_list <- lapply(sample_ids, function(sid) {
 # Subset object
 obj_sub <- subset(merged_obj, subset = sample_id == sid & !is.na(animal))
 # Plot UMAP colored by animal
 DimPlot(obj_sub, group.by = "animal", reduction = "umap") +
   ggtitle(paste("Sample:", sid)) +
   theme(plot.title = element_text(hjust = 0.5))
})
# Combine plots using patchwork
combined_plot <- wrap_plots(umap_list, ncol = 1)
combined_plot

# Save the combined plot
output_file <- file.path(pdf_dir, "combined_umap_by_animal.pdf")
ggsave(filename = output_file, plot = combined_plot, width = 10, height = 15, dpi = 600, bg = "transparent")
cat("✅ Combined UMAP by animal saved to", output_file, "\n")

