# Load libraries
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(patchwork)
library(cowplot)
library(RColorBrewer)
library(pheatmap)
 
# Define file paths and sample names
file_paths <- list(
  "/home/rawdata/total_nk/animal25_totalnk_filtered_feature_bc_matrix.h5",
  "/home/rawdata/total_nk/animal27_totalnk_filtered_feature_bc_matrix.h5",
  "/home/rawdata/total_nk/animal51_totalnk_filtered_feature_bc_matrix.h5",
  "/home/rawdata/total_nk/animal26_totalnk_filtered_feature_bc_matrix.h5",
  "/home/rawdata/total_nk/animal28_totalnk_filtered_feature_bc_matrix.h5",
  "/home/rawdata/total_nk/animal52_totalnk_filtered_feature_bc_matrix.h5"
)
sample_names <- c("Animal25", "Animal27", "Animal51", "Animal26", "Animal28", "Animal52")
 
# 1. Load individual H5 files into Seurat objects
seurat_objects <- mapply(function(file, name) {
  data <- Read10X_h5(file)
  seurat_obj <- CreateSeuratObject(counts = data, project = name, min.cells = 3, min.features = 200)
  seurat_obj$sample <- name  # Add sample metadata
  return(seurat_obj)
}, file_paths, sample_names, SIMPLIFY = FALSE)
 
# Assign meaningful names to the list
names(seurat_objects) <- sample_names
 
# Print sample information with details
for (name in names(seurat_objects)) {
  cat("Sample:", name, "\n")
  print(seurat_objects[[name]])
  cat("Number of cells:", ncol(seurat_objects[[name]]), "\n")
  cat("Number of genes:", nrow(seurat_objects[[name]]), "\n\n")
}
 
# Optional: Save the list for later use
saveRDS(seurat_objects, "seurat_objects_list.rds")
 
# 2. quality control and pre-processing
# Load the saved RDS file
seurat_objects <- readRDS("seurat_objects_list.rds")
 
# Verify the loaded objects (optional)
cat("Number of samples loaded:", length(seurat_objects), "\n")
lapply(names(seurat_objects), function(name) {
  cat("Sample:", name, " - Cells:", ncol(seurat_objects[[name]]), " Genes:", nrow(seurat_objects[[name]]), "\n")
})

# to see the gene names
head(rownames(seurat_objects[[1]]), 20)

# Create a list to store all combined plots
combined_plots <- list()

# Perform QC filtering and generate violin plots
for (i in 1:length(seurat_objects)) {
  # Define animal name
  animal_name <- names(seurat_objects)[i]
  
  # Apply QC filtering based on animal name
  if (animal_name == "Animal27") {
    seurat_objects[[i]] <- subset(seurat_objects[[i]], 
                                  subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & 
                                           nCount_RNA > 250 & nCount_RNA < 7500)
  } else if (animal_name == "Animal28") {
    seurat_objects[[i]] <- subset(seurat_objects[[i]], 
                                  subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & 
                                           nCount_RNA > 250 & nCount_RNA < 10000)
  } else if (animal_name == "Animal51") {
    seurat_objects[[i]] <- subset(seurat_objects[[i]], 
                                  subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & 
                                           nCount_RNA > 250 & nCount_RNA < 5000)
  } else if (animal_name == "Animal52") {
    seurat_objects[[i]] <- subset(seurat_objects[[i]], 
                                  subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & 
                                           nCount_RNA > 250 & nCount_RNA < 7500)
  } else {
    # Default threshold for other animals (e.g., Animal25, Animal26)
    seurat_objects[[i]] <- subset(seurat_objects[[i]], 
                                  subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & 
                                           nCount_RNA > 500 & nCount_RNA < 10000)
  }
  
  # Create separate violin plots for each feature, removing the legend
  p1 <- VlnPlot(seurat_objects[[i]], features = "nFeature_RNA", ncol = 1) +
        ggtitle("Gene counts QC") +
        theme(plot.title = element_text(hjust = 0.5),  # Center the title
              legend.position = "none")               # Remove legend
  
  p2 <- VlnPlot(seurat_objects[[i]], features = "nCount_RNA", ncol = 1) +
        ggtitle("UMI counts QC") +
        theme(plot.title = element_text(hjust = 0.5),  # Center the title
              legend.position = "none")               # Remove legend
  
  # Combine plots for this sample using patchwork
  combined_plot <- p1 + p2 + plot_layout(ncol = 2)
  
  # Add an overall title for the combined plot
  combined_plot <- combined_plot + plot_annotation(title = paste("QC for", animal_name))
  
  # Store the combined plot in the list
  combined_plots[[animal_name]] <- combined_plot
  
  # Print status message for this sample
  cat("✅ Plot generated for", animal_name, "\n")
}

# Stack all plots vertically using cowplot
stacked_plot <- cowplot::plot_grid(plotlist = combined_plots, ncol = 1)

# Save the stacked plot as a single PNG
output_file <- "qc_violin_stacked.png"
ggsave(output_file, plot = stacked_plot, width = 10, height = 6 * length(seurat_objects), dpi = 600, units = "in")

# Print final message
cat("✅ All QC plots combined into", output_file, "\n")

# Normalize with SCTransform
seurat_objects <- lapply(seurat_objects, SCTransform, verbose = FALSE)
 
# Save the QC'd and normalized objects
saveRDS(seurat_objects, "seurat_objects_qc_sct.rds")
 
# Final completion message
cat("QC and figure generation process is finished.\n")


# 3. Batch correction and integration
# Load the saved RDS file
seurat_objects <- readRDS("seurat_objects_qc_sct.rds")

# Step 1: Run PCA on the first sample to determine dims
# Create a list to store all elbow plots
combined_elbow_plots <- list()

for (i in 1:length(seurat_objects)) {
  # Define animal name
  animal_name <- names(seurat_objects)[i]
  
  # Run PCA on the current sample
  seurat_objects[[i]] <- RunPCA(seurat_objects[[i]], verbose = FALSE)
  
  # Generate and customize the elbow plot
  elbow_plot <- ElbowPlot(seurat_objects[[i]], ndims = 50) +
                ggtitle(paste("Elbow Plot for PCA -", animal_name)) +  # Dynamic title
                theme(plot.title = element_text(hjust = 0.5))  # Center the title
  
  # Store the plot in the list
  combined_elbow_plots[[animal_name]] <- elbow_plot
  
  # Print status message
  cat("✅ Elbow plot generated for", animal_name, "\n")
}

stacked_elbow_plot <- cowplot::plot_grid(plotlist = combined_elbow_plots, ncol = 1)

output_file <- "combined_elbow_plots.png"
ggsave(output_file, plot = stacked_elbow_plot, 
       width = 10, height = 6 * length(seurat_objects), dpi = 600, units = "in")

# Print final message
cat("✅ All elbow plots combined into", output_file, "\n")
# Save the plot
ggsave("elbow_plot.png", plot = elbow_plot, width = 10, height = 6, dpi = 600, units = "in")

# Optional: Display the plot
# print(elbow_plot)

# Step 2: Determine the number of dimensions for integration
dims_to_use <- 1:20  # Replace with your chosen dims

# Step 3: Select integration features
features <- SelectIntegrationFeatures(object.list = seurat_objects, nfeatures = 3000)

# Step 4: Find integration anchors
anchors <- FindIntegrationAnchors(object.list = seurat_objects, anchor.features = features, dims = dims_to_use)

# Step 5: Integrate data
integrated_data <- IntegrateData(anchorset = anchors, dims = dims_to_use)

# 4. Clustering and visualizing the integration
DefaultAssay(integrated_data) <- "integrated"
integrated_data <- ScaleData(integrated_data)
integrated_data <- RunPCA(integrated_data)
integrated_data <- RunUMAP(integrated_data, dims = dims_to_use)

# Save the integrated data object
saveRDS(integrated_data, "integrated_data.rds")

# Generate and customize the DimPlot
dim_plot <- DimPlot(integrated_data, group.by = "sample", label = TRUE, label = FALSE, pt.size = 0.5) +
            ggtitle("UMAP - Integration by Sample") +
            theme(plot.title = element_text(hjust = 0.5),
                  legend.position = "right") +  # Adjust legend position (or "none" to remove)
            scale_color_brewer(palette = "Set2")  # Optional: Change color palette

# Save the DimPlot as a PNG
ggsave("dimplot_integration.png", plot = dim_plot, width = 10, height = 8, dpi = 600, units = "in")


# 5. Clustering to Identify Cell Groups
# Load your integrated Seurat object
integrated_data <- readRDS("integrated_data.rds")

# Set default assay now that all six samples are integrated with batch correction
DefaultAssay(integrated_data) <- "integrated"

# Now identifying clusters with the integrated data and then identify DEGs for each cluster
integrated_data <- FindNeighbors(integrated_data, dims = dims_to_use)
integrated_data <- FindClusters(integrated_data, resolution = 0.5)

# Save the integrated data object with multiple clusters identified across six samples
saveRDS(integrated_data, "integrated_data_clusters.rds")

# Visualize the clusters
cluster_plot <- DimPlot(integrated_data, label = TRUE, group.by = "seurat_clusters")
ggsave("cluster_plot.png", plot = cluster_plot, width = 10, height = 8, dpi = 600, units = "in")

# Validate the integrated clusters with well-known Markers
feature_plot <- FeaturePlot(integrated_data, features = c("NCR1", "CD3E"), ncol = 2)
ggsave("feature_plot_markers.png", plot = feature_plot, width = 12, height = 6, dpi = 600, units = "in")




# 6. Differential Expression Analysis for finding commonly expressed genes per cluster
# Load the integrated data object with clusters
integrated_data <- readRDS("integrated_data_clusters.rds")
# Switching to the RNA assay for DEG analysis
DefaultAssay(integrated_data) <- "RNA"

# Join layers in the RNA assay
integrated_data <- JoinLayers(integrated_data, assay = "RNA")
# Normalize a log-normalization on the integrated data for DEG analysis
integrated_data <- NormalizeData(integrated_data, assay = "RNA")

# Find DEGs in all clusters across six samples
degs_all_clusters <- FindAllMarkers(integrated_data, 
                                    assay = "RNA",
                                    min.pct = 0.25, 
                                    logfc.threshold = 0.25, 
                                    test.use = "wilcox")
# Filter DEGs based on adjusted p-value and log2 fold change
degs_filtered <- degs_all_clusters %>%
                filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)

# Identify shared DEGs across all samples
# Get the list of clusters
clusters <- unique(degs_filtered$cluster)

# Initialize a list to store shared DEGs
shared_degs_list <- list()

# Loop over each cluster
for (cl in clusters) {
  # Subset DEGs for this cluster
  degs_cluster <- degs_filtered %>% filter(cluster == cl)
  genes <- unique(degs_cluster$gene)
  
  # Subset cells in this cluster
  cells_in_cluster <- WhichCells(integrated_data, idents = cl)
  cluster_data <- subset(integrated_data, cells = cells_in_cluster)
  
  # Get average expression by sample
  avg_exp <- AverageExpression(cluster_data, assays = "RNA", features = genes, group.by = "sample")$RNA
  
  # Check if the DEG is consistently up- or down-regulated
  shared_genes <- c()
  for (gene in genes) {
    gene_logfc <- degs_cluster %>% filter(gene == !!gene) %>% pull(avg_log2FC)
    direction <- ifelse(gene_logfc > 0, "up", "down")
    
    # Get expression values across samples
    exp_values <- avg_exp[gene, ]
    
    # Check consistency: all samples should have the same direction as the overall logFC
    if (direction == "up") {
      consistent <- all(exp_values > 0)  # All samples should have positive expression relative to mean
    } else {
      consistent <- all(exp_values < 0)  # All samples should have negative expression relative to mean
    }
    
    if (consistent) {
      shared_genes <- c(shared_genes, gene)
    }
  }
  
  # Store shared DEGs for this cluster
  shared_degs_list[[as.character(cl)]] <- degs_cluster %>% filter(gene %in% shared_genes)
}

# Combine shared DEGs across all clusters
shared_degs <- do.call(rbind, shared_degs_list)

# a. visualzing the shared DEGs
# Select top 5 DEGs per cluster (if too many)
top_degs <- shared_degs %>%
            group_by(cluster) %>%
            top_n(5, wt = abs(avg_log2FC)) %>%
            ungroup() %>%
            pull(gene) %>%
            unique()

# Ensure the RNA assay is scaled
integrated_data <- ScaleData(integrated_data, assay = "RNA")

# Generate dot plot with a custom color palette
dot_plot <- DotPlot(integrated_data, 
                    features = top_degs, 
                    group.by = "seurat_clusters", 
                    split.by = "sample", 
                    assay = "RNA",
                    cols = RColorBrewer::brewer.pal(6, "Set2")) +  # Use a palette with 6 colors
            ggtitle("Dot Plot of Shared DEGs Across All Samples") +
            theme(plot.title = element_text(hjust = 0.5))
# Save the dot plot as a PNG
ggsave("dotplot_shared_degs.png", plot = dot_plot, width = 14, height = 8, dpi = 600, units = "in")

# generateing violin plot
# Select a few DEGs for visualization
top3_degs <- head(top_degs, 3)

vln_plot <- VlnPlot(integrated_data, features = top3_degs, group.by = "seurat_clusters", split.by = "sample", assay = "RNA") +
            ggtitle("Violin Plot of Example Shared DEGs") +
            theme(plot.title = element_text(hjust = 0.5))
ggsave("vlnplot_shared_degs.png", plot = vln_plot, width = 12, height = 8, dpi = 600, units = "in")

# b. exporting the shared DEGs
write.csv(shared_degs, "shared_degs_across_samples.csv", row.names = FALSE)

avg_exp_all <- AverageExpression(integrated_data, assays = "RNA", features = unique(shared_degs$gene), group.by = c("seurat_clusters", "sample"))$RNA
write.csv(avg_exp_all, "shared_degs_avg_expression.csv")

# save the shared GEX analysis
saveRDS(shared_degs, "shared_degs_gex_analysis.rds")




# 6. EDIT: Differential Expression Analysis for finding commonly expressed genes per sample (not cluster)
# # Load the integrated data object with clusters
# integrated_data <- readRDS("integrated_data_clusters.rds")
# # Switching to the RNA assay for DEG analysis
# DefaultAssay(integrated_data) <- "RNA"

# # Join layers in the RNA assay
# integrated_data <- JoinLayers(integrated_data, assay = "RNA")
# # Normalize a log-normalization on the integrated data for DEG analysis
# integrated_data <- NormalizeData(integrated_data, assay = "RNA")

# step a. Set sample as the identity for comparison
Idents(integrated_data) <- "sample"

# step b. Find DEGs across all samples by comparing each sample to the others
# We'll store DEGs for each sample vs. others in a list
sample_degs <- list()
# Get unique sample names
samples <- unique(integrated_data@meta.data$sample)

# Loop over each sample to find DEGs
for (s in samples) {
  degs <- FindMarkers(integrated_data, 
                      ident.1 = s,  # Sample of interest
                      ident.2 = NULL,  # Compare to all other samples
                      assay = "RNA",
                      min.pct = 0.25, 
                      logfc.threshold = 0.25, 
                      test.use = "wilcox")
  degs$gene <- rownames(degs)
  degs$sample <- s
  sample_degs[[s]] <- degs
}

# Combine all DEGs into a single data frame
all_sample_degs <- do.call(rbind, sample_degs)

# step c. Rank DEGs globally by adjusted p-value and log2FC
# We'll select DEGs based on the smallest adjusted p-value and largest |log2FC|
global_degs_df <- all_sample_degs %>%
                  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5) %>%
                  group_by(gene) %>%
                  summarise(mean_log2FC = mean(abs(avg_log2FC)), 
                            min_p_val_adj = min(p_val_adj))

print(summary(global_degs_df$min_p_val_adj))
print(table(global_degs_df$min_p_val_adj == 0))  # Check how many genes have p_val_adj = 0

# Select exactly 20 genes using slice_head()
global_degs_top20 <- global_degs_df %>%
               arrange(min_p_val_adj, desc(mean_log2FC)) %>%  # Sort by p-value, then by log2FC
               slice_head(n = 20) %>%  # Strictly take the top 20
               pull(gene) %>%
               unique()
print(global_degs_top20)

# Compute average expression of global_degs_top20 genes across samples
avg_exp_top20 <- AverageExpression(integrated_data, 
                             assays = "RNA", 
                             features = global_degs_top20, 
                             group.by = "sample")$RNA                             

# Select exactly 100 genes using slice_head()
global_degs_top100 <- global_degs_df %>%
               arrange(min_p_val_adj, desc(mean_log2FC)) %>%  # Sort by p-value, then by log2FC
               slice_head(n = 100) %>%  # Strictly take the top 100
               pull(gene) %>%
               unique()

# Compute average expression of global_degs_top100 across samples
avg_exp_top100 <- AverageExpression(integrated_data, 
                             assays = "RNA", 
                             features = global_degs_top100, 
                             group.by = "sample")$RNA
# Ensure the matrix is correctly formatted (genes as columns, samples as rows)
write.csv(avg_exp_top100, "shared_degs_top100_avg_expression.csv")
print(dim(avg_exp_top100))  # Should be 6 rows (samples) x 20 columns (genes)

# step d. Visualize the global DEGs
# Generate heatmap with a dark purple to light purple gradient
heatmap_plot <- pheatmap(avg_exp_top20, 
                         scale = "row", 
                         cluster_rows = FALSE, 
                         cluster_cols = FALSE, 
                         show_rownames = TRUE, 
                         show_colnames = TRUE, 
                         color = colorRampPalette(c("#00a2fa", "#53008e"))(50), 
                         main = "Heatmap of Top 20 Global DEGs Across Samples",
                         fontsize = 10,
                         border_color = "white")

# Save the heatmap as a PNG
png("6_heatmap_global_degs_top20.png", width = 14, height = 8, units = "in", res = 600)
print(heatmap_plot)
dev.off()

# Generate dot plot with a custom color palette
dot_plot <- DotPlot(integrated_data, 
                    features = global_degs_top20, 
                    group.by = "sample", 
                    assay = "RNA") +
            ggtitle("Dot Plot of Top 20 Global DEGs Across Samples") +
            theme(plot.title = element_text(hjust = 0.5),
                  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))  # Tilt x-axis labels by 45 degrees

# Save the dot plot as a PNG
ggsave("new_dotplot_global_degs_top20.png", plot = dot_plot, width = 14, height = 8, dpi = 600, units = "in")