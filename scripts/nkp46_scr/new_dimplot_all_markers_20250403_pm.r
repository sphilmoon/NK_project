# Load libraries
library(Seurat)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(cowplot)

# Define directories
data_dir <- "/home/rawdata/nkp46_pos_neg/"
dge_output_dir <- "/home/outputs/nkp46_outputs"

# Create output directory if it doesn’t exist
dir.create(dge_output_dir, recursive = TRUE, showWarnings = FALSE)

# List of animals and conditions
animals <- c("25", "26", "27", "28", "52")
conditions <- c("nkp46+", "nkp46-")

# Initialize a list to store Seurat objects
seurat_list <- list()

# Load each .h5 file and create a Seurat object
for (condition in conditions) {
  for (animal in animals) {
    file_name <- paste0(data_dir, condition, "_animal", animal, ".h5")
    if (!file.exists(file_name)) {
      warning(paste("File not found:", file_name))
      next
    }
    data <- Read10X_h5(file_name)
    seurat_obj <- CreateSeuratObject(counts = data, 
                                     project = paste(condition, "animal", animal, sep = "_"),
                                     min.cells = 3, min.features = 200)
    seurat_obj$condition <- condition
    seurat_obj$animal <- paste0("Animal", animal)
    seurat_list[[paste(condition, animal, sep = "_")]] <- seurat_obj
  }
}
print(names(seurat_list))

# Preprocess each Seurat object
seurat_list <- lapply(seurat_list, function(seurat_obj) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj <- subset(seurat_obj, 
                       subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  return(seurat_obj)
})
lapply(seurat_list, function(x) ncol(x))

# Select features for integration
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 2000)

# Scale Data, Run PCA, and Generate Elbow Plots
combined_elbow_plots <- list()
for (i in 1:length(seurat_list)) {
  sample_name <- names(seurat_list)[i]
  seurat_list[[i]] <- ScaleData(seurat_list[[i]], features = features)
  seurat_list[[i]] <- RunPCA(seurat_list[[i]], features = features)
  elbow_plot <- ElbowPlot(seurat_list[[i]], ndims = 50) +
    ggtitle(paste("Elbow Plot for PCA -", sample_name)) +
    theme(plot.title = element_text(hjust = 0.5))
  combined_elbow_plots[[sample_name]] <- elbow_plot
  cat("✅ Elbow plot generated for", sample_name, "\n")
}
stacked_elbow_plot <- plot_grid(plotlist = combined_elbow_plots, ncol = 1)
output_file <- file.path(dge_output_dir, "combined_elbow_plots.pdf")
ggsave(output_file, plot = stacked_elbow_plot, 
       width = 10, height = 6 * length(seurat_list), dpi = 600, units = "in", 
       limitsize = FALSE)
cat("✅ Combined elbow plots saved to", output_file, "\n")

# Find Integration Anchors
dims_to_use <- 1:22
anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features, 
                                  reduction = "cca", dims = dims_to_use)

# Integrate the Data
integrated_data <- IntegrateData(anchorset = anchors, dims = dims_to_use)

# Dimensionality Reduction and Clustering
DefaultAssay(integrated_data) <- "integrated"
integrated_data <- ScaleData(integrated_data)
integrated_data <- RunPCA(integrated_data, npcs = 30)
integrated_data <- RunUMAP(integrated_data, dims = dims_to_use)
integrated_data <- FindNeighbors(integrated_data, dims = dims_to_use)
integrated_data <- FindClusters(integrated_data, resolution = 0.5)

# Save the integrated data object
saveRDS(integrated_data, file.path(dge_output_dir, "nkp46_integrated_data.rds"))

# Prepare for Dot Plots
DefaultAssay(integrated_data) <- "RNA"
integrated_data <- JoinLayers(integrated_data, assay = "RNA")

# Define animals and NK subset markers
animals <- paste0("Animal", c("25", "26", "27", "28", "52"))
chemotaxis_markers <- c("S1PR1", "CXCR3", "CXCR4", "S1PR5", "CX3CR1", "CXCR2", "S1PR4")
inhibiting_receptors_markers <- c("KLRC1", "CD300A", "TIGIT", "KIR3DL2", "KIR3DL1", "KIR2DL3", "KIR2DL1", "KLRB1", "SIGLEC7", "HAVCR2", "TIM-3")
activating_receptors_markers <- c("KLRK1", "KLRC2", "NCR1", "CD160", "NCR3", "SLAMF6", "SLAMF7")
cytokine_receptors_markers <- c("TGFBR1", "IL12RB1", "IL2RG", "TGFBR3", "TGFBR2", "IL10RA", "IL18R1", "IL18RAP", "IL2RB", "IL10RB")
activation_markers <- c("CD69", "TNFRSF18")
maturation_markers <- c("ITGAM", "ITGAX", "NCAM1", "FCGR3A")
nk1_markers <- c("PLAC8", "CD38", "CCL4", "CHST2", "FCER1G", "IGFBP7", "AKR1C3", "SPON2")
nk2_markers <- c("CD2", "GZMK", "XCL2", "CMC1", "SELL", "XCL1", "KLRC1", "IL7R", "TCF7", "CD44", "GPR183", "TPT1", "IL2RB", "COTL1", "CLDND1")
nk3_markers <- c("CD52", "S100A4", "LGALS1", "PTMS", "LINC01871", "VIM", "S100A6", "ITGB1", "IL32", "CCL5", "GZMH", "CD3E", "KLRC2")
adaptive_nk_markers <- c("KLRC2", "CD16", "FCER1G", "SH2D1B", "ZBTB16", "PRF1", "GZMA", "GZMB", "CD52", "NKG2C", "CD3", "CCL5", "NCAM1", "B3GAT1")
immature_nk_markers <- c("CD34", "CD7", "CD45RA", "CD244", "KIT", "IL2RB", "IL1R1", "NKG2D", "NCR1", "NCR3", "KLRB1", "GATA2", "EOMES", "CD16")
transitional_nk_markers <- c("CD7", "CD244", "KIT", "IL2RB", "IL1R1", "NKG2D", "NCR1", "NCR3", "NKG2A", "KLRF1", "KLRB1", "NCAM1", "TBX21", "GATA2", "EOMES")
mature_nk_markers <- c("CD7", "IL2RB", "CD244", "KIT", "IL1R1", "NKG2D", "NCR1", "NCR3", "NKG2A", "KLRF1", "KLRB1", "CD16", "KIR", "NCAM1", "TBX21", "CX3CR1")
terminally_mature_nk_markers <- c("CD7", "IL2RB", "CD244", "B3GAT1", "IL1R1", "NKG2D", "NCR1", "NCR3", "NKG2A", "KLRF1", "KLRB1", "CD16", "KIR", "NCAM1", "TMEM107", "LINC00910", "C12orf57", "RP11-386I14.4", "PTCH2", "NKTR", "DDX17", "CX3CR1", "MYBL1", "HAVCR2", "ZEB2")
active_nk_markers <- c("FOS", "FOSB", "JUN", "JUNB")
inhibiting_receptors_mice_markers <- c("Ly49A", "Ly49C", "Ly49I", "Ly49P", "NKR-P1B", "NKG2A", "CD244")
activating_receptors_mice_markers <- c("Ly49D", "NKG2D", "NKG2C", "Ly49H", "NKp46", "NKR-P1C", "NKR-P1G", "NKR-P1F", "CD226")
immature_nk_mice_markers <- c("CD27", "IL2RB", "CD244", "NKG2D", "KLRB1", "NCR1", "LY49", "ITGA2", "Runx3", "Id2", "Gata3", "Tox1/2", "Ets2", "IRF2", "Eomes", "SPN", "SELL", "CD226", "KLRB1", "ITGAM")
mature_nk_mice_markers <- c("IL2RB", "CD244", "NKG2D", "KLRB1", "NKG2A/C", "NCR1", "LY49", "ITGA2", "ITGAV", "ITAGX", "ITGAM", "SPN", "CD27", "Runx3", "Id2", "Gata3", "Tox1/2", "Ets2", "IRF2", "Eomes", "IKZF3", "TBX21", "SELL", "CD226", "KLRB1", "NKp46")
nk_mice_markers <- c("IL2RB", "CD244", "NKG2D", "KLRB1", "NKG2A/C", "NCR1", "LY49", "ITGA2", "ITGAV", "ITAGX", "ITGAM", "SPN", "CD27", "KLRG1", "IKZF3", "TBX21")
canonical_nk_markers <- c("NCR1", "CD2", "CD16", "GZMA", "PRF1", "KLRK1")
canonical_non_nk_markers <- c("CD40", "CD3a", "CD3b", "CD4", "CD8a", "CD14", "CD68", "A4GALT", "CD9", "CD1", "FCGR2A", "CD163L1")

# Check and define NKp46 condition if not present
if (!"condition" %in% colnames(integrated_data@meta.data)) {
  ncr1_exp <- GetAssayData(integrated_data, assay = "RNA", slot = "data")["NCR1", ]
  integrated_data$condition <- ifelse(ncr1_exp > median(ncr1_exp, na.rm = TRUE), "nkp46+", "nkp46-")
  cat("🔧 Defined NKp46+/- conditions based on NCR1 expression median.\n")
}

# Ensure only desired animals are included
integrated_data <- subset(integrated_data, subset = animal %in% animals)
cat("🌟 Subset data to include only Animals 25, 26, 27, 28, 52.\n")

# Generate Dot Plots (One Marker per Figure)
cat("🌟 Generating dot plots for NKp46+ and NKp46- across all animals...\n")
for (marker in nk1_markers) {  # Example: using nk1_markers; replace with desired set
  dot_plot <- DotPlot(integrated_data, 
                      features = marker, 
                      group.by = "condition",  # X-axis: NKp46+, NKp46-
                      split.by = "animal",     # Y-axis: Animal25, Animal26, etc.
                      dot.scale = 8) +
              ggtitle(paste("Expression of", marker, "in NKp46+ and NKp46- Across Animals")) +
              theme(plot.title = element_text(hjust = 0.5),
                    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                    axis.text.y = element_text(size = 10),
                    axis.title.x = element_text(size = 12),
                    axis.title.y = element_text(size = 12)) +
              scale_color_gradient(low = "lightblue", high = "darkblue") +  # Color gradient for expression
              xlab("Condition") +
              ylab("Animal")
  
  # Save the plot
  output_file <- file.path(dge_output_dir, paste0("dotplot_", marker, "_nkp46_all_animals.pdf"))
  ggsave(filename = output_file, plot = dot_plot, width = 8, height = 6, dpi = 600)
  cat(sprintf("📊 Dot plot for %s saved to %s \n", marker, output_file))
}