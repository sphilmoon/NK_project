# Load libraries
library(Seurat)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(cowplot)

# Define directories
data_dir <- "/home/rawdata/nkp46_pos_neg/"
dge_output_dir <- "/home/outputs/nkp46_outputs"

# Create output directory if it doesn't exist
dir.create(dge_output_dir, recursive = TRUE, showWarnings = FALSE)

# Load your integrated Seurat object
integrated_data <- readRDS(file.path(dge_output_dir, "nkp46_integrated_data.rds"))

# Step 7: Prepare for Dot Plots
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

# Combine all markers into a single vector
all_markers <- c(
  chemotaxis_markers,
  inhibiting_receptors_markers,
  activating_receptors_markers,
  cytokine_receptors_markers,
  activation_markers,
  maturation_markers,
  nk1_markers,
  nk2_markers,
  nk3_markers,
  adaptive_nk_markers,
  immature_nk_markers,
  transitional_nk_markers,
  mature_nk_markers,
  terminally_mature_nk_markers,
  active_nk_markers,
  inhibiting_receptors_mice_markers,
  activating_receptors_mice_markers,
  immature_nk_mice_markers,
  mature_nk_mice_markers,
  nk_mice_markers,
  canonical_nk_markers,
  canonical_non_nk_markers
)

# Step 8: Generate Custom Dot Plots (One Marker per Figure)
cat("🌟 Generating custom dot plots for NKp46+ and NKp46- across all animals...\n")
for (marker in all_markers) {
  # Extract expression data
  expr_data <- GetAssayData(integrated_data, assay = "RNA", slot = "data")[marker, , drop = FALSE]
  
  # Check if the marker exists in the data
  if (nrow(expr_data) == 0) {
    cat(sprintf("⚠️ Marker %s not found in the data, skipping...\n", marker))
    next
  }
  
  # Create a data frame with metadata
  plot_data <- data.frame(
    animal = integrated_data$animal,
    condition = integrated_data$condition,
    expression = as.vector(expr_data[marker, ])
  )
  
  # Calculate average expression and percent expressed per group
  summary_data <- plot_data %>%
    group_by(animal, condition) %>%
    summarise(
      avg_expression = mean(expression, na.rm = TRUE),
      percent_expressed = mean(expression > 0, na.rm = TRUE) * 100,
      .groups = "drop"
    )
  
  # Create the dot plot using ggplot2
  dot_plot <- ggplot(summary_data, aes(x = condition, y = animal, size = percent_expressed, color = avg_expression)) +
    geom_point() +
    scale_size_continuous(range = c(2, 8), breaks = c(20, 40, 60)) +
    scale_color_gradient(low = "lightgrey", high = "red") +
    ggtitle(paste("Expression of", marker, "in NKp46+ and NKp46- Across Animals")) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12, margin = margin(t = 10)),
      axis.title.y = element_text(size = 12, margin = margin(r = 10)),
      legend.position = "right"
    ) +
    xlab("Condition") +
    ylab("Animal")
  
  # Save the plot
  output_file <- file.path(dge_output_dir, paste0("dotplot_", marker, "_nkp46_all_animals.pdf"))
  ggsave(filename = output_file, plot = dot_plot, width = 8, height = 10, dpi = 600)
  cat(sprintf("📊 Dot plot for %s saved to %s \n", marker, output_file))
}