# # Load libraries
# library(Seurat)
# library(SeuratDisk)
# library(dplyr)
# library(ggplot2)
# library(DESeq2)

# # Define file paths and sample names
# file_paths <- list(
#   "/home/rawdata/total_nk/animal25_totalnk_filtered_feature_bc_matrix.h5",
#   "/home/rawdata/total_nk/animal27_totalnk_filtered_feature_bc_matrix.h5",
#   "/home/rawdata/total_nk/animal51_totalnk_filtered_feature_bc_matrix.h5",
#   "/home/rawdata/total_nk/animal26_totalnk_filtered_feature_bc_matrix.h5",
#   "/home/rawdata/total_nk/animal28_totalnk_filtered_feature_bc_matrix.h5",
#   "/home/rawdata/total_nk/animal52_totalnk_filtered_feature_bc_matrix.h5"
# )
# sample_names <- c("Animal25", "Animal27", "Animal51", "Animal26", "Animal28", "Animal52")

# # step 1: Load individual H5 files into Seurat objects
# seurat_objects <- mapply(function(file, name) {
#   data <- Read10X_h5(file)
#   seurat_obj <- CreateSeuratObject(counts = data, project = name, min.cells = 3, min.features = 200)
#   seurat_obj$sample <- name  # Add sample metadata
#   return(seurat_obj)
# }, file_paths, sample_names, SIMPLIFY = FALSE)

# # Assign meaningful names to the list
# names(seurat_objects) <- sample_names

# # Print sample information with details
# for (name in names(seurat_objects)) {
#   cat("Sample:", name, "\n")
#   print(seurat_objects[[name]])
#   cat("Number of cells:", ncol(seurat_objects[[name]]), "\n")
#   cat("Number of genes:", nrow(seurat_objects[[name]]), "\n\n")
# }

# # Optional: Save the list for later use
# saveRDS(seurat_objects, "seurat_objects_list.rds")





# step 2: quality control and pre-processing

# Load libraries
library(Seurat)
library(SeuratDisk)  # For H5 file import if needed later
library(dplyr)
library(ggplot2)

# Load the saved RDS file
seurat_objects <- readRDS("seurat_objects_list.rds")

# Verify the loaded objects (optional)
cat("Number of samples loaded:", length(seurat_objects), "\n")
lapply(names(seurat_objects), function(name) {
  cat("Sample:", name, " - Cells:", ncol(seurat_objects[[name]]), " Genes:", nrow(seurat_objects[[name]]), "\n")
})

# Perform QC filtering and generate violin plots
for (i in 1:length(seurat_objects)) {
  # Calculate mitochondrial percentage
  seurat_objects[[i]][["percent.mt"]] <- PercentageFeatureSet(seurat_objects[[i]], pattern = "^MT-")
  
  # Filter based on QC thresholds
  seurat_objects[[i]] <- subset(seurat_objects[[i]], 
                                subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
  
  # Create violin plot
  p <- VlnPlot(seurat_objects[[i]], 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 3) +
       ggtitle(paste("QC for", names(seurat_objects)[i]))
  
  # Define output file name
  file_name <- paste0("QC_violin_", names(seurat_objects)[i], ".png")
  
  # Save plot as PNG with 600 DPI
  ggsave(file_name, plot = p, width = 10, height = 6, dpi = 600, units = "in")
  
  # Print status message
  cat("Figure generation done for", names(seurat_objects)[i], "- Saved as", file_name, "\n")
}

# Normalize with SCTransform
seurat_objects <- lapply(seurat_objects, SCTransform, verbose = FALSE)

# Save the QC'd and normalized objects
saveRDS(seurat_objects, "seurat_objects_qc_sct.rds")

# Final completion message
cat("QC and figure generation process is finished.\n")
