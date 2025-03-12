# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# set working directory
setwd("~/ephemeral/soumya_NK/rawdata/total_nk")

# step 1: pre-processing the individual samples
# creating a Seurat object for each sample
samples <- list()

# Loop through the six samples
for (i in 1:6) {
  # Load H5 file
  data <- Read10X_h5(paste0("sample", i, "_filtered_feature_bc_matrix.h5"))
  
  # Create Seurat object
  samples[[i]] <- CreateSeuratObject(counts = data, 
                                     project = paste0("Sample", i), 
                                     min.cells = 3, 
                                     min.features = 200)
  
  # Add sample identifier
  samples[[i]]$sample <- paste0("Sample", i)
  
  # Calculate mitochondrial percentage (adjust pattern if cattle MT genes differ)
  samples[[i]][["percent.mt"]] <- PercentageFeatureSet(samples[[i]], pattern = "^MT-")
  
  # Visualize QC metrics (optional, comment out if not needed)
  # VlnPlot(samples[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  # Filter cells based on QC thresholds (adjust as needed)
  samples[[i]] <- subset(samples[[i]], 
                         subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  
  # Normalize and find variable features
  samples[[i]] <- NormalizeData(samples[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  samples[[i]] <- FindVariableFeatures(samples[[i]], selection.method = "vst", nfeatures = 2000)
}