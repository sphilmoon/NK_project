# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(dplyr)
library(ggplot2)
# ------------------------- #
# Define Output Directories
# ------------------------- #
output_dir <- "/home/outputs/totalNK_outputs/3_gex/wo_51_52"
pdf_dir <- file.path(output_dir, "pdf")
dge_dir <- file.path(output_dir, "dge")
rds_dir <- file.path(output_dir, "rds")
for (dir in list(pdf_dir, dge_dir, rds_dir)) {
 if (!dir.exists(dir)) {
   cat("📂 Creating directory:", dir, "\n")
   dir.create(dir, recursive = TRUE)
 }
}

# ------------------------- #
# Step 1: Load Merged Seurat Object
# ------------------------- #
merged_obj <- readRDS("/home/outputs/totalNK_outputs/2_umap/wo_51_52/rds/integrated_data_dims25_res0.3_genecounts.rds")

# table(merged_obj$orig.ident)
# table(merged_obj$sample_id)

# ------------------------- #
# Step 2: Split by Sample
# ------------------------- #
seurat_list <- SplitObject(merged_obj, split.by = "orig.ident")
# Step 3: SCTransform per Sample (retain full model)
seurat_list <- lapply(seurat_list, function(x) {
 SCTransform(x, verbose = FALSE, return.only.var.genes = FALSE)
})

# Step 4: Integration
features <- SelectIntegrationFeatures(seurat_list, nfeatures = 3000)
seurat_list <- PrepSCTIntegration(seurat_list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", anchor.features = features)
seurat_obj <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# ------------------------- #
# Step 5: Dimensionality Reduction and Clustering
# ------------------------- #
DefaultAssay(seurat_obj) <- "integrated"
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, npcs = 30)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:25)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:25)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)
# Store cluster identities
seurat_obj$seurat_clusters <- Idents(seurat_obj)



genes_to_check <- c("CD1D", "TCRbeta", "CD3E", "CD3D", "CD3G", "TRDC", "TRBC1", "TRBC2")
# Check presence
present_genes <- genes_to_check[genes_to_check %in% rownames(seurat_obj)]
# Report results
if (length(present_genes) > 0) {
 cat("✅ Found genes in dataset:", paste(present_genes, collapse = ", "), "\n")
} else {
 cat("❌ None of the specified genes were found in the dataset.\n")
}






# ------------------------- #
# Generate and Save UMAP Plot
# ------------------------- #
umap_plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) +
 ggtitle("UMAP of Integrated Data (dims=25, res=0.3)") +
 theme_minimal()
umap_pdf_file <- file.path(pdf_dir, "check_umap_integrated_data_dims25_res0.3_genecounts.pdf")
pdf(umap_pdf_file, width = 8, height = 6)
print(umap_plot)
dev.off()
cat("✅ UMAP plot saved to", umap_pdf_file, "\n")

# ------------------------- #
# Save Final Seurat Object
# ------------------------- #
rds_file_out <- file.path(rds_dir, "sct_integrated_data_dims25_res0.3_genecounts.rds")
saveRDS(seurat_obj, file = rds_file_out)
cat("✅ Integrated Seurat object saved to", rds_file_out, "\n")





# ------------------------- #
# Step 6: Prepare for DE
# ------------------------- #
seurat_obj <- PrepSCTFindMarkers(seurat_obj)
# ------------------------- #
# Step 7: Run DE Analysis
# ------------------------- #
markers <- FindAllMarkers(
 seurat_obj,
 only.pos = TRUE,
 min.pct = 0.25,
 logfc.threshold = 0.25
)
dge_file_out <- file.path(dge_dir, "markers_all_clusters_dims25_res0.3.csv")
write.csv(markers, dge_file_out, row.names = FALSE)
cat("✅ DEG results saved to", dge_file_out, "\n")