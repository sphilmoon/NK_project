# Script: umap_visualizations_all_conditions.R
# Purpose: Generate UMAP visualizations for NKp46-, NKp46+, Total NK, merged conditions, and grid by animal and condition

# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)

# ------------------------- #
# Define Directory Structure
# ------------------------- #
nkp_neg_rds <- "/home/outputs/nkp46_outputs/chat/neg/rds/nkp46_neg.rds"
nkp_pos_rds <- "/home/outputs/nkp46_outputs/chat/pos/rds/nkp46_pos.rds"
totalNK_rds <- "/home/outputs/totalNK_outputs/2_umap/wo_51_52/rds/integrated_data_dims25_res0.3_genecounts.rds"

output_dir <- "/home/outputs/all_merged_TotalNK_nkp46/chat"

all_merged_rds <- file.path(output_dir, "rds")
pdf_dir <- file.path(output_dir, "pdf")

# Create directories
dir.create(all_merged_rds, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------- #
# 1. Load Seurat Objects for Each Condition
# ------------------------- #
# Function to safely load and validate Seurat objects
load_seurat_object <- function(rds_path, condition_name) {
  obj <- tryCatch({
    readRDS(rds_path)
  }, error = function(e) {
    stop("❌ Failed to load Seurat object from ", rds_path, ": ", e$message)
  })
  
  # Validate that it's a Seurat object
  if (!inherits(obj, "Seurat")) {
    stop("❌ Loaded object from ", rds_path, " is not a valid Seurat object")
  }
  
  # Check number of cells
  n_cells <- ncol(obj)
  cat("📊 ", condition_name, " has ", n_cells, " cells\n")
  if (n_cells == 0) {
    stop("❌ ", condition_name, " has 0 cells. Cannot proceed.")
  }
  
  # Print metadata columns and sample data
  cat("📋 Metadata columns in ", condition_name, ":\n")
  print(colnames(obj@meta.data))
  cat("📄 First few rows of metadata in ", condition_name, ":\n")
  print(head(obj@meta.data))
  
  return(obj)
}

totalNK <- load_seurat_object(totalNK_rds, "totalNK")
nkp46_neg <- load_seurat_object(nkp_neg_rds, "nkp46_neg")
nkp46_pos <- load_seurat_object(nkp_pos_rds, "nkp46_pos")

# ------------------------- #
# 2. Assign Sample Identity (if not already present)
# ------------------------- #
# Function to safely assign metadata
assign_sample_id <- function(obj, condition_name, value) {
  n_cells <- ncol(obj)
  cat("🔍 Assigning 'sample_id' to ", condition_name, " with ", n_cells, " cells\n")
  metadata_vec <- rep(value, n_cells)
  cat("🔍 Metadata vector length:", length(metadata_vec), "\n")
  
  obj <- AddMetaData(obj, metadata = metadata_vec, col.name = "sample_id")
  cat("✅ Successfully assigned 'sample_id' to ", condition_name, "\n")
  return(obj)
}

# Assign sample_id with individual checks
totalNK <- assign_sample_id(totalNK, "totalNK", "TotalNK")
nkp46_neg <- assign_sample_id(nkp46_neg, "nkp46_neg", "NKp46neg")
nkp46_pos <- assign_sample_id(nkp46_pos, "nkp46_pos", "NKp46pos")

# Verify that 'sample' exists in each object
objects <- list(totalNK = totalNK, nkp46_neg = nkp46_neg, nkp46_pos = nkp46_pos)
for (name in names(objects)) {
  obj <- objects[[name]]
  cat("🔍 Checking 'sample' field in ", name, "\n")
  cat("📄 Metadata columns in ", name, " during check:\n")
  print(colnames(obj@meta.data))
  if (!"sample" %in% colnames(obj@meta.data)) {
    cat("⚠️ 'sample' field not found in ", name, "\n")
    tryCatch({
      if ("orig.ident" %in% colnames(obj@meta.data)) {
        cat("🔍 Found 'orig.ident' in ", name, ", inferring 'sample'\n")
        objects[[name]] <- AddMetaData(obj, metadata = obj$orig.ident, col.name = "sample")
        cat("✅ Inferred 'sample' from 'orig.ident' for ", name, "\n")
      } else {
        cat("❌ 'orig.ident' not found in ", name, ". Available columns:\n")
        print(colnames(obj@meta.data))
        stop("Cannot proceed without 'sample' field. Please check metadata.")
      }
    }, error = function(e) {
      cat("❌ Error inferring 'sample' from 'orig.ident' in ", name, ": ", e$message, "\n")
      cat("📄 Available metadata in ", name, ":\n")
      print(head(obj@meta.data))
      stop("Metadata inference failed. Check the object.")
    })
  } else {
    cat("✅ 'sample' field exists in ", name, "\n")
    cat("📄 First few 'sample' values in ", name, ":\n")
    print(head(obj$sample))
  }
}

# Update the objects after verification
totalNK <- objects$totalNK
nkp46_neg <- objects$nkp46_neg
nkp46_pos <- objects$nkp46_pos

# ------------------------- #
# 3. Merge All Three Conditions into One Integrated Object
# ------------------------- #

# Merge the objects
merged_obj <- merge(totalNK,
                    y = c(nkp46_pos, nkp46_neg),
                    add.cell.ids = c("TotalNK", "NKp46pos", "NKp46neg"))
cat("✅ Merged all three Seurat objects\n")

# Debug: Print metadata column names for merged object
cat("📋 Metadata columns in merged_obj:\n")
print(colnames(merged_obj@meta.data))

# Ensure 'sample' and 'sample_id' are preserved after merging
if (!"sample" %in% colnames(merged_obj@meta.data) || !"sample_id" %in% colnames(merged_obj@meta.data)) {
  cat("⚠️ 'sample' or 'sample_id' missing in merged object. Attempting to reassign...\n")
  # Reassign sample_id based on add.cell.ids
  merged_obj$sample_id <- factor(gsub("_.*", "", colnames(merged_obj)),
                                 levels = c("TotalNK", "NKp46pos", "NKp46neg"))
  # Attempt to infer 'sample' if possible
  if ("orig.ident" %in% colnames(merged_obj@meta.data)) {
    merged_obj$sample <- merged_obj$orig.ident
  } else {
    stop("Cannot infer 'sample' field in merged object. Please check metadata.")
  }
  cat("✅ Reassigned 'sample' and 'sample_id' in merged object\n")
}

# Re-normalize and re-scale the merged object
merged_obj <- NormalizeData(merged_obj, verbose = FALSE)
merged_obj <- FindVariableFeatures(merged_obj, selection.method = "vst", nfeatures = 2000)
merged_obj <- ScaleData(merged_obj, verbose = FALSE)
cat("✅ Re-normalized and scaled merged object\n")

# Run PCA
merged_obj <- RunPCA(merged_obj, npcs = 50, verbose = FALSE)
cat("✅ Ran PCA on merged object\n")

# Generate and save ElbowPlot
elbow_plot <- ElbowPlot(merged_obj, ndims = 50) +
              ggtitle("Elbow Plot for Merged Conditions") +
              theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(pdf_dir, "elbow_plot_merged.pdf"),
       plot = elbow_plot, width = 8, height = 6, dpi = 600)
cat("✅ Saved ElbowPlot to", file.path(pdf_dir, "elbow_plot_merged.pdf"), "\n")

# Proceed with dims = 1:25 as specified
dims_to_use <- 1:25
merged_obj <- FindNeighbors(merged_obj, dims = dims_to_use)
merged_obj <- FindClusters(merged_obj, resolution = 0.5)
merged_obj <- RunUMAP(merged_obj, dims = dims_to_use, verbose = FALSE)
cat("✅ Processed merged object with dims = 1:25 and resolution = 0.5\n")

# Save the merged object
saveRDS(merged_obj, file.path(all_merged_rds, "merged_all_conditions_dims25_res0.5.rds"))
cat("✅ Saved merged Seurat object to", file.path(all_merged_rds, "merged_all_conditions_dims25_res0.5.rds"), "\n")

# ------------------------- #
# 4. Generate Combined UMAP for All Conditions
# ------------------------- #
umap_all_combined <- DimPlot(merged_obj,
                             group.by = "sample_id", label = TRUE, label.size = 5,
                             repel = TRUE, pt.size = 0.3) +
                     ggtitle("UMAP - All Conditions Combined") +
                     theme(plot.title = element_text(hjust = 0.5), legend.position = "right") +
                     scale_color_brewer(palette = "Set2")
ggsave(file.path(pdf_dir, "umap_all_conditions_combined.pdf"),
       plot = umap_all_combined, width = 10, height = 8, dpi = 600)
cat("✅ Saved combined UMAP for all conditions to", file.path(pdf_dir, "umap_all_conditions_combined.pdf"), "\n")

# ------------------------- #
# 5. Generate UMAP Grid: Animals (Rows) x Conditions (Columns)
# ------------------------- #

# Define animals and conditions
animals <- c("Animal25", "Animal26", "Animal27", "Animal28")
conditions <- c("NKp46neg", "NKp46pos", "TotalNK")  # Match sample_id values

# Create a list to store UMAP plots
umap_grid_plots <- list()

# Generate UMAP for each animal-condition combination
for (animal in animals) {
  for (condition in conditions) {
    # Subset the merged object for the specific animal and condition
    subset_cells <- merged_obj$sample == animal & merged_obj$sample_id == condition
    if (sum(subset_cells) == 0) {
      cat("⚠️ No cells found for", animal, "and", condition, "- skipping\n")
      umap_grid_plots[[paste(animal, condition, sep = "_")]] <- ggplot() + theme_void() + ggtitle("No Data")
      next
    }
    
    subset_obj <- subset(merged_obj, subset = sample == animal & sample_id == condition)
    
    # Generate UMAP plot (no grouping since it's a single animal-condition pair)
    umap_plot <- DimPlot(subset_obj, group.by = "sample_id", pt.size = 0.3) +
                 ggtitle(paste(animal, condition, sep = " - ")) +
                 theme(plot.title = element_text(hjust = 0.5, size = 10),
                       axis.text = element_text(size = 8),
                       axis.title = element_text(size = 10),
                       legend.position = "none",
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
    
    umap_grid_plots[[paste(animal, condition, sep = "_")]] <- umap_plot
  }
}

# Arrange plots in a 4x3 grid (rows: animals, columns: conditions)
umap_grid <- wrap_plots(umap_grid_plots, nrow = length(animals), ncol = length(conditions))

# Save the grid as a PDF
ggsave(file.path(pdf_dir, "umap_grid_animals_conditions.pdf"),
       plot = umap_grid, width = 15, height = 12, dpi = 600)
cat("✅ Saved UMAP grid (animals x conditions) to", file.path(pdf_dir, "umap_grid_animals_conditions.pdf"), "\n")

# ------------------------- #
# Final Message
# ------------------------- #
cat("🎉 Processing complete. Outputs saved in:", output_dir, "\n")