# ------------------------- #
# Load Required Libraries
# ------------------------- #
library(Seurat)
library(ggplot2)
library(cowplot)
library(hdf5r)

# ------------------------- #
# Define Directory Structure
# ------------------------- #
output_dir <- "/home/outputs/nkp46_outputs/chat"

pos_output_dir <- file.path(output_dir, "pos")
neg_output_dir <- file.path(output_dir, "neg")

pos_rds_output_dir <- file.path(pos_output_dir, "rds")
neg_rds_output_dir <- file.path(neg_output_dir, "rds")

pos_pdf_output_dir <- file.path(pos_output_dir, "pdf")
neg_pdf_output_dir <- file.path(neg_output_dir, "pdf")

dir.create(pos_rds_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(neg_rds_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pos_pdf_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(neg_pdf_output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------- #
# Define Input Directories
# ------------------------- #
nkp46_neg_dir <- "/home/rawdata/neg_nkp46"
nkp46_pos_dir <- "/home/rawdata/pos_nkp46"
sample_names <- c("animal25", "animal26", "animal27", "animal28")

# ------------------------- #
# Build File Paths (Validated)
# ------------------------- #
build_h5_paths <- function(base_dir, condition_prefix) {
    fnames <- paste0(tolower(sample_names), "_", condition_prefix, "_filtered_feature_bc_matrix.h5")
    full_paths <- file.path(base_dir, fnames)

    missing <- !file.exists(full_paths)
    if (any(missing)) {
        cat("❌ Missing files:\n")
        print(full_paths[missing])
        stop("Some .h5 files are missing. Please check the paths or filenames.")
    }

    return(full_paths)
}

# Build paths and validate existence
neg_file_paths <- build_h5_paths(nkp46_neg_dir, "ncr1neg")
pos_file_paths <- build_h5_paths(nkp46_pos_dir, "ncr1pos")

# ------------------------- #
# Function to Process Each Condition
# ------------------------- #
process_condition <- function(h5_files, sample_names, condition_label, rds_dir, pdf_dir) {
    message("📦 Processing condition: ", condition_label)

    # Load and annotate each sample
    seurat_objects <- mapply(function(file, name) {
        data <- Read10X_h5(file)
        obj <- CreateSeuratObject(counts = data, project = name, min.cells = 3, min.features = 200)
        obj$sample <- name
        obj$condition <- condition_label
        return(obj)
    }, h5_files, sample_names, SIMPLIFY = FALSE)

    names(seurat_objects) <- sample_names

    # Quality control filtering
    for (i in seq_along(seurat_objects)) {
        obj <- seurat_objects[[i]]
        obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 &
            nCount_RNA > 200 & nCount_RNA < 10000)
        seurat_objects[[i]] <- obj
    }

    # Generate Violin Plots
    plots <- lapply(seurat_objects, function(obj) {
        p1 <- VlnPlot(obj, features = "nFeature_RNA") + NoLegend() + ggtitle("nFeature_RNA")
        p2 <- VlnPlot(obj, features = "nCount_RNA") + NoLegend() + ggtitle("nCount_RNA")
        plot_grid(p1, p2, ncol = 2)
    })

    names(plots) <- sample_names
    qc_plot <- plot_grid(plotlist = plots, ncol = 1)
    ggsave(file.path(pdf_dir, paste0("nkp46_", condition_label, "_qc_violin.pdf")),
        plot = qc_plot, width = 10, height = 6 * length(plots), dpi = 600
    )

    # Normalize and integrate
    seurat_objects <- lapply(seurat_objects, SCTransform, verbose = FALSE)
    features <- SelectIntegrationFeatures(seurat_objects, nfeatures = 3000)
    anchors <- FindIntegrationAnchors(seurat_objects, anchor.features = features, dims = 1:25)
    integrated <- IntegrateData(anchors, dims = 1:25)

    DefaultAssay(integrated) <- "integrated"
    integrated <- ScaleData(integrated)
    integrated <- RunPCA(integrated, npcs = 25)
    integrated <- FindNeighbors(integrated, dims = 1:25)
    integrated <- FindClusters(integrated, resolution = 0.3)
    integrated <- RunUMAP(integrated, dims = 1:25)

    # Save integrated object
    saveRDS(integrated, file.path(rds_dir, paste0("nkp46_", condition_label, ".rds")))
    message("💾 Saved integrated object for ", condition_label)

    # UMAP Plots
    umap_by_sample <- DimPlot(integrated, group.by = "sample", pt.size = 0.3) +
        ggtitle(paste0("UMAP by Sample - ", condition_label)) +
        theme(plot.title = element_text(hjust = 0.5))

    ggsave(file.path(pdf_dir, paste0("nkp46_", condition_label, "_umap_by_sample.pdf")),
        plot = umap_by_sample, width = 10, height = 8, dpi = 600
    )

    umap_combined <- DimPlot(integrated, group.by = "seurat_clusters", pt.size = 0.3) +
        ggtitle(paste0("UMAP (Clusters) - ", condition_label)) +
        theme(plot.title = element_text(hjust = 0.5))

    ggsave(file.path(pdf_dir, paste0("nkp46_", condition_label, "_umap_combined.pdf")),
        plot = umap_combined, width = 10, height = 8, dpi = 600
    )
}

# ------------------------- #
# Process Both Conditions
# ------------------------- #
process_condition(neg_file_paths, sample_names, "neg", neg_rds_output_dir, neg_pdf_output_dir)
process_condition(pos_file_paths, sample_names, "pos", pos_rds_output_dir, pos_pdf_output_dir)
