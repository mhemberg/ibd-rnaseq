################################################################################
# Process Validation Dataset (GSE233063)
# IBD Transcriptomics Analysis
#
# Purpose: Load and process cytokine treatment validation scRNA-seq data
# Input: GSE233063 raw data (10X format)
# Output: Merged, QC'd, and processed Seurat object
#
# Original file: other_scRNAseq/data.docx (date: 20240418)
################################################################################

# Load configuration and utilities
source("scripts/config.R")
source("scripts/utils_functions.R")

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# Script-specific parameters
OUTPUT_PREFIX <- "validation_gse233063"

################################################################################
# Main Analysis
################################################################################

message("==========================================")
message("Processing Validation Dataset (GSE233063)")
message("==========================================\n")

# Define treatment groups
message("Treatment groups: ", paste(TREATMENT_GROUPS, collapse = ", "))
message("True order: ", paste(TREATMENT_ORDER, collapse = ", "))

# Initialize list to store Seurat objects
seurat_list <- list()

# Loop through each treatment group
for (i in seq_along(TREATMENT_GROUPS)) {
  
  group_name <- TREATMENT_GROUPS[i]
  true_group <- TREATMENT_ORDER[i]
  
  message("\n--- Processing group: ", group_name, " ---")
  
  # Construct path to 10X data
  data_path <- file.path(VALIDATION_DATA_DIR, group_name)
  
  # Check if directory exists
  if (!dir.exists(data_path)) {
    warning("Data directory not found: ", data_path)
    warning("Skipping group: ", group_name)
    next
  }
  
  # Read 10X data
  message("Reading 10X data from: ", data_path)
  count_matrix <- Read10X(data.dir = data_path)
  
  # Create Seurat object
  message("Creating Seurat object...")
  seurat_obj <- CreateSeuratObject(counts = count_matrix)
  
  # Calculate QC metrics
  message("Calculating QC metrics...")
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^MT-"
  )
  
  # Apply QC filters
  n_cells_before <- ncol(seurat_obj)
  message("Cells before QC: ", n_cells_before)
  
  seurat_obj <- subset(
    seurat_obj,
    subset = percent.mt <= QC_MT_THRESHOLD & nFeature_RNA >= QC_MIN_FEATURES
  )
  
  n_cells_after <- ncol(seurat_obj)
  message("Cells after QC: ", n_cells_after)
  message("Cells removed: ", n_cells_before - n_cells_after)
  
  # Add group metadata
  seurat_obj$group <- true_group
  seurat_obj$treatment <- group_name
  
  # Store in list
  seurat_list[[i]] <- seurat_obj
}

# Check if any objects were created
if (length(seurat_list) == 0) {
  stop("No Seurat objects were created. Check data paths.")
}

# Remove control group for merging
message("\n--- Merging datasets ---")
message("Removing control group from merge list...")
seurat_list_no_ctrl <- seurat_list[-1]

# Merge all objects
message("Merging ", length(seurat_list_no_ctrl) + 1, " Seurat objects...")
merged_seurat <- merge(
  seurat_list[[1]],  # Control group as base
  y = seurat_list_no_ctrl
)

message("Total cells after merging: ", ncol(merged_seurat))

# Process merged object
message("\n--- Processing merged object ---")

# Set RNA assay as active
merged_seurat@active.assay <- "RNA"

# Normalize data
message("Normalizing data...")
merged_seurat <- NormalizeData(merged_seurat)

# Find variable features
message("Finding variable features...")
merged_seurat <- FindVariableFeatures(merged_seurat)

# Scale data
message("Scaling data...")
merged_seurat <- ScaleData(merged_seurat)

# Run PCA
message("Running PCA (", PCA_DIMS, " components)...")
merged_seurat <- RunPCA(
  merged_seurat,
  verbose = FALSE,
  npcs = PCA_DIMS
)

# Run UMAP
message("Running UMAP...")
merged_seurat <- RunUMAP(merged_seurat, dims = 1:PCA_DIMS)

# Find neighbors and clusters
message("Finding neighbors...")
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:PCA_DIMS)

message("Finding clusters (resolution = ", CLUSTER_RESOLUTION, ")...")
merged_seurat <- FindClusters(merged_seurat, resolution = CLUSTER_RESOLUTION)

# Save processed object
output_file <- file.path(PROCESSED_DATA_DIR, paste0(OUTPUT_PREFIX, "_processed.rds"))
message("\n--- Saving processed object ---")
message("Output file: ", output_file)
saveRDS(merged_seurat, file = output_file)

################################################################################
# Generate QC Plots
################################################################################

message("\n--- Generating visualizations ---")

# UMAP by group
plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_umap_by_group.png"))
message("Saving UMAP by group: ", plot_file)

png(plot_file, width = 800, height = 600, res = 150)
print(DimPlot(merged_seurat, group.by = "group"))
dev.off()

# UMAP split by group
plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_umap_split_by_group.png"))
message("Saving UMAP split by group: ", plot_file)

png(plot_file, width = 1200, height = 600, res = 150)
print(DimPlot(merged_seurat, group.by = "group", split.by = "group"))
dev.off()

# UMAP by cluster
plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_umap_by_cluster.png"))
message("Saving UMAP by cluster: ", plot_file)

png(plot_file, width = 800, height = 600, res = 150)
print(DimPlot(merged_seurat, pt.size = 0.1, label = TRUE, raster = FALSE))
dev.off()

################################################################################
# Find Cluster Markers
################################################################################

message("\n--- Finding cluster markers ---")
merged_seurat <- JoinLayers(merged_seurat)
deg_results <- FindAllMarkers(
  merged_seurat,
  only.pos = TRUE,
  logfc.threshold = DEG_LOGFC_THRESHOLD,
  min.pct = DEG_MIN_PCT
)

# Filter by adjusted p-value
deg_filtered <- deg_results[deg_results$p_val_adj < DEG_PADJ_CUTOFF, ]

# Save results
deg_output <- file.path(DEG_DIR, paste0(OUTPUT_PREFIX, "_cluster_markers.txt"))
message("Saving cluster markers: ", deg_output)
save_table_standard(deg_filtered, deg_output)

################################################################################
# Summary Statistics
################################################################################

message("\n--- Analysis Summary ---")
message("Total cells: ", ncol(merged_seurat))
message("Total genes: ", nrow(merged_seurat))
message("Number of clusters: ", length(unique(merged_seurat$seurat_clusters)))
message("Number of treatments: ", length(unique(merged_seurat$treatment)))

# Cells per treatment
cells_per_treatment <- table(merged_seurat$treatment)
message("\nCells per treatment:")
for (treatment in names(cells_per_treatment)) {
  message("  ", treatment, ": ", cells_per_treatment[treatment])
}

# Save session info
save_session_info(file.path(RESULTS_DIR, paste0(OUTPUT_PREFIX, "_session_info.txt")))

message("\n==========================================")
message("Validation data processing complete!")
message("==========================================")
