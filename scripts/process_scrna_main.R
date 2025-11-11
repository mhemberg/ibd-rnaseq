################################################################################
# Main scRNA-seq Processing Pipeline
# IBD Transcriptomics Analysis
#
# Purpose: Complete preprocessing pipeline for IBD scRNA-seq data
# Input: Raw 10X data from multiple patients
# Output: Quality-controlled, batch-corrected Seurat objects by cell type
#
# Note: This is a computationally intensive pipeline
################################################################################

# Load configuration and utilities
source("scripts/config.R")
source("scripts/utils_functions.R")

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
})

# Script-specific parameters
OUTPUT_PREFIX <- "scrna_main"
SAMPLE_INFO_FILE <- file.path(RAW_DATA_DIR, "sample_info.txt")

message("==========================================")
message("Main scRNA-seq Processing Pipeline")
message("==========================================\n")

message("\nNOTE: This is a computationally intensive pipeline.")
message("Ensure you have:")
message("  - Sufficient RAM (>32GB recommended)")
message("  - Raw 10X data in: ", file.path(RAW_DATA_DIR, "cellranger_data"))
message("  - Sample metadata file: ", SAMPLE_INFO_FILE)

################################################################################
# Load Sample Metadata
################################################################################

message("\n--- Loading sample metadata ---")

if (!file.exists(SAMPLE_INFO_FILE)) {
  stop("Sample info file not found: ", SAMPLE_INFO_FILE)
}

sample_info <- read.table(SAMPLE_INFO_FILE, header = TRUE, row.names = 1)
message("Samples to process: ", nrow(sample_info))

################################################################################
# Load and Merge Data
################################################################################

message("\n--- Loading 10X data from all samples ---")

seurat_list <- list()

for (i in 1:nrow(sample_info)) {
  
  sample_id <- rownames(sample_info)[i]
  message("\nProcessing sample ", i, "/", nrow(sample_info), ": ", sample_id)
  
  # Construct path to 10X data
  #data_path <- file.path(RAW_DATA_DIR, "cellranger_data", sample_id, "filtered_feature_bc_matrix")
  data_path <- file.path(RAW_DATA_DIR, "cellranger_data", "geo_processed", sample_id)
  
  if (!dir.exists(data_path)) {
    warning("Data path not found: ", data_path)
    next
  }
  tryCatch({
  # Load 10X data
  counts <- Read10X(data.dir = data_path)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    min.cells = 0,
    min.features = 1
  )
  
  # Add metadata
  seurat_obj$Group <- sample_info$Group[i]
  seurat_obj$Sex <- sample_info$Sex[i]
  seurat_obj$Age <- sample_info$Age[i]
  seurat_obj$CD_location <- sample_info$CD_location[i]
  seurat_obj$pat <- strsplit(sample_id, "-")[[1]][2]
  seurat_obj$time <- strsplit(sample_id, "-")[[1]][3]
  
  # Calculate mitochondrial percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Apply QC filters
  seurat_obj <- subset(
    seurat_obj,
    subset = percent.mt <= QC_MT_THRESHOLD & nFeature_RNA >= QC_MIN_FEATURES
  )
  
  message("  Cells after QC: ", ncol(seurat_obj))
  
  seurat_list[[i]] <- seurat_obj
  }, error = function(e) {
  message("Error reading data for ", sample_id)
  })
}

# Remove NULL entries
seurat_list <- seurat_list[!sapply(seurat_list, is.null)]

message("\n--- Merging ", length(seurat_list), " samples ---")

# Merge all objects
total_seurat <- merge(seurat_list[[1]], y = seurat_list[-1])
message("Total cells after merging: ", ncol(total_seurat))

# Clean up
rm(seurat_list)
gc()

################################################################################
# Initial Processing
################################################################################

message("\n--- Initial normalization and clustering ---")

# Normalize
total_seurat <- NormalizeData(total_seurat)

# Find variable features per patient for consensus
message("Finding variable features per patient...")

pat_list <- unique(total_seurat@meta.data$pat)
pat_list <- pat_list[!is.na(pat_list)]  # Remove NAs if any
message("Number of unique patients: ", length(pat_list))
var_genes_list <- list()

for (i in seq_along(pat_list)) {
  current_pat <- pat_list[i]
  message("  Patient ", i, "/", length(pat_list), ": ", current_pat)
  
  # Subset to current patient using explicit metadata access
  pat_indices <- which(total_seurat@meta.data$pat == current_pat)
  message("    Cells for this patient: ", length(pat_indices))
  
  if (length(pat_indices) == 0) {
    warning("No cells found for patient ", current_pat)
    next
  }
  
  if (length(pat_indices) < 50) {
    warning("Too few cells for patient ", current_pat, ", skipping")
    next
  }
  
  pat_subset <- total_seurat[, pat_indices]
  message("    Finding variable features...")
  pat_subset <- FindVariableFeatures(pat_subset, nfeatures = 1500, verbose = FALSE)
  
  var_genes_list[[current_pat]] <- VariableFeatures(pat_subset)
  message("    Found ", length(var_genes_list[[current_pat]]), " variable features")
}

for (i in seq_along(pat_list)) {
  current_pat <- pat_list[i]
  message("  Patient ", i, "/", length(pat_list), ": ", current_pat)
  
  # Subset to current patient
  pat_indices <- which(total_seurat$pat == current_pat)
  pat_subset <- total_seurat[, pat_indices]
  
  pat_subset <- FindVariableFeatures(pat_subset, nfeatures = 1500)
  var_genes_list[[current_pat]] <- VariableFeatures(pat_subset)
}

# Get consensus variable genes (appear in multiple patients)
all_var_genes <- unlist(var_genes_list)
var_gene_counts <- table(all_var_genes)
consensus_genes <- names(var_gene_counts)[var_gene_counts >= 3]  # In at least 3 patients

message("Consensus variable genes: ", length(consensus_genes))

# Filter out ribosomal, mitochondrial, and IG genes
genes_to_remove <- c(
  grep("^MT-", consensus_genes, value = TRUE),
  grep("^RPS", consensus_genes, value = TRUE),
  grep("^RPL", consensus_genes, value = TRUE),
  grep("^IGH", consensus_genes, value = TRUE),
  grep("^IGK", consensus_genes, value = TRUE),
  grep("^IGL", consensus_genes, value = TRUE)
)

consensus_genes_filtered <- setdiff(consensus_genes, genes_to_remove)
message("After filtering: ", length(consensus_genes_filtered), " genes")

# Set variable features
VariableFeatures(total_seurat) <- consensus_genes_filtered

# Scale and run PCA
message("Scaling data...")
total_seurat <- ScaleData(total_seurat)

message("Running PCA...")
total_seurat <- RunPCA(total_seurat, verbose = FALSE, npcs = 50)

# Run UMAP
message("Running UMAP...")
total_seurat <- RunUMAP(total_seurat, dims = 1:30)

# Find neighbors and clusters
total_seurat <- FindNeighbors(total_seurat, dims = 1:30)
total_seurat <- FindClusters(total_seurat, resolution = 0.06)

# Save initial object
initial_file <- file.path(PROCESSED_DATA_DIR, "total_initial.rds")
saveRDS(total_seurat, initial_file)
message("Initial object saved: ", initial_file)

################################################################################
# Identify and Separate Major Cell Types
################################################################################

message("\n--- Identifying major cell types ---")

# Find all markers to identify cell types
all_markers <- FindAllMarkers(
  total_seurat,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Save markers
marker_file <- file.path(DEG_DIR, "total_initial_markers.txt")
save_table_standard(all_markers, marker_file)

# Separate by major cell types (based on marker genes)
# Epithelial: clusters with EPCAM, KRT8, KRT18
# T cells: clusters with CD3D, CD3E
# B cells: clusters with CD79A, MS4A1
# Stroma: clusters with COL1A2, VWF
# Myeloid: clusters with CD14, CD68

message("\nManual annotation required:")
message("  Review markers in: ", marker_file)
message("  Annotate clusters as: Epithelium, T_cell, BandPcell, Stroma, Myeloid")
message("  Save annotations to metadata$maincelltype")

# Apply MT filtering by cell type (example - adjust as needed)
message("\n--- Applying cell-type specific MT filtering ---")

# For this example, using cluster-based separation
# In practice, you would use marker-based annotation

epi_clusters <- c("0", "3", "8", "9")  # Example epithelial clusters
tcell_clusters <- c("1")               # Example T cell clusters
bcell_clusters <- c("2", "6")         # Example B cell clusters
stroma_clusters <- c("4", "7")        # Example stroma clusters
myeloid_clusters <- c("5", "10")      # Example myeloid clusters

# Subset by cell type
epi <- subset(total_seurat, subset = seurat_clusters %in% epi_clusters)
tcell <- subset(total_seurat, subset = seurat_clusters %in% tcell_clusters)
bcell <- subset(total_seurat, subset = seurat_clusters %in% bcell_clusters)
stroma <- subset(total_seurat, subset = seurat_clusters %in% stroma_clusters)
myeloid <- subset(total_seurat, subset = seurat_clusters %in% myeloid_clusters)

# Apply differential MT cutoffs
epi_filtered <- subset(epi, subset = percent.mt < 65)
tcell_filtered <- subset(tcell, subset = percent.mt < 25)
bcell_filtered <- subset(bcell, subset = percent.mt < 25)
stroma_filtered <- subset(stroma, subset = percent.mt < 25)
myeloid_filtered <- subset(myeloid, subset = percent.mt < 25)

message("Cells after MT filtering:")
message("  Epithelial: ", ncol(epi_filtered))
message("  T cells: ", ncol(tcell_filtered))
message("  B cells: ", ncol(bcell_filtered))
message("  Stroma: ", ncol(stroma_filtered))
message("  Myeloid: ", ncol(myeloid_filtered))

################################################################################
# Process Each Cell Type with Harmony
################################################################################

#' Process Cell Type with Harmony Batch Correction
#' 
#' @param seurat_obj Seurat object for one cell type
#' @param name Cell type name
#' @param harmony_dims Number of harmony dimensions
#' @return Processed Seurat object
process_celltype_harmony <- function(seurat_obj, name, harmony_dims = 30) {
  
  message("\n=== Processing ", name, " ===")
  
  # Normalize
  message("  Normalizing...")
  seurat_obj <- NormalizeData(seurat_obj)
  
  message("  Finding variable features per patient...")
  pat_list <- unique(seurat_obj@meta.data$pat)
  pat_list <- pat_list[!is.na(pat_list)]  # Remove NAs
  message("    Number of unique patients: ", length(pat_list))
  var_genes_list <- list()
  
  for (i in seq_along(pat_list)) {
    current_pat <- pat_list[i]
    message("    Patient ", i, "/", length(pat_list), ": ", current_pat)
    
    # Subset to current patient
    pat_indices <- which(seurat_obj@meta.data$pat == current_pat)
    
    if (length(pat_indices) == 0 || length(pat_indices) < 50) {
      message("      Skipping (insufficient cells)")
      next
    }
    
    message("      Cells: ", length(pat_indices))
    pat_subset <- seurat_obj[, pat_indices]
    pat_subset <- FindVariableFeatures(pat_subset, nfeatures = 1500, verbose = FALSE)
    
    var_genes_list[[current_pat]] <- VariableFeatures(pat_subset)
  }
  # Consensus variable genes
  all_var_genes <- unlist(var_genes_list)
  var_gene_counts <- table(all_var_genes)
  consensus_genes <- names(var_gene_counts)[var_gene_counts >= 3]
  
  # Filter unwanted genes
  genes_to_remove <- c(
    grep("^MT-", consensus_genes, value = TRUE),
    grep("^RPS", consensus_genes, value = TRUE),
    grep("^RPL", consensus_genes, value = TRUE),
    grep("^IGH", consensus_genes, value = TRUE),
    grep("^IGK", consensus_genes, value = TRUE),
    grep("^IGL", consensus_genes, value = TRUE),
    grep("^HLA", consensus_genes, value = TRUE)
  )
  
  consensus_genes_filtered <- setdiff(consensus_genes, genes_to_remove)
  message("  Variable genes: ", length(consensus_genes_filtered))
  
  # Set variable features
  VariableFeatures(seurat_obj) <- consensus_genes_filtered
  
  # Scale data
  message("  Scaling...")
  seurat_obj <- ScaleData(seurat_obj)
  
  # Run PCA
  message("  Running PCA...")
  seurat_obj <- RunPCA(seurat_obj, npcs = 100, verbose = FALSE)
  
  # Run Harmony for batch correction
  message("  Running Harmony...")
  seurat_obj <- RunHarmony(
    seurat_obj,
    group.by.vars = "pat",
    dims.use = 1:harmony_dims,
    plot_convergence = FALSE
  )
  
  # Run UMAP
  message("  Running UMAP...")
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:harmony_dims, reduction = "harmony")
  
  # Find neighbors and clusters at multiple resolutions
  message("  Clustering...")
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:harmony_dims)
  
  # Test multiple resolutions
  resolutions <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5)
  for (res in resolutions) {
    seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
  }
  
  return(seurat_obj)
}

# Process each cell type
message("\n--- Processing cell types with Harmony ---")

# Epithelial (50 dims)
epi_processed <- process_celltype_harmony(epi_filtered, "Epithelial", harmony_dims = 50)
save(epi_processed, file = file.path(PROCESSED_DATA_DIR, "epi2_badrm_harmony50"))

# T cells (35 dims)
tcell_processed <- process_celltype_harmony(tcell_filtered, "T_cell", harmony_dims = 35)
save(tcell_processed, file = file.path(PROCESSED_DATA_DIR, "tcell2_badrm_harmony35"))

# B cells (20 dims)
bcell_processed <- process_celltype_harmony(bcell_filtered, "B_cell", harmony_dims = 20)
save(bcell_processed, file = file.path(PROCESSED_DATA_DIR, "bcell2_badrm_harmony20"))

# Stroma (30 dims)
stroma_processed <- process_celltype_harmony(stroma_filtered, "Stroma", harmony_dims = 30)
save(stroma_processed, file = file.path(PROCESSED_DATA_DIR, "stroma2_badrm_harmony30"))

# Myeloid (25 dims)
myeloid_processed <- process_celltype_harmony(myeloid_filtered, "Myeloid", harmony_dims = 25)
save(myeloid_processed, file = file.path(PROCESSED_DATA_DIR, "myeloid2_badrm_harmony25"))

################################################################################
# Generate Visualizations
################################################################################

message("\n--- Generating visualizations for each cell type ---")

for (celltype_info in list(
  list(obj = epi_processed, name = "epithelial"),
  list(obj = tcell_processed, name = "tcell"),
  list(obj = bcell_processed, name = "bcell"),
  list(obj = stroma_processed, name = "stroma"),
  list(obj = myeloid_processed, name = "myeloid")
)) {
  
  obj <- celltype_info$obj
  name <- celltype_info$name
  
  message("\nGenerating plots for ", name)
  
  # UMAP by patient
  plot_file <- file.path(FIGURES_DIR, paste0(name, "_umap_by_patient.png"))
  png(plot_file, width = 800, height = 600)
  print(DimPlot(obj, group.by = "pat", reduction = "umap"))
  dev.off()
  
  # UMAP by group
  plot_file <- file.path(FIGURES_DIR, paste0(name, "_umap_by_group.png"))
  png(plot_file, width = 800, height = 600)
  print(DimPlot(obj, group.by = "Group", reduction = "umap"))
  dev.off()
  
  # UMAP split by group
  plot_file <- file.path(FIGURES_DIR, paste0(name, "_umap_split.png"))
  png(plot_file, width = 2000, height = 600)
  print(DimPlot(obj, group.by = "Group", split.by = "Group", reduction = "umap"))
  dev.off()
}

################################################################################
# Find Markers for Each Cell Type
################################################################################

message("\n--- Finding markers for each cell type ---")

# This is computationally intensive - may want to run separately
for (celltype_info in list(
  list(obj = epi_processed, name = "epithelial", res = "RNA_snn_res.1.3"),
  list(obj = tcell_processed, name = "tcell", res = "RNA_snn_res.1.3"),
  list(obj = bcell_processed, name = "bcell", res = "RNA_snn_res.0.7"),
  list(obj = stroma_processed, name = "stroma", res = "RNA_snn_res.1.5"),
  list(obj = myeloid_processed, name = "myeloid", res = "RNA_snn_res.1.3")
)) {
  
  obj <- celltype_info$obj
  name <- celltype_info$name
  resolution <- celltype_info$res
  
  message("\nFinding markers for ", name, " at resolution ", resolution)

# Set identity - extract as vector
if (resolution %in% colnames(obj@meta.data)) {
  Idents(obj) <- obj@meta.data[[resolution]]
} else {
  warning("Resolution ", resolution, " not found in metadata. Available: ", 
          paste(grep("RNA_snn_res", colnames(obj@meta.data), value = TRUE), collapse = ", "))
  next
}
  
  # Find all markers
  markers <- FindAllMarkers(
    obj,
    only.pos = FALSE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    verbose = FALSE
  )
  
  # Save markers
  marker_file <- file.path(DEG_DIR, paste0(name, "_markers.txt"))
  save_table_standard(markers, marker_file)
  
  message("  Markers saved: ", marker_file)
}

################################################################################
# Summary Report
################################################################################

message("\n==========================================")
message("Analysis Summary")
message("==========================================")
message("Total samples processed: ", length(pat_list))
message("\nProcessed cell types:")
message("  Epithelial: ", ncol(epi_processed), " cells")
message("  T cells: ", ncol(tcell_processed), " cells")
message("  B cells: ", ncol(bcell_processed), " cells")
message("  Stroma: ", ncol(stroma_processed), " cells")
message("  Myeloid: ", ncol(myeloid_processed), " cells")
message("\nOutput files:")
message("  Processed objects: ", PROCESSED_DATA_DIR)
message("  Markers: ", DEG_DIR)
message("  Figures: ", FIGURES_DIR)

# Save session info
save_session_info(file.path(RESULTS_DIR, paste0(OUTPUT_PREFIX, "_session_info.txt")))

