################################################################################
# Project Configuration File
# IBD Transcriptomics Analysis
# 
# Purpose: Centralized configuration for all analysis scripts
# Usage: source("config.R") at the beginning of each script
################################################################################

# Project Metadata ----
PROJECT_NAME <- "ibd-transcriptomics"
AUTHOR <- "Your Name"
DATE_CREATED <- "2025-10-22"

# Directory Structure ----
# Root directory - modify this to match your system
ROOT_DIR <- "/Users/mh1015/JaewonBarcelonaFiles/ClaudeScripts/githubtest/ibd-transcriptomics"  # CHANGE THIS TO YOUR PROJECT ROOT

# Data directories
DATA_DIR <- file.path(ROOT_DIR, "data")
RAW_DATA_DIR <- file.path(DATA_DIR, "raw")
PROCESSED_DATA_DIR <- file.path(DATA_DIR, "processed")
EXTERNAL_DATA_DIR <- file.path(DATA_DIR, "external")

# Results directories
RESULTS_DIR <- file.path(ROOT_DIR, "results")
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
TABLES_DIR <- file.path(RESULTS_DIR, "tables")
DEG_DIR <- file.path(RESULTS_DIR, "deg_results")

# Reference data directories
REF_DIR <- file.path(ROOT_DIR, "references")
GENOME_DIR <- file.path(REF_DIR, "genome")
ANNOTATION_DIR <- file.path(REF_DIR, "annotations")

# Create directories if they don't exist
dirs_to_create <- c(
  DATA_DIR, RAW_DATA_DIR, PROCESSED_DATA_DIR, EXTERNAL_DATA_DIR,
  RESULTS_DIR, FIGURES_DIR, TABLES_DIR, DEG_DIR,
  REF_DIR, GENOME_DIR, ANNOTATION_DIR
)

for (dir in dirs_to_create) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    message("Created directory: ", dir)
  }
}

# External Data Paths ----
# scRNA-seq data
SCRNA_DATA_PATH <- file.path(RAW_DATA_DIR, "tcell2_badrm_harmony35")
STROMA_DATA_PATH <- file.path(RAW_DATA_DIR, "stroma2_badrm_harmony30")
IBD_DEG_DATA_PATH <- file.path(EXTERNAL_DATA_DIR, "DE_data_new1.RDS")

# Other scRNA-seq validation data
VALIDATION_DATA_DIR <- file.path(RAW_DATA_DIR, "GSE233063")

# Bulk RNA-seq paths
BULK_STAR_DIR <- file.path(RAW_DATA_DIR, "bulk_rna", "star")
BULK_ANNOTATION <- file.path(GENOME_DIR, "gencode.v41.annotation.gtf")

# Microarray data
IL13_ARRAY_PATH <- file.path(RAW_DATA_DIR, "GSE190705_read_counts.txt")

# Reference Genome Paths ----
GENOME_FASTA <- file.path(GENOME_DIR, "GRCh38.p13.genome.fa")
GENOME_GTF <- file.path(GENOME_DIR, "gencode.v41.annotation.gtf")
STAR_INDEX_DIR <- file.path(GENOME_DIR, "STAR-2.7.10a")

# Analysis Parameters ----
# Quality control thresholds
QC_MT_THRESHOLD <- 25        # Maximum mitochondrial percentage
QC_MIN_FEATURES <- 200       # Minimum features per cell
QC_MIN_CELLS <- 3            # Minimum cells per feature

# Clustering parameters
PCA_DIMS <- 30               # Number of PCA dimensions to use
CLUSTER_RESOLUTION <- 0.1    # Clustering resolution

# Differential expression parameters
DEG_LOGFC_THRESHOLD <- 0.25  # Log2 fold change threshold
DEG_MIN_PCT <- 0.25          # Minimum percentage of cells expressing
DEG_PADJ_CUTOFF <- 0.05      # Adjusted p-value cutoff
DEG_FC_THRESHOLDS <- c(1, 1.4, 2)  # Fold change thresholds for analysis

# Fisher exact test parameters
FISHER_TOTAL_GENES <- 33538  # Total number of genes for Fisher test

# Sample Groups ----
# Treatment groups for validation dataset
TREATMENT_GROUPS <- c("ctrl", "il1b", "il36", "osm", "tgfb", "tl1a", "tnfa")
TREATMENT_ORDER <- c("ctrl", "il1b", "osm", "tgfb", "tl1a", "il36", "tnfa")

# Patient groups (coding: 1=CDinfl, 2=CDuninfl, 3=CDinfl, 4=PFDinfl, 5=CDuninfl, 6=PFDuninfl, 7=HC)
PATIENT_GROUPS <- list(
  CD_UNINFL = 5,
  PFD_UNINFL = 6,
  CD_INFL = 3,
  PFD_INFL = 4,
  HC = 7,
  CD_UNINFL_ALT = 2,
  CD_INFL_ALT = 1
)

# Comparison sets for analysis
COMPARISON_SETS <- list(
  c(6, 5),   # PFD uninflamed vs CD uninflamed
  c(3, 4),   # CD inflamed vs PFD inflamed
  c(4, 5),   # PFD inflamed vs CD uninflamed
  c(6, 3),   # PFD uninflamed vs CD inflamed
  c(6, 7),   # PFD uninflamed vs Healthy Control
  c(7, 5),   # Healthy Control vs CD uninflamed
  c(1, 2),   # CD inflamed alt vs CD uninflamed alt
  c(5, 2),   # CD uninflamed vs CD uninflamed alt
  c(5, 1)    # CD uninflamed vs CD inflamed alt
)

# Cell Types ----
# IL-22 expressing cell types
IL22_CELLTYPES <- c(
  "CD4_effector_memory",
  "CD4_effector_TNF",
  "Th17",
  "ILC3"
)

# Bulk RNA-seq sample groups
BULK_GROUPS <- c("ltb", "tgfb", "tnf")

# Plot Parameters ----
PLOT_WIDTH <- 7
PLOT_HEIGHT <- 5
PLOT_DPI <- 300
CAIRO_TYPE <- "cairo"  # For better rendering

# Color schemes
COLORS_TREATMENT <- c(
  ctrl = "#808080",
  il1b = "#E41A1C",
  il36 = "#377EB8",
  osm = "#4DAF4A",
  tgfb = "#984EA3",
  tl1a = "#FF7F00",
  tnfa = "#FFFF33"
)

# Required Libraries ----
REQUIRED_PACKAGES <- c(
  "Seurat",
  "DESeq2",
  "edgeR",
  "biomaRt",
  "Rsubread",
  "ggplot2",
  "dplyr",
  "tidyr",
  "pheatmap"
)

# Function to check and install required packages
check_packages <- function() {
  missing_pkgs <- REQUIRED_PACKAGES[!(REQUIRED_PACKAGES %in% installed.packages()[,"Package"])]
  
  if (length(missing_pkgs) > 0) {
    message("The following packages are missing: ", paste(missing_pkgs, collapse = ", "))
    message("Please install them using BiocManager::install() or install.packages()")
    return(FALSE)
  }
  
  message("All required packages are installed.")
  return(TRUE)
}

# Session Information ----
save_session_info <- function(output_file = NULL) {
  if (is.null(output_file)) {
    output_file <- file.path(RESULTS_DIR, paste0("session_info_", Sys.Date(), ".txt"))
  }
  
  sink(output_file)
  cat("Analysis Date:", as.character(Sys.Date()), "\n\n")
  cat("R Version Information:\n")
  print(sessionInfo())
  sink()
  
  message("Session info saved to: ", output_file)
}

# Utility Functions ----
# Create timestamped output filename
create_output_name <- function(base_name, extension = ".txt", timestamp = TRUE) {
  if (timestamp) {
    timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
    return(paste0(base_name, "_", timestamp_str, extension))
  } else {
    return(paste0(base_name, extension))
  }
}

# Print configuration summary
print_config <- function() {
  cat("================================\n")
  cat("IBD Transcriptomics Analysis Configuration\n")
  cat("================================\n")
  cat("Project:", PROJECT_NAME, "\n")
  cat("Root Directory:", ROOT_DIR, "\n")
  cat("Data Directory:", DATA_DIR, "\n")
  cat("Results Directory:", RESULTS_DIR, "\n")
  cat("PCA Dimensions:", PCA_DIMS, "\n")
  cat("DEG Log FC Threshold:", DEG_LOGFC_THRESHOLD, "\n")
  cat("================================\n")
}

# Initialization message
message("Configuration loaded successfully!")
message("Root directory: ", ROOT_DIR)
message("Run print_config() to see full configuration")
