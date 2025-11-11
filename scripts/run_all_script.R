################################################################################
# Complete Analysis Pipeline
# IBD Transcriptomics Analysis
#
# Purpose: Master script to run the complete analysis pipeline
# Usage: source("run_all_scripts.R")
#
# Note: This script runs all analyses in sequence. Individual scripts can
#       also be run independently if only specific analyses are needed.
################################################################################

# Set options
options(width = 120)
options(warn = 1)  # Print warnings as they occur

# Record start time
start_time <- Sys.time()

cat("\n")
cat("================================================================================\n")
cat("  IBD TRANSCRIPTOMICS COMPLETE ANALYSIS PIPELINE\n")
cat("================================================================================\n")
cat("Start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================================\n\n")

# Load configuration
source("scripts/config.R")

# Print configuration
print_config()

# Check packages
if (!check_packages()) {
  stop("Required packages are missing. Please install them and try again.")
}

cat("\nSetup complete!\n")


# 1.1: Validation dataset (GSE233063)
cat("\n[1.1] Processing validation dataset...\n")
tryCatch({
  source("scripts/process_validation.R")
  cat("[1.1] ✓ Validation data processing complete\n")
}, error = function(e) {
  cat("[1.1] ✗ Error:", conditionMessage(e), "\n")
  cat("[1.1] Continuing with remaining analyses...\n")
})

# 1.2: IL-13 microarray
cat("\n[1.2] Processing IL-13 microarray data...\n")
tryCatch({
  source("scripts/il13_preprocess.R")
  cat("[1.2] ✓ IL-13 microarray processing complete\n")
}, error = function(e) {
  cat("[1.2] ✗ Error:", conditionMessage(e), "\n")
  cat("[1.2] Continuing with remaining analyses...\n")
})

# 1.3: Bulk RNA-seq (optional - requires aligned data)
cat("\n[1.3] Bulk RNA-seq processing (optional)...\n")
if (file.exists("scripts/process_bulk_rnaseq.R")) {
  cat("[1.3] Skipped - run separately if you have aligned bulk RNA-seq data\n")
} else {
  cat("[1.3] Script not yet created\n")
}


# 2.1: IL-22 T cell analysis
cat("\n[2.1] IL-22 T cell differential expression...\n")
tryCatch({
  source("scripts/deg_il22_script.R")
  cat("[2.1] ✓ IL-22 T cell analysis complete\n")
}, error = function(e) {
  cat("[2.1] ✗ Error:", conditionMessage(e), "\n")
  cat("[2.1] Continuing with remaining analyses...\n")
})

# 2.2: TL1A and IL22RA1 epithelial analysis
cat("\n[2.2] TL1A and IL22RA1 epithelial analysis...\n")
tryCatch({
  source("scripts/tl1a_il22ra1_analysis.R")
  cat("[2.2] ✓ TL1A/IL22RA1 epithelial analysis complete\n")
}, error = function(e) {
  cat("[2.2] ✗ Error:", conditionMessage(e), "\n")
  cat("[2.2] Continuing with remaining analyses...\n")
})

# 3.1: Fisher enrichment
cat("\n[3.1] Fisher exact test enrichment...\n")
tryCatch({
  source("scripts/fisher_enrichment.R")
  cat("[3.1] ✓ Fisher enrichment analysis complete\n")
}, error = function(e) {
  cat("[3.1] ✗ Error:", conditionMessage(e), "\n")
  cat("[3.1] Continuing with remaining analyses...\n")
})

# 3.2: IL-13 module scoring
cat("\n[3.2] IL-13 module scoring...\n")
tryCatch({
  source("scripts/il13_module_score.R")
  cat("[3.2] ✓ IL-13 module scoring complete\n")
}, error = function(e) {
  cat("[3.2] ✗ Error:", conditionMessage(e), "\n")
  cat("[3.2] Continuing with remaining analyses...\n")
})


# 4.1: IL-13 signature heatmaps
cat("\n[4.1] IL-13 signature heatmaps...\n")
tryCatch({
  source("scripts/il13_heatmap.R")
  cat("[4.1] ✓ IL-13 heatmap generation complete\n")
}, error = function(e) {
  cat("[4.1] ✗ Error:", conditionMessage(e), "\n")
  cat("[4.1] Continuing with remaining analyses...\n")
})

################################################################################
# Final Summary
################################################################################

end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "mins")

cat("\n")
cat("================================================================================\n")
cat("  ANALYSIS PIPELINE COMPLETE\n")
cat("================================================================================\n")
cat("Start time:  ", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("End time:    ", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Elapsed time:", round(elapsed_time, 2), "minutes\n")
cat("================================================================================\n")

# Count output files
n_figures <- length(list.files(FIGURES_DIR, pattern = "\\.(pdf|png)$"))
n_tables <- length(list.files(TABLES_DIR, pattern = "\\.txt$"))
n_deg <- length(list.files(DEG_DIR, pattern = "\\.txt$"))

cat("\nOutput Summary:\n")
cat("  Figures:", n_figures, "files in", FIGURES_DIR, "\n")
cat("  Tables:", n_tables, "files in", TABLES_DIR, "\n")
cat("  DEG results:", n_deg, "files in", DEG_DIR, "\n")


# Save pipeline completion info
pipeline_log <- file.path(RESULTS_DIR, "pipeline_completion.txt")
sink(pipeline_log)
cat("IBD Transcriptomics Analysis Pipeline\n")
cat("======================================\n\n")
cat("Completion Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total Runtime:", round(elapsed_time, 2), "minutes\n\n")
cat("System Information:\n")
print(sessionInfo())
sink()

cat("Pipeline log saved to:", pipeline_log, "\n\n")
