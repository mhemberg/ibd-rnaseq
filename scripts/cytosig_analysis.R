################################################################################
# CytoSig Cytokine Signature Analysis
# IBD Transcriptomics Analysis
#
# Purpose: Infer cytokine activity from gene expression using CytoSig
# Input: Seurat object with expression data
# Output: Cytokine activity scores per cell, comparisons across groups
#
# Reference: CytoSig (Jiang et al. 2021)
# Note: Requires CytoSig to be run externally or signature matrix
################################################################################

# Load configuration and utilities
source("scripts/config.R")
source("scripts/utils_functions.R")

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(viridis)
  library(data.table)
})

# Script-specific parameters
OUTPUT_PREFIX <- "cytosig"
CYTOSIG_DIR <- file.path(RAW_DATA_DIR, "cytosig")
CYTOSIG_RESULTS <- file.path(CYTOSIG_DIR, "ibd_total.Zscore")

################################################################################
# Custom Functions
################################################################################

#' Compare CytoSig Scores Across Groups
#'
#' @param seurat_obj Seurat object with cell metadata
#' @param compare_set Comparison set (e.g., c(6, 5))
#' @param cytosig_scores Matrix of cytosig scores (cytokines x cells)
#' @param annotation_level Annotation level ("group" or "intermediate")
#' @return Results data frame
cytosig_compare <- function(seurat_obj, compare_set, cytosig_scores, annotation_level) {
  
  # Determine which group variable to use
  if ((compare_set[1] %in% c(1, 2)) && (compare_set[2] %in% c(1, 2))) {
    group_order <- seurat_obj$group
  } else {
    group_order <- seurat_obj$new_group
  }
  
  set1 <- compare_set[1]
  set2 <- compare_set[2]
  
  message("\n=== Comparison: Group ", set1, " vs Group ", set2, " ===")
  
  # Get unique cell types
  if (annotation_level == "intermediate") {
    celltype_unique <- unique(seurat_obj$anotation_intermediate)
  } else if (annotation_level == "group") {
    celltype_unique <- unique(seurat_obj$anotation_group)
  } else {
    stop("annotation_level must be 'intermediate' or 'group'")
  }
  
  # Initialize results
  results <- data.frame()
  
  # Test each cell type
  for (celltype in celltype_unique) {
    
    message("  Testing cell type: ", celltype)
    
    # Get cell indices
    if (annotation_level == "intermediate") {
      set1_index <- which(seurat_obj$anotation_intermediate == celltype & group_order == set1)
      set2_index <- which(seurat_obj$anotation_intermediate == celltype & group_order == set2)
    } else {
      set1_index <- which(seurat_obj$anotation_group == celltype & group_order == set1)
      set2_index <- which(seurat_obj$anotation_group == celltype & group_order == set2)
    }
    
    # Skip if insufficient cells
    if (length(set1_index) == 0 || length(set2_index) == 0) {
      message("    Skipping - insufficient cells")
      next
    }
    
    # Get main cell type
    maincelltype <- as.character(seurat_obj$maincelltype[set1_index[1]])
    
    # Test each cytokine signature
    for (cytokine_idx in 1:nrow(cytosig_scores)) {
      
      cytokine <- rownames(cytosig_scores)[cytokine_idx]
      
      # Extract scores for each group
      group1_scores <- cytosig_scores[cytokine_idx, set1_index]
      group2_scores <- cytosig_scores[cytokine_idx, set2_index]
      
      # Calculate means
      mean1 <- mean(group1_scores)
      mean2 <- mean(group2_scores)
      
      # Wilcoxon test
      test_result <- tryCatch({
        wilcox.test(group1_scores, group2_scores, alternative = "two.sided")
      }, error = function(e) {
        list(p.value = NA)
      })
      
      pval <- test_result$p.value
      
      # Store results
      results <- rbind(results, data.frame(
        cytokine = cytokine,
        pval = pval,
        mean_group1 = mean1,
        mean_group2 = mean2,
        cell = celltype,
        maincelltype = maincelltype,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Adjust p-values
  results$qval <- p.adjust(results$pval, method = "BH")
  
  # Reorder columns
  results <- results[, c("cytokine", "pval", "qval", 
                        paste0("mean_group", set1), paste0("mean_group", set2),
                        "cell", "maincelltype")]
  
  # Rename mean columns
  colnames(results)[4:5] <- c(as.character(set1), as.character(set2))
  
  return(results)
}

################################################################################
# Main Analysis
################################################################################

message("==========================================")
message("CytoSig Cytokine Signature Analysis")
message("==========================================\n")

# Check if CytoSig results exist
if (!file.exists(CYTOSIG_RESULTS)) {
  message("CytoSig results not found: ", CYTOSIG_RESULTS)
  stop("Please generate CytoSig results first")
}

################################################################################
# Load Data
################################################################################

message("--- Loading CytoSig results ---")

# Load CytoSig scores
cytosig_data <- fread(CYTOSIG_RESULTS, check.names = FALSE, sep = "\t")
message("Loaded CytoSig scores: ", nrow(cytosig_data), " cytokines x ", 
        ncol(cytosig_data), " cells")

# Convert to matrix format
cytosig_data <- as.data.frame(cytosig_data)

# Remove duplicate columns
dup_cols <- duplicated(colnames(cytosig_data))
cytosig_data <- cytosig_data[, !dup_cols]

# Set row names (first column is cytokine names)
cytokine_names <- cytosig_data[, 1]
cytosig_matrix <- as.matrix(cytosig_data[, -1])
rownames(cytosig_matrix) <- cytokine_names

message("Unique cells: ", ncol(cytosig_matrix))
message("Cytokines analyzed: ", nrow(cytosig_matrix))

# Load Seurat object
message("\n--- Loading Seurat object ---")
total_file <- file.path(RAW_DATA_DIR, "total_no_umap.rdata")
validate_file_exists(total_file)
total3 <- get(load(total_file))

message("Total cells in Seurat object: ", ncol(total3))

# Verify cell order matches
if (ncol(cytosig_matrix) != ncol(total3)) {
  warning("Cell count mismatch between CytoSig and Seurat!")
  message("CytoSig cells: ", ncol(cytosig_matrix))
  message("Seurat cells: ", ncol(total3))
}

################################################################################
# Compare Across Patient Groups
################################################################################

message("\n--- Comparing cytokine activity across patient groups ---")

# Define comparisons
compare_sets <- list(c(6,5), c(3,4), c(4,5), c(6,3), c(6,7), c(7,5), 
                    c(1,2), c(5,2), c(5,1))

# Annotation levels to test
annotation_levels <- c("group", "intermediate")

for (annot_level in annotation_levels) {
  
  message("\n=== Annotation level: ", annot_level, " ===")
  
  for (compare_set in compare_sets) {
    
    set1 <- compare_set[1]
    set2 <- compare_set[2]
    
    # Perform comparison
    results <- tryCatch({
      cytosig_compare(total3, compare_set, cytosig_matrix, annot_level)
    }, error = function(e) {
      message("Error in comparison: ", e$message)
      return(NULL)
    })
    
    if (is.null(results) || nrow(results) == 0) {
      message("No results for comparison ", set1, " vs ", set2)
      next
    }
    
    # Save results
    output_file <- file.path(
      TABLES_DIR,
      paste0("signature_comp_", annot_level, "_", annot_level, "_", set1, "_", set2, ".txt")
    )
    save_table_standard(results, output_file)
    
    # Report significant results
    sig_results <- results[results$qval < 0.05, ]
    message("  Significant results (q < 0.05): ", nrow(sig_results))
  }
}

################################################################################
# Generate Dotplot for Key Cell Type
################################################################################

message("\n--- Generating dotplot visualization ---")

# Example: Undifferentiated epithelial cells, comparison 6 vs 5
comp_file <- file.path(
  TABLES_DIR,
  "signature_comp_group_group_6_5.txt"
)

if (file.exists(comp_file)) {
  
  message("Creating dotplot for Undifferentiated cells...")
  
  # Read with automatic row name detection
  first_lines <- readLines(comp_file, n = 2)
  header_cols <- length(strsplit(first_lines[1], "\t")[[1]])
  data_cols <- length(strsplit(first_lines[2], "\t")[[1]])
  
  if (data_cols > header_cols) {
    comp_data <- read.table(comp_file, check.names = FALSE, sep = "\t", 
                           row.names = 1, header = TRUE)
  } else {
    comp_data <- read.table(comp_file, check.names = FALSE, sep = "\t", header = TRUE)
  }
  
  # Ensure cytokine column exists
  if (!"cytokine" %in% colnames(comp_data)) {
    if (ncol(comp_data) >= 1) {
      colnames(comp_data)[1] <- "cytokine"
    }
  }
  
  # Filter for specific cell type
  if ("cell" %in% colnames(comp_data)) {
    cell_data <- comp_data[comp_data$cell == "Undifferentiated", ]
    
    if (nrow(cell_data) > 0) {
      
      # Order by q-value
      cyto_order <- cell_data$cytokine[order(cell_data$qval)]
      cell_data$cytokine <- factor(cell_data$cytokine, levels = rev(cyto_order))
      
      # Prepare data for plotting
      plot_data_6 <- data.frame(
        cytokine = cell_data$cytokine,
        qval = cell_data$qval,
        score = cell_data[, "6"],
        group = "PFDun",
        stringsAsFactors = FALSE
      )
      
      plot_data_5 <- data.frame(
        cytokine = cell_data$cytokine,
        qval = cell_data$qval,
        score = cell_data[, "5"],
        group = "CDun",
        stringsAsFactors = FALSE
      )
      
      plot_data <- rbind(plot_data_6, plot_data_5)
      plot_data$group <- factor(plot_data$group, levels = c("CDun", "PFDun"))
      plot_data$qval <- -log10(plot_data$qval)
      
      # Create dotplot
      plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_undiff_dotplot.pdf"))
      pdf(plot_file, width = 3, height = 9)
      
      p <- ggplot(plot_data, aes(x = group, y = cytokine, color = score, size = qval)) +
        geom_point() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        scale_color_viridis(discrete = FALSE, option = "inferno") +
        labs(
          x = "Group",
          y = "Cytokine",
          color = "Score",
          size = "-log10(q-value)"
        )
      
      print(p)
      dev.off()
      
      message("Dotplot saved: ", plot_file)
    }
  }
}

################################################################################
# Summary Statistics
################################################################################

message("\n--- Generating summary statistics ---")

# Initialize summary data frame
summary_stats <- data.frame(
  comparison = character(),
  annotation = character(),
  n_tests = integer(),
  n_significant = integer(),
  stringsAsFactors = FALSE
)

# Count significant results per comparison
for (annot_level in annotation_levels) {
  for (compare_set in compare_sets) {
    
    set1 <- compare_set[1]
    set2 <- compare_set[2]
    
    comp_file <- file.path(
      TABLES_DIR,
      paste0("signature_comp_", annot_level, "_", annot_level, "_", set1, "_", set2, ".txt")
    )
    
    if (file.exists(comp_file)) {
      
      # Read first few lines to detect format
      first_lines <- readLines(comp_file, n = 2)
      header_cols <- length(strsplit(first_lines[1], "\t")[[1]])
      data_cols <- length(strsplit(first_lines[2], "\t")[[1]])
      
      # Read appropriately based on column counts
      if (data_cols > header_cols) {
        # Data has row names
        comp_data <- read.table(comp_file, check.names = FALSE, header = TRUE, 
                               row.names = 1, sep = "\t")
      } else {
        # No row names
        comp_data <- read.table(comp_file, check.names = FALSE, header = TRUE, 
                               sep = "\t")
      }
      
      summary_stats <- rbind(summary_stats, data.frame(
        comparison = paste0(set1, "_vs_", set2),
        annotation = annot_level,
        n_tests = nrow(comp_data),
        n_significant = sum(comp_data$qval < 0.05, na.rm = TRUE),
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Save summary
summary_file <- file.path(TABLES_DIR, paste0(OUTPUT_PREFIX, "_summary.txt"))
save_table_standard(summary_stats, summary_file)

message("\n--- Summary Statistics ---")
print(summary_stats)

################################################################################
# Summary Report
################################################################################

message("\n==========================================")
message("Analysis Summary")
message("==========================================")
message("Cytokines analyzed: ", nrow(cytosig_matrix))
message("Cells analyzed: ", ncol(cytosig_matrix))
message("Annotation levels: ", paste(annotation_levels, collapse = ", "))
message("Comparisons performed: ", length(compare_sets))
message("\nTotal significant results (q < 0.05): ", sum(summary_stats$n_significant))
message("\nOutput files:")
message("  Comparison results: ", TABLES_DIR)
message("  Visualizations: ", FIGURES_DIR)

# Save session info
save_session_info(file.path(RESULTS_DIR, paste0(OUTPUT_PREFIX, "_session_info.txt")))

