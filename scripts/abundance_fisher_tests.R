################################################################################
# Cell Abundance Fisher Exact Tests
# IBD Transcriptomics Analysis
#
# Purpose: Test for differential cell type abundance between conditions
# Input: Seurat object with cell type annotations
# Output: Fisher test results, fold changes, visualizations
#
# Comparisons:
# - PFD vs CD (overall)
# - PFD uninflamed vs CD uninflamed (6 vs 5)
# - CD inflamed vs PFD inflamed (3 vs 4)
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
OUTPUT_PREFIX <- "abundance_fisher"

################################################################################
# Custom Functions
################################################################################

#' Perform Fisher Exact Test for Cell Type Abundance
#'
#' @param seurat_obj Seurat object
#' @param name Annotation level ("maincelltype", "anotation_group", etc.)
#' @param group_var Group variable name (e.g., "pfd", "new_group")
#' @return Data frame with Fisher test results
fisher_abundance <- function(seurat_obj, name, group_var = "pfd") {
  
  message("\n=== Testing abundance: ", name, " ===")
  
  # Extract metadata as data frame
  metadata <- seurat_obj@meta.data
  
  # Debug information
  message("Metadata dimensions: ", nrow(metadata), " x ", ncol(metadata))
  message("Group variable: ", group_var)
  message("Group variable in colnames: ", group_var %in% colnames(metadata))
  
  # Extract vectors using pull-like operation
  group_vec <- metadata[[group_var]]
  
  # Check what we got
  message("Extracted group_vec class: ", paste(class(group_vec), collapse = ", "))
  message("Is data.frame: ", is.data.frame(group_vec))
  message("Is list: ", is.list(group_vec))
  
  # Convert to simple vector
  if (is.data.frame(group_vec)) {
    group_vec <- group_vec[[1]]
  }
  if (is.list(group_vec) && !is.data.frame(group_vec)) {
    group_vec <- unlist(group_vec, use.names = FALSE)
  }
  if (is.factor(group_vec)) {
    group_vec <- as.character(group_vec)
  }
  
  message("After conversion - class: ", class(group_vec), ", length: ", length(group_vec))
  
  # Extract cell type annotation
  if (name == "maincelltype") {
    celltype_vec <- metadata[["maincelltype"]]
  } else if (name == "anotation_group") {
    celltype_vec <- metadata[["anotation_group"]]
  } else if (name == "anotation_intermediate") {
    celltype_vec <- metadata[["anotation_intermediate"]]
  } else if (name == "anotation_refined") {
    celltype_vec <- metadata[["anotation_refined"]]
  } else {
    stop("Unknown annotation level: ", name)
  }
  
  # Convert to simple vector
  if (is.data.frame(celltype_vec)) {
    celltype_vec <- celltype_vec[[1]]
  }
  if (is.list(celltype_vec) && !is.data.frame(celltype_vec)) {
    celltype_vec <- unlist(celltype_vec, use.names = FALSE)
  }
  if (is.factor(celltype_vec)) {
    celltype_vec <- as.character(celltype_vec)
  }
  
  message("Cell type vector class: ", class(celltype_vec), ", length: ", length(celltype_vec))
  
  # Verify lengths match
  if (length(group_vec) != length(celltype_vec)) {
    stop("Length mismatch: group_vec = ", length(group_vec), ", celltype_vec = ", length(celltype_vec))
  }
  
  # Create contingency table
  message("Creating contingency table...")
  contingency <- table(group_vec, celltype_vec)
  message("Contingency table dimensions: ", paste(dim(contingency), collapse = " x "))
  
  # Calculate group totals
  group_totals <- rowSums(contingency)
  cell_types <- colnames(contingency)
  
  # Initialize results
  fc_results <- numeric(length(cell_types))
  pval_results <- numeric(length(cell_types))
  
  # Test each cell type
  for (i in seq_along(cell_types)) {
    
    celltype <- cell_types[i]
    
    # Create 2x2 table: this cell type vs all others
    cells_in_type <- contingency[, i]
    cells_not_in_type <- group_totals - cells_in_type
    
    test_table <- cbind(cells_in_type, cells_not_in_type)
    
    # Calculate fold change (ratio of proportions)
    proportions <- cells_in_type / group_totals
    fc <- proportions[2] / proportions[1]  # Group 2 / Group 1
    
    # Fisher exact test
    test_result <- fisher.test(as.matrix(test_table), alternative = "two.sided")
    
    fc_results[i] <- fc
    pval_results[i] <- test_result$p.value
  }
  
  # Adjust p-values
  qval_results <- p.adjust(pval_results, method = "BH")
  
  # Compile results
  results <- data.frame(
    fold_change = fc_results,
    p_value = pval_results,
    q_value = qval_results,
    row.names = cell_types
  )
  
  # Sort by p-value
  results <- results[order(results$p_value), ]
  
  return(results)
}

################################################################################
# Main Analysis
################################################################################

message("==========================================")
message("Cell Abundance Fisher Exact Tests")
message("==========================================\n")

################################################################################
# Part 1: PFD vs CD (Overall)
################################################################################

message("\n=== PART 1: PFD vs CD (Overall) ===\n")

# Load data
message("Loading data...")
total_file <- file.path(RAW_DATA_DIR, "total_no_umap.rdata")
validate_file_exists(total_file)

# Load the RDS file
total3 <- get(load(total_file))
message("Loaded object(s): ", paste(loaded_obj, collapse = ", "))


message("Total cells: ", ncol(total3))

# Create PFD vs CD grouping
message("\nCreating PFD vs CD groups...")

# Initialize in metadata
total3@meta.data$pfd <- "PFD"
total3@meta.data$pfd[total3@meta.data$new_group %in% c(3, 5)] <- "CD"

message("PFD cells: ", sum(total3@meta.data$pfd == "PFD"))
message("CD cells: ", sum(total3@meta.data$pfd == "CD"))

# Test all annotation levels
annotation_levels <- c("maincelltype", "anotation_group", 
                      "anotation_intermediate", "anotation_refined")

for (annot in annotation_levels) {
  
  results <- fisher_abundance(total3, annot, group_var = "pfd")
  
  # Save results
  output_file <- file.path(
    TABLES_DIR,
    paste0("pfd_cd_abundance_fisher_", annot, ".txt")
  )
  save_table_standard(results, output_file)
  
  # Report significant results
  sig_results <- results[results$q_value < 0.05, ]
  message("\n  Significant cell types (q < 0.05): ", nrow(sig_results))
  if (nrow(sig_results) > 0) {
    for (i in 1:min(5, nrow(sig_results))) {
      message("    ", rownames(sig_results)[i], 
              ": FC = ", round(sig_results$fold_change[i], 3),
              ", q = ", format(sig_results$q_value[i], scientific = TRUE, digits = 3))
    }
  }
}

################################################################################
# Part 2: PFD Uninflamed vs CD Uninflamed (6 vs 5)
################################################################################

message("\n=== PART 2: PFD Uninflamed vs CD Uninflamed (6 vs 5) ===\n")

# Subset to groups 6 and 5
total3_65 <- subset(total3, subset = new_group %in% c(6, 5))
message("Cells in comparison: ", ncol(total3_65))

# Set factor levels
total3_65@meta.data$new_group <- factor(total3_65@meta.data$new_group, levels = c("5", "6"))

message("Group 5 (CD uninflamed): ", sum(total3_65@meta.data$new_group == "5"))
message("Group 6 (PFD uninflamed): ", sum(total3_65@meta.data$new_group == "6"))

# Check for cell types with very few cells
message("\nChecking cell type counts...")
for (annot in c("anotation_intermediate", "anotation_group")) {
  counts <- table(total3_65@meta.data[[annot]])
  low_count <- counts[counts <= 10]
  if (length(low_count) > 0) {
    message("  ", annot, " - cell types with ≤10 cells: ", length(low_count))
    message("    ", paste(names(low_count), collapse = ", "))
  }
}

# Test all annotation levels
for (annot in annotation_levels) {
  
  results <- fisher_abundance(total3_65, annot, group_var = "new_group")
  
  # Save results
  output_file <- file.path(
    TABLES_DIR,
    paste0("6_5_abundance_fisher_", annot, ".txt")
  )
  save_table_standard(results, output_file)
  
  # Report significant results
  sig_results <- results[results$q_value < 0.05, ]
  message("\n  Significant cell types (q < 0.05): ", nrow(sig_results))
  if (nrow(sig_results) > 0) {
    for (i in 1:min(5, nrow(sig_results))) {
      message("    ", rownames(sig_results)[i], 
              ": FC = ", round(sig_results$fold_change[i], 3),
              ", q = ", format(sig_results$q_value[i], scientific = TRUE, digits = 3))
    }
  }
}

################################################################################
# Part 3: CD Inflamed vs PFD Inflamed (3 vs 4)
################################################################################

message("\n=== PART 3: CD Inflamed vs PFD Inflamed (3 vs 4) ===\n")

# Subset to groups 3 and 4
total3_34 <- subset(total3, subset = new_group %in% c(3, 4))
message("Cells in comparison: ", ncol(total3_34))

# Set factor levels
total3_34@meta.data$new_group <- factor(total3_34@meta.data$new_group, levels = c("4", "3"))

message("Group 4 (PFD inflamed): ", sum(total3_34@meta.data$new_group == "4"))
message("Group 3 (CD inflamed): ", sum(total3_34@meta.data$new_group == "3"))

# Check for cell types with very few cells
message("\nChecking cell type counts...")
for (annot in c("anotation_intermediate", "anotation_group")) {
  counts <- table(total3_34@meta.data[[annot]])
  low_count <- counts[counts <= 10]
  if (length(low_count) > 0) {
    message("  ", annot, " - cell types with ≤10 cells: ", length(low_count))
  }
}

# Test all annotation levels
for (annot in annotation_levels) {
  
  results <- fisher_abundance(total3_34, annot, group_var = "new_group")
  
  # Save results
  output_file <- file.path(
    TABLES_DIR,
    paste0("3_4_abundance_fisher_", annot, ".txt")
  )
  save_table_standard(results, output_file)
  
  # Report significant results
  sig_results <- results[results$q_value < 0.05, ]
  message("\n  Significant cell types (q < 0.05): ", nrow(sig_results))
  if (nrow(sig_results) > 0) {
    for (i in 1:min(5, nrow(sig_results))) {
      message("    ", rownames(sig_results)[i], 
              ": FC = ", round(sig_results$fold_change[i], 3),
              ", q = ", format(sig_results$q_value[i], scientific = TRUE, digits = 3))
    }
  }
}

################################################################################
# Generate Visualizations
################################################################################

message("\n--- Generating visualizations ---")

#' Create Barplot of Log2 Fold Changes
#'
#' @param results_file Path to results file
#' @param output_file Path to output PDF
#' @param title Plot title
#' @param seurat_obj Seurat object for ordering
create_fc_barplot <- function(results_file, output_file, title, seurat_obj) {
  
  if (!file.exists(results_file)) return()
  
  results <- read.table(results_file, check.names = FALSE, row.names = 1)
  
  # Calculate log2 FC
  log2fc <- log2(results$fold_change)
  names(log2fc) <- rownames(results)
  
  # Handle infinite values
  log2fc[is.infinite(log2fc)] <- ifelse(
    log2fc[is.infinite(log2fc)] > 0,
    max(log2fc[!is.infinite(log2fc)]),
    min(log2fc[!is.infinite(log2fc)])
  )
  
  # Order by cell type hierarchy if possible
  if ("anotation_group" %in% names(results_file)) {
    # Order by main cell type
    metadata <- seurat_obj@meta.data
    epi <- unique(metadata$anotation_group[metadata$maincelltype == "Epi"])
    stroma <- unique(metadata$anotation_group[metadata$maincelltype == "Stroma"])
    tcell <- unique(metadata$anotation_group[metadata$maincelltype == "Tcell"])
    bcell <- unique(metadata$anotation_group[metadata$maincelltype == "Bcell"])
    myeloid <- unique(metadata$anotation_group[metadata$maincelltype == "Myeloid"])
    
    cell_order <- c(epi, stroma, tcell, bcell, myeloid)
    cell_order <- intersect(cell_order, names(log2fc))
    
    if (length(cell_order) > 0) {
      log2fc <- log2fc[cell_order]
    }
  }
  
  # Create plot
  pdf(output_file, width = 10, height = 6)
  par(mar = c(10, 5, 4, 2))
  
  barplot(
    log2fc,
    las = 2,
    ylab = "Log2 Fold Change",
    main = title,
    col = ifelse(log2fc > 0, "red", "blue"),
    cex.names = 0.7
  )
  abline(h = 0, lty = 2)
  
  dev.off()
  
  message("  Saved: ", output_file)
}

# Generate plots for each comparison and annotation level
message("\nGenerating fold change barplots...")

# PFD vs CD
for (annot in annotation_levels) {
  results_file <- file.path(TABLES_DIR, paste0("pfd_cd_abundance_fisher_", annot, ".txt"))
  plot_file <- file.path(FIGURES_DIR, paste0("pfd_cd_abundance_fisher_", annot, ".pdf"))
  create_fc_barplot(results_file, plot_file, paste0("PFD vs CD: ", annot), total3)
}

# 6 vs 5
for (annot in annotation_levels) {
  results_file <- file.path(TABLES_DIR, paste0("6_5_abundance_fisher_", annot, ".txt"))
  plot_file <- file.path(FIGURES_DIR, paste0("6_5_abundance_fisher_", annot, ".pdf"))
  create_fc_barplot(results_file, plot_file, 
                   paste0("PFD Uninflamed vs CD Uninflamed (6 vs 5): ", annot), total3_65)
}

# 3 vs 4
for (annot in annotation_levels) {
  results_file <- file.path(TABLES_DIR, paste0("3_4_abundance_fisher_", annot, ".txt"))
  plot_file <- file.path(FIGURES_DIR, paste0("3_4_abundance_fisher_", annot, ".pdf"))
  create_fc_barplot(results_file, plot_file, 
                   paste0("CD Inflamed vs PFD Inflamed (3 vs 4): ", annot), total3_34)
}

################################################################################
# Generate Summary Table
################################################################################

message("\n--- Generating summary table ---")

summary_data <- data.frame(
  comparison = character(),
  annotation = character(),
  n_cell_types = integer(),
  n_significant = integer(),
  top_enriched = character(),
  top_depleted = character(),
  stringsAsFactors = FALSE
)

comparisons <- list(
  list(name = "pfd_cd", title = "PFD vs CD"),
  list(name = "6_5", title = "6 vs 5"),
  list(name = "3_4", title = "3 vs 4")
)

for (comp_info in comparisons) {
  comp_name <- comp_info$name
  comp_title <- comp_info$title
  
  for (annot in annotation_levels) {
    
    results_file <- file.path(
      TABLES_DIR,
      paste0(comp_name, "_abundance_fisher_", annot, ".txt")
    )
    
    if (!file.exists(results_file)) next
    
    results <- read.table(results_file, check.names = FALSE, row.names = 1)
    sig_results <- results[results$q_value < 0.05, ]
    
    # Get top enriched and depleted
    if (nrow(sig_results) > 0) {
      enriched <- sig_results[sig_results$fold_change > 1, ]
      depleted <- sig_results[sig_results$fold_change < 1, ]
      
      top_enriched <- if (nrow(enriched) > 0) rownames(enriched)[1] else "None"
      top_depleted <- if (nrow(depleted) > 0) rownames(depleted)[1] else "None"
    } else {
      top_enriched <- "None"
      top_depleted <- "None"
    }
    
    summary_data <- rbind(summary_data, data.frame(
      comparison = comp_title,
      annotation = annot,
      n_cell_types = nrow(results),
      n_significant = nrow(sig_results),
      top_enriched = top_enriched,
      top_depleted = top_depleted,
      stringsAsFactors = FALSE
    ))
  }
}

# Save summary
summary_file <- file.path(TABLES_DIR, paste0(OUTPUT_PREFIX, "_summary.txt"))
save_table_standard(summary_data, summary_file)

message("\n--- Summary Table ---")
print(summary_data)

################################################################################
# Summary Report
################################################################################

message("\n==========================================")
message("Analysis Summary")
message("==========================================")
message("Comparisons performed:")
message("  1. PFD vs CD (overall)")
message("  2. PFD uninflamed vs CD uninflamed (6 vs 5)")
message("  3. CD inflamed vs PFD inflamed (3 vs 4)")
message("\nAnnotation levels tested:")
for (annot in annotation_levels) {
  message("  - ", annot)
}
message("\nTotal significant results (q < 0.05): ", sum(summary_data$n_significant))
message("\nOutput files:")
message("  Results tables: ", TABLES_DIR)
message("  Barplots: ", FIGURES_DIR)

# Save session info
save_session_info(file.path(RESULTS_DIR, paste0(OUTPUT_PREFIX, "_session_info.txt")))

