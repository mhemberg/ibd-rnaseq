################################################################################
# CellPhoneDB Cell-Cell Interaction Analysis
# IBD Transcriptomics Analysis
#
# Purpose: Analyze cell-cell interactions using ligand-receptor pairs
# Input: CellPhoneDB results, scRNA-seq data
# Output: Interaction matrices, heatmaps, and network visualizations
#
# Note: This script assumes CellPhoneDB has been run externally
# See: https://www.cellphonedb.org/
################################################################################

# Load configuration and utilities
source("scripts/config.R")
source("scripts/utils_functions.R")

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(ggplot2)
  library(RColorBrewer)
})

# Script-specific parameters
OUTPUT_PREFIX <- "cellphonedb"
CELLPHONEDB_DIR <- file.path(RAW_DATA_DIR, "cellphonedb")

################################################################################
# Custom Functions
################################################################################

#' Create Ligand-Receptor Matrix for T cells
#'
#' @param cellphone_data CellPhoneDB results data frame
#' @param interacting_pair_list List of interacting pairs to analyze
#' @return Matrix of interaction scores (rows: source, columns: target)
lr_mat_tcell <- function(cellphone_data, interacting_pair_list) {
  
  # Get unique cell types
  cluster_table <- table(cellphone_data$Cluster_1, cellphone_data$Subset_1)
  cluster_subset <- apply(cluster_table, 1, function(x) colnames(cluster_table)[which(x != 0)])
  cell_list <- names(sort(cluster_subset))
  
  # Initialize result matrix
  result_matrix <- matrix(0, nrow = length(cell_list), ncol = length(cell_list))
  rownames(result_matrix) <- cell_list
  colnames(result_matrix) <- cell_list
  
  # Fill matrix
  for (i in seq_along(cell_list)) {
    cell1 <- cell_list[i]
    
    for (j in seq_along(cell_list)) {
      cell2 <- cell_list[j]
      
      # Calculate mean interaction score for specified pairs
      interaction_values <- c()
      
      for (int_pair in interacting_pair_list) {
        # Filter data for this cell pair
        tmp_data <- cellphone_data[
          cellphone_data$Cluster_1 == cell1 & cellphone_data$Cluster_2 == cell2,
        ]
        
        # Check if interaction pair exists
        if (int_pair %in% tmp_data$interacting_pair) {
          val <- tmp_data$mean[tmp_data$interacting_pair == int_pair]
        } else {
          val <- 0
        }
        
        interaction_values <- c(interaction_values, val)
      }
      
      # Store mean value
      result_matrix[i, j] <- mean(interaction_values)
    }
  }
  
  return(result_matrix)
}

#' Calculate Maximum Value Across All Groups
#'
#' @param interacting_pair_list List of interacting pairs
#' @return Maximum interaction value
max_eval <- function(interacting_pair_list) {
  
  group_list <- c(3, 4, 5, 6, 7)
  max_values <- c()
  
  for (grp in group_list) {
    
    # Load CellPhoneDB results
    cp_file <- file.path(
      CELLPHONEDB_DIR,
      paste0("group", grp, "_toda_info.csv")
    )
    
    if (!file.exists(cp_file)) {
      warning("CellPhoneDB file not found: ", cp_file)
      next
    }
    
    cellphone_data <- read.table(cp_file, sep = ";", dec = ",", header = TRUE)
    
    # Calculate interaction matrix
    lr_matrix <- lr_mat_tcell(cellphone_data, interacting_pair_list)
    
    max_values <- c(max_values, max(lr_matrix))
  }
  
  return(max(max_values))
}

################################################################################
# Main Analysis
################################################################################

message("==========================================")
message("CellPhoneDB Interaction Analysis")
message("==========================================\n")

# Check if CellPhoneDB directory exists
if (!dir.exists(CELLPHONEDB_DIR)) {
  message("CellPhoneDB directory not found: ", CELLPHONEDB_DIR)
  message("\nSkipping CellPhoneDB analysis...")
  quit(save = "no")
}

################################################################################
# LTB (Lymphotoxin Beta) Analysis
################################################################################

message("\n--- Analyzing LTB interactions ---")

# Define interaction pairs
ltb_pairs <- c("LTBR_LTB")

# Calculate maximum value for scaling
message("Calculating maximum interaction value...")
max_val <- max_eval(ltb_pairs)
message("Maximum value: ", round(max_val, 4))

# Process each group
group_list <- c(3, 4, 5, 6, 7)

for (grp in group_list) {
  
  message("\n=== Processing Group ", grp, " ===")
  
  # Load CellPhoneDB results
  cp_file <- file.path(
    CELLPHONEDB_DIR,
    paste0("group", grp, "_toda_info.csv")
  )
  
  if (!file.exists(cp_file)) {
    message("File not found, skipping...")
    next
  }
  
  cellphone_data <- read.table(cp_file, sep = ";", dec = ",", header = TRUE)
  
  # Load corresponding CellChat object for visualization
  cellchat_file <- file.path(
    PROCESSED_DATA_DIR,
    paste0("cellchat_intermediate_anot_group", grp)
  )
  
  if (!file.exists(cellchat_file)) {
    message("CellChat object not found, skipping visualization...")
    next
  }
  
  load(cellchat_file)
  
  # Calculate interaction matrix
  lr_matrix <- lr_mat_tcell(cellphone_data, ltb_pairs)
  
  # Scale to maximum value for comparison
  lr_matrix[1, 1] <- max_val
  
  # Update CellChat object
  cellchat@net$weight <- lr_matrix
  
  # Generate heatmap
  plot_file <- file.path(
    FIGURES_DIR,
    paste0(OUTPUT_PREFIX, "_ltb_heatmap_group", grp, ".pdf")
  )
  
  pdf(plot_file, width = 12, height = 12)
  
  tryCatch({
    heatmap_plot <- netVisual_heatmap(
      cellchat,
      measure = "weight",
      color.heatmap = "Reds"
    )
    print(heatmap_plot)
  }, error = function(e) {
    message("Error generating heatmap: ", e$message)
  })
  
  dev.off()
  
  message("Heatmap saved: ", plot_file)
  
  # Save interaction matrix
  matrix_file <- file.path(
    TABLES_DIR,
    paste0(OUTPUT_PREFIX, "_ltb_matrix_group", grp, ".txt")
  )
  save_table_standard(lr_matrix, matrix_file)
}

################################################################################
# Additional Interaction Pairs
################################################################################

message("\n--- Analyzing additional interactions ---")

# Define additional interaction pairs of interest
additional_pairs <- list(
  IL22 = c("IL22RA1_IL22", "IL22RA2_IL22"),
  TNFSF15 = c("TNFRSF25_TNFSF15"),
  IL13 = c("IL13RA1_IL13", "IL13RA2_IL13")
)

for (pair_name in names(additional_pairs)) {
  
  message("\n--- Processing ", pair_name, " interactions ---")
  
  pairs <- additional_pairs[[pair_name]]
  
  for (grp in group_list) {
    
    cp_file <- file.path(
      CELLPHONEDB_DIR,
      paste0("group", grp, "_toda_info.csv")
    )
    
    if (!file.exists(cp_file)) next
    
    cellphone_data <- read.table(cp_file, sep = ";", dec = ",", header = TRUE)
    
    # Calculate interaction matrix
    lr_matrix <- lr_mat_tcell(cellphone_data, pairs)
    
    # Save matrix
    matrix_file <- file.path(
      TABLES_DIR,
      paste0(OUTPUT_PREFIX, "_", tolower(pair_name), "_matrix_group", grp, ".txt")
    )
    save_table_standard(lr_matrix, matrix_file)
    
    message("Group ", grp, " ", pair_name, " matrix saved")
  }
}

################################################################################
# Compare Interactions Across Groups
################################################################################

message("\n--- Comparing interactions across groups ---")

# For each interaction type, compare across groups
for (pair_name in c("ltb", names(additional_pairs))) {
  
  message("\n=== Comparing ", pair_name, " across groups ===")
  
  # Collect matrices from all groups
  matrices <- list()
  
  for (grp in group_list) {
    
    matrix_file <- file.path(
      TABLES_DIR,
      paste0(OUTPUT_PREFIX, "_", tolower(pair_name), "_matrix_group", grp, ".txt")
    )
    
    if (file.exists(matrix_file)) {
      matrices[[as.character(grp)]] <- read.table(matrix_file, check.names = FALSE)
    }
  }
  
  if (length(matrices) < 2) {
    message("Insufficient data for comparison")
    next
  }

# Calculate differences (e.g., group 6 vs group 5)
  if ("6" %in% names(matrices) && "5" %in% names(matrices)) {
    
    mat6 <- as.matrix(matrices[["6"]])
    mat5 <- as.matrix(matrices[["5"]])
    
    # Check if dimensions match
    if (!identical(dim(mat6), dim(mat5))) {
      message("Warning: Matrix dimensions don't match")
      message("  Group 6: ", nrow(mat6), " x ", ncol(mat6))
      message("  Group 5: ", nrow(mat5), " x ", ncol(mat5))
      
      # Find common cell types
      common_rows <- intersect(rownames(mat6), rownames(mat5))
      common_cols <- intersect(colnames(mat6), colnames(mat5))
      
      if (length(common_rows) == 0 || length(common_cols) == 0) {
        message("  No common cell types, skipping difference calculation")
        next
      }
      
      message("  Using ", length(common_rows), " common cell types")
      
      # Subset to common cell types
      mat6 <- mat6[common_rows, common_cols, drop = FALSE]
      mat5 <- mat5[common_rows, common_cols, drop = FALSE]
    }
    
    # Calculate difference
    diff_matrix <- mat6 - mat5
    
    # Save difference matrix
    diff_file <- file.path(
      TABLES_DIR,
      paste0(OUTPUT_PREFIX, "_", tolower(pair_name), "_diff_6vs5.txt")
    )
    save_table_standard(diff_matrix, diff_file)
    
    message("Difference matrix (6 vs 5) saved: ", diff_file)
  }
}

################################################################################
# Summary Report
################################################################################

message("\n==========================================")
message("Analysis Summary")
message("==========================================")
message("Interaction pairs analyzed:")
message("  - LTB interactions")
for (name in names(additional_pairs)) {
  message("  - ", name, " interactions")
}
message("\nOutput files:")
message("  Interaction matrices: ", TABLES_DIR)
message("  Heatmaps: ", FIGURES_DIR)

# Save session info
save_session_info(file.path(RESULTS_DIR, paste0(OUTPUT_PREFIX, "_session_info.txt")))


