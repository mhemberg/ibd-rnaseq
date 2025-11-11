################################################################################
# IL-13 Signature Heatmap Visualization
# IBD Transcriptomics Analysis
#
# Purpose: Generate heatmaps of IL-13 response signatures across epithelial cell types
# Input: IL-13 signature comparison results
# Output: Heatmaps showing signature enrichment and differences
################################################################################

# Load configuration and utilities
source("scripts/config.R")
source("scripts/utils_functions.R")

# Load required libraries
suppressPackageStartupMessages({
  library(gplots)
  library(Seurat)
})

# Script-specific parameters
OUTPUT_PREFIX <- "il13_signature_heatmap"
EPI_DATA_PATH <- file.path(RAW_DATA_DIR, "epi2_badrm_harmony50.rds")
COMPARISON_NAME <- "6_5"  # PFD uninflamed vs CD uninflamed

################################################################################
# Load Data
################################################################################

message("==========================================")
message("IL-13 Signature Heatmap Generation")
message("==========================================\n")

# Load epithelial data to get cell type order
message("Loading epithelial data for cell type ordering...")
validate_file_exists(EPI_DATA_PATH)
load(EPI_DATA_PATH)

epi_data <- d
rm(d)

# Get cell type order from data
message("Determining cell type order...")
celltype_table <- table(epi_data$anotation_group, epi_data$anotation_intermediate)

# Extract ordered cell types (non-zero entries)
cell_order <- character(0)
for (i in 1:nrow(celltype_table)) {
  non_zero_cells <- celltype_table[i, ][celltype_table[i, ] != 0]
  cell_order <- c(cell_order, names(non_zero_cells))
}

message("Cell type order (", length(cell_order), " types):")
for (cell in cell_order) {
  message("  ", cell)
}

################################################################################
# Load Signature Comparison Results
################################################################################

message("\n--- Loading signature comparison results ---")

# Load IL-13 signature scores
sig_file <- file.path(
  TABLES_DIR,
  paste0("signature_comp_il13_response_intermediate_", COMPARISON_NAME, ".txt")
)

message("Loading from: ", sig_file)
validate_file_exists(sig_file)

signature_data <- read.table(sig_file, check.names = FALSE, header = TRUE)
message("Loaded ", nrow(signature_data), " rows")

################################################################################
# Prepare Heatmap Data
################################################################################

message("\n--- Preparing heatmap data ---")


sig_data_filtered <- signature_data[which(signature_data[,1]=="il13_micro_qval_0.01"),]

message("Rows after filtering: ", nrow(sig_data_filtered))

# Set cell types as rownames
rownames(sig_data_filtered) <- sig_data_filtered$cell

# Reorder rows to match cell type order
sig_data_ordered <- sig_data_filtered[cell_order, ]

# Remove RBHI goblet cells (row 12 in original code)
# Identify which row this is
rbhi_idx <- which(rownames(sig_data_ordered) == "RBHI_goblet")
if (length(rbhi_idx) > 0) {
  sig_data_ordered <- sig_data_ordered[-rbhi_idx, ]
  message("Removed RBHI_goblet cells")
}

# Select columns for heatmap (groups 6 and 5 mean scores)
# Columns 4 and 3 in original correspond to group means
heatmap_data <- sig_data_ordered[, c(4, 3)]


# Verify data is valid
if (nrow(heatmap_data) == 0 || ncol(heatmap_data) == 0) {
  stop("No valid data for heatmap after cleaning")
}

message("Heatmap dimensions: ", nrow(heatmap_data), " x ", ncol(heatmap_data))

################################################################################
# Generate Heatmap 1: IL-13 Signature Scores by Group
################################################################################

message("\n--- Generating signature score heatmap ---")

plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_scores_heatmap.pdf"))
pdf(plot_file, width = 7, height = 10)

heatmap.2(
  as.matrix(heatmap_data),
  col = bluered(100),
  trace = "none",
  Rowv = FALSE,           # No row clustering
  Colv = FALSE,           # No column clustering
  dendrogram = "none",    # No dendrogram
  margins = c(15, 10),    # Increased margins for labels
  scale = "none",         # No scaling
  main = "IL-13 Response Signature\nPFD Uninflamed vs CD Uninflamed",
  xlab = "Patient Group",
  ylab = "Cell Type",
  key.title = "Score",
  key.xlab = "Module Score",
  cexRow = 0.8,
  cexCol = 1.2,
  labCol = c("PFD Uninfl (6)", "CD Uninfl (5)")
)

dev.off()
message("Signature score heatmap saved to: ", plot_file)

################################################################################
# Generate Heatmap 2: Difference Between Groups
################################################################################

message("\n--- Generating difference heatmap ---")

# Calculate difference (PFD uninflamed - CD uninflamed)
score_difference <- as.numeric(heatmap_data[, "6"]) - as.numeric(heatmap_data[, "5"])
names(score_difference) <- rownames(heatmap_data)

# Create matrix for heatmap (duplicate column for visualization)
diff_matrix <- cbind(score_difference, score_difference)
rownames(diff_matrix) <- names(score_difference)
colnames(diff_matrix) <- c("Difference", "Difference")

# Report statistics
message("\nDifference statistics:")
message("  Mean difference: ", round(mean(score_difference), 4))
message("  Range: ", round(min(score_difference), 4), " to ", 
        round(max(score_difference), 4))

# Cell types with largest positive difference
top_positive <- head(sort(score_difference, decreasing = TRUE), 5)
message("\nTop 5 cell types (higher in PFD uninflamed):")
for (i in 1:length(top_positive)) {
  message("  ", names(top_positive)[i], ": ", round(top_positive[i], 4))
}

# Cell types with largest negative difference
top_negative <- head(sort(score_difference, decreasing = FALSE), 5)
message("\nTop 5 cell types (higher in CD uninflamed):")
for (i in 1:length(top_negative)) {
  message("  ", names(top_negative)[i], ": ", round(top_negative[i], 4))
}

# Generate heatmap
plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_difference_heatmap.pdf"))
pdf(plot_file, width = 5, height = 10)

heatmap.2(
  as.matrix(diff_matrix),
  col = bluered(100),
  trace = "none",
  Rowv = FALSE,
  Colv = FALSE,
  dendrogram = "none",
  margins = c(15, 10),
  scale = "none",
  main = "IL-13 Response Difference\n(PFD Uninfl - CD Uninfl)",
  xlab = "",
  ylab = "Cell Type",
  key.title = "Difference",
  key.xlab = "Score Difference",
  cexRow = 0.8,
  cexCol = 1.2,
  labCol = c("PFD vs CD", "")
)

dev.off()
message("Difference heatmap saved to: ", plot_file)

################################################################################
# Generate Barplot of Differences
################################################################################

message("\n--- Generating difference barplot ---")

# Sort by difference
diff_sorted <- sort(score_difference)

plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_difference_barplot.pdf"))
pdf(plot_file, width = 8, height = 10)

par(mar = c(5, 12, 4, 2))  # Increase left margin for labels

barplot(
  diff_sorted,
  horiz = TRUE,
  las = 1,
  main = "IL-13 Response Difference by Cell Type",
  xlab = "Score Difference (PFD Uninfl - CD Uninfl)",
  col = ifelse(diff_sorted > 0, "red", "blue"),
  border = NA
)
abline(v = 0, lty = 2, col = "black")

dev.off()
message("Difference barplot saved to: ", plot_file)

################################################################################
# Save Processed Data
################################################################################

message("\n--- Saving processed data ---")

# Save ordered signature data
output_file <- file.path(TABLES_DIR, paste0(OUTPUT_PREFIX, "_ordered_scores.txt"))
save_table_standard(heatmap_data, output_file)

# Save difference values
diff_df <- data.frame(
  cell_type = names(diff_sorted),
  score_difference = diff_sorted,
  row.names = names(diff_sorted)
)
output_file <- file.path(TABLES_DIR, paste0(OUTPUT_PREFIX, "_score_differences.txt"))
save_table_standard(diff_df, output_file)

################################################################################
# Additional Comparisons (Optional)
################################################################################

message("\n--- Generating heatmaps for additional comparisons ---")

# Additional comparison: CD inflamed vs PFD inflamed (3 vs 4)
additional_comps <- c("3_4", "4_5")

for (comp in additional_comps) {
  
  sig_file_comp <- file.path(
    TABLES_DIR,
    paste0("signature_comp_il13_response_intermediate_", comp, ".txt")
  )
  
  if (!file.exists(sig_file_comp)) {
    message("  Skipping comparison ", comp, " (file not found)")
    next
  }
  
  message("\n  Processing comparison: ", comp)
  
  # Load data
  sig_data_comp <- read.table(sig_file_comp, check.names = FALSE, 
                              row.names = NULL, header = TRUE)
  
  # Filter and process
  sig_filtered_comp <- sig_data_comp[
    sig_data_comp$row.names == "il13_micro_qval_0.01",
  ]
  
  if (nrow(sig_filtered_comp) == 0) {
    message("    No data for this comparison")
    next
  }
  
  rownames(sig_filtered_comp) <- sig_filtered_comp$cell
  
  # Get group numbers
  groups <- strsplit(comp, "_")[[1]]
  group1 <- groups[1]
  group2 <- groups[2]
  
  # Create heatmap
  if (group1 %in% colnames(sig_filtered_comp) && 
      group2 %in% colnames(sig_filtered_comp)) {
    
    heatmap_data_comp <- sig_filtered_comp[, c(group1, group2)]
    
    plot_file_comp <- file.path(
      FIGURES_DIR,
      paste0(OUTPUT_PREFIX, "_scores_heatmap_", comp, ".pdf")
    )
    
    pdf(plot_file_comp, width = 7, height = 10)
    
    heatmap.2(
      as.matrix(heatmap_data_comp),
      col = bluered(100),
      trace = "none",
      Rowv = FALSE,
      Colv = FALSE,
      dendrogram = "none",
      margins = c(15, 10),
      scale = "none",
      main = paste0("IL-13 Response Signature\nGroup ", group1, " vs Group ", group2),
      cexRow = 0.8,
      cexCol = 1.2
    )
    
    dev.off()
    message("    Heatmap saved to: ", plot_file_comp)
  }
}

################################################################################
# Summary Report
################################################################################

message("\n==========================================")
message("Analysis Summary")
message("==========================================")
message("Cell types visualized: ", nrow(heatmap_data))
message("Primary comparison: PFD uninflamed vs CD uninflamed (6 vs 5)")
message("Mean IL-13 score difference: ", round(mean(score_difference), 4))
message("\nOutput files:")
message("  Heatmaps: ", FIGURES_DIR)
message("  Data tables: ", TABLES_DIR)

# Save session info
save_session_info(file.path(RESULTS_DIR, paste0(OUTPUT_PREFIX, "_session_info.txt")))

