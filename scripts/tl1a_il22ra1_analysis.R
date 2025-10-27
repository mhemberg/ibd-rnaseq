################################################################################
# TNFSF15 (TL1A) and IL22RA1 Expression Analysis in Epithelial Cells
# IBD Transcriptomics Analysis
#
# Purpose: Analyze TL1A and IL-22 receptor expression in epithelial subsets
# Input: Epithelial scRNA-seq data
# Output: Expression statistics, DEG results, and visualizations
#
# Original files: 
#   - scRNAseq/il22/epi_il22_tnfsf15.docx (date: 20241231, 20250110)
#   - il22_microarray/epi_il22_tnfsf15.docx
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
OUTPUT_PREFIX <- "tl1a_il22ra1_epithelial"
EPI_DATA_PATH <- file.path(RAW_DATA_DIR, "epi2_badrm_harmony50.rds")

# Genes of interest
TL1A_GENE <- "TNFSF15"
IL22R_GENE <- "IL22RA1"

# TL1A-expressing cell types
TL1A_CELLTYPES <- c("PLCG2_enterocytes", "TA", "M_cells")
IL22R_CELLTYPES <- c("Entero_AQP8pos")

################################################################################
# Custom Functions
################################################################################

#' Perform Differential Expression with Zero/Non-zero Analysis
#' 
#' Tests expression differences between two groups, including non-zero cell analysis
#' 
#' @param seurat_obj Seurat object
#' @param cell_types Cell types to include
#' @param gene Gene name to test
#' @param group1_id Group 1 identifier (e.g., 6 for PFD uninflamed)
#' @param group2_id Group 2 identifier (e.g., 5 for CD uninflamed)
#' @return Data frame with comprehensive statistics
deg_special <- function(seurat_obj, cell_types, gene, group1_id, group2_id) {
  
  # Subset to selected cell types
  epi_subset <- subset(seurat_obj, subset = anotation_intermediate %in% cell_types)
  
  # Get expression data
  expr_matrix <- GetAssayData(epi_subset)
  gene_expr <- expr_matrix[gene, ]
  
  # Extract expression by group
  group1_expr <- gene_expr[which(epi_subset$new_group == group1_id)]
  group2_expr <- gene_expr[which(epi_subset$new_group == group2_id)]
  
  # Overall statistics
  mean_group1 <- mean(group1_expr)
  mean_group2 <- mean(group2_expr)
  
  # Statistical test (all cells)
  test_all <- wilcox.test(group1_expr, group2_expr, alternative = "two.sided")
  pval_all <- test_all$p.value
  
  # Non-zero cells only
  group1_nonzero <- group1_expr[group1_expr != 0]
  group2_nonzero <- group2_expr[group2_expr != 0]
  
  # Non-zero statistics
  pct_nonzero_group1 <- length(group1_nonzero) / length(group1_expr)
  pct_nonzero_group2 <- length(group2_nonzero) / length(group2_expr)
  mean_nonzero_group1 <- ifelse(length(group1_nonzero) > 0, mean(group1_nonzero), 0)
  mean_nonzero_group2 <- ifelse(length(group2_nonzero) > 0, mean(group2_nonzero), 0)
  
  # Statistical test (non-zero cells only)
  if (length(group1_nonzero) > 0 && length(group2_nonzero) > 0) {
    test_nonzero <- wilcox.test(group1_nonzero, group2_nonzero, alternative = "two.sided")
    pval_nonzero <- test_nonzero$p.value
  } else {
    pval_nonzero <- NA
  }
  
  # Compile results
  results <- c(
    mean_group1,
    mean_group2,
    pval_all,
    pct_nonzero_group1,
    pct_nonzero_group2,
    mean_nonzero_group1,
    mean_nonzero_group2,
    pval_nonzero
  )
  
  names(results) <- c(
    "mean6", "mean5", "pval",
    "nonzero6", "nonzero5",
    "nz_mean6", "nz_mean5", "nz_pval"
  )
  
  return(results)
}

################################################################################
# Load Data
################################################################################

message("==========================================")
message("TL1A and IL-22R Epithelial Expression Analysis")
message("==========================================\n")

# Load epithelial scRNA-seq data
message("Loading epithelial data from: ", EPI_DATA_PATH)
validate_file_exists(EPI_DATA_PATH)
load(EPI_DATA_PATH)

epi_data <- d
rm(d)

message("Total epithelial cells: ", ncol(epi_data))
message("Total genes: ", nrow(epi_data))

################################################################################
# Part 1: TNFSF15 (TL1A) Expression Analysis
################################################################################

message("\n==========================================")
message("TNFSF15 (TL1A) Expression Analysis")
message("==========================================\n")

# Get TL1A expression data
expr_matrix <- GetAssayData(epi_data)
tl1a_expr <- expr_matrix[TL1A_GENE, ]
names(tl1a_expr) <- epi_data$anotation_intermediate

# Calculate expression statistics by cell type
message("--- TL1A expression by cell type ---")

cell_types_unique <- unique(epi_data$anotation_intermediate)
nonzero_ratios <- numeric(length(cell_types_unique))
names(nonzero_ratios) <- cell_types_unique

for (celltype in cell_types_unique) {
  celltype_expr <- tl1a_expr[names(tl1a_expr) == celltype]
  nonzero_ratio <- length(which(celltype_expr != 0)) / length(celltype_expr)
  nonzero_ratios[celltype] <- nonzero_ratio
}

# Sort by ratio
nonzero_sorted <- sort(nonzero_ratios, decreasing = TRUE)

message("Top 10 cell types by TL1A+ cell percentage:")
for (i in 1:min(10, length(nonzero_sorted))) {
  message("  ", names(nonzero_sorted)[i], ": ", 
          round(100 * nonzero_sorted[i], 2), "%")
}

# Calculate pseudobulk expression
tl1a_matrix <- expr_matrix[TL1A_GENE, , drop = FALSE]
tl1a_pseudobulk <- get_pseudobulk(
  tl1a_matrix,
  as.character(epi_data$anotation_intermediate)
)

tl1a_pseudobulk_vec <- as.numeric(tl1a_pseudobulk)
names(tl1a_pseudobulk_vec) <- colnames(tl1a_pseudobulk)
tl1a_pseudobulk_sorted <- sort(tl1a_pseudobulk_vec)

# Save results
output_file <- file.path(TABLES_DIR, paste0(OUTPUT_PREFIX, "_tl1a_nonzero_ratio.txt"))
save_table_standard(
  data.frame(cell_type = names(nonzero_sorted), ratio = nonzero_sorted),
  output_file
)

output_file <- file.path(TABLES_DIR, paste0(OUTPUT_PREFIX, "_tl1a_mean_expression.txt"))
save_table_standard(
  data.frame(cell_type = names(tl1a_pseudobulk_sorted), 
             mean_expression = tl1a_pseudobulk_sorted),
  output_file
)

# Barplot of mean expression
plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_tl1a_mean_barplot.pdf"))
pdf(plot_file, width = 10, height = 6)

par(mar = c(8, 5, 4, 2))
barplot(
  tl1a_pseudobulk_sorted,
  las = 2,
  ylab = "Mean TL1A Expression",
  main = "TNFSF15 (TL1A) Expression by Cell Type",
  col = "steelblue",
  cex.names = 0.7
)

dev.off()
message("\nTL1A barplot saved to: ", plot_file)

################################################################################
# Part 2: TL1A Differential Expression Analysis
################################################################################

message("\n--- TL1A differential expression analysis ---")

# Subset to PFD uninflamed (6) and CD uninflamed (5) for comparison
epi_subset <- subset(epi_data, subset = new_group %in% c(6, 5))

# Test different cell type combinations
test_configs <- list(
  list(name = "PLCG2_enterocytes", celltypes = c("PLCG2_enterocytes")),
  list(name = "TA", celltypes = c("TA")),
  list(name = "PLCG2_TA_combined", celltypes = c("PLCG2_enterocytes", "TA")),
  list(name = "M_cells", celltypes = c("M_cells")),
  list(name = "All_TL1A_high", celltypes = c("PLCG2_enterocytes", "TA", "M_cells"))
)

deg_results_all <- list()

for (config in test_configs) {
  
  message("\nTesting: ", config$name)
  
  result <- deg_special(
    epi_subset,
    config$celltypes,
    TL1A_GENE,
    group1_id = 6,
    group2_id = 5
  )
  
  deg_results_all[[config$name]] <- result
  
  # Print results
  message("  Mean expression - PFD uninfl (6): ", round(result["mean6"], 4))
  message("  Mean expression - CD uninfl (5): ", round(result["mean5"], 4))
  message("  P-value (all cells): ", format(result["pval"], scientific = TRUE, digits = 3))
  message("  % TL1A+ cells - PFD uninfl: ", round(100 * result["nonzero6"], 2), "%")
  message("  % TL1A+ cells - CD uninfl: ", round(100 * result["nonzero5"], 2), "%")
  
  if (!is.na(result["nz_pval"])) {
    message("  P-value (TL1A+ cells only): ", 
            format(result["nz_pval"], scientific = TRUE, digits = 3))
  }
}

# Combine and save results
deg_results_df <- do.call(rbind, deg_results_all)
output_file <- file.path(DEG_DIR, paste0(OUTPUT_PREFIX, "_tl1a_deg_results.txt"))
save_table_standard(deg_results_df, output_file)

################################################################################
# Part 3: IL22RA1 Expression Analysis
################################################################################

message("\n==========================================")
message("IL22RA1 Expression Analysis")
message("==========================================\n")

# Calculate IL22RA1 pseudobulk expression
il22r_matrix <- expr_matrix[IL22R_GENE, , drop = FALSE]
il22r_pseudobulk <- get_pseudobulk(
  il22r_matrix,
  as.character(epi_data$anotation_intermediate)
)

il22r_pseudobulk_vec <- as.numeric(il22r_pseudobulk)
names(il22r_pseudobulk_vec) <- colnames(il22r_pseudobulk)
il22r_pseudobulk_sorted <- sort(il22r_pseudobulk_vec)

message("Top 10 cell types by IL22RA1 expression:")
top10_il22r <- tail(il22r_pseudobulk_sorted, 10)
for (i in length(top10_il22r):1) {
  message("  ", names(top10_il22r)[i], ": ", round(top10_il22r[i], 4))
}

# Save IL22RA1 results
output_file <- file.path(TABLES_DIR, paste0(OUTPUT_PREFIX, "_il22ra1_mean_expression.txt"))
save_table_standard(
  data.frame(cell_type = names(il22r_pseudobulk_sorted), 
             mean_expression = il22r_pseudobulk_sorted),
  output_file
)

# Barplot
plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_il22ra1_mean_barplot.pdf"))
pdf(plot_file, width = 10, height = 6)

par(mar = c(8, 5, 4, 2))
barplot(
  il22r_pseudobulk_sorted,
  las = 2,
  ylab = "Mean IL22RA1 Expression",
  main = "IL22RA1 Expression by Cell Type",
  col = "darkgreen",
  cex.names = 0.7
)

dev.off()
message("IL22RA1 barplot saved to: ", plot_file)

# IL22RA1 differential expression
message("\n--- IL22RA1 differential expression analysis ---")

il22r_result <- deg_special(
  epi_subset,
  IL22R_CELLTYPES,
  IL22R_GENE,
  group1_id = 6,
  group2_id = 5
)

message("Testing: ", paste(IL22R_CELLTYPES, collapse = ", "))
message("  Mean expression - PFD uninfl (6): ", round(il22r_result["mean6"], 4))
message("  Mean expression - CD uninfl (5): ", round(il22r_result["mean5"], 4))
message("  P-value: ", format(il22r_result["pval"], scientific = TRUE, digits = 3))

# Save IL22RA1 DEG results
output_file <- file.path(DEG_DIR, paste0(OUTPUT_PREFIX, "_il22ra1_deg_results.txt"))
save_table_standard(
  data.frame(t(il22r_result)),
  output_file
)

################################################################################
# Part 4: TL1A+ Cell Marker Analysis
################################################################################

message("\n--- Finding markers for TL1A+ epithelial cells ---")

# Create TL1A expression group
epi_data$tl1a_group <- "other"
epi_data$tl1a_group[epi_data$anotation_intermediate %in% TL1A_CELLTYPES] <- "tl1a_high"

# Set identity
Idents(epi_data) <- as.factor(epi_data$tl1a_group)

# Find markers
message("Finding markers (this may take a few minutes)...")
tl1a_markers <- FindAllMarkers(epi_data, only.pos = TRUE)

# Save markers
output_file <- file.path(DEG_DIR, paste0(OUTPUT_PREFIX, "_tl1a_high_cell_markers.txt"))
save_table_standard(tl1a_markers, output_file)

# Report TNFSF15 statistics
tl1a_stats <- tl1a_markers[tl1a_markers$gene == TL1A_GENE, ]
if (nrow(tl1a_stats) > 0) {
  message("\nTNFSF15 marker statistics:")
  message("  avg_log2FC: ", round(tl1a_stats$avg_log2FC, 2))
  message("  p_val_adj: ", format(tl1a_stats$p_val_adj, scientific = TRUE, digits = 3))
}

################################################################################
# Part 5: Visualizations
################################################################################

message("\n--- Generating visualizations ---")

# Boxplot: TL1A expression by patient group (TL1A-high cells)
tl1a_high_cells <- subset(epi_data, subset = anotation_intermediate %in% TL1A_CELLTYPES)
tl1a_high_expr <- GetAssayData(tl1a_high_cells)
tl1a_high_values <- tl1a_high_expr[TL1A_GENE, ]

plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_tl1a_boxplot.pdf"))
pdf(plot_file, width = 4.8, height = 5)

boxplot(
  tl1a_high_values[tl1a_high_cells$new_group == "5"],
  tl1a_high_values[tl1a_high_cells$new_group == "6"],
  tl1a_high_values[tl1a_high_cells$new_group == "3"],
  tl1a_high_values[tl1a_high_cells$new_group == "4"],
  names = c("CD uninfl", "PFD uninfl", "CD infl", "PFD infl"),
  ylab = "TL1A Expression",
  main = "TL1A in High-Expressing Epithelial Cells",
  pch = 20,
  col = c("lightblue", "lightcoral", "lightblue", "lightcoral")
)

dev.off()
message("TL1A boxplot saved to: ", plot_file)

# Violin plot: TL1A expression
Idents(tl1a_high_cells) <- factor(tl1a_high_cells$new_group, levels = c(5, 6, 3, 4, 7))

plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_tl1a_violin.pdf"))
pdf(plot_file, width = 7, height = 5)
tryCatch({
	print(VlnPlot(tl1a_high_cells, features = TL1A_GENE, pt.size = 0.5))

dev.off()
message("TL1A violin plot saved to: ", plot_file)
}, error = function(e) {
# TL1A boxplot
plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_tl1a_violin.pdf"))
pdf(plot_file, width = 7, height = 5)

# Extract expression
expr_data <- GetAssayData(tl1a_high_cells, slot = "data")
tl1a_values <- as.numeric(expr_data[TL1A_GENE, ])
groups <- factor(tl1a_high_cells$new_group, levels = c(5, 6, 3, 4, 7))

# Boxplot
boxplot(tl1a_values ~ groups,
        col = c("lightblue", "lightcoral", "lightblue", "lightcoral", "lightgreen"),
        main = paste0(TL1A_GENE, " Expression"),
        xlab = "Patient Group",
        ylab = "Expression Level",
        names = c("CD uninfl", "PFD uninfl", "CD infl", "PFD infl", "Healthy"))

dev.off()
})

# Boxplot: IL22RA1 expression
il22r_high_cells <- subset(epi_data, subset = anotation_intermediate %in% IL22R_CELLTYPES)

if (ncol(il22r_high_cells) > 0) {
  il22r_high_expr <- GetAssayData(il22r_high_cells)
  il22r_high_values <- il22r_high_expr[IL22R_GENE, ]
  
  plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_il22ra1_boxplot.pdf"))
  pdf(plot_file, width = 4.8, height = 5)
  
  boxplot(
    il22r_high_values[il22r_high_cells$new_group == "5"],
    il22r_high_values[il22r_high_cells$new_group == "6"],
    il22r_high_values[il22r_high_cells$new_group == "3"],
    il22r_high_values[il22r_high_cells$new_group == "4"],
    names = c("CD uninfl", "PFD uninfl", "CD infl", "PFD infl"),
    ylab = "IL22RA1 Expression",
    main = "IL22RA1 in Entero_AQP8pos Cells",
    pch = 20,
    col = c("lightblue", "lightcoral", "lightblue", "lightcoral")
  )
  
  dev.off()
  message("IL22RA1 boxplot saved to: ", plot_file)
  tryCatch({ # Violin plot

  Idents(il22r_high_cells) <- factor(il22r_high_cells$new_group, levels = c(5, 6, 3, 4, 7))
  
  plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_il22ra1_violin.pdf"))
  pdf(plot_file, width = 7, height = 5)
  
  print(VlnPlot(il22r_high_cells, features = IL22R_GENE, pt.size = 0.5))
  
  dev.off()
  message("IL22RA1 violin plot saved to: ", plot_file)
}, error = function(e) {
plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_il22ra1_violin.pdf"))
pdf(plot_file, width = 7, height = 5)

# Extract expression
expr_data <- GetAssayData(tl1a_high_cells, slot = "data")
il22r_values <- as.numeric(expr_data[IL22R_GENE, ])
groups <- factor(tl1a_high_cells$new_group, levels = c(5, 6, 3, 4, 7))

# Boxplot
boxplot(il22r_values ~ groups,
        col = c("lightblue", "lightcoral", "lightblue", "lightcoral", "lightgreen"),
        main = paste0(IL22R_GENE, " Expression"),
        xlab = "Patient Group",
        ylab = "Expression Level",
        names = c("CD uninfl", "PFD uninfl", "CD infl", "PFD infl", "Healthy"))

dev.off()
})

}

################################################################################
# Summary Report
################################################################################

message("\n==========================================")
message("Analysis Summary")
message("==========================================")
message("Total epithelial cells: ", ncol(epi_data))
message("\nTL1A (TNFSF15):")
message("  High-expressing cell types: ", paste(TL1A_CELLTYPES, collapse = ", "))
message("  Cells in TL1A-high types: ", ncol(tl1a_high_cells))
message("\nIL22RA1:")
message("  High-expressing cell types: ", paste(IL22R_CELLTYPES, collapse = ", "))
message("  Cells in IL22RA1-high types: ", ncol(il22r_high_cells))
message("\nOutput files:")
message("  Expression statistics: ", TABLES_DIR)
message("  DEG results: ", DEG_DIR)
message("  Visualizations: ", FIGURES_DIR)

# Save session info
save_session_info(file.path(RESULTS_DIR, paste0(OUTPUT_PREFIX, "_session_info.txt")))

message("\n==========================================")
message("TL1A and IL22RA1 analysis complete!")
message("==========================================")
