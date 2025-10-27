################################################################################
# IL-13 Module Scoring and Signature Analysis
# IBD Transcriptomics Analysis
#
# Purpose: Score IL-13 response signatures across epithelial cell types
# Input: Epithelial scRNA-seq data, IL-13 microarray DEG results
# Output: Module scores per cell, signature comparisons across patient groups
#
# Original file: IL13_microarray/addmodulescore.docx (date: 20250302)
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
OUTPUT_PREFIX <- "il13_module_scoring"
EPI_DATA_PATH <- file.path(RAW_DATA_DIR, "epi2_badrm_harmony50.rds")

# IL-13 DEG thresholds
IL13_FC_THRESHOLD_1 <- 1
IL13_QVAL_THRESHOLD_1 <- 0.05
IL13_QVAL_THRESHOLD_2 <- 0.01

################################################################################
# Load Data
################################################################################

message("==========================================")
message("IL-13 Module Scoring Analysis")
message("==========================================\n")

# Load epithelial scRNA-seq data
message("Loading epithelial data from: ", EPI_DATA_PATH)
validate_file_exists(EPI_DATA_PATH)
load(EPI_DATA_PATH)

epi_data <- d
rm(d)

message("Total epithelial cells: ", ncol(epi_data))
message("Total genes: ", nrow(epi_data))

# Load IL-13 microarray DEG results
deg_file <- file.path(DEG_DIR, "il13_microarray_deseq2_results.txt")
message("\nLoading IL-13 microarray DEG results from: ", deg_file)
validate_file_exists(deg_file)

deg_results <- read.table(deg_file, check.names = FALSE, header = TRUE)
message("Loaded ", nrow(deg_results), " genes from IL-13 DEG results")

################################################################################
# Extract IL-13 Response Gene Sets
################################################################################

message("\n--- Extracting IL-13 response gene signatures ---")

# Get epithelial genes (for filtering)
epi_genes <- rownames(epi_data)

# Gene set 1: FC > 1, q-value < 0.05
genes_qval05 <- rownames(deg_results)[
  deg_results$log2FoldChange > IL13_FC_THRESHOLD_1 &
  deg_results$padj < IL13_QVAL_THRESHOLD_1
]
genes_qval05 <- na.omit(genes_qval05)
genes_qval05_filtered <- intersect(epi_genes, genes_qval05)
message("IL-13 genes (FC>1, q<0.05): ", length(genes_qval05))
message("  Present in epithelial data: ", length(genes_qval05_filtered))

# Gene set 2: FC > 1, q-value < 0.01 (more stringent)
genes_qval01 <- rownames(deg_results)[
  deg_results$log2FoldChange > IL13_FC_THRESHOLD_1 &
  deg_results$padj < IL13_QVAL_THRESHOLD_2
]
genes_qval01 <- na.omit(genes_qval01)
genes_qval01_filtered <- intersect(epi_genes, genes_qval01)
message("IL-13 genes (FC>1, q<0.01): ", length(genes_qval01))
message("  Present in epithelial data: ", length(genes_qval01_filtered))

# Save gene lists
genelist_file_1 <- file.path(TABLES_DIR, paste0(OUTPUT_PREFIX, "_fc1_qval0.05_genes.txt"))
save_gene_list(genes_qval05, genelist_file_1)

genelist_file_2 <- file.path(TABLES_DIR, paste0(OUTPUT_PREFIX, "_fc1_qval0.01_genes.txt"))
save_gene_list(genes_qval01, genelist_file_2)

################################################################################
# Calculate Module Scores
################################################################################

message("\n--- Calculating module scores ---")

# Initialize results matrix
module_scores <- matrix(nrow = 2, ncol = ncol(epi_data))
colnames(module_scores) <- colnames(epi_data)
rownames(module_scores) <- c("il13_micro_qval_0.05", "il13_micro_qval_0.01")

# Module score 1: q-value < 0.05
message("\nCalculating module score for q < 0.05 signature...")
if (length(genes_qval05_filtered) >= 5) {
  n_metadata_cols <- ncol(epi_data@meta.data)
  
  epi_scored <- AddModuleScore(
    epi_data,
    features = list(genes_qval05_filtered),
    ctrl = 5,
    name = "IL13_Response"
  )
  
  # Extract scores (take mean if multiple columns added)
  score_cols <- (n_metadata_cols + 1):ncol(epi_scored@meta.data)
  scores <- apply(epi_scored@meta.data[, score_cols, drop = FALSE], 1, mean)
  module_scores[1, ] <- scores
  
  message("  Mean score: ", round(mean(scores), 4))
  message("  Score range: ", round(min(scores), 4), " to ", round(max(scores), 4))
} else {
  warning("Too few genes for q < 0.05 signature (", length(genes_qval05_filtered), " genes)")
}

# Module score 2: q-value < 0.01
message("\nCalculating module score for q < 0.01 signature...")
if (length(genes_qval01_filtered) >= 5) {
  n_metadata_cols <- ncol(epi_data@meta.data)
  
  epi_scored <- AddModuleScore(
    epi_data,
    features = list(genes_qval01_filtered),
    ctrl = 5,
    name = "IL13_Response_Strict"
  )
  
  # Extract scores
  score_cols <- (n_metadata_cols + 1):ncol(epi_scored@meta.data)
  scores <- apply(epi_scored@meta.data[, score_cols, drop = FALSE], 1, mean)
  module_scores[2, ] <- scores
  
  message("  Mean score: ", round(mean(scores), 4))
  message("  Score range: ", round(min(scores), 4), " to ", round(max(scores), 4))
} else {
  warning("Too few genes for q < 0.01 signature (", length(genes_qval01_filtered), " genes)")
}

# Save module scores
output_file <- file.path(TABLES_DIR, paste0(OUTPUT_PREFIX, "_scores.txt"))
save_table_standard(module_scores, output_file)

################################################################################
# Compare Signatures Across Patient Groups
################################################################################

message("\n--- Comparing IL-13 signatures across patient groups ---")

# Prepare GMT list for signature comparison
gmt_list <- list()
gmt_list[[1]] <- module_scores
names(gmt_list) <- "il13_response"

# Get unique cell types
epi_data$anotation_intermediate <- as.character(epi_data$anotation_intermediate)
celltype_unique <- unique(epi_data$anotation_intermediate)
message("Cell types to analyze: ", length(celltype_unique))

annotation_level <- "intermediate"

# Loop through comparison sets
for (i in seq_along(COMPARISON_SETS)) {
  
  comp_set <- COMPARISON_SETS[[i]]
  group1 <- comp_set[1]
  group2 <- comp_set[2]
  
  message("\n=== Comparison: Group ", group1, " vs Group ", group2, " ===")
  
  # Determine which group variable to use
  if ((comp_set[1] == 1 || comp_set[1] == 2) && 
      (comp_set[2] == 1 || comp_set[2] == 2)) {
    group_order <- epi_data$group
  } else {
    group_order <- epi_data$new_group
  }
  
  # Process each GMT (signature set)
  for (w in seq_along(gmt_list)) {
    
    gmt_name <- names(gmt_list)[w]
    gmt_data <- gmt_list[[w]]
    signature_names <- rownames(gmt_data)
    
    message("  Signature set: ", gmt_name)
    
    # Initialize results matrix
    results <- matrix(nrow = 0, ncol = 5)
    
    # Test each cell type
    for (celltype in celltype_unique) {
      
      # Get cell indices for each group
      group1_indices <- which(
        epi_data$anotation_intermediate == celltype &
        group_order == group1
      )
      group2_indices <- which(
        epi_data$anotation_intermediate == celltype &
        group_order == group2
      )
      
      # Skip if insufficient cells
      if (length(group1_indices) == 0 || length(group2_indices) == 0) {
        next
      }
      
      # Get main cell type annotation
      main_celltype <- as.character(epi_data$maincelltype[group1_indices[1]]) 
      
      # Test each signature
      for (sig_idx in 1:nrow(gmt_data)) {
        
        sig_scores <- gmt_data[sig_idx, ]
        
        # Extract scores for each group
        group1_scores <- sig_scores[group1_indices]
        group2_scores <- sig_scores[group2_indices]
        
        # Calculate means
        if (length(group1_indices) == 1) {
          mean_group1 <- group1_scores
        } else {
          mean_group1 <- mean(group1_scores)
        }
        
        if (length(group2_indices) == 1) {
          mean_group2 <- group2_scores
        } else {
          mean_group2 <- mean(group2_scores)
        }
        
        # Wilcoxon test
        test_result <- wilcox.test(
          group1_scores,
          group2_scores,
          alternative = "two.sided"
        )
        
        # Store results
	if (length(group1_indices) > 0) {
  main_celltype <- as.character(epi_data$maincelltype[group1_indices[1]])
} else {
  main_celltype <- "Unknown"
}

# Then check length
cat("Creating result_row with", length(c(sig_name, test_result$p.value, mean_group1, mean_group2, celltype, main_celltype)), "elements\n")

	sig_name <- rownames(gmt_data)[sig_idx]
result_row <- c(
  sig_name,
  test_result$p.value,
  mean_group1,
  mean_group2,
  celltype,
  as.character(main_celltype)  # ← Force to character
)

	results <- rbind(results, result_row)
      }
    }
    
    # Check if we have results
    if (nrow(results) == 0) {
      message("    No cell types with sufficient cells")
      next
    }
    
    # Format results
    results_df <- as.data.frame(results, stringsAsFactors = FALSE)
    colnames(results_df) <- c("row.names", "pval", group1, group2, "cell", "maincelltype")
    
    # Convert numeric columns - NOW column 2 is pval, 3 is group1, 4 is group2
    results_df$pval <- as.numeric(results_df$pval)
    results_df[[as.character(group1)]] <- as.numeric(results_df[[as.character(group1)]])
    results_df[[as.character(group2)]] <- as.numeric(results_df[[as.character(group2)]])
    results_df$pval <- as.numeric(results_df$pval)
    
    # Save results
    output_file <- file.path(
      TABLES_DIR,
      paste0("signature_comp_", gmt_name, "_", annotation_level, "_", 
             group1, "_", group2, ".txt")
    )
    save_table_standard(results_df, output_file)
    
    # Report significant results
    sig_results <- results_df[results_df$pval < 0.05, ]
    message("    Cell types with p < 0.05: ", nrow(sig_results))
    
    if (nrow(sig_results) > 0) {
      for (j in 1:min(5, nrow(sig_results))) {
        message("      ", sig_results$cell[j], ": p = ", 
                format(sig_results$pval[j], scientific = TRUE, digits = 3))
      }
    }
  }
}

################################################################################
# Visualization - Module Score Distribution
################################################################################

message("\n--- Generating module score visualizations ---")

# Add module scores to Seurat object for visualization
epi_data$IL13_Score_Q05 <- module_scores["il13_micro_qval_0.05", ]
epi_data$IL13_Score_Q01 <- module_scores["il13_micro_qval_0.01", ]

# Feature plot
plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_featureplot.pdf"))
pdf(plot_file, width = 12, height = 5)

print(FeaturePlot(
  epi_data,
  features = c("IL13_Score_Q05", "IL13_Score_Q01"),
  ncol = 2
))

dev.off()
message("Feature plot saved to: ", plot_file)

# Violin plot by patient group
plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_violin_by_group.pdf"))
pdf(plot_file, width = 10, height = 5)
tryCatch({
	Idents(epi_data) <- factor(epi_data$new_group, levels = c(5, 6, 4, 3, 7))

	print(VlnPlot(
	  epi_data,
	  features = c("IL13_Score_Q05", "IL13_Score_Q01"),
	  ncol = 2,
  	  pt.size = 0
	  ))

	  dev.off()
	  message("Violin plot saved to: ", plot_file)
}, error = function(e) {
  # Violin plot by patient group (using boxplots)

# Extract data
plot_data <- data.frame(
  IL13_Score_Q05 = epi_data$IL13_Score_Q05,
  IL13_Score_Q01 = epi_data$IL13_Score_Q01,
  group = factor(epi_data$new_group, levels = c(5, 6, 4, 3, 7))
)

# Two-panel plot
par(mfrow = c(1, 2))

boxplot(IL13_Score_Q05 ~ group, data = plot_data,
        col = "lightblue",
        main = "IL13 Score Q05",
        xlab = "Patient Group",
        ylab = "Module Score",
        names = c("CD\nuninfl", "PFD\nuninfl", "PFD\ninfl", "CD\ninfl", "HC"))

boxplot(IL13_Score_Q01 ~ group, data = plot_data,
        col = "lightcoral",
        main = "IL13 Score Q01",
        xlab = "Patient Group",
        ylab = "Module Score",
        names = c("CD\nuninfl", "PFD\nuninfl", "PFD\ninfl", "CD\ninfl", "HC"))

	dev.off()
	message("Failed to save violin plot, instead made boxplot and saved to: ", plot_file)
}) 

################################################################################
# Summary Report
################################################################################

message("\n==========================================")
message("Analysis Summary")
message("==========================================")
message("Epithelial cells analyzed: ", ncol(epi_data))
message("IL-13 signature genes (q<0.05): ", length(genes_qval05_filtered))
message("IL-13 signature genes (q<0.01): ", length(genes_qval01_filtered))
message("Cell types analyzed: ", length(celltype_unique))
message("Patient group comparisons: ", length(COMPARISON_SETS))
message("\nOutput files:")
message("  Module scores: ", TABLES_DIR)
message("  Signature comparisons: ", TABLES_DIR)
message("  Visualizations: ", FIGURES_DIR)

# Save session info
save_session_info(file.path(RESULTS_DIR, paste0(OUTPUT_PREFIX, "_session_info.txt")))

message("\n==========================================")
message("IL-13 module scoring analysis complete!")
message("==========================================")
