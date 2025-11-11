################################################################################
# Validation Dataset Differential Expression Analysis
# IBD Transcriptomics Analysis
#
# Purpose: Perform DEG analysis on GSE233063 validation data for enrichment tests
# Input: Processed validation Seurat object
# Output: DEG results for each treatment vs control
################################################################################

# Load configuration and utilities
source("scripts/config.R")
source("scripts/utils_functions.R")

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# Script-specific parameters
OUTPUT_PREFIX <- "validation_gse233063"
FC_THRESHOLD <- 1.4

################################################################################
# Main Analysis
################################################################################

message("==========================================")
message("Validation Dataset DEG Analysis")
message("==========================================\n")

# Load processed validation data
processed_file <- file.path(PROCESSED_DATA_DIR, paste0(OUTPUT_PREFIX, "_processed.rds"))
message("Loading processed data from: ", processed_file)
validate_file_exists(processed_file)

validation_obj <- readRDS(processed_file)
message("Loaded object with ", ncol(validation_obj), " cells")

# Set identity to treatment group
Idents(validation_obj) <- as.factor(validation_obj$treatment)

################################################################################
# DEG Analysis: Each Treatment vs Control
################################################################################

message("\n--- Performing DEG analysis for each treatment ---")

# Get treatment groups (excluding control)
treatments <- setdiff(unique(validation_obj$treatment), "ctrl")

for (treatment in treatments) {
  
  message("\n=== Analyzing: ", treatment, " vs ctrl ===")
  
  # Perform differential expression
  deg_results <- tryCatch({
    FindMarkers(
      validation_obj,
      ident.1 = treatment,
      ident.2 = "ctrl",
      logfc.threshold = 0.01,  # Low threshold for comprehensive analysis
      min.pct = 0.1,
      test.use = "wilcox",
      verbose = FALSE
    )
  }, error = function(e) {
    message("Error in DEG analysis for ", treatment, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(deg_results)) next
  
  # Add gene names as column
  deg_results$gene <- rownames(deg_results)
  
  # Save raw DEG results
  output_file <- file.path(DEG_DIR, paste0(OUTPUT_PREFIX, "_deg_raw_", treatment))
  save_table_standard(deg_results, output_file)
  
  message("  Total DEGs: ", nrow(deg_results))
  message("  Saved to: ", output_file)
  
  # Filter significant DEGs
  deg_sig <- deg_results[deg_results$p_val_adj < DEG_PADJ_CUTOFF, ]
  message("  Significant (padj < ", DEG_PADJ_CUTOFF, "): ", nrow(deg_sig))
  
  # Save significant DEGs
  if (nrow(deg_sig) > 0) {
    output_file_sig <- file.path(DEG_DIR, paste0(OUTPUT_PREFIX, "_deg_sig_", treatment))
    save_table_standard(deg_sig, output_file_sig)
  }
  
  # Extract genes by fold change threshold
  for (fc_thresh in c(1, FC_THRESHOLD, 2)) {
    
    # Upregulated genes
    genes_up <- rownames(deg_sig)[deg_sig$avg_log2FC > fc_thresh]
    message("  FC > ", fc_thresh, ": ", length(genes_up), " up")
    
    if (length(genes_up) > 0) {
      gene_file <- file.path(
        TABLES_DIR,
        paste0(OUTPUT_PREFIX, "_", treatment, "_up_fc", fc_thresh, "_genes.txt")
      )
      save_gene_list(genes_up, gene_file)
    }
    
    # Downregulated genes
    genes_down <- rownames(deg_sig)[deg_sig$avg_log2FC < -fc_thresh]
    message("  FC < -", fc_thresh, ": ", length(genes_down), " down")
    
    if (length(genes_down) > 0) {
      gene_file <- file.path(
        TABLES_DIR,
        paste0(OUTPUT_PREFIX, "_", treatment, "_down_fc", fc_thresh, "_genes.txt")
      )
      save_gene_list(genes_down, gene_file)
    }
  }
}

################################################################################
# Generate Summary Statistics
################################################################################

message("\n--- Generating summary statistics ---")

summary_stats <- data.frame(
  treatment = character(),
  total_degs = integer(),
  sig_degs = integer(),
  up_fc1.4 = integer(),
  down_fc1.4 = integer(),
  stringsAsFactors = FALSE
)

for (treatment in treatments) {
  
  deg_file <- file.path(DEG_DIR, paste0(OUTPUT_PREFIX, "_deg_raw_", treatment))
  
  if (!file.exists(deg_file)) next
  
  deg_data <- read.table(deg_file, check.names = FALSE, header = TRUE)
  deg_sig <- deg_data[deg_data$p_val_adj < DEG_PADJ_CUTOFF, ]
  
  n_up <- sum(deg_sig$avg_log2FC > FC_THRESHOLD, na.rm = TRUE)
  n_down <- sum(deg_sig$avg_log2FC < -FC_THRESHOLD, na.rm = TRUE)
  
  summary_stats <- rbind(summary_stats, data.frame(
    treatment = treatment,
    total_degs = nrow(deg_data),
    sig_degs = nrow(deg_sig),
    up_fc1.4 = n_up,
    down_fc1.4 = n_down,
    stringsAsFactors = FALSE
  ))
}

# Save summary
summary_file <- file.path(TABLES_DIR, paste0(OUTPUT_PREFIX, "_deg_summary.txt"))
save_table_standard(summary_stats, summary_file)

message("\n--- DEG Summary ---")
print(summary_stats)

################################################################################
# Visualization
################################################################################

message("\n--- Generating visualizations ---")

# Volcano plots for each treatment
for (treatment in treatments) {
  
  deg_file <- file.path(DEG_DIR, paste0(OUTPUT_PREFIX, "_deg_raw_", treatment))
  
  if (!file.exists(deg_file)) next
  
  deg_data <- read.table(deg_file, check.names = FALSE, header = TRUE)
  
  # Prepare data for volcano plot
  deg_data$significant <- deg_data$p_val_adj < DEG_PADJ_CUTOFF & 
                          abs(deg_data$avg_log2FC) > FC_THRESHOLD
  deg_data$log10padj <- -log10(deg_data$p_val_adj + 1e-300)
  
  # Create volcano plot
  plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_volcano_", treatment, ".pdf"))
  pdf(plot_file, width = 7, height = 6)
  
  plot(
    deg_data$avg_log2FC,
    deg_data$log10padj,
    pch = 20,
    cex = 0.5,
    col = ifelse(deg_data$significant, "red", "gray"),
    xlab = "Log2 Fold Change",
    ylab = "-Log10(Adjusted P-value)",
    main = paste0(treatment, " vs Control")
  )
  
  # Add threshold lines
  abline(h = -log10(DEG_PADJ_CUTOFF), col = "blue", lty = 2)
  abline(v = c(-FC_THRESHOLD, FC_THRESHOLD), col = "blue", lty = 2)
  
  # Add legend
  legend("topright", 
         legend = c("Significant", "Not significant"),
         col = c("red", "gray"),
         pch = 20,
         cex = 0.8)
  
  dev.off()
  message("Volcano plot saved: ", plot_file)
}

# Barplot of DEG counts
plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_deg_counts.pdf"))
pdf(plot_file, width = 8, height = 6)

par(mar = c(8, 5, 4, 2))

# Create matrix for barplot
plot_data <- as.matrix(t(summary_stats[, c("up_fc1.4", "down_fc1.4")]))
colnames(plot_data) <- summary_stats$treatment

barplot(
  plot_data,
  beside = TRUE,
  las = 2,
  col = c("red", "blue"),
  ylab = "Number of DEGs",
  main = paste0("Differentially Expressed Genes (FC > ", FC_THRESHOLD, ")"),
  legend.text = c("Upregulated", "Downregulated"),
  args.legend = list(x = "topright")
)

dev.off()
message("DEG counts barplot saved: ", plot_file)

################################################################################
# Summary Report
################################################################################

message("\n==========================================")
message("Analysis Summary")
message("==========================================")
message("Treatments analyzed: ", length(treatments))
message("Total comparisons: ", nrow(summary_stats))
message("\nDEG Statistics (FC > ", FC_THRESHOLD, "):")

for (i in 1:nrow(summary_stats)) {
  message("  ", summary_stats$treatment[i], ": ",
          summary_stats$up_fc1.4[i], " up, ",
          summary_stats$down_fc1.4[i], " down")
}

message("\nOutput files:")
message("  DEG results: ", DEG_DIR)
message("  Gene lists: ", TABLES_DIR)
message("  Figures: ", FIGURES_DIR)

# Save session info
save_session_info(file.path(RESULTS_DIR, paste0(OUTPUT_PREFIX, "_deg_session_info.txt")))

