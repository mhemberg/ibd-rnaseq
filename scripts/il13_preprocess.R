################################################################################
# IL-13 Microarray Data Preprocessing
# IBD Transcriptomics Analysis
#
# Purpose: Download, preprocess, and analyze IL-13 treatment microarray data
# Input: GSE190705 data from GEO
# Output: Processed count matrix and DESeq2 results
################################################################################

# Load configuration and utilities
source("scripts/config.R")
source("scripts/utils_functions.R")

# Load required libraries
suppressPackageStartupMessages({
  library(biomaRt)
  library(DESeq2)
})

# Script-specific parameters
OUTPUT_PREFIX <- "il13_microarray"
GEO_ACCESSION <- "GSE190705"
DATA_URL <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190705/suppl/GSE190705_read_counts.txt.gz"

################################################################################
# Download Data
################################################################################

message("==========================================")
message("IL-13 Microarray Data Preprocessing")
message("==========================================\n")

# Set download path
download_file <- file.path(RAW_DATA_DIR, "GSE190705_read_counts.txt.gz")
extracted_file <- file.path(RAW_DATA_DIR, "GSE190705_read_counts.txt")

# Download if not already present
if (!file.exists(extracted_file)) {
  message("Downloading data from GEO...")
  message("URL: ", DATA_URL)
  message("Destination: ", download_file)
  
  download.file(DATA_URL, download_file, method = "auto")
  message("Download complete")
  
  # Extract gz file
  message("Extracting file...")
  system(paste("gunzip", download_file))
  message("Extraction complete")
} else {
  message("Data file already exists: ", extracted_file)
}

################################################################################
# Load and Parse Data
################################################################################

message("\n--- Loading count data ---")

# Read count data
count_data <- read.table(extracted_file, check.names = FALSE, header = TRUE)
message("Dimensions: ", nrow(count_data), " genes x ", ncol(count_data), " columns")

# Extract Ensembl IDs
ensembl_ids <- count_data$geneID
message("Total genes: ", length(ensembl_ids))

################################################################################
# Convert Ensembl IDs to Gene Symbols
################################################################################

message("\n--- Converting Ensembl IDs to gene symbols ---")

# Use biomaRt for conversion
gene_info <- convert_ensembl_to_symbol(ensembl_ids)

message("Successfully mapped genes: ", nrow(gene_info))

################################################################################
# Process Count Matrix
################################################################################

message("\n--- Processing count matrix ---")

# Find overlap between mapped genes and count data
overlap_ids <- intersect(gene_info$ensembl_gene_id, count_data$geneID)
message("Overlapping genes: ", length(overlap_ids))

# Subset to overlapping genes
rownames(gene_info) <- gene_info$ensembl_gene_id
rownames(count_data) <- count_data$geneID

count_overlap <- count_data[overlap_ids, ]
gene_info_overlap <- gene_info[overlap_ids, ]

# Replace Ensembl IDs with gene symbols
count_overlap$geneID <- gene_info_overlap$external_gene_name

# Remove duplicates (keep first occurrence)
dup_mask <- duplicated(count_overlap$geneID)
n_duplicates <- sum(dup_mask)
message("Duplicate genes found: ", n_duplicates)

count_unique <- count_overlap[!dup_mask, ]
message("Unique genes: ", nrow(count_unique))

# Select sample columns (based on GSE190705 metadata)
# Columns 2,3,4 are control samples; 8,9,10 are IL-13 treated samples
sample_columns <- c(2, 3, 4, 8, 9, 10)
count_matrix <- count_unique[, sample_columns]
rownames(count_matrix) <- count_unique$geneID

message("Final count matrix: ", nrow(count_matrix), " genes x ", ncol(count_matrix), " samples")

# Remove genes with zero counts across all samples
gene_sums <- rowSums(count_matrix)
count_filtered <- count_matrix[gene_sums != 0, ]
message("After filtering zero-expression genes: ", nrow(count_filtered), " genes")

# Save parsed table
output_file <- file.path(PROCESSED_DATA_DIR, paste0(OUTPUT_PREFIX, "_parsed_counts.txt"))
save_table_standard(count_filtered, output_file)

################################################################################
# Create Sample Annotation
################################################################################

message("\n--- Creating sample annotation ---")

# Create annotation data frame (adjust based on actual GSE190705 metadata)
# Sample IDs from column names
sample_ids <- colnames(count_filtered)

# Create treatment and sample ID annotations
# Assuming first 3 are control, last 3 are IL-13 treated
# And they represent paired samples (1, 2, 3)
annotation <- data.frame(
  sample_id = factor(c(1, 2, 3, 1, 2, 3)),
  treatment = factor(c("control", "control", "control", "IL13", "IL13", "IL13"))
)
rownames(annotation) <- sample_ids

message("Sample annotation:")
print(annotation)

# Save annotation
annot_file <- file.path(PROCESSED_DATA_DIR, paste0(OUTPUT_PREFIX, "_annotation.txt"))
save_table_standard(annotation, annot_file)

################################################################################
# DESeq2 Analysis
################################################################################

message("\n--- Running DESeq2 analysis ---")

# Create DESeq2 object with paired design
# Design: ~sample_id + treatment
# This accounts for paired samples (treatment is the variable of interest)
dds <- DESeqDataSetFromMatrix(
  countData = count_filtered,
  colData = annotation,
  design = ~ sample_id + treatment
)

message("DESeq2 design: ~ sample_id + treatment")
message("  This design accounts for paired samples")
message("  Treatment effect is tested after controlling for sample_id")

# Run DESeq2
message("\nRunning differential expression analysis...")
dds <- DESeq(dds)

# Extract results
results_deseq <- results(dds)
message("DESeq2 analysis complete")

# Summary of results
message("\nResults summary:")
print(summary(results_deseq))

# Save DESeq2 results
deseq_output <- file.path(DEG_DIR, paste0(OUTPUT_PREFIX, "_deseq2_results.txt"))
save_table_standard(results_deseq, deseq_output)

################################################################################
# Filter and Export Significant DEGs
################################################################################

message("\n--- Filtering significant DEGs ---")

# Convert to data frame and remove NAs
results_df <- as.data.frame(results_deseq)
results_clean <- na.omit(results_df)

# Filter by significance
sig_genes <- results_clean[results_clean$padj < DEG_PADJ_CUTOFF, ]
message("Genes with padj < ", DEG_PADJ_CUTOFF, ": ", nrow(sig_genes))

# Filter by fold change thresholds
for (fc_threshold in DEG_FC_THRESHOLDS) {
  
  # Upregulated genes
  up_genes <- sig_genes[sig_genes$log2FoldChange > fc_threshold, ]
  message("  Upregulated (FC > ", fc_threshold, "): ", nrow(up_genes))
  
  # Save upregulated genes
  if (nrow(up_genes) > 0) {
    up_file <- file.path(
      TABLES_DIR,
      paste0(OUTPUT_PREFIX, "_upregulated_fc", fc_threshold, ".txt")
    )
    save_table_standard(up_genes, up_file)
    
    # Save gene list
    up_genelist_file <- file.path(
      TABLES_DIR,
      paste0(OUTPUT_PREFIX, "_upregulated_fc", fc_threshold, "_genes.txt")
    )
    save_gene_list(rownames(up_genes), up_genelist_file)
  }
  
  # Downregulated genes
  down_genes <- sig_genes[sig_genes$log2FoldChange < -fc_threshold, ]
  message("  Downregulated (FC < -", fc_threshold, "): ", nrow(down_genes))
  
  # Save downregulated genes
  if (nrow(down_genes) > 0) {
    down_file <- file.path(
      TABLES_DIR,
      paste0(OUTPUT_PREFIX, "_downregulated_fc", fc_threshold, ".txt")
    )
    save_table_standard(down_genes, down_file)
    
    # Save gene list
    down_genelist_file <- file.path(
      TABLES_DIR,
      paste0(OUTPUT_PREFIX, "_downregulated_fc", fc_threshold, "_genes.txt")
    )
    save_gene_list(rownames(down_genes), down_genelist_file)
  }
}

################################################################################
# Generate Visualizations
################################################################################

message("\n--- Generating visualizations ---")

# MA plot
plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_ma_plot.pdf"))
pdf(plot_file, width = 7, height = 5)

plotMA(results_deseq, main = "IL-13 Treatment MA Plot", ylim = c(-5, 5))

dev.off()
message("MA plot saved to: ", plot_file)

# Volcano plot
plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_volcano_plot.pdf"))
pdf(plot_file, width = 7, height = 5)

# Prepare data for volcano plot
plot_data <- results_clean
plot_data$significant <- plot_data$padj < DEG_PADJ_CUTOFF
plot_data$log10padj <- -log10(plot_data$padj)

# Create volcano plot
plot(
  plot_data$log2FoldChange,
  plot_data$log10padj,
  pch = 20,
  cex = 0.5,
  col = ifelse(plot_data$significant, "red", "gray"),
  xlab = "Log2 Fold Change",
  ylab = "-Log10(Adjusted P-value)",
  main = "IL-13 Treatment Volcano Plot"
)

# Add threshold lines
abline(h = -log10(DEG_PADJ_CUTOFF), col = "blue", lty = 2)
abline(v = c(-1, 1), col = "blue", lty = 2)

dev.off()
message("Volcano plot saved to: ", plot_file)

# Top genes table
n_top_genes <- 50
top_genes <- head(sig_genes[order(sig_genes$padj), ], n_top_genes)

message("\nTop ", n_top_genes, " significant genes (by adjusted p-value):")
print(head(top_genes[, c("log2FoldChange", "padj")], 10))

# Save top genes
top_genes_file <- file.path(TABLES_DIR, paste0(OUTPUT_PREFIX, "_top", n_top_genes, "_genes.txt"))
save_table_standard(top_genes, top_genes_file)

################################################################################
# Expression Heatmap Data Preparation
################################################################################

message("\n--- Preparing heatmap data ---")

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Select top variable genes for heatmap
top_var_genes <- head(order(apply(normalized_counts, 1, var), decreasing = TRUE), 100)
heatmap_data <- normalized_counts[top_var_genes, ]

# Log transform for visualization
heatmap_data_log <- log2(heatmap_data + 1)

# Save heatmap data
heatmap_file <- file.path(TABLES_DIR, paste0(OUTPUT_PREFIX, "_heatmap_data.txt"))
save_table_standard(heatmap_data_log, heatmap_file)

message("Heatmap data saved (top 100 variable genes)")

################################################################################
# Quality Control Plots
################################################################################

message("\n--- Generating QC plots ---")

# Sample distance heatmap
plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_sample_distances.pdf"))
pdf(plot_file, width = 6, height = 5)

sample_dists <- dist(t(assay(dds)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- paste(
  annotation$treatment,
  annotation$sample_id,
  sep = "_"
)
colnames(sample_dist_matrix) <- rownames(sample_dist_matrix)

# Use pheatmap if available
if (require(pheatmap, quietly = TRUE)) {
  pheatmap(
    sample_dist_matrix,
    clustering_distance_rows = sample_dists,
    clustering_distance_cols = sample_dists,
    main = "Sample-to-Sample Distances"
  )
} else {
  heatmap(sample_dist_matrix, main = "Sample-to-Sample Distances")
}

dev.off()
message("Sample distance heatmap saved to: ", plot_file)

# PCA plot
plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_pca.pdf"))
pdf(plot_file, width = 7, height = 5)

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# PCA plot
plotPCA(vsd, intgroup = "treatment") +
  ggtitle("PCA - IL-13 Treatment")

dev.off()
message("PCA plot saved to: ", plot_file)

# Dispersion plot
plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_dispersion.pdf"))
pdf(plot_file, width = 7, height = 5)

plotDispEsts(dds, main = "Dispersion Estimates")

dev.off()
message("Dispersion plot saved to: ", plot_file)

################################################################################
# Summary Report
################################################################################

message("\n==========================================")
message("Analysis Summary")
message("==========================================")
message("Dataset: ", GEO_ACCESSION)
message("Total genes analyzed: ", nrow(count_filtered))
message("Samples: ", ncol(count_filtered))
message("Design: Paired (control vs IL-13 treatment)")
message("\nSignificant DEGs (padj < ", DEG_PADJ_CUTOFF, "): ", nrow(sig_genes))

# Report by fold change threshold
for (fc_threshold in DEG_FC_THRESHOLDS) {
  n_up <- sum(sig_genes$log2FoldChange > fc_threshold)
  n_down <- sum(sig_genes$log2FoldChange < -fc_threshold)
  message("  FC > ", fc_threshold, ": ", n_up, " up, ", n_down, " down")
}

message("\nOutput files:")
message("  Processed counts: ", PROCESSED_DATA_DIR)
message("  DEG results: ", DEG_DIR)
message("  Gene lists: ", TABLES_DIR)
message("  Figures: ", FIGURES_DIR)

# Top upregulated genes
top_up <- head(sig_genes[order(sig_genes$log2FoldChange, decreasing = TRUE), ], 5)
message("\nTop 5 upregulated genes:")
for (i in 1:min(5, nrow(top_up))) {
  gene <- rownames(top_up)[i]
  fc <- top_up$log2FoldChange[i]
  padj <- top_up$padj[i]
  message("  ", gene, ": log2FC = ", round(fc, 2), 
          ", padj = ", format(padj, scientific = TRUE, digits = 3))
}

# Top downregulated genes
top_down <- head(sig_genes[order(sig_genes$log2FoldChange, decreasing = FALSE), ], 5)
message("\nTop 5 downregulated genes:")
for (i in 1:min(5, nrow(top_down))) {
  gene <- rownames(top_down)[i]
  fc <- top_down$log2FoldChange[i]
  padj <- top_down$padj[i]
  message("  ", gene, ": log2FC = ", round(fc, 2), 
          ", padj = ", format(padj, scientific = TRUE, digits = 3))
}

# Save session info
save_session_info(file.path(RESULTS_DIR, paste0(OUTPUT_PREFIX, "_session_info.txt")))

