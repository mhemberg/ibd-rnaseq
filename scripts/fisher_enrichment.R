################################################################################
# Fisher Exact Test Enrichment Analysis
# IBD Transcriptomics Analysis
#
# Purpose: Perform Fisher exact test for gene set enrichment across cell types
# Input: DEG results from validation dataset, reference IBD DEG data
# Output: Fisher test results and overlapping gene lists
#
################################################################################

# Load configuration and utilities
source("scripts/config.R")
source("scripts/utils_functions.R")

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
})

# Script-specific parameters
OUTPUT_PREFIX <- "fisher_enrichment"
FC_THRESHOLD <- 1.4  # Fold change threshold for defining DEGs

################################################################################
# Load Reference Data
################################################################################

message("==========================================")
message("Fisher Exact Test Enrichment Analysis")
message("==========================================\n")

# Load IBD reference DEG data
message("Loading IBD reference DEG data from: ", IBD_DEG_DATA_PATH)
validate_file_exists(IBD_DEG_DATA_PATH)
ibd_deg <- readRDS(IBD_DEG_DATA_PATH)

# Filter for stroma subset with group annotation
stroma_deg <- ibd_deg[
  ibd_deg$subset == "stroma" & ibd_deg$annotation == "anotation_group",
]

message("Loaded ", nrow(stroma_deg), " stroma DEG records")

################################################################################
# Extract Positive DEGs from Validation Data
################################################################################

message("\n--- Extracting DEGs from validation dataset ---")

# Process each treatment group
for (group_name in BULK_GROUPS) {
  
  message("\nProcessing group: ", group_name)
  
  # Construct path to DEG file
  deg_file <- file.path(
    DEG_DIR,
    paste0("validation_gse233063_deg_raw_", group_name)
  )
  
  # Check if file exists
  if (!file.exists(deg_file)) {
    warning("DEG file not found: ", deg_file)
    warning("Skipping group: ", group_name)
    next
  }
  
  # Read DEG results
  deg_data <- read.table(deg_file, check.names = FALSE, header = TRUE)
  message("  Loaded ", nrow(deg_data), " genes")
  
  # Filter by fold change threshold
  positive_genes <- rownames(deg_data)[deg_data$avg_log2FC > FC_THRESHOLD]
  message("  Genes with FC > ", FC_THRESHOLD, ": ", length(positive_genes))
  
  # Save gene list
  output_file <- file.path(
    TABLES_DIR,
    paste0(OUTPUT_PREFIX, "_", group_name, "_positive_genes.txt")
  )
  save_gene_list(positive_genes, output_file)
}

################################################################################
# Fisher Enrichment Function
################################################################################

#' Perform Fisher Enrichment for Stromal Cell Types
#' 
#' Tests enrichment of query genes in stromal cell type signatures
#' 
#' @param query_genes Character vector of query gene names
#' @param stroma_data Data frame with stromal DEG data
#' @param comparison Comparison to test (e.g., "6_vs_5")
#' @param direction Direction ("UPP" for upregulated, "DWW" for downregulated)
#' @param group_name Name of treatment group (for output naming)
#' @return List with results data frame and S1 overlap genes
perform_fisher_stroma <- function(query_genes,
                                  stroma_data,
                                  comparison,
                                  direction,
                                  group_name) {
  
  # Get unique cell types
  cell_types <- unique(stroma_data$cluster)
  
  # Initialize results vectors
  p_values <- numeric(0)
  q_values <- numeric(0)
  cell_names <- character(0)
  s1_overlap <- "none"
  
  # Test each cell type
  for (cell_type in cell_types) {
    
    message("  Testing cell type: ", cell_type)
    
    # Extract genes for this cell type, comparison, and direction
    cell_genes <- stroma_data$gene[
      stroma_data$cluster == cell_type &
      stroma_data$sign == direction &
      stroma_data$comp == comparison
    ]
    
    # Skip if no genes
    if (length(cell_genes) == 0) {
      message("    No genes found, skipping")
      next
    }
    
    # Perform Fisher exact test
    test_result <- fisher_test_enrichment(
      query_genes = query_genes,
      reference_genes = cell_genes,
      universe_size = FISHER_TOTAL_GENES
    )
    
    # Store overlap for S1 cell type
    if (cell_type == "S1" && test_result$n_overlap > 0) {
      s1_overlap <- test_result$overlap
    }
    
    # Store results
    p_values <- c(p_values, test_result$p_value)
    cell_names <- c(cell_names, cell_type)
    
    message("    Overlap: ", test_result$n_overlap, " genes")
    message("    P-value: ", format(test_result$p_value, scientific = TRUE, digits = 3))
  }
  
  # Calculate adjusted p-values
  q_values <- p.adjust(p_values, method = "BH")
  
  # Create results data frame
  results_df <- data.frame(
    cell_type = cell_names,
    p_value = p_values,
    q_value = q_values
  )
  rownames(results_df) <- cell_names
  
  # Save results
  output_file <- file.path(
    TABLES_DIR,
    paste0(OUTPUT_PREFIX, "_", group_name, "_", direction, "_fisher_", comparison, ".txt")
  )
  save_table_standard(results_df, output_file)
  
  # Save S1 overlap genes
  if (length(s1_overlap) > 0 && s1_overlap[1] != "none") {
    overlap_file <- file.path(
      TABLES_DIR,
      paste0(OUTPUT_PREFIX, "_", group_name, "_", direction, "_s1_overlap_", comparison, ".txt")
    )
    save_gene_list(s1_overlap, overlap_file)
  }
  
  return(list(
    results = results_df,
    s1_overlap = s1_overlap
  ))
}

################################################################################
# Run Fisher Enrichment for All Groups
################################################################################

message("\n--- Running Fisher enrichment analysis ---")

# Comparison to analyze
comparison <- "6_vs_5"  # PFD uninflamed vs CD uninflamed

# Loop through treatment groups
for (group_name in BULK_GROUPS) {
  
  message("\n=== Analyzing group: ", group_name, " ===")
  
  # Load positive genes
  gene_file <- file.path(
    TABLES_DIR,
    paste0(OUTPUT_PREFIX, "_", group_name, "_positive_genes.txt")
  )
  
  if (!file.exists(gene_file)) {
    warning("Gene list file not found: ", gene_file)
    warning("Skipping group: ", group_name)
    next
  }
  
  query_genes <- read_gene_list(gene_file)
  
  # Test upregulated genes (UPP)
  message("\n--- Testing upregulated stroma genes (UPP) ---")
  upp_results <- perform_fisher_stroma(
    query_genes = query_genes,
    stroma_data = stroma_deg,
    comparison = comparison,
    direction = "UPP",
    group_name = group_name
  )
  
  # Test downregulated genes (DWW)
  message("\n--- Testing downregulated stroma genes (DWW) ---")
  dww_results <- perform_fisher_stroma(
    query_genes = query_genes,
    stroma_data = stroma_deg,
    comparison = comparison,
    direction = "DWW",
    group_name = group_name
  )
  
  # Report significant enrichments
  message("\n--- Significant enrichments (q < 0.05) ---")
  
  sig_upp <- upp_results$results[upp_results$results$q_value < 0.05, ]
  if (nrow(sig_upp) > 0) {
    message("Upregulated (UPP):")
    for (i in 1:nrow(sig_upp)) {
      message("  ", rownames(sig_upp)[i], ": q = ", 
              format(sig_upp$q_value[i], scientific = TRUE, digits = 3))
    }
  } else {
    message("Upregulated (UPP): None")
  }
  
  sig_dww <- dww_results$results[dww_results$results$q_value < 0.05, ]
  if (nrow(sig_dww) > 0) {
    message("Downregulated (DWW):")
    for (i in 1:nrow(sig_dww)) {
      message("  ", rownames(sig_dww)[i], ": q = ", 
              format(sig_dww$q_value[i], scientific = TRUE, digits = 3))
    }
  } else {
    message("Downregulated (DWW): None")
  }
}

################################################################################
# Multi-Comparison Analysis (Optional)
################################################################################

message("\n--- Running enrichment for additional comparisons ---")

# Additional comparisons to test
additional_comparisons <- list(
  c(3, 4),  # CD inflamed vs PFD inflamed
  c(4, 5)   # PFD inflamed vs CD uninflamed
)

for (comp in additional_comparisons) {
  
  comp_name <- paste(comp[1], comp[2], sep = "_vs_")
  message("\n=== Comparison: ", comp_name, " ===")
  
  for (group_name in BULK_GROUPS) {
    
    message("\nGroup: ", group_name)
    
    # Load positive genes
    gene_file <- file.path(
      TABLES_DIR,
      paste0(OUTPUT_PREFIX, "_", group_name, "_positive_genes.txt")
    )
    
    if (!file.exists(gene_file)) next
    
    query_genes <- read_gene_list(gene_file)
    
    # Test both directions
    for (direction in c("UPP", "DWW")) {
      
      message("  Direction: ", direction)
      
      results <- perform_fisher_stroma(
        query_genes = query_genes,
        stroma_data = stroma_deg,
        comparison = comp_name,
        direction = direction,
        group_name = paste0(group_name, "_", comp_name)
      )
      
      # Count significant results
      n_sig <- sum(results$results$q_value < 0.05)
      message("    Significant cell types: ", n_sig)
    }
  }
}

################################################################################
# Summary Report
################################################################################

message("\n==========================================")
message("Analysis Summary")
message("==========================================")
message("Treatment groups analyzed: ", paste(BULK_GROUPS, collapse = ", "))
message("Fold change threshold: ", FC_THRESHOLD)
message("Primary comparison: 6_vs_5 (PFD uninfl vs CD uninfl)")
message("Total genes in universe: ", FISHER_TOTAL_GENES)
message("\nOutput files:")
message("  Gene lists: ", TABLES_DIR)
message("  Fisher results: ", TABLES_DIR)

# Save session info
save_session_info(file.path(RESULTS_DIR, paste0(OUTPUT_PREFIX, "_session_info.txt")))

