################################################################################
# Pathway Module Scoring with GO/KEGG
# IBD Transcriptomics Analysis
#
# Purpose: Calculate module scores for GO/KEGG pathways across cell types
# Input: Seurat object, GO/KEGG gene sets
# Output: Module scores per cell, pathway comparisons across groups
################################################################################

# Load configuration and utilities
source("scripts/config.R")
source("scripts/utils_functions.R")

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(GSEABase)
  library(dplyr)
})

# Script-specific parameters
OUTPUT_PREFIX <- "pathway_modules"
GMT_DIR <- file.path(REF_DIR, "genesets")

################################################################################
# Custom Functions
################################################################################

#' Compare Module Scores Across Patient Groups
#'
#' @param seurat_obj Seurat object
#' @param module_scores Matrix of module scores (pathways x cells)
#' @param compare_set Comparison (e.g., c(6, 5))
#' @param annotation_level "group" or "intermediate"
#' @return Results data frame
compare_pathway_scores <- function(seurat_obj, module_scores, compare_set, annotation_level) {
  
  # Determine group variable
  if ((compare_set[1] %in% c(1, 2)) && (compare_set[2] %in% c(1, 2))) {
    group_order <- seurat_obj$group
  } else {
    group_order <- seurat_obj$new_group
  }
  
  set1 <- compare_set[1]
  set2 <- compare_set[2]
  
  # Get cell types
  if (annotation_level == "group") {
    celltype_unique <- unique(seurat_obj$anotation_group)
  } else {
    celltype_unique <- unique(seurat_obj$anotation_intermediate)
  }
  
  results <- data.frame()
  
  # Test each cell type
  for (celltype in celltype_unique) {
    
    # Get cell indices
    if (annotation_level == "group") {
      set1_idx <- which(seurat_obj$anotation_group == celltype & group_order == set1)
      set2_idx <- which(seurat_obj$anotation_group == celltype & group_order == set2)
    } else {
      set1_idx <- which(seurat_obj$anotation_intermediate == celltype & group_order == set1)
      set2_idx <- which(seurat_obj$anotation_intermediate == celltype & group_order == set2)
    }
    
    if (length(set1_idx) == 0 || length(set2_idx) == 0) next
    
    # Get main cell type
    maincelltype <- as.character(seurat_obj$maincelltype[set1_idx[1]])
    
    # Test each pathway
    for (pathway_idx in 1:nrow(module_scores)) {
      
      pathway <- rownames(module_scores)[pathway_idx]
      
      # Extract scores
      scores1 <- module_scores[pathway_idx, set1_idx]
      scores2 <- module_scores[pathway_idx, set2_idx]
      
      # Skip if all negative (artifacts)
      mean1 <- mean(scores1)
      mean2 <- mean(scores2)
      
      # Wilcoxon test
      test_res <- wilcox.test(scores1, scores2, alternative = "two.sided")
      
      results <- rbind(results, data.frame(
        pathway = pathway,
        pval = test_res$p.value,
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
  
  # Filter significant only
  results_sig <- results[results$pval < 0.05, ]
  
  # Reorder columns
  colnames(results_sig)[3:4] <- c(as.character(set1), as.character(set2))
  
  return(results_sig)
}

################################################################################
# Main Analysis
################################################################################

message("==========================================")
message("Pathway Module Scoring Analysis")
message("==========================================\n")

# Create GMT directory if needed
if (!dir.exists(GMT_DIR)) {
  dir.create(GMT_DIR, recursive = TRUE)
}

################################################################################
# Download/Load Gene Sets
################################################################################

message("--- Loading gene sets ---")

# GO Component
go_component_file <- file.path(GMT_DIR, "gnsymbol2go.human.Component")
if (!file.exists(go_component_file)) {
  message("GO Component GMT not found.")
  go_component_file <- NULL
}

# GO Function
go_function_file <- file.path(GMT_DIR, "gnsymbol2go.human.Function")
if (!file.exists(go_function_file)) {
  go_function_file <- file.path(REF_DIR, "gnsymbol2go.human.Function")
  if (!file.exists(go_function_file)) {
    message("GO Function gene set not found, will skip")
    go_function_file <- NULL
  }
}

# KEGG
kegg_file <- file.path(GMT_DIR, "kegg.human_gmt")
if (!file.exists(kegg_file)) {
  message("KEGG GMT not found. Will skip KEGG analysis")
  kegg_file <- NULL
}

# Reactome
reactome_file <- file.path(GMT_DIR, "ReactomePathways.gmt")
if (!file.exists(reactome_file)) {
  message("Reactome GMT not found. Will skip Reactome analysis")
  reactome_file <- NULL
}

# Load gene sets
gmt_list <- list()

if (file.exists(go_component_file)) {
  message("Loading GO Component...")
  gmt1 <- getGmt(go_component_file)
  gmt_list[["go_component"]] <- gmt1
}

if (!is.null(go_function_file) && file.exists(go_function_file)) {
  message("Loading GO Function...")
  gmt2 <- getGmt(go_function_file)
  gmt_list[["go_function"]] <- gmt2
}

if (!is.null(kegg_file) && file.exists(kegg_file)) {
  message("Loading KEGG...")
  gmt3 <- getGmt(kegg_file)
  gmt_list[["kegg"]] <- gmt3
}

if (!is.null(reactome_file) && file.exists(reactome_file)) {
  message("Loading Reactome...")
  gmt4 <- getGmt(reactome_file)
  gmt_list[["reactome"]] <- gmt4
}

if (length(gmt_list) == 0) {
  stop("No gene sets loaded. Please provide GMT files.")
}

message("Loaded ", length(gmt_list), " gene set collections")

################################################################################
# Load Seurat Object
################################################################################

message("\n--- Loading Seurat object ---")
total_file <- file.path(RAW_DATA_DIR, "total_no_umap.rdata")
validate_file_exists(total_file)

total3 <- get(load(total_file))

message("Total cells: ", ncol(total3))

# Get gene universe
total_genes <- rownames(total3)
message("Total genes in dataset: ", length(total_genes))

################################################################################
# Calculate Module Scores
################################################################################

message("\n--- Calculating module scores ---")

# Store all module scores
all_module_scores <- list()

for (gmt_name in names(gmt_list)) {
  
  message("\n=== Processing: ", gmt_name, " ===")
  
  gmt <- gmt_list[[gmt_name]]
  n_pathways <- length(gmt@.Data)
  
  message("Number of pathways: ", n_pathways)
  
  # Initialize result matrix
  module_scores <- matrix(nrow = 0, ncol = ncol(total3))
  colnames(module_scores) <- colnames(total3)
  
  # Process each pathway
  for (i in 1:n_pathways) {
    
    pathway_genes <- gmt@.Data[[i]]@geneIds
    pathway_name <- gmt@.Data[[i]]@setName
    
    if (i %% 100 == 0) {
      message("  Processing pathway ", i, "/", n_pathways)
    }
    
    # Filter to genes in dataset
    pathway_genes_filtered <- intersect(pathway_genes, total_genes)
    
    # Need at least 5 genes
    if (length(pathway_genes_filtered) < 5) next
    
    # Add module score
    temp_obj <- AddModuleScore(
      total3,
      features = list(pathway_genes_filtered),
      ctrl = 5,
      name = "temp_score"
    )
    
    # Extract scores
    n_meta_before <- ncol(total3@meta.data)
    n_meta_after <- ncol(temp_obj@meta.data)
    score_cols <- (n_meta_before + 1):n_meta_after
    
    if (length(score_cols) == 1) {
      scores <- temp_obj@meta.data[, score_cols]
    } else {
      scores <- apply(temp_obj@meta.data[, score_cols], 1, mean)
    }
    
    # Add to result
    module_scores <- rbind(module_scores, scores)
    rownames(module_scores)[nrow(module_scores)] <- pathway_name
  }
  
  message("Calculated scores for ", nrow(module_scores), " pathways")
  
  # Save module scores
  score_file <- file.path(TABLES_DIR, paste0(gmt_name, "_addmodules.txt"))
  save_table_standard(module_scores, score_file)
  
  # Store in list
  all_module_scores[[gmt_name]] <- module_scores
}

################################################################################
# Compare Across Patient Groups
################################################################################

message("\n--- Comparing pathway scores across patient groups ---")

# Comparisons
compare_sets <- list(c(6,5), c(3,4), c(4,5), c(6,3), c(6,7), c(7,5), 
                    c(1,2), c(5,2), c(5,1))

# Annotation levels
annotation_levels <- c("group", "intermediate")

for (annot_level in annotation_levels) {
  
  message("\n=== Annotation level: ", annot_level, " ===")
  
  for (compare_set in compare_sets) {
    
    set1 <- compare_set[1]
    set2 <- compare_set[2]
    
    message("\n--- Comparison: ", set1, " vs ", set2, " ---")
    
    # Process each gene set collection
    for (gmt_name in names(all_module_scores)) {
      
      message("  Gene set: ", gmt_name)
      
      module_scores <- all_module_scores[[gmt_name]]
      
      # Compare
      results <- tryCatch({
        compare_pathway_scores(total3, module_scores, compare_set, annot_level)
      }, error = function(e) {
        message("    Error: ", e$message)
        return(NULL)
      })
      
      if (is.null(results) || nrow(results) == 0) {
        message("    No significant results")
        next
      }
      
      # Order by p-value
      results <- results[order(results$pval), ]
      
      # Save
      output_file <- file.path(
        TABLES_DIR,
        paste0("signature_comp_", gmt_name, "_", annot_level, "_", set1, "_", set2, ".txt")
      )
      save_table_standard(results, output_file)
      
      message("    Significant pathways: ", nrow(results))
    }
  }
}

################################################################################
# Generate Top Pathway Summary
################################################################################

message("\n--- Generating summary of top pathways ---")

summary_list <- list()

for (gmt_name in names(all_module_scores)) {
  
  # Look at main comparison (6 vs 5)
  comp_file <- file.path(
    TABLES_DIR,
    paste0("signature_comp_", gmt_name, "_group_6_5.txt")
  )
  
  if (file.exists(comp_file)) {
    comp_data <- read.table(comp_file, check.names = FALSE, header = TRUE, sep = '\t')
    
    # Top 10 pathways
    top_pathways <- head(comp_data[order(comp_data$pval), ], 10)
    
    summary_list[[gmt_name]] <- top_pathways
    
    message("\nTop 10 pathways for ", gmt_name, " (6 vs 5):")
    for (i in 1:min(10, nrow(top_pathways))) {
      message("  ", i, ". ", top_pathways$pathway[i], 
              " (q = ", format(top_pathways$qval[i], scientific = TRUE, digits = 3), ")")
    }
  }
}

################################################################################
# Summary Report
################################################################################

message("\n==========================================")
message("Analysis Summary")
message("==========================================")
message("Gene set collections analyzed: ", paste(names(gmt_list), collapse = ", "))
message("Cells analyzed: ", ncol(total3))
message("Comparisons performed: ", length(compare_sets))
message("\nOutput files:")
message("  Module scores: ", TABLES_DIR, "/*_addmodules.txt")
message("  Comparison results: ", TABLES_DIR, "/signature_comp_*.txt")

# Save session info
save_session_info(file.path(RESULTS_DIR, paste0(OUTPUT_PREFIX, "_session_info.txt")))

