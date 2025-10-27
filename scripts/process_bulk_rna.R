################################################################################
# Bulk RNA-seq Processing and Analysis
# IBD Transcriptomics Analysis
#
# Purpose: Process bulk RNA-seq data from cytokine-treated organoids
# Input: STAR-aligned BAM/SAM files
# Output: Count matrices, DEG results, module scores
#
# Original file: bulk_seq.docx
#
# Treatments: LTB (Lymphotoxin beta), TGF-β, TNF-α
# Design: Paired samples (treated vs control from same donor)
################################################################################

# Load configuration and utilities
source("scripts/config.R")
source("scripts/utils_functions.R")

# Load required libraries
suppressPackageStartupMessages({
  library(Rsubread)
  library(edgeR)
  library(DESeq2)
  library(Seurat)
})

# Script-specific parameters
OUTPUT_PREFIX <- "bulk_rnaseq"
STAR_OUTPUT_DIR <- file.path(RAW_DATA_DIR, "bulk_rna", "star")
ANNOTATION_GTF <- file.path(GENOME_DIR, "gencode.v41.annotation.gtf")

# Treatment groups
TREATMENT_GROUPS <- c("ltb", "tgfb", "tnf")

# Fold change thresholds for analysis
FC_THRESHOLDS <- c(1, 1.4, 2)

################################################################################
# PART 1: STAR Alignment Instructions
################################################################################

message("==========================================")
message("Bulk RNA-seq Processing")
message("==========================================\n")

message("NOTE: STAR alignment must be run separately on command line")
message("See comments below for STAR commands\n")

# STAR commands for reference (run these in terminal first):
# 
# # 1. Build genome index (one time only)
# STAR --runThreadN 10 \
#      --runMode genomeGenerate \
#      --genomeDir /path/to/STAR_index \
#      --genomeFastaFiles /path/to/GRCh38.p13.genome.fa \
#      --sjdbGTFfile /path/to/gencode.v41.annotation.gtf
# 
# # 2. Align samples (run for each sample)
# STAR --outReadsUnmapped Fastx \
#      --outSAMstrandField intronMotif \
#      --genomeDir /path/to/STAR_index \
#      --readFilesIn sample_R1.fastq.gz sample_R2.fastq.gz \
#      --runThreadN 6 \
#      --readFilesCommand zcat \
#      --outFileNamePrefix sample_
# 
# This generates: sample_Aligned.out.sam

################################################################################
# PART 2: Feature Counting
################################################################################

message("\n--- Feature Counting with featureCounts ---")

# Check if STAR output directory exists
if (!dir.exists(STAR_OUTPUT_DIR)) {
  warning("STAR output directory not found: ", STAR_OUTPUT_DIR)
  message("Please run STAR alignment first or update STAR_OUTPUT_DIR in config.R")
  message("Skipping feature counting...")
} else {
  
  # Get sample list
  sample_list_file <- file.path(RAW_DATA_DIR, "bulk_rna", "sample_list")
  
  if (!file.exists(sample_list_file)) {
    # Create sample list from directory contents
    sample_dirs <- list.dirs(STAR_OUTPUT_DIR, recursive = FALSE, full.names = FALSE)
    sample_list <- sample_dirs
    writeLines(sample_list, sample_list_file)
    message("Created sample list from STAR output directories")
  } else {
    sample_list <- readLines(sample_list_file)
  }
  
  message("Processing ", length(sample_list), " samples")
  
  # Initialize result matrices
  rpkm_res <- NULL
  count_res <- NULL
  
  # Process each sample
  for (j in seq_along(sample_list)) {
    
    sample_id <- sample_list[j]
    message("\nProcessing sample ", j, "/", length(sample_list), ": ", sample_id)
    
    # Construct path to aligned SAM file
    sam_file <- file.path(STAR_OUTPUT_DIR, sample_id, "Aligned.out.sam")
    
    if (!file.exists(sam_file)) {
      warning("SAM file not found: ", sam_file)
      next
    }
    
    # Run featureCounts
    message("  Running featureCounts...")
    count <- featureCounts(
      files = sam_file,
      annot.ext = ANNOTATION_GTF,
      isGTFAnnotationFile = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = "gene_name",
      nthreads = 8,
      isPairedEnd = TRUE
    )
    
    # Extract counts
    sample_counts <- count$counts
    
    # Calculate RPKM
    sample_rpkm <- rpkm(
      count$counts,
      gene.length = count$annotation$Length
    )
    
    # Add to result matrices
    if (is.null(count_res)) {
      count_res <- sample_counts
      rpkm_res <- sample_rpkm
    } else {
      count_res <- cbind(count_res, sample_counts)
      rpkm_res <- cbind(rpkm_res, sample_rpkm)
    }
  }
  
  # Set column names
  colnames(count_res) <- sample_list
  colnames(rpkm_res) <- sample_list
  
  # Fix gene name: CCN2 -> CTGF (legacy naming)
  ccn2_idx <- which(rownames(count_res) == "CCN2")
  if (length(ccn2_idx) > 0) {
    rownames(count_res)[ccn2_idx] <- "CTGF"
    rownames(rpkm_res)[ccn2_idx] <- "CTGF"
    count$annotation$GeneID[ccn2_idx] <- "CTGF"
    message("Renamed CCN2 to CTGF")
  }
  
  # Save annotation
  annot_file <- file.path(PROCESSED_DATA_DIR, paste0(OUTPUT_PREFIX, "_annotation_info.txt"))
  save_table_standard(count$annotation, annot_file)
  
  # Save count matrices
  count_file <- file.path(PROCESSED_DATA_DIR, paste0(OUTPUT_PREFIX, "_counts.txt"))
  save_table_standard(count_res, count_file)
  
  rpkm_file <- file.path(PROCESSED_DATA_DIR, paste0(OUTPUT_PREFIX, "_rpkm.txt"))
  save_table_standard(rpkm_res, rpkm_file)
  
  message("\n✓ Feature counting complete")
}

################################################################################
# PART 3: Filter to Single-Cell Gene Universe
################################################################################

message("\n--- Filtering to single-cell gene universe ---")

# Check if we have count data
count_file <- file.path(PROCESSED_DATA_DIR, paste0(OUTPUT_PREFIX, "_counts.txt"))

if (!file.exists(count_file)) {
  warning("Count file not found. Skipping filtering and downstream analysis.")
  message("Run STAR alignment and feature counting first.")
} else {
  
  # Load count data
  count_res <- read.table(count_file, check.names = FALSE, header = TRUE)
  rpkm_file <- file.path(PROCESSED_DATA_DIR, paste0(OUTPUT_PREFIX, "_rpkm.txt"))
  rpkm_res <- read.table(rpkm_file, check.names = FALSE, header = TRUE)
  
  message("Total genes in bulk RNA-seq: ", nrow(count_res))
  
  # Load single-cell total gene list (if available)
  sc_gene_file <- file.path(RAW_DATA_DIR, "single_cell_total_gn")
  
  if (file.exists(sc_gene_file)) {
    total_gn <- readLines(sc_gene_file)
    message("Single-cell gene universe: ", length(total_gn))
    
    # Find overlap
    gn <- rownames(count_res)
    overlap <- intersect(total_gn, gn)
    message("Overlapping genes: ", length(overlap))
    
    # Filter to overlap
    count_res_filtered <- count_res[overlap, ]
    rpkm_res_filtered <- rpkm_res[overlap, ]
    
  } else {
    message("Single-cell gene list not found. Using all genes.")
    count_res_filtered <- count_res
    rpkm_res_filtered <- rpkm_res
  }
  
  # Remove genes with zero counts across all samples
  gene_sums <- rowSums(count_res_filtered)
  zero_genes <- which(gene_sums == 0)
  
  if (length(zero_genes) > 0) {
    count_res_final <- count_res_filtered[-zero_genes, ]
    rpkm_res_final <- rpkm_res_filtered[-zero_genes, ]
    message("Removed ", length(zero_genes), " genes with zero counts")
  } else {
    count_res_final <- count_res_filtered
    rpkm_res_final <- rpkm_res_filtered
  }
  
  message("Final gene count: ", nrow(count_res_final))
  
  # Save filtered matrices
  count_overlap_file <- file.path(PROCESSED_DATA_DIR, paste0(OUTPUT_PREFIX, "_count_overlap_gn.txt"))
  save_table_standard(count_res_final, count_overlap_file)
  
  rpkm_overlap_file <- file.path(PROCESSED_DATA_DIR, paste0(OUTPUT_PREFIX, "_rpkm_overlap_gn.txt"))
  save_table_standard(rpkm_res_final, rpkm_overlap_file)
}

################################################################################
# PART 4: Differential Expression Analysis with DESeq2
################################################################################

message("\n--- DESeq2 Differential Expression Analysis ---")

#' Run DESeq2 with Paired Design
#' 
#' @param count Count matrix
#' @param treatment_name Treatment name (ltb, tgfb, tnf)
run_deseq2_paired <- function(count, treatment_name) {
  
  message("\nAnalyzing treatment: ", treatment_name)
  
  # Load comparison metadata
  comp_file <- file.path(RAW_DATA_DIR, "bulk_rna", paste0(treatment_name, "_comp"))
  
  if (!file.exists(comp_file)) {
    warning("Comparison file not found: ", comp_file)
    message("Expected format: sample_id, treatment columns")
    message("Creating template file...")
    
    # Create template
    sample_names <- colnames(count)
    template <- data.frame(
      treatment = rep(c("control", "treated"), length.out = length(sample_names)),
      sample_id = rep(1:(length(sample_names)/2), each = 2),
      row.names = sample_names
    )
    write.table(template, comp_file, sep = "\t", quote = FALSE)
    message("Template created. Please edit and re-run.")
    return(NULL)
  }
  
  # Load comparison info
  col_comp <- read.table(comp_file, header = TRUE, row.names = 1)
  col_comp$treatment <- as.factor(col_comp$treatment)
  col_comp$sample_id <- as.factor(col_comp$sample_id)
  
  message("  Samples: ", nrow(col_comp))
  message("  Design: ~sample_id + treatment (paired)")
  
  # Subset count matrix to match samples
  count_subset <- count[, rownames(col_comp)]
  
  # Create DESeq2 object
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = count_subset,
    colData = col_comp,
    design = ~ sample_id + treatment
  )
  
  # Run DESeq2
  message("  Running DESeq2...")
  dds <- DESeq(dds)
  
  # Extract results
  res <- results(dds)
  
  # Save results
  output_file <- file.path(DEG_DIR, paste0(treatment_name, "_deseq2.txt"))
  save_table_standard(res, output_file)
  
  message("  ✓ Results saved: ", output_file)
  
  # Summary statistics
  sig_genes <- sum(res$padj < 0.05, na.rm = TRUE)
  message("  Significant genes (padj < 0.05): ", sig_genes)
  
  return(res)
}

# Load filtered count data
count_overlap_file <- file.path(PROCESSED_DATA_DIR, paste0(OUTPUT_PREFIX, "_count_overlap_gn.txt"))

if (file.exists(count_overlap_file)) {
  count_data <- read.table(count_overlap_file, check.names = FALSE, header = TRUE)
  
  # Run DESeq2 for each treatment
  deseq_results <- list()
  for (treatment in TREATMENT_GROUPS) {
    deseq_results[[treatment]] <- run_deseq2_paired(count_data, treatment)
  }
  
  message("\n✓ DESeq2 analysis complete for all treatments")
} else {
  warning("Filtered count data not found. Skipping DESeq2 analysis.")
}

################################################################################
# PART 5: Extract DEG Gene Lists
################################################################################

message("\n--- Extracting DEG gene lists ---")

for (treatment in TREATMENT_GROUPS) {
  
  message("\nProcessing ", treatment, " DEGs...")
  
  deg_file <- file.path(DEG_DIR, paste0(treatment, "_deseq2.txt"))
  
  if (!file.exists(deg_file)) {
    warning("DEG file not found: ", deg_file)
    next
  }
  
  # Load DEG results
  deg_data <- read.table(deg_file, check.names = FALSE, header = TRUE)
  deg_data <- na.omit(deg_data)
  
  # Filter by significance
  deg_sig <- deg_data[deg_data$padj < 0.05, ]
  
  message("  Total significant genes: ", nrow(deg_sig))
  
  # Extract genes by fold change thresholds
  for (fc_threshold in FC_THRESHOLDS) {
    
    # Upregulated genes
    pos_genes <- rownames(deg_sig)[deg_sig$log2FoldChange > fc_threshold]
    
    if (length(pos_genes) > 0) {
      pos_file <- file.path(
        TABLES_DIR,
        paste0(treatment, "_deg_gn")  # Match original naming
      )
      save_gene_list(pos_genes, pos_file)
    }
    
    # Downregulated genes
    neg_genes <- rownames(deg_sig)[deg_sig$log2FoldChange < -fc_threshold]
    
    if (length(neg_genes) > 0) {
      neg_file <- file.path(
        TABLES_DIR,
        paste0(treatment, "_neg_gn")  # Match original naming
      )
      save_gene_list(neg_genes, neg_file)
    }
    
    message("  FC > ", fc_threshold, ": ", length(pos_genes), " up, ", 
            length(neg_genes), " down")
  }
  
  # Save combined DEG table (up and down at fc 1.4)
  all_genes <- c(
    rownames(deg_sig)[deg_sig$log2FoldChange > 1.4],
    rownames(deg_sig)[deg_sig$log2FoldChange < -1.4]
  )
  
  if (length(all_genes) > 0) {
    deg_table <- deg_data[all_genes, ]
    deg_table_file <- file.path(DEG_DIR, paste0(treatment, "_deg_table.txt"))
    save_table_standard(deg_table, deg_table_file)
  }
}

################################################################################
# PART 6: Module Scoring in Stromal Cells
################################################################################

message("\n--- Module scoring in stromal cells ---")

# Check if stromal data is available
stroma_file <- file.path(RAW_DATA_DIR, "stroma2_badrm_harmony30")

if (!file.exists(stroma_file)) {
  message("Stromal scRNA-seq data not found. Skipping module scoring.")
  message("Expected file: ", stroma_file)
} else {
  
  message("Loading stromal scRNA-seq data...")
  load(stroma_file)
  stroma <- d
  rm(d)
  
  message("Stromal cells: ", ncol(stroma))
  
  # Initialize results matrix
  module_scores <- matrix(nrow = 0, ncol = ncol(stroma))
  colnames(module_scores) <- colnames(stroma)
  
  # Process each treatment
  for (treatment in TREATMENT_GROUPS) {
    
    message("\nModule scoring for ", treatment, "...")
    
    deg_file <- file.path(DEG_DIR, paste0(treatment, "_deseq2.txt"))
    
    if (!file.exists(deg_file)) {
      next
    }
    
    # Load DEG results
    deg_data <- read.table(deg_file, check.names = FALSE, header = TRUE)
    deg_data <- na.omit(deg_data)
    deg_data <- deg_data[deg_data$padj < 0.05, ]
    
    # Test different fold change thresholds
    for (fc_threshold in FC_THRESHOLDS) {
      
      # Get upregulated genes
      pos_genes <- rownames(deg_data)[deg_data$log2FoldChange > fc_threshold]
      
      signature_name <- paste0(treatment, "_", fc_threshold)
      message("  Scoring ", signature_name, " (", length(pos_genes), " genes)")
      
      # Need at least 5 genes for module scoring
      if (length(pos_genes) < 5) {
        message("    Skipping - too few genes")
        next
      }
      
      # Get current number of metadata columns
      n_meta_before <- ncol(stroma@meta.data)
      
      # Add module score
      stroma_scored <- AddModuleScore(
        stroma,
        features = list(pos_genes),
        ctrl = 5,
        name = signature_name
      )
      
      # Extract scores
      n_meta_after <- ncol(stroma_scored@meta.data)
      score_cols <- (n_meta_before + 1):n_meta_after
      
      if (length(score_cols) == 1) {
        scores <- stroma_scored@meta.data[, score_cols]
      } else {
        scores <- apply(stroma_scored@meta.data[, score_cols], 1, mean)
      }
      
      # Add to results
      module_scores <- rbind(module_scores, scores)
      rownames(module_scores)[nrow(module_scores)] <- signature_name
      
      # Update stroma object for next iteration
      stroma <- stroma_scored
    }
  }
  
  # Save module scores
  if (nrow(module_scores) > 0) {
    module_file <- file.path(TABLES_DIR, "stroma_bulkrna_addmodules.txt")
    save_table_standard(module_scores, module_file)
    message("\n✓ Module scores saved: ", module_file)
  }
}

################################################################################
# PART 7: Signature Comparison Across Patient Groups
################################################################################

message("\n--- Comparing signatures across patient groups ---")

if (exists("stroma") && exists("module_scores") && nrow(module_scores) > 0) {
  
  # Create GMT list (gene module scores)
  gmt_list <- list()
  gmt_list[[1]] <- module_scores
  names(gmt_list) <- "cytokine_response"
  
  # Get unique cell types
  celltype_unique <- unique(stroma$anotation_group)
  annotation_level <- "group"
  
  message("Cell types to analyze: ", length(celltype_unique))
  
  # Comparison sets
  compare_sets <- list(
    c(6, 5), c(3, 4), c(4, 5), c(6, 3), c(6, 7), c(7, 5),
    c(1, 2), c(5, 2), c(5, 1)
  )
  
  # Process each comparison
  for (comp_idx in seq_along(compare_sets)) {
    
    comp_set <- compare_sets[[comp_idx]]
    set1 <- comp_set[1]
    set2 <- comp_set[2]
    
    message("\nComparison: Group ", set1, " vs Group ", set2)
    
    # Determine which group variable to use
    if ((comp_set[1] == 1 || comp_set[1] == 2) && 
        (comp_set[2] == 1 || comp_set[2] == 2)) {
      group_order <- stroma$group
    } else {
      group_order <- stroma$new_group
    }
    
    # Process each GMT (signature set)
    for (gmt_idx in seq_along(gmt_list)) {
      
      gmt_name <- names(gmt_list)[gmt_idx]
      gmt_data <- gmt_list[[gmt_idx]]
      
      message("  Signature set: ", gmt_name)
      
      # Initialize results
      results <- matrix(nrow = 0, ncol = 5)
      
      # Test each cell type
      for (celltype in celltype_unique) {
        
        # Get cell indices
        set1_indices <- which(stroma$anotation_group == celltype & group_order == set1)
        set2_indices <- which(stroma$anotation_group == celltype & group_order == set2)
        
        if (length(set1_indices) == 0 || length(set2_indices) == 0) {
          next
        }
        
        # Get main cell type
        maincelltype <- as.character(stroma$maincelltype[set1_indices[1]])
        
        # Test each signature row
        for (sig_idx in 1:nrow(gmt_data)) {
          
          sig_scores <- gmt_data[sig_idx, ]
          
          # Extract scores
          group1_scores <- sig_scores[set1_indices]
          group2_scores <- sig_scores[set2_indices]
          
          # Calculate means
          mean1 <- if (length(set1_indices) == 1) group1_scores else mean(group1_scores)
          mean2 <- if (length(set2_indices) == 1) group2_scores else mean(group2_scores)
          
          # Wilcoxon test
          test_result <- wilcox.test(group1_scores, group2_scores, alternative = "two.sided")
          
          # Store results
          result_row <- c(test_result$p.value, mean1, mean2, celltype, maincelltype)
          results <- rbind(results, result_row)
        }
      }
      
      if (nrow(results) == 0) {
        message("    No valid comparisons")
        next
      }
      
      # Format results
      results_df <- as.data.frame(results, stringsAsFactors = FALSE)
      colnames(results_df) <- c("pval", set1, set2, "cell", "maincelltype")
      
      # Convert numeric columns
      results_df$pval <- as.numeric(results_df$pval)
      results_df[[as.character(set1)]] <- as.numeric(results_df[[as.character(set1)]])
      results_df[[as.character(set2)]] <- as.numeric(results_df[[as.character(set2)]])
      
      # Save results
      output_file <- file.path(
        TABLES_DIR,
        paste0("signature_comp_", gmt_name, "_", annotation_level, "_", set1, "_", set2, ".txt")
      )
      save_table_standard(results_df, output_file)
      
      # Report significant results
      sig_count <- sum(results_df$pval < 0.05)
      message("    Significant (p < 0.05): ", sig_count)
    }
  }
  
  message("\n✓ Signature comparison complete")
}

################################################################################
# Summary Report
################################################################################

message("\n==========================================")
message("Analysis Summary")
message("==========================================")

message("\nTreatment groups analyzed: ", paste(TREATMENT_GROUPS, collapse = ", "))

for (treatment in TREATMENT_GROUPS) {
  deg_file <- file.path(DEG_DIR, paste0(treatment, "_deseq2.txt"))
  if (file.exists(deg_file)) {
    deg <- read.table(deg_file, check.names = FALSE, header = TRUE)
    deg_sig <- na.omit(deg)
    deg_sig <- deg_sig[deg_sig$padj < 0.05, ]
    n_up <- sum(deg_sig$log2FoldChange > 1.4)
    n_down <- sum(deg_sig$log2FoldChange < -1.4)
    message("  ", treatment, ": ", n_up, " up, ", n_down, " down (FC > 1.4)")
  }
}

message("\nOutput files:")
message("  Processed counts: ", PROCESSED_DATA_DIR)
message("  DEG results: ", DEG_DIR)
message("  Gene lists: ", TABLES_DIR)
message("  Module scores: ", TABLES_DIR)

# Save session info
save_session_info(file.path(RESULTS_DIR, paste0(OUTPUT_PREFIX, "_session_info.txt")))

message("\n==========================================")
message("Bulk RNA-seq processing complete!")
message("==========================================")
