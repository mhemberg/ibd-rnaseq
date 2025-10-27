################################################################################
# Common Utility Functions
# IBD Transcriptomics Analysis
#
# Purpose: Shared functions used across multiple analysis scripts
# Usage: source("scripts/utils/common_functions.R")
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(edgeR)
})

# Pseudobulk Analysis ----

#' Calculate Pseudobulk Expression
#' 
#' Aggregates single-cell expression data by cell type
#' 
#' @param expression_matrix Expression matrix (genes x cells)
#' @param celltype_labels Vector of cell type labels for each cell
#' @return Matrix of pseudobulk expression (genes x cell types)
#' 
#' @examples
#' pseudobulk <- get_pseudobulk(seurat_obj[["RNA"]]@data, seurat_obj$celltype)
get_pseudobulk <- function(expression_matrix, celltype_labels) {
  
  # Get unique cell types
  cell_types <- unique(celltype_labels)
  
  # Initialize result matrix
  result <- matrix(nrow = nrow(expression_matrix), ncol = length(cell_types))
  colnames(result) <- cell_types
  rownames(result) <- rownames(expression_matrix)
  
  # Calculate mean expression for each cell type
  for (i in seq_along(cell_types)) {
    cell_type <- cell_types[i]
    cell_indices <- which(celltype_labels == cell_type)
    
    # Handle single cell case
    if (length(cell_indices) == 1) {
      result[, i] <- expression_matrix[, cell_indices]
    } else {
      result[, i] <- rowMeans(expression_matrix[, cell_indices, drop = FALSE])
    }
  }
  
  return(result)
}


# Statistical Tests ----

#' Fisher Exact Test for Gene Set Enrichment
#' 
#' Performs Fisher exact test to assess enrichment of query genes in reference set
#' 
#' @param query_genes Character vector of query gene names
#' @param reference_genes Character vector of reference gene names
#' @param universe_size Total number of genes in background (default: 33538)
#' @param alternative Direction of test ("greater", "less", "two.sided")
#' @return List containing p-value, overlap genes, and contingency table
#' 
#' @examples
#' result <- fisher_test_enrichment(deg_genes, pathway_genes)
fisher_test_enrichment <- function(query_genes, 
                                   reference_genes, 
                                   universe_size = 33538,
                                   alternative = "greater") {
  
  # Find overlap
  overlap_genes <- intersect(query_genes, reference_genes)
  
  # Construct 2x2 contingency table
  # | In Reference | Not in Reference |
  # |--------------+------------------|
  # | In Query     | m1               | m3               |
  # | Not in Query | m2               | m4               |
  
  m1 <- length(overlap_genes)                      # In both
  m2 <- length(reference_genes) - m1               # In reference only
  m3 <- length(query_genes) - m1                   # In query only
  m4 <- universe_size - m1 - m2 - m3              # In neither
  
  contingency_matrix <- matrix(c(m1, m2, m3, m4), nrow = 2, ncol = 2)
  
  # Perform Fisher exact test
  test_result <- fisher.test(contingency_matrix, alternative = alternative)
  
  return(list(
    p_value = test_result$p.value,
    overlap = overlap_genes,
    n_overlap = m1,
    n_query = length(query_genes),
    n_reference = length(reference_genes),
    contingency_matrix = contingency_matrix,
    odds_ratio = test_result$estimate
  ))
}


#' Batch Fisher Enrichment Analysis
#' 
#' Performs Fisher exact test across multiple cell types/groups
#' 
#' @param query_genes Character vector of query genes
#' @param reference_data Data frame with columns: gene, cluster, sign, comp
#' @param comparison Comparison name to filter (e.g., "6_vs_5")
#' @param direction Direction to test ("UPP" or "DWW")
#' @param universe_size Total genes in background
#' @return Data frame with p-values and adjusted p-values per cell type
fisher_enrichment_batch <- function(query_genes,
                                   reference_data,
                                   comparison,
                                   direction = "UPP",
                                   universe_size = 33538) {
  
  # Filter reference data
  ref_filtered <- reference_data[
    reference_data$comp == comparison & reference_data$sign == direction,
  ]
  
  # Get unique cell types
  cell_types <- unique(ref_filtered$cluster)
  
  # Initialize results
  p_values <- numeric(0)
  q_values <- numeric(0)
  cell_names <- character(0)
  overlap_genes <- list()
  
  # Test each cell type
  for (cell in cell_types) {
    cell_genes <- ref_filtered$gene[ref_filtered$cluster == cell]
    
    # Skip if no genes
    if (length(cell_genes) == 0) next
    
    # Perform test
    test_result <- fisher_test_enrichment(
      query_genes = query_genes,
      reference_genes = cell_genes,
      universe_size = universe_size
    )
    
    # Store results
    p_values <- c(p_values, test_result$p_value)
    cell_names <- c(cell_names, cell)
    overlap_genes[[cell]] <- test_result$overlap
  }
  
  # Adjust p-values
  q_values <- p.adjust(p_values, method = "BH")
  
  # Combine results
  results <- data.frame(
    cell_type = cell_names,
    p_value = p_values,
    q_value = q_values,
    row.names = cell_names
  )
  
  return(list(
    results = results,
    overlap_genes = overlap_genes
  ))
}


# DESeq2 Helpers ----

#' Run DESeq2 with Paired Design
#' 
#' Performs differential expression analysis using DESeq2 with paired samples
#' 
#' @param count_matrix Count matrix (genes x samples)
#' @param metadata Data frame with sample information
#' @param design_formula Formula for DESeq2 design (default: ~sample_id + treatment)
#' @return DESeq2 results object
#' 
#' @examples
#' results <- run_deseq2_paired(counts, metadata)
run_deseq2_paired <- function(count_matrix,
                              metadata,
                              design_formula = ~ sample_id + treatment) {
  
  # Ensure factors
  metadata$treatment <- as.factor(metadata$treatment)
  metadata$sample_id <- as.factor(metadata$sample_id)
  
  # Reorder count matrix to match metadata
  count_matrix <- count_matrix[, rownames(metadata)]
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = metadata,
    design = design_formula
  )
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Extract results
  res <- results(dds)
  
  return(res)
}


#' Filter DEGs by Thresholds
#' 
#' Filters differential expression results by fold change and p-value
#' 
#' @param deg_results DESeq2 results object or data frame
#' @param log2fc_threshold Log2 fold change threshold (default: 1)
#' @param padj_threshold Adjusted p-value threshold (default: 0.05)
#' @param direction Filter direction: "up", "down", or "both" (default: "both")
#' @return Filtered results data frame
filter_degs <- function(deg_results,
                       log2fc_threshold = 1,
                       padj_threshold = 0.05,
                       direction = "both") {
  
  # Convert to data frame if needed
  if (class(deg_results)[1] == "DESeqResults") {
    deg_results <- as.data.frame(deg_results)
  }
  
  # Remove NA values
  deg_results <- na.omit(deg_results)
  
  # Filter by adjusted p-value
  deg_filtered <- deg_results[deg_results$padj < padj_threshold, ]
  
  # Filter by fold change and direction
  if (direction == "up") {
    deg_filtered <- deg_filtered[deg_filtered$log2FoldChange > log2fc_threshold, ]
  } else if (direction == "down") {
    deg_filtered <- deg_filtered[deg_filtered$log2FoldChange < -log2fc_threshold, ]
  } else if (direction == "both") {
    deg_filtered <- deg_filtered[abs(deg_filtered$log2FoldChange) > log2fc_threshold, ]
  }
  
  return(deg_filtered)
}


# Seurat Helpers ----

#' Load and QC Seurat Object
#' 
#' Loads Seurat object and applies quality control filters
#' 
#' @param file_path Path to Seurat object (.rds file)
#' @param mt_threshold Maximum mitochondrial percentage (default: 25)
#' @param min_features Minimum features per cell (default: 200)
#' @return Filtered Seurat object
load_and_qc_seurat <- function(file_path,
                               mt_threshold = 25,
                               min_features = 200) {
  
  message("Loading Seurat object from: ", file_path)
  seurat_obj <- readRDS(file_path)
  
  # Calculate mitochondrial percentage if not present
  if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
    message("Calculating mitochondrial percentage...")
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  }
  
  # Apply QC filters
  message("Applying QC filters...")
  n_cells_before <- ncol(seurat_obj)
  
  seurat_obj <- subset(
    seurat_obj,
    subset = percent.mt <= mt_threshold & nFeature_RNA >= min_features
  )
  
  n_cells_after <- ncol(seurat_obj)
  message("Cells before QC: ", n_cells_before)
  message("Cells after QC: ", n_cells_after)
  message("Cells removed: ", n_cells_before - n_cells_after)
  
  return(seurat_obj)
}


#' Standard Seurat Processing Pipeline
#' 
#' Runs normalization, scaling, PCA, and UMAP on Seurat object
#' 
#' @param seurat_obj Seurat object
#' @param n_pcs Number of PCs to compute and use (default: 30)
#' @param resolution Clustering resolution (default: 0.1)
#' @return Processed Seurat object
process_seurat <- function(seurat_obj,
                           n_pcs = 30,
                           resolution = 0.1) {
  
  message("Normalizing data...")
  seurat_obj <- NormalizeData(seurat_obj)
  
  message("Finding variable features...")
  seurat_obj <- FindVariableFeatures(seurat_obj)
  
  message("Scaling data...")
  seurat_obj <- ScaleData(seurat_obj)
  
  message("Running PCA...")
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE, npcs = n_pcs)
  
  message("Running UMAP...")
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_pcs)
  
  message("Finding neighbors...")
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_pcs)
  
  message("Finding clusters...")
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  
  return(seurat_obj)
}


#' Perform Differential Expression in Seurat
#' 
#' Wrapper for FindMarkers with common parameters
#' 
#' @param seurat_obj Seurat object
#' @param ident.1 Identity class 1
#' @param ident.2 Identity class 2
#' @param group.by Metadata column to group by
#' @param logfc.threshold Log fold change threshold (default: 0.25)
#' @param min.pct Minimum percentage threshold (default: 0.25)
#' @param test.use Statistical test (default: "wilcox")
#' @return Data frame with differential expression results
perform_deg_analysis <- function(seurat_obj,
                                 ident.1,
                                 ident.2 = NULL,
                                 group.by = NULL,
                                 logfc.threshold = 0.25,
                                 min.pct = 0.25,
                                 test.use = "wilcox") {
  
  # Set identity if group.by specified
  if (!is.null(group.by)) {
    Idents(seurat_obj) <- seurat_obj[[group.by]]
  }
  
  message("Running differential expression: ", ident.1, " vs ", 
          ifelse(is.null(ident.2), "rest", ident.2))
  
  deg_results <- FindMarkers(
    seurat_obj,
    ident.1 = ident.1,
    ident.2 = ident.2,
    logfc.threshold = logfc.threshold,
    min.pct = min.pct,
    test.use = test.use,
    verbose = FALSE
  )
  
  return(deg_results)
}


#' Extract IL-22 Expressing Cells
#' 
#' Identifies and extracts cells expressing a specific gene
#' 
#' @param seurat_obj Seurat object
#' @param gene Gene name to check (default: "IL22")
#' @param assay Assay to use (default: "RNA")
#' @return Subset Seurat object with only expressing cells
extract_expressing_cells <- function(seurat_obj,
                                     gene = "IL22",
                                     assay = "RNA") {
  
  # Get expression data
  expr_data <- GetAssayData(seurat_obj, assay = assay)
  
  # Check if gene exists
  if (!gene %in% rownames(expr_data)) {
    stop("Gene ", gene, " not found in ", assay, " assay")
  }
  
  # Find expressing cells
  gene_expr <- expr_data[gene, ]
  expressing_cells <- names(gene_expr[gene_expr != 0])
  
  message("Total cells: ", ncol(seurat_obj))
  message(gene, " expressing cells: ", length(expressing_cells))
  message("Percentage: ", round(100 * length(expressing_cells) / ncol(seurat_obj), 2), "%")
  
  # Subset object
  seurat_subset <- seurat_obj[, expressing_cells]
  
  return(seurat_subset)
}


# Data Processing Helpers ----

#' Convert Ensembl IDs to Gene Symbols
#' 
#' Uses biomaRt to convert Ensembl gene IDs to gene symbols
#' 
#' @param ensembl_ids Character vector of Ensembl gene IDs
#' @param species Species for biomaRt (default: "hsapiens_gene_ensembl")
#' @return Data frame with ensembl_gene_id and external_gene_name columns
convert_ensembl_to_symbol <- function(ensembl_ids,
                                      species = "hsapiens_gene_ensembl") {
  
  suppressPackageStartupMessages(library(biomaRt))
  
  message("Connecting to Ensembl biomaRt...")
  ensembl <- useEnsembl(biomart = "ensembl", dataset = species)
  
  message("Retrieving gene symbols for ", length(ensembl_ids), " genes...")
  gene_info <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
  )
  
  # Remove empty gene names
  gene_info <- gene_info[gene_info$external_gene_name != "", ]
  
  message("Successfully mapped ", nrow(gene_info), " genes")
  
  return(gene_info)
}


#' Remove Duplicate Genes
#' 
#' Handles duplicate gene symbols by keeping first occurrence
#' 
#' @param data_matrix Data matrix with gene symbols as rownames
#' @return Matrix with duplicates removed
remove_duplicate_genes <- function(data_matrix) {
  
  n_genes_before <- nrow(data_matrix)
  duplicated_mask <- duplicated(rownames(data_matrix))
  data_matrix <- data_matrix[!duplicated_mask, ]
  n_genes_after <- nrow(data_matrix)
  
  if (n_genes_before > n_genes_after) {
    message("Removed ", n_genes_before - n_genes_after, " duplicate genes")
  }
  
  return(data_matrix)
}


#' Filter Zero Expression Genes
#' 
#' Removes genes with zero expression across all samples
#' 
#' @param count_matrix Count matrix (genes x samples)
#' @return Filtered count matrix
filter_zero_genes <- function(count_matrix) {
  
  n_genes_before <- nrow(count_matrix)
  gene_sums <- rowSums(count_matrix)
  count_matrix <- count_matrix[gene_sums != 0, ]
  n_genes_after <- nrow(count_matrix)
  
  message("Removed ", n_genes_before - n_genes_after, " genes with zero expression")
  
  return(count_matrix)
}


# Plotting Helpers ----

#' Save Plot with Standard Settings
#' 
#' Wrapper for saving plots with consistent parameters
#' 
#' @param filename Output filename
#' @param plot_function Function that generates the plot
#' @param width Plot width in inches (default: 7)
#' @param height Plot height in inches (default: 5)
#' @param device Device to use ("pdf", "png", default: "pdf")
#' @param dpi DPI for raster devices (default: 300)
#' @param cairo Use cairo for better rendering (default: TRUE)
save_plot_standard <- function(filename,
                               plot_function,
                               width = 7,
                               height = 5,
                               device = "pdf",
                               dpi = 300,
                               cairo = TRUE) {
  
  # Add device extension if not present
  if (!grepl(paste0("\\.", device, "$"), filename)) {
    filename <- paste0(filename, ".", device)
  }
  
  # Open device
  if (device == "pdf") {
    pdf(filename, width = width, height = height)
  } else if (device == "png") {
    if (cairo) {
      png(filename, width = width, height = height, units = "in", 
          res = dpi, type = "cairo")
    } else {
      png(filename, width = width, height = height, units = "in", res = dpi)
    }
  }
  
  # Execute plot function
  plot_function()
  
  # Close device
  dev.off()
  
  message("Plot saved to: ", filename)
}


# File I/O Helpers ----

#' Save Table with Standard Format
#' 
#' Wrapper for write.table with consistent parameters
#' 
#' @param data Data frame or matrix to save
#' @param filename Output filename
#' @param sep Column separator (default: tab)
#' @param quote Quote strings (default: FALSE)
#' @param row.names Include row names (default: TRUE)
save_table_standard <- function(data,
                               filename,
                               sep = "\t",
                               quote = FALSE,
                               row.names = TRUE) {
  
  write.table(data, filename, sep = sep, quote = quote, row.names = row.names)
  message("Table saved to: ", filename)
}


#' Save Gene List to File
#' 
#' Saves a character vector of gene names to file
#' 
#' @param genes Character vector of gene names
#' @param filename Output filename
save_gene_list <- function(genes, filename) {
  writeLines(genes, filename)
  message("Saved ", length(genes), " genes to: ", filename)
}


#' Read Gene List from File
#' 
#' Reads a gene list from file
#' 
#' @param filename Input filename
#' @return Character vector of gene names
read_gene_list <- function(filename) {
  genes <- readLines(filename)
  message("Loaded ", length(genes), " genes from: ", filename)
  return(genes)
}


# Validation Functions ----

#' Check Required Columns
#' 
#' Validates that a data frame contains required columns
#' 
#' @param data Data frame to check
#' @param required_cols Character vector of required column names
#' @param data_name Name of data frame for error message
check_required_columns <- function(data, required_cols, data_name = "data") {
  
  missing_cols <- setdiff(required_cols, colnames(data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", data_name, ": ", 
         paste(missing_cols, collapse = ", "))
  }
  
  return(TRUE)
}


#' Validate File Exists
#' 
#' Checks if file exists, stops with error if not
#' 
#' @param file_path Path to file
validate_file_exists <- function(file_path) {
  
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  return(TRUE)
}


# Summary Statistics ----

#' Calculate Gene Expression Summary Statistics
#' 
#' Computes mean, median, and percentage expressing for a gene
#' 
#' @param expression_vector Numeric vector of expression values
#' @return Named list with summary statistics
expression_summary <- function(expression_vector) {
  
  expressing_cells <- expression_vector[expression_vector != 0]
  
  return(list(
    n_total = length(expression_vector),
    n_expressing = length(expressing_cells),
    pct_expressing = 100 * length(expressing_cells) / length(expression_vector),
    mean_all = mean(expression_vector),
    mean_expressing = ifelse(length(expressing_cells) > 0, 
                             mean(expressing_cells), 
                             0),
    median_expressing = ifelse(length(expressing_cells) > 0,
                               median(expressing_cells),
                               0)
  ))
}


# Module Scoring ----

#' Add Module Score to Seurat Object
#' 
#' Wrapper for AddModuleScore with standard parameters
#' 
#' @param seurat_obj Seurat object
#' @param gene_list List of gene vectors for scoring
#' @param ctrl Number of control genes (default: 5)
#' @param name Module name prefix (default: "Module")
#' @return Seurat object with module scores added
add_module_score_wrapper <- function(seurat_obj,
                                     gene_list,
                                     ctrl = 5,
                                     name = "Module") {
  
  message("Adding module score for ", length(gene_list), " gene sets")
  
  seurat_obj <- AddModuleScore(
    seurat_obj,
    features = gene_list,
    ctrl = ctrl,
    name = name
  )
  
  return(seurat_obj)
}


# Print utility message
message("Common utility functions loaded successfully!")
message("Available functions:")
message("  - get_pseudobulk()")
message("  - fisher_test_enrichment()")
message("  - run_deseq2_paired()")
message("  - load_and_qc_seurat()")
message("  - process_seurat()")
message("  - perform_deg_analysis()")
message("  - extract_expressing_cells()")
message("  - convert_ensembl_to_symbol()")
message("  - save_plot_standard()")
message("  - save_table_standard()")
message("  See function documentation for usage details")