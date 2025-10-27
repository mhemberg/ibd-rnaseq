################################################################################
# IL-22 Differential Expression Analysis in T Cells
# IBD Transcriptomics Analysis
#
# Purpose: Analyze IL-22 expression and differential expression in T cell subsets
# Input: T cell scRNA-seq data (tcell2_badrm_harmony35)
# Output: DEG results, IL-22 expression statistics, and visualizations
#
# Original file: scRNAseq/il22/il22_deg.docx (dates: 20240916, 20240917, 20250113)
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
OUTPUT_PREFIX <- "il22_tcells"
IL22_GENE <- "IL22"

################################################################################
# Load Data
################################################################################

message("==========================================")
message("IL-22 T Cell Differential Expression Analysis")
message("==========================================\n")

# Load T cell data
message("Loading T cell data from: ", SCRNA_DATA_PATH)
validate_file_exists(SCRNA_DATA_PATH)
load(SCRNA_DATA_PATH)

# Store in consistent variable name
tcell_data <- d
rm(d)

message("Total cells loaded: ", ncol(tcell_data))
message("Total genes: ", nrow(tcell_data))

################################################################################
# Subset IL-22 Expressing Cell Types
################################################################################

message("\n--- Subsetting IL-22 expressing cell types ---")
message("Cell types: ", paste(IL22_CELLTYPES, collapse = ", "))

# Subset to IL-22 expressing cell types
tcell_subset <- subset(
  tcell_data,
  subset = anotation_intermediate %in% IL22_CELLTYPES
)

message("Cells after subsetting: ", ncol(tcell_subset))

# Calculate IL-22 expression statistics
message("\n--- IL-22 Expression Statistics ---")

# Overall statistics
expr_data <- GetAssayData(tcell_subset)
il22_expr <- expr_data[IL22_GENE, ]
n_il22_positive <- sum(il22_expr != 0)
pct_il22_positive <- 100 * n_il22_positive / ncol(tcell_subset)

message("IL-22+ cells in subset: ", n_il22_positive, " / ", ncol(tcell_subset))
message("Percentage: ", round(pct_il22_positive, 2), "%")

# Statistics by cell type
for (celltype in IL22_CELLTYPES) {
  celltype_cells <- subset(tcell_subset, subset = anotation_intermediate == celltype)
  celltype_expr <- GetAssayData(celltype_cells)
  celltype_il22 <- celltype_expr[IL22_GENE, ]
  n_positive <- sum(celltype_il22 != 0)
  pct_positive <- 100 * n_positive / ncol(celltype_cells)
  
  message(celltype, ": ", n_positive, " / ", ncol(celltype_cells), 
          " (", round(pct_positive, 2), "%)")
}

################################################################################
# Differential Expression Analysis
################################################################################

message("\n--- Running Differential Expression Analysis ---")

# Set active identity to patient group
Idents(tcell_subset) <- as.factor(tcell_subset$new_group)

# Perform DEG analysis for each comparison
for (i in seq_along(COMPARISON_SETS)) {
  
  comp_set <- COMPARISON_SETS[[i]]
  group1 <- comp_set[1]
  group2 <- comp_set[2]
  
  message("\nComparison: Group ", group1, " vs Group ", group2)
  tryCatch({  
    # Subset to only these two groups
    subset_comp <- subset(
    tcell_subset,
    subset = new_group %in% c(group1, group2)
    )

    if (ncol(subset_comp) < 10) {
       message("  Skipping - too few cells")
       next
    }
  
  # Run Wilcoxon test
  deg_results <- FindMarkers(
    subset_comp,
    ident.1 = group1,
    ident.2 = group2,
    logfc.threshold = 0.01,  # Lower threshold for exploratory analysis
    verbose = FALSE,
    test.use = "wilcox"
  )
  
  # Save results
  output_file <- file.path(
    DEG_DIR,
    paste0(OUTPUT_PREFIX, "_deg_", group1, "_vs_", group2, ".txt")
  )
  save_table_standard(deg_results, output_file)
  
  message("  DEGs found: ", nrow(deg_results))
  message("  Saved to: ", output_file)
  }, error = function(e) {
  message("Error comparing groups, skipping this comparison.\n")
})
}

################################################################################
# IL-22 Expression Analysis by Group
################################################################################

message("\n--- Analyzing IL-22 expression by patient group ---")

# Extract IL-22 expressing cells only
il22_positive_indices <- which(il22_expr != 0)
tcell_il22_pos <- tcell_subset[, il22_positive_indices]

message("IL-22+ cells for analysis: ", ncol(tcell_il22_pos))

# Compile IL-22 expression data by group
group_order <- c(5, 6, 4, 3, 7)  # CD uninfl, PFD uninfl, PFD infl, CD infl, HC
il22_expr_data <- data.frame()

for (group_id in group_order) {
  tryCatch({  
  
  group_cells <- subset(tcell_il22_pos, subset = new_group == group_id)
  
  if (ncol(group_cells) == 0) {
    message("Group ", group_id, ": No IL-22+ cells")
    next
  }
  
  group_expr <- GetAssayData(group_cells)
  group_il22 <- group_expr[IL22_GENE, ]
  group_il22_nonzero <- group_il22[group_il22 != 0]
  
  # Store data
  temp_df <- data.frame(
    IL22 = group_il22_nonzero,
    Group = group_id
  )
  il22_expr_data <- rbind(il22_expr_data, temp_df)
  
  # Calculate mean
  mean_expr <- mean(group_il22_nonzero)
  message("Group ", group_id, " - Mean IL-22 expression: ", round(mean_expr, 3))
  }, error = function(e) {
  message("Error evaluating group ", group_id, ", skipping.\n")
})
}

# Convert group to factor with correct order
il22_expr_data$Group <- factor(il22_expr_data$Group, levels = group_order)

# Statistical tests
message("\n--- Statistical comparisons ---")

# PFD uninfl vs CD uninfl (6 vs 5)
if (6 %in% il22_expr_data$Group && 5 %in% il22_expr_data$Group) {
  test_result <- wilcox.test(
    il22_expr_data$IL22[il22_expr_data$Group == 5],
    il22_expr_data$IL22[il22_expr_data$Group == 6],
    alternative = "two.sided"
  )
  message("CD uninfl vs PFD uninfl (5 vs 6): p-value = ", 
          format(test_result$p.value, scientific = TRUE, digits = 4))
}

# CD infl vs PFD infl (3 vs 4)
if (3 %in% il22_expr_data$Group && 4 %in% il22_expr_data$Group) {
  test_result <- wilcox.test(
    il22_expr_data$IL22[il22_expr_data$Group == 3],
    il22_expr_data$IL22[il22_expr_data$Group == 4],
    alternative = "two.sided"
  )
  message("CD infl vs PFD infl (3 vs 4): p-value = ", 
          format(test_result$p.value, scientific = TRUE, digits = 4))
}

################################################################################
# Visualization - Boxplot
################################################################################

message("\n--- Generating boxplot ---")

# Create group labels
group_labels <- c("5" = "CD uninfl", "6" = "PFD uninfl", 
                  "3" = "CD infl", "4" = "PFD infl", "7" = "HC")

# Boxplot of IL-22 expression
plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_expression_boxplot.pdf"))
pdf(plot_file, width = 4.8, height = 5)

boxplot(
  IL22 ~ Group,
  data = il22_expr_data,
  names = group_labels[levels(il22_expr_data$Group)],
  ylab = "IL-22 Expression",
  main = "IL-22 Expression in IL-22+ T Cells",
  las = 2,
  pch = 20
)

dev.off()
message("Boxplot saved to: ", plot_file)

################################################################################
# Visualization - Violin Plots
################################################################################

message("\n--- Generating violin plots ---")

# Extract IL-22 expressing cells with factor levels
tcell_il22_pos$group_factor <- factor(
  tcell_il22_pos$new_group,
  levels = c(5, 6, 4, 3, 7)
)
Idents(tcell_il22_pos) <- tcell_il22_pos$group_factor

tryCatch({# Violin plot - IL-22 expressing cells only
plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_expression_violin.pdf"))
pdf(plot_file, width = 7, height = 5)

print(VlnPlot(tcell_il22_pos, features = IL22_GENE, pt.size = 0))

dev.off()
message("Violin plot (IL-22+ cells) saved to: ", plot_file)
}, error = function(e) {
# Boxplot - IL-22 expressing cells only
plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_expression_violin.pdf"))
pdf(plot_file, width = 7, height = 5)

# Extract data
expr_data <- GetAssayData(tcell_il22_pos, slot = "data")
il22_values <- as.numeric(expr_data[IL22_GENE, ])
groups <- factor(tcell_il22_pos$new_group, levels = c(5, 6, 4, 3, 7))

# Simple boxplot
boxplot(il22_values ~ groups,
        col = c("lightblue", "lightcoral", "lightblue", "lightcoral", "lightgreen"),
        main = paste0(IL22_GENE, " Expression in IL-22+ T Cells"),
        xlab = "Patient Group",
        ylab = "Expression Level",
        names = c("CD uninfl", "PFD uninfl", "CD infl", "PFD infl", "Healthy"))

dev.off()
message("Boxplot (IL-22+ cells) saved to: ", plot_file)
})
tryCatch({# Violin plot - all T cells
tcell_data_il22 <- extract_expressing_cells(tcell_data, gene = IL22_GENE)
tcell_data_il22$group_factor <- factor(
  tcell_data_il22$new_group,
  levels = c(5, 6, 4, 3, 7)
)
Idents(tcell_data_il22) <- tcell_data_il22$group_factor

plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_all_tcells_expression_violin.pdf"))
pdf(plot_file, width = 7, height = 5)

print(VlnPlot(tcell_data_il22, features = IL22_GENE, pt.size = 0))

dev.off()
message("Violin plot (all IL-22+ T cells) saved to: ", plot_file)
}, error = function(e) {
# Boxplot - IL-22 expressing cells only
plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_all_tcells_expression_violin.pdf"))
pdf(plot_file, width = 7, height = 5)

# Extract data
expr_data <- GetAssayData(tcell_data_il22, slot = "data")
il22_values <- as.numeric(expr_data[IL22_GENE, ])
groups <- factor(tcell_data_il22$new_group, levels = c(5, 6, 4, 3, 7))

# Simple boxplot
boxplot(il22_values ~ groups,
        col = c("lightblue", "lightcoral", "lightblue", "lightcoral", "lightgreen"),
        main = paste0(IL22_GENE, " Expression in T Cells"),
        xlab = "Patient Group",
        ylab = "Expression Level",
        names = c("CD uninfl", "PFD uninfl", "CD infl", "PFD infl", "Healthy"))

dev.off()
message("Boxplot (All T cells) saved to: ", plot_file)
})

################################################################################
# Pseudobulk Analysis by Cell Type
################################################################################

message("\n--- Pseudobulk analysis by cell type ---")

# Get expression data
expr_matrix <- GetAssayData(tcell_data)
il22_row <- expr_matrix[IL22_GENE, , drop = FALSE]

# Calculate pseudobulk
pseudobulk_expr <- get_pseudobulk(
  il22_row,
  as.character(tcell_data$anotation_intermediate)
)

# Convert to numeric vector
pseudobulk_vec <- as.numeric(pseudobulk_expr)
names(pseudobulk_vec) <- colnames(pseudobulk_expr)

# Sort by expression
pseudobulk_sorted <- sort(pseudobulk_vec, decreasing = FALSE)

# Save results
output_file <- file.path(TABLES_DIR, paste0(OUTPUT_PREFIX, "_pseudobulk_by_celltype.txt"))
save_table_standard(
  data.frame(celltype = names(pseudobulk_sorted), 
             mean_IL22 = pseudobulk_sorted),
  output_file
)

# Barplot
plot_file <- file.path(FIGURES_DIR, paste0(OUTPUT_PREFIX, "_pseudobulk_barplot.pdf"))
pdf(plot_file, width = 8, height = 6)

barplot(
  pseudobulk_sorted,
  las = 2,
  ylab = "Mean IL-22 Expression",
  main = "IL-22 Expression by T Cell Subset",
  cex.names = 0.7
)

dev.off()
message("Pseudobulk barplot saved to: ", plot_file)

################################################################################
# IL-22+ vs Other T Cells DEG Analysis
################################################################################

message("\n--- IL-22+ vs other T cells DEG analysis ---")

# Create new grouping variable
tcell_data$il22_group <- "other"
tcell_data$il22_group[tcell_data$anotation_intermediate == "ILC3"] <- "il22_exp"
tcell_data$il22_group[tcell_data$anotation_intermediate == "Th17"] <- "il22_exp"
tcell_data$il22_group[tcell_data$anotation_intermediate == "CD4_effector_memory"] <- "il22_exp"
tcell_data$il22_group[tcell_data$anotation_intermediate == "CD4_effector_TNF"] <- "il22_exp"

# Set identity
Idents(tcell_data) <- as.factor(tcell_data$il22_group)

# Find markers
message("Finding markers for IL-22 expressing cell types...")
deg_il22_exp <- FindAllMarkers(tcell_data, only.pos = TRUE)

# Save results
output_file <- file.path(DEG_DIR, paste0(OUTPUT_PREFIX, "_il22exp_vs_other_markers.txt"))
save_table_standard(deg_il22_exp, output_file)

# Report IL-22 statistics
il22_stats <- deg_il22_exp[deg_il22_exp$gene == IL22_GENE, ]
if (nrow(il22_stats) > 0) {
  message("\nIL-22 gene statistics:")
  message("  avg_log2FC: ", round(il22_stats$avg_log2FC, 2))
  message("  p_val_adj: ", format(il22_stats$p_val_adj, scientific = TRUE, digits = 3))
}

################################################################################
# Summary Report
################################################################################

message("\n==========================================")
message("Analysis Summary")
message("==========================================")
message("Total T cells analyzed: ", ncol(tcell_data))
message("IL-22 expressing T cells: ", n_il22_positive, " (", round(pct_il22_positive, 2), "%)")
message("Number of comparisons performed: ", length(COMPARISON_SETS))
message("\nOutput files:")
message("  DEG results: ", DEG_DIR)
message("  Figures: ", FIGURES_DIR)
message("  Tables: ", TABLES_DIR)

# Save session info
save_session_info(file.path(RESULTS_DIR, paste0(OUTPUT_PREFIX, "_session_info.txt")))

message("\n==========================================")
message("IL-22 T cell analysis complete!")
message("==========================================")
