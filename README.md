# IBD Transcriptomics Analysis

Comprehensive analysis of transcriptomic data (scRNA-seq, bulk RNA-seq, and microarray) related to inflammatory bowel disease (IBD), focusing on cytokine signaling (IL-22, IL-13, TGF-β, TNF-α, TL1A) in intestinal cells.

## Overview

This repository contains analysis pipelines for:
- Single-cell RNA sequencing (scRNA-seq) data
- Bulk RNA sequencing data
- Microarray data (IL-22 and IL-13 studies)
- Cell-cell interaction analysis
- Differential gene expression analysis
- Functional enrichment and pathway analysis

## Project Structure

```
ibd-transcriptomics/
├── README.md                           # This file
├── config.R                            # Central configuration file
├── data/                               # Data directory (see Data section)
│   ├── raw/                           # Raw data files
│   ├── processed/                     # Processed data objects
│   └── external/                      # External reference datasets
├── scripts/
│   ├── 01_preprocessing/              # Data preprocessing scripts
│   │   ├── process_validation_data.R
│   │   ├── process_il13_microarray.R
│   │   └── process_bulk_rnaseq.R
│   ├── 02_quality_control/            # QC and filtering
│   │   └── qc_scrna.R
│   ├── 03_differential_expression/     # DEG analysis
│   │   ├── deg_il22_tcells.R
│   │   ├── deg_bulk_rnaseq.R
│   │   └── deg_validation.R
│   ├── 04_functional_analysis/        # Pathway and enrichment analysis
│   │   ├── fisher_enrichment.R
│   │   ├── module_scoring.R
│   │   └── cytokine_signatures.R
│   ├── 05_cell_interactions/          # Cell-cell communication
│   │   ├── cellphonedb_analysis.R
│   │   └── ligand_receptor_heatmap.R
│   ├── 06_visualization/              # Plotting scripts
│   │   └── generate_figures.R
│   └── utils/                         # Utility functions
│       ├── seurat_helpers.R
│       ├── deseq_helpers.R
│       ├── pseudobulk.R
│       └── statistical_tests.R
├── results/
│   ├── figures/                       # Output plots
│   ├── tables/                        # Output tables
│   └── deg_results/                   # DEG results
├── references/                        # Reference files
│   ├── genome/                        # Genome files
│   └── annotations/                   # Annotation files
└── docs/                             # Additional documentation
    ├── analysis_workflow.md
    └── data_sources.md
```

## Installation

### Prerequisites

- R version ≥ 4.2.0
- RStudio (recommended)
- Sufficient disk space for genomic data (~50GB)

### Required R Packages

```r
# CRAN packages
install.packages(c("ggplot2", "dplyr", "tidyr", "pheatmap"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "Seurat",
  "DESeq2",
  "edgeR",
  "biomaRt",
  "Rsubread"
))
```

### External Tools (for bulk RNA-seq)

- **STAR aligner** v2.7.10a or later
- **featureCounts** (part of Rsubread package)

## Quick Start

### 1. Clone the Repository

```bash
git clone https://github.com/mhemberg/ibd-transcriptomics.git
cd ibd-transcriptomics
```

### 2. Configure Paths

Edit `config.R` to set your local paths:

```r
# Open config.R and modify:
ROOT_DIR <- "/path/to/your/project"  # Change this line
```

### 3. Download Data

See `docs/data_sources.md` for data download instructions. Key datasets:

- **GSE233063**: Validation cytokine treatment scRNA-seq data
- **GSE190705**: IL-13 microarray data
- **Custom datasets**: Contact authors for access to proprietary datasets

### 4. Run Analysis Pipeline

```r
# Start R session
source("config.R")

# Check package installation
check_packages()

# Run preprocessing
source("scripts/01_preprocessing/process_validation_data.R")

# Continue with subsequent analysis steps...
```

## Analysis Workflows

### Workflow 1: scRNA-seq IL-22 Analysis

Analyzes IL-22 expressing T cells and ILC3s in IBD patients.

```r
source("config.R")
source("scripts/utils/seurat_helpers.R")
source("scripts/03_differential_expression/deg_il22_tcells.R")
```

**Key steps:**
1. Load scRNA-seq data
2. Subset IL-22+ cell types
3. Differential expression between patient groups
4. Pseudobulk analysis
5. Generate visualizations

### Workflow 2: Bulk RNA-seq Cytokine Response

Analyzes bulk RNA-seq from cytokine-treated intestinal organoids.

```r
source("config.R")
source("scripts/01_preprocessing/process_bulk_rnaseq.R")
source("scripts/03_differential_expression/deg_bulk_rnaseq.R")
```

**Key steps:**
1. STAR alignment (command line)
2. Feature counting
3. DESeq2 differential expression
4. Signature module scoring

### Workflow 3: Validation Dataset Analysis

Analyzes GSE233063 cytokine treatment validation data.

```r
source("config.R")
source("scripts/01_preprocessing/process_validation_data.R")
source("scripts/04_functional_analysis/fisher_enrichment.R")
```

**Key steps:**
1. Load and merge treatment groups
2. Quality control filtering
3. Clustering and cell type annotation
4. Fisher exact test enrichment

### Workflow 4: IL-13 Microarray Analysis

Processes and analyzes IL-13 treatment microarray data.

```r
source("config.R")
source("scripts/01_preprocessing/process_il13_microarray.R")
source("scripts/04_functional_analysis/il13_module_scoring.R")
source("scripts/06_visualization/il13_signature_heatmap.R")
```

**Key steps:**
1. Download data from GEO
2. Gene ID conversion (Ensembl to Gene Symbol)
3. DESeq2 analysis with paired design
4. Module scoring in epithelial cells
5. Generate signature heatmaps

### Workflow 5: TL1A and IL-22 Receptor Analysis

Analyzes TNFSF15 (TL1A) and IL22RA1 expression in epithelial cells.

```r
source("config.R")
source("scripts/03_differential_expression/deg_tl1a_il22ra1_epithelial.R")
```

**Key steps:**
1. Calculate expression statistics by cell type
2. Identify TL1A-expressing cell populations
3. Differential expression between patient groups
4. Non-zero cell analysis
5. Find markers for TL1A-high cells

### Complete Pipeline

To run the entire analysis pipeline:

```r
source("run_complete_analysis.R")
```

This will execute all preprocessing, DEG, functional analysis, and visualization steps in sequence.

## Key Functions

### Utility Functions (`scripts/utils/`)

#### `get_pseudobulk()`
Calculates pseudobulk expression by averaging across cell types.

```r
pseudobulk_expr <- get_pseudobulk(
  expression_matrix = seurat_obj[["RNA"]]@data,
  celltype_labels = seurat_obj$celltype
)
```

#### `fisher_test_enrichment()`
Performs Fisher exact test for gene set enrichment.

```r
results <- fisher_test_enrichment(
  query_genes = deg_genes,
  reference_genes = pathway_genes,
  universe_size = 33538
)
```

#### `run_deseq2_paired()`
Runs DESeq2 with paired sample design.

```r
results <- run_deseq2_paired(
  count_matrix = counts,
  metadata = sample_info,
  design_formula = ~ sample_id + treatment
)
```

#### `deg_special()`
Performs comprehensive DEG analysis including non-zero cell statistics.

```r
results <- deg_special(
  seurat_obj = epi_data,
  cell_types = c("PLCG2_enterocytes", "TA"),
  gene = "TNFSF15",
  group1_id = 6,
  group2_id = 5
)
```

## Data Description

### scRNA-seq Data

- **tcell2_badrm_harmony35**: T cell subset, Harmony batch-corrected
- **stroma2_badrm_harmony30**: Stromal cell subset
- **epi2_badrm_harmony50**: Epithelial cell subset
- **Patient groups**:
  - CD uninflamed (group 5)
  - PFD uninflamed (group 6)
  - CD inflamed (group 3)
  - PFD inflamed (group 4)
  - Healthy control (group 7)

### Key Cell Types

**T Cells:**
- ILC3: Innate lymphoid cells type 3
- Th17: T helper 17 cells
- CD4_effector_memory: CD4+ effector memory T cells
- CD4_effector_TNF: CD4+ TNF-producing effector cells

**Epithelial Cells:**
- PLCG2_enterocytes: PLCG2+ enterocytes
- TA: Transit-amplifying cells
- M_cells: Microfold cells
- Entero_AQP8pos: AQP8+ enterocytes
- RBHI_goblet: RBHI+ goblet cells

**Stromal Cells:**
- S1-S4: Stromal fibroblast subtypes

### Bulk RNA-seq

- Intestinal organoids treated with:
  - Lymphotoxin beta (LTB)
  - TGF-β
  - TNF-α
- Paired design (treated vs control from same donor)

### Microarray

- **IL-22**: Epithelial responses (custom data)
- **IL-13**: Epithelial responses (GSE190705)

## Comparison Sets

Standard comparisons across analyses:
1. PFD uninflamed vs CD uninflamed (6 vs 5)
2. CD inflamed vs PFD inflamed (3 vs 4)
3. PFD inflamed vs CD uninflamed (4 vs 5)
4. PFD uninflamed vs CD inflamed (6 vs 3)
5. PFD uninflamed vs Healthy (6 vs 7)
6. Healthy vs CD uninflamed (7 vs 5)

## Output Files

### Figures
- UMAP plots: `results/figures/*_umap.pdf`
- Violin plots: `results/figures/*_violin.pdf`
- Heatmaps: `results/figures/*_heatmap.pdf`
- Barplots: `results/figures/*_barplot.pdf`
- Boxplots: `results/figures/*_boxplot.pdf`
- MA plots and volcano plots: `results/figures/*_ma_plot.pdf`, `*_volcano_plot.pdf`

### Tables
- DEG results: `results/deg_results/*_deg.txt`
- Fisher test results: `results/tables/*_fisher.txt`
- Module scores: `results/tables/*_module_scores.txt`
- Gene lists: `results/tables/*_genes.txt`
- Expression statistics: `results/tables/*_expression.txt`

### Analysis-Specific Outputs

**IL-22 T Cell Analysis:**
- `il22_tcells_deg_6_vs_5.txt`: DEG results
- `il22_tcells_expression_violin.pdf`: Expression violin plots
- `il22_tcells_pseudobulk_by_celltype.txt`: Pseudobulk expression

**IL-13 Analysis:**
- `il13_microarray_deseq2_results.txt`: DESeq2 results
- `il13_module_scoring_scores.txt`: Module scores per cell
- `il13_signature_heatmap_scores_heatmap.pdf`: Signature heatmaps

**TL1A/IL22RA1 Analysis:**
- `tl1a_il22ra1_epithelial_tl1a_deg_results.txt`: TL1A DEG results
- `tl1a_il22ra1_epithelial_tl1a_mean_barplot.pdf`: Expression barplots
- `tl1a_il22ra1_epithelial_tl1a_high_cell_markers.txt`: TL1A-high cell markers

## Analysis Scripts Summary

### Preprocessing Scripts (scripts/01_preprocessing/)
1. **process_validation_data.R**: GSE233063 cytokine treatment scRNA-seq data
   - Loads 10X data for 7 treatment conditions
   - QC filtering, normalization, PCA, UMAP
   - Clustering and marker identification
2. **process_il13_microarray.R**: IL-13 microarray (GSE190705)
   - Auto-downloads from GEO
   - Ensembl ID to gene symbol conversion
   - DESeq2 with paired design
3. **process_bulk_rnaseq.R**: Bulk RNA-seq (LTB, TGF-β, TNF-α)
   - STAR alignment workflow
   - featureCounts quantification
   - DESeq2 and module scoring

### Differential Expression Scripts (scripts/03_differential_expression/)
1. **deg_il22_tcells.R**: IL-22 expression in T cell subsets
   - Analyzes ILC3, Th17, CD4+ effector cells
   - Pseudobulk analysis
   - Wilcoxon tests across patient groups
2. **deg_tl1a_il22ra1_epithelial.R**: TL1A (TNFSF15) and IL-22 receptor
   - Expression patterns in epithelial subtypes
   - Non-zero cell analysis
   - Marker identification
3. **deg_bulk_rnaseq.R**: Cytokine-treated organoid DEG
   - DESeq2 with donor-paired design
   - Multiple fold-change thresholds

### Functional Analysis Scripts (scripts/04_functional_analysis/)
1. **fisher_enrichment.R**: Gene set enrichment via Fisher exact test
   - Tests overlap with stromal cell signatures
   - Multiple comparison correction
   - Batch processing across treatment groups
2. **il13_module_scoring.R**: IL-13 response signature scoring
   - AddModuleScore in epithelial cells
   - Statistical comparisons across groups
   - Multiple q-value thresholds
3. **cell_abundance_fisher.R**: Cell type abundance testing
   - Tests differential abundance between groups
   - Multiple annotation levels (main, group, intermediate, refined)
   - Handles rare cell types
4. **module_scoring.R**: General pathway analysis
   - KEGG, Reactome, GO pathways
   - Custom gene set scoring

### Visualization Scripts (scripts/06_visualization/)
1. **il13_signature_heatmap.R**: IL-13 signature heatmaps
   - Cell type ordering by main category
   - Handles infinite values
   - Difference calculations
2. **generate_figures.R**: Publication-quality figure generation
   - Standardized color schemes
   - Consistent formatting

### Utility Scripts (scripts/utils/)
1. **common_functions.R**: 30+ utility functions
   - `get_pseudobulk()`: Pseudobulk expression
   - `fisher_test_enrichment()`: Fisher exact testing
   - `run_deseq2_paired()`: DESeq2 wrapper
   - `deg_special()`: Comprehensive DEG with non-zero analysis
   - `load_and_qc_seurat()`: Load and QC Seurat objects
   - `process_seurat()`: Standard Seurat pipeline
   - `convert_ensembl_to_symbol()`: ID conversion
   - `save_plot_standard()`: Consistent plotting
   - And many more...

### Python Utilities (scripts/utils/)
1. **pat_var_gn_consensus_rand.py**: Consensus variable gene selection
   - Identifies genes variable across patients
   - Random tie-breaking for consensus
2. **deg_threshold.py**: DEG result filtering
   - Applies q-value and fold-change thresholds
   - Cell ratio filtering

## Citation

If you use this code, please cite:

Gudino et al, [TL1A-activated T cells remodel the rectal mucosa in Crohn’s disease patients with perianal fistulizing disease](https://www.biorxiv.org/content/10.1101/2025.06.26.657455.full.pdf)

## Notes

- Ensure sufficient memory (≥32GB RAM recommended for full scRNA-seq analysis)
- Processing times vary: preprocessing (~1-2 hours), DEG analysis (~30 min - 2 hours)
- All hardcoded paths have been parameterized in `config.R`
- Korean comments in original code have been translated to English
- Analysis reproducibility ensured through session info logging


## Data Description

### scRNA-seq Data

- **tcell2_badrm_harmony35**: T cell subset, Harmony batch-corrected
- **stroma2_badrm_harmony30**: Stromal cell subset
- **Patient groups**:
  - CD uninflamed (group 5)
  - PFD uninflamed (group 6)
  - CD inflamed (group 3)
  - PFD inflamed (group 4)
  - Healthy control (group 7)

### Bulk RNA-seq

- Intestinal organoids treated with:
  - Lymphotoxin beta (LTB)
  - TGF-β
  - TNF-α
- Paired design (treated vs control from same donor)

### Microarray

- **IL-22**: Epithelial responses
- **IL-13**: Epithelial responses (GSE190705)

## Comparison Sets

Standard comparisons across analyses:
1. PFD uninflamed vs CD uninflamed (6 vs 5)
2. CD inflamed vs PFD inflamed (3 vs 4)
3. PFD inflamed vs CD uninflamed (4 vs 5)
4. PFD uninflamed vs CD inflamed (6 vs 3)
5. PFD uninflamed vs Healthy (6 vs 7)
6. Healthy vs CD uninflamed (7 vs 5)

## Output Files

### Figures
- UMAP plots: `results/figures/*_umap.pdf`
- Violin plots: `results/figures/*_violin.pdf`
- Heatmaps: `results/figures/*_heatmap.pdf`

### Tables
- DEG results: `results/deg_results/*_deg.txt`
- Fisher test results: `results/tables/*_fisher.txt`
- Module scores: `results/tables/*_module_scores.txt`

