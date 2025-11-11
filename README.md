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
ibd-rnaseq/
├── README.md                           # This file
├── data/                               # Data directory (see Data section)
│   ├── raw/                           # Raw data files
│   │	├── GSE233063/                 
│   │	├── bulk_rna/                 
│   │	    └── star/                 
│   │	├── cellphonedb/                 
│   │	├── cytosig/                 
│   │	└── cellranger_data/                 
│   │	    └── geo_processed/                 
│   └── processed/                     # Processed data objects
├── scripts/
│   ├── abundance_fisher_tests.R
│   ├── process_validation.R
│   ├── process_il13_microarray.R
│   ├── process_bulk_rna.R
│   ├── process_scrna_main.R
│   ├── deg_il22_script.R
│   ├── deg_validation.R
│   ├── fisher_enrichment.R
│   ├── il13_heatmap.R
│   ├── il13_preprocess.R
│   ├── il13_module_score.R
│   ├── pathway_module_scoring.R
│   ├── cytosig_analysis.R
│   ├── cellphonedb_analysis.R
│   ├── run_all_script.R
│   ├── tl1a_il22ra1_analysis.R
│   ├── config.R
│   └── utils_functions.R
├── results/
│   ├── figures/                       # Output plots
│   ├── tables/                        # Output tables
│   └── deg_results/                   # DEG results
```

## Installation

### Prerequisites

- R version ≥ 4.2.0
- Sufficient disk space for genomic data

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
git clone https://github.com/mhemberg/ibd-rnaseq.git
cd ibd-rnaseq
```

### 2. Configure Paths

Edit `config.R` to set your local paths:

```r
# Open config.R and modify:
ROOT_DIR <- "/path/to/your/project"  # Change this line
```

### 3. Run Analysis Pipeline

```r
# Start R session
source("config.R")

# Check package installation
check_packages()

# Continue with subsequent analysis steps...
```

## Data Description

### scRNA-seq Data generated as part of this project, can be downloaded from GEO

Run process_scrna_main.R to process the outputs from Cell Ranger which should be stored in data/raw/cellranger_data/geo_processed (one subdirectory for each patient). You will also need the sample_info.txt file (part of this repository). This will result in datasets for the three main compartments of interest:

- **tcell2_badrm_harmony35**: T cell subset, Harmony batch-corrected
- **stroma2_badrm_harmony30**: Stromal cell subset
- **epi2_badrm_harmony50**: Epithelial cell subset
- **Patient groups**:
  - CD uninflamed (group 5)
  - PFD uninflamed (group 6)
  - CD inflamed (group 3)
  - PFD inflamed (group 4)
  - Healthy control (group 7)

### Cellphonedb data

To run cellphonedb_analysis.R you need processed data which can be downloaded from ....

### Cytosig data

For the cytosig_analysis.R script you need the ibd_total.Zscore file from ...

### Pathway analysis

To run pathway_module_scoring.R you need to download GO, Reactome, and KEGG files from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb)

### scRNA-seq Data from external sources

The process_validation.R script requires data from [GSE233063](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233063).

### Bulk RNA-seq

To run process_bulk_rna.R you need the bulk RNA seq data aligned with STAR from ...

### Microarray **GSE190705**

For process_il13_microarray.R we need data from both samples treated with IL-22 [GSE190705](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190705) and IL-13 (custom data).

## Comparison Sets

For most analyses, we consider the following groups:

1. PFD uninflamed vs CD uninflamed (6 vs 5)
2. CD inflamed vs PFD inflamed (3 vs 4)
3. PFD inflamed vs CD uninflamed (4 vs 5)
4. PFD uninflamed vs CD inflamed (6 vs 3)
5. PFD uninflamed vs Healthy (6 vs 7)
6. Healthy vs CD uninflamed (7 vs 5)


## Citation

If you use this code, please cite:

Gudino et al, [TL1A-activated T cells remodel the rectal mucosa in Crohn’s disease patients with perianal fistulizing disease](https://www.biorxiv.org/content/10.1101/2025.06.26.657455.full.pdf)



