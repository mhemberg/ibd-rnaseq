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
│   │   ├── process_validation_data.R
│   │   ├── process_il13_microarray.R
│   │   └── process_bulk_rnaseq.R
│   │   └── qc_scrna.R
│   │   ├── deg_il22_tcells.R
│   │   ├── deg_bulk_rnaseq.R
│   │   └── deg_validation.R
│   │   ├── fisher_enrichment.R
│   │   ├── module_scoring.R
│   │   └── cytokine_signatures.R
│   │   ├── cellphonedb_analysis.R
│   │   └── ligand_receptor_heatmap.R
│   │   └── config.R
│       └── utils_functions.R
├── results/
│   ├── figures/                       # Output plots
│   ├── tables/                        # Output tables
│   └── deg_results/                   # DEG results
├── references/                        # Reference files
│   ├── genome/                        # Genome files
│   └── annotations/                   # Annotation files
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

# Continue with subsequent analysis steps...
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


## Citation

If you use this code, please cite:

Gudino et al, [TL1A-activated T cells remodel the rectal mucosa in Crohn’s disease patients with perianal fistulizing disease](https://www.biorxiv.org/content/10.1101/2025.06.26.657455.full.pdf)



