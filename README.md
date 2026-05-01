# STE-DC2I

## Overview
This repository contains the R implementation of **STE-DC2I** for identifying subtype-specific driver genes and dark-causal dependencies from pseudotime-resolved single-cell RNA-seq data.

This repository now provides:
1. a **fully reproducible example workflow** for **GSE146771 subtype C14**;
2. all core R scripts used in the manuscript;
3. an **R environment lock file** (`renv.lock`);
4. a **Dockerfile** for containerized execution;
5. an example dataset and expected output files.

---

## Repository structure
```text
STE-DC2I/
├── README.md
├── renv.lock
├── Dockerfile
├── reproduce_GSE146771_C14.Rmd
├── run_reproduce.R
├── 1_data_preprocessing.R
├── 2_target_marker_gene.R
├── 3_dark_causal_infer_analysis.R
├── 4_comparision_validation_analysis.R
├── 5_enrichment_analysis.R
├── dark_causal_infer_fun.R
├── data/
│   └── [original or processed input files]
├── example_data/
│   ├── GSE146771_C14_demo.rds
│   └── metadata.tsv
├── output_example/
│   ├── driver_rankings.tsv
│   ├── dark_causal_pairs.tsv
│   ├── enrichment_results.tsv
│   └── figures/
├── Fig.3/
├── Fig.4/
├── Fig.5/
└── sessionInfo.txt
```

## Workflow description
Step 1. Data preprocessing
Script: 1_data_preprocessing.R
Functions:
quality control
normalization
variable gene filtering
malignant cell subset preparation

Step 2. Target marker identification
Script: 2_target_marker_gene.R
Functions:
subtype-specific marker extraction
candidate target marker definition

Step 3. Dark causality inference

Script: 3_dark_causal_infer_analysis.R
Core function file:
dark_causal_infer_fun.R
Functions:
pseudotime trajectory processing
state-space reconstruction
symbolic trend encoding
history projection / nearest-neighbor matching
dark causality scoring
driver gene prioritization

Step 4. Comparison and validation
Script: 4_comparision_validation_analysis.R
Functions:
benchmarking
overlap analysis
validation metric calculation

Step 5. Enrichment analysis
Script: 5_enrichment_analysis.R
Functions:
GO / KEGG / Reactome / DO enrichment
driver-specific functional annotation

## Software requirements
R: 4.4.3
RStudio: optional
Operating systems tested: Linux
Package management: renv

Required R packages
Main dependencies include, but are not limited to:
Seurat
monocle3
tidyverse
data.table
Matrix
ggplot2
clusterProfiler / enrichplot / DOSE
org.Hs.eg.db
patchwork
reshape2
pheatmap
igraph
patchwork
ggrepel
tidyr
dplyr

## Quick start (recommended: renv)
1. Clone the repository 
git clone https://github.com/linxi159/STE-DC2I.git
cd STE-DC2I

2. Install renv
Open R and run:
install.packages("renv")
renv::restore()

3. Run the reproducible example
rmarkdown::render("reproduce_GSE146771_C14.Rmd")
OR:
source("run_reproduce.R")

## Docker-based reproducibility
A Docker environment is provided for containerized execution.
Build image

docker build -t ste-dc2i-r .
Run example
docker run --rm -it -v $(pwd):/work ste-dc2i-r Rscript run_reproduce.R






