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
