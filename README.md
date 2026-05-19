# ACS and NFAT Signaling: Integrative Multiomics Analysis

This repository contains R scripts for the analysis presented in the following paper:

> [Authors]. [Title]. *[Journal]*. [Year]. DOI: [DOI]

## Scripts

| Script | Description |
|--------|-------------|
| `src/00_setup.R` | Package loading and common settings |
| `src/01_baseline_table.R` | Baseline characteristics table |
| `src/02_vsr_analysis.R` | VSR analysis |
| `src/03_1_rnaseq_deseq2_umap.R` | RNA-seq: DESeq2 and UMAP |
| `src/03_2_rnaseq_decoupler.R` | RNA-seq: decoupleR pathway analysis |
| `src/03_3_rnaseq_enrichr.R` | RNA-seq: Enrichr gene set enrichment |
| `src/03_4_bisque.R` | Cell type deconvolution (BisqueRNA) |
| `src/04_metabolomics.R` | Metabolomics analysis |
| `src/05_wgcna_metabolomics.R` | WGCNA for metabolomics |
| `src/06_cell_rnaseq.R` | Single-cell RNA-seq analysis |
| `src/07_lipidomics.R` | Lipidomics analysis |
| `src/08_wgcna_lipidomics.R` | WGCNA for lipidomics |
| `src/09_mr_analysis.R` | Mendelian randomization analysis |
