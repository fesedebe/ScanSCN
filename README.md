# ScanSCN

**ScanSCN** is a toolkit for identifing and analyzing small cell neuroendocrine (SCN) molecular features in transcriptomics datasets. 

## Overview
It supports R, Python and Snakemake workflows for both bulk and single-cell RNA-seq data and is applicable across a range of cancer types.

## Features

### Bulk RNA-seq
- **Signature Enrichment**: Enrichment against curated resistance gene sets
- **SCN Scoring**: PCA projection of query expression onto SCN reference
- **Similarity Assessment**: Compare DEGs, expression, and pathway trends vs known SCN profiles

### Single-cell RNA-seq
- **Preprocessing**: Mitochondrial filtering, doublet detection (scDblFinder), ambient RNA correction (DecontX)
- **Normalization & Integration**: SCTransform with optional cell cycle regression, RPCA or Harmony integration
- **Clustering & DE**: Leiden clustering, MAST-based marker identification
- **Trajectory Analysis**: Monocle 3 pseudotime inference
- **Signature Scoring**: AddModuleScore for SCN and other resistance signatures
- **SCN Gene Prioritization**: Cell classification → pseudobulk aggregation → DE feature selection → Ensemble ML model training for feature importance ranking