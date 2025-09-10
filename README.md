# ScanSCN
A predictive modeling pipeline to quantify rare tumor states and identify their actionable gene drivers

## Overview
Tumors can shift into a small-cell neuroendocrine (SCN) state to escape therapy. Detecting these states in transcriptomic data can be challenging: they are rare, heterogeneous, and defined by subtle transcriptional programs. ScanSCN provides a predictive modeling approach to quantify SCN states and identify the gene drivers underlying them. It takes in both bulk and single-cell RNA-seq data, producing:

- **SCN scores** → continuous values indicating how “SCN-like” a sample/subpopulation is.
- **Gene driver rankings** → interpretable feature importance lists that highlight genes contributing most to SCN variability.

## Features
- Linear projection model: Applies PCA-derived gene weights from reference datasets to calculate SCN scores in new cohorts.
- Single-cell workflow: Identifies candidate SCN subpopulations, aggregates pseudobulks, and applies predictive models (RF, XGBoost) to quantify SCN-ness.
- Feature importance analysis: Extracts interpretable gene rankings (Gini/gain importance) as candidate biomarkers.

## Installation
Clone and install:

```bash
git clone https://github.com/fesedebe/ScanSCN.git
cd scan_scn
pip install -e .
```
## Repo Structure
```
scan_scn/bulk                           # Signature Enrichment, SCN Scoring & Similarity Assessment
scan_scn/singlecell/preprocessing       # Mitochondrial Filtering, Doublet Detection, Ambient RNA Correction
scan_scn/singlecell/cell_state          # Normalization & Integration, Clustering & Marker ID, Trajectory Analysis
scan_scn/singlecell/scn_gene_importance # Pseudobulk, Ensemble Modeling & Feature Importance for Gene Drivers
workflows                               # Snakemake Workflows
```