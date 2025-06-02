#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(scater)
  library(dplyr)
  source("scan_scn/singlecell/preprocessing/load_and_annotate.R")
  source("scan_scn/singlecell/preprocessing/decontX_and_doublet_detection.R")
  source("scan_scn/singlecell/preprocessing/qc_metrics.R")
  source("scan_scn/singlecell/preprocessing/filter_cells.R")
})

run_qc_pipeline <- function(matrix_dir, aggr_csv, output_path) {
  sobj <- load_raw_10X_matrix(matrix_dir, project = "ucla_ov")
  anno_df <- create_barcode_sample_annotation(sobj, aggr_csv)
  sobj <- add_sample_metadata(sobj, anno_df)

  sce <- as.SingleCellExperiment(sobj)
  sce <- run_decontX_wrapper(sce, sample_col = "sample")

  sce_clean <- apply_decontaminated_counts(sce, sce)
  sce_clean <- run_doublet_detection(sce_clean, sample_col = "sample")

  sobj_clean <- as.Seurat(sce_clean)
  sobj_clean <- add_qc_metrics(sobj_clean)

  sobj_clean$MT.3MADs <- flag_outliers_mad(sobj_clean$percent.mt, sobj_clean$sample)
  sobj_clean$feature.mad.hi <- flag_outliers_mad(sobj_clean$nFeature_RNA, sobj_clean$sample)
  sobj_clean$umi.mad.hi <- flag_outliers_mad(sobj_clean$nCount_RNA, sobj_clean$sample, nmads = 4)

  sobj_filt <- filter_cells(sobj_clean)
  saveRDS(sobj_filt, file = output_path)

  message("Final Seurat object saved to: ", output_path)
}

if (interactive() || identical(Sys.getenv("R_SCRIPT_DEBUG"), "TRUE")) {
  run_qc_pipeline(
    matrix_dir = "data/input/filtered_feature_bc_matrix",
    aggr_csv = "data/internal/aggregation_1.csv",
    output_path = "res/ucla_qc_filtered.rds"
  )
} else {
  matrix_dir <- snakemake@input[["matrix_dir"]]
  aggr_csv <- snakemake@input[["aggr_csv"]]
  output_path <- snakemake@output[["output"]]

  run_qc_pipeline(matrix_dir, aggr_csv, output_path)
}