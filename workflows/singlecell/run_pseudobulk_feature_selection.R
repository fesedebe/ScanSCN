#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  source("scan_scn/singlecell/scn_gene_importance/pseudobulk_feature_selection.R") # nolint
})

# ---- Wrapper Function ----
run_pseudobulk_feature_selection <- function(seurat_path, output_path = NULL, threshold = 0.03) {
  sobj <- readRDS(seurat_path)

  results <- pseudobulk_feature_selection::run_deseq2_on_pseudobulk(
    sobj = sobj,
    scn_score_col = "SCNBalanis_AUC",
    value_thresh = threshold
  )

  if (is.null(output_path)) {
    output_path <- sub("\\.rds$", "_scn_DE_results.txt", seurat_path)
    output_path <- gsub("data/", "res/", output_path)
  }

  write.table(results, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Results written to: ", output_path)
}

# ---- Execution block ----
if (interactive() || identical(Sys.getenv("R_SCRIPT_DEBUG"), "TRUE")) {
  run_pseudobulk_feature_selection(
    seurat_path = "data/ucla.sobj.dx.filt.epi.pat.sct_ccr_integ.srpca_leidenclust_sigscores.rds",
    output_path = "res/scn_DE_results_thresh003.txt",
    threshold = 0.03
  )
} else {
  seurat_path <- snakemake@input[[1]]
  output_path <- snakemake@output[[1]]
  threshold <- as.numeric(snakemake@params[["threshold"]])

  run_pseudobulk_feature_selection(
    seurat_path = seurat_path,
    output_path = output_path,
    threshold = threshold
  )
}