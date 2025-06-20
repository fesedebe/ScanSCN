#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  source("scan_scn/singlecell/cell_state/signature_scoring.R")
})

run_score_signatures <- function(input_rds, output_rds) {
  obj <- readRDS(input_rds)
  gene_sets <- load_default_signature_sets()
  obj <- score_signatures(obj, gene_sets)
  saveRDS(obj, output_rds)
}

# Entry point
if (interactive() || identical(Sys.getenv("R_SCRIPT_DEBUG"), "TRUE")) {
  run_score_signatures(
    input_rds = "res/ucla_clustered.rds",
    output_rds = "res/ucla_clustered_with_signatures.rds"
  )
} else {
  run_score_signatures(
    input_rds = snakemake@input[[1]],
    output_rds = snakemake@output[[1]]
  )
}
