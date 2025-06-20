#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  source("scan_scn/singlecell/cell_state/normalize_and_integrate.R")
})

normalize_and_save <- function(input_path, output_path, method = "rpca") {
  sobj <- readRDS(input_path)

  integrated <- normalize_and_integrate(
    sobj = sobj,
    method = method,
    vars.to.regress = c("S.Score", "G2M.Score")
  )

  saveRDS(integrated, output_path)
  message("Saved integrated Seurat object to: ", output_path)
}

if (interactive() || identical(Sys.getenv("R_SCRIPT_DEBUG"), "TRUE")) {
  normalize_and_save(
    input_path = "res/ucla_qc_filtered.rds",
    output_path = "res/ucla_integrated.rds",
    method = "rpca"
  )
} else {
  input_path <- snakemake@input[[1]]
  output_path <- snakemake@output[[1]]
  normalize_and_save(input_path, output_path)
}