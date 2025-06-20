#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  source("scan_scn/singlecell/cell_state/cluster_and_score.R")
})

run_pipeline <- function(input_rds, output_rds) {
  sobj <- readRDS(input_rds)

  sobj <- cluster_and_score_cells(
    sobj = sobj,
    resolution = 0.4,
    reduction = "pca",
    save_path = output_rds
  )

  message("Clustering + scoring complete. Saved to ", output_rds)
}

# Main block
if (interactive() || identical(Sys.getenv("R_SCRIPT_DEBUG"), "TRUE")) {
  run_pipeline(
    input_rds = "res/ucla_integrated.rds",
    output_rds = "res/ucla_clustered_scored.rds"
  )
} else {
  run_pipeline(
    input_rds = snakemake@input[[1]],
    output_rds = snakemake@output[[1]]
  )
}