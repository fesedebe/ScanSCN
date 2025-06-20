#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(ggplot2)
  source("scan_scn/singlecell/cell_state/trajectory_inference.R")
})

infer_pseudotime_main <- function(input_rds, output_rds, output_plot) {
  sobj <- readRDS(input_rds)

  result <- infer_pseudotime(sobj, method = "monocle", save_path = output_rds)

  ggsave(output_plot, result$plots$pseudotime, width = 6, height = 5)
  message("✅ Pseudotime plot and RDS saved.")
}

if (interactive() || identical(Sys.getenv("R_SCRIPT_DEBUG"), "TRUE")) {
  infer_pseudotime_main(
    input_rds = "data/ucla_integrated.rds",
    output_rds = "res/ucla_monocle_pseudotime.rds",
    output_plot = "res/ucla_monocle_pseudotime.pdf"
  )
} else {
  infer_pseudotime_main(
    snakemake@input[[1]],
    snakemake@output[[1]],
    snakemake@output[[2]]
  )
}