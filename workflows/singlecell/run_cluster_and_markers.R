#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  source("scan_scn/singlecell/cell_state/cluster_and_identify_markers.R")
})

# Main block
run_cluster_and_markers <- function(input_rds, out_rds, out_markers) {
  sobj <- readRDS(input_rds)

  clustered <- run_leiden_clustering(sobj)
  saveRDS(clustered, out_rds)

  markers <- find_cluster_markers(clustered)
  write.table(markers, out_markers, sep = "\t", quote = FALSE, row.names = FALSE)
}

if (interactive() || identical(Sys.getenv("R_SCRIPT_DEBUG"), "TRUE")) {
  run_cluster_and_markers(
    input_rds = "data/ucla_integrated.rds",
    out_rds = "res/ucla_clustered.rds",
    out_markers = "res/ucla_cluster_markers.tsv"
  )
} else {
  run_cluster_and_markers(
    snakemake@input[[1]],
    snakemake@output[[1]],
    snakemake@output[[2]]
  )
}
