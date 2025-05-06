
source("scan_scn/bulk/enrichment_fgsea.R")
source("scan_scn/bulk/enrichment_plotting.R")

run_enrichment_and_plot <- function(deg_path, gmt_path, output_plot, w = 9, h = 5) {
  suppressPackageStartupMessages({
    library(data.table)
    library(fgsea)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
  })

  resis_genesets <- fgsea::gmtPathways(gmt_path)
  deg_df <- data.table::fread(deg_path)

  res <- run_fgsea(
    pathways = resis_genesets,
    deg_df = deg_df,
    deg_df_slpval = "sign_log_p"
  )

  plot <- plot_fgsea_bar(res, title = "Enrichment of Resistant Pathways")
  ggsave(output_plot, plot = plot, device = "pdf", width = w, height = h)
}

if (interactive() || identical(Sys.getenv("R_SCRIPT_DEBUG"), "TRUE")) {
  run_enrichment_and_plot(
    deg_path = "data/input/OV_recur_DESeq.txt",
    gmt_path = "data/signature_sets/resistance_signatures.gmt",
    output_plot = "results/fgsea_resistance_plot.pdf"
  )
} else {
  deg_path <- snakemake@input[[1]]
  gmt_path <- snakemake@input[[2]]
  output_plot <- snakemake@output[[1]]

  run_enrichment_and_plot(deg_path, gmt_path, output_plot)
}