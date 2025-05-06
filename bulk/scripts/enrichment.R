args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript enrichment.R <DEG_table> <geneset_gmt> <output_plot>")
}

deg_path <- args[1]
gmt_path <- args[2]
output_plot <- args[3]

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

print("✅ Starting enrichment.R")

# Load functions & input files
source("bulk/scripts/fgsea_utils.R")
source("bulk/scripts/fgsea_plotting.R")

resis_genesets <- fgsea::gmtPathways(gmt_path)
deg_df <- data.table::fread(deg_path)

# Run enrichment
res <- run_fgsea(
  pathways = resis_genesets,
  deg_df = deg_df,
  deg_df_slpval = "sign_log_p"
)

# Plot and save
plot <- plot_fgsea_bar(fgsea_results, title = "Enrichment of Resistant Pathways")
ggsave(filename = output_plot, plot = plot, device = "pdf", width = 7.5, height = 5.8)

