
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  source("scan_scn/bulk/visualization/plot_paired_boxplot.R")
})

run_plot_scn_scores <- function(scores_path, output_plot, w = 6, h = 5) {
  df <- read.csv(scores_path, sep = "\t", header = TRUE)

  df$Treatment_Status <- as.factor(df$Treatment_Status)
  df$Dataset <- as.factor(df$Dataset)

  p <- create_paired_boxplot(
    data = df,
    x = "Treatment_Status",
    y = "SCN_Score",
    title = "Change in SCN Score",
    fill = "Treatment_Status",
    facet_by = "Dataset"
  )

  ggsave(output_plot, plot = p, width = w, height = h)
  message(paste("Plot saved to", output_plot))
}

# interactive or Snakemake
if (interactive() || identical(Sys.getenv("R_SCRIPT_DEBUG"), "TRUE")) {
  run_plot_scn_scores(
    scores_path = "data/input/OV_SCNscores.txt",
    output_plot = "res/test_scn_score_plot.pdf"
  )
} else {
  scores_path <- snakemake@input[[1]]
  output_plot <- snakemake@output[[1]]
  run_plot_scn_scores(scores_path, output_plot)
}
