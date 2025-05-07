#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  source("scan_scn/bulk/visualization/plot_paired_boxplot.R")
})

# ---- Main block ----
if (interactive() || identical(Sys.getenv("R_SCRIPT_DEBUG"), "TRUE") || !interactive()) {
  
  # Load test data
  df <- read.csv("data/input/OV_SCNscores.txt", sep = "\t", header = TRUE)

  # Optional: ensure correct column types
  df$Treatment_Status <- as.factor(df$Treatment_Status)
  df$Dataset <- as.factor(df$Dataset)

  # Call plot function
  p <- create_paired_boxplot(
    data = df,
    x = "Treatment_Status",
    y = "SCN_Score",
    title = "Change in SCN Score",
    fill = "Treatment_Status",
    facet_by = "Dataset"
  )

  # Save the plot
  ggsave("res/test_scn_score_plot.pdf", p, width = 6, height = 5)
  message("Plot saved to res/test_scn_score_plot.pdf")
}