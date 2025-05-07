
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(cowplot)
})

utils::globalVariables(c(
  "label", "p_value"
))

# ---- Internal utility functions ----

combine_pvalues_fisher <- function(p_values) {
  chisq_stat <- -2 * sum(log(p_values))
  pchisq(chisq_stat, df = 2 * length(p_values), lower.tail = FALSE)
}

build_theme <- function(text_size, legend.position, legend.justification) {
  cowplot::theme_cowplot() +
    ggplot2::theme(
      text = element_text(size = text_size),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(size = 18),
      strip.text = element_text(size = 16),
      legend.position = legend.position,
      legend.justification = legend.justification,
      legend.title = element_text(size = 15, face = "bold"),
      legend.text = element_text(size = text_size),
      axis.text.y = element_text(size = text_size),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )
}

add_combined_pvalue_label <- function(p, data, x, y, facet_by, p_anno_x, p_anno_y, p_anno_lab, p_anno_n, p_paired) {
  p_values <- data |>
    group_by(across(all_of(facet_by))) |>
    summarize(p_value = wilcox.test(.data[[y]] ~ .data[[x]], paired = p_paired)$p.value, .groups = "drop") |>
    pull(p_value) # nolint

  combined_p_value <- format(combine_pvalues_fisher(p_values), scientific = TRUE, digits = 3)

  if (is.null(p_anno_lab)) {
    p_anno_lab <- unique(data[[facet_by]])
  }

  anno_data <- data.frame(
    Dataset = p_anno_lab,
    label = c(paste0("Combined p = ", combined_p_value), rep(" ", p_anno_n))
  )

  p + geom_text(data = anno_data, mapping = aes(x = p_anno_x, y = p_anno_y, label = label), size = 4) # nolint
}

add_plot_tag <- function(p, plot_tag) {
  p +
    labs(tag = plot_tag) +
    theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0.01, 0.98)
    )
}

# ---- Main exported function ----

#' Create a Paired Boxplot with Optional Statistical Annotations
#' @export
create_paired_boxplot <- function(
  data, x, y,
  xlab = NULL,
  ylab = "Signature Score",
  title = NULL,
  facet_by = "Dataset",
  color = "black",
  fill,
  line_color = "gray",
  line_size = 0.4,
  text_size = 14.5,
  legend.position = "top",
  legend.justification = "center",
  p_anno_lab = NULL,
  p_anno_n = 1,
  p_anno_x = 1.1,
  p_anno_y = 75,
  p_label_y = NULL,
  ggtheme = NULL,
  palette = c("#0072B2", "#D55E00"),
  y_limits = NULL,
  y_breaks = NULL,
  p_method = "wilcox.test",
  p_label = "p.signif",
  p_label_x = 1.3,
  p_paired = TRUE,
  plot_tag = NULL
) {
  if (is.null(ggtheme)) {
    ggtheme <- build_theme(text_size, legend.position, legend.justification)
  }

  p <- ggpubr::ggpaired(
    data = data,
    x = x,
    y = y,
    xlab = xlab,
    ylab = ylab,
    title = title,
    facet.by = facet_by,
    color = color,
    fill = fill,
    line.color = line_color,
    line.size = line_size,
    ggtheme = ggtheme,
    palette = palette
  ) + labs(fill = "Treatment Status")

  if (p_method != "combined") {
    p <- p + ggpubr::stat_compare_means(
      paired = p_paired,
      label = p_label,
      label.x = p_label_x,
      label.y = p_label_y,
      method = p_method
    )
  } else {
    p <- add_combined_pvalue_label(p, data, x, y, facet_by, p_anno_x, p_anno_y, p_anno_lab, p_anno_n, p_paired)
  }

  if (!is.null(y_limits) || !is.null(y_breaks)) {
    p <- p + scale_y_continuous(limits = y_limits, breaks = y_breaks)
  }

  if (!is.null(plot_tag)) {
    p <- add_plot_tag(p, plot_tag)
  }

  return(p)
}