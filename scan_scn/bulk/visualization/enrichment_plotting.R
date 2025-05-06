
#' Barplot of NES for enriched pathways
plot_fgsea_bar <- function(fgsea_df, topn = 25, bottomn = 25, sort.val = NULL, title = "GSEA Enrichment") {
  fgsea_df_top <- fgsea_df[c(tail(order(fgsea_df$NES), topn), head(order(fgsea_df$NES), bottomn)), ]
  fgsea_df_top$logp <- abs(fgsea_df_top$signed_logp)

  ggpubr::ggbarplot(
    fgsea_df_top,
    x = "pathway",
    y = "NES",
    fill = "logp",
    color = NA,
    sort.by.groups = FALSE,
    sort.val = sort.val,
    title = title,
    ylab = "Normalized Enrichment Score (NES)",
    xlab = "Pathway",
    rotate = TRUE
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 12),
      axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.text = ggplot2::element_text(size = 12),
      legend.position = "right"
    ) +
    ggplot2::labs(fill = "-log10(pval)")
}

plot_fgsea_rubrary <- function(fgsea_df, NES_cutoff = 2, sig_cutoff = c("pval", 0.05), title = "GSEA") {
  Rubrary::plot_GSEA_barplot(
    gsea_res = fgsea_df,
    gsea_pws = fgsea_df$pathway,
    NES_cutoff = NES_cutoff,
    sig_cutoff = sig_cutoff,
    pw_size = 5,
    colors = c("firebrick", "darkblue"),
    title = title
  )
}
