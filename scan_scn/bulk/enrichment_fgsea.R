
#' Run FGSEA (Fast Gene Set Enrichment Analysis)
#'
#' @param pathways A named list of gene sets
#' @param deg_df Data frame with at least gene and signed log p columns
#' @param deg_df_slpval Column name in deg_df with signed log p-values
#' @param deg_df_gene Column name with gene symbols (default: "gene")
#' @param subset Logical: subset by cluster? (default FALSE)
#' @param cluster_id Optional cluster label if subset is TRUE
#' @param minSize, maxSize GSEA parameters
#' @param nPermSimple Number of permutations
#'
#' @return A fgsea result data.frame
#' @export
run_fgsea <- function(
  pathways,
  deg_df,
  deg_df_slpval,
  deg_df_gene = "gene",
  subset = FALSE,
  cluster_id = NULL,
  minSize = 15,
  maxSize = 500,
  nPermSimple = 1000
) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("The 'fgsea' package is required but not installed.")
  }

  if (subset) {
    if (is.null(cluster_id)) {
      stop("You must provide 'cluster_id' when 'subset = TRUE'")
    }
    deg_df <- deg_df[deg_df$cluster == cluster_id, ]
  }

  stats <- setNames(
    deg_df[[deg_df_slpval]],
    deg_df[[deg_df_gene]]
  )

  set.seed(999)
  fgseaRes <- fgsea::fgsea(
    pathways = pathways,
    stats = stats,
    eps = 0.0,
    minSize = minSize,
    maxSize = maxSize,
    nPermSimple = nPermSimple
  ) |>
    dplyr::arrange(dplyr::desc(.data$NES)) |> # nolint
    dplyr::mutate(signed_logp = sign(.data$NES) * -log10(.data$pval)) |>
    tidyr::drop_na()

  fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, paste, collapse = ",")
  fgseaRes$NAME <- fgseaRes$pathway

  return(fgseaRes)
}

#' Write a gene set list to .gmt format
#'
#' @param gene_sets Named list of gene vectors
#' @param description_df Data frame with Signature and Publication columns
#' @param file Output file path
#' @export
write_gmt_with_desc <- function(gene_sets, description_df, file) {
  description_lookup <- stats::setNames(description_df$Publication, description_df$Signature)
  gmt_lines <- lapply(names(gene_sets), function(name) {
    genes <- gene_sets[[name]]
    desc <- description_lookup[[name]]
    if (is.null(desc)) desc <- "NA"
    line <- paste(c(name, desc, genes), collapse = "\t")
    return(line)
  })
  writeLines(unlist(gmt_lines), con = file)
}
