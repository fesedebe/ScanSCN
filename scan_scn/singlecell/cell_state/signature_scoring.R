#' Score gene signatures using Seurat::AddModuleScore
#'
#' @param sobj A Seurat object (with normalized data)
#' @param gene_sets A named list of character vectors (gene sets)
#' @param ctrl Number of control genes used in AddModuleScore (default: 75)
#' @return Seurat object with added module score columns
#' @export
score_signatures <- function(sobj, gene_sets, ctrl = 75) {
  stopifnot(is.list(gene_sets))
  
  sobj <- Seurat::AddModuleScore(
    object = sobj,
    features = gene_sets,
    ctrl = ctrl,
    search = TRUE,
    name = names(gene_sets)
  )
  
  return(sobj)
}

#' Load default SCN-related gene sets for scoring
#'
#' @return Named list of gene vectors
#' @export
load_default_signature_sets <- function() {
  list(
    RBLoss_Malorni = read.csv("data/internal/RBLoss_Malorni.csv", header = TRUE)[, 1],
    RBLoss_Chen = read.delim("data/internal/RBLoss_Chen.txt", header = TRUE)[, 1],
    SCN_Balanis = read.delim("data/internal/SCN_Balanis.txt", header = TRUE) |>
      dplyr::arrange(dplyr::desc(PC1.v)) |>
      dplyr::slice_head(n = 500) |>
      dplyr::pull(Gene)
  )
}