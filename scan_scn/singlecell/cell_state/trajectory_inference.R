#' Convert Seurat object to CellDataSet (Monocle3)
#'
#' @param sobj A Seurat object
#' @param assay Assay to use (default: "SCT")
#' @return A Monocle3 CellDataSet object
#' @export
convert_to_cds <- function(sobj, assay = "SCT") {
  SeuratWrappers::as.cell_data_set(sobj, assay = assay)
}

#' Run Monocle3 trajectory inference
#'
#' @param cds A CellDataSet object
#' @return A CellDataSet object with pseudotime and trajectory graph
#' @export
run_monocle3_trajectory <- function(cds) {
  cds %>%
    monocle3::cluster_cells() %>%
    monocle3::learn_graph() %>%
    monocle3::order_cells()
}

#' Plot Monocle3 trajectory results
#'
#' @param cds A CellDataSet object with trajectory results
#' @return Named list of ggplot objects
#' @export
plot_monocle3_results <- function(cds) {
  list(
    umap = monocle3::plot_cells(cds, show_trajectory_graph = FALSE),
    graph = monocle3::plot_cells(
      cds,
      label_groups_by_cluster = FALSE,
      label_leaves = FALSE,
      label_branch_points = FALSE
    ) + ggplot2::labs(title = "Monocle3 Graph"),
    pseudotime = monocle3::plot_cells(
      cds,
      color_cells_by = "pseudotime",
      label_cell_groups = FALSE,
      label_leaves = FALSE,
      label_branch_points = TRUE,
      graph_label_size = 2
    ) + ggplot2::labs(title = "Monocle3 Pseudotime")
  )
}

#' Wrapper to run pseudotime inference using Monocle3 (default)
#'
#' @param sobj Seurat object
#' @param method Method to use: "monocle" (default)
#' @param save_path Path to save the RDS (optional)
#' @return List with CellDataSet and plots
#' @export
infer_pseudotime <- function(sobj, method = "monocle", save_path = NULL) {
  if (method != "monocle") stop("Only 'monocle' method supported for now.")

  cds <- convert_to_cds(sobj)
  cds <- run_monocle3_trajectory(cds)
  plots <- plot_monocle3_results(cds)

  if (!is.null(save_path)) {
    saveRDS(cds, save_path)
  }

  return(list(cds = cds, plots = plots))
}