#' Label SCN-high cells within each cluster
#'
#' Adds a logical metadata column "SCNhicells" to the Seurat object based on SCN scores.
#'
#' @param seurat_obj A Seurat object with SCN scores in metadata (e.g. SCNBalanis_AUC)
#' @param score_col Metadata column with SCN score (default: "SCNBalanis_AUC")
#' @param threshold_mode One of "percentile", "MAD", or "value"
#' @param value_thresh Fixed threshold if using "value" mode
#' @param thresh_percentile Top percentile if using "percentile" mode
#'
#' @return Seurat object with updated "SCNhicells" column
label_scnhigh_cells <- function(seurat_obj, 
                                score_col = "SCNBalanis_AUC",
                                threshold_mode = "value",
                                value_thresh = 0.03,
                                thresh_percentile = 10) {
  seurat_obj$SCNhicells <- FALSE
  Idents(seurat_obj) <- "seurat_clusters"
  clusters <- levels(Idents(seurat_obj))

  for (clust in clusters) {
    cells <- WhichCells(seurat_obj, idents = clust)
    scores <- seurat_obj[[score_col]][cells, 1]

    threshold <- switch(
      threshold_mode,
      percentile = quantile(scores, probs = 1 - thresh_percentile / 100),
      MAD        = median(scores) + 4 * mad(scores),
      value      = value_thresh,
      stop("Invalid threshold_mode: use 'percentile', 'MAD', or 'value'")
    )

    scn_hi <- cells[scores > threshold]
    seurat_obj$SCNhicells[scn_hi] <- TRUE
  }

  return(seurat_obj)
}


#' Aggregate pseudobulk counts by group (cluster + SCN status)
#'
#' @param seurat_obj A Seurat object with RNA assay
#' @param group_col Column to define groups (e.g., "group" with cluster + SCN status)
#' 
#' @return List with: matrix of pseudobulk counts, sample_table
aggregate_pseudobulk_counts <- function(seurat_obj, group_col = "group") {
  seurat_obj$group <- paste0(seurat_obj$seurat_clusters, "_", seurat_obj$SCNhicells)

  metadata <- seurat_obj@meta.data
  counts <- Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  groups <- unique(metadata[[group_col]])

  pseudobulk_mat <- sapply(groups, function(grp) {
    cells <- rownames(metadata[metadata[[group_col]] == grp, ])
    Matrix::rowSums(counts[, cells, drop = FALSE])
  })

  sample_table <- data.frame(
    group = groups,
    cluster = sapply(strsplit(groups, "_"), `[`, 1),
    SCNhicells = factor(sapply(strsplit(groups, "_"), `[`, 2), levels = c("FALSE", "TRUE"))
  )

  return(list(counts = pseudobulk_mat, metadata = sample_table))
}


#' Run DESeq2 to identify genes differentially expressed by SCN-high status
#'
#' @param pseudobulk_counts Matrix of pseudobulk counts
#' @param sample_table Metadata dataframe with 'SCNhicells' and 'cluster'
#'
#' @return DESeq2 results as a data.frame with logFC and signed log10 p-adj
run_deseq2_on_pseudobulk <- function(pseudobulk_counts, sample_table) {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = pseudobulk_counts,
    colData = sample_table,
    design = ~ cluster + SCNhicells
  )

  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = c("SCNhicells", "TRUE", "FALSE"))

  result_df <- as.data.frame(res) %>%
    dplyr::mutate(
      gene = rownames(.),
      signedlogpadj = -log10(padj + 1e-300) * sign(log2FoldChange)
    ) %>%
    dplyr::arrange(desc(signedlogpadj))

  return(result_df)
}