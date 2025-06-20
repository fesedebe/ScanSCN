#' Run clustering on an integrated Seurat object
#' @export
run_leiden_clustering <- function(seurat_obj, norm_method = "SCT", resolution = 0.8,
                                  reduction = c("pca", "harmony"), n_pcs = 30, k_param = 20,
                                  run_umap = TRUE, run_tsne = FALSE, verbose = TRUE) {

  # Scale data if Log Normalization was used                                  
  if(norm_method == "LogNormalize") {
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  }

  if (is.null(seurat_obj@reductions[[reduction]])) {
    seurat_obj <- RunPCA(seurat_obj, npcs = n_pcs, verbose = verbose)
  }

  seurat_obj <- FindNeighbors(seurat_obj, reduction = reduction, dims = 1:n_pcs, k.param = k_param)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution, algorithm = 4, method = "igraph")

  if (run_umap) {
    seurat_obj <- RunUMAP(seurat_obj, reduction = reduction, dims = 1:n_pcs)
  }
  if (run_tsne) {
    seurat_obj <- RunTSNE(seurat_obj, reduction = reduction, dims = 1:n_pcs)
  }

  return(seurat_obj)
}

#' Identify cluster or group markers using MAST
#' @export
find_cluster_markers <- function(seurat_obj, group.by = "seurat_clusters", assay = "SCT") {
  Idents(seurat_obj) <- group.by

  all_markers <- FindAllMarkers(
    object = PrepSCTFindMarkers(sobj),
    test.use = "MAST",
    min.pct = min.pct,
    assay = assay
  ) %>%
    mutate(p_val = ifelse(p_val == 0, .Machine$double.xmin, p_val),
           p_val_adj = ifelse(p_val_adj == 0, .Machine$double.xmin, p_val_adj),
           signed_log_p_val = sign(avg_log2FC) * -log10(p_val),
           signed_log_p_val_adj = sign(avg_log2FC) * -log10(p_val_adj))

  filtered <- all_markers %>%
    filter(p_val_adj < padj.threshold, abs(avg_log2FC) >= logfc.threshold)

  return(list(full = all_markers, filtered = filtered))
}

#' Wrapper: Cluster and Score
cluster_markers_and_score <- function(seurat_obj) {
  obj <- run_leiden_clustering(seurat_obj)
  markers <- find_cluster_markers(obj)
  sigs <- load_default_signatures()
  obj <- score_cells_with_signatures(obj, sigs)
  return(list(obj = obj, markers = markers))
}