#' Run Leiden Clustering on Seurat Object
run_leiden_clustering <- function(seurat_obj, norm_method = "SCT", resolution = 0.8,
                                  reduction = c("pca", "harmony"), n_pcs = 30, k_param = 20,
                                  run_umap = TRUE, run_tsne = FALSE, verbose = TRUE) {
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

#' Score Cells with Gene Signatures
load_default_signatures <- function() {
  list(
    SCN_Balanis = read.delim("data/internal/SCN_Balanis.txt")[, 1],
    RBLoss_Chen = read.delim("data/internal/RBLoss_Chen.txt")[, 1],
    RBLoss_Malorni = read.csv("data/internal/Malorni_RBsig_genes.csv")[, 1]
  )
}

score_cells_with_signatures <- function(seurat_obj, gene_sets, ctrl = 75) {
  for (gs_name in names(gene_sets)) {
    genes <- gene_sets[[gs_name]]
    seurat_obj <- AddModuleScore(
      object = seurat_obj,
      features = list(genes),
      ctrl = ctrl,
      search = TRUE,
      name = gs_name
    )
  }
  return(seurat_obj)
}

#' Wrapper: Cluster and Score
cluster_and_score_cells <- function(
  sobj,
  gene_sets = list(),
  resolution = 0.4,
  reduction = "pca",
  save_path = NULL
) {
  sobj <- run_leiden_clustering(sobj, resolution = resolution, reduction = reduction)
  if (length(gene_sets) == 0) {
  gene_sets <- load_default_signatures()
}
  sobj <- score_cells_with_signatures(sobj, gene_sets)

  if (!is.null(save_path)) {
    saveRDS(sobj, file = save_path)
    message("✅ Saved clustered + scored object to: ", save_path)
  }
  return(sobj)
}
