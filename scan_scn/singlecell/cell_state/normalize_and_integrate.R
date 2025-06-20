# normalize_and_integrate.R
# Contains functions for scRNA-seq normalization, cell cycle scoring, and integration

#' Score Cell Cycle
score_cell_cycle <- function(sobj) {
  Seurat::CellCycleScoring(
    sobj,
    s.features = Seurat::cc.genes$s.genes,
    g2m.features = Seurat::cc.genes$g2m.genes
  )
}

#' Split and SCTransform Normalize Each Sample
sctransform_by_sample <- function(sobj, vars.to.regress) {
  sobj.list <- Seurat::SplitObject(sobj, split.by = "sample")
  sobj.list <- lapply(sobj.list, function(x) {
    SCTransform(x, vars.to.regress = vars.to.regress, verbose = FALSE)
  })
  return(sobj.list)
}

#' Run RPCA Integration
integrate_with_rpca <- function(sobj.list, nfeatures = 2000) {
  features <- SelectIntegrationFeatures(sobj.list, nfeatures = nfeatures)
  sobj.list <- PrepSCTIntegration(sobj.list, anchor.features = features)
  sobj.list <- lapply(sobj.list, RunPCA, features = features)
  anchors <- FindIntegrationAnchors(
    object.list = sobj.list,
    anchor.features = features,
    normalization.method = "SCT",
    reduction = "rpca"
  )
  integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  DefaultAssay(integrated) <- "integrated"
  return(integrated)
}

#' Run Harmony Integration
integrate_with_harmony <- function(sobj, vars.to.regress) {
  sobj.list <- Seurat::SplitObject(sobj, split.by = "sample")
  merged <- Reduce(function(x, y) merge(x, y), sobj.list)
  merged <- SCTransform(merged, vars.to.regress = vars.to.regress, verbose = FALSE)
  merged <- RunPCA(merged)
  merged <- harmony::RunHarmony(merged, group.by.vars = "sample")
  return(merged)
}

#' Normalize and Integrate scRNA-seq Data (Wrapper)
#'
#' @export
normalize_and_integrate <- function(
  sobj,
  method = "rpca",
  vars.to.regress = c("S.Score", "G2M.Score"),
  nfeatures = 2000
) {
  stopifnot("sample" %in% colnames(sobj@meta.data))

  sobj <- score_cell_cycle(sobj)

  if (method == "rpca") {
    sobj.list <- sctransform_by_sample(sobj, vars.to.regress)
    integrated <- integrate_with_rpca(sobj.list, nfeatures)
    return(integrated)

  } else if (method == "harmony") {
    integrated <- integrate_with_harmony(sobj, vars.to.regress)
    return(integrated)

  } else {
    stop("Unsupported method: choose 'rpca' or 'harmony'")
  }
}