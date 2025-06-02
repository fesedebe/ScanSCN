# Calculate nCount, nFeature, and percent.mt
add_qc_metrics <- function(seurat_obj) {
  seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  return(seurat_obj)
}

# Outlier flagging using MAD
flag_outliers_mad <- function(metric, batch, nmads = 3, type = "higher") {
  isOutlier(
    metric = metric,
    batch = batch,
    nmads = nmads,
    type = type,
    log = TRUE
  )
}