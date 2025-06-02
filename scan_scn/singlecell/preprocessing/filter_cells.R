# Filter cells based on multiple QC metrics
filter_cells <- function(seurat_obj) {
  subset(
    seurat_obj,
    subset = scDblFinder.class == "singlet" &
             MT.3MADs == FALSE &
             percent.mt < 25 &
             nCount_RNA > 500 &
             umi.mad.hi == FALSE &
             nFeature_RNA > 200 &
             feature.mad.hi == FALSE
  )
}

# Save Seurat object
save_seurat_object <- function(seurat_obj, path) {
  saveRDS(seurat_obj, file = path)
}