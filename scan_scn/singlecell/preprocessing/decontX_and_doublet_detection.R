# Run DecontX
run_decontX_wrapper <- function(sce, sample_col) {
  sce <- celda::decontX(sce, batch = colData(sce)[[sample_col]])
  return(sce)
}

# Run scDblFinder on both raw and decontaminated SCEs
run_doublet_detection <- function(sce, sample_col) {
  sce <- scDblFinder::scDblFinder(sce, samples = sample_col)
  return(sce)
}

# Apply decontaminated counts to a new SCE
apply_decontaminated_counts <- function(sce, decontaminated_sce) {
  counts(sce) <- SummarizedExperiment::assay(decontaminated_sce, "decontXcounts")
  return(sce)
}