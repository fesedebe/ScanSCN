# Load 10X matrix into Seurat object
load_raw_10X_matrix <- function(matrix_dir, project, min_cells = 3, min_features = 200) {
  sobj <- Seurat::CreateSeuratObject(
    counts = Seurat::Read10X(matrix_dir),
    project = project,
    min.cells = min_cells,
    min.features = min_features
  )
  return(sobj)
}

# Generate annotation from CellRanger aggregation CSV
create_barcode_sample_annotation <- function(seurat_obj, aggr_csv) {
  aggr_anno <- read.csv(aggr_csv)
  barcodes <- Seurat::Cells(seurat_obj)
  barcode_parts <- stringr::str_split(barcodes, "-", simplify = TRUE)
  barcode_df <- data.frame(
    barcode = barcodes,
    id = barcode_parts[, 2],
    read = barcode_parts[, 1],
    sample = aggr_anno$sample_id[match(barcode_parts[, 2], rownames(aggr_anno))],
    stringsAsFactors = FALSE
  )
  return(barcode_df)
}

# Add sample metadata to Seurat object
add_sample_metadata <- function(seurat_obj, barcode_df) {
  rownames(barcode_df) <- barcode_df$barcode
  seurat_obj <- Seurat::AddMetaData(seurat_obj, barcode_df[, c("sample", "read")])
  return(seurat_obj)
}