#' @title Data Loading Functions
#'
#' @description Functions for loading Seurat objects and extracting metadata
#'   for the RIPPLE spatial analysis pipeline.
#'
#' @name data_loading
NULL

#' Load just metadata and coordinates (fast, no expression matrix)
#'
#' Extracts cell metadata and spatial coordinates from a Seurat object without
#' keeping the full expression matrix in memory. This is much faster and uses
#' less RAM when only spatial information is needed.
#'
#' @param path Character. Path to Seurat RDS file.
#' @param celltype_column Character. Name of the metadata column containing cell
#'   type annotations (e.g., \code{"cell_type"}).
#' @param sample_column Character. Name of the metadata column containing
#'   sample/replicate IDs (default: \code{"sample_id"}).
#' @param x_column Character or NULL. Name of X coordinate column. If NULL,
#'   auto-detected via \code{\link{get_coord_columns}}.
#' @param y_column Character or NULL. Name of Y coordinate column. If NULL,
#'   auto-detected via \code{\link{get_coord_columns}}.
#'
#' @return A \code{data.table} with cell metadata including columns \code{barcode},
#'   \code{x}, \code{y}, and the cell type and sample columns.
#'
#' @details Loads the full Seurat object, extracts metadata and coordinates,
#'   then discards the object to free memory. Coordinate columns are resolved
#'   via \code{\link{get_coord_columns}} and aliased to \code{x} and \code{y}.
#'
#' @examples
#' \dontrun{
#' meta <- load_metadata_only(
#'   path = "/path/to/seurat.rds",
#'   celltype_column = "cell_type",
#'   sample_column = "sample_id"
#' )
#' }
#'
#' @importFrom data.table as.data.table setnames
#' @export
load_metadata_only <- function(path, celltype_column = "cell_type",
                               sample_column = "sample_id",
                               x_column = NULL, y_column = NULL) {
  message("Loading metadata from Seurat object...")

  if (!file.exists(path)) {
    stop("Seurat file not found at: ", path)
  }
  obj <- readRDS(path)

  # Extract metadata with barcodes

  meta <- data.table::as.data.table(obj@meta.data, keep.rownames = "barcode")

  # Create cell_id from barcode if not present (backward compatibility)
  if (!"cell_id" %in% names(meta)) {
    meta[, cell_id := barcode]
  }

  # Resolve coordinate columns and create x/y aliases
  coord_cols <- get_coord_columns(meta, x_col = x_column, y_col = y_column)
  if (coord_cols[1] != "x") meta[, x := get(coord_cols[1])]
  if (coord_cols[2] != "y") meta[, y := get(coord_cols[2])]

  # Verify cell type column exists
  if (!celltype_column %in% names(meta)) {
    stop(sprintf("Column '%s' not found in Seurat metadata.\n  Available columns: %s",
                 celltype_column,
                 paste(head(names(meta), 20), collapse = ", ")))
  }

  message(sprintf("\nUsing %s column for cell types", celltype_column))
  message(sprintf("Loaded metadata for %s cells", nrow(meta)))

  # Print cell type summary
  message(sprintf("\nCell type summary (using %s):", celltype_column))
  ct_summary <- meta[, .N, by = c(celltype_column)][order(-N)]
  print(ct_summary[1:min(15, nrow(ct_summary))])

  # Clean up Seurat object to free memory
  rm(obj)
  gc(verbose = FALSE)

  return(meta)
}


#' Check data exists and print summary
#'
#' Verifies that a Seurat RDS file exists at the given path and prints
#' a summary of its contents including cell count, gene count, and available
#' metadata columns.
#'
#' @param path Character. Path to the Seurat RDS file.
#' @param celltype_column Character. Name of the cell type annotation column to
#'   check for (default: \code{"cell_type"}).
#' @param sample_column Character. Name of the sample ID column to check for
#'   (default: \code{"sample_id"}).
#'
#' @return Invisible \code{TRUE} if file exists and loads correctly,
#'   invisible \code{FALSE} otherwise.
#'
#' @examples
#' \dontrun{
#' check_data("/path/to/seurat.rds", celltype_column = "cell_type")
#' }
#'
#' @export
check_data <- function(path, celltype_column = "cell_type",
                       sample_column = "sample_id") {
  message("Checking Seurat object at: ", path)
  message("")

  if (!file.exists(path)) {
    message("  File NOT FOUND")
    return(invisible(FALSE))
  }

  size_mb <- file.info(path)$size / (1024 * 1024)
  message(sprintf("  Found (%.1f MB)", size_mb))

  # Quick load to check contents
  message("  Loading to check contents...")
  obj <- readRDS(path)

  message(sprintf("  Cells: %s", format(ncol(obj), big.mark = ",")))
  message(sprintf("  Genes: %s", format(nrow(obj), big.mark = ",")))
  message("")

  message("  Metadata columns:")
  key_cols <- unique(c(celltype_column, sample_column))

  # Check common coordinate column names
  for (pair in list(c("spatial_x", "spatial_y"), c("x", "y"),
                    c("x_centroid", "y_centroid"))) {
    if (all(pair %in% colnames(obj@meta.data))) {
      key_cols <- c(key_cols, pair)
      break
    }
  }

  for (col in key_cols) {
    if (col %in% colnames(obj@meta.data)) {
      message(sprintf("    + %s", col))
    } else {
      message(sprintf("    - %s (missing)", col))
    }
  }

  rm(obj)
  gc(verbose = FALSE)

  return(invisible(TRUE))
}
