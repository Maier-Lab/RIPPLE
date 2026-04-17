#' @title Data Loading Functions
#'
#' @description Functions for loading Seurat objects and extracting metadata
#'   for the RIPPLE spatial analysis pipeline.
#'
#' @name data_loading
NULL

#' Load just metadata and coordinates (fast, no expression matrix)
#'
#' Extracts cell metadata and spatial coordinates from a RIPPLE-compatible
#' input (Seurat, SingleCellExperiment, SpatialExperiment, or an `.rds` file
#' containing one of these) without keeping the full expression matrix in
#' memory.
#'
#' @param input Either a path to an `.rds` file or an in-memory Seurat,
#'   SingleCellExperiment, or SpatialExperiment object.
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
#' @examples
#' \dontrun{
#' meta <- load_metadata_only(
#'   input = "/path/to/seurat.rds",
#'   celltype_column = "cell_type",
#'   sample_column = "sample_id"
#' )
#' }
#'
#' @importFrom data.table as.data.table setnames
#' @export
load_metadata_only <- function(input, celltype_column = "cell_type",
                               sample_column = "sample_id",
                               x_column = NULL, y_column = NULL) {
  message("Loading metadata...")

  data <- .resolve_input(input, require_expr = FALSE, verbose = FALSE)
  meta <- data$meta
  rm(data)

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
    stop(sprintf(
      "Column '%s' not found in metadata.\n  Available columns: %s",
      celltype_column,
      paste(head(names(meta), 20), collapse = ", ")
    ))
  }

  message(sprintf("\nUsing %s column for cell types", celltype_column))
  message(sprintf("Loaded metadata for %s cells", nrow(meta)))

  message(sprintf("\nCell type summary (using %s):", celltype_column))
  ct_summary <- meta[, .N, by = c(celltype_column)][order(-N)]
  print(ct_summary[1:min(15, nrow(ct_summary))])

  gc(verbose = FALSE)
  return(meta)
}


#' Check input data and print summary
#'
#' Loads a RIPPLE-compatible input (Seurat, SCE, SpatialExperiment, or an
#' `.rds` file containing one of these) and prints a summary of its contents
#' including cell/gene counts and availability of key metadata columns.
#'
#' @param input Either a path to an `.rds` file or an in-memory object.
#' @param celltype_column Character. Name of the cell type annotation column to
#'   check for (default: \code{"cell_type"}).
#' @param sample_column Character. Name of the sample ID column to check for
#'   (default: \code{"sample_id"}).
#'
#' @return Invisible \code{TRUE} if loadable, invisible \code{FALSE} otherwise.
#'
#' @examples
#' \dontrun{
#' check_data("/path/to/data.rds", celltype_column = "cell_type")
#' check_data(my_spe, celltype_column = "cell_type")
#' }
#'
#' @export
check_data <- function(input, celltype_column = "cell_type",
                       sample_column = "sample_id") {
  if (is.character(input) && length(input) == 1) {
    message("Checking input at: ", input)
    if (!file.exists(input)) {
      message("  File NOT FOUND")
      return(invisible(FALSE))
    }
    size_mb <- file.info(input)$size / (1024 * 1024)
    message(sprintf("  Found (%.1f MB)", size_mb))
  } else {
    message(
      "Checking in-memory input of class: ",
      paste(class(input), collapse = ", ")
    )
  }

  data <- tryCatch(
    .resolve_input(input, require_expr = FALSE, verbose = FALSE),
    error = function(e) {
      message("  Failed to load: ", e$message)
      NULL
    }
  )
  if (is.null(data)) {
    return(invisible(FALSE))
  }

  message(sprintf("  Cells: %s", format(ncol(data$counts), big.mark = ",")))
  message(sprintf("  Genes: %s", format(nrow(data$counts), big.mark = ",")))
  message("")

  message("  Metadata columns:")
  meta_names <- names(data$meta)
  key_cols <- unique(c(celltype_column, sample_column))

  # Check common coordinate column names
  for (pair in list(
    c("spatial_x", "spatial_y"), c("x", "y"),
    c("x_centroid", "y_centroid")
  )) {
    if (all(pair %in% meta_names)) {
      key_cols <- c(key_cols, pair)
      break
    }
  }

  for (col in key_cols) {
    if (col %in% meta_names) {
      message(sprintf("    + %s", col))
    } else {
      message(sprintf("    - %s (missing)", col))
    }
  }

  gc(verbose = FALSE)
  return(invisible(TRUE))
}
