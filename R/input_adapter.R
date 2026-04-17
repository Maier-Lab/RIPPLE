#' @title Input Adapter
#'
#' @description Internal function for normalizing input data from multiple
#'   canonical classes (Seurat, SingleCellExperiment, SpatialExperiment, or
#'   an `.rds` file containing one of these) into a common structure used by
#'   the RIPPLE pipeline.
#'
#' @name input_adapter
#' @keywords internal
NULL

#' Resolve RIPPLE input to normalized list
#'
#' Accepts any supported input format and returns a normalized list with
#' `counts`, `meta`, and optional `expr` components.
#'
#' Supported inputs:
#' \itemize{
#'   \item Path to an \code{.rds} file containing a Seurat, SingleCellExperiment,
#'     or SpatialExperiment object.
#'   \item In-memory Seurat object (with `@meta.data` containing coordinates).
#'   \item In-memory SingleCellExperiment (with coordinates in `colData`).
#'   \item In-memory SpatialExperiment (with coordinates in `spatialCoords()`).
#' }
#'
#' For CSV or plain matrix input, use \code{\link{make_ripple_input}} first to
#' build a canonical Seurat or SpatialExperiment object.
#'
#' @param input One of the supported input formats (see Details).
#' @param require_expr Logical. If TRUE, also resolves a normalized expression
#'   matrix (used by L-R integration). If the source object doesn't provide one,
#'   counts are log-normalized on the fly.
#' @param verbose Logical. Print progress messages (default: TRUE).
#'
#' @return A named list with components:
#' \describe{
#'   \item{counts}{Sparse or dense matrix (genes x cells) with raw integer counts.}
#'   \item{meta}{A \code{data.table} with cell metadata. Column \code{barcode}
#'     contains cell identifiers; other columns are user metadata fields.}
#'   \item{expr}{Normalized expression matrix (only if \code{require_expr = TRUE}).}
#' }
#'
#' @keywords internal
.resolve_input <- function(input, require_expr = FALSE, verbose = TRUE) {
  .msg <- function(...) if (isTRUE(verbose)) message(...)

  # ------------------------------------------------------------------
  # Case 1: Character path
  # ------------------------------------------------------------------
  if (is.character(input) && length(input) == 1) {
    if (!file.exists(input)) {
      stop("Input path does not exist: ", input, call. = FALSE)
    }
    .msg("Loading from file: ", input)
    obj <- readRDS(input)
    return(.resolve_input(obj, require_expr = require_expr, verbose = verbose))
  }

  # ------------------------------------------------------------------
  # Case 2: Seurat object
  # ------------------------------------------------------------------
  if (inherits(input, "Seurat")) {
    .msg("Detected Seurat object")
    counts <- .get_seurat_layer(input, "counts")
    meta <- data.table::as.data.table(input@meta.data, keep.rownames = "barcode")

    result <- list(counts = counts, meta = meta)

    if (require_expr) {
      # Suppress Seurat's "empty layer" warning; we fall back to log-normalizing
      expr <- suppressWarnings(tryCatch(
        .get_seurat_layer(input, "data"),
        error = function(e) NULL
      ))
      if (is.null(expr) || identical(dim(expr), c(0L, 0L)) ||
        (inherits(expr, "Matrix") && length(expr@x) == 0)) {
        .msg("No 'data' layer found; log-normalizing counts")
        expr <- .lognormalize(counts)
      }
      result$expr <- expr
    }

    return(result)
  }

  # ------------------------------------------------------------------
  # Case 3: SpatialExperiment (must check before SCE since it inherits)
  # ------------------------------------------------------------------
  if (inherits(input, "SpatialExperiment")) {
    if (!requireNamespace("SpatialExperiment", quietly = TRUE) ||
      !requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("Packages 'SpatialExperiment' and 'SummarizedExperiment' are ",
        "required to use SpatialExperiment input.\n",
        "Install with: BiocManager::install(c('SpatialExperiment', ",
        "'SummarizedExperiment'))",
        call. = FALSE
      )
    }
    .msg("Detected SpatialExperiment object")

    counts <- .get_sce_assay(input, "counts")
    meta <- .build_meta_from_sce(input, include_spatial_coords = TRUE)

    result <- list(counts = counts, meta = meta)

    if (require_expr) {
      assay_names <- SummarizedExperiment::assayNames(input)
      if ("logcounts" %in% assay_names) {
        result$expr <- SummarizedExperiment::assay(input, "logcounts")
      } else {
        .msg("No 'logcounts' assay found; log-normalizing counts")
        result$expr <- .lognormalize(counts)
      }
    }

    return(result)
  }

  # ------------------------------------------------------------------
  # Case 4: SingleCellExperiment (or any SummarizedExperiment)
  # ------------------------------------------------------------------
  if (inherits(input, "SingleCellExperiment") ||
    inherits(input, "SummarizedExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("Package 'SummarizedExperiment' is required to use SCE input.\n",
        "Install with: BiocManager::install('SummarizedExperiment')",
        call. = FALSE
      )
    }
    .msg("Detected SingleCellExperiment/SummarizedExperiment object")

    counts <- .get_sce_assay(input, "counts")
    meta <- .build_meta_from_sce(input, include_spatial_coords = FALSE)

    result <- list(counts = counts, meta = meta)

    if (require_expr) {
      assay_names <- SummarizedExperiment::assayNames(input)
      if ("logcounts" %in% assay_names) {
        result$expr <- SummarizedExperiment::assay(input, "logcounts")
      } else {
        .msg("No 'logcounts' assay found; log-normalizing counts")
        result$expr <- .lognormalize(counts)
      }
    }

    return(result)
  }

  # ------------------------------------------------------------------
  # Unsupported
  # ------------------------------------------------------------------
  stop("Unsupported input type: ", paste(class(input), collapse = ", "),
    "\nExpected one of: file path (.rds), Seurat, SingleCellExperiment, ",
    "or SpatialExperiment.\n",
    "For plain matrices/CSVs, use make_ripple_input() first.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------

#' Get a layer from a Seurat object, handling v4/v5 differences
#' @keywords internal
#' @noRd
.get_seurat_layer <- function(obj, layer) {
  tryCatch(
    Seurat::GetAssayData(obj, layer = layer),
    error = function(e) {
      tryCatch(
        Seurat::GetAssayData(obj, slot = layer),
        error = function(e2) {
          stop("Could not extract '", layer, "' from Seurat object. ",
            "Error: ", e$message,
            call. = FALSE
          )
        }
      )
    }
  )
}

#' Extract an assay from an SCE/SPE, preferring 'counts' if the name is missing
#' @keywords internal
#' @noRd
.get_sce_assay <- function(sce, preferred) {
  assay_names <- SummarizedExperiment::assayNames(sce)
  if (preferred %in% assay_names) {
    return(SummarizedExperiment::assay(sce, preferred))
  }
  # Fall back to first assay with a message
  message(
    "No '", preferred, "' assay found; using first assay: ",
    assay_names[1]
  )
  SummarizedExperiment::assay(sce, 1)
}

#' Build a metadata data.table from colData (and optionally spatialCoords)
#' @keywords internal
#' @noRd
.build_meta_from_sce <- function(sce, include_spatial_coords = FALSE) {
  col_data <- as.data.frame(SummarizedExperiment::colData(sce))

  # Ensure rownames are cell barcodes
  if (is.null(rownames(col_data)) ||
    all(rownames(col_data) == as.character(seq_len(nrow(col_data))))) {
    rownames(col_data) <- colnames(sce)
  }

  # Optionally fold in spatialCoords
  if (isTRUE(include_spatial_coords)) {
    sc <- tryCatch(
      SpatialExperiment::spatialCoords(sce),
      error = function(e) NULL
    )
    if (!is.null(sc) && ncol(sc) >= 2) {
      coord_df <- as.data.frame(sc)
      # Standardize column names if needed
      if (!any(c("x", "spatial_x", "x_centroid") %in% names(coord_df))) {
        names(coord_df)[1:2] <- c("x", "y")
      }
      col_data <- cbind(col_data, coord_df)
    }
  }

  meta <- data.table::as.data.table(col_data, keep.rownames = "barcode")
  meta
}

#' Log-normalize a counts matrix (equivalent to Seurat's LogNormalize)
#' @keywords internal
#' @noRd
.lognormalize <- function(counts, scale_factor = 1e4) {
  # Coerce Seurat Assay-like objects to a plain (sparse) matrix
  if (!is.matrix(counts) && !inherits(counts, "Matrix")) {
    counts <- tryCatch(
      methods::as(counts, "CsparseMatrix"),
      error = function(e) as.matrix(counts)
    )
  }
  total_per_cell <- Matrix::colSums(counts)
  total_per_cell[total_per_cell == 0] <- 1
  normed <- Matrix::t(Matrix::t(counts) / total_per_cell) * scale_factor
  log1p(normed)
}
