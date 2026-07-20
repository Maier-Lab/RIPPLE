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
.resolve_input <- function(input, require_expr = FALSE, verbose = TRUE,
                           assay = NULL) {
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
    return(.resolve_input(
      obj, require_expr = require_expr, verbose = verbose, assay = assay
    ))
  }

  # ------------------------------------------------------------------
  # Case 2: Seurat object
  # ------------------------------------------------------------------
  if (inherits(input, "Seurat")) {
    .msg("Detected Seurat object")
    chosen_assay <- .resolve_seurat_assay(input, assay = assay, verbose = verbose)
    counts <- .get_seurat_layer(input, "counts", assay = chosen_assay)
    .check_integer_counts(counts, chosen_assay)
    meta_df <- .inject_seurat_fov_coords(input, verbose = verbose)
    meta <- data.table::as.data.table(meta_df, keep.rownames = "barcode")

    result <- list(counts = counts, meta = meta)

    if (require_expr) {
      # Suppress Seurat's "empty layer" warning; we fall back to log-normalizing
      expr <- suppressWarnings(tryCatch(
        .get_seurat_layer(input, "data", assay = chosen_assay),
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
.get_seurat_layer <- function(obj, layer, assay = NULL) {
  args_v5 <- list(object = obj, layer = layer)
  args_v4 <- list(object = obj, slot = layer)
  if (!is.null(assay)) {
    args_v5$assay <- assay
    args_v4$assay <- assay
  }
  tryCatch(
    do.call(Seurat::GetAssayData, args_v5),
    error = function(e) {
      tryCatch(
        do.call(Seurat::GetAssayData, args_v4),
        error = function(e2) {
          stop("Could not extract '", layer, "' from Seurat object",
            if (!is.null(assay)) paste0(" (assay '", assay, "')") else "",
            ". Error: ", e$message,
            call. = FALSE
          )
        }
      )
    }
  )
}

#' Pick which Seurat assay to use for RIPPLE
#'
#' Resolution order:
#' 1. If \code{assay} is explicitly provided, validate it exists.
#' 2. Otherwise try known raw-count assay names in order: Xenium, RNA, Spatial.
#'    First match wins.
#' 3. Otherwise fall back to the object's active assay, with a warning.
#'
#' This avoids the common trap where SCTransform (or another normalization)
#' has taken over as the active assay, and RIPPLE silently gets normalized
#' values instead of raw counts.
#'
#' @keywords internal
#' @noRd
.resolve_seurat_assay <- function(obj, assay = NULL, verbose = TRUE) {
  .msg <- function(...) if (isTRUE(verbose)) message(...)
  available <- Seurat::Assays(obj)

  if (!is.null(assay)) {
    if (!assay %in% available) {
      stop(
        "Assay '", assay, "' not found. Available: ",
        paste(available, collapse = ", "),
        call. = FALSE
      )
    }
    .msg("  Using assay: '", assay, "' (user-specified)")
    return(assay)
  }

  priority <- c("Xenium", "RNA", "Spatial")
  match <- priority[priority %in% available]
  if (length(match) > 0) {
    chosen <- match[1]
    .msg("  Using assay: '", chosen, "' (auto-selected raw-count assay)")
    return(chosen)
  }

  active <- Seurat::DefaultAssay(obj)
  warning(
    "No standard raw-count assay found (looked for Xenium, RNA, Spatial). ",
    "Falling back to the active assay '", active, "'. ",
    "If this is normalized data (e.g. SCT), RIPPLE's Poisson GLM will ",
    "misbehave. Pass assay = '<name>' to run_ripple(), or set ",
    "DefaultAssay(obj) <- '<name>' before calling. ",
    "Available assays: ", paste(available, collapse = ", "),
    call. = FALSE
  )
  active
}

#' Extract centroids from Seurat FOV objects into @meta.data if missing
#'
#' Xenium (and MERSCOPE / CosMx via Seurat's FOV class) data loaded via
#' \code{LoadXenium()} keeps per-cell spatial centroids inside FOV objects
#' at \code{obj@@images$fov}, not in \code{@@meta.data}. RIPPLE's
#' coordinate resolver (\code{\link{get_coord_columns}}) only reads
#' \code{@@meta.data}, so out of the box a Xenium Seurat object errors
#' with "Could not find spatial coordinate columns".
#'
#' This helper closes that gap. If the metadata already has any of the
#' recognised coordinate column names (x/y, x_centroid/y_centroid,
#' spatial_x/spatial_y), it does nothing. Otherwise it iterates every
#' FOV-class image, calls \code{Seurat::GetTissueCoordinates()} on each,
#' row-binds the results, aligns rows to the metadata's cell order, and
#' adds \code{x_centroid} / \code{y_centroid} columns. A one-line
#' message reports what was done.
#'
#' Never overwrites existing metadata coordinate columns.
#'
#' @keywords internal
#' @noRd
.inject_seurat_fov_coords <- function(obj, verbose = TRUE) {
  meta <- obj@meta.data

  # Already has coordinates in meta.data — respect the user's setup.
  known <- c("x", "y", "x_centroid", "y_centroid", "spatial_x", "spatial_y")
  if (any(known %in% names(meta))) {
    return(meta)
  }

  img_names <- tryCatch(Seurat::Images(obj), error = function(e) character())
  if (!length(img_names)) return(meta)
  fov_names <- Filter(function(nm) inherits(obj[[nm]], "FOV"), img_names)
  if (!length(fov_names)) return(meta)

  coord_list <- lapply(fov_names, function(nm) {
    tryCatch(
      Seurat::GetTissueCoordinates(obj[[nm]]),
      error = function(e) NULL
    )
  })
  coord_list <- Filter(Negate(is.null), coord_list)
  if (!length(coord_list)) return(meta)

  all_coords <- do.call(rbind, coord_list)

  # GetTissueCoordinates returns a data.frame with x/y numeric cols and
  # (Seurat >= 4) a "cell" column carrying the barcode. Older versions
  # use rownames instead.
  cell_col <- if ("cell" %in% names(all_coords)) {
    as.character(all_coords$cell)
  } else {
    rownames(all_coords)
  }

  # Align to meta.data cell order; drop cells with no FOV entry.
  match_idx <- match(rownames(meta), cell_col)
  if (all(is.na(match_idx))) {
    warning(
      "FOV centroids present but none of the barcodes matched @meta.data. ",
      "Skipping FOV coordinate extraction. Add x_centroid / y_centroid to ",
      "@meta.data manually.",
      call. = FALSE
    )
    return(meta)
  }
  meta$x_centroid <- all_coords$x[match_idx]
  meta$y_centroid <- all_coords$y[match_idx]

  n_missing <- sum(is.na(meta$x_centroid))
  if (isTRUE(verbose)) {
    message(
      "  Extracted centroids from ", length(fov_names),
      " FOV(s) (", paste(fov_names, collapse = ", "),
      ") into x_centroid / y_centroid",
      if (n_missing) paste0(" (", n_missing, " cells missing coords)") else ""
    )
  }
  meta
}

#' Warn if the counts matrix contains non-integer values
#'
#' RIPPLE fits a Poisson GLM with an offset for cell size, which assumes
#' raw integer counts. SCTransform, normalization, and imputation all
#' produce non-integer values that break the assumption silently.
#'
#' Samples the first \code{n_sample} non-zero values rather than the whole
#' matrix so the check is O(1) even on 100M-cell inputs.
#'
#' @keywords internal
#' @noRd
.check_integer_counts <- function(counts, assay_name = NULL,
                                  n_sample = 1000) {
  if (!nrow(counts) || !ncol(counts)) return(invisible(NULL))
  vals <- if (inherits(counts, "Matrix") && .hasSlot(counts, "x")) {
    utils::head(counts@x, n_sample)
  } else {
    utils::head(as.numeric(counts), n_sample)
  }
  vals <- vals[is.finite(vals) & vals != 0]
  if (!length(vals)) return(invisible(NULL))
  if (!isTRUE(all.equal(vals, round(vals), tolerance = 1e-8))) {
    warning(
      "Counts matrix ",
      if (!is.null(assay_name)) paste0("(assay '", assay_name, "') "),
      "contains non-integer values. RIPPLE fits a Poisson GLM and expects ",
      "raw integer counts. Common cause: SCTransform, LogNormalize, or ",
      "another normalization has replaced the counts. Pass assay = '<raw ",
      "assay>' to run_ripple(), or set DefaultAssay(obj) <- '<raw assay>' ",
      "before calling.",
      call. = FALSE
    )
  }
  invisible(NULL)
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
