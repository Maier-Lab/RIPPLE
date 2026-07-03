#' @title Input Constructors
#'
#' @description User-facing constructors for building a RIPPLE-compatible
#'   canonical object (Seurat or SpatialExperiment) from raw matrices or CSVs.
#'
#' @name make_input
NULL

#' Build a RIPPLE input object from raw matrices
#'
#' Converts a counts matrix plus metadata (and optional coordinates) into a
#' canonical Seurat or SpatialExperiment object ready to be passed to
#' \code{\link{run_ripple}}.
#'
#' @param counts Matrix (sparse or dense). Genes in rows, cells in columns.
#'   Must have rownames (gene names) and colnames (cell barcodes).
#' @param metadata data.frame or data.table with one row per cell. Rownames
#'   (or a \code{barcode}/\code{cell_id} column) must match \code{colnames(counts)}.
#' @param coords Optional. Matrix or data.frame of spatial coordinates with
#'   one row per cell. Must have two columns (x, y). Rownames must match cells,
#'   or the order must match \code{colnames(counts)}.
#' @param output_class Character. Either \code{"SpatialExperiment"} (default)
#'   or \code{"Seurat"}.
#' @param x_column Character. Name to assign to the x coordinate column in
#'   the output (default: \code{"x"}).
#' @param y_column Character. Name to assign to the y coordinate column in
#'   the output (default: \code{"y"}).
#'
#' @return A \code{SpatialExperiment} or \code{Seurat} object containing the
#'   counts, metadata, and spatial coordinates.
#'
#' @details
#' For \code{output_class = "SpatialExperiment"}: the coordinates are stored
#' in the \code{spatialCoords()} slot; non-coordinate metadata is stored in
#' \code{colData()}. Requires the \pkg{SpatialExperiment} Bioconductor package.
#'
#' For \code{output_class = "Seurat"}: the coordinates are added as metadata
#' columns (default names \code{x} and \code{y}); other metadata becomes the
#' Seurat object's \code{@meta.data}.
#'
#' @examples
#' \dontrun{
#' # Build a SpatialExperiment from raw matrices
#' spe <- make_ripple_input(
#'   counts = my_counts, # matrix or dgCMatrix
#'   metadata = my_metadata, # data.frame with sample_id, cell_type, etc.
#'   coords = my_coords, # two-column matrix (x, y)
#'   output_class = "SpatialExperiment"
#' )
#'
#' # Now run RIPPLE
#' results <- run_ripple(
#'   input           = spe,
#'   query_celltype  = "Tumor",
#'   celltype_column = "cell_type"
#' )
#' }
#'
#' @importFrom data.table as.data.table setnames setDT
#' @export
make_ripple_input <- function(counts,
                              metadata,
                              coords = NULL,
                              output_class = c("SpatialExperiment", "Seurat"),
                              x_column = "x",
                              y_column = "y") {
  output_class <- match.arg(output_class)

  # ------------------------------------------------------------------
  # 1. Validate counts
  # ------------------------------------------------------------------
  if (is.null(rownames(counts))) {
    stop("counts must have rownames (gene names)", call. = FALSE)
  }
  if (is.null(colnames(counts))) {
    stop("counts must have colnames (cell barcodes)", call. = FALSE)
  }

  # ------------------------------------------------------------------
  # 2. Normalize metadata (ensure rownames are cell barcodes)
  # ------------------------------------------------------------------
  if (inherits(metadata, "data.table")) {
    meta_df <- as.data.frame(metadata)
    if ("barcode" %in% names(meta_df)) {
      rownames(meta_df) <- meta_df$barcode
      meta_df$barcode <- NULL
    } else if ("cell_id" %in% names(meta_df)) {
      rownames(meta_df) <- meta_df$cell_id
      meta_df$cell_id <- NULL
    }
  } else if (is.data.frame(metadata)) {
    meta_df <- metadata
  } else {
    stop("metadata must be a data.frame or data.table", call. = FALSE)
  }

  if (is.null(rownames(meta_df))) {
    stop("metadata must have rownames (or a 'barcode'/'cell_id' column) ",
      "matching counts column names",
      call. = FALSE
    )
  }

  # ------------------------------------------------------------------
  # 3. Align cells between counts and metadata
  # ------------------------------------------------------------------
  common <- intersect(colnames(counts), rownames(meta_df))
  if (length(common) == 0) {
    stop("No matching barcodes between counts columns and metadata rownames",
      call. = FALSE
    )
  }
  if (length(common) < ncol(counts)) {
    message(
      "Subsetting counts to ", length(common),
      " cells present in metadata"
    )
    counts <- counts[, common, drop = FALSE]
  }
  meta_df <- meta_df[common, , drop = FALSE]

  # ------------------------------------------------------------------
  # 4. Process coordinates (if provided)
  # ------------------------------------------------------------------
  coord_matrix <- NULL
  if (!is.null(coords)) {
    if (is.data.frame(coords) || is.matrix(coords)) {
      if (ncol(coords) < 2) {
        stop("coords must have at least 2 columns (x, y)", call. = FALSE)
      }
      coord_matrix <- as.matrix(coords[, 1:2, drop = FALSE])

      # Align by rownames if available, otherwise assume order matches
      if (!is.null(rownames(coord_matrix))) {
        coord_matrix <- coord_matrix[common, , drop = FALSE]
      } else {
        if (nrow(coord_matrix) != length(common)) {
          stop("coords has no rownames and length does not match cells: ",
            "got ", nrow(coord_matrix), " rows, expected ",
            length(common),
            call. = FALSE
          )
        }
        rownames(coord_matrix) <- common
      }

      colnames(coord_matrix) <- c(x_column, y_column)
    } else {
      stop("coords must be a matrix or data.frame", call. = FALSE)
    }
  }

  # ------------------------------------------------------------------
  # 5. Build canonical object
  # ------------------------------------------------------------------
  if (output_class == "SpatialExperiment") {
    missing_pkgs <- c("SpatialExperiment", "S4Vectors")[
      !vapply(c("SpatialExperiment", "S4Vectors"), requireNamespace,
        logical(1), quietly = TRUE)
    ]
    if (length(missing_pkgs)) {
      stop("Package(s) ", paste(missing_pkgs, collapse = ", "),
        " required for output_class = 'SpatialExperiment'.\n",
        "Install with: BiocManager::install(c(",
        paste(sprintf("'%s'", missing_pkgs), collapse = ", "), "))",
        call. = FALSE
      )
    }
    if (is.null(coord_matrix)) {
      stop("coords is required when output_class = 'SpatialExperiment'",
        call. = FALSE
      )
    }

    spe <- SpatialExperiment::SpatialExperiment(
      assays = list(counts = counts),
      colData = S4Vectors::DataFrame(meta_df),
      spatialCoords = coord_matrix
    )
    return(spe)
  }

  if (output_class == "Seurat") {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package 'Seurat' is required for output_class = 'Seurat'",
        call. = FALSE
      )
    }

    # Merge coordinates into metadata for Seurat
    if (!is.null(coord_matrix)) {
      meta_df[[x_column]] <- coord_matrix[, 1]
      meta_df[[y_column]] <- coord_matrix[, 2]
    }

    # Coerce dense matrix to sparse to silence Seurat's coercion warning
    if (is.matrix(counts)) {
      counts <- methods::as(counts, "CsparseMatrix")
    }
    obj <- Seurat::CreateSeuratObject(
      counts    = counts,
      meta.data = meta_df
    )
    return(obj)
  }
}

#' Read a RIPPLE input from a CSV directory
#'
#' Convenience wrapper that loads counts, metadata, and coordinates from CSV
#' files in a directory, then calls \code{\link{make_ripple_input}} to build
#' a canonical object.
#'
#' @param dir_path Path to a directory containing:
#' \describe{
#'   \item{counts.csv}{First column = gene names, remaining columns = cells.}
#'   \item{metadata.csv}{First column = cell barcodes, remaining columns =
#'     metadata fields.}
#'   \item{coords.csv (optional)}{First column = cell barcodes, two additional
#'     columns for x and y coordinates.}
#' }
#' @param output_class Character. Passed to \code{\link{make_ripple_input}}.
#' @param ... Additional arguments forwarded to \code{\link{make_ripple_input}}.
#'
#' @return A \code{SpatialExperiment} or \code{Seurat} object.
#'
#' @examples
#' \dontrun{
#' spe <- read_ripple_csv("/path/to/csv_dir", output_class = "SpatialExperiment")
#' results <- run_ripple(spe,
#'   query_celltype = "Tumor",
#'   celltype_column = "cell_type"
#' )
#' }
#'
#' @importFrom data.table fread setnames
#' @export
read_ripple_csv <- function(dir_path,
                            output_class = c("SpatialExperiment", "Seurat"),
                            ...) {
  output_class <- match.arg(output_class)

  counts_file <- file.path(dir_path, "counts.csv")
  meta_file <- file.path(dir_path, "metadata.csv")
  coords_file <- file.path(dir_path, "coords.csv")

  if (!file.exists(counts_file)) {
    stop("counts.csv not found in: ", dir_path, call. = FALSE)
  }
  if (!file.exists(meta_file)) {
    stop("metadata.csv not found in: ", dir_path, call. = FALSE)
  }

  message("Loading counts from: ", counts_file)
  counts_df <- data.table::fread(counts_file)
  gene_names <- counts_df[[1]]
  counts <- as.matrix(counts_df[, -1, with = FALSE])
  rownames(counts) <- gene_names
  counts <- methods::as(counts, "CsparseMatrix")

  message("Loading metadata from: ", meta_file)
  meta_dt <- data.table::fread(meta_file)
  data.table::setnames(meta_dt, 1, "barcode")

  coords_mat <- NULL
  if (file.exists(coords_file)) {
    message("Loading coordinates from: ", coords_file)
    coords_dt <- data.table::fread(coords_file)
    data.table::setnames(coords_dt, 1, "barcode")
    coord_names <- setdiff(names(coords_dt), "barcode")
    if (length(coord_names) < 2) {
      stop("coords.csv must have at least 2 coordinate columns", call. = FALSE)
    }
    coords_mat <- as.matrix(coords_dt[, coord_names[1:2], with = FALSE])
    rownames(coords_mat) <- coords_dt$barcode
  }

  make_ripple_input(
    counts       = counts,
    metadata     = meta_dt,
    coords       = coords_mat,
    output_class = output_class,
    ...
  )
}
