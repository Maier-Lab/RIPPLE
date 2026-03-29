#' =============================================================================
#' Load data for RIPPLE spatial analysis
#'
#' Loads Seurat objects and extracts metadata for the pipeline.
#' Sources utils.R (which sources config.R) for shared configuration.
#'
#' Usage:
#'   source("scripts/load_data.R")
#'   obj <- load_seurat()
#'   data <- load_metadata_only()  # Fast: just coords + cell types
#' =============================================================================

# Load utilities - get script directory (works with Rscript and source())
get_script_dir <- function() {
  # Try commandArgs for Rscript
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  # Try sys.frame for source()
  for (i in seq_len(sys.nframe())) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  }
  # Fallback: current working directory
  getwd()
}
script_dir <- get_script_dir()
source(file.path(script_dir, "utils.R"))

# =============================================================================
# Loading Functions
# =============================================================================

#' Load just metadata and coordinates (fast, for spatial-only analysis)
#'
#' Extracts cell_id, cell types, x, y, and key metadata from the Seurat object.
#' Uses CELLTYPE_COLUMN from utils.R to determine which annotation to use.
#' For HyMy annotation: merges from CSV. For L1: uses existing column.
#'
#' @param path Path to Seurat RDS file (optional)
#' @return data.table with cell metadata and coordinates
load_metadata_only <- function(path = NULL) {
  message("Loading metadata from Seurat object...")
  message(sprintf("  Query cell type: %s (column: %s)", QUERY_CELLTYPE, CELLTYPE_COLUMN))

  obj <- load_seurat(path)

  # Extract metadata - use "barcode" for rownames (cell_id column may already exist)
  meta <- as.data.table(obj@meta.data, keep.rownames = "barcode")

  # Create cell_id from barcode if not present (for backward compatibility)
  if (!"cell_id" %in% names(meta)) {
    meta[, cell_id := barcode]
  }

  # Resolve coordinate columns and create x/y aliases

  coord_cols <- get_coord_columns(meta)
  if (coord_cols[1] != "x") meta[, x := get(coord_cols[1])]
  if (coord_cols[2] != "y") meta[, y := get(coord_cols[2])]

  # Handle cell type annotation based on configuration
  if (USE_HYMY_ANNOTATION && nchar(INPUT_PATH) == 0) {
    # HyMy annotation: Merge from CSV (legacy CeMM mode only)
    if (!"cell_type_with_HyMy" %in% names(meta) ||
        !any(meta$cell_type_with_HyMy == "HyMy_GMM", na.rm = TRUE)) {
      message("\nMerging HyMy annotations from CSV...")
      hymy <- load_hymy_annotations()
      # Rename cell_id to barcode for merge
      setnames(hymy, "cell_id", "barcode")
      meta <- merge(meta, hymy[, .(barcode, cell_type_with_HyMy)],
                    by = "barcode", all.x = TRUE)
    }
  } else if (ANNOTATION_LEVEL == "L1") {
    # L1 annotation: Use cell_type_assignment_L1 directly
    # Create cell_type_with_HyMy as alias for compatibility with existing scripts
    if (CELLTYPE_COLUMN %in% names(meta)) {
      meta[, cell_type_with_HyMy := get(CELLTYPE_COLUMN)]
      message(sprintf("\nUsing %s column for cell types", CELLTYPE_COLUMN))
    } else {
      stop(sprintf("Column '%s' not found in Seurat metadata", CELLTYPE_COLUMN))
    }
  } else {
    # User-specified cell type column: verify it exists, no alias needed
    if (CELLTYPE_COLUMN %in% names(meta)) {
      message(sprintf("\nUsing %s column for cell types", CELLTYPE_COLUMN))
    } else {
      stop(sprintf("Column '%s' not found in Seurat metadata.\n  Available columns: %s",
                    CELLTYPE_COLUMN,
                    paste(head(names(meta), 20), collapse = ", ")))
    }
  }

  message(sprintf("\nLoaded metadata for %s cells", nrow(meta)))

  # Print summary using the appropriate column
  message(sprintf("\nCell type summary (using %s):", CELLTYPE_COLUMN))
  ct_summary <- meta[, .N, by = c(CELLTYPE_COLUMN)][order(-N)]
  print(ct_summary[1:min(15, nrow(ct_summary))])

  # Print query cell type count
  query_count <- sum(meta[[CELLTYPE_COLUMN]] == QUERY_CELLTYPE, na.rm = TRUE)
  message(sprintf("\n  Query cell type (%s): %d cells", QUERY_CELLTYPE, query_count))

  # Clean up Seurat object to free memory
  rm(obj)
  gc(verbose = FALSE)

  return(meta)
}


#' Load full Seurat object (for expression-based analyses)
#'
#' This is just an alias for load_seurat() from utils.R
#' Use this when you need gene expression data.
#'
#' @param path Path to Seurat RDS file (optional)
#' @return Seurat object
load_full_data <- function(path = NULL) {
  load_seurat(path)
}


#' Quick check that data exists and print summary
check_data <- function(path = NULL) {
  path <- if (is.null(path)) SEURAT_PATH else path

  message("Checking Seurat object at: ", path)
  message("")

  if (!file.exists(path)) {
    message("  ✗ File NOT FOUND")
    message("")
    message("To create the Seurat object, run:")
    message("  Rscript scripts/annotate_HyMy_with_UCell.R")
    return(invisible(FALSE))
  }

  size_mb <- file.info(path)$size / (1024 * 1024)
  message(sprintf("  ✓ Found (%.1f MB)", size_mb))

  # Quick load to check contents
  message("  Loading to check contents...")
  obj <- readRDS(path)

  message(sprintf("  Cells: %s", format(ncol(obj), big.mark = ",")))
  message(sprintf("  Genes: %s", format(nrow(obj), big.mark = ",")))
  message("")

  message("  Metadata columns:")
  # Build dynamic key columns list based on configuration
  key_cols <- unique(c(CELLTYPE_COLUMN, SAMPLE_COL))
  # Add coordinate columns (user-specified or common candidates)
  if (nchar(X_COL_ENV) > 0 && nchar(Y_COL_ENV) > 0) {
    key_cols <- c(key_cols, X_COL_ENV, Y_COL_ENV)
  } else {
    # Check common coordinate column names
    for (pair in list(c("spatial_x", "spatial_y"), c("x", "y"),
                      c("x_centroid", "y_centroid"))) {
      if (all(pair %in% colnames(obj@meta.data))) {
        key_cols <- c(key_cols, pair)
        break
      }
    }
  }
  # Add condition column if set
  if (nchar(CONDITION_COL) > 0) {
    key_cols <- c(key_cols, CONDITION_COL)
  }
  for (col in key_cols) {
    if (col %in% colnames(obj@meta.data)) {
      message(sprintf("    ✓ %s", col))
    } else {
      message(sprintf("    ✗ %s (missing)", col))
    }
  }

  rm(obj)
  gc(verbose = FALSE)

  return(invisible(TRUE))
}


# =============================================================================
# Run check if sourced directly
# =============================================================================

if (sys.nframe() == 0) {
  message("Checking data availability...\n")
  check_data()
}
