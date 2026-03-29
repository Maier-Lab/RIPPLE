#' =============================================================================
#' RIPPLE Pipeline Configuration
#' =============================================================================
#'
#' Lightweight config module: env var resolution + platform detection.
#' No package loads, no function definitions.
#'
#' Sourced by:
#'   - utils.R (and therefore all heavy scripts via load_data.R)
#'   - Lightweight scripts that only need config (merge_*, recompute_*)
#'
#' Query cell type resolution (3-tier):
#'   Priority 1: QUERY_CELLTYPE + CELLTYPE_COLUMN env vars (both must be set)
#'   Priority 2: ANNOTATION_LEVEL="L1" (legacy, IL1B_myeloid)
#'   Priority 3: Default (legacy, HyMy_GMM)
#'
#' Usage:
#'   source("config.R")
#'   # or from another directory:
#'   source(file.path(script_dir, "config.R"))
#' =============================================================================

# =============================================================================
# Platform Detection
# =============================================================================

if (.Platform$OS.type == "windows") {
  BASE_PATH <- "N:/lab_maier/Projects/mXenium"
} else {
  BASE_PATH <- "/nobackup/lab_maier/Projects/mXenium"
}

PROJECT_ROOT <- file.path(BASE_PATH, "CMM")
MXENIUM_ROOT <- BASE_PATH  # Parent project root

# =============================================================================
# Query Cell Type Configuration (3-tier resolution)
# =============================================================================

QUERY_CELLTYPE <- Sys.getenv("QUERY_CELLTYPE", unset = "")
CELLTYPE_COLUMN <- Sys.getenv("CELLTYPE_COLUMN", unset = "")
ANNOTATION_LEVEL <- Sys.getenv("ANNOTATION_LEVEL", unset = "")

if (nchar(QUERY_CELLTYPE) > 0 && nchar(CELLTYPE_COLUMN) > 0) {
  OUTPUT_SUFFIX <- paste0("_", QUERY_CELLTYPE)
  USE_HYMY_ANNOTATION <- as.logical(Sys.getenv("MERGE_ANNOTATION_CSV", unset = "FALSE"))
  message(sprintf(">>> Query cell type: %s (column: %s)", QUERY_CELLTYPE, CELLTYPE_COLUMN))

} else if (nchar(ANNOTATION_LEVEL) > 0 && ANNOTATION_LEVEL == "L1") {
  CELLTYPE_COLUMN <- "cell_type_assignment_L1"
  QUERY_CELLTYPE <- "IL1B_myeloid"
  OUTPUT_SUFFIX <- "_L1"
  USE_HYMY_ANNOTATION <- FALSE
  message(">>> Using L1 annotation (IL1B_myeloid as query type)")

} else {
  # Default: HyMy (backward compatible)
  CELLTYPE_COLUMN <- "cell_type_with_HyMy"
  QUERY_CELLTYPE <- "HyMy_GMM"
  OUTPUT_SUFFIX <- ""
  USE_HYMY_ANNOTATION <- TRUE
  ANNOTATION_LEVEL <- "HyMy"
  message(">>> Using HyMy annotation (HyMy_GMM as query type)")
}

# Alias (some scripts use CELLTYPE_COL)
CELLTYPE_COL <- CELLTYPE_COLUMN

# Display label for plots (defaults to QUERY_CELLTYPE)
QUERY_LABEL <- Sys.getenv("QUERY_LABEL", unset = QUERY_CELLTYPE)

# =============================================================================
# Analysis Name
# =============================================================================

ANALYSIS_NAME <- Sys.getenv("ANALYSIS_NAME", unset = "hymy_distance_correlation_v2")

# =============================================================================
# Output Paths
# =============================================================================

OUTPUT_ROOT <- file.path(PROJECT_ROOT, "results", paste0("spatial_analysis", OUTPUT_SUFFIX))
