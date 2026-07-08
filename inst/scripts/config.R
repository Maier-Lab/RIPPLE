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
# Data Paths (env var resolution)
# =============================================================================

INPUT_PATH  <- Sys.getenv("INPUT_PATH",  unset = "")
OUTPUT_DIR  <- Sys.getenv("OUTPUT_DIR",  unset = "")

# Column names
SAMPLE_COL    <- Sys.getenv("SAMPLE_COLUMN",    unset = "sample_id")
CONDITION_COL <- Sys.getenv("CONDITION_COLUMN",  unset = "")
CONDITION_VAL <- Sys.getenv("CONDITION_VALUE",   unset = "")

# Coordinate columns (empty = auto-detect)
X_COL_ENV <- Sys.getenv("X_COLUMN", unset = "")
Y_COL_ENV <- Sys.getenv("Y_COLUMN", unset = "")

# =============================================================================
# Platform Detection (legacy — only needed when INPUT_PATH is not set)
# =============================================================================

if (nchar(INPUT_PATH) == 0) {
  # No INPUT_PATH set: fall back to a base path from the RIPPLE_BASE_PATH
  # environment variable (point this at your project root).
  BASE_PATH <- Sys.getenv("RIPPLE_BASE_PATH", unset = "/path/to/project")
  PROJECT_ROOT <- file.path(BASE_PATH, "analysis")
  MXENIUM_ROOT <- BASE_PATH  # Parent project root
} else {
  # External mode: derive PROJECT_ROOT from OUTPUT_DIR or working directory
  BASE_PATH    <- dirname(INPUT_PATH)
  MXENIUM_ROOT <- BASE_PATH
  PROJECT_ROOT <- if (nchar(OUTPUT_DIR) > 0) dirname(OUTPUT_DIR) else getwd()
}

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

if (nchar(OUTPUT_DIR) > 0) {
  OUTPUT_ROOT <- file.path(OUTPUT_DIR, paste0("spatial_analysis", OUTPUT_SUFFIX))
} else {
  OUTPUT_ROOT <- file.path(PROJECT_ROOT, "results", paste0("spatial_analysis", OUTPUT_SUFFIX))
}

# =============================================================================
# Configuration Summary
# =============================================================================

message("--- RIPPLE Configuration ---")
if (nchar(INPUT_PATH) > 0) {
  message(sprintf("  INPUT_PATH:       %s", INPUT_PATH))
} else {
  message(sprintf("  INPUT_PATH:       (not set; using RIPPLE_BASE_PATH fallback)"))
}
message(sprintf("  OUTPUT_ROOT:      %s", OUTPUT_ROOT))
message(sprintf("  SAMPLE_COL:       %s", SAMPLE_COL))
if (nchar(CONDITION_COL) > 0) {
  message(sprintf("  CONDITION_COL:    %s = %s",
                  CONDITION_COL,
                  if (nchar(CONDITION_VAL) > 0) CONDITION_VAL else "(all)"))
}
if (nchar(X_COL_ENV) > 0 && nchar(Y_COL_ENV) > 0) {
  message(sprintf("  Coordinates:      %s, %s", X_COL_ENV, Y_COL_ENV))
} else {
  message("  Coordinates:      (auto-detect)")
}
message(sprintf("  CELLTYPE_COL:     %s", CELLTYPE_COL))
message(sprintf("  QUERY_CELLTYPE:   %s", QUERY_CELLTYPE))
message("----------------------------")
