#' =============================================================================
#' RIPPLE Package Configuration
#' =============================================================================
#'
#' Package options system using R's options() with ripple. prefix.
#' Environment variables are checked as fallback for backward compatibility
#' with SLURM scripts.
#'
#' @name ripple-config
#' @keywords internal
NULL


# ============================================================================
# .onLoad — initialize package options from env vars
# ============================================================================

.onLoad <- function(libname, pkgname) {
  op <- options()
  ripple_defaults <- list(
    ripple.sample_column       = Sys.getenv("SAMPLE_COLUMN", "sample_id"),
    ripple.condition_column    = Sys.getenv("CONDITION_COLUMN", ""),
    ripple.condition_value     = Sys.getenv("CONDITION_VALUE", ""),
    ripple.k_neighbors         = as.integer(Sys.getenv("K_NEIGHBORS", "1")),
    ripple.max_distance_um     = as.numeric(Sys.getenv("MAX_DISTANCE_UM", "200")),
    ripple.fdr_threshold       = as.numeric(Sys.getenv("FDR_THRESHOLD", "0.05")),
    ripple.min_cells_per_sample = as.integer(Sys.getenv("MIN_CELLS_PER_SAMPLE", "30")),
    ripple.min_expr_cells      = as.integer(Sys.getenv("MIN_EXPR_CELLS", "25")),
    ripple.min_expr_pct        = as.numeric(Sys.getenv("MIN_EXPR_PCT", "0.01")),
    ripple.sign_consistency    = as.numeric(Sys.getenv("SIGN_CONSISTENCY_THRESHOLD", "1.0")),
    ripple.verbose             = TRUE
  )
  # Only set options that aren't already set
  toset <- !(names(ripple_defaults) %in% names(op))
  if (any(toset)) options(ripple_defaults[toset])
  invisible()
}


#' Get or set RIPPLE configuration options
#'
#' RIPPLE uses R's \code{options()} system for configuration. All options are
#' prefixed with \code{ripple.}. Environment variables are checked as fallback
#' when an option is not explicitly set (handled at package load time).
#'
#' @param ... Named arguments to set options, or a single unnamed character
#'   string to retrieve an option value. When setting, provide name-value pairs
#'   (e.g., \code{ripple_config(k_neighbors = 5)}). When getting, provide the
#'   option name without the \code{ripple.} prefix (e.g.,
#'   \code{ripple_config("k_neighbors")}).
#'
#' @return When getting (single character argument), returns the option value.
#'   When setting (named arguments), returns previous values invisibly as a
#'   named list.
#'
#' @details
#' Available options (shown without the \code{ripple.} prefix):
#' \describe{
#'   \item{sample_column}{Metadata column for sample/replicate IDs
#'     (default: \code{"sample_id"})}
#'   \item{condition_column}{Metadata column for condition filtering
#'     (default: \code{""} = no filter)}
#'   \item{condition_value}{Which condition to analyze
#'     (default: \code{""} = all)}
#'   \item{k_neighbors}{Number of nearest query cells for distance
#'     (default: \code{1})}
#'   \item{max_distance_um}{Maximum distance in micrometers
#'     (default: \code{200})}
#'   \item{fdr_threshold}{FDR significance cutoff
#'     (default: \code{0.05})}
#'   \item{min_cells_per_sample}{Minimum cells of target type per sample
#'     (default: \code{30})}
#'   \item{min_expr_cells}{Absolute floor for expressing cells per sample
#'     (default: \code{25})}
#'   \item{min_expr_pct}{Minimum fraction of cells expressing per sample
#'     (default: \code{0.01})}
#'   \item{sign_consistency}{Required fraction of samples agreeing on beta
#'     direction (default: \code{1.0})}
#'   \item{verbose}{Print progress messages (default: \code{TRUE})}
#' }
#'
#' @examples
#' \dontrun{
#' # Get a single option
#' ripple_config("k_neighbors")
#'
#' # Set one or more options
#' ripple_config(k_neighbors = 5, fdr_threshold = 0.1)
#'
#' # List all current RIPPLE options
#' ripple_config()
#' }
#'
#' @export
ripple_config <- function(...) {
  args <- list(...)

  # No arguments: return all ripple options
  if (length(args) == 0) {
    all_opts <- options()
    ripple_opts <- all_opts[grepl("^ripple\\.", names(all_opts))]
    # Strip the ripple. prefix for display
    names(ripple_opts) <- sub("^ripple\\.", "", names(ripple_opts))
    return(ripple_opts)
  }

  # Single unnamed character: get a specific option
  if (is.null(names(args)) && length(args) == 1 && is.character(args[[1]])) {
    opt_name <- paste0("ripple.", args[[1]])
    return(getOption(opt_name))
  }

  # Named arguments: set options
  if (!is.null(names(args)) && all(nzchar(names(args)))) {
    old_vals <- list()
    new_opts <- list()
    for (nm in names(args)) {
      full_name <- paste0("ripple.", nm)
      old_vals[[nm]] <- getOption(full_name)
      new_opts[[full_name]] <- args[[nm]]
    }
    do.call(options, new_opts)
    return(invisible(old_vals))
  }

  stop("Usage: ripple_config('option_name') to get, ",
       "ripple_config(option_name = value) to set.",
       call. = FALSE)
}
