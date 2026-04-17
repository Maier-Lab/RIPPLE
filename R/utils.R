#' @title General Utility Functions
#'
#' @description Miscellaneous helper functions for statistical calculations,
#'   gene scoring, neighborhood analysis, and enrichment analysis.
#'
#' @name utils
NULL

#' Calculate enrichment scores
#'
#' Computes enrichment statistics by comparing observed counts to expected
#' counts across categories.
#'
#' @param observed Named numeric vector of observed counts.
#' @param expected Named numeric vector of expected counts. Names must
#'   overlap with \code{observed}.
#'
#' @return A \code{data.table} with columns:
#' \describe{
#'   \item{\code{category}}{Character. Category names.}
#'   \item{\code{observed}}{Numeric. Observed counts.}
#'   \item{\code{expected}}{Numeric. Expected counts.}
#'   \item{\code{enrichment}}{Numeric. Observed / expected ratio.}
#'   \item{\code{log2_enrichment}}{Numeric. Log2 of the enrichment ratio.}
#' }
#'
#' @details Only categories present in both \code{observed} and \code{expected}
#'   are included. Categories with zero expected counts get \code{NA} for
#'   enrichment and log2_enrichment.
#'
#' @examples
#' \dontrun{
#' obs <- c(A = 50, B = 10, C = 40)
#' exp <- c(A = 33, B = 33, C = 33)
#' calculate_enrichment(obs, exp)
#' }
#'
#' @importFrom data.table data.table
#' @export
calculate_enrichment <- function(observed, expected) {
  # Ensure same order
  common_names <- intersect(names(observed), names(expected))
  observed <- observed[common_names]
  expected <- expected[common_names]

  # Avoid division by zero
  expected_safe <- ifelse(expected == 0, NA_real_, expected)

  data.table::data.table(
    category = common_names,
    observed = observed,
    expected = expected,
    enrichment = observed / expected_safe,
    log2_enrichment = log2(observed / expected_safe)
  )
}


#' P-value from permutation null distribution
#'
#' Computes an empirical p-value by comparing an observed statistic to a
#' null distribution from permutations.
#'
#' @param observed Numeric. The observed test statistic.
#' @param null_distribution Numeric vector. Null statistics from permutations.
#' @param alternative Character. Type of test: \code{"two.sided"} (default),
#'   \code{"greater"}, or \code{"less"}.
#'
#' @return Numeric. The empirical p-value, corrected for finite permutations
#'   (minimum: \code{1 / (n_perms + 1)}).
#'
#' @examples
#' \dontrun{
#' null <- rnorm(1000)
#' permutation_pvalue(2.5, null) # two-sided
#' permutation_pvalue(2.5, null, "greater") # one-sided
#' }
#'
#' @export
permutation_pvalue <- function(observed, null_distribution, alternative = "two.sided") {
  n_perms <- length(null_distribution)

  p <- switch(alternative,
    "two.sided" = sum(abs(null_distribution) >= abs(observed)) / n_perms,
    "greater" = sum(null_distribution >= observed) / n_perms,
    "less" = sum(null_distribution <= observed) / n_perms,
    stop("Unknown alternative: ", alternative)
  )

  # Correct for finite permutations
  p <- max(p, 1 / (n_perms + 1))

  return(p)
}
