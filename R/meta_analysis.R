#' @title Meta-Analysis Functions
#'
#' @description Fisher's combined p-value with a sign-consistency gate for
#'   combining per-sample results across biological replicates. This is the
#'   only cross-replicate combination step in RIPPLE; per-gene effect sizes
#'   are summarised as the median of per-sample log-rate coefficients.
#'
#' @name meta_analysis
NULL


#' Fisher's combined p-value with sign consistency gate - RECOMMENDED APPROACH
#'
#' Combines per-sample p-values using Fisher's method, with an additional
#' requirement that a minimum fraction of samples agree on the direction
#' of the effect. This prevents combining p-values when samples disagree
#' on whether a gene is induced or repressed.
#'
#' @param pvals Numeric vector of per-sample p-values.
#' @param coefs Numeric vector of per-sample coefficients (for sign checking).
#' @param min_samples Integer. Minimum number of valid samples required
#'   (default: 2).
#' @param sign_threshold Numeric. Required fraction of samples agreeing on
#'   effect direction (default: 1.0, i.e., all samples must agree).
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{fisher_pval}}{Numeric. Fisher's combined p-value
#'     (set to 1.0 if sign consistency is below threshold).}
#'   \item{\code{fisher_stat}}{Numeric. Fisher's test statistic
#'     (-2 * sum(log(p_i))).}
#'   \item{\code{median_coef}}{Numeric. Median coefficient across samples
#'     (equal weighting per replicate).}
#'   \item{\code{sign_consistency}}{Numeric. Fraction of samples agreeing
#'     on direction.}
#'   \item{\code{n_valid}}{Integer. Number of valid samples used.}
#' }
#'
#' @details Fisher's method combines independent p-values:
#'
#'   \code{X^2 = -2 * sum(log(p_i)), df = 2k}
#'
#'   where k is the number of p-values. The sign consistency gate ensures
#'   that the combined result is only significant when replicates agree on
#'   the direction of the effect. Genes with contradictory directions across
#'   samples get \code{fisher_pval = 1} regardless of their individual
#'   p-values.
#'
#'   P-values are clamped to a minimum of \code{1e-15} to avoid
#'   \code{log(0)}.
#'
#' @examples
#' \dontrun{
#' result <- compute_fisher_pval(
#'   pvals = c(0.01, 0.001, 0.05),
#'   coefs = c(-0.005, -0.003, -0.007)
#' )
#' cat("Fisher p-value:", result$fisher_pval, "\n")
#' cat("Median coefficient:", result$median_coef, "\n")
#' }
#'
#' @importFrom stats pchisq
#' @export
compute_fisher_pval <- function(pvals, coefs, min_samples = 2,
                                sign_threshold = 1.0) {
  # Filter to valid entries
  valid <- !is.na(pvals) & !is.na(coefs) & pvals > 0
  valid_pvals <- pvals[valid]
  valid_coefs <- coefs[valid]
  n_valid <- length(valid_pvals)

  na_result <- list(
    fisher_pval = NA_real_,
    fisher_stat = NA_real_,
    median_coef = NA_real_,
    sign_consistency = NA_real_,
    n_valid = n_valid
  )

  if (n_valid < min_samples) {
    return(na_result)
  }

  # Median coefficient (equal replicate weighting)
  median_coef <- median(valid_coefs)

  # Sign consistency check
  # Exact zeros are treated as agreeing with both directions so they don't
  # artificially drag the consistency below threshold. If all coefficients
  # are zero, consistency is defined as 1 (nothing disagrees).
  n_pos <- sum(valid_coefs > 0)
  n_neg <- sum(valid_coefs < 0)
  n_nonzero <- n_pos + n_neg
  if (n_nonzero == 0) {
    sign_consistency <- 1.0
  } else {
    sign_consistency <- max(n_pos, n_neg) / n_nonzero
  }

  if (sign_consistency < sign_threshold) {
    # Contradictory directions -- not significant
    return(list(
      fisher_pval = 1.0,
      fisher_stat = NA_real_,
      median_coef = median_coef,
      sign_consistency = sign_consistency,
      n_valid = n_valid
    ))
  }

  # Fisher's method: X^2 = -2 * sum(log(p_i)), df = 2k
  clamped_pvals <- pmax(valid_pvals, 1e-15) # avoid log(0)
  fisher_stat <- -2 * sum(log(clamped_pvals))
  fisher_pval <- stats::pchisq(fisher_stat, df = 2 * n_valid, lower.tail = FALSE)

  list(
    fisher_pval = fisher_pval,
    fisher_stat = fisher_stat,
    median_coef = median_coef,
    sign_consistency = sign_consistency,
    n_valid = n_valid
  )
}
