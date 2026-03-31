#' @title Meta-Analysis Functions
#'
#' @description Functions for combining per-sample results across biological
#'   replicates. Includes random-effects meta-analysis and Fisher's combined
#'   p-value with sign consistency gating.
#'
#' @name meta_analysis
NULL

#' Random-effects meta-analysis across samples
#'
#' Combines per-sample log-rate coefficients using a random-effects
#' meta-analysis model via \code{meta::metagen}. This accounts for both
#' within-sample and between-sample variability. This is not the recommended solution for most cases.
#'
#' @param coefs Numeric vector of per-sample coefficient estimates.
#' @param ses Numeric vector of per-sample standard errors.
#' @param sample_ids Character vector of sample identifiers.
#' @param method Character. Method for estimating between-study variance
#'   (default: \code{"REML"}). Passed to \code{meta::metagen}.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{combined_coef}}{Numeric. Random-effects combined coefficient.}
#'   \item{\code{combined_se}}{Numeric. Standard error of the combined estimate.}
#'   \item{\code{pval}}{Numeric. P-value for the combined estimate.}
#'   \item{\code{i2}}{Numeric. I-squared heterogeneity statistic (0-1).}
#'   \item{\code{n_samples}}{Integer. Number of valid samples used.}
#' }
#' All values are \code{NA_real_} if fewer than 2 valid samples.
#'
#' @details Uses the REML estimator for between-study variance (tau-squared)
#'   by default, which is for small numbers of studies. Samples
#'   with NA coefficients, NA standard errors, or zero standard errors are
#'   excluded.
#'
#' @examples
#' \dontrun{
#' result <- run_meta_analysis(
#'   coefs = c(-0.005, -0.003, -0.007),
#'   ses = c(0.001, 0.002, 0.001),
#'   sample_ids = c("mouse1", "mouse2", "mouse3")
#' )
#' cat("Combined coefficient:", result$combined_coef, "\n")
#' cat("P-value:", result$pval, "\n")
#' }
#'
#' @export
run_meta_analysis <- function(coefs, ses, sample_ids, method = "REML") {
  # Remove NAs
  valid_idx <- !is.na(coefs) & !is.na(ses) & ses > 0

  na_result <- list(
    combined_coef = NA_real_,
    combined_se = NA_real_,
    pval = NA_real_,
    i2 = NA_real_,
    n_samples = sum(valid_idx)
  )

  if (sum(valid_idx) < 2) {
    return(na_result)
  }

  coefs <- coefs[valid_idx]
  ses <- ses[valid_idx]
  sample_ids <- sample_ids[valid_idx]

  meta_result <- tryCatch({
    meta::metagen(
      TE = coefs,
      seTE = ses,
      studlab = sample_ids,
      random = TRUE,
      method.tau = method
    )
  }, error = function(e) NULL)

  if (is.null(meta_result)) {
    return(na_result)
  }

  list(
    combined_coef = meta_result$TE.random,
    combined_se = meta_result$seTE.random,
    pval = meta_result$pval.random,
    i2 = meta_result$I2,
    n_samples = length(coefs)
  )
}


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
  n_pos <- sum(valid_coefs > 0)
  n_neg <- sum(valid_coefs < 0)
  sign_consistency <- max(n_pos, n_neg) / n_valid

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
  clamped_pvals <- pmax(valid_pvals, 1e-15)  # avoid log(0)
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
