#' @title GLM Functions for Distance-Dependent Expression Analysis
#'
#' @description Core statistical functions for the RIPPLE pipeline. Fits Poisson
#'   generalized linear models to test whether gene expression changes as a
#'   function of distance from query cells.
#'
#' @name glm_functions
NULL

#' Fit Poisson GLM for distance-dependent expression
#'
#' Fits a Poisson generalized linear model to test whether gene expression
#' changes as a function of distance from query cells, with a cell-size offset.
#'
#' @param counts Integer vector of raw transcript counts for one gene.
#' @param distances Numeric vector of distances to nearest query cell (um).
#' @param total_counts Numeric vector of total counts per cell (library size).
#' @param min_cells Integer. Minimum non-zero cells required (default: 25).
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{beta}}{Numeric. Log-rate change per um of distance.}
#'   \item{\code{se}}{Numeric. Standard error of the coefficient.}
#'   \item{\code{pval}}{Numeric. Wald z-test p-value.}
#'   \item{\code{dispersion}}{Numeric. Residual deviance / residual df
#'     (values >> 1 suggest overdispersion).}
#'   \item{\code{n_cells}}{Integer. Number of valid cells used in the fit.}
#' }
#' Returns \code{NULL} if too few expressing cells.
#'
#' @details The model is:
#'
#'   \code{log(E[Y]) = alpha + beta * distance + log(total_counts)}
#'
#'   A negative beta indicates expression increases near query cells (induced).
#'   A positive beta indicates expression decreases near query cells (repressed).
#'
#'   The \code{offset(log(total_counts))} converts the model from counts to rates,
#'   controlling for differences in cell size, ambient RNA, and segmentation
#'   quality. This is critical, without it, cells near query cells may appear
#'   to express more of everything due to technical artifacts!
#'
#'   Invalid values (NA, non-finite, zero total_counts) are automatically removed
#'   before fitting.
#'
#' @examples
#' \dontrun{
#' result <- fit_poisson(
#'   counts = rpois(100, lambda = 5),
#'   distances = runif(100, 0, 200),
#'   total_counts = rpois(100, lambda = 5000)
#' )
#' if (!is.null(result)) {
#'   cat("Beta:", result$beta, "P-value:", result$pval, "\n")
#' }
#' }
#'
#' @importFrom stats glm poisson
#' @export
fit_poisson <- function(counts, distances, total_counts, min_cells = 25) {

  # Remove NAs and invalid values
  valid_idx <- !is.na(counts) & !is.na(distances) &
               is.finite(counts) & is.finite(distances) &
               !is.na(total_counts) & total_counts > 0
  counts <- counts[valid_idx]
  distances <- distances[valid_idx]
  log_total <- log(total_counts[valid_idx])

  # Need enough cells
  if (length(counts) < min_cells) return(NULL)

  # Need some non-zero counts for the model to be meaningful
  if (sum(counts > 0) < min_cells) return(NULL)

  # Poisson GLM with offset for cell size with tryCatch for robustness
  fit <- tryCatch({
    suppressWarnings(stats::glm(counts ~ distances + offset(log_total),
                                family = stats::poisson()))
  }, error = function(e) NULL)

  if (is.null(fit) || !fit$converged) return(NULL)

  coef_summary <- summary(fit)$coefficients
  if (!"distances" %in% rownames(coef_summary)) return(NULL)

  # Overdispersion diagnostic: residual deviance / residual df
  # Values >> 1 suggest negative binomial would be more appropriate. shouldn't really happen for Xenium at least.
  dispersion <- fit$deviance / fit$df.residual

  list(
    beta = coef_summary["distances", "Estimate"],       # log-rate change per um
    se = coef_summary["distances", "Std. Error"],
    pval = coef_summary["distances", "Pr(>|z|)"],       # Wald z-test
    dispersion = dispersion,
    n_cells = length(counts)
  )
}


#' Fit bivariate Poisson GLM with confounder control
#'
#' Fits a bivariate Poisson GLM with two distance predictors, one to the
#' query cell type and one to a control cell type, to isolate query-specific
#' effects from general niche effects. ONLY SUGGESTIVE, don't take it at face value for interpretation.
#'
#' @param counts Integer vector of raw transcript counts for one gene.
#' @param dist_query Numeric vector of distances to nearest query cell (um).
#' @param dist_control Numeric vector of distances to nearest control cell (um).
#' @param total_counts Numeric vector of total counts per cell (library size).
#' @param min_cells Integer. Minimum non-zero cells required (default: 25).
#'
#' @return A list with the query-specific coefficient:
#' \describe{
#'   \item{\code{beta}}{Numeric. Partial log-rate change per um for query distance,
#'     controlling for control distance.}
#'   \item{\code{se}}{Numeric. Standard error.}
#'   \item{\code{pval}}{Numeric. Wald z-test p-value.}
#'   \item{\code{dispersion}}{Numeric. Overdispersion estimate.}
#'   \item{\code{n_cells}}{Integer. Number of valid cells used.}
#' }
#' Returns \code{NULL} if too few expressing cells.
#'
#' @details The model is:
#'
#'   \code{log(E[Y]) = alpha + beta_query * dist_query + beta_control * dist_control
#'     + log(total_counts)}
#'
#'   The returned \code{beta} is for the query distance term only (partial effect).
#'   If this is significant, the gene's expression gradient near query cells
#'   persists even after accounting for proximity to the control cell type,
#'   suggesting a query-specific effect rather than a general niche effect.
#'
#' @examples
#' \dontrun{
#' result <- fit_poisson_controlled(
#'   counts = rpois(100, lambda = 5),
#'   dist_query = runif(100, 0, 200),
#'   dist_control = runif(100, 0, 200),
#'   total_counts = rpois(100, lambda = 5000)
#' )
#' }
#'
#' @importFrom stats glm poisson
#' @export
fit_poisson_controlled <- function(counts, dist_query, dist_control, total_counts,
                                   min_cells = 25) {
  # Remove NAs and invalid values
  valid_idx <- !is.na(counts) & !is.na(dist_query) & !is.na(dist_control) &
               is.finite(counts) & is.finite(dist_query) & is.finite(dist_control) &
               !is.na(total_counts) & total_counts > 0
  counts <- counts[valid_idx]
  dist_query <- dist_query[valid_idx]
  dist_control <- dist_control[valid_idx]
  log_total <- log(total_counts[valid_idx])

  if (length(counts) < min_cells) return(NULL)

  # Need some non-zero counts for the model to be meaningful
  if (sum(counts > 0) < min_cells) return(NULL)

  # Fit bivariate Poisson GLM with cell-size offset
  fit <- tryCatch({
    suppressWarnings(stats::glm(counts ~ dist_query + dist_control + offset(log_total),
                                family = stats::poisson()))
  }, error = function(e) NULL)

  if (is.null(fit) || !fit$converged) return(NULL)

  coef_summary <- summary(fit)$coefficients
  if (!"dist_query" %in% rownames(coef_summary)) return(NULL)

  # Overdispersion diagnostic
  dispersion <- fit$deviance / fit$df.residual

  list(
    beta = coef_summary["dist_query", "Estimate"],       # log-rate change per um
    se = coef_summary["dist_query", "Std. Error"],
    pval = coef_summary["dist_query", "Pr(>|z|)"],       # Wald z-test
    dispersion = dispersion,
    n_cells = length(counts)
  )
}


#' Classify expression decay pattern - use at your own risk
#'
#' Fits multiple Poisson regression models to classify the shape of the
#' expression-distance relationship. Compares linear, step (at various
#' thresholds), and exponential decay models using AIC.
#'
#' @param counts Integer vector of raw transcript counts.
#' @param distances Numeric vector of distances to query cells (um).
#' @param total_counts Numeric vector of total counts per cell (library size).
#'
#' @return A character string indicating the best-fitting pattern:
#'   \code{"linear"}, \code{"exponential"}, \code{"step_10um"},
#'   \code{"step_25um"}, \code{"step_50um"}, \code{"none"},
#'   \code{"insufficient_data"}, or \code{"no_variation"}.
#'
#' @details Fits the following models (all with Poisson family and cell-size offset):
#' \itemize{
#'   \item Linear: \code{counts ~ distance + offset(log(total_counts))}
#'   \item Step at 10/25/50 um: \code{counts ~ I(distance < threshold) + offset(log(total_counts))}
#'   \item Exponential: fitted via NLS on binned mean rates
#' }
#'
#' The exponential model gets a +10 AIC penalty since it uses a different fitting
#' method (NLS on binned data vs. GLM on individual cells).
#'
#' @examples
#' \dontrun{
#' pattern <- classify_decay_pattern(
#'   counts = rpois(500, lambda = 3),
#'   distances = runif(500, 0, 200),
#'   total_counts = rpois(500, lambda = 5000)
#' )
#' }
#'
#' @importFrom stats glm poisson AIC nls
#' @export
classify_decay_pattern <- function(counts, distances, total_counts) {

  valid_idx <- !is.na(counts) & !is.na(distances) &
               is.finite(counts) & is.finite(distances) &
               !is.na(total_counts) & total_counts > 0
  counts <- counts[valid_idx]
  distances <- distances[valid_idx]
  log_total <- log(total_counts[valid_idx])

  if (length(counts) < 30) return("insufficient_data")
  if (sum(counts > 0) < 5) return("no_variation")

  # Fit linear Poisson model with offset
  fit_linear <- tryCatch({
    fit <- suppressWarnings(stats::glm(counts ~ distances + offset(log_total),
                                       family = stats::poisson()))
    if (fit$converged) fit else NULL
  }, error = function(e) NULL)

  # Fit step models at different thresholds with offset
  fit_step_10 <- tryCatch({
    fit <- suppressWarnings(stats::glm(counts ~ I(distances < 10) + offset(log_total),
                                       family = stats::poisson()))
    if (fit$converged) fit else NULL
  }, error = function(e) NULL)

  fit_step_25 <- tryCatch({
    fit <- suppressWarnings(stats::glm(counts ~ I(distances < 25) + offset(log_total),
                                       family = stats::poisson()))
    if (fit$converged) fit else NULL
  }, error = function(e) NULL)

  fit_step_50 <- tryCatch({
    fit <- suppressWarnings(stats::glm(counts ~ I(distances < 50) + offset(log_total),
                                       family = stats::poisson()))
    if (fit$converged) fit else NULL
  }, error = function(e) NULL)

  # Fit exponential decay (approximate via binned mean rates)
  fit_exp <- tryCatch({
    bins <- cut(distances, breaks = seq(0, max(distances) + 10, by = 20),
                include.lowest = TRUE)
    bin_rates <- tapply(counts / exp(log_total), bins, mean)
    bin_mids <- tapply(distances, bins, mean)

    valid_bins <- !is.na(bin_rates) & !is.na(bin_mids) & bin_rates > 0
    if (sum(valid_bins) < 3) return(NULL)

    bin_rates <- bin_rates[valid_bins]
    bin_mids <- bin_mids[valid_bins]
    bin_n <- tapply(counts, bins, length)[valid_bins]

    log_rates <- log(bin_rates)

    stats::nls(log_rates ~ a * exp(-b * bin_mids) + c,
               start = list(a = max(log_rates) - min(log_rates), b = 0.02,
                            c = min(log_rates)),
               weights = bin_n,
               control = list(maxiter = 100, warnOnly = TRUE))
  }, error = function(e) NULL, warning = function(w) NULL)

  aics <- c(
    linear = if (!is.null(fit_linear)) stats::AIC(fit_linear) else Inf,
    exponential = if (!is.null(fit_exp)) stats::AIC(fit_exp) + 10 else Inf,
    step_10um = if (!is.null(fit_step_10)) stats::AIC(fit_step_10) else Inf,
    step_25um = if (!is.null(fit_step_25)) stats::AIC(fit_step_25) else Inf,
    step_50um = if (!is.null(fit_step_50)) stats::AIC(fit_step_50) else Inf
  )

  if (all(is.infinite(aics))) return("none")
  names(aics)[which.min(aics)]
}
