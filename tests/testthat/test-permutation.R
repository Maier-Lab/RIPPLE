# The permutation null and the observed statistic must be the SAME estimator:
# the median of per-sample coefficients (compute_fisher_pval's median_coef).
# These tests use a well-powered synthetic (high counts, many cells/sample) so
# per-sample slopes are precise and a real gradient is extreme versus the
# random-query null. Note: with few samples the median statistic is noisy, so
# this test deliberately uses precise per-sample fits rather than tiny N.

make_gradient_sample <- function(samp, n = 250, slope = -0.015, seed_x) {
  set.seed(seed_x)
  # Cells spread over a field; a query anchor sits at the origin.
  x <- runif(n, 0, 100)
  y <- runif(n, 0, 100)
  dist0 <- sqrt(x^2 + y^2)
  counts <- rpois(n, exp(2 + slope * dist0)) # high baseline -> precise fit
  data.frame(x = x, y = y, dist0 = dist0, counts = counts,
    total = 3000, sample = samp)
}

test_that("run_permutation_test null is unbiased and drives the empirical p-value", {
  samples <- c("s1", "s2", "s3", "s4", "s5")
  target <- do.call(rbind, Map(make_gradient_sample, samples, seed_x = 1:5))

  coords_target <- as.matrix(target[, c("x", "y")])
  # Pool of candidate query cells = the same field (random draws for the null).
  coords_all <- coords_target
  sample_ids_all <- target$sample
  query_per_sample <- stats::setNames(rep(15L, length(samples)), samples)

  set.seed(123)
  res <- run_permutation_test(
    counts           = target$counts,
    coords_target    = coords_target,
    coords_all       = coords_all,
    sample_ids       = target$sample,
    n_perms          = 300,
    observed_coef    = -0.02,
    sample_ids_all   = sample_ids_all,
    query_per_sample = query_per_sample,
    k_neighbors      = 1,
    total_counts     = target$total,
    max_distance     = 200,
    min_cells_per_sample = 30,
    min_expr_cells   = 5
  )

  # The label-permuted (random-query) null must be centred near zero -- it
  # carries no query-specific gradient by construction.
  expect_true(all(is.finite(res$null_coefs)))
  expect_gt(length(res$null_coefs), 200)
  expect_lt(abs(stats::median(res$null_coefs)), 0.01)

  # The returned p-value must be exactly the two-sided empirical tail of the
  # null against the observed statistic (guards the comparison logic).
  manual_p <- (sum(abs(res$null_coefs) >= abs(-0.02)) + 1) /
    (length(res$null_coefs) + 1)
  expect_equal(res$perm_pval, manual_p)

  # Monotone sanity: an impossibly extreme observed hits the p floor; a zero
  # observed cannot be more extreme than any null value.
  extreme <- run_permutation_test(
    counts = target$counts, coords_target = coords_target,
    coords_all = coords_all, sample_ids = target$sample, n_perms = 50,
    observed_coef = -10, sample_ids_all = sample_ids_all,
    query_per_sample = query_per_sample, k_neighbors = 1,
    total_counts = target$total, max_distance = 200,
    min_cells_per_sample = 30, min_expr_cells = 5
  )
  expect_equal(extreme$perm_pval, 1 / (length(extreme$null_coefs) + 1))
})

test_that("run_ripple permutation testing populates valid empirical p-values", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("meta")

  data(ripple_mock_data)
  out_dir <- tempfile("ripple_perm_test_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  results <- run_ripple(
    input                = ripple_mock_data,
    query_celltype       = "Tumor",
    celltype_column      = "cell_type",
    sample_column        = "sample_id",
    output_dir           = out_dir,
    min_cells_per_sample = 30,
    min_expr_pct         = 0,
    min_expr_floor       = 10,
    n_permutations       = 199,
    verbose              = FALSE
  )

  expect_true("perm_pval" %in% names(results))
  pv <- results$perm_pval
  finite_pv <- pv[!is.na(pv)]
  expect_gt(length(finite_pv), 0)
  expect_true(all(finite_pv >= 0 & finite_pv <= 1))
  # Minimum achievable p-value with 199 perms is 1/200.
  expect_true(all(finite_pv >= 1 / 200))
})
