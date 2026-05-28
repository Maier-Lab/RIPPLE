test_that("compute_fisher_pval combines consistent results", {
  # 4 samples all showing negative effect
  pvals <- c(0.01, 0.02, 0.05, 0.03)
  coefs <- c(-0.005, -0.003, -0.004, -0.006)

  result <- compute_fisher_pval(pvals, coefs)

  expect_type(result, "list")
  expect_true(result$fisher_pval < 0.01)
  expect_equal(result$sign_consistency, 1.0)
  expect_equal(result$median_coef, median(coefs))
})

test_that("compute_fisher_pval gates inconsistent signs", {
  # 4 samples with mixed signs
  pvals <- c(0.01, 0.02, 0.05, 0.03)
  coefs <- c(-0.005, 0.003, -0.004, 0.006) # mixed signs

  result <- compute_fisher_pval(pvals, coefs)

  expect_equal(result$fisher_pval, 1.0) # Gated
  expect_true(result$sign_consistency < 1.0)
})

test_that("compute_fisher_pval handles single sample", {
  result <- compute_fisher_pval(c(0.01), c(-0.005), min_samples = 1)
  expect_type(result, "list")
})

# run_meta_analysis() was removed when the REML inverse-variance step was
# dropped from the hot path. gradient_score is now the median of per-sample
# coefficients (see test below).

test_that("median_coef equals median of per-sample coefficients", {
  coefs <- c(-0.005, -0.003, -0.004, -0.006)
  pvals <- c(0.01, 0.02, 0.05, 0.03)

  result <- compute_fisher_pval(pvals, coefs)
  expect_equal(result$median_coef, median(coefs))
})
