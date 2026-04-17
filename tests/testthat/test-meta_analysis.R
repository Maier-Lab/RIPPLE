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

test_that("run_meta_analysis returns correct structure", {
  coefs <- c(-0.005, -0.003, -0.004, -0.006)
  ses <- c(0.001, 0.002, 0.001, 0.002)
  ids <- paste0("sample_", 1:4)

  result <- run_meta_analysis(coefs, ses, ids)

  expect_type(result, "list")
  expect_true("combined_coef" %in% names(result))
  expect_true("combined_pval" %in% names(result))
  expect_true(result$combined_coef < 0)
})
