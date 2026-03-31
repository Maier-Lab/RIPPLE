test_that("permutation_pvalue is bounded", {
  null <- rnorm(1000, mean = 0)
  # Extreme observed value
  pval <- permutation_pvalue(5.0, null)
  expect_true(pval >= 0 && pval <= 1)
  expect_true(pval < 0.01)

  # Typical value
  pval2 <- permutation_pvalue(0.0, null)
  expect_true(pval2 > 0.5)
})
