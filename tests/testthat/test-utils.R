test_that("shannon_entropy works correctly", {
  # Equal distribution: maximum entropy
  counts_equal <- c(25, 25, 25, 25)
  # Skewed distribution: lower entropy
  counts_skewed <- c(97, 1, 1, 1)

  h_equal <- shannon_entropy(counts_equal)
  h_skewed <- shannon_entropy(counts_skewed)

  expect_true(h_equal > h_skewed)
  expect_equal(shannon_entropy(c(100, 0, 0, 0)), 0)  # single type = 0 entropy
})

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
