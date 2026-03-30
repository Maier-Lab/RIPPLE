test_that("fit_poisson returns correct structure", {
  # Create synthetic data: 100 cells, expression decreases with distance
  set.seed(42)
  n <- 100
  distances <- runif(n, 0, 200)
  total_counts <- rpois(n, 5000)
  # Gene with planted negative gradient (induced near query)
  lambda <- exp(-3 + (-0.005) * distances + log(total_counts))
  counts <- rpois(n, lambda)

  result <- fit_poisson(counts, distances, total_counts)

  expect_type(result, "list")
  expect_named(result, c("beta", "se", "pval", "dispersion", "n_cells"))
  expect_true(result$beta < 0)  # Should detect negative gradient
  expect_true(result$pval < 0.05)
})

test_that("fit_poisson returns NULL with too few cells", {
  result <- fit_poisson(
    counts = c(rep(0, 95), rep(1, 5)),
    distances = runif(100, 0, 200),
    total_counts = rep(5000, 100),
    min_cells = 25
  )
  expect_null(result)
})

test_that("fit_poisson_controlled returns query-specific coefficient", {
  set.seed(42)
  n <- 200
  dist_query <- runif(n, 0, 200)
  dist_control <- runif(n, 0, 200)
  total_counts <- rpois(n, 5000)
  # Gene driven by query distance only
  lambda <- exp(-3 + (-0.005) * dist_query + 0 * dist_control + log(total_counts))
  counts <- rpois(n, lambda)

  result <- fit_poisson_controlled(counts, dist_query, dist_control, total_counts)

  expect_type(result, "list")
  expect_true(result$beta < 0)
})
