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
  expect_true(result$beta < 0) # Should detect negative gradient
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

test_that("classify_decay_pattern always returns a character scalar", {
  # Regression guard: an earlier bug used `return(NULL)` inside the fit_exp
  # tryCatch block, which exited the whole function (not just the tryCatch)
  # when there were fewer than 3 valid bins. sapply() over samples then
  # produced a list with mixed types and crashed downstream table() calls.
  # The function must ALWAYS return a character scalar for any non-trivial
  # input that passes the initial n-cells/expression checks.

  # Normal case with clear gradient -> should return a meaningful label.
  set.seed(123)
  n <- 200
  distances <- runif(n, 0, 200)
  total_counts <- rpois(n, 5000)
  lambda <- exp(-3 + (-0.01) * distances + log(total_counts))
  counts <- rpois(n, lambda)

  res <- classify_decay_pattern(counts, distances, total_counts)
  expect_type(res, "character")
  expect_length(res, 1L)
  expect_false(is.na(res))

  # Degenerate case: distances all in one narrow window (< 3 bins after
  # cut at 20 um). The fit_exp tryCatch should return NULL internally
  # but classify_decay_pattern must still return a character scalar.
  set.seed(456)
  narrow_distances <- runif(n, 0, 15) # all fit in one 20 um bin
  narrow_counts <- rpois(n, 3)
  narrow_totals <- rpois(n, 5000)
  res_narrow <- classify_decay_pattern(narrow_counts, narrow_distances,
                                       narrow_totals)
  expect_type(res_narrow, "character")
  expect_length(res_narrow, 1L)

  # Edge case: flat data (no gradient) -> should still return a character,
  # whether "linear" (flat slope) or similar, never NULL.
  set.seed(789)
  flat_counts <- rpois(n, 5)
  flat_distances <- runif(n, 0, 200)
  flat_totals <- rpois(n, 5000)
  res_flat <- classify_decay_pattern(flat_counts, flat_distances, flat_totals)
  expect_type(res_flat, "character")
  expect_length(res_flat, 1L)
})
