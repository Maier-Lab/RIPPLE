# Smoke test for run_ripple_atlas: the largest orchestrator in the package,
# previously exercised by no test (which is how the msigdbr break stayed
# invisible). This drives the non-fGSEA panels end-to-end on real run_ripple
# output and asserts the wiring holds -- files are produced and the documented
# return value is a character vector of paths.

test_that("run_ripple_atlas produces panels from run_ripple output", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("meta")

  data(ripple_mock_data)
  out_dir <- tempfile("ripple_atlas_test_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  # Produce a real results directory (summary/all_genes_results.csv).
  run_ripple(
    input                = ripple_mock_data,
    query_celltype       = "Tumor",
    celltype_column      = "cell_type",
    sample_column        = "sample_id",
    output_dir           = out_dir,
    min_cells_per_sample = 30,
    min_expr_pct         = 0,
    min_expr_floor       = 10,
    verbose              = FALSE
  )
  results_dir <- file.path(out_dir, "ripple")
  expect_true(file.exists(
    file.path(results_dir, "summary", "all_genes_results.csv")
  ))

  # Run the atlas without fGSEA (no network / msigdbr dependency).
  saved <- suppressWarnings(run_ripple_atlas(
    results_dir = results_dir,
    query_label = "Tumor",
    run_fgsea   = FALSE,
    verbose     = FALSE
  ))

  # Documented return value: a character vector of generated file paths.
  expect_type(saved, "character")
  expect_gt(length(saved), 0)
  expect_true(all(file.exists(saved)))

  # Plots land in the default results_dir/plots directory.
  plots_dir <- file.path(results_dir, "plots")
  expect_true(dir.exists(plots_dir))
  expect_gt(length(list.files(plots_dir)), 0)
})

test_that("run_ripple_atlas errors clearly when the summary file is missing", {
  bad <- tempfile("ripple_atlas_bad_")
  dir.create(file.path(bad, "summary"), recursive = TRUE)
  on.exit(unlink(bad, recursive = TRUE), add = TRUE)

  expect_error(
    run_ripple_atlas(results_dir = bad, run_fgsea = FALSE, verbose = FALSE),
    "all_genes_results.csv"
  )
})
