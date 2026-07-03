test_that("ripple_mock_data is the expected SpatialExperiment", {
  skip_if_not_installed("SpatialExperiment")

  data(ripple_mock_data)

  expect_s4_class(ripple_mock_data, "SpatialExperiment")
  expect_equal(nrow(ripple_mock_data), 50)
  expect_equal(ncol(ripple_mock_data), 600)
  expect_true(all(c("cell_type", "sample_id") %in%
    names(SummarizedExperiment::colData(ripple_mock_data))))

  expected_types <- c("Fibroblast", "T_cell", "Tumor")
  expect_setequal(unique(ripple_mock_data$cell_type), expected_types)
  expect_equal(length(unique(ripple_mock_data$sample_id)), 3)
})

test_that("run_ripple recovers the planted gradient on ripple_mock_data", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("meta")

  data(ripple_mock_data)
  out_dir <- tempfile("ripple_mock_test_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE))

  results <- run_ripple(
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

  expect_true(is.data.frame(results))
  expect_true(all(c(
    "gene", "cell_type", "median_coef", "fisher_fdr",
    "sign_consistency"
  ) %in% names(results)))

  # Planted INDUCED genes in T_cell: should all be significant + negative coef
  induced <- results[grepl("^INDUCED", gene) & cell_type == "T_cell"]
  expect_equal(nrow(induced), 5)
  expect_true(all(induced$median_coef < 0),
    info = "INDUCED genes should have negative coefficients"
  )
  expect_true(all(induced$fisher_fdr < 0.01),
    info = "INDUCED genes should be highly significant"
  )
  expect_true(all(induced$sign_consistency == 1.0),
    info = "INDUCED genes should have full sign consistency"
  )

  # Planted REPRESSED genes in T_cell: significant + positive coef
  repressed <- results[grepl("^REPRESSED", gene) & cell_type == "T_cell"]
  expect_equal(nrow(repressed), 5)
  expect_true(all(repressed$median_coef > 0),
    info = "REPRESSED genes should have positive coefficients"
  )
  expect_true(all(repressed$fisher_fdr < 0.01),
    info = "REPRESSED genes should be highly significant"
  )

  # Background genes in T_cell: type I error rate should be low
  bg <- results[grepl("^BG_", gene) & cell_type == "T_cell"]
  fp_rate <- sum(bg$fisher_fdr < 0.05, na.rm = TRUE) / nrow(bg)
  expect_lt(fp_rate, 0.10,
    label = "Background gene false positive rate (T_cell)"
  )
})

test_that("run_ripple warns (not just messages) when cell types are skipped", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("meta")

  data(ripple_mock_data)
  out_dir <- tempfile("ripple_skip_test_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  # Fibroblast has 70 cells/sample; a 90-cell floor leaves it with 0 valid
  # samples, so it must be skipped -- and that must surface as a warning even
  # with verbose = FALSE (the batch/SLURM case).
  expect_warning(
    run_ripple(
      input                = ripple_mock_data,
      query_celltype       = "Tumor",
      celltype_column      = "cell_type",
      sample_column        = "sample_id",
      output_dir           = out_dir,
      target_celltypes     = c("Fibroblast", "T_cell"),
      min_cells_per_sample = 90,
      min_expr_pct         = 0,
      min_expr_floor       = 10,
      verbose              = FALSE
    ),
    "cell type\\(s\\) skipped"
  )
})

test_that("run_ripple warns and returns empty when no cell type qualifies", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("meta")

  data(ripple_mock_data)
  out_dir <- tempfile("ripple_empty_test_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  suppressWarnings(
    res <- expect_warning(
      run_ripple(
        input                = ripple_mock_data,
        query_celltype       = "Tumor",
        celltype_column      = "cell_type",
        sample_column        = "sample_id",
        output_dir           = out_dir,
        min_cells_per_sample = 500,
        min_expr_pct         = 0,
        min_expr_floor       = 10,
        verbose              = FALSE
      ),
      "No cell types had sufficient data"
    )
  )
  expect_true(is.data.frame(res))
  expect_equal(nrow(res), 0)
})

test_that("run_ripple warns on NA cell types but still recovers the gradient", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("meta")

  data(ripple_mock_data)
  spe <- ripple_mock_data

  # Inject NA cell types into a handful of Fibroblasts (not query/target).
  fib_idx <- which(spe$cell_type == "Fibroblast")[1:5]
  spe$cell_type[fib_idx] <- NA

  out_dir <- tempfile("ripple_na_test_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  expect_warning(
    results <- run_ripple(
      input                = spe,
      query_celltype       = "Tumor",
      celltype_column      = "cell_type",
      sample_column        = "sample_id",
      output_dir           = out_dir,
      min_cells_per_sample = 30,
      min_expr_pct         = 0,
      min_expr_floor       = 10,
      verbose              = FALSE
    ),
    "NA in cell-type column"
  )

  # Distances must be finite (no NA-row corruption) and gradient still found.
  induced <- results[grepl("^INDUCED", gene) & cell_type == "T_cell"]
  expect_equal(nrow(induced), 5)
  expect_true(all(induced$median_coef < 0))
  expect_true(all(induced$fisher_fdr < 0.01))
})
