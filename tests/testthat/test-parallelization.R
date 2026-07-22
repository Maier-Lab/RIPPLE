# The parallelization vignette recommends fanning out over target cell
# types with future_lapply(), passing a per-celltype subset of the input
# and target_celltypes = <one type>. This test asserts that the pattern
# produces results identical to a single-call run_ripple() with the same
# arguments. If it drifts, the vignette is wrong.

test_that("per-celltype fan-out via future_lapply matches a single run_ripple call", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("future.apply")
  skip_on_cran()

  data(ripple_mock_data)
  spe <- ripple_mock_data

  query   <- "Tumor"
  targets <- setdiff(unique(spe$cell_type), query)

  common_args <- list(
    query_celltype       = query,
    celltype_column      = "cell_type",
    sample_column        = "sample_id",
    min_cells_per_sample = 30,
    min_expr_pct         = 0,
    min_expr_floor       = 10,
    verbose              = FALSE
  )

  # --- Serial baseline ---
  serial_dir <- withr::local_tempdir()
  serial_result <- suppressMessages(do.call(
    run_ripple,
    c(list(input = spe, output_dir = serial_dir), common_args)
  ))

  # --- Parallel fan-out (multisession, 2 workers) ---
  # Wrap in withr so we restore the plan even if the test errors.
  old_plan <- future::plan(future::multisession, workers = 2)
  withr::defer(future::plan(old_plan))

  par_root <- withr::local_tempdir()
  par_list <- future.apply::future_lapply(targets, function(ct) {
    keep <- spe$cell_type %in% c(query, ct)
    spe_sub <- spe[, keep]
    suppressMessages(do.call(ripple::run_ripple, c(
      list(
        input            = spe_sub,
        output_dir       = file.path(par_root, ct),
        target_celltypes = ct,
        analysis_name    = "ripple"
      ),
      common_args
    )))
  }, future.seed = TRUE)
  names(par_list) <- targets
  par_combined <- data.table::rbindlist(par_list, fill = TRUE)

  # --- Compare on gene x cell_type primary key ---
  key <- c("cell_type", "gene")
  data.table::setorderv(serial_result, key)
  data.table::setorderv(par_combined,  key)

  expect_equal(nrow(par_combined), nrow(serial_result))
  expect_identical(par_combined$gene, serial_result$gene)
  expect_identical(
    as.character(par_combined$cell_type),
    as.character(serial_result$cell_type)
  )

  # Numeric outputs must match to a tight tolerance across all cell types.
  for (col in c("median_coef", "fisher_pval", "fisher_fdr",
                "sign_consistency", "gradient_score")) {
    if (col %in% names(serial_result)) {
      expect_equal(
        par_combined[[col]], serial_result[[col]],
        tolerance = 1e-8,
        info = paste0("Column '", col, "' drifted between serial and parallel")
      )
    }
  }
})
