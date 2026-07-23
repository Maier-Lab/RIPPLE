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

  # --- consolidate_parallel_ripple round-trip ---
  # The worker output_dir was file.path(par_root, ct); run_ripple wrote to
  # <output_dir>/<analysis_name>/ = <par_root>/<ct>/ripple/
  worker_dirs <- file.path(par_root, targets, "ripple")
  expect_true(all(dir.exists(worker_dirs)))

  combined_dir <- file.path(par_root, "combined")
  suppressMessages(
    consolidate_parallel_ripple(
      worker_dirs    = worker_dirs,
      output_dir     = combined_dir,
      query_celltype = query
    )
  )

  # Consolidated tree has the shape ripple_plot_qc expects
  expect_true(file.exists(file.path(combined_dir, "summary",
                                    "all_genes_results.csv")))
  expect_true(file.exists(file.path(combined_dir, "qc",
                                    "cell_distances.csv.gz")))
  for (ct in targets) {
    expect_true(dir.exists(file.path(combined_dir, "per_celltype", ct)))
  }

  # Query rows should appear ONCE per (sample, dist) after consolidation,
  # not N_workers times. Compare to any single worker's query row count.
  merged_dist <- data.table::fread(
    file.path(combined_dir, "qc", "cell_distances.csv.gz")
  )
  one_worker_dist <- data.table::fread(
    file.path(worker_dirs[1], "qc", "cell_distances.csv.gz")
  )
  n_query_merged <- sum(merged_dist$cell_type == query)
  n_query_one    <- sum(one_worker_dist$cell_type == query)
  expect_equal(n_query_merged, n_query_one,
    info = "Query cells duplicated across workers instead of deduped"
  )

  # Consolidated summary matches serial baseline on gene x celltype x FDR
  consolidated <- data.table::fread(
    file.path(combined_dir, "summary", "all_genes_results.csv")
  )
  data.table::setorderv(consolidated, key)
  expect_equal(nrow(consolidated), nrow(serial_result))
  expect_equal(
    consolidated$fisher_fdr, serial_result$fisher_fdr, tolerance = 1e-8
  )
})
