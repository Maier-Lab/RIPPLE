test_that("ripple_config gets and sets options with the ripple. prefix", {
  old <- ripple_config(k_neighbors = 7L)
  on.exit(do.call(ripple_config, old), add = TRUE)

  expect_equal(ripple_config("k_neighbors"), 7L)
  expect_equal(getOption("ripple.k_neighbors"), 7L)

  all_opts <- ripple_config()
  expect_true("k_neighbors" %in% names(all_opts))
  expect_false(any(grepl("^ripple\\.", names(all_opts))))
})

test_that("ripple_config warns on unknown (misspelled) option names", {
  expect_warning(
    old <- ripple_config(max_distanceum = 500),
    "Unknown RIPPLE option"
  )
  on.exit(options(ripple.max_distanceum = NULL), add = TRUE)
})

test_that(".env_num returns default for unset or malformed values", {
  # Unset -> default, no warning
  expect_silent(v1 <- ripple:::.env_num("RIPPLE_NO_SUCH_VAR_XYZ", 42))
  expect_equal(v1, 42)

  withr::local_envvar(RIPPLE_TEST_NUM = "not_a_number")
  expect_warning(
    v2 <- ripple:::.env_num("RIPPLE_TEST_NUM", 200),
    "not a valid number"
  )
  expect_equal(v2, 200)

  withr::local_envvar(RIPPLE_TEST_NUM = "3.5")
  expect_equal(ripple:::.env_num("RIPPLE_TEST_NUM", 200), 3.5)
})

test_that(".resolve consults the package option when no explicit value is passed", {
  withr::local_options(ripple.max_distance_um = 123)
  # Explicit value wins over option
  expect_equal(ripple:::.resolve(50, "max_distance_um", 200), 50)
  # NULL -> option value, not the hardcoded default
  expect_equal(ripple:::.resolve(NULL, "max_distance_um", 200), 123)
})

test_that("run_ripple honors options() for an omitted argument (regression)", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("meta")

  data(ripple_mock_data)

  run_one <- function(...) {
    out_dir <- tempfile("ripple_cfg_test_")
    dir.create(out_dir)
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
    run_ripple(
      input                = ripple_mock_data,
      query_celltype       = "Tumor",
      celltype_column      = "cell_type",
      sample_column        = "sample_id",
      output_dir           = out_dir,
      min_cells_per_sample = 30,
      min_expr_pct         = 0,
      min_expr_floor       = 10,
      verbose              = FALSE,
      ...
    )
  }

  # Argument omitted, but option set to 60: must match passing 60 explicitly,
  # and must differ from the built-in default (200). Before the fix, the
  # omitted-argument run silently used 200 and matched the default instead.
  res_explicit_60 <- run_one(max_distance_um = 60)
  res_default_200 <- run_one(max_distance_um = 200)

  withr::local_options(ripple.max_distance_um = 60)
  res_option_60 <- run_one()

  key <- c("gene", "cell_type", "median_coef")
  expect_equal(res_option_60[, ..key], res_explicit_60[, ..key])
  expect_false(isTRUE(all.equal(
    res_option_60[, ..key], res_default_200[, ..key]
  )))
})
