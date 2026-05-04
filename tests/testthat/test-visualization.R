test_that("plot_gene_category_dotplot builds a ggplot from curated panel", {
  skip_if_not_installed("ggplot2")

  results <- data.table::data.table(
    gene = c("G1", "G2", "G3", "G4", "G5", "G6"),
    median_coef = c(-0.005, -0.001, 0.003, 0.002, -0.004, 0.001),
    fisher_fdr = c(1e-10, 0.2, 1e-5, 0.8, 1e-20, 0.03)
  )

  panel <- list(
    "Category A" = c("G1", "G2"),
    "Category B" = c("G3", "G4"),
    "Category C" = c("G5", "G6")
  )

  p <- plot_gene_category_dotplot(results, panel, query_label = "Test")

  expect_s3_class(p, "ggplot")
  expect_true("category" %in% names(p$data))
  expect_setequal(
    as.character(levels(p$data$category)),
    c("Category A", "Category B", "Category C")
  )
  expect_equal(nrow(p$data), 6)
  expect_true(all(p$data$is_sig == (p$data$fisher_fdr < 0.05)))
})

test_that("plot_gene_category_dotplot drops missing genes and warns when empty", {
  results <- data.table::data.table(
    gene = c("G1", "G2"),
    median_coef = c(-0.001, 0.001),
    fisher_fdr = c(0.01, 0.5)
  )

  # Genes not in results — should error
  expect_error(
    plot_gene_category_dotplot(
      results,
      list("Missing" = c("Z1", "Z2"))
    ),
    "None of the genes"
  )

  # Partial match — keeps overlap only
  p <- plot_gene_category_dotplot(
    results,
    list("A" = c("G1", "Z_absent"), "B" = c("G2"))
  )
  expect_equal(nrow(p$data), 2)
})

test_that("plot_gene_category_dotplot validates inputs", {
  results <- data.table::data.table(
    gene = "G1", median_coef = 0, fisher_fdr = 0.5
  )

  expect_error(
    plot_gene_category_dotplot(results, c("G1", "G2")),
    "named list"
  )
  expect_error(
    plot_gene_category_dotplot(results, list(c("G1"))),
    "named list"
  )
  expect_error(
    plot_gene_category_dotplot(
      data.table::data.table(gene = "G1"),
      list("A" = "G1")
    ),
    "Required column"
  )
})


test_that("ripple_plot_qc assembles a dashboard from run_ripple output", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("meta")
  skip_if_not_installed("patchwork")

  data(ripple_mock_data)
  out_dir <- tempfile("ripple_qc_test_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE))

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

  # Pipeline should write the new per-cell distance file
  expect_true(file.exists(file.path(
    results_dir, "qc", "cell_distances.csv.gz"
  )))

  # Without query markers — bleed-through panel renders without highlight
  qc <- ripple_plot_qc(
    results_dir, query_label = "Tumor", top_n_bleed = 5
  )
  expect_s3_class(qc, "patchwork")

  # With query markers — bleed-through panel highlights matching genes
  qc2 <- ripple_plot_qc(
    results_dir,
    query_signature_genes = c("INDUCED_1", "INDUCED_2"),
    query_label = "Tumor",
    top_n_bleed = 5
  )
  expect_s3_class(qc2, "patchwork")
})


test_that("plot_gradient_curve renders pooled mode", {
  bin_stats <- data.table::data.table(
    dist_mid        = seq(5, 195, by = 10),
    prop_expressing = seq(0.5, 0.05, length.out = 20),
    n_cells         = round(seq(200, 50, length.out = 20)),
    se              = rep(0.02, 20)
  )
  p <- plot_gradient_curve(
    bin_stats,
    gene_name      = "Cxcl12", cell_type = "T_cell",
    gradient_score = -0.005,    fdr = 1e-3
  )
  expect_s3_class(p, "ggplot")
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomLine"   %in% geoms)
  expect_true("GeomRibbon" %in% geoms)
  expect_equal(sum(geoms == "GeomLine"), 1)
})

test_that("plot_gradient_curve per-sample mode draws bold mean + 95% CI ribbon", {
  set.seed(7)
  samples <- c("S1", "S2", "S3", "S4")
  bin_per_sample <- data.table::rbindlist(lapply(samples, function(s) {
    data.table::data.table(
      sample_id = s,
      bin_mid   = seq(5, 195, by = 10),
      mean_rate = pmax(
        0,
        seq(0.004, 0.001, length.out = 20) + rnorm(20, 0, 5e-4)
      ),
      n_cells   = sample(50:300, 20)
    )
  }))

  p <- plot_gradient_curve(
    bin_per_sample,
    gene_name      = "Cxcl12", cell_type = "T_cell",
    gradient_score = -0.003,    fdr = 1e-4,
    sample_col     = "sample_id",
    x_col          = "bin_mid",
    y_col          = "mean_rate",
    y_lab          = "Mean expression rate"
  )

  expect_s3_class(p, "ggplot")
  # 2 geom_line (per-sample faint + bold mean), 1 ribbon (95% CI), 1 point
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_equal(sum(geoms == "GeomLine"),   2)
  expect_equal(sum(geoms == "GeomRibbon"), 1)
  expect_equal(sum(geoms == "GeomPoint"),  1)

  # The ribbon should encode the 95% CI = mean +/- 1.96 * SE across samples.
  # Compute expected widths and compare with the layer's mapped data.
  expected <- bin_per_sample[, .(
    mean_val = mean(mean_rate),
    se_val   = stats::sd(mean_rate) / sqrt(.N)
  ), by = bin_mid][order(bin_mid)]
  expected[, ymin := pmax(0, mean_val - 1.96 * se_val)]
  expected[, ymax := mean_val + 1.96 * se_val]

  # Locate the ribbon layer and force its computed positions
  ribbon_layer <- p$layers[[which(geoms == "GeomRibbon")]]
  built <- ggplot2::ggplot_build(p)$data[[which(geoms == "GeomRibbon")]]
  built <- built[order(built$x), ]
  expect_equal(built$ymin, expected$ymin, tolerance = 1e-9)
  expect_equal(built$ymax, expected$ymax, tolerance = 1e-9)
})

test_that("plot_gradient_curve errors when sample_col is missing", {
  bs <- data.table::data.table(
    dist_mid = 1:5, prop_expressing = c(0.4, 0.3, 0.2, 0.15, 0.1)
  )
  expect_error(
    plot_gradient_curve(bs, "G", "CT", sample_col = "nonexistent"),
    "Missing required columns"
  )
})

test_that("plot_prop_curve is a thin wrapper with proportion-expressing defaults", {
  set.seed(11)
  samples <- c("S1", "S2", "S3")
  bin_per_sample <- data.table::rbindlist(lapply(samples, function(s) {
    data.table::data.table(
      sample_id       = s,
      dist_mid        = seq(5, 195, by = 10),
      prop_expressing = pmin(
        1,
        pmax(0, seq(0.4, 0.05, length.out = 20) + rnorm(20, 0, 0.02))
      ),
      n_cells         = sample(50:300, 20)
    )
  }))

  p <- plot_prop_curve(
    bin_per_sample,
    gene_name      = "Cxcl12", cell_type = "T_cell",
    gradient_score = -0.005,    fdr = 1e-4,
    sample_col     = "sample_id"
  )
  expect_s3_class(p, "ggplot")

  # Defaults from the wrapper should appear on the built plot
  built <- ggplot2::ggplot_build(p)
  expect_equal(p$labels$y, "Proportion expressing")
  # Should have 2 lines (per-sample faint + bold mean), 1 ribbon, 1 point
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_equal(sum(geoms == "GeomLine"),   2)
  expect_equal(sum(geoms == "GeomRibbon"), 1)
  expect_equal(sum(geoms == "GeomPoint"),  1)
})

test_that("plot_decay_curve still works as a deprecated alias", {
  bs <- data.table::data.table(
    dist_mid        = seq(5, 195, by = 10),
    prop_expressing = seq(0.5, 0.05, length.out = 20),
    n_cells         = round(seq(200, 50, length.out = 20)),
    se              = rep(0.02, 20)
  )
  expect_warning(
    p <- plot_decay_curve(
      bs,
      gene_name      = "X", cell_type = "T_cell",
      gradient_score = -0.005, fdr = 1e-3
    ),
    "deprecated", ignore.case = TRUE
  )
  expect_s3_class(p, "ggplot")
})

test_that("ripple_plot_qc errors when summary file is missing", {
  bad <- tempfile("ripple_qc_bad_")
  dir.create(bad)
  on.exit(unlink(bad, recursive = TRUE))
  dir.create(file.path(bad, "summary"))
  expect_error(
    ripple_plot_qc(bad),
    "Missing required file"
  )
})
