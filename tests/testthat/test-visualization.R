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

  # Mock data has only 2 target cell types, so the default broad_threshold (4)
  # trips the unreachable warning (issue #14); orthogonal to dashboard assembly.

  # Without query markers — bleed-through panel renders without highlight
  qc <- suppressWarnings(ripple_plot_qc(
    results_dir, query_label = "Tumor", top_n_bleed = 5
  ))
  expect_s3_class(qc, "patchwork")

  # With query markers — bleed-through panel highlights matching genes
  qc2 <- suppressWarnings(ripple_plot_qc(
    results_dir,
    query_signature_genes = c("INDUCED_1", "INDUCED_2"),
    query_label = "Tumor",
    top_n_bleed = 5
  ))
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

test_that("plot_gradient_curve auto-detects per-sample mode via sample_id column", {
  set.seed(7)
  bs <- data.table::rbindlist(lapply(c("S1", "S2", "S3"), function(s) {
    data.table::data.table(
      sample_id       = s,
      bin_center      = seq(5, 195, by = 10),
      prop_expressing = pmax(0, seq(0.5, 0.05, length.out = 20) +
                               rnorm(20, 0, 0.02)),
      n_cells         = sample(50:300, 20)
    )
  }))
  p <- plot_gradient_curve(
    bs,
    gene_name = "Cxcl12", cell_type = "T_cell",
    gradient_score = -0.005, fdr = 1e-3
  )
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  # Per-sample mode is signalled by 2 line layers (faint per-sample + bold mean).
  expect_equal(sum(geoms == "GeomLine"), 2)
  expect_equal(sum(geoms == "GeomRibbon"), 1)
})

test_that("plot_gradient_curve drops bins below min_cells_per_bin", {
  bs <- data.table::rbindlist(lapply(c("S1", "S2", "S3"), function(s) {
    data.table::data.table(
      sample_id       = s,
      bin_center      = c(5, 15, 25),
      prop_expressing = c(0.5, 0.3, 0.1),
      n_cells         = c(50, 5, 50)   # middle bin is below the default 10
    )
  }))
  p <- plot_gradient_curve(
    bs, gene_name = "X", cell_type = "Y",
    gradient_score = -0.005, fdr = 1e-3
  )
  ribbon_idx <- which(vapply(p$layers,
                             function(l) class(l$geom)[1], character(1))
                      == "GeomRibbon")
  built <- ggplot2::ggplot_build(p)$data[[ribbon_idx]]
  # Only bins 5 and 25 should survive the filter; bin 15 (n_cells = 5) drops.
  expect_setequal(built$x, c(5, 25))
})

test_that("plot_gradient_curve drops bins below min_samples_per_bin", {
  bs <- data.table::data.table(
    sample_id       = c("S1", "S2", "S1"),
    bin_center      = c(5, 5, 25),       # bin 25 has only 1 sample
    prop_expressing = c(0.5, 0.4, 0.1),
    n_cells         = c(50, 50, 50)
  )
  p <- plot_gradient_curve(
    bs, gene_name = "X", cell_type = "Y",
    gradient_score = -0.005, fdr = 1e-3
  )
  ribbon_idx <- which(vapply(p$layers,
                             function(l) class(l$geom)[1], character(1))
                      == "GeomRibbon")
  built <- ggplot2::ggplot_build(p)$data[[ribbon_idx]]
  expect_setequal(built$x, 5)
})

test_that("bin_decay_data with sample_ids returns per-sample table", {
  set.seed(42)
  counts <- rpois(1000, lambda = 2)
  distances <- runif(1000, 0, 200)
  ids <- rep(c("A", "B"), each = 500)
  out <- bin_decay_data(counts, distances, n_bins = 10, sample_ids = ids)
  expect_true("sample_id" %in% names(out))
  expect_setequal(unique(out$sample_id), c("A", "B"))
  # No bin should have fewer than the default min_cells_per_bin
  expect_true(all(out$n_cells >= 10))
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

test_that("plot_confounder_bar renders + canonicalises legacy class names", {
  s4 <- data.table::data.table(
    gene           = paste0("g", 1:8),
    cell_type      = c("CT1","CT1","CT1","CT2","CT2","CT3","CT3","CT3"),
    classification = c("TRC-Ccl21a_specific", "enhanced", "niche_driven",
                       "underpowered", "no_stage2_result",
                       "TRC-Ccl21a_specific", "confounder_specific",
                       "confounder_driven")
  )
  p <- plot_confounder_bar(s4, query_label = "TRC-Ccl21a",
                           control_label = "TRC-Cxcl12")
  expect_s3_class(p, "ggplot")
  # All raw classes (legacy + dynamic + already-canonical) collapsed to the
  # canonical 5; reversed not present
  cls <- as.character(unique(p$data$classification))
  expect_true(all(cls %in% c("confounder_specific", "enhanced",
                             "confounder_driven", "underpowered",
                             "no_conf_result")))
  # Cell types ordered by total descending: CT1 (3) and CT3 (3) tie before CT2 (2)
  ct_order <- levels(p$data$cell_type)
  expect_equal(ct_order[length(ct_order)], "CT2")
  # Subtitle reflects query/control labels
  expect_match(p$labels$subtitle, "TRC-Ccl21a query")
  expect_match(p$labels$subtitle, "TRC-Cxcl12 control")
})

test_that("plot_confounder_scatter labels requested genes via ggrepel", {
  s4 <- data.table::data.table(
    gene = c("Ccr7", "Sell", "Bg1", "Bg2", "Bg3"),
    stage1_coef        = c(-0.005, -0.004,  0.001,  0.002, -0.003),
    stage2_median_coef = c(-0.004, -0.0035, 0.0001, 0.0015, -0.0001),
    classification     = c("TRC-Ccl21a_specific", "enhanced",
                           "underpowered", "niche_driven",
                           "no_stage2_result"),
    cell_type          = "T_cell_all"
  )
  p <- plot_confounder_scatter(
    s4,
    label_genes   = c("Ccr7", "Sell"),
    query_label   = "TRC-Ccl21a",
    control_label = "TRC-Cxcl12"
  )
  expect_s3_class(p, "ggplot")
  # Layers: 3 ref-lines (abline, hline, vline), 1 point, 1 text-repel
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomTextRepel" %in% geoms)
  txt_idx <- which(geoms == "GeomTextRepel")
  expect_equal(nrow(p$layers[[txt_idx]]$data), 2L)
  # Coloured by canonicalised class
  cls <- as.character(unique(p$data$classification))
  expect_true(all(cls %in% c("confounder_specific", "enhanced",
                             "confounder_driven", "underpowered",
                             "no_conf_result")))
  # No labels when label_genes = NULL
  p2 <- plot_confounder_scatter(s4)
  geoms2 <- vapply(p2$layers, function(l) class(l$geom)[1], character(1))
  expect_false("GeomTextRepel" %in% geoms2)
})

test_that("plot_confounder_* error on missing required columns", {
  expect_error(
    plot_confounder_bar(data.table::data.table(gene = "g1")),
    "Missing required columns"
  )
  expect_error(
    plot_confounder_scatter(data.table::data.table(gene = "g1")),
    "Missing required columns"
  )
  expect_error(
    plot_confounder_ratio(data.table::data.table(gene = "g1")),
    "Missing required columns"
  )
})

test_that("plot_confounder_ratio computes stage2/stage1 and draws thresholds", {
  s4 <- data.table::data.table(
    gene = c("Ccr7", "Sell", "Dropme", "NoStage2", "Lef1"),
    stage1_coef        = c(-0.005, -0.004, 0,       -0.003,  -0.002),
    stage2_median_coef = c(-0.004, -0.0035, -0.001, NA,      0.0025),
    classification     = c("TRC-Ccl21a_specific", "enhanced",
                           "TRC-Ccl21a_specific", "no_stage2_result",
                           "reversed"),
    cell_type          = "T_cell_all"
  )

  expect_message(
    p <- plot_confounder_ratio(
      s4,
      label_genes   = c("Ccr7", "Sell"),
      query_label   = "TRC-Ccl21a",
      control_label = "TRC-Cxcl12"
    ),
    "dropped 2 gene"
  )

  expect_s3_class(p, "ggplot")

  # Layers: 4 hlines + 1 vline + 1 point + 1 text-repel = 7
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_equal(sum(geoms == "GeomHline"),     4)
  expect_equal(sum(geoms == "GeomVline"),     1)
  expect_equal(sum(geoms == "GeomPoint"),     1)
  expect_equal(sum(geoms == "GeomTextRepel"), 1)

  # Ratio computed correctly for the 3 plottable rows
  pt_data <- p$data
  expect_equal(sort(round(pt_data$ratio, 4)),
               sort(round(c(-0.004 / -0.005,
                            -0.0035 / -0.004,
                             0.0025 / -0.002), 4)))

  # Subtitle reports all 5 genes (counts include the dropped ones)
  expect_match(p$labels$subtitle, "TRC-Ccl21a")

  # ggrepel layer should carry exactly the 2 requested label rows
  txt_idx <- which(geoms == "GeomTextRepel")
  expect_equal(nrow(p$layers[[txt_idx]]$data), 2L)

  # Class colour palette is canonicalised (no raw "TRC-Ccl21a_specific")
  cls <- as.character(unique(pt_data$classification))
  expect_true(all(cls %in% c("confounder_specific", "enhanced",
                             "confounder_driven", "underpowered",
                             "reversed", "no_conf_result")))
})

test_that("plot_confounder_ratio errors when no plottable rows remain", {
  s4 <- data.table::data.table(
    gene = c("a", "b"),
    stage1_coef        = c(0, 0),
    stage2_median_coef = c(NA, NA),
    classification     = c("no_stage2_result", "no_stage2_result"),
    cell_type          = "CT"
  )
  expect_error(
    suppressMessages(plot_confounder_ratio(s4)),
    "No genes remain"
  )
})

test_that("plot_k_diagnostics measures distances within each sample", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("S4Vectors")

  # Two samples occupying the SAME coordinate frame. Sample A: query at x=0,
  # its target at x=5. Sample B: query at x=5, its target at x=0. Pooled nn2
  # would match A's target (5,0) to B's query at (5,0) -> distance 0 (wrong);
  # per-sample it must find its own query at distance exactly 5.
  qA <- cbind(0, c(0, 1, 2)) # 3 query cells at x=0 (sample A)
  tA <- cbind(5, 0) # target at x=5 (sample A)
  qB <- cbind(5, c(0, 1, 2)) # 3 query cells at x=5 (sample B)
  tB <- cbind(0, 0) # target at x=0 (sample B)

  coords <- rbind(qA, tA, qB, tB)
  ct <- c(rep("query", 3), "target", rep("query", 3), "target")
  samp <- c(rep("A", 4), rep("B", 4))
  colnames(coords) <- c("x", "y")
  rownames(coords) <- paste0("cell_", seq_len(nrow(coords)))

  counts <- matrix(1L, nrow = 2, ncol = nrow(coords),
    dimnames = list(c("g1", "g2"), rownames(coords)))
  counts <- methods::as(counts, "CsparseMatrix")

  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(cell_type = ct, patient = samp,
      row.names = rownames(coords)),
    spatialCoords = coords
  )

  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  res <- plot_k_diagnostics(
    input = spe, query_celltype = "query", celltype_column = "cell_type",
    sample_column = "patient", k_range = 1:3, verbose = FALSE
  )

  # Each target's nearest within-sample query is exactly 5 um away. A pooled
  # search would instead report 0 (matching the other sample's query).
  target_k1 <- res[cell_type == "target" & k == 1]
  expect_equal(target_k1$mean_dist, 5, tolerance = 1e-6)
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
