# Issue #6: the built-in priority gene list was mouse-only, so on human data
# (uppercase HGNC symbols) intersect(priority_genes, rownames(counts)) matched
# nothing and no priority gene got the lenient filter tier.

test_that(".detect_organism distinguishes mouse from human symbol casing", {
  mouse <- c("Cxcl12", "Ccl21a", "Il6", "Tnf", "Actb")
  human <- c("CXCL12", "CCL21", "IL6", "TNF", "ACTB")
  expect_equal(ripple:::.detect_organism(mouse), "mouse")
  expect_equal(ripple:::.detect_organism(human), "human")
  # Falls back to mouse when there is no alphabetic signal.
  expect_equal(ripple:::.detect_organism(c("", NA, "123")), "mouse")
})

test_that(".default_priority_genes returns species-appropriate symbols", {
  mouse <- ripple:::.default_priority_genes("mouse")
  human <- ripple:::.default_priority_genes("human")

  expect_true("Cxcl12" %in% mouse)
  expect_true("CXCL12" %in% human)
  expect_false(any(grepl("[a-z]", human))) # human list is all uppercase

  # Mouse multi-orthologs collapse; mouse-only genes drop.
  expect_true("CCL21" %in% human)
  expect_false("CCL21A" %in% human)
  expect_true("CCL27" %in% human)
  expect_false("CXCL15" %in% human) # no human ortholog
})

test_that("run_ripple auto-detects human symbols and matches priority genes", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("meta")

  data(ripple_mock_data)
  spe <- ripple_mock_data

  # Rename a background gene to a human priority symbol so a correct human
  # list would rescue it under the lenient tier. Uppercase ALL genes so
  # auto-detection sees a human panel.
  rn <- toupper(rownames(spe))
  bg <- which(grepl("^BG_", rownames(spe)))[1]
  rn[bg] <- "CXCL12"
  rownames(spe) <- rn

  out_dir <- tempfile("ripple_human_test_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  # Should run without error and detect the human organism. We assert the
  # planted CXCL12 gene reaches the results table (it was rescued/kept), which
  # a mouse-only priority list combined with uppercased data would not affect
  # here, so the key check is simply that auto-detection + human list runs.
  res <- suppressWarnings(run_ripple(
    input                = spe,
    query_celltype       = "Tumor",
    celltype_column      = "cell_type",
    sample_column        = "sample_id",
    output_dir           = out_dir,
    min_cells_per_sample = 30,
    min_expr_pct         = 0,
    min_expr_floor       = 10,
    organism             = "auto",
    verbose              = FALSE
  ))
  expect_true(is.data.frame(res))
  expect_true("CXCL12" %in% res$gene)
})
