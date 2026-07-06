# Integration tests for downstream pipeline stages:
#   - Stage 3: merge_ripple_results()
#   - Stage 4: run_ripple_confounder()
#   - Stage 5: classify_gene_specificity(), run_ripple_fgsea()
#   - Stage 6: classify_lr_artifacts()

test_that("classify_gene_specificity classifies genes by cell-type breadth", {
  results <- data.table::data.table(
    gene = c(
      "SPECIFIC_A", "SPECIFIC_A",    # sig in 1 ct only
      "MODERATE_B", "MODERATE_B", "MODERATE_B",
      "WIDE_C", "WIDE_C", "WIDE_C",
      "BROAD_D", "BROAD_D", "BROAD_D", "BROAD_D", "BROAD_D"
    ),
    cell_type = c(
      "CT1", "CT2",
      "CT1", "CT2", "CT3",
      "CT1", "CT2", "CT3",
      "CT1", "CT2", "CT3", "CT4", "CT5"
    ),
    fisher_fdr = c(
      0.001, 0.5,                      # SPECIFIC_A: sig in CT1 only
      0.001, 0.001, 0.02,              # MODERATE_B: sig in 3 ct
      0.001, 0.001, 0.001,             # WIDE_C: sig in 3 ct (boundary)
      0.001, 0.001, 0.001, 0.001, 0.01 # BROAD_D: sig in 5 ct
    )
  )

  spec <- classify_gene_specificity(
    results,
    fdr_threshold = 0.05, broad_threshold = 4
  )

  expect_s3_class(spec, "data.table")
  # Only three classes exist now; "ubiquitous" was removed.
  expect_true(all(spec$specificity_class %in% c("specific", "moderate", "broad")))
  expect_equal(
    spec[gene == "SPECIFIC_A"]$specificity_class, "specific"
  )
  expect_equal(
    spec[gene == "MODERATE_B"]$specificity_class, "moderate"
  )
  # WIDE_C is sig in 3 ct — below broad_threshold, so "moderate"
  expect_equal(
    spec[gene == "WIDE_C"]$specificity_class, "moderate"
  )
  expect_equal(
    spec[gene == "BROAD_D"]$specificity_class, "broad"
  )
  expect_equal(spec[gene == "BROAD_D"]$n_celltypes, 5L)
})

test_that("raising broad_threshold demotes genes from broad to moderate", {
  # WIDE is significant in exactly 4 cell types.
  results <- data.table::data.table(
    gene = rep("WIDE", 4),
    cell_type = c("CT1", "CT2", "CT3", "CT4"),
    fisher_fdr = rep(0.001, 4)
  )

  # broad_threshold = 4 -> WIDE (4 cell types) is broad.
  s4 <- classify_gene_specificity(results, broad_threshold = 4)
  expect_equal(s4[gene == "WIDE"]$specificity_class, "broad")

  # Raise broad_threshold to 6 -> 4 cell types no longer clears the bar, so
  # WIDE is demoted to moderate. broad_threshold is the boundary (its name).
  # (6 > the 4 cell types present trips the #14 unreachable warning, which is
  # expected here and orthogonal to the demotion behaviour under test.)
  s6 <- suppressWarnings(classify_gene_specificity(results, broad_threshold = 6))
  expect_equal(s6[gene == "WIDE"]$specificity_class, "moderate")
  # And no "ubiquitous" phantom tier appears in between.
  expect_false("ubiquitous" %in% s6$specificity_class)
})

test_that("classify_gene_specificity respects broad_sig_threshold", {
  # Gene SHARED is significant in 4 cell types at FDR < 0.05 (default
  # broad-class call), but only at strong significance (FDR < 0.001) in
  # ONE of those 4 cell types. Tightening the significance bar should
  # therefore demote it from "broad" to a non-broad class.
  results <- data.table::data.table(
    gene = c(
      "SHARED", "SHARED", "SHARED", "SHARED",
      "SPECIFIC_A"
    ),
    cell_type = c("CT1", "CT2", "CT3", "CT4", "CT1"),
    fisher_fdr = c(
      1e-5,   # CT1: ** + *** + **** all hit
      0.02,   # CT2: only * hits
      0.03,   # CT3: only * hits
      0.04,   # CT4: only * hits
      1e-6    # specific control
    )
  )

  # Default behaviour (any * counts) -> SHARED hits 4 cell types -> broad
  spec_default <- classify_gene_specificity(
    results, fdr_threshold = 0.05, broad_threshold = 4
  )
  expect_equal(spec_default[gene == "SHARED"]$specificity_class,
               "broad")
  expect_equal(spec_default[gene == "SHARED"]$n_celltypes, 4L)
  expect_equal(spec_default[gene == "SHARED"]$n_celltypes_strict, 4L)

  # Tighten to ** (FDR < 0.01) -> SHARED only hits CT1 strictly -> NOT
  # broad. With "ubiquitous" removed, a widely-but-weakly significant gene
  # that misses the strict broad bar is classified "moderate".
  spec_strict <- classify_gene_specificity(
    results, fdr_threshold = 0.05,
    broad_threshold     = 4,
    broad_sig_threshold = "**"
  )
  expect_equal(spec_strict[gene == "SHARED"]$specificity_class,
               "moderate")
  expect_equal(spec_strict[gene == "SHARED"]$n_celltypes,        4L)
  expect_equal(spec_strict[gene == "SHARED"]$n_celltypes_strict, 1L)

  # Numeric threshold matches the star convention
  spec_numeric <- classify_gene_specificity(
    results, fdr_threshold = 0.05,
    broad_threshold     = 4,
    broad_sig_threshold = 0.01
  )
  expect_equal(spec_numeric[gene == "SHARED"]$specificity_class,
               spec_strict[gene == "SHARED"]$specificity_class)

  # Bad string -> error
  expect_error(
    classify_gene_specificity(results,
                              broad_sig_threshold = "bogus"),
    "must be NULL"
  )

  # Looser sig_threshold than fdr_threshold -> warning
  expect_warning(
    classify_gene_specificity(results, fdr_threshold = 0.01,
                              broad_sig_threshold = 0.05),
    "MORE aggressive"
  )
})

test_that("classify_gene_specificity handles empty input gracefully", {
  empty <- data.table::data.table(
    gene = character(), cell_type = character(), fisher_fdr = numeric()
  )
  spec <- classify_gene_specificity(empty)
  expect_s3_class(spec, "data.table")
  expect_equal(nrow(spec), 0)
})

test_that("classify_gene_specificity warns when broad_threshold is unreachable", {
  # Only 2 cell types present but broad_threshold = 4 -> no gene can ever be
  # classified 'broad', so the contamination heuristic is inactive (issue #14).
  results <- data.table::data.table(
    gene      = c("G1", "G1", "G2", "G2"),
    cell_type = c("CT1", "CT2", "CT1", "CT2"),
    fisher_fdr = c(0.001, 0.001, 0.4, 0.001)
  )
  expect_warning(
    spec <- classify_gene_specificity(results, broad_threshold = 4),
    "broad_threshold \\(4\\) exceeds the number of cell types"
  )
  # No gene classified broad; function still returns a valid table.
  expect_false("broad" %in% spec$specificity_class)

  # Reachable threshold -> no warning.
  expect_no_warning(
    classify_gene_specificity(results, broad_threshold = 2)
  )
})

test_that("classify_lr_artifacts flags LR pairs by rule", {
  lr <- data.table::data.table(
    ligand = c("Ccl5", "Il6", "Csf3", "Tnf", "Cxcl12"),
    receptor = c("Ccr5", "Il6r", "Csf3r", "Tnfrsf1a", "Cxcr4"),
    cell_type = c(
      "T_cell",          # receptor Ccr5 — clean
      "Macrophage",      # receptor Il6r on myeloid — clean
      "T_cell",          # receptor Csf3r on non-myeloid — ARTIFACT
      "Fibroblast",      # receptor Tnfrsf1a; low receptor_pct — low_confidence
      "T_cell"           # receptor in contam list, non-myeloid, low pct — suspect
    ),
    receptor_pct_target = c(0.30, 0.40, 0.15, 0.01, 0.03)
  )

  out <- classify_lr_artifacts(
    lr,
    query_signature = c("Csf3r", "Ly6g"),
    contamination_genes = c("Cxcr4", "Igkc"),
    myeloid_celltypes = c("Macrophage", "Monocyte", "cDC1"),
    low_expr_threshold = 0.02
  )

  expect_s3_class(out, "data.table")
  expect_equal(
    out[ligand == "Ccl5"]$artifact_flag, "clean"
  )
  expect_equal(
    out[ligand == "Il6"]$artifact_flag, "clean"
  )
  expect_equal(
    out[ligand == "Csf3"]$artifact_flag, "artifact"
  )
  expect_equal(
    out[ligand == "Tnf"]$artifact_flag, "low_confidence"
  )
  expect_equal(
    out[ligand == "Cxcl12"]$artifact_flag, "suspect"
  )
})

test_that("classify_lr_artifacts returns empty on empty input", {
  empty <- data.table::data.table(
    ligand = character(), receptor = character(), cell_type = character()
  )
  out <- classify_lr_artifacts(empty)
  expect_equal(nrow(out), 0)
  expect_true("artifact_flag" %in% names(out))
})

test_that("run_ripple_fgsea runs with a custom gene_sets list", {
  skip_if_not_installed("fgsea")

  # Use a custom named list so we don't need msigdbr
  set.seed(1)
  results <- data.table::data.table(
    gene = c(paste0("INDUCED_", 1:5), paste0("REPRESSED_", 1:5),
             paste0("BG_", 1:20)),
    cell_type = "T_cell",
    # Jitter so fgsea's preranked stat has no ties
    median_coef = c(
      -0.01 + rnorm(5, 0, 5e-4),
      0.01 + rnorm(5, 0, 5e-4),
      rnorm(20, 0, 5e-4)
    ),
    fisher_fdr = c(rep(1e-5, 10), runif(20, 0.2, 1))
  )

  gene_sets <- list(
    INDUCED_SET = paste0("INDUCED_", 1:5),
    REPRESSED_SET = paste0("REPRESSED_", 1:5)
  )

  fgsea_res <- suppressMessages(run_ripple_fgsea(
    results,
    gene_sets = gene_sets,
    coef_col = "median_coef",
    min_size = 3, max_size = 100,
    min_genes = 10,
    n_perm = 1000,
    seed = 42
  ))

  expect_s3_class(fgsea_res, "data.table")
  expect_true(all(
    c("pathway", "pval", "padj", "NES", "cell_type") %in% names(fgsea_res)
  ))
  expect_equal(
    sort(unique(fgsea_res$pathway)),
    c("INDUCED_SET", "REPRESSED_SET")
  )
  # INDUCED_SET has negative coefficients — negative NES (depleted at top of
  # sorted stat, enriched at bottom); REPRESSED_SET the opposite.
  induced_nes <- fgsea_res[pathway == "INDUCED_SET"]$NES
  repressed_nes <- fgsea_res[pathway == "REPRESSED_SET"]$NES
  expect_true(induced_nes < 0)
  expect_true(repressed_nes > 0)
})


test_that("run_ripple_fgsea exclude_broad drops cross-cell-type genes", {
  skip_if_not_installed("fgsea")

  # 5 cell types; BROAD_* are significant in all 5 (>= threshold = 4).
  # SPECIFIC_* are significant only in T_cell.
  set.seed(2)
  ct_levels <- c("T_cell", "Macrophage", "Fibroblast", "B_cell", "DC")
  results <- data.table::rbindlist(lapply(ct_levels, function(ct) {
    data.table::data.table(
      gene        = c(paste0("BROAD_",    1:5),
                      paste0("SPECIFIC_", 1:5),
                      paste0("BG_",       1:20)),
      cell_type   = ct,
      median_coef = c(
        -0.01 + stats::rnorm(5, 0, 5e-4),
        if (ct == "T_cell") -0.01 + stats::rnorm(5, 0, 5e-4)
                            else stats::rnorm(5, 0, 5e-4),
        stats::rnorm(20, 0, 5e-4)
      ),
      fisher_fdr  = c(
        rep(1e-5, 5),                                  # BROAD_*: sig everywhere
        if (ct == "T_cell") rep(1e-5, 5) else runif(5, 0.2, 1),
        runif(20, 0.2, 1)
      )
    )
  }))

  gene_sets <- list(
    BROAD_SET    = paste0("BROAD_",    1:5),
    SPECIFIC_SET = paste0("SPECIFIC_", 1:5)
  )

  # Without filter
  unfiltered <- suppressMessages(run_ripple_fgsea(
    results, gene_sets = gene_sets, coef_col = "median_coef",
    min_size = 3, max_size = 100, min_genes = 10,
    n_perm = 1000, seed = 42
  ))

  # With filter — BROAD_* should be dropped
  filtered <- suppressMessages(run_ripple_fgsea(
    results, gene_sets = gene_sets, coef_col = "median_coef",
    min_size = 3, max_size = 100, min_genes = 10,
    n_perm = 1000, seed = 42,
    exclude_broad   = TRUE,
    broad_threshold = 4L
  ))

  # Unfiltered: BROAD_SET tested in every cell type
  expect_setequal(unique(unfiltered$pathway), c("BROAD_SET", "SPECIFIC_SET"))

  # Filtered: BROAD_SET should NOT be present (no genes left for that set)
  expect_false("BROAD_SET" %in% unique(filtered$pathway))
  expect_true("SPECIFIC_SET" %in% unique(filtered$pathway))
})


# -----------------------------------------------------------------------------
# End-to-end: run Stage 1 once, then exercise merge + confounder on that output.
# Grouped into one test_that block so Stage 1 only runs once.
# -----------------------------------------------------------------------------

test_that("downstream stages chain off a single run_ripple run", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("meta")

  data(ripple_mock_data)
  out_dir <- tempfile("ripple_stages_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE))

  stage1 <- run_ripple(
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
  expect_true(dir.exists(results_dir))
  expect_true(file.exists(file.path(
    results_dir, "summary", "all_genes_results.csv"
  )))

  # --- merge_ripple_results ---
  merged <- suppressMessages(merge_ripple_results(
    results_dir,
    recompute_fisher = TRUE,
    verbose = FALSE
  ))
  expect_s3_class(merged, "data.table")
  expect_true(all(c("gene", "cell_type", "fisher_fdr") %in% names(merged)))
  expect_setequal(unique(merged$cell_type), c("Fibroblast", "T_cell"))
  # Planted INDUCED / REPRESSED should survive merge for T_cell
  induced <- merged[grepl("^INDUCED", gene) & cell_type == "T_cell"]
  expect_equal(nrow(induced), 5)
  expect_true(all(induced$fisher_fdr < 0.01))

  # --- classify_gene_specificity on merged output ---
  # Mock data has only 2 target cell types, so the default broad_threshold (4)
  # legitimately triggers the "unreachable" warning (issue #14); suppress it
  # here since it is orthogonal to the classification behaviour under test.
  spec <- suppressWarnings(
    classify_gene_specificity(merged, fdr_threshold = 0.05)
  )
  expect_s3_class(spec, "data.table")
  # Planted signal is T_cell-specific — the 10 planted genes should appear
  # as either specific or moderate (not broad)
  planted <- spec[grepl("^INDUCED|^REPRESSED", gene)]
  expect_true(nrow(planted) >= 1)
  expect_true(!any(planted$specificity_class == "broad"))

  # --- run_ripple_confounder with Fibroblast as control ---
  stage2 <- suppressMessages(run_ripple_confounder(
    input            = ripple_mock_data,
    results_dir      = results_dir,
    query_celltype   = "Tumor",
    celltype_column  = "cell_type",
    control_celltype = "Fibroblast",
    sample_column    = "sample_id",
    target_celltypes = "T_cell",
    min_cells_per_sample = 30,
    min_control_cells = 20,
    verbose          = FALSE
  ))
  expect_s3_class(stage2, "data.table")
  expect_true("gene" %in% names(stage2))
  # Stage 2 writes to "<results_dir>_stage2" by default
  stage2_sum <- file.path(
    paste0(results_dir, "_stage2"), "summary", "stage2_all_results.csv"
  )
  expect_true(file.exists(stage2_sum))
})
