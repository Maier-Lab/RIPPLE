# ============================================================================
# Method comparison benchmark: RIPPLE vs MISTy vs nnSVG
# ============================================================================
# Runs MISTy and nnSVG on the same synthetic datasets used by bench_power.R
# and compares their outputs to RIPPLE.
#
# Claims tested:
#   1. nnSVG detects the planted gradient genes (positive control) but also
#      flags many additional spatially variable genes unrelated to the query
#      anchor. This reflects a specificity gap: an SVG test picks up any
#      spatial structure, while RIPPLE conditions on distance to a specific
#      cell type.
#   2. MISTy variable importance (celltype-aware para-view anchored on the
#      query type) correlates with RIPPLE |beta|. MISTy and RIPPLE largely
#      agree on which genes have a signal but MISTy does not provide
#      direction, per-sample replication, or a calibrated FDR.
#
# Design: fixed effect size (beta = -0.01) x 5 samples x 10 iterations.
# A single config is used instead of the full bench_power grid because
# MISTy and nnSVG are considerably more expensive to run than RIPPLE.
#
# Run with:
#   Rscript data-raw/benchmarks/bench_comparison.R
#
# Output:
#   data-raw/benchmarks/results/bench_comparison_results.rds
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(SpatialExperiment)
  devtools::load_all(quiet = TRUE)
})
source("data-raw/benchmarks/benchmark_helpers.R")

have_misty <- requireNamespace("mistyR", quietly = TRUE)
have_nnsvg <- requireNamespace("nnSVG", quietly = TRUE)
if (!have_misty && !have_nnsvg) {
  stop("Neither mistyR nor nnSVG installed. ",
       "Run data-raw/benchmarks/_install_comparators.R first.")
}
cat("mistyR available:", have_misty, " | nnSVG available:", have_nnsvg, "\n")

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
beta          <- -0.01
n_samples     <- 5
n_iterations  <- 10
n_gradient    <- 5
n_background  <- 45
fdr_threshold <- 0.05
base_seed     <- 9000

total_runs <- n_iterations
cat(sprintf("=== Comparison Benchmark ===\n"))
cat(sprintf("  beta = %.3f, N = %d, iterations = %d\n",
            beta, n_samples, n_iterations))
cat(sprintf("  %d gradient genes + %d background genes per run\n",
            n_gradient, n_background))

# ---------------------------------------------------------------------------
# Helper: run nnSVG per sample, return gene x sample log10(padj) matrix
# ---------------------------------------------------------------------------
run_nnsvg_persample <- function(spe) {
  if (!have_nnsvg) return(NULL)
  suppressPackageStartupMessages({
    library(nnSVG)
    library(SpatialExperiment)
  })

  samples <- unique(as.character(colData(spe)$sample_id))
  per_sample_padj <- list()

  for (samp in samples) {
    keep <- colData(spe)$sample_id == samp
    spe_s <- spe[, keep]
    # nnSVG requires log-normalized counts in 'logcounts' assay
    counts_mat <- as.matrix(assay(spe_s, "counts"))
    lib_sizes  <- pmax(colSums(counts_mat), 1)
    logc       <- log1p(t(t(counts_mat) / lib_sizes * 1e4))
    assay(spe_s, "logcounts") <- logc

    # Run nnSVG. assay_name = "logcounts" is the default; we set it
    # explicitly. nnSVG fits a nearest-neighbour Gaussian process per gene.
    res <- tryCatch(
      nnSVG::nnSVG(spe_s, assay_name = "logcounts", verbose = FALSE),
      error = function(e) {
        warning(sprintf("nnSVG failed for sample %s: %s", samp, e$message))
        NULL
      }
    )
    if (!is.null(res)) {
      rd <- rowData(res)
      per_sample_padj[[samp]] <- data.table(
        gene   = rownames(res),
        pval   = rd$pval,
        padj   = rd$padj
      )
    }
  }
  per_sample_padj
}

# ---------------------------------------------------------------------------
# Helper: aggregate nnSVG per-sample padj across samples (min-p pooling)
# ---------------------------------------------------------------------------
aggregate_nnsvg <- function(per_sample_padj, fdr = 0.05) {
  if (length(per_sample_padj) == 0) return(NULL)
  dt <- rbindlist(lapply(names(per_sample_padj), function(s) {
    x <- per_sample_padj[[s]]
    x[, sample_id := s]
    x
  }))
  agg <- dt[, .(
    min_padj  = min(padj, na.rm = TRUE),
    n_sig     = sum(padj < fdr, na.rm = TRUE),
    n_samples = .N
  ), by = gene]
  # A gene is called an SVG if it's significant in a majority of samples
  agg[, is_svg := n_sig >= ceiling(n_samples / 2)]
  agg
}

# ---------------------------------------------------------------------------
# Helper: run MISTy with a celltype-anchored para-view per sample,
# return variable importance aggregated across samples
# ---------------------------------------------------------------------------
run_misty_persample <- function(spe, query_celltype = "Tumor",
                                target_celltype = "T_cell") {
  if (!have_misty) return(NULL)
  suppressPackageStartupMessages({
    library(mistyR)
    library(SpatialExperiment)
  })

  samples <- unique(as.character(colData(spe)$sample_id))
  importance_list <- list()

  for (samp in samples) {
    keep <- colData(spe)$sample_id == samp
    spe_s <- spe[, keep]

    is_target <- colData(spe_s)$cell_type == target_celltype
    is_query  <- colData(spe_s)$cell_type == query_celltype
    if (sum(is_target) < 30 || sum(is_query) < 10) next

    target_expr <- t(as.matrix(assay(spe_s, "counts")[, is_target]))
    target_expr <- log1p(target_expr / pmax(rowSums(target_expr), 1) * 1e4)
    target_xy   <- as.data.frame(spatialCoords(spe_s)[is_target, , drop = FALSE])
    colnames(target_xy) <- c("x", "y")

    # Para-view: one-hot indicator of the query cell type across neighbours
    # within a kernel. Use Gaussian kernel with radius = 50 um.
    query_xy <- as.data.frame(spatialCoords(spe_s)[is_query, , drop = FALSE])
    colnames(query_xy) <- c("x", "y")
    # For each target cell: minimum distance to any query cell
    # (proxy for "query neighbour density"); MISTy's paraview uses kernel
    # sums, but here we use a simple distance-derived feature vector so that
    # the per-gene importance can be interpreted as "how much does proximity
    # to query predict this gene's expression".
    d <- sqrt(outer(target_xy$x, query_xy$x, "-")^2 +
              outer(target_xy$y, query_xy$y, "-")^2)
    min_d <- apply(d, 1, min)
    # Para-view table: one "feature" = exp(-min_d / 50) (Gaussian kernel)
    paraview <- data.frame(query_proximity = exp(-min_d / 50))

    misty_views <- tryCatch(
      {
        v <- mistyR::create_initial_view(as.data.frame(target_expr))
        v <- mistyR::add_views(v,
          mistyR::create_view("paraview_query", paraview, "query_proximity")
        )
        v
      },
      error = function(e) {
        warning(sprintf("MISTy view creation failed (%s): %s", samp, e$message))
        NULL
      }
    )
    if (is.null(misty_views)) next

    misty_out <- tryCatch(
      mistyR::run_misty(misty_views, results.folder = tempfile("misty_")),
      error = function(e) {
        warning(sprintf("MISTy run failed (%s): %s", samp, e$message))
        NULL
      }
    )
    if (is.null(misty_out)) next

    imp <- mistyR::collect_results(misty_out)$importances
    # Filter to the para-view and take importances of query_proximity
    imp_dt <- as.data.table(imp)
    if (!"view" %in% names(imp_dt) || !"Predictor" %in% names(imp_dt)) next
    imp_q <- imp_dt[view == "paraview_query" &
                    Predictor == "query_proximity",
                    .(gene = Target, importance = Importance,
                      sample_id = samp)]
    importance_list[[samp]] <- imp_q
  }
  rbindlist(importance_list, fill = TRUE)
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
all_results <- list()

for (iter in seq_len(n_iterations)) {
  seed <- base_seed + iter
  cat(sprintf("\n[%d/%d] iteration %d (seed=%d)\n",
              iter, n_iterations, iter, seed))

  spe <- generate_benchmark_data(
    n_samples      = n_samples,
    n_gradient_neg = n_gradient,
    n_gradient_pos = 0,
    n_background   = n_background,
    beta           = beta,
    seed           = seed
  )
  # Add cell_type column alias for MISTy (expects "cell_type")
  colData(spe)$cell_type <- colData(spe)$cell_type

  gradient_genes <- paste0("GRAD_NEG_", seq_len(n_gradient))
  bg_genes       <- paste0("BG_", sprintf("%02d", seq_len(n_background)))
  all_gene_names <- c(gradient_genes, bg_genes)

  # ---- RIPPLE ---------------------------------------------------------
  cat("  RIPPLE...\n")
  ripple_res <- tryCatch(
    run_ripple_quiet(spe),
    error = function(e) { warning("RIPPLE failed: ", e$message); NULL }
  )
  ripple_dt <- if (!is.null(ripple_res)) {
    r <- ripple_res[cell_type == "T_cell",
                    .(gene, median_coef, fisher_fdr, sign_consistency)]
    r[, is_planted   := gene %in% gradient_genes]
    r[, ripple_sig   := fisher_fdr < fdr_threshold]
    r
  } else NULL

  # ---- nnSVG -----------------------------------------------------------
  nnsvg_agg <- NULL
  if (have_nnsvg) {
    cat("  nnSVG (per-sample)...\n")
    per_s <- run_nnsvg_persample(spe)
    nnsvg_agg <- aggregate_nnsvg(per_s, fdr = fdr_threshold)
    if (!is.null(nnsvg_agg)) {
      nnsvg_agg[, is_planted := gene %in% gradient_genes]
    }
  }

  # ---- MISTy -----------------------------------------------------------
  misty_agg <- NULL
  if (have_misty) {
    cat("  MISTy (per-sample)...\n")
    misty_raw <- tryCatch(
      run_misty_persample(spe, query_celltype = "Tumor",
                          target_celltype = "T_cell"),
      error = function(e) { warning("MISTy failed: ", e$message); NULL }
    )
    if (!is.null(misty_raw) && nrow(misty_raw) > 0) {
      misty_agg <- misty_raw[, .(
        mean_importance = mean(importance, na.rm = TRUE),
        n_samples_misty = .N
      ), by = gene]
      misty_agg[, is_planted := gene %in% gradient_genes]
    }
  }

  # ---- Store per-iteration summary ------------------------------------
  iter_summary <- list(
    iter       = iter,
    seed       = seed,
    ripple     = ripple_dt,
    nnsvg      = nnsvg_agg,
    misty      = misty_agg
  )

  # Aggregate single-row metrics for cross-iteration summary
  metrics <- data.table(iter = iter, seed = seed)
  if (!is.null(ripple_dt)) {
    metrics[, ripple_tp := sum(ripple_dt$ripple_sig & ripple_dt$is_planted)]
    metrics[, ripple_fp := sum(ripple_dt$ripple_sig & !ripple_dt$is_planted)]
    metrics[, ripple_n_tested := nrow(ripple_dt)]
  }
  if (!is.null(nnsvg_agg)) {
    metrics[, nnsvg_tp := sum(nnsvg_agg$is_svg & nnsvg_agg$is_planted)]
    metrics[, nnsvg_fp := sum(nnsvg_agg$is_svg & !nnsvg_agg$is_planted)]
    metrics[, nnsvg_n_tested := nrow(nnsvg_agg)]
  }
  if (!is.null(misty_agg) && !is.null(ripple_dt)) {
    joined <- merge(
      misty_agg[, .(gene, mean_importance)],
      ripple_dt[, .(gene, median_coef, fisher_fdr)],
      by = "gene", all.x = TRUE
    )
    metrics[, misty_ripple_cor := if (nrow(joined) >= 5) {
      stats::cor(abs(joined$median_coef), joined$mean_importance,
                 use = "pairwise.complete.obs", method = "spearman")
    } else NA_real_]
  }
  iter_summary$metrics <- metrics

  all_results[[iter]] <- iter_summary
}

# ---------------------------------------------------------------------------
# Summarize
# ---------------------------------------------------------------------------
cat("\n=== Cross-iteration summary ===\n")
summary_dt <- rbindlist(lapply(all_results, `[[`, "metrics"), fill = TRUE)
print(summary_dt)

# Aggregate
agg <- data.table(
  n_iter              = nrow(summary_dt),
  ripple_mean_tp      = mean(summary_dt$ripple_tp, na.rm = TRUE),
  ripple_mean_fp      = mean(summary_dt$ripple_fp, na.rm = TRUE),
  nnsvg_mean_tp       = mean(summary_dt$nnsvg_tp, na.rm = TRUE),
  nnsvg_mean_fp       = mean(summary_dt$nnsvg_fp, na.rm = TRUE),
  misty_ripple_rho    = mean(summary_dt$misty_ripple_cor, na.rm = TRUE)
)
cat("\nAggregated:\n")
print(agg)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out_path <- "data-raw/benchmarks/results/bench_comparison_results.rds"
saveRDS(list(per_iteration = all_results, summary = summary_dt,
             aggregate = agg), file = out_path)
cat("\nSaved:", out_path, "\n")
