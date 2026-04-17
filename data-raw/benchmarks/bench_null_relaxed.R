# ============================================================================
# FDR null calibration benchmark — relaxed sign-consistency gate (0.75)
# ============================================================================
# Same design as bench_null.R but with sign_consistency = 0.75 instead of 1.0.
# Quantifies the FDR tradeoff when relaxing the requirement that ALL replicates
# must agree on coefficient direction.
#
# Run with:
#   Rscript data-raw/benchmarks/bench_null_relaxed.R
#
# Output:
#   data-raw/benchmarks/results/bench_null_relaxed_results.rds
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  devtools::load_all(quiet = TRUE)
})
source("data-raw/benchmarks/benchmark_helpers.R")

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
n_iterations <- 50
sample_sizes <- c(3, 5, 10)
n_background <- 50
base_seed <- 2026 # same seeds as bench_null.R for direct comparison
sign_consistency <- 0.75

cat("=== FDR Null Calibration — Relaxed Sign Consistency (0.75) ===\n")
cat(sprintf(
  "  %d iterations x %d sample sizes = %d runs\n",
  n_iterations, length(sample_sizes),
  n_iterations * length(sample_sizes)
))

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
all_results <- list()
counter <- 0
total_runs <- n_iterations * length(sample_sizes)

for (n_samp in sample_sizes) {
  for (iter in seq_len(n_iterations)) {
    counter <- counter + 1
    seed <- base_seed * 1000 + n_samp * 100 + iter

    if (counter %% 10 == 1 || counter == total_runs) {
      cat(sprintf(
        "[%d/%d] N=%d, iter=%d (seed=%d)\n",
        counter, total_runs, n_samp, iter, seed
      ))
    }

    spe <- generate_benchmark_data(
      n_samples      = n_samp,
      n_gradient_neg = 0,
      n_gradient_pos = 0,
      n_background   = n_background,
      seed           = seed
    )

    res <- tryCatch(
      run_ripple_quiet(spe, sign_consistency = sign_consistency),
      error = function(e) {
        warning(sprintf("N=%d iter=%d failed: %s", n_samp, iter, e$message))
        NULL
      }
    )

    if (!is.null(res)) {
      tcell_res <- res[cell_type == "T_cell"]
      n_tested <- nrow(tcell_res)
      n_sig_fdr <- sum(tcell_res$fisher_fdr < 0.05, na.rm = TRUE)
      n_sig_pval <- sum(tcell_res$fisher_pval < 0.05, na.rm = TRUE)

      all_results[[counter]] <- data.table(
        n_samples = n_samp,
        iteration = iter,
        seed = seed,
        sign_consistency = sign_consistency,
        n_genes_tested = n_tested,
        n_sig_fdr = n_sig_fdr,
        n_sig_pval = n_sig_pval,
        empirical_fdr = n_sig_fdr / max(n_tested, 1),
        empirical_fwer = as.integer(n_sig_fdr > 0)
      )
    }
  }
}

results_dt <- rbindlist(all_results)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat("\n=== Results (sign_consistency = 0.75) ===\n")
summary_dt <- results_dt[, .(
  n_runs = .N,
  mean_genes_tested = mean(n_genes_tested),
  mean_fdr = mean(empirical_fdr),
  sd_fdr = sd(empirical_fdr),
  max_fdr = max(empirical_fdr),
  mean_fwer = mean(empirical_fwer),
  total_sig = sum(n_sig_fdr),
  total_tested = sum(n_genes_tested)
), by = n_samples]

summary_dt[, pooled_fdr := total_sig / total_tested]
print(summary_dt)

cat("\nTarget: empirical FDR ≤ 0.075 (1.5x nominal 0.05)\n")
cat("Pass:", all(summary_dt$pooled_fdr <= 0.075), "\n")

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out_path <- "data-raw/benchmarks/results/bench_null_relaxed_results.rds"
saveRDS(list(per_run = results_dt, summary = summary_dt), file = out_path)
cat("\nSaved:", out_path, "\n")
