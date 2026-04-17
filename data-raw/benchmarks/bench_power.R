# ============================================================================
# Power benchmark
# ============================================================================
# Varies gradient effect size (beta) and sample count to compute power curves.
# For each combination, generates datasets with 5 gradient genes + 45 background
# genes. Power = fraction of gradient genes recovered at FDR < 0.05.
#
# Design: 4 effect sizes × 3 sample sizes × 30 iterations.
#
# Run with:
#   Rscript data-raw/benchmarks/bench_power.R
#
# Output:
#   data-raw/benchmarks/results/bench_power_results.rds
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  devtools::load_all(quiet = TRUE)
})
source("data-raw/benchmarks/benchmark_helpers.R")

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
betas <- c(-0.002, -0.005, -0.01, -0.02)
sample_sizes <- c(3, 5, 10)
n_iterations <- 30
n_gradient <- 5
n_background <- 45
base_seed <- 7777

cat("=== Power Benchmark ===\n")
total_runs <- length(betas) * length(sample_sizes) * n_iterations
cat(sprintf(
  "  %d betas x %d sample sizes x %d iterations = %d runs\n",
  length(betas), length(sample_sizes), n_iterations, total_runs
))

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
all_results <- list()
counter <- 0

for (b in betas) {
  for (n_samp in sample_sizes) {
    for (iter in seq_len(n_iterations)) {
      counter <- counter + 1
      seed <- base_seed * 1000 + abs(b * 1e5) + n_samp * 100 + iter

      if (counter %% 20 == 1 || counter == total_runs) {
        cat(sprintf(
          "[%d/%d] beta=%.3f, N=%d, iter=%d\n",
          counter, total_runs, b, n_samp, iter
        ))
      }

      spe <- generate_benchmark_data(
        n_samples      = n_samp,
        n_gradient_neg = n_gradient,
        n_gradient_pos = 0,
        n_background   = n_background,
        beta           = b,
        seed           = seed
      )

      res <- tryCatch(
        run_ripple_quiet(spe),
        error = function(e) {
          warning(sprintf(
            "beta=%.3f N=%d iter=%d failed: %s",
            b, n_samp, iter, e$message
          ))
          NULL
        }
      )

      if (!is.null(res)) {
        tcell_res <- res[cell_type == "T_cell"]

        gradient_genes <- paste0("GRAD_NEG_", seq_len(n_gradient))
        bg_genes <- paste0("BG_", sprintf("%02d", seq_len(n_background)))

        grad_tested <- tcell_res[gene %in% gradient_genes]
        bg_tested <- tcell_res[gene %in% bg_genes]

        tp <- sum(grad_tested$fisher_fdr < 0.05, na.rm = TRUE)
        fn <- nrow(grad_tested) - tp
        fp <- sum(bg_tested$fisher_fdr < 0.05, na.rm = TRUE)
        tn <- nrow(bg_tested) - fp

        all_results[[counter]] <- data.table(
          beta = b,
          n_samples = n_samp,
          iteration = iter,
          seed = seed,
          n_grad_tested = nrow(grad_tested),
          n_bg_tested = nrow(bg_tested),
          tp = tp, fn = fn, fp = fp, tn = tn,
          sensitivity = tp / max(tp + fn, 1),
          specificity = tn / max(tn + fp, 1),
          fdr_obs = fp / max(tp + fp, 1)
        )
      }
    }
  }
}

results_dt <- rbindlist(all_results)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat("\n=== Power Summary ===\n")
summary_dt <- results_dt[, .(
  n_runs       = .N,
  mean_power   = mean(sensitivity),
  sd_power     = sd(sensitivity),
  mean_fdr     = mean(fdr_obs),
  mean_spec    = mean(specificity)
), by = .(beta, n_samples)]

print(summary_dt[order(beta, n_samples)])

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out_path <- "data-raw/benchmarks/results/bench_power_results.rds"
saveRDS(list(per_run = results_dt, summary = summary_dt), file = out_path)
cat("\nSaved:", out_path, "\n")
