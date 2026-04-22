# ============================================================================
# Runtime and scalability benchmark
# ============================================================================
# Measures wall-clock time (via bench::mark) and peak R memory (via gc()
# max-used Mb column) for run_ripple across three dataset sizes.
#
# Dataset sizes reflect realistic spatial transcriptomics workloads:
#   Small:  ~3k   cells x 100 genes x 3 samples
#   Medium: ~50k  cells x 300 genes x 5 samples
#   Large:  ~250k cells x 500 genes x 5 samples
#
# bench::mark is used for timing with a single iteration per configuration
# (results are per-config; loop provides cross-size comparison). An outer
# rep loop (n_reps = 3) gives cross-run variability for the summary.
#
# Run with:
#   Rscript data-raw/benchmarks/bench_runtime.R
#
# Output:
#   data-raw/benchmarks/results/bench_runtime_results.rds
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(bench)
  devtools::load_all(quiet = TRUE)
})
source("data-raw/benchmarks/benchmark_helpers.R")

# ---------------------------------------------------------------------------
# Size presets
# ---------------------------------------------------------------------------
size_configs <- list(
  small = list(
    cells = list(Tumor = 150, T_cell = 600, Fibroblast = 250),
    n_genes = 100, n_samples = 3
  ),
  medium = list(
    cells = list(Tumor = 1500, T_cell = 6000, Fibroblast = 2500),
    n_genes = 300, n_samples = 5
  ),
  large = list(
    cells = list(Tumor = 7500, T_cell = 30000, Fibroblast = 12500),
    n_genes = 500, n_samples = 5
  )
)

n_reps <- 3
base_seed <- 12345

cat("=== Runtime Benchmark (bench::mark + gc Mb) ===\n")
for (s in names(size_configs)) {
  cfg <- size_configs[[s]]
  total_cells <- sum(unlist(cfg$cells)) * cfg$n_samples
  cat(sprintf(
    "  %s:  %s cells total x %d genes x %d samples\n",
    s, format(total_cells, big.mark = ","),
    cfg$n_genes, cfg$n_samples
  ))
}
cat(sprintf("  %d reps per size = %d total runs\n\n",
            n_reps, n_reps * length(size_configs)))

# ---------------------------------------------------------------------------
# gc_peak_mb: sum of max-used Mb across Ncells and Vcells
# ---------------------------------------------------------------------------
gc_peak_mb <- function() {
  g <- gc()
  # gc() returns 6 columns: used, (Mb), gc trigger, (Mb), max used, (Mb)
  # We want the last column (max used in Mb) summed across Ncells + Vcells rows.
  sum(g[, ncol(g)])
}

# ---------------------------------------------------------------------------
# measure_one_run: uses bench::mark for time (1 iteration, just to get its
# precise timing + memory allocation tracking) and gc() for peak resident.
# ---------------------------------------------------------------------------
measure_one_run <- function(cfg, seed) {
  gc(reset = TRUE, full = TRUE)

  # Data generation time measured separately (not the main benchmark target)
  t_gen <- system.time({
    spe <- generate_benchmark_data(
      n_samples        = cfg$n_samples,
      n_gradient_neg   = 0,
      n_gradient_pos   = 0,
      n_background     = cfg$n_genes,
      cells_per_sample = cfg$cells,
      seed             = seed
    )
  })[3]

  # Reset memory tracking before the RIPPLE run so peak reflects RIPPLE only
  gc(reset = TRUE, full = TRUE)

  # bench::mark with a single iteration. memory = FALSE because run_ripple
  # internally parallelises (bench's memory profiler cannot track that);
  # peak memory is captured via gc() immediately after the run instead.
  bm <- bench::mark(
    run_ripple_quiet(spe),
    iterations = 1,
    check      = FALSE,
    memory     = FALSE,
    time_unit  = "s",
    filter_gc  = FALSE
  )

  t_ripple <- as.numeric(bm$median)
  peak_mb  <- gc_peak_mb()

  total_cells <- sum(unlist(cfg$cells)) * cfg$n_samples

  list(
    n_cells    = total_cells,
    n_genes    = cfg$n_genes,
    n_samples  = cfg$n_samples,
    t_generate = as.numeric(t_gen),
    t_ripple   = t_ripple,
    peak_mb    = peak_mb
  )
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
all_results <- list()
counter <- 0
total_runs <- n_reps * length(size_configs)

for (size in names(size_configs)) {
  cfg <- size_configs[[size]]
  for (rep in seq_len(n_reps)) {
    counter <- counter + 1
    seed <- base_seed + match(size, names(size_configs)) * 100 + rep

    cat(sprintf("[%d/%d] size=%s rep=%d... ", counter, total_runs, size, rep))

    result <- tryCatch(
      measure_one_run(cfg, seed),
      error = function(e) {
        cat("ERROR:", e$message, "\n")
        NULL
      }
    )

    if (!is.null(result)) {
      cat(sprintf("RIPPLE %.1fs, peak %.0f MB\n",
                  result$t_ripple, result$peak_mb))
      result$size <- size
      result$rep  <- rep
      result$seed <- seed
      all_results[[counter]] <- as.data.table(result)
    }

    rm(result)
    gc(reset = TRUE, full = TRUE)
  }
}

results_dt <- rbindlist(all_results)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat("\n=== Summary ===\n")
summary_dt <- results_dt[, .(
  n_reps        = .N,
  n_cells       = mean(n_cells),
  n_genes       = mean(n_genes),
  n_samples     = mean(n_samples),
  mean_t_ripple = mean(t_ripple),
  sd_t_ripple   = sd(t_ripple),
  mean_peak_mb  = mean(peak_mb)
), by = size]
print(summary_dt)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out_path <- "data-raw/benchmarks/results/bench_runtime_results.rds"
saveRDS(list(per_run = results_dt, summary = summary_dt), file = out_path)
cat("\nSaved:", out_path, "\n")
