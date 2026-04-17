# ============================================================================
# Runtime and scalability benchmark
# ============================================================================
# Measures wall-clock time and peak memory for run_ripple across dataset sizes.
#
# Dataset sizes reflect realistic spatial transcriptomics workloads:
#   Small:  ~3k  cells  x 100 genes x 3 samples
#   Medium: ~50k cells  x 300 genes x 5 samples
#   Large:  ~250k cells x 500 genes x 5 samples
#
# Each size is run 3 times to get a stable mean / SD.
#
# Run with:
#   Rscript data-raw/benchmarks/bench_runtime.R
#
# Output:
#   data-raw/benchmarks/results/bench_runtime_results.rds
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
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

cat("=== Runtime Benchmark ===\n")
for (s in names(size_configs)) {
  cfg <- size_configs[[s]]
  total_cells <- sum(unlist(cfg$cells)) * cfg$n_samples
  cat(sprintf(
    "  %s:  %s cells total x %d genes x %d samples\n",
    s, format(total_cells, big.mark = ","),
    cfg$n_genes, cfg$n_samples
  ))
}
cat(sprintf(
  "  %d reps per size = %d total runs\n\n",
  n_reps, n_reps * length(size_configs)
))

# ---------------------------------------------------------------------------
# Benchmark function
# ---------------------------------------------------------------------------
measure_one_run <- function(cfg, seed) {
  # Reset memory tracking
  gc(reset = TRUE, full = TRUE)
  mem_before <- sum(gc()[, "used"]) # MB (approx, ncells column is in KB-ish units)

  t_gen_start <- proc.time()
  spe <- generate_benchmark_data(
    n_samples        = cfg$n_samples,
    n_gradient_neg   = 0,
    n_gradient_pos   = 0,
    n_background     = cfg$n_genes,
    cells_per_sample = cfg$cells,
    seed             = seed
  )
  t_gen <- (proc.time() - t_gen_start)[3]

  t_run_start <- proc.time()
  res <- run_ripple_quiet(spe)
  t_run <- (proc.time() - t_run_start)[3]

  # Peak memory from gc() max-used tracking
  gc_info <- gc()
  peak_mb <- sum(gc_info[, "max used"] * c(56, 8) / 1024 / 1024)

  # n_cells is the total cells across all samples
  total_cells <- sum(unlist(cfg$cells)) * cfg$n_samples

  list(
    n_cells     = total_cells,
    n_genes     = cfg$n_genes,
    n_samples   = cfg$n_samples,
    t_generate  = as.numeric(t_gen),
    t_ripple    = as.numeric(t_run),
    t_total     = as.numeric(t_gen + t_run),
    peak_mb     = peak_mb,
    n_results   = nrow(res)
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
      cat(sprintf(
        "%.1fs (RIPPLE: %.1fs), %.0f MB peak\n",
        result$t_total, result$t_ripple, result$peak_mb
      ))
      result$size <- size
      result$rep <- rep
      result$seed <- seed
      all_results[[counter]] <- as.data.table(result)
    }

    # Force garbage collection between runs
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
  n_reps = .N,
  n_cells = mean(n_cells),
  n_genes = mean(n_genes),
  n_samples = mean(n_samples),
  mean_t_ripple = mean(t_ripple),
  sd_t_ripple = sd(t_ripple),
  mean_peak_mb = mean(peak_mb)
), by = size]
print(summary_dt)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out_path <- "data-raw/benchmarks/results/bench_runtime_results.rds"
saveRDS(list(per_run = results_dt, summary = summary_dt), file = out_path)
cat("\nSaved:", out_path, "\n")
