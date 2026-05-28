# ripple 0.2.0 (development)

## Breaking changes — hot-path simplification

* `run_meta_analysis()` (REML inverse-variance random-effects meta-analysis
  via `meta::metagen`) has been **removed**. Profiling showed it accounted
  for ~60% of `run_ripple()` wall-clock while populating columns that
  duplicated the Fisher path.
* `gradient_score` is now defined as `median_coef` (median of per-sample
  Poisson GLM coefficients), not as the REML inverse-variance combined
  coefficient. Per-gene shifts vs the previous value are small for typical
  N = 3–10 replicates.
* The following columns are no longer written by `run_ripple()`:
  `combined_coef`, `combined_se`, `pval`, `i2`, `fdr`. Downstream
  functions (`run_ripple_atlas()`, `run_ripple_lr()`, `merge_ripple_results()`,
  `run_ripple_confounder()`) retain legacy fallbacks so older CSV results
  on disk still load.
* `meta` is no longer an `Imports` dependency.
* `create_forest_plot()` summary row is now drawn at the median of
  per-sample coefficients (no horizontal CI bar). The subtitle reports only
  the sign-consistency count.

## Speedup

* Medium synthetic benchmark (5 samples × 300 genes × 2 target types):
  32.5 s → 14.5 s (~2.2× faster end-to-end).
* Real-data runs on full imaging panels with many target cell types should
  see proportional gains. Parallelisation across target cell types via
  SLURM array templates in `inst/slurm/` remains the largest additional
  lever.

---

# ripple 0.1.0

Initial public release of the RIPPLE package.

## Core method

* Per-replicate Poisson GLM with cell-size offset for distance-conditioned
  gene expression (`fit_poisson()`, `fit_poisson_controlled()`).
* Cross-replicate inference via Fisher's combined p-value with sign-consistency
  gating (`compute_fisher_pval()`).
* Two-tier expression filtering (strict for regular genes, lenient for
  user-supplied priority genes) to rescue sparse but biologically important
  transcripts.
* Optional confounder control via bivariate GLM with a second cell type
  (`run_ripple_confounder()`), with classification of genes as
  query-specific / enhanced / niche-driven / underpowered.

## Pipeline

* `run_ripple()` is the main entry point; `merge_ripple_results()`
  aggregates per-celltype output.
* `run_ripple_atlas()` produces publication-style figures
  (volcano, decay curves, dotplot, heatmap, fGSEA panels, contamination
  flagging).
* `run_ripple_fgsea()` performs reproducible pathway enrichment with a
  user-supplied seed; `run_ripple_lr()` integrates results with
  ligand–receptor databases via NicheNet.
* CPU permutation via `run_permutation_tests()`; GPU permutation script
  shipped under `inst/python/run_permutation_gpu.py`.

## Inputs

* Accepts in-memory or `.rds`-stored `Seurat`, `SingleCellExperiment`, and
  `SpatialExperiment` objects via a unified input adapter.
* `make_ripple_input()` builds a canonical object from raw counts, metadata,
  and coordinates; `read_ripple_csv()` loads from a directory of CSVs.
* All entry points perform input validation with informative error messages
  before any compute begins.

## Data

* `ripple_mock_data` — synthetic 50-gene × 600-cell × 3-sample dataset with
  a planted distance-dependent gradient in T cells; ships with the package
  for examples, tests, and tutorials.

## Reproducibility

* fGSEA results are deterministic given a user seed (`fgsea_seed`).
* Sign-consistency gating handles zero-coefficient samples consistently.

## Documentation

* Comprehensive README covering pipeline stages, statistical model,
  configuration, and troubleshooting.
* Roxygen2 documentation for all 37 exported functions.

## Packaging

* Installable as a standalone R package; legacy standalone scripts have
  been removed from the main code path (still available in
  [HyMy-distance-correlation-analysis](https://github.com/CMangana/HyMy-distance-correlation-analysis)
  for SLURM-driven workflows).
