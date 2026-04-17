# ripple 0.1.0

Initial public release of the RIPPLE package.

## Core method

* Per-replicate Poisson GLM with cell-size offset for distance-conditioned
  gene expression (`fit_poisson()`, `fit_poisson_controlled()`).
* Cross-replicate inference via Fisher's combined p-value with sign-consistency
  gating (`run_meta_analysis()`, `compute_fisher_pval()`).
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
