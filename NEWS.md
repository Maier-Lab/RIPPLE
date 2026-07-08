# ripple 0.1.0

Initial release of the RIPPLE package, accompanying the preprint. Version
numbering starts counting from the first public release; all pre-release
development is collected here under 0.1.0.

## Core method

* Per-replicate Poisson GLM with cell-size offset for distance-conditioned
  gene expression (`fit_poisson()`, `fit_poisson_controlled()`).
* Cross-replicate inference via Fisher's combined p-value with sign-consistency
  gating (`compute_fisher_pval()`).
* `gradient_score` is defined as `median_coef`, the median of the per-sample
  Poisson GLM coefficients (equal weight per replicate).
* Two-tier expression filtering (strict for regular genes, lenient for
  a curated set of priority genes) to rescue sparse but biologically important
  transcripts. The built-in priority list (chemokines, cytokines, interleukins,
  interferons, and receptors) is curated from MGI, NCBI Gene, and Zlotnik &
  Yoshie (2012); a species-matched human list is derived automatically.
* Optional confounder control via bivariate GLM with a second cell type
  (`run_ripple_confounder()`), with classification of genes as
  query-specific / enhanced / niche-driven / underpowered.

## Pipeline

* `run_ripple()` is the main entry point; a single call already combines
  across samples (Fisher) and across cell types, and writes
  `summary/all_genes_results.csv`. `merge_ripple_results()` stitches together
  cell types that were run as separate jobs.
* `run_ripple_atlas()` produces publication-style figures
  (volcano, decay curves, dotplot, heatmap, fGSEA panels, contamination
  flagging).
* `run_ripple_fgsea()` performs reproducible pathway enrichment with a
  user-supplied seed; `run_ripple_lr()` integrates results with
  ligand-receptor databases via NicheNet.
* CPU permutation via `run_permutation_tests()`; GPU permutation script
  shipped under `inst/python/run_permutation_gpu.py`.

## Inputs

* Accepts in-memory or `.rds`-stored `Seurat`, `SingleCellExperiment`, and
  `SpatialExperiment` objects via a unified input adapter.
* `make_ripple_input()` builds a canonical object from raw counts, metadata,
  and coordinates; `read_ripple_csv()` loads from a directory of CSVs.
* All entry points perform input validation with informative error messages
  before any compute begins.
* Configuration is resolved as explicit argument, then `options(ripple.*)`,
  then environment variable, then a built-in default, so SLURM/env-driven
  runs honour the documented options.

## Gene specificity

* `classify_gene_specificity()` labels genes `specific` (1 cell type),
  `moderate` (2 up to `broad_threshold` - 1), or `broad`
  (>= `broad_threshold`). `broad_threshold` is the single boundary that
  defines the broad class; raising it flags fewer genes.

## Decay curves

* `bin_decay_data()` gains a `sample_ids` argument. When supplied, it returns
  a per-(sample, bin) table with a `sample_id` column instead of pooling cells
  across samples, making the recommended per-sample workflow a one-liner. It
  also gains `min_cells_per_bin` (default `10L`); bins with fewer cells are
  dropped to avoid unstable proportions.
* `plot_gradient_curve()` auto-detects per-sample mode when the input has a
  `sample_id` column, and gains `min_cells_per_bin` (default `10L`) and
  `min_samples_per_bin` (default `2L`) filters. Pooled mode remains available.
  Its docstring flags that pooled mode overstates precision when replicates
  disagree, and recommends per-sample mode for manuscript figures.
* `plot_k_diagnostics()` computes distances within each sample (never pooled
  across samples), so overlapping per-sample coordinate frames cannot produce
  meaningless cross-sample distances.

## Diagnostics and warnings

* Data-quality caveats are raised via `warning()` (not verbose-only messages),
  so they surface in batch/SLURM runs: distance-cap saturation, cell types
  skipped for too few cells/samples/genes, only two valid samples, high
  collinearity in the confounder model, empty results, and running on a
  subset of target cell types (which weakens the cross-cell-type
  contamination check).

## Performance

* Medium synthetic benchmark (5 samples x 300 genes x 2 target types):
  32.5 s -> 14.5 s (~2.2x faster end-to-end). Real-data runs on full imaging
  panels with many target cell types should see proportional gains;
  parallelising across target cell types via the `inst/slurm/` array templates
  remains the largest additional lever.

## Data

* `ripple_mock_data`, a synthetic 50-gene x 600-cell x 3-sample dataset with
  a planted distance-dependent gradient in T cells; ships with the package
  for examples, tests, and tutorials.

## Reproducibility

* fGSEA results are deterministic given a user seed (`fgsea_seed`).
* Sign-consistency gating handles zero-coefficient samples consistently.
* Permutation testing compares the observed statistic against a null built
  from the same statistic (the median of per-sample coefficients).

## Documentation

* Comprehensive README covering pipeline stages, statistical model,
  configuration, and troubleshooting.
* Roxygen2 documentation for all exported functions; four vignettes
  (getting started, CosMx NSCLC walkthrough, benchmarks, methods positioning).

## Packaging

* Installable as a standalone R package; legacy standalone scripts have
  been removed from the main code path
* Optional dependencies (Bioconductor input classes, fgsea/msigdbr, spdep,
  nichenetr, pheatmap) are guarded at runtime with actionable install
  messages.
