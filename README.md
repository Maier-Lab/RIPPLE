# RIPPLE

**Replicate-Aware Inference of Paracrine Profiles via Likelihood Estimation**

RIPPLE is an R package that detects distance-dependent gene expression gradients from a chosen query cell type in spatial transcriptomics data.

---

## The question

> *"Which genes in cell type B change expression as a function of physical distance from cell type A, reproducibly across biological replicates?"*

Given a spatial transcriptomics dataset with cell type annotations and biological replicates, RIPPLE fits a per-sample Poisson GLM for each gene in each target cell type, using Euclidean distance to the nearest query cell as the predictor and `log(total_counts)` as an offset. Per-sample coefficients are combined across replicates via Fisher's combined p-value with a sign-consistency gate. The result is a ranked list of genes with signed, interpretable gradient coefficients, calibrated FDR, and per-sample reproducibility.

**Supported platforms:** Xenium, MERFISH, CosMx, CODEX -- any imaging-based platform with single-cell resolved coordinates and raw integer counts.

For a standalone SLURM-driven script version, see [HyMy-distance-correlation-analysis](https://github.com/CMangana/HyMy-distance-correlation-analysis).

---

## Validation

RIPPLE has been benchmarked on synthetic data. Scripts live in `data-raw/benchmarks/` and rendered results are in `vignettes/benchmarks.Rmd`.

- **FDR calibration.** 150 null simulations (50 genes each, no planted gradient) across N = 3, 5, 10 replicates yielded 1 false positive out of 7,500 tests (empirical false discovery proportion < 0.05%). The framework is conservative by design; the sign-consistency gate is the primary driver.
- **Power.** 100% recovery of planted gradient genes for |beta| >= 0.005 at all tested sample sizes (N = 3, 5, 10). At the weakest effect (|beta| = 0.002), power scales from 55% (N = 3) to 85% (N = 10).
- **Sign-consistency tradeoff.** Relaxing `sign_consistency` from 1.0 to 0.75 does not meaningfully increase the false discovery proportion (still < 0.1%) and is a reasonable choice for studies with N >= 6 where one discordant sample is plausible.

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("CMangana/RIPPLE")
```

Installing from GitHub pulls in all required dependencies automatically. Optional functionality requires additional packages:

| Optional feature | Extra packages |
|------------------|----------------|
| Multi-panel figures and vignettes | `knitr`, `rmarkdown` |
| Pathway enrichment (`run_ripple_fgsea`) | `fgsea`, `msigdbr` (Bioconductor) |
| Ligand-receptor integration (`run_ripple_lr`) | `nichenetr` (GitHub: `saeyslab/nichenetr`) |
| Input from `SingleCellExperiment`/`SpatialExperiment` | `SingleCellExperiment`, `SpatialExperiment`, `SummarizedExperiment`, `S4Vectors` (Bioconductor) |
| GPU permutation testing (Stage 2) | Python 3, PyTorch (CUDA), NumPy < 2, SciPy |

Install Bioconductor extras with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("fgsea", "msigdbr", "SpatialExperiment", "SingleCellExperiment"))
```

---

## Quick start

The package ships with `ripple_mock_data`, a small synthetic `SpatialExperiment` (50 genes, 600 cells, 3 samples) containing a planted distance-dependent gradient. You can run the full pipeline on it without any external data:

```r
library(ripple)
data(ripple_mock_data)

results <- run_ripple(
  input           = ripple_mock_data,
  query_celltype  = "Tumor",
  celltype_column = "cell_type",
  sample_column   = "sample_id",
  output_dir      = tempfile("ripple_demo_")
)

# The planted INDUCED_* and REPRESSED_* genes should top the list
head(results[order(fisher_fdr)], 12)
```

On real data, `input` can be a `Seurat`, `SingleCellExperiment`, or `SpatialExperiment` object (in memory or as an `.rds` path). For raw matrices or CSV inputs, use `make_ripple_input()` or `read_ripple_csv()` to build a canonical object first. See `?run_ripple` for the full argument list.

---

## Vignettes

Three long-form vignettes ship with the package. They cover, respectively, a real-data walkthrough, positioning vs. other spatial tools, and benchmark details.

| Vignette | Description |
|----------|-------------|
| `cosmx_nsclc_walkthrough` | End-to-end example on the public CosMx NSCLC dataset (He et al., 2022). Loads cached results from `inst/extdata/`, so the vignette renders without re-running the full pipeline. |
| `methods_positioning` | Landscape table comparing RIPPLE with Hotspot, nnSVG, MISTy, COMMOT, and BANKSY. Clarifies the distinction between co-localization tests and continuous gradient detection. |
| `benchmarks` | FDR calibration and power curves from the synthetic benchmark suite. |

Browse locally:

```r
browseVignettes("ripple")

# Or build the full pkgdown site (shipped _pkgdown.yml):
pkgdown::build_site()
```

---

## Pipeline overview

Each stage is optional except Stage 1.

| Stage | Function(s) | Purpose |
|-------|-------------|---------|
| 1. Distance correlation | `run_ripple()` | Per-sample Poisson GLM + Fisher combined p-value with sign-consistency gate |
| 2. Permutation validation | `run_permutation_tests()` (R) or `inst/python/run_permutation_gpu.py` (GPU) | Validates query specificity via label permutation |
| 3. Merge and summarize | `merge_ripple_results()`, `compute_fisher_pval()` | Combines per-celltype results, recomputes Fisher p-values |
| 4. Confounder control | `run_ripple_confounder()` | Bivariate GLM isolating query-specific from shared-niche effects |
| 5. Atlas figures | `run_ripple_atlas()`, `run_ripple_fgsea()`, `plot_gradient_volcano()`, `plot_decay_curve()` | Multi-panel figures, pathway enrichment, contamination flagging |
| 6. Ligand-receptor integration | `run_ripple_lr()`, `classify_lr_artifacts()` | Matches gradient genes to L-R pairs via NicheNet |

Diagnostics:

- `check_spatial_autocorrelation()` computes Moran's I on Poisson residuals for selected genes, flagging cases where the independence assumption may be violated.
- `check_data()` and `load_metadata_only()` give fast metadata access without loading the full expression matrix.

---

## Core statistical model

For each gene $g$ in target cell type $t$, per biological replicate $s$:

```r
glm(counts ~ distance_to_query + offset(log(total_counts)), family = poisson)
```

- `distance_to_query`: Euclidean distance (um) to the nearest query cell (default k = 1)
- `offset(log(total_counts))`: cell-size correction -- converts raw counts to rates, controls for ambient RNA and segmentation differences that co-vary with cell size
- **Coefficient (beta)**: log-rate change per um. Negative = expression increases near query cells (induced). Positive = expression decreases (repressed).

Per-sample coefficients are combined via Fisher's combined p-value. The sign-consistency gate requires all replicates to agree on the direction of the effect (`sign_consistency = 1.0` by default; relax to 0.75 for N >= 6). `fisher_fdr` is the primary significance metric for all downstream analyses.

---

## Input expectations

| Component | Expected |
|-----------|----------|
| Counts | Raw integer counts in `assays(spe, "counts")` for a `SpatialExperiment`, or `obj[["RNA"]]$counts` for Seurat. The Poisson model handles normalization internally via the offset -- pre-normalized data will produce incorrect results. |
| Spatial coordinates | `spatialCoords()` for `SpatialExperiment`, or X/Y columns in cell metadata for Seurat/SCE. Common column names (`x_centroid`/`y_centroid`, `spatial_x`/`spatial_y`, `x`/`y`) are auto-detected; override via `x_column` / `y_column`. |
| Cell types | Metadata column named by `celltype_column`, containing the query population. |
| Replicate ID | Metadata column named by `sample_column` (default `"sample_id"`). Minimum 3 biological replicates; 4+ recommended. |
| Condition (optional) | Metadata column named by `condition_column`, with the target value in `condition_value`. |

---

## Output structure

`run_ripple()` writes to `{output_dir}/{analysis_name}/` (with a `_k{n}` suffix if `k_neighbors > 1`):

```
{output_dir}/{analysis_name}/
  per_celltype/
    {CellType}/
      meta_analysis_results.csv   # Per-gene meta-analysis (coefficient, FDR)
      coef_per_sample.csv         # Per-sample GLM coefficients
      gradient_scores.csv         # Gradient magnitude scores
      decay_classification.csv    # Decay pattern classification
      gradient_volcano.pdf        # Per-celltype volcano
      coefficient_strips.pdf      # Per-sample coefficient strips
      forest_plots/               # Per-gene forest plots
  summary/
    all_genes_results.csv         # Merged across all cell types
    top_gradient_genes.csv        # Top 50 per cell type by FDR
    decay_pattern_summary.csv     # Decay pattern counts
  qc/
    distance_distribution.pdf     # QC: distance-to-query distribution
```

`run_ripple_atlas()` adds a `plots/` (or `atlas/`) subdirectory with multi-panel figures and fGSEA output.

---

## Recommended QC workflow

After running `run_ripple()` and `merge_ripple_results()`, always check for contamination before interpreting individual genes:

```r
# 1. Classify gene specificity — do this FIRST
specificity <- classify_gene_specificity(all_results, fdr_threshold = 0.05)
table(specificity$specificity_class)
#   specific  moderate  ubiquitous  contamination
#       412       187          45             23

# 2. Exclude contamination genes (significant in too many cell types)
contam_genes <- specificity[specificity_class == "contamination"]$gene
clean_results <- all_results[!gene %in% contam_genes]

# 3. THEN look at individual genes
plot_gradient_volcano(clean_results[cell_type == "CD8_T"], query_label = "Tumor")
```

Genes flagged as "contamination" are significant in many cell types simultaneously, which usually indicates segmentation artifacts (query cell transcripts leaking into neighboring cells) rather than genuine biology.

For k-selection, use `plot_k_diagnostics()` before running the full pipeline to choose an appropriate `k_neighbors` value:

```r
plot_k_diagnostics(my_spe, query_celltype = "Tumor", celltype_column = "cell_type")
```

---

## Citation

A bibentry is shipped in `inst/CITATION` and can be retrieved with `citation("ripple")`. The placeholder entry will be replaced once the accompanying manuscript is posted on bioRxiv.

> *Citation forthcoming (bioRxiv preprint in preparation).*

---

## License

MIT. See `LICENSE`.
