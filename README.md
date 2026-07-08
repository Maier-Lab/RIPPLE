# RIPPLE


<p align="center">
    <picture align="center">
    <img width="200" alt="ripple logo" src="https://github.com/user-attachments/assets/5aaab8a9-bb66-4117-b378-3f33dc13ea3c" />
    </picture>
</p>
<p align="center">
  <b>Replicate-Aware Inference of Paracrine Profiles via Likelihood Estimation</b>
</p>

<p align="center">
  <a href="https://github.com/Maier-Lab/RIPPLE/actions/workflows/R-CMD-check.yaml">
    <img alt="R-CMD-check status" src="https://github.com/Maier-Lab/RIPPLE/actions/workflows/R-CMD-check.yaml/badge.svg">
  </a>
</p>

RIPPLE is an R package that detects distance-dependent gene expression gradients from a chosen query cell type in spatial transcriptomics data.

---

## The question

> *"Which genes in cell type B change expression as a function of physical distance from cell type A, reproducibly across biological replicates?"*

For each gene in each target cell type, RIPPLE fits a per-sample Poisson GLM with Euclidean distance to the nearest query cell as the predictor and `log(total_counts)` as an offset. Per-sample coefficients are then combined across biological replicates via Fisher's combined p-value with a sign-consistency gate. The output is a ranked list of genes with signed gradient scores, calibrated FDR, and per-sample reproducibility.

**Supported platforms:** Xenium, CosMx, MERFISH, etc. Suited to any imaging-based platform with single-cell resolved coordinates and integer counts. Not designed for spot-resolution platforms (e.g. Visium without deconvolution) where one spot mixes multiple cell types.

---

## Installation

To build the vignettes at install time (they are skipped by default), install
the packages above first, then:

```r
devtools::install_github("Maier-Lab/RIPPLE", build_vignettes = TRUE)
```

It is recommended you install this as well, because `SpatialExperiment` is needed even for the quick start below. Install the Bioconductor
packages first so the object loads and the vignettes build:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "SpatialExperiment", "SingleCellExperiment", "SummarizedExperiment",
  "S4Vectors", "fgsea", "msigdbr"
))
install.packages(c("R.utils", "knitr", "rmarkdown"))
```
Then:
```r
# install.packages("devtools")
devtools::install_github("Maier-Lab/RIPPLE")
```

Optional functionality requires additional packages:

| Optional feature | Extra packages |
|------------------|----------------|
| Bundled data + quick start (`ripple_mock_data` is a `SpatialExperiment`) | `SpatialExperiment`, `SummarizedExperiment`, `S4Vectors` (Bioconductor) |
| QC dashboard (`ripple_plot_qc`, reads `cell_distances.csv.gz`) | `R.utils` |
| Multi-panel figures and vignettes | `knitr`, `rmarkdown` |
| Pathway enrichment (`run_ripple_fgsea`) | `fgsea`, `msigdbr` (Bioconductor) |
| Ligand-receptor integration (`run_ripple_lr`) | `nichenetr` (GitHub: `saeyslab/nichenetr`) |
| Input from `SingleCellExperiment` | `SingleCellExperiment` (Bioconductor) |
| GPU permutation testing (Stage 2) | Python 3, PyTorch (CUDA), NumPy < 2, SciPy |

**Note:** building the vignettes may pull in the `magick`
R package, which needs ImageMagick installed at the system level. If you hit a
`magick`/ImageMagick error, install it via your OS package manager (e.g.
`conda install -c conda-forge r-magick`) and
reinstall.

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

Four vignettes ship with the package:

| Vignette | Description |
|----------|-------------|
| `getting_started` | 5-minute end-to-end run on the bundled synthetic dataset. The fastest way to see what RIPPLE does. |
| `cosmx_nsclc_walkthrough` | Applied walkthrough on the public CosMx NSCLC dataset (He et al., 2022), reproducing paper Figure 4 panel by panel. Loads cached results from `inst/extdata/`, so the vignette renders without re-running the full pipeline. |
| `methods_positioning` | Landscape table comparing RIPPLE with nnSVG / SpatialDE / SPARK, STdiff (spatialGE), MISTy, COMMOT / CellChat, NicheNet, and BANKSY / CellCharter. Clarifies the distinction between co-localization tests and continuous gradient detection. |
| `benchmarks` | FDR calibration, power curves, and runtime measurements from the synthetic benchmark suite. |

Browse locally:

```r
browseVignettes("ripple")

# Or build the full pkgdown site (shipped _pkgdown.yml):
pkgdown::build_site()
```

### Interactive explainer

A standalone HTML page walks through how RIPPLE works, with diagrams.
It is not a vignette; open it directly in a browser:

[`vignettes/interactive_ripple_explainer.html`](vignettes/interactive_ripple_explainer.html)

GitHub shows the raw source rather than the rendered page for standalone
HTML files. To view it rendered, clone the repo and open the file locally:

```sh
git clone https://github.com/Maier-Lab/RIPPLE.git
open RIPPLE/vignettes/interactive_ripple_explainer.html   # macOS
# or:  xdg-open / start, depending on your OS
```

---

## Pipeline overview

Each stage is optional except Stage 1.

| Stage | Function(s) | Purpose |
|-------|-------------|---------|
| 1. Distance correlation | `run_ripple()` | Per-sample Poisson GLM + Fisher combined p-value with sign-consistency gate |
| 2. Merge and summarize | `merge_ripple_results()`, `compute_fisher_pval()` | Combines per-celltype results, recomputes Fisher p-values. IF you run_ripple(), you will get these out too, but you can use these functions if your run gets interrupted, for ex. |
| 3. Permutation validation | `run_permutation_tests()` (R) or `inst/python/run_permutation_gpu.py` (GPU) | Validates query specificity via label permutation. Warning: running without GPU is very slow, you may only want to do it later for genes of interest. |
| 4. Confounder control | `run_ripple_confounder()` | Bivariate GLM isolating query-specific from shared-niche effects |
| 5. Atlas figures | `run_ripple_atlas()`, `run_ripple_fgsea()`, `plot_gradient_volcano()`, `plot_gradient_curve()` | Multi-panel figures, pathway enrichment, contamination flagging |
| 6. Ligand-receptor integration | `run_ripple_lr()`, `classify_lr_artifacts()` | Matches gradient genes to L-R pairs via NicheNet |

Diagnostics:

- `ripple_plot_qc()` builds a one-glance multi-panel dashboard from the artefacts `run_ripple()` writes (per-sample distance density, cell composition, sign consistency, dispersion, specificity breakdown, and query-marker bleed-through).
- `classify_gene_specificity()` flags broad-expression genes (significant in many cell types at once). Three candidate sources to combine: this cross-cell-type flag, the `find_ambient_family_genes()` blocklist for Ig / J-chain / ribosomal / mitochondrial, and a user-supplied `query_signature_genes` list. None of these is a contamination measurement. You decide using domain knowledge.
- `plot_k_diagnostics()` helps pick `k_neighbors` before running the full pipeline.
- `check_spatial_autocorrelation()` computes Moran's I on Poisson residuals for selected genes, flagging cases where the independence assumption may be violated.
- `check_data()` and `load_metadata_only()` give fast metadata access without loading the full expression matrix.

---

## Core statistical model

For each gene in each target cell type, per biological replicate:

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

After running `run_ripple()` and `merge_ripple_results()`, run the
specificity check before interpreting individual genes:

```r
# 1. Classify gene specificity. Do this FIRST.
specificity <- classify_gene_specificity(all_results, fdr_threshold = 0.05)
table(specificity$specificity_class)
#   broad  moderate  specific
#      68       187       412
# Three classes: "specific" (1 cell type), "moderate" (2 up to
# broad_threshold-1), and "broad" (>= broad_threshold cell types).
# broad_threshold (default 4) is the single boundary -- raise it to flag
# fewer genes as broad. See ?classify_gene_specificity.

# 2. Pull the broad-class candidates aside.
broad_genes <- specificity[specificity_class == "broad"]$gene
clean_results <- all_results[!gene %in% broad_genes]

# 3. THEN look at individual genes.
plot_gradient_volcano(clean_results[cell_type == "CD8_T"], query_label = "Tumor")
```

Genes flagged as `broad` are significant in many cell types at once.
Common reasons: ambient RNA (query transcripts leaking into neighbouring
cells), housekeeping genes, or real shared biology (cytokines, MHC II,
stress programs). The flag is a heuristic, not a verdict. You decide
what is actually contamination using domain knowledge. The
`getting_started` vignette walks through the full curation workflow.

For a one-glance summary of the whole run, `ripple_plot_qc()` builds a
multi-panel QC dashboard from the artefacts `run_ripple()` writes:

```r
ripple_plot_qc(
  results_dir           = "ripple_output/tumor_ripple",
  query_signature_genes = c("EPCAM", "KRT8", "KRT19"),  # your query markers
  query_label           = "Tumor"
)
```

For k-selection, use `plot_k_diagnostics()` before the full pipeline:

```r
plot_k_diagnostics(my_spe, query_celltype = "Tumor", celltype_column = "cell_type")
```

---

RIPPLE was developed with Claude Opus 4.6/8 on Claude Code by Anthropic.

---
## Citation

A bibentry ships in `inst/CITATION` and can be retrieved with
`citation("ripple")`. It will be replaced once the accompanying
manuscript is posted on bioRxiv.

> *Citation (bioRxiv preprint in preparation).*


