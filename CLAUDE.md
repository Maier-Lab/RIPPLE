# RIPPLE -- Replicate-Aware Inference of Paracrine Profiles via Likelihood Estimation

## Overview

RIPPLE is an R package that detects distance-dependent gene expression gradients from any **query cell type** in spatial transcriptomics data (Xenium, MERFISH, CosMx, seqFISH, CODEX, Visium with deconvolution).

**Core question:** *"Which genes in cell type B change expression as a function of physical distance from cell type A?"*

For the standalone script version (env-var-driven, SLURM-compatible), see [HyMy-distance-correlation-analysis](https://github.com/CMangana/HyMy-distance-correlation-analysis).

---

## Package Structure

```
RIPPLE/
├── DESCRIPTION              # Package metadata + dependencies
├── NAMESPACE                # Exported functions + imports
├── LICENSE
├── R/                       # Package source (functions with roxygen2 docs)
│   ├── config.R             # Package options system (.onLoad, ripple_config)
│   ├── data_loading.R       # load_seurat, load_metadata_only, check_data
│   ├── spatial.R            # get_coord_columns, build_knn_graph, distance_to_type
│   ├── glm.R                # fit_poisson, fit_poisson_controlled, classify_decay
│   ├── meta_analysis.R      # run_meta_analysis, compute_fisher_pval
│   ├── permutation.R        # run_permutation_test, merge_permutation_results
│   ├── pipeline.R           # run_ripple, run_ripple_confounder, merge_ripple_results
│   ├── visualization.R      # plot_gradient_volcano, plot_decay_curve, spatial plots
│   └── utils.R              # entropy, enrichment, gene scoring helpers
├── tests/testthat/          # Unit tests
├── inst/
│   ├── scripts/             # Copies of standalone scripts (installed with package)
│   ├── slurm/               # SLURM job templates
│   └── python/              # GPU permutation script
└── README.md
```

---

## Exported Functions (28 total)

| Module | Functions | Purpose |
|--------|-----------|---------|
| `pipeline.R` | `run_ripple`, `run_ripple_confounder`, `merge_ripple_results` | Main entry points |
| `config.R` | `ripple_config` | Get/set package options |
| `data_loading.R` | `load_seurat`, `load_metadata_only`, `check_data` | Data loading |
| `glm.R` | `fit_poisson`, `fit_poisson_controlled`, `classify_decay_pattern` | Core statistical models |
| `meta_analysis.R` | `run_meta_analysis`, `compute_fisher_pval` | Cross-replicate inference |
| `permutation.R` | `run_permutation_test`, `merge_permutation_results` | Null distribution |
| `spatial.R` | `get_coord_columns`, `build_knn_graph`, `build_radius_graph`, `calculate_distance_to_type`, `get_neighbor_cell_types`, `calculate_neighbor_composition` | Spatial utilities |
| `visualization.R` | `plot_spatial_scatter`, `plot_spatial_single`, `plot_spatial_by_sample`, `plot_gradient_volcano`, `plot_decay_curve`, `plot_violin` | Plotting |
| `utils.R` | `shannon_entropy`, `calculate_enrichment`, `score_gene_signature`, `score_multiple_modules`, `permutation_pvalue`, `calculate_neighborhood_entropy` | Helpers |

---

## Usage

### Quick start

```r
library(ripple)
results <- run_ripple(
  input_path      = "my_seurat.rds",
  query_celltype  = "Tumor",
  celltype_column = "cell_type",
  output_dir      = "./results"
)
```

### With condition filtering

```r
results <- run_ripple(
  input_path       = "my_seurat.rds",
  query_celltype   = "Tumor",
  celltype_column  = "cell_type",
  output_dir       = "./results",
  condition_column = "treatment",
  condition_value  = "treated"
)
```

### Merge results across cell types

```r
merged <- merge_ripple_results(
  results_dir    = "./results/spatial_analysis_Tumor/tumor_ripple",
  recompute_fisher = TRUE
)
```

### Confounder control

```r
stage2 <- run_ripple_confounder(
  input_path       = "my_seurat.rds",
  results_dir      = "./results/spatial_analysis_Tumor/tumor_ripple",
  query_celltype   = "Tumor",
  celltype_column  = "cell_type",
  control_celltype = "CAF"
)
```

### Atlas figures

```r
run_ripple_atlas(
  results_dir = "./results/spatial_analysis_Tumor/tumor_ripple",
  output_dir  = "./results/spatial_analysis_Tumor/tumor_ripple/atlas",
  query_label = "Tumor",
  run_fgsea   = TRUE,
  organism    = "human"
)
```

---

## Core Statistical Model

For each gene in each target cell type, per biological replicate:

```r
glm(counts ~ distance_to_query + offset(log(total_counts)), family = poisson)
```

- **`distance_to_query`**: Euclidean distance (um) to nearest query cell (k=1 NN)
- **`offset(log(total_counts))`**: cell-size correction -- converts to rates, controls for ambient RNA / segmentation differences
- **Coefficient (beta)**: log-rate change per um. Negative = expression increases near query cells (induced). Positive = decreases (repressed).

Per-sample coefficients are combined via **Fisher's combined p-value** with sign consistency gate (all replicates must agree on direction). `fisher_fdr` is the **primary significance metric** for all downstream analyses.

### Why This Design

- **Raw counts + Poisson**: Xenium data is integer counts; the offset handles normalization internally
- **Per-sample fitting**: Avoids pseudoreplication (true N = mice/patients, not cells)
- **Fisher's + sign gate**: Equal weight per replicate; requires reproducibility across all biological replicates
- **Cell-size offset**: Critical -- without it, cells near query may appear to express more of everything due to ambient RNA or segmentation artifacts

---

## Pipeline Stages

| Stage | Type | R Package Function(s) | Description |
|-------|------|-----------------------|-------------|
| **1. Distance Correlation** | Core | `run_ripple()` | Per-sample Poisson GLM: gene expression ~ distance to query, with cell-size offset |
| **2. GPU Permutation** | Optional | `run_permutation_test()` | GPU-accelerated null distribution (label permutation) -- validates query specificity |
| **3. Merge & Summarize** | Core | `merge_ripple_results()`, `compute_fisher_pval()` | Integrate permutation p-values, compute Fisher's combined p-value, merge across cell types |
| **4. Confounder Control** | Optional | `run_ripple_confounder()` | Bivariate GLM adding distance-to-control-cell-type -- isolates query-specific effects from shared niche |
| **5. Visualization** | Optional | `run_ripple_atlas()`, `plot_gradient_volcano()`, `plot_decay_curve()` | Volcanos, decay curves, heatmaps, dotplots, fGSEA pathway enrichment |
| **6. L-R Integration** | Optional | `run_ripple_lr()`, `classify_lr_artifacts()` | Match gradient genes to ligand-receptor pairs via NicheNet, artifact classification |

---

## Configuration

RIPPLE can be configured via environment variables (read by `ripple_config()` and `.onLoad()`) or by passing arguments directly to R functions. Function arguments take precedence over environment variables.

### Key parameters

| Parameter / Env Var | Default | Description |
|---------------------|---------|-------------|
| `input_path` / `INPUT_PATH` | -- | Path to Seurat object (`.rds`) with raw counts |
| `query_celltype` / `QUERY_CELLTYPE` | -- | Cell type label for the source population |
| `celltype_column` / `CELLTYPE_COLUMN` | -- | Metadata column containing cell type annotations |
| `output_dir` / `OUTPUT_DIR` | `./results` | Output directory |
| `sample_column` / `SAMPLE_COLUMN` | `sample_id` | Metadata column for sample/replicate IDs |
| `condition_column` / `CONDITION_COLUMN` | -- | Metadata column for condition/group filtering |
| `condition_value` / `CONDITION_VALUE` | -- | Which condition to analyze |
| `target_celltypes` / `TARGET_CELLTYPES` | auto-detect all | Comma-separated target cell types |
| `query_label` / `QUERY_LABEL` | Same as query_celltype | Display name for plots |
| `query_signature_genes` / `QUERY_SIGNATURE_GENES` | -- | Comma-separated markers for contamination check |
| `control_celltype` / `CONTROL_CELLTYPE` | -- | Control cell type for Stage 4 |
| `k_neighbors` / `K_NEIGHBORS` | `1` | k for kNN distance calculation |

### Internal / Advanced

| Parameter | Default | Description | When to Change |
|-----------|---------|-------------|----------------|
| `MAX_DISTANCE_UM` | 200 | Maximum distance (um) to consider | Reduce for contact-dependent (~50); increase for morphogens (~500) |
| `MIN_EXPR_CELLS` | 25 | Minimum non-zero cells per sample per gene | Lower for rare cell types |
| `MIN_EXPR_PCT` | 0.01 | Minimum fraction expressing | |
| `FDR_THRESHOLD` | 0.05 | Fisher FDR cutoff | 0.1 for discovery; 0.01 for high-confidence |
| `SIGN_CONSISTENCY_THRESHOLD` | 1.0 | Required fraction of samples agreeing on beta direction | Relax to 0.75 with N>=6 samples |

---

## Spatial Statistics (CRITICAL)

### N = Samples, Not Cells!

**True N = number of biological replicates**, not number of cells.

```r
# WRONG - pseudoreplication
wilcox.test(data[condition == "control"]$entropy,
            data[condition == "treated"]$entropy)

# RIGHT - sample-level aggregation
sample_stats <- data[, .(median_entropy = median(entropy)),
                     by = .(sample_id, condition)]
wilcox.test(sample_stats[condition == "control"]$median_entropy,
            sample_stats[condition == "treated"]$median_entropy,
            exact = FALSE)
```

### Statistical Testing Summary

| Analysis Type | Unit | Test |
|--------------|------|------|
| Condition comparison | Sample median | Wilcoxon rank-sum |
| Paired comparison | Per-sample difference | Wilcoxon signed-rank |
| Spatial enrichment | Per-sample log2 enrichment | One-sample Wilcoxon (H0: 0) |
| Permutation test | Stratify by sample | Within-sample shuffling |

### Always:
1. Set `set.seed()` for reproducibility
2. Aggregate to sample-level before testing
3. Use Wilcoxon (not t-test) for small N
4. Apply FDR correction across tests
5. Report effect sizes (mean +/- SD, fold change)
6. Be explicit about N in reports

### Expression Filtering (Two-Tier)

- **Tier 1 (strict, regular genes):** >= max(1% of cells, 25) expressing in ALL valid samples
- **Tier 2 (lenient, priority genes):** floor threshold (25 cells), pass in >= 2 samples. Rescues biologically important but sparse genes.

### Confounder Control (Stage 4)

Choose a control cell type that (1) co-localizes with query cells and (2) is biologically distinct:
- Query = inflammatory myeloid -> Control = monocytes (shared niche, different function)
- Query = tumor cells -> Control = CAFs (both in tumor core)
- Query = Tregs -> Control = CD4 T cells (shared T cell zones)

If no obvious confounder exists, skip Stage 4 and rely on permutation testing.

**Gene classification after Stage 4:**

| Classification | Criteria |
|---------------|----------|
| `query_specific` | Fisher FDR < 0.05 in Stage 4, same sign as Stage 1 |
| `enhanced` | query_specific + abs(stage4_coef) > abs(stage1_coef) x 1.1 |
| `niche_driven` | FDR >= 0.05 AND coefficient attenuated >50% |
| `underpowered` | FDR >= 0.05 AND coefficient preserved >=50% |

---

## Development Workflow

```r
devtools::load_all()    # Load package for development (no install)
devtools::test()        # Run unit tests
devtools::document()    # Regenerate NAMESPACE + man/ from roxygen2
devtools::check()       # Full R CMD check (needs Rtools on Windows)
devtools::install()     # Install locally
```

---

## Dependencies

### R Packages

```r
# Core (required)
library(Seurat)       # Data loading, count extraction
library(data.table)   # Core data manipulation
library(ggplot2)      # All plotting
library(patchwork)    # Multi-panel figure assembly
library(RANN)         # Fast kNN
library(meta)         # Meta-analysis across replicates
library(ggrepel)      # Volcano label placement
library(Matrix)       # Sparse matrix operations
library(scales)       # Axis formatting
library(dplyr)        # Data wrangling
library(tidyr)        # Reshaping

# Visualization (Stage 5)
library(pheatmap)     # Heatmaps
library(fgsea)        # Pathway enrichment
library(msigdbr)      # Gene set collections

# L-R Integration (Stage 6)
library(nichenetr)    # Ligand-receptor database + activity prediction
```

### Python (GPU Permutation, Stage 2)

- PyTorch with CUDA support
- NumPy 1.x (not 2.x -- incompatible with PyTorch 2.1)
- SciPy

---

## Implementation Notes

### Seurat Metadata
- Use `keep.rownames = "barcode"` when converting to data.table (avoids duplicate `cell_id` column errors)
- `cell_id` column is truncated -- use full barcode for merging

### Spatial Visualization
- Use `coord_fixed()` for spatial data
- Use individual plots + `wrap_plots()` instead of `facet_wrap(scales="free")` + `coord_fixed()` (ggplot2 bug)

### Contamination Check
If top "induced" genes in a target cell type are known markers of the query cell type, suspect ambient RNA / segmentation artifacts. Provide query markers via `query_signature_genes` to enable automatic flagging.

### Overdispersion
Median dispersion 0.3-0.6 across cell types confirms Poisson is appropriate. If consistently >2, consider quasi-Poisson or negative binomial.

---

## Troubleshooting

| Problem | Fix |
|---------|-----|
| Duplicate column error (`cell_id`) | Use `keep.rownames = "barcode"` |
| facet_wrap + coord_fixed error | Use `plot_spatial_by_sample()` |
| GPU OOM in permutation | Reduce `batch_size` in GPU permutation script |
| Missing query cells in some samples | Handled gracefully; samples with <30 query cells excluded |
| rbindlist column mismatch | Use `rbindlist(results, fill = TRUE)` |
