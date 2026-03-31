# RIPPLE -- Replicate-Aware Inference of Paracrine Profiles via Likelihood Estimation

## Overview

RIPPLE detects distance-dependent gene expression gradients from any **query cell type** in spatial transcriptomics data (Xenium, MERFISH, CosMx, seqFISH, CODEX, Visium with deconvolution).

**Core question:** *"Which genes in cell type B change expression as a function of physical distance from cell type A?"*

**Origin:** Developed for HyMy (Hyperinflammatory Myeloid) cells in mouse lymph node Xenium data (CMM/mXenium project). Now generalized into a reusable, configurable pipeline driven entirely by environment variables.

---

## Pipeline Stages

| Stage | Type | Script(s) | Description |
|-------|------|-----------|-------------|
| **1. Distance Correlation** | Core | `hymy_distance_correlation_v2.R` | Per-sample Poisson GLM: gene expression ~ distance to query, with cell-size offset |
| **2. GPU Permutation** | Optional | `run_permutation_gpu.py` | GPU-accelerated null distribution (label permutation) -- validates query specificity |
| **3. Merge & Summarize** | Core | `merge_permutation_pvals.R`, `recompute_meta_summary.R`, `merge_distance_correlation_results.R` | Integrate permutation p-values, compute Fisher's combined p-value, merge across cell types |
| **4. Confounder Control** | Optional | `hymy_distance_correlation_stage2.R` | Bivariate GLM adding distance-to-control-cell-type -- isolates query-specific effects from shared niche. Requires `CONTROL_CELLTYPE`. |
| **5. Visualization** | Optional | `distance_correlation_atlas.R`, `hypothesis_visualizations.R`, `plot_decay_curves.R` | Volcanos, decay curves, heatmaps, dotplots, fGSEA pathway enrichment |
| **6. L-R Integration** | Optional | `gradient_lr_integration.R`, `gradient_lr_atlas.R`, `gradient_lr_biology_figure.R` | Match gradient genes to ligand-receptor pairs via NicheNet, artifact classification, curated biology figures |

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

## Generalization Status

All items from the original hardcoded HyMy/mXenium pipeline have been generalized. The pipeline is now fully configurable via environment variables.

| Element | Previous (Hardcoded) | Current (Generalized) | Status |
|---------|----------------------|-----------------------|--------|
| Input data path | Fixed path to HyMy Seurat RDS | `INPUT_PATH` env var -- user provides any Seurat `.rds` | Done |
| Query cell type | `HyMy_GMM` / `IL1B_myeloid` | `QUERY_CELLTYPE` env var | Done |
| Cell type column | `cell_type_with_HyMy` | `CELLTYPE_COLUMN` env var | Done |
| Control cell type (Stage 4) | `Monocyte` hardcoded | `CONTROL_CELLTYPE` env var | Done |
| Target cell types | 14 hardcoded types | Auto-detected from data, or `TARGET_CELLTYPES` env var | Done |
| Sample ID column | Hardcoded `sample_id` | `SAMPLE_COLUMN` env var (default: `sample_id`) | Done |
| Spatial coordinates | Hardcoded `spatial_x`/`spatial_y` | `X_COLUMN`/`Y_COLUMN` env vars with auto-detection | Done |
| Condition filtering | Hardcoded `group` with "naive"/"TDLN" | `CONDITION_COLUMN` + `CONDITION_VALUE` env vars | Done |
| Output directory | Fixed CeMM paths | `OUTPUT_DIR` env var (default: `./results`) | Done |
| Base paths | `N:/lab_maier/...` / `/nobackup/...` | Removed -- all paths from env vars | Done |
| Analysis subdirectory | Hardcoded name | `ANALYSIS_NAME` env var | Done |
| Display label | Hardcoded "HyMy" | `QUERY_LABEL` env var | Done |
| Contamination genes | HyMy signature genes | `QUERY_SIGNATURE_GENES` env var (comma-separated) | Done |
| AnnData for GPU | Fixed path | `ADATA_PATH` env var | Done |
| Annotation level routing | `ANNOTATION_LEVEL` env var with HyMy/L1 logic | Removed -- single `QUERY_CELLTYPE` parameter | Done |
| SLURM config | CeMM cluster-specific | Portable via `CONDA_SETUP`, `CONDA_ENV`, `CUDA_MODULE` env vars | Done |
| Priority genes | ~200 curated chemokines/cytokines/receptors | Optional user list or omit | Done |
| Positive controls | CSF3, IL33, CXCL12 | Optional user list or omit | Done |

---

## Environment Variable Reference

### Required

| Variable | Description |
|----------|-------------|
| `INPUT_PATH` | Path to Seurat object (`.rds`) with raw counts in the RNA assay |
| `QUERY_CELLTYPE` | Cell type label for the source population (e.g., `"Tumor"`) |
| `CELLTYPE_COLUMN` | Metadata column containing cell type annotations |

### Data Configuration (Optional)

| Variable | Default | Description |
|----------|---------|-------------|
| `OUTPUT_DIR` | `./results` | Output directory |
| `SAMPLE_COLUMN` | `sample_id` | Metadata column for sample/replicate IDs |
| `CONDITION_COLUMN` | -- | Metadata column for condition/group filtering |
| `CONDITION_VALUE` | -- | Which condition to analyze (if unset, all samples used) |
| `X_COLUMN` | auto-detect | X spatial coordinate column |
| `Y_COLUMN` | auto-detect | Y spatial coordinate column |
| `TARGET_CELLTYPES` | auto-detect all | Comma-separated target cell types |
| `ADATA_PATH` | -- | AnnData `.h5ad` file for GPU permutation |

### Analysis Configuration (Optional)

| Variable | Default | Description |
|----------|---------|-------------|
| `QUERY_LABEL` | Same as `QUERY_CELLTYPE` | Display name for plots |
| `QUERY_SIGNATURE_GENES` | -- | Comma-separated markers for contamination check |
| `CONTROL_CELLTYPE` | -- | Control cell type for Stage 4 confounder analysis |
| `ANALYSIS_NAME` | `hymy_distance_correlation_v2` | Subdirectory name for output |
| `K_NEIGHBORS` | `1` | k for kNN distance calculation |
| `N_PERMUTATIONS` | `0` | Label permutations (0 = skip; use GPU Stage 2 for production) |

### HPC / SLURM (Optional)

| Variable | Default | Description |
|----------|---------|-------------|
| `CONDA_SETUP` | -- | Path to conda setup script |
| `CONDA_ENV` | -- | Conda environment name |
| `CUDA_MODULE` | -- | CUDA module to load for GPU permutation |
| `CELLTYPE_INDEX` | -- | 1-based index for SLURM array jobs |

### Internal / Advanced

| Parameter | Default | Description | When to Change |
|-----------|---------|-------------|----------------|
| `MAX_DISTANCE_UM` | 200 | Maximum distance (um) to consider | Reduce for contact-dependent (~50); increase for morphogens (~500) |
| `MIN_EXPR_CELLS` | 25 | Minimum non-zero cells per sample per gene | Lower for rare cell types |
| `MIN_EXPR_PCT` | 0.01 | Minimum fraction expressing | |
| `FDR_THRESHOLD` | 0.05 | Fisher FDR cutoff | 0.1 for discovery; 0.01 for high-confidence |
| `SIGN_CONSISTENCY_THRESHOLD` | 1.0 | Required fraction of samples agreeing on beta direction | Relax to 0.75 with N>=6 samples |

---

## Usage

### Minimal (discovery)

```bash
export INPUT_PATH="/path/to/seurat.rds"
export QUERY_CELLTYPE="Tumor"
export CELLTYPE_COLUMN="cell_type"

Rscript scripts/hymy_distance_correlation_v2.R
Rscript scripts/recompute_meta_summary.R
Rscript scripts/merge_distance_correlation_results.R
```

Results land in `./results/spatial_analysis_Tumor/`.

### With condition filtering

```bash
export CONDITION_COLUMN="treatment"
export CONDITION_VALUE="treated"
```

### Full pipeline (publication-ready)

```bash
export INPUT_PATH="/path/to/seurat.rds"
export QUERY_CELLTYPE="Tumor"
export CELLTYPE_COLUMN="cell_type"
export OUTPUT_DIR="./results"
export CONDITION_COLUMN="treatment"
export CONDITION_VALUE="treated"

# Stage 1: Core analysis
N_PERMUTATIONS=0 Rscript scripts/hymy_distance_correlation_v2.R

# Stage 2: GPU permutation (optional)
ADATA_PATH="/path/to/adata.h5ad" python scripts/run_permutation_gpu.py

# Stage 3: Merge + Fisher's p-value
Rscript scripts/merge_permutation_pvals.R
Rscript scripts/recompute_meta_summary.R
Rscript scripts/merge_distance_correlation_results.R

# Stage 4: Confounder control (optional)
CONTROL_CELLTYPE="Macrophage" Rscript scripts/hymy_distance_correlation_stage2.R

# Stage 5: Visualization
Rscript scripts/distance_correlation_atlas.R
Rscript scripts/hypothesis_visualizations.R

# Stage 6: L-R integration (requires NicheNet)
Rscript scripts/gradient_lr_integration.R
Rscript scripts/gradient_lr_atlas.R
```

### SLURM cluster (parallel across cell types)

```bash
# Array job: one job per target cell type (auto-detected)
sbatch --array=1-20 scripts/run_hymy_distance_correlation_v2.sh
```

### Input data requirements

The Seurat `.rds` object must have:
- **Raw counts** in `RNA@counts` (not normalized)
- **Spatial coordinates** in metadata (any column names -- auto-detected, or set `X_COLUMN`/`Y_COLUMN`)
- **Cell type annotations** in a metadata column matching `CELLTYPE_COLUMN`, including the `QUERY_CELLTYPE` label
- **Sample/replicate IDs** in a metadata column matching `SAMPLE_COLUMN` (default: `sample_id`) -- minimum 3 replicates

---

## R Package Structure

RIPPLE is an installable R package. The core analysis functions live in `R/`, while the original standalone scripts are preserved in `scripts/` (and `inst/scripts/`, `inst/slurm/`).

### Two ways to use RIPPLE

**As an R package** (recommended):
```r
library(ripple)
results <- run_ripple(
  input_path = "my_seurat.rds",
  query_celltype = "Tumor",
  celltype_column = "cell_type"
)
```

**As standalone scripts** (for SLURM / backward compat):
```bash
export INPUT_PATH="my_seurat.rds"
export QUERY_CELLTYPE="Tumor"
export CELLTYPE_COLUMN="cell_type"
Rscript scripts/hymy_distance_correlation_v2.R
```

### Package directory layout

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
├── tests/testthat/          # Unit tests (31 tests, all passing)
├── scripts/                 # Standalone scripts (env var interface, for SLURM)
├── inst/
│   ├── scripts/             # Copies of standalone scripts (installed with package)
│   ├── slurm/               # SLURM job templates
│   └── python/              # GPU permutation script
├── README.md
└── RIPPLE_user_guide.html   # Statistical methodology guide
```

### Exported functions (28 total)

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

### Development workflow

```r
devtools::load_all()    # Load package for development (no install)
devtools::test()        # Run unit tests
devtools::document()    # Regenerate NAMESPACE + man/ from roxygen2
devtools::check()       # Full R CMD check (needs Rtools on Windows)
devtools::install()     # Install locally
```

### Standalone scripts (scripts/)

These are the original pipeline scripts, configured via environment variables. They are independent of the R package and can be run directly with `Rscript`. Useful for SLURM array jobs.

| Script | Stage | Role |
|--------|-------|------|
| `hymy_distance_correlation_v2.R` | 1 (Core) | Poisson GLM distance correlation |
| `run_permutation_gpu.py` | 2 (Optional) | GPU-accelerated permutation testing |
| `merge_permutation_pvals.R` | 3 (Core) | Integrate GPU permutation p-values |
| `recompute_meta_summary.R` | 3 (Core) | Fisher's combined p-value + sign gate |
| `merge_distance_correlation_results.R` | 3 (Core) | Merge per-celltype results |
| `hymy_distance_correlation_stage2.R` | 4 (Optional) | Bivariate GLM confounder control |
| `distance_correlation_atlas.R` | 5 (Optional) | Volcanos, heatmaps, fGSEA |
| `hypothesis_visualizations.R` | 5 (Optional) | Per-sample decay curves |
| `plot_decay_curves.R` | 5 (Optional) | Decay curve plots |
| `gradient_lr_integration.R` | 6 (Optional) | L-R pair matching via NicheNet |
| `gradient_lr_atlas.R` | 6 (Optional) | L-R summary figures + artifact classification |
| `gradient_lr_biology_figure.R` | 6 (Optional) | Curated L-R biology figure |
| `config.R` | All | Env var resolution (sourced by all scripts) |
| `utils.R` | All | Shared utilities |
| `load_data.R` | All | Data loading functions |
| `plot_binary_assumption.R` | QC | Verify Poisson assumption |

---

## Downstream Analysis Details

### Atlas Visualization (Stage 5)

`distance_correlation_atlas.R` produces:

| Panel | Output | Description |
|-------|--------|-------------|
| P1 | `gene_counts_by_celltype.pdf` | Stacked bar: induced vs repressed, by specificity |
| P2 | `dotplot_specific_genes.pdf` | Top 5 cell-type-specific genes per cell type |
| P3 | `multi_volcano.pdf` | Multi-panel volcano with contamination flagging |
| P5 | `specificity_distribution.pdf` | Gene specificity breakdown |
| P6 | `specific_genes_heatmap.pdf` | Heatmap of top cell-type-specific genes |
| P7 | `contamination_candidates.csv` | Genes significant in >=4 cell types (suspect artifacts) |
| P8-9 | `fgsea_hallmark_*.pdf` | fGSEA pathway enrichment (Hallmark gene sets via msigdbr) |
| D2 | `sample_contribution_heatmap.pdf` | Per-sample coefficient patterns |
| E1 | `ligand_receptor_pairs.pdf` | Known L-R pairs in spatial gradients |
| E3 | `query_signature_leakage_check.pdf` | Query cell markers detected in other cell types (contamination) |
| P10 | `stage2_classification_breakdown.pdf` | Query-specific vs niche-driven by cell type |

### Decay Curves (Stage 5)

`hypothesis_visualizations.R` creates per-sample reproducibility plots:
- Thin lines = individual replicates
- Thick line + 95% CI ribbon = pooled mean
- Red = induced near query, Blue = repressed, Gray = niche-driven
- Meta-analysis coefficient + FDR annotated per panel

### L-R Integration (Stage 6)

`gradient_lr_integration.R` runs three analyses per cell type:

| Part | Method | Question |
|------|--------|----------|
| 1: Direct L-R Mapping | Match gradient receptors to query-expressed ligands (per-replicate) | Which L-R pairs have gradient + expression support? |
| 2: NicheNet Activity | `predict_ligand_activities()` on query-induced genes | Which query ligands best predict the gradient pattern? |
| 3: Enrichment | Fisher's exact on NicheNet target predictions | Are predicted targets enriched in gradient genes? |

**Directions:**
- **A (Query -> Target)**: Full pipeline (Parts 1-3 + combined scoring)
- **B (Target -> Query)**: Direct L-R mapping only

**Artifact Classification** (`gradient_lr_atlas.R`): 4-tier system (artifact / suspect / low_confidence / clean) based on query signature leakage, contamination candidates, lineage-mismatched markers, and expression thresholds. All downstream figures filter to `clean` only.

**Biology Figure** (`gradient_lr_biology_figure.R`): Curated side-by-side Direction A + B layout with dotplot, gradient lollipop, and reproducibility panels.

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

### Condition-Specific Cell Types

Query cell types may be rare in one condition but abundant in another. Use `CONDITION_COLUMN` and `CONDITION_VALUE` to restrict analysis to the condition where the query is abundant:

```bash
export CONDITION_COLUMN="treatment"
export CONDITION_VALUE="treated"
```

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

### Spatial Confounding

A positive proximity result could mean:
- Biological interaction, OR
- Both cell types independently localize to the same region

**Solutions:**
1. Compare within same spatial domain
2. Use control cell types (Stage 4)
3. Permutation testing (Stage 2)

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
If top "induced" genes in a target cell type are known markers of the query cell type, suspect ambient RNA / segmentation artifacts. Provide query markers via `QUERY_SIGNATURE_GENES` to enable automatic flagging.

### Overdispersion
Median dispersion 0.3-0.6 across cell types confirms Poisson is appropriate. If consistently >2, consider quasi-Poisson or negative binomial.

---

## Troubleshooting

| Problem | Fix |
|---------|-----|
| DOS line endings on HPC | `sed -i 's/\r$//' *.sh *.R` |
| Duplicate column error (`cell_id`) | Use `keep.rownames = "barcode"` |
| facet_wrap + coord_fixed error | Use `plot_spatial_by_sample()` from `utils.R` |
| GPU OOM in permutation | Reduce `batch_size` in `run_permutation_gpu.py` |
| Missing query cells in some samples | Script handles gracefully; samples with <30 query cells excluded |
| rbindlist column mismatch | Use `rbindlist(results, fill = TRUE)` |
