# RIPPLE

**Replicate-Aware Inference of Paracrine Profiles via Likelihood Estimation**

RIPPLE detects distance-dependent gene expression gradients from any query cell type in spatial transcriptomics data.

---

## Overview

RIPPLE tries to answer a common question in spatial biology: *"Which genes in cell type B change expression as a function of physical distance from cell type A?"*

Given a spatial transcriptomics dataset with cell type annotations and biological replicates, RIPPLE fits per-sample Poisson GLMs to model gene expression as a function of distance to a query cell population, then combines evidence across replicates using Fisher's combined p-value with a sign consistency gate. The result is a ranked list of genes with distance-dependent expression gradients, along with effect sizes, significance measures, and confounder controls.

**Supported platforms:** Xenium, MERFISH, CosMx, CODEX. In theory, any platform that provides single-cell resolved spatial coordinates and count data.

**Use cases:**
- Tumor border effects on neighboring immune and stromal cells
- Inflammatory myeloid cell influence on tissue microenvironment
- CAF-mediated remodeling of the tumor niche
- Treg-mediated immunosuppression in tissue

---

## Requirements

- **Spatial transcriptomics data** with single-cell resolution (coordinates + counts)
- **Cell type annotations** including the query population of interest
- **Biological replicates** (minimum 3, ideally 4+)
- **Raw count matrix** (not normalized, the model handles normalization internally via cell-size offset)

---

## Installation

### R Dependencies

```r
# Core packages (CRAN)
install.packages(c(
  "data.table", "ggplot2", "patchwork", "RANN", "meta",
  "ggrepel", "Matrix", "scales", "rstatix", "dplyr", "tidyr"
))

# Seurat
install.packages("Seurat")

# Visualization
install.packages("pheatmap")

# Bioconductor packages (pathway enrichment)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("fgsea", "msigdbr"))

# NicheNet (for L-R integration, Stage 6)
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("saeyslab/nichenetr")
```

**Full R package list:**

| Package | Source | Used in |
|---------|--------|---------|
| `Seurat` | CRAN | Data loading, count extraction |
| `data.table` | CRAN | Core data manipulation |
| `ggplot2` | CRAN | All plotting |
| `patchwork` | CRAN | Multi-panel figure assembly |
| `RANN` | CRAN | Fast k-nearest neighbor search |
| `meta` | CRAN | Meta-analysis across replicates |
| `ggrepel` | CRAN | Non-overlapping volcano labels |
| `Matrix` | CRAN | Sparse matrix operations |
| `scales` | CRAN | Axis formatting |
| `rstatix` | CRAN | Tidy statistical tests |
| `dplyr` | CRAN | Data wrangling |
| `tidyr` | CRAN | Reshaping (pivot_wider/longer) |
| `pheatmap` | CRAN | Heatmap visualization |
| `fgsea` | Bioconductor | Fast gene set enrichment analysis |
| `msigdbr` | Bioconductor | MSigDB gene set collections |
| `nichenetr` | GitHub | Ligand-receptor database and activity prediction |

### Python Dependencies (GPU Permutation Testing)

Required only for Stage 2 (GPU-accelerated permutation testing):

- **PyTorch** with CUDA support
- **NumPy** 1.x (not 2.x -- incompatible with PyTorch 2.1)
- **SciPy**

```bash
pip install torch numpy scipy
```

---

## Minimal Example

The absolute minimum to run RIPPLE -- three environment variables and one command:

```bash
export INPUT_PATH="/path/to/my/seurat.rds"
export QUERY_CELLTYPE="Tumor"
export CELLTYPE_COLUMN="cell_type"

Rscript scripts/hymy_distance_correlation_v2.R
```

This will auto-detect spatial coordinates, sample IDs, and all non-query cell types, then fit per-sample Poisson GLMs for every gene in every target cell type.

---

## Quick Start

```bash
export INPUT_PATH="/path/to/my/seurat.rds"
export QUERY_CELLTYPE="Tumor"
export CELLTYPE_COLUMN="cell_type"

# Run core analysis (Stage 1)
Rscript scripts/hymy_distance_correlation_v2.R

# Merge results (Stage 3)
Rscript scripts/recompute_meta_summary.R
Rscript scripts/merge_distance_correlation_results.R
```

---

## Configuration

All configuration is via environment variables. Set them before running any script.

### Required

| Variable | Description |
|----------|-------------|
| `INPUT_PATH` | Path to a Seurat object (`.rds`) with raw counts |
| `QUERY_CELLTYPE` | Cell type label for the source population (e.g., `"Tumor"`) |
| `CELLTYPE_COLUMN` | Metadata column containing cell type annotations |

### Data Configuration (Optional)

| Variable | Default | Description |
|----------|---------|-------------|
| `OUTPUT_DIR` | `./results` | Output directory for all results |
| `SAMPLE_COLUMN` | `sample_id` | Metadata column containing sample/replicate IDs |
| `CONDITION_COLUMN` | -- | Metadata column for condition/group filtering (e.g., `"treatment"`) |
| `CONDITION_VALUE` | -- | Which condition to analyze (e.g., `"treated"`). If unset, all samples are used |
| `X_COLUMN` | auto-detect | Metadata column for X spatial coordinates |
| `Y_COLUMN` | auto-detect | Metadata column for Y spatial coordinates |
| `TARGET_CELLTYPES` | auto-detect all | Comma-separated list of target cell types to analyze (e.g., `"CD8_T,Macrophage"`) |
| `ADATA_PATH` | -- | Path to AnnData `.h5ad` file for GPU permutation (Stage 2) |

### Analysis Configuration (Optional)

| Variable | Default | Description |
|----------|---------|-------------|
| `QUERY_LABEL` | Same as `QUERY_CELLTYPE` | Display name for the query cell type in plot titles |
| `QUERY_SIGNATURE_GENES` | -- | Comma-separated marker genes for contamination check (e.g., `"CD274,PDCD1"`) |
| `CONTROL_CELLTYPE` | -- | Control cell type for Stage 4 confounder analysis |
| `ANALYSIS_NAME` | `hymy_distance_correlation_v2` | Subdirectory name for output organization |
| `K_NEIGHBORS` | `1` | Number of nearest query cells for distance calculation |
| `N_PERMUTATIONS` | `0` | Number of label permutations (0 = skip; use GPU Stage 2 instead for production) |

### HPC / SLURM (Optional)

| Variable | Default | Description |
|----------|---------|-------------|
| `CONDA_SETUP` | -- | Path to conda setup script (e.g., `/path/to/conda.sh`) |
| `CONDA_ENV` | -- | Conda environment name for R dependencies |
| `CUDA_MODULE` | -- | CUDA module to load for GPU permutation (e.g., `"cuda/12.2"`) |
| `CELLTYPE_INDEX` | -- | 1-based index for SLURM array jobs (run one target cell type per job) |

---

## Pipeline Stages

### Stage 1 (Core): Distance Correlation

**Script:** `hymy_distance_correlation_v2.R`

For each gene in each non-query cell type, fits a per-sample Poisson GLM:

```r
glm(counts ~ distance_to_query + offset(log(total_counts)), family = poisson)
```

- `distance_to_query`: Euclidean distance (um) to nearest query cell (k-NN, default k=1)
- `offset(log(total_counts))`: cell-size correction, converts raw counts to rates
- **Coefficient (beta)**: log-rate change per um. Negative = expression increases near query cells (induced). Positive = expression decreases near query (repressed).

Per-sample coefficients are combined via Fisher's combined p-value. A sign consistency gate requires all replicates to agree on the direction of the effect.

### Stage 2 (Optional): GPU Permutation Testing

**Script:** `run_permutation_gpu.py`

GPU-accelerated null distribution via label permutation. Validates that the observed distance-expression relationship is specific to the query cell type and not a spatial artifact. Shuffles query cell labels across all cells and re-fits the GLM to build a null distribution.

### Stage 3 (Core): Merge and Summarize

**Scripts:** `merge_permutation_pvals.R`, `recompute_meta_summary.R`, `merge_distance_correlation_results.R`

Integrates GPU permutation p-values into the R results, recomputes Fisher's combined p-values, and merges per-cell-type results into summary tables.

### Stage 4 (Optional): Confounder Control

**Script:** `hymy_distance_correlation_stage2.R`

Requires `CONTROL_CELLTYPE` to be set. Bivariate GLM adding distance to a control cell type as a covariate. Isolates query-specific effects from shared spatial niche effects. Classifies genes as:

| Classification | Criteria |
|----------------|----------|
| `query_specific` | Significant in bivariate model, same sign as Stage 1 |
| `enhanced` | Query-specific and effect size increased >10% vs Stage 1 |
| `niche_driven` | Not significant and coefficient attenuated >50% |
| `underpowered` | Not significant but coefficient preserved (>=50%) |

### Stage 5 (Optional): Visualization

**Scripts:** `distance_correlation_atlas.R`, `hypothesis_visualizations.R`, `plot_decay_curves.R`

Produces publication-quality figures: multi-panel volcanos, dotplots of cell-type-specific genes, heatmaps, contamination flagging, fGSEA pathway enrichment (Hallmark gene sets), per-sample decay curves with reproducibility ribbons, and diagnostic panels.

### Stage 6 (Optional): Ligand-Receptor Integration

**Scripts:** `gradient_lr_integration.R`, `gradient_lr_atlas.R`, `gradient_lr_biology_figure.R`

Matches gradient genes to ligand-receptor pairs using NicheNet. Runs three analyses per cell type: (1) direct L-R mapping of gradient receptors to query-expressed ligands, (2) NicheNet ligand activity prediction, and (3) Fisher's exact enrichment of NicheNet target predictions among gradient genes. Includes a 4-tier artifact classification system.

---

## Full Pipeline Example

```bash
# --- Required ---
export INPUT_PATH="/path/to/my/seurat.rds"
export QUERY_CELLTYPE="Tumor"
export CELLTYPE_COLUMN="cell_type"

# --- Optional data configuration ---
export OUTPUT_DIR="./results"
export SAMPLE_COLUMN="patient_id"
export QUERY_LABEL="Tumor"
export QUERY_SIGNATURE_GENES="EPCAM,KRT8,KRT18"

# Stage 1 (Core): Distance correlation
Rscript scripts/hymy_distance_correlation_v2.R

# Stage 2 (Optional): GPU permutation testing
export ADATA_PATH="/path/to/adata.h5ad"
python scripts/run_permutation_gpu.py

# Stage 3 (Core): Merge results
Rscript scripts/merge_permutation_pvals.R
Rscript scripts/recompute_meta_summary.R
Rscript scripts/merge_distance_correlation_results.R

# Stage 4 (Optional): Confounder control -- requires a suitable control cell type
export CONTROL_CELLTYPE="CAF"
Rscript scripts/hymy_distance_correlation_stage2.R

# Stage 5 (Optional): Visualization
Rscript scripts/distance_correlation_atlas.R
Rscript scripts/hypothesis_visualizations.R
Rscript scripts/plot_decay_curves.R

# Stage 6 (Optional): L-R integration -- requires NicheNet
Rscript scripts/gradient_lr_integration.R
Rscript scripts/gradient_lr_atlas.R
Rscript scripts/gradient_lr_biology_figure.R
```

---

## SLURM Cluster Usage

Stage 1 can be parallelized across target cell types using SLURM array jobs. The target cell types are auto-detected at runtime (or set via `TARGET_CELLTYPES`). Each array task processes one cell type.

```bash
# Set environment variables in your SLURM script or export them before submission
export INPUT_PATH="/path/to/seurat.rds"
export QUERY_CELLTYPE="Tumor"
export CELLTYPE_COLUMN="cell_type"

# Submit array job -- set the range to match the number of target cell types
# (e.g., if you have 12 non-query cell types, use --array=1-12)
sbatch --array=1-12 scripts/run_hymy_distance_correlation_v2.sh
```

Inside the SLURM script, `CELLTYPE_INDEX` is set from `$SLURM_ARRAY_TASK_ID` so each job processes a different target cell type.

For GPU permutation (Stage 2), submit with a GPU partition:

```bash
sbatch scripts/run_permutation_gpu.sh
```

Configure HPC-specific settings via environment variables:

```bash
export CONDA_SETUP="/path/to/conda.sh"
export CONDA_ENV="my_R_env"
export CUDA_MODULE="cuda/12.2"
```

---

## Advanced Examples

### Condition filtering

Analyze only samples from a specific condition:

```bash
export INPUT_PATH="/path/to/seurat.rds"
export QUERY_CELLTYPE="Tumor"
export CELLTYPE_COLUMN="cell_type"
export CONDITION_COLUMN="treatment"
export CONDITION_VALUE="treated"

Rscript scripts/hymy_distance_correlation_v2.R
```

### Custom target cell types

Restrict analysis to specific target cell types instead of auto-detecting all:

```bash
export INPUT_PATH="/path/to/seurat.rds"
export QUERY_CELLTYPE="Tumor"
export CELLTYPE_COLUMN="cell_type"
export TARGET_CELLTYPES="CD8_T,Macrophage,CAF,Endothelial"

Rscript scripts/hymy_distance_correlation_v2.R
```

### Stage 4 with custom control

Isolate tumor-specific effects from shared niche effects driven by CAFs:

```bash
export INPUT_PATH="/path/to/seurat.rds"
export QUERY_CELLTYPE="Tumor"
export CELLTYPE_COLUMN="cell_type"
export CONTROL_CELLTYPE="CAF"

Rscript scripts/hymy_distance_correlation_stage2.R
```

### Full publication pipeline

End-to-end analysis with all stages, permutation testing, confounder control, and visualization:

```bash
export INPUT_PATH="/path/to/seurat.rds"
export QUERY_CELLTYPE="Tumor"
export CELLTYPE_COLUMN="cell_type"
export SAMPLE_COLUMN="patient_id"
export OUTPUT_DIR="./results"
export QUERY_LABEL="Tumor"
export QUERY_SIGNATURE_GENES="EPCAM,KRT8,KRT18"
export CONTROL_CELLTYPE="CAF"
export ADATA_PATH="/path/to/adata.h5ad"
export ANALYSIS_NAME="tumor_gradient_analysis"

# Stage 1
Rscript scripts/hymy_distance_correlation_v2.R

# Stage 2
python scripts/run_permutation_gpu.py

# Stage 3
Rscript scripts/merge_permutation_pvals.R
Rscript scripts/recompute_meta_summary.R
Rscript scripts/merge_distance_correlation_results.R

# Stage 4
Rscript scripts/hymy_distance_correlation_stage2.R

# Stage 5
Rscript scripts/distance_correlation_atlas.R
Rscript scripts/hypothesis_visualizations.R
Rscript scripts/plot_decay_curves.R

# Stage 6
Rscript scripts/gradient_lr_integration.R
Rscript scripts/gradient_lr_atlas.R
Rscript scripts/gradient_lr_biology_figure.R
```

---

## Input Data Format

RIPPLE expects a Seurat object (`.rds`) specified by `INPUT_PATH` with the following structure:

| Component | Description |
|-----------|-------------|
| **Raw counts** | Unnormalized counts in the `RNA` assay, `counts` slot. Essential -- the Poisson GLM with cell-size offset handles normalization internally. Using pre-normalized data will produce incorrect results. |
| **Spatial coordinates** | X and Y coordinate columns in cell metadata. Column names are auto-detected (common names like `x_centroid`/`y_centroid`, `spatial_x`/`spatial_y`, `x`/`y` are recognized). Override with `X_COLUMN` and `Y_COLUMN` if needed. |
| **Cell type annotations** | A metadata column whose name matches `CELLTYPE_COLUMN`. |
| **Sample/replicate ID** | A metadata column for biological replicate identity. Defaults to `sample_id`; override with `SAMPLE_COLUMN`. |
| **Condition/group** | Optional. A metadata column for filtering to a subset of samples. Set `CONDITION_COLUMN` and `CONDITION_VALUE` to use. |

---

## Output Structure

Results are organized under the analysis output directory:

```
{OUTPUT_DIR}/spatial_analysis/{ANALYSIS_NAME}/
  per_celltype/
    {CellType}/
      meta_analysis_results.csv    # Per-gene summary with coefficients, p-values, FDR
      per_sample_results.csv       # Per-sample GLM coefficients
      permutation_pvals.csv        # Permutation test results (if run)
  summary/
    all_celltypes_results.csv      # Merged results across all cell types
    significant_genes_summary.csv  # Filtered to significant hits
  plots/
    multi_volcano.pdf              # Multi-panel volcano plot
    dotplot_specific_genes.pdf     # Cell-type-specific gene dotplot
    specific_genes_heatmap.pdf     # Heatmap of top hits
    fgsea_hallmark_*.pdf           # Pathway enrichment results
    forest_plots/                  # Per-gene forest plots
```

---

## Statistical Methodology

RIPPLE uses a per-sample Poisson GLM framework designed for count-based spatial transcriptomics data:

1. **Per-sample fitting** avoids pseudoreplication -- the true N is the number of biological replicates, not the number of cells
2. **Poisson GLM with offset** models raw counts directly, using `log(total_counts)` as an offset to account for cell size, ambient RNA, and segmentation differences
3. **Fisher's combined p-value** aggregates per-sample evidence with equal weight per replicate
4. **Sign consistency gate** requires all replicates to agree on the direction of the effect, ensuring reproducibility
5. **Two-tier expression filtering** applies strict thresholds for regular genes and lenient thresholds for curated priority genes (e.g., cytokines, chemokines)
6. **Permutation testing** (Stage 2) validates specificity by shuffling query cell labels
7. **Confounder control** (Stage 4) uses a bivariate model to separate query-specific effects from shared spatial niche effects

For a complete description of the statistical methodology, parameter choices, and interpretation guidelines, see `RIPPLE_user_guide.html`.

---

## Citation

If you use RIPPLE in your research, please cite:

> *Citation forthcoming.*

---

## License

*License forthcoming.*
