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

## Quick Start

```bash
# Set your query cell type and the metadata column containing cell type labels
export QUERY_CELLTYPE="Tumor"
export CELLTYPE_COLUMN="cell_type"

# Stage 1: Distance correlation (skip permutations for a fast first pass)
N_PERMUTATIONS=0 Rscript scripts/hymy_distance_correlation_v2.R
```

This will fit per-sample Poisson GLMs for every gene in every non-query cell type, testing whether expression changes as a function of distance to the nearest query cell.

---

## Configuration

All configuration is via environment variables. Set them before running any script.

| Variable | Required | Default | Description |
|----------|----------|---------|-------------|
| `QUERY_CELLTYPE` | Yes | -- | Cell type label for the source population (e.g., `"Tumor"`, `"HyMy_GMM"`) |
| `CELLTYPE_COLUMN` | Yes | -- | Metadata column containing cell type annotations |
| `QUERY_LABEL` | No | Same as `QUERY_CELLTYPE` | Short display name used in plot titles |
| `QUERY_SIGNATURE_GENES` | No | Empty | Comma-separated marker genes for contamination check |
| `CELLTYPE_INDEX` | No | -- | 1-based index for SLURM array jobs (run one cell type per job) |
| `K_NEIGHBORS` | No | `1` | Number of nearest query cells for distance calculation |
| `N_PERMUTATIONS` | No | `500` | Number of label permutations (set 0 to skip, 1000+ for publication) |
| `ANALYSIS_NAME` | No | `hymy_distance_correlation_v2` | Analysis subdirectory name for output |
| `ANNOTATION_LEVEL` | No | -- | Legacy shorthand: `"HyMy"` or `"L1"` (sets query type and column automatically) |
| `GRADIENT_SOURCE` | No | `hymy_distance_correlation` | Which Stage 1 results to use for downstream stages |

---

## Pipeline Stages

### Stage 1: Distance Correlation

**Script:** `hymy_distance_correlation_v2.R`

For each gene in each non-query cell type, fits a per-sample Poisson GLM:

```r
glm(counts ~ distance_to_query + offset(log(total_counts)), family = poisson)
```

- `distance_to_query`: Euclidean distance (um) to nearest query cell (k-NN, default k=1)
- `offset(log(total_counts))`: cell-size correction, converts raw counts to rates
- **Coefficient (beta)**: log-rate change per um. Negative = expression increases near query cells (induced). Positive = expression decreases near query (repressed).

Per-sample coefficients are combined via Fisher's combined p-value. A sign consistency gate requires all replicates to agree on the direction of the effect.

### Stage 2: GPU Permutation Testing

**Script:** `run_permutation_gpu.py`

GPU-accelerated null distribution via label permutation. Validates that the observed distance-expression relationship is specific to the query cell type and not a spatial artifact. Shuffles query cell labels across all cells and re-fits the GLM to build a null distribution.

### Stage 3: Merge and Summarize

**Scripts:** `merge_permutation_pvals.R`, `recompute_meta_summary.R`, `merge_distance_correlation_results.R`

Integrates GPU permutation p-values into the R results, recomputes Fisher's combined p-values, and merges per-cell-type results into summary tables.

### Stage 4: Confounder Control

**Script:** `hymy_distance_correlation_stage2.R`

Bivariate GLM adding distance to a control cell type as a covariate. Isolates query-specific effects from shared spatial niche effects. Classifies genes as:

| Classification | Criteria |
|----------------|----------|
| `query_specific` | Significant in bivariate model, same sign as Stage 1 |
| `enhanced` | Query-specific and effect size increased >10% vs Stage 1 |
| `niche_driven` | Not significant and coefficient attenuated >50% |
| `underpowered` | Not significant but coefficient preserved (>=50%) |

### Stage 5: Visualization

**Scripts:** `distance_correlation_atlas.R`, `hypothesis_visualizations.R`

Produces publication-quality figures: multi-panel volcanos, dotplots of cell-type-specific genes, heatmaps, contamination flagging, fGSEA pathway enrichment (Hallmark gene sets), per-sample decay curves with reproducibility ribbons, and diagnostic panels.

### Stage 6: Ligand-Receptor Integration

**Scripts:** `gradient_lr_integration.R`, `gradient_lr_atlas.R`

Matches gradient genes to ligand-receptor pairs using NicheNet. Runs three analyses per cell type: (1) direct L-R mapping of gradient receptors to query-expressed ligands, (2) NicheNet ligand activity prediction, and (3) Fisher's exact enrichment of NicheNet target predictions among gradient genes. Includes a 4-tier artifact classification system.

---

## Full Pipeline Example

```bash
export QUERY_CELLTYPE="Tumor"
export CELLTYPE_COLUMN="cell_type"
export QUERY_LABEL="Tumor"

# Stage 1: Run distance correlation (skip permutations for speed)
N_PERMUTATIONS=0 Rscript scripts/hymy_distance_correlation_v2.R

# Stage 2: GPU permutation testing
python scripts/run_permutation_gpu.py

# Stage 3: Merge results
Rscript scripts/merge_permutation_pvals.R
Rscript scripts/recompute_meta_summary.R
Rscript scripts/merge_distance_correlation_results.R

# Stage 4: Confounder control (optional -- requires a suitable control cell type)
Rscript scripts/hymy_distance_correlation_stage2.R

# Stage 5: Visualization
Rscript scripts/distance_correlation_atlas.R

# Stage 6: L-R integration (requires NicheNet)
Rscript scripts/gradient_lr_integration.R
Rscript scripts/gradient_lr_atlas.R
```

---

## Input Data Format

RIPPLE expects a Seurat object (`.rds`) with the following structure:

| Component | Description |
|-----------|-------------|
| **Counts** | Raw (unnormalized) counts in the `RNA` assay, `counts` layer |
| **Spatial coordinates** | `spatial_x` and `spatial_y` (or `x` and `y`) columns in cell metadata |
| **Cell type annotations** | A metadata column matching the value of `CELLTYPE_COLUMN` |
| **Sample/replicate ID** | `sample_id` column in metadata |
| **Condition/group** | `condition` or `group` column in metadata (optional, for condition-specific analyses) |

Raw counts are essential -- the Poisson GLM with cell-size offset handles normalization internally. Using pre-normalized data will produce incorrect results.

---

## Output Structure

Results are organized under the analysis output directory:

```
results/spatial_analysis/{ANALYSIS_NAME}/
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
