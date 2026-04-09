# RIPPLE

**Replicate-Aware Inference of Paracrine Profiles via Likelihood Estimation**

RIPPLE is an R package that detects distance-dependent gene expression gradients from any query cell type in spatial transcriptomics data.

---

## Overview

RIPPLE tries to answer a common question in spatial biology: *"Which genes in cell type B change expression as a function of physical distance from cell type A?"*

Given a spatial transcriptomics dataset with cell type annotations and biological replicates, RIPPLE fits per-sample Poisson GLMs to model gene expression as a function of distance to a query cell population, then combines evidence across replicates using Fisher's combined p-value with a sign consistency gate. The result is a ranked list of genes with distance-dependent expression gradients, along with effect sizes, significance measures, and confounder controls.

For the standalone script version (environment-variable-driven, SLURM-compatible), see [HyMy-distance-correlation-analysis](https://github.com/CMangana/HyMy-distance-correlation-analysis).

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

```r
# Install from GitHub
devtools::install_github("CMangana/RIPPLE")
```

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

```r
library(ripple)
results <- run_ripple(
  input           = "my_seurat.rds",
  query_celltype  = "Tumor",
  celltype_column = "cell_type",
  output_dir      = "./results"
)
```

This will auto-detect spatial coordinates, sample IDs, and all non-query cell types, then fit per-sample Poisson GLMs for every gene in every target cell type.

### Try it on the bundled mock dataset

RIPPLE ships with `ripple_mock_data`, a small synthetic SpatialExperiment (50 genes x 600 cells, 3 samples) containing a planted distance-dependent gradient. You can run the full pipeline on it without any external data:

```r
library(ripple)
data(ripple_mock_data)

results <- run_ripple(
  input           = ripple_mock_data,
  query_celltype  = "Tumor",
  celltype_column = "cell_type",
  output_dir      = tempfile("ripple_demo_"),
  min_expr_pct    = 0,
  min_expr_floor  = 10
)

# The top hits should be the planted INDUCED_* and REPRESSED_* genes in T_cell
head(results[order(fisher_fdr)], 12)
```

The 5 INDUCED genes have expression that decays away from Tumor cells (negative `median_coef`), and the 5 REPRESSED genes show the opposite pattern (positive `median_coef`). The 40 background genes should mostly be non-significant.

---

## Configuration

All configuration can be passed as function arguments or via environment variables (read by `ripple_config()` and `.onLoad()`). Function arguments take precedence over environment variables.

### Required

| Argument / Env Var | Description |
|--------------------|-------------|
| `input` / `INPUT_PATH` | Seurat, SingleCellExperiment, or SpatialExperiment object (or path to an `.rds` containing one). Must have raw counts. |
| `query_celltype` / `QUERY_CELLTYPE` | Cell type label for the source population (e.g., `"Tumor"`) |
| `celltype_column` / `CELLTYPE_COLUMN` | Metadata column containing cell type annotations |

### Data Configuration (Optional)

| Argument / Env Var | Default | Description |
|--------------------|---------|-------------|
| `output_dir` / `OUTPUT_DIR` | `./results` | Output directory for all results |
| `sample_column` / `SAMPLE_COLUMN` | `sample_id` | Metadata column containing sample/replicate IDs |
| `condition_column` / `CONDITION_COLUMN` | -- | Metadata column for condition/group filtering (e.g., `"treatment"`) |
| `condition_value` / `CONDITION_VALUE` | -- | Which condition to analyze (e.g., `"treated"`). If unset, all samples are used |
| `x_column` / `X_COLUMN` | auto-detect | Metadata column for X spatial coordinates |
| `y_column` / `Y_COLUMN` | auto-detect | Metadata column for Y spatial coordinates |
| `target_celltypes` / `TARGET_CELLTYPES` | auto-detect all | Target cell types to analyze (character vector or comma-separated string) |

### Analysis Configuration (Optional)

| Argument / Env Var | Default | Description |
|--------------------|---------|-------------|
| `query_label` / `QUERY_LABEL` | Same as `query_celltype` | Display name for the query cell type in plot titles |
| `query_signature_genes` / `QUERY_SIGNATURE_GENES` | -- | Marker genes for contamination check |
| `control_celltype` / `CONTROL_CELLTYPE` | -- | Control cell type for Stage 4 confounder analysis |
| `analysis_name` / `ANALYSIS_NAME` | `hymy_distance_correlation_v2` | Subdirectory name for output organization |
| `k_neighbors` / `K_NEIGHBORS` | `1` | Number of nearest query cells for distance calculation |
| `n_permutations` / `N_PERMUTATIONS` | `0` | Number of label permutations (0 = skip; use GPU Stage 2 instead for production) |

---

## Pipeline Stages

### Stage 1 (Core): Distance Correlation

**Function:** `run_ripple()`

For each gene in each non-query cell type, fits a per-sample Poisson GLM:

```r
glm(counts ~ distance_to_query + offset(log(total_counts)), family = poisson)
```

- `distance_to_query`: Euclidean distance (um) to nearest query cell (k-NN, default k=1)
- `offset(log(total_counts))`: cell-size correction, converts raw counts to rates
- **Coefficient (beta)**: log-rate change per um. Negative = expression increases near query cells (induced). Positive = expression decreases near query (repressed).

Per-sample coefficients are combined via Fisher's combined p-value. A sign consistency gate requires all replicates to agree on the direction of the effect.

### Stage 2 (Optional): GPU Permutation Testing

**Function:** `run_permutation_test()` (R) or GPU script via `system.file("python", "run_permutation_gpu.py", package = "ripple")`

GPU-accelerated null distribution via label permutation. Validates that the observed distance-expression relationship is specific to the query cell type and not a spatial artifact. Shuffles query cell labels across all cells and re-fits the GLM to build a null distribution.

### Stage 3 (Core): Merge and Summarize

**Functions:** `merge_ripple_results()`, `compute_fisher_pval()`

Integrates GPU permutation p-values into the R results, recomputes Fisher's combined p-values, and merges per-cell-type results into summary tables.

### Stage 4 (Optional): Confounder Control

**Function:** `run_ripple_confounder()`

Requires `control_celltype` to be set. Bivariate GLM adding distance to a control cell type as a covariate. Isolates query-specific effects from shared spatial niche effects. Classifies genes as:

| Classification | Criteria |
|----------------|----------|
| `query_specific` | Significant in bivariate model, same sign as Stage 1 |
| `enhanced` | Query-specific and effect size increased >10% vs Stage 1 |
| `niche_driven` | Not significant and coefficient attenuated >50% |
| `underpowered` | Not significant but coefficient preserved (>=50%) |

### Stage 5 (Optional): Visualization

**Functions:** `run_ripple_atlas()`, `plot_gradient_volcano()`, `plot_decay_curve()`

Produces publication-quality figures: multi-panel volcanos, dotplots of cell-type-specific genes, heatmaps, contamination flagging, fGSEA pathway enrichment (Hallmark gene sets), per-sample decay curves with reproducibility ribbons, and diagnostic panels.

### Stage 6 (Optional): Ligand-Receptor Integration

**Functions:** `run_ripple_lr()`, `classify_lr_artifacts()`

Matches gradient genes to ligand-receptor pairs using NicheNet. Runs three analyses per cell type: (1) direct L-R mapping of gradient receptors to query-expressed ligands, (2) NicheNet ligand activity prediction, and (3) Fisher's exact enrichment of NicheNet target predictions among gradient genes. Includes a 4-tier artifact classification system.

---

## Tutorial: Step-by-Step Analysis

This tutorial walks through a complete RIPPLE analysis using the R package functions. We assume you have a Seurat object with Xenium data, cell type annotations, and at least 3 biological replicates.

### Step 0: Check the binary assumption

RIPPLE uses a Poisson GLM on raw transcript counts. This works well when most detected genes have 1 transcript per cell (the typical regime for Xenium/MERFISH). Before running the analysis, verify this assumption by plotting the relationship between number of detected features and total counts per cell.

```r
library(Seurat)
library(ggplot2)

obj <- readRDS("my_seurat.rds")
meta <- obj@meta.data
meta$nCount <- colSums(obj[["RNA"]]$counts)
meta$nFeature <- colSums(obj[["RNA"]]$counts > 0)

ggplot(meta, aes(x = nFeature, y = nCount)) +
  geom_point(size = 0.1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Binary assumption check",
    subtitle = "Points near the red line = most genes have 1 transcript (Poisson appropriate)",
    x = "Number of detected genes (nFeature)",
    y = "Total transcript counts (nCount)"
  ) +
  theme_minimal()
```

If points cluster tightly along the `y = x` line, the Poisson model with a cell-size offset is appropriate. If there is substantial spread above the line (many genes with 2+ transcripts per cell), consider a negative binomial model instead.

### Step 1: Load data and inspect

```r
library(ripple)

# Quick check that the data has the expected columns
check_data("my_seurat.rds", celltype_column = "cell_type", sample_column = "sample_id")

# Load metadata only (fast, no expression matrix in memory)
meta <- load_metadata_only(
  "my_seurat.rds",
  celltype_column = "cell_type",
  sample_column = "sample_id"
)

# Verify query cell type is present and has enough cells
table(meta[["cell_type"]])
```

### Step 2: Run the core analysis

```r
results <- run_ripple(
  input            = "my_seurat.rds",
  query_celltype   = "Tumor",
  celltype_column  = "cell_type",
  output_dir       = "./results",
  sample_column    = "sample_id",
  analysis_name    = "tumor_ripple"
)

# Results is a data.table with one row per gene x cell type
head(results[order(fisher_fdr)])
```

`run_ripple()` auto-detects all non-query cell types as targets, fits per-sample Poisson GLMs, combines across replicates with Fisher's method, and writes per-cell-type CSVs to `./results/spatial_analysis_Tumor/tumor_ripple/`.

To restrict to a specific condition (e.g., only treated samples):

```r
results <- run_ripple(
  input            = "my_seurat.rds",
  query_celltype   = "Tumor",
  celltype_column  = "cell_type",
  output_dir       = "./results",
  condition_column = "treatment",
  condition_value  = "treated"
)
```

### Step 3: Merge results across cell types

```r
merged <- merge_ripple_results(
  results_dir = "./results/spatial_analysis_Tumor/tumor_ripple",
  recompute_fisher = TRUE
)

# Top hits
merged[fisher_fdr < 0.05][order(fisher_fdr)]
```

### Step 4: Visualize the results

#### Volcano plot

```r
library(data.table)

# Load merged results
all_results <- fread("./results/spatial_analysis_Tumor/tumor_ripple/summary/all_genes_results.csv")

# Volcano for one cell type
cd8_results <- all_results[cell_type == "CD8_T"]
plot_gradient_volcano(
  cd8_results,
  query_label = "Tumor",
  title = "CD8 T cells: distance-dependent genes"
)
```

#### Decay curves

```r
# Bin expression by distance for a specific gene
obj <- readRDS("my_seurat.rds")
counts_matrix <- obj[["RNA"]]$counts

# Get the gene of interest
gene <- "PDCD1"
gene_counts <- as.numeric(counts_matrix[gene, ])

# You need the distances (saved during run_ripple in per_celltype output)
# Or compute them directly:
meta <- as.data.table(obj@meta.data, keep.rownames = "barcode")
coord_cols <- get_coord_columns(meta)
query_coords <- as.matrix(meta[cell_type == "Tumor", ..coord_cols])
target_mask <- meta$cell_type == "CD8_T"
target_coords <- as.matrix(meta[target_mask, ..coord_cols])

nn <- RANN::nn2(query_coords, target_coords, k = 1)
distances <- as.vector(nn$nn.dists)

# Bin and plot
bins <- bin_decay_data(
  gene_counts[target_mask],
  distances,
  n_bins = 20,
  max_distance = 200
)
plot_decay_curve(bins, gene_name = "PDCD1", cell_type = "CD8_T")
```

#### Gene specificity classification

```r
# Which genes are specific to one cell type vs contamination?
specificity <- classify_gene_specificity(all_results, fdr_threshold = 0.05)
table(specificity$specificity_class)
#   specific  moderate  ubiquitous  contamination
#       412       187          45             23
```

### Step 5: Confounder control (optional)

If your query cell type co-localizes with another cell type, validate that the gradients are query-specific:

```r
stage2 <- run_ripple_confounder(
  input            = "my_seurat.rds",
  results_dir      = "./results/spatial_analysis_Tumor/tumor_ripple",
  query_celltype   = "Tumor",
  celltype_column  = "cell_type",
  control_celltype = "CAF"   # cells that share the same niche
)

# Classification breakdown
table(stage2$classification)
#   query_specific  enhanced  niche_driven  underpowered
#             892       134            67            45
```

### Step 6: Generate atlas figures (optional)

```r
run_ripple_atlas(
  results_dir = "./results/spatial_analysis_Tumor/tumor_ripple",
  output_dir  = "./results/spatial_analysis_Tumor/tumor_ripple/atlas",
  query_label = "Tumor",
  run_fgsea   = TRUE,
  organism    = "human"
)
```

This generates volcanos, dotplots, heatmaps, fGSEA enrichment panels, and contamination reports.

---

## Tutorial: Downstream Integrations

After running the core RIPPLE analysis, you can feed the results into pathway enrichment and ligand-receptor analysis.

### Pathway enrichment with fGSEA

RIPPLE provides `run_ripple_fgsea()` which ranks genes by their gradient coefficient and runs fast gene set enrichment analysis per cell type.

```r
library(ripple)
library(data.table)

# Load merged results
results <- fread("./results/spatial_analysis_Tumor/tumor_ripple/summary/all_genes_results.csv")

# Run fGSEA with Hallmark gene sets
gsea <- run_ripple_fgsea(
  results,
  gene_sets = "hallmark",    # or "kegg", "reactome", "go_bp"
  organism  = "human",       # or "mouse"
  coef_col  = "median_coef",
  fdr_col   = "fisher_fdr"
)

# Top enriched pathways across all cell types
gsea[padj < 0.05][order(pval)]

# Filter to one cell type
gsea[cell_type == "CD8_T" & padj < 0.05][order(NES)]
```

You can also pass your own custom gene sets:

```r
my_gene_sets <- list(
  exhaustion = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "ENTPD1"),
  cytotoxicity = c("GZMA", "GZMB", "PRF1", "GNLY", "NKG7"),
  stemness = c("TCF7", "SELL", "IL7R", "LEF1", "CCR7")
)

gsea_custom <- run_ripple_fgsea(results, gene_sets = my_gene_sets)
```

### Ligand-receptor integration with NicheNet

RIPPLE provides `run_ripple_lr()` which matches gradient genes to ligand-receptor pairs using the NicheNet database. It scores L-R pairs by combining direct expression matching, NicheNet activity prediction, and downstream target enrichment.

```r
# Requires nichenetr: devtools::install_github("saeyslab/nichenetr")

lr_results <- run_ripple_lr(
  results_dir      = "./results/spatial_analysis_Tumor/tumor_ripple",
  input            = "my_seurat.rds",
  query_celltype   = "Tumor",
  celltype_column  = "cell_type",
  organism         = "human",
  output_dir       = "./results/lr_integration"
)

# Top L-R pairs (Direction A: Query -> Target)
lr_results[artifact_flag == "clean"][order(-combined_score)][1:20]
```

#### Artifact classification

L-R results can include false positives from segmentation artifacts. Use `classify_lr_artifacts()` to flag them:

```r
# Provide query marker genes to detect leakage
lr_clean <- classify_lr_artifacts(
  lr_results,
  query_signature = c("EPCAM", "KRT8", "KRT18"),
  contamination_genes = c("JCHAIN", "IGKC"),  # genes in many cell types
  low_expr_threshold = 0.02
)

# See classification breakdown
table(lr_clean$artifact_flag)
#   clean  suspect  low_confidence  artifact
#     891       42              31        18

# Use only clean pairs for downstream analysis
lr_clean[artifact_flag == "clean"][order(-combined_score)]
```

#### Combining with Stage 4 results

If you ran the confounder control, you can cross-reference to focus on query-specific L-R pairs:

```r
# Load Stage 4 classifications
stage2 <- fread("./results/spatial_analysis_Tumor/tumor_ripple_stage2/summary/stage2_all_results.csv")

# Keep only L-R pairs where the receptor gene is query-specific (not niche-driven)
query_specific_genes <- stage2[classification == "query_specific"]$gene
lr_validated <- lr_clean[
  artifact_flag == "clean" &
  receptor %in% query_specific_genes
][order(-combined_score)]

# These are high-confidence L-R pairs: clean artifacts + query-specific gradients
lr_validated[1:10, .(ligand, receptor, cell_type, combined_score)]
```

---

## Input Data Format

RIPPLE expects a Seurat object (`.rds`) with the following structure:

| Component | Description |
|-----------|-------------|
| **Raw counts** | Unnormalized counts in the `RNA` assay, `counts` slot. Essential -- the Poisson GLM with cell-size offset handles normalization internally. Using pre-normalized data will produce incorrect results. |
| **Spatial coordinates** | X and Y coordinate columns in cell metadata. Column names are auto-detected (common names like `x_centroid`/`y_centroid`, `spatial_x`/`spatial_y`, `x`/`y` are recognized). Override with `x_column` and `y_column` arguments if needed. |
| **Cell type annotations** | A metadata column whose name matches `celltype_column`. |
| **Sample/replicate ID** | A metadata column for biological replicate identity. Defaults to `sample_id`; override with `sample_column`. |
| **Condition/group** | Optional. A metadata column for filtering to a subset of samples. Set `condition_column` and `condition_value` to use. |

---

## Output Structure

Results are organized under the analysis output directory:

```
{output_dir}/spatial_analysis/{analysis_name}/
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

---

## Citation

If you use RIPPLE in your research, please cite:

> *Citation forthcoming.*

---

## License

*License forthcoming.*
