#!/usr/bin/env python3
"""
RIPPLE Stage 2: GPU-Accelerated Permutation Testing
====================================================

Replaces the CPU-bound permutation step of the distance correlation R scripts
with GPU-accelerated kNN using PyTorch CUDA tensors.

Supports both v1 (logistic) and v2 (Poisson GLM with cell size offset).
Model selection is automatic based on ANALYSIS_NAME:
  - v1 (default): binary logistic regression on expression detection
  - v2 (ANALYSIS_NAME contains "v2"): Poisson GLM on raw counts with offset

Reads:
  - h5ad file with expression, coordinates, metadata
  - meta_analysis_results.csv from the R script (run with N_PERMUTATIONS=0)

Writes:
  - permutation_pvals.csv per cell type (same format as R script)

Usage:
  python run_permutation_gpu.py --celltype LEC
  python run_permutation_gpu.py --celltype CD8_T_cells --annotation-level L1
  CELLTYPE_INDEX=5 python run_permutation_gpu.py  # Array job mode
  QUERY_CELLTYPE=MyType CELLTYPE_COLUMN=my_col python run_permutation_gpu.py

Requirements: PyTorch with CUDA, anndata, numpy, pandas, scipy
Environment: harpy (/nobackup/lab_maier/envs/harpy)
"""

import os
import sys
import time
import argparse
import numpy as np
import pandas as pd
import torch
from scipy import sparse
from pathlib import Path

# PyTorch 2.1 is ABI-incompatible with NumPy 2.x. Fail early with an actionable
# message instead of a cryptic C-level ABI error deep inside a tensor op.
if int(np.__version__.split(".")[0]) >= 2:
    sys.exit(
        f"NumPy {np.__version__} detected, but this script requires numpy<2 "
        "(PyTorch 2.1 is not ABI-compatible with NumPy 2.x). "
        "Install a compatible version with: pip install 'numpy<2'"
    )

# =============================================================================
# Configuration
# =============================================================================

# Env var resolution (matches config.R)
INPUT_PATH = os.environ.get("INPUT_PATH", "")
ADATA_PATH_ENV = os.environ.get("ADATA_PATH", "")
OUTPUT_DIR = os.environ.get("OUTPUT_DIR", "")
SAMPLE_COL = os.environ.get("SAMPLE_COLUMN", "sample_id")
CONDITION_COL = os.environ.get("CONDITION_COLUMN", "")
CONDITION_VAL = os.environ.get("CONDITION_VALUE", "")

# Resolve ADATA_PATH and PROJECT_ROOT
if ADATA_PATH_ENV:
    ADATA_PATH = Path(ADATA_PATH_ENV)
elif INPUT_PATH:
    # Try to find .h5ad alongside the .rds
    rds_dir = Path(INPUT_PATH).parent
    candidates = list(rds_dir.glob("*.h5ad"))
    ADATA_PATH = candidates[0] if candidates else None
else:
    ADATA_PATH = None  # Will be set from legacy paths below

# Platform detection (legacy — only needed when INPUT_PATH is not set)
if not INPUT_PATH:
    if sys.platform == "win32":
        BASE_PATH = Path("N:/lab_maier/Projects/mXenium")
    else:
        BASE_PATH = Path("/nobackup/lab_maier/Projects/mXenium")
    PROJECT_ROOT = BASE_PATH / "CMM"
    if ADATA_PATH is None:
        ADATA_PATH = BASE_PATH / "results" / "cell_type_assignment" / "adata_incl_hymy.h5ad"
else:
    PROJECT_ROOT = Path(OUTPUT_DIR).parent if OUTPUT_DIR else Path.cwd()

# Statistical parameters (match R script)
MIN_CELLS_PER_SAMPLE = 30
MAX_DISTANCE_UM = 200.0
N_PERMUTATIONS = 500
PERM_TOP_N = 100
FDR_THRESHOLD = 0.05

# Priority genes for permutation testing (same curated list as R script)
PERM_PRIORITY_GENES = [
    # CC Chemokines
    "Ccl1", "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl6", "Ccl7", "Ccl8", "Ccl9",
    "Ccl11", "Ccl12", "Ccl17", "Ccl19", "Ccl20", "Ccl21a", "Ccl21b", "Ccl21c",
    "Ccl22", "Ccl24", "Ccl25", "Ccl27a", "Ccl27b", "Ccl28",
    # CXC Chemokines
    "Cxcl1", "Cxcl2", "Cxcl3", "Cxcl4", "Cxcl5", "Cxcl7", "Cxcl9", "Cxcl10",
    "Cxcl11", "Cxcl12", "Cxcl13", "Cxcl14", "Cxcl15", "Cxcl16", "Cxcl17",
    # CX3C / XC
    "Cx3cl1", "Xcl1",
    # CC Receptors
    "Ccr1", "Ccr2", "Ccr3", "Ccr4", "Ccr5", "Ccr6", "Ccr7", "Ccr8", "Ccr9", "Ccr10",
    # CXC Receptors
    "Cxcr1", "Cxcr2", "Cxcr3", "Cxcr4", "Cxcr5", "Cxcr6",
    # CX3C / XC / Atypical
    "Cx3cr1", "Xcr1", "Ackr1", "Ackr2", "Ackr3", "Ackr4",
    # Interleukins
    "Il1a", "Il1b", "Il1rn", "Il2", "Il3", "Il4", "Il5", "Il6", "Il7", "Il9",
    "Il10", "Il11", "Il12a", "Il12b", "Il13", "Il14", "Il15", "Il16",
    "Il17a", "Il17b", "Il17c", "Il17d", "Il17f",
    "Il18", "Il19", "Il20", "Il21", "Il22", "Il23a", "Il24", "Il25",
    "Il27", "Il31", "Il33", "Il34",
    "Il36a", "Il36b", "Il36g", "Il36rn",
    "Il1f5", "Il1f9", "Il1f10",
    # IL Receptors
    "Il1r1", "Il1r2", "Il1rl1", "Il1rl2", "Il1rap",
    "Il2ra", "Il2rb", "Il2rg",
    "Il3ra", "Il4ra", "Il5ra", "Il6ra", "Il6st", "Il7r", "Il9r",
    "Il10ra", "Il10rb", "Il11ra1", "Il12rb1", "Il12rb2",
    "Il13ra1", "Il13ra2", "Il15ra",
    "Il17ra", "Il17rb", "Il17rc", "Il17rd", "Il17re",
    "Il18r1", "Il18rap", "Il20ra", "Il20rb", "Il21r",
    "Il22ra1", "Il22ra2", "Il23r", "Il27ra", "Il31ra",
    # Interferons
    "Ifna1", "Ifna2", "Ifna4", "Ifna5", "Ifna6", "Ifna7", "Ifna9",
    "Ifna11", "Ifna12", "Ifna13", "Ifna14", "Ifnab",
    "Ifnb1", "Ifne", "Ifnk", "Ifng", "Ifnl2", "Ifnl3",
    # IFN Receptors
    "Ifnar1", "Ifnar2", "Ifngr1", "Ifngr2", "Ifnlr1",
    # TNF Superfamily
    "Tnf", "Lta", "Ltb", "Fasl", "Cd40lg",
    "Tnfsf4", "Tnfsf8", "Tnfsf9", "Tnfsf10", "Tnfsf11", "Tnfsf12",
    "Tnfsf13", "Tnfsf13b", "Tnfsf14", "Tnfsf15", "Tnfsf18",
    # TNF Receptors
    "Tnfrsf1a", "Tnfrsf1b", "Tnfrsf4", "Fas", "Cd40",
    "Tnfrsf8", "Tnfrsf9", "Tnfrsf10b", "Tnfrsf11a", "Tnfrsf11b",
    "Tnfrsf12a", "Tnfrsf13b", "Tnfrsf13c", "Tnfrsf14",
    "Tnfrsf17", "Tnfrsf18", "Tnfrsf19", "Tnfrsf21", "Tnfrsf25",
    # CSF
    "Csf1", "Csf2", "Csf3", "Csf1r", "Csf2ra", "Csf2rb", "Csf2rb2", "Csf3r",
    # TGF-beta
    "Tgfb1", "Tgfb2", "Tgfb3", "Tgfbr1", "Tgfbr2", "Tgfbr3",
    # VEGF
    "Vegfa", "Vegfb", "Vegfc", "Vegfd", "Flt1", "Kdr", "Flt4", "Nrp1", "Nrp2",
    # gp130
    "Lif", "Osm", "Cntf", "Lifr", "Osmr", "Cntfr",
    # Other
    "Tslp", "Mif", "Spp1", "Kitl", "Kit", "Epo", "Epor", "Thpo", "Mpl",
]

# Target cell type definitions (match R script)
TARGET_CELLTYPES = {
    "LEC": ["LEC"],
    "FRC": ["FRC"],
    "BEC": ["BEC"],
    "CD4_T_cells": ["Naive_CD4", "Tfh", "Treg"],
    "CD8_T_cells": ["Naive_CD8", "Activated_CD8", "Cytotoxic_CD8", "Tpex"],
    "gdT_cells": ["gdT_cell"],
    "Macrophages": ["Macrophages"],
    "Monocyte": ["Monocyte"],
    "Fibroblasts_mac": ["Fibroblasts_mac"],
    "cDC1": ["cDC1"],
    "cDC2": ["cDC2"],
    "mature_migDC": ["mature_migDC"],
    "B_cells": ["B_cell", "Follicular_B"],
    "Plasma_cell": ["Plasma_cell"],
}

# Cell type index mapping (for SLURM array jobs)
CELLTYPE_INDEX_MAP = {
    1: "LEC", 2: "FRC", 3: "BEC", 4: "CD4_T_cells", 5: "CD8_T_cells",
    6: "gdT_cells", 7: "Macrophages", 8: "Monocyte", 9: "Fibroblasts_mac",
    10: "cDC1", 11: "cDC2", 12: "mature_migDC",
    13: "B_cells", 14: "Plasma_cell",
}


# =============================================================================
# GPU kNN Function
# =============================================================================

def gpu_knn_distances(query_coords: torch.Tensor, target_coords: torch.Tensor,
                      k: int = 1, batch_size: int = 10000) -> torch.Tensor:
    """
    Compute k-nearest-neighbor distances from target to query points on GPU.

    Uses batched pairwise distance computation to avoid OOM on large datasets.
    For k=1, returns the minimum distance from each target point to any query point.

    Args:
        query_coords: (n_query, 2) tensor of query cell coordinates
        target_coords: (n_target, 2) tensor of target cell coordinates
        k: number of nearest neighbors (default 1)
        batch_size: target cells per batch to control GPU memory

    Returns:
        (n_target,) tensor of distances to nearest query cell
    """
    n_target = target_coords.shape[0]
    min_dists = torch.empty(n_target, device=target_coords.device)

    for start in range(0, n_target, batch_size):
        end = min(start + batch_size, n_target)
        batch = target_coords[start:end]  # (batch, 2)

        # Pairwise squared distances: ||a - b||^2 = ||a||^2 + ||b||^2 - 2*a.b
        dists_sq = torch.cdist(batch, query_coords, p=2.0)  # (batch, n_query)

        if k == 1:
            min_dists[start:end] = dists_sq.min(dim=1).values
        else:
            # Mean of k nearest distances (not just k-th)
            min_dists[start:end] = dists_sq.topk(k, dim=1, largest=False).values.mean(dim=1)

    return min_dists


# =============================================================================
# Logistic Regression (CPU, per-sample)
# =============================================================================

def fit_logistic_sample(expressing: np.ndarray, distances: np.ndarray,
                        min_cells: int = 5):
    """
    Fit logistic regression: P(expressing) ~ distance.
    Returns (coef, se) or (nan, nan) on failure.

    Uses statsmodels-free implementation via scipy for speed.
    Falls back to simple IRLS for logistic regression.
    """
    valid = np.isfinite(expressing) & np.isfinite(distances)
    y = expressing[valid].astype(np.float64)
    x = distances[valid].astype(np.float64)

    if len(y) < min_cells or y.var() == 0:
        return np.nan, np.nan

    # Logistic regression via IRLS (iteratively reweighted least squares)
    # Model: logit(p) = beta0 + beta1 * x
    n = len(y)
    X = np.column_stack([np.ones(n), x])

    # Initialize with OLS on logit scale (approximate)
    beta = np.zeros(2)

    for iteration in range(25):  # Max IRLS iterations
        eta = X @ beta
        # Clip to avoid overflow in exp
        eta = np.clip(eta, -20, 20)
        mu = 1.0 / (1.0 + np.exp(-eta))
        # Ensure mu is bounded away from 0/1
        mu = np.clip(mu, 1e-10, 1 - 1e-10)

        W = mu * (1 - mu)
        z = eta + (y - mu) / W  # Working response

        # Weighted least squares: (X'WX)^-1 X'Wz
        XtWX = X.T @ (X * W[:, None])
        XtWz = X.T @ (W * z)

        try:
            beta_new = np.linalg.solve(XtWX, XtWz)
        except np.linalg.LinAlgError:
            return np.nan, np.nan

        # Check convergence
        if np.max(np.abs(beta_new - beta)) < 1e-8:
            beta = beta_new
            break
        beta = beta_new
    else:
        # Did not converge
        return np.nan, np.nan

    # Standard error from Fisher information
    try:
        cov = np.linalg.inv(XtWX)
        se = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        return np.nan, np.nan

    if not np.isfinite(beta[1]) or not np.isfinite(se[1]) or se[1] <= 0:
        return np.nan, np.nan

    return beta[1], se[1]  # coefficient and SE for distance term


# =============================================================================
# Poisson Regression (CPU, per-sample) — v2
# =============================================================================

def fit_poisson_sample(counts: np.ndarray, distances: np.ndarray,
                       log_total_counts: np.ndarray, min_cells: int = 5):
    """
    Fit Poisson GLM: counts ~ distance + offset(log(total_counts)).
    Returns (coef, se) or (nan, nan) on failure.

    Uses custom IRLS for Poisson with log link.
    The offset enters through mu = exp(X@beta + offset) but the working
    response uses X@beta (without offset) for the WLS solve.
    """
    valid = np.isfinite(counts) & np.isfinite(distances) & np.isfinite(log_total_counts)
    y = counts[valid].astype(np.float64)
    x = distances[valid].astype(np.float64)
    offset = log_total_counts[valid].astype(np.float64)

    if len(y) < min_cells or y.sum() == 0:
        return np.nan, np.nan

    n = len(y)
    X = np.column_stack([np.ones(n), x])

    # Initialize beta from log(mean rate)
    mean_rate = np.maximum(y.sum() / np.exp(offset).sum(), 1e-10)
    beta = np.array([np.log(mean_rate), 0.0])

    for iteration in range(25):
        lin_pred = X @ beta  # Linear predictor WITHOUT offset
        eta = lin_pred + offset  # Full predictor WITH offset
        # Clip to avoid overflow
        eta = np.clip(eta, -20, 20)
        mu = np.exp(eta)
        mu = np.maximum(mu, 1e-10)

        # Poisson IRLS: weights = mu, working response = lin_pred + (y - mu) / mu
        W = mu
        z = lin_pred + (y - mu) / mu

        # Weighted least squares: (X'WX)^{-1} X'Wz
        XtWX = X.T @ (X * W[:, None])
        XtWz = X.T @ (W * z)

        try:
            beta_new = np.linalg.solve(XtWX, XtWz)
        except np.linalg.LinAlgError:
            return np.nan, np.nan

        if np.max(np.abs(beta_new - beta)) < 1e-8:
            beta = beta_new
            break
        beta = beta_new
    else:
        return np.nan, np.nan

    # Standard error from Fisher information
    try:
        cov = np.linalg.inv(XtWX)
        se = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        return np.nan, np.nan

    if not np.isfinite(beta[1]) or not np.isfinite(se[1]) or se[1] <= 0:
        return np.nan, np.nan

    return beta[1], se[1]  # coefficient and SE for distance term


# =============================================================================
# Permutation Test (GPU-accelerated)
# =============================================================================

def run_permutation_test_gpu(
    gene_data: np.ndarray,
    target_coords_gpu: torch.Tensor,
    all_coords_gpu: torch.Tensor,
    sample_ids_target: np.ndarray,
    sample_ids_all: np.ndarray,
    query_per_sample: dict,
    observed_coef: float,
    n_perms: int,
    device: torch.device,
    rng: np.random.Generator,
    k_neighbors: int = 1,
    use_poisson: bool = False,
    log_total_counts_target: np.ndarray = None,
) -> float:
    """
    GPU-accelerated permutation test for a single gene.

    For each permutation:
    1. Sample pseudo-query cells (stratified by sample) on CPU
    2. Compute kNN distances on GPU (the bottleneck that's now fast)
    3. Fit per-sample regression on CPU (logistic or Poisson)
    4. Combine via inverse-variance weighting

    Args:
        gene_data: Binary expressing (v1) or raw counts (v2) for target cells
        use_poisson: If True, use Poisson GLM with offset instead of logistic
        log_total_counts_target: Required if use_poisson=True; log(total_counts)
        k_neighbors: Number of nearest neighbors for distance computation

    Returns empirical two-sided p-value.
    """
    unique_samples = list(query_per_sample.keys())
    null_coefs = np.empty(n_perms)
    null_coefs[:] = np.nan

    # Pre-compute sample masks for all cells (for stratified sampling)
    sample_masks_all = {s: np.where(sample_ids_all == s)[0] for s in unique_samples}
    sample_masks_target = {s: np.where(sample_ids_target == s)[0] for s in unique_samples}

    target_coords_cpu = target_coords_gpu.cpu().numpy()

    for i in range(n_perms):
        # 1. Stratified sampling of pseudo-query cells
        pseudo_query_list = []
        for samp in unique_samples:
            samp_indices = sample_masks_all[samp]
            n_to_sample = query_per_sample[samp]
            if n_to_sample > 0 and len(samp_indices) >= n_to_sample:
                chosen = rng.choice(samp_indices, size=n_to_sample, replace=False)
                pseudo_query_list.append(chosen)

        if len(pseudo_query_list) == 0:
            continue
        pseudo_query_idx = np.concatenate(pseudo_query_list)
        if len(pseudo_query_idx) < 5:
            continue

        # 2. GPU kNN: distances from target cells to pseudo-query cells
        pseudo_query_coords = all_coords_gpu[pseudo_query_idx]
        perm_distances = gpu_knn_distances(pseudo_query_coords, target_coords_gpu,
                                           k=k_neighbors)
        perm_distances = torch.clamp(perm_distances, max=MAX_DISTANCE_UM)
        perm_distances_cpu = perm_distances.cpu().numpy()

        # 3. Per-sample regression
        coefs = []
        ses = []
        for samp in unique_samples:
            idx = sample_masks_target[samp]
            if len(idx) < MIN_CELLS_PER_SAMPLE:
                continue

            samp_data = gene_data[idx]
            samp_dist = perm_distances_cpu[idx]

            if use_poisson:
                if samp_data.sum() == 0:
                    continue
                samp_offset = log_total_counts_target[idx]
                c, s = fit_poisson_sample(samp_data, samp_dist, samp_offset)
            else:
                if samp_data.var() == 0:
                    continue
                c, s = fit_logistic_sample(samp_data, samp_dist)

            if np.isfinite(c) and np.isfinite(s) and s > 0:
                coefs.append(c)
                ses.append(s)

        # 4. Inverse-variance weighted mean
        if len(coefs) >= 2:
            coefs = np.array(coefs)
            ses = np.array(ses)
            weights = 1.0 / (ses ** 2)
            null_coefs[i] = np.sum(weights * coefs) / np.sum(weights)

    # Empirical two-sided p-value
    valid_null = null_coefs[np.isfinite(null_coefs)]
    if len(valid_null) < 10:
        return np.nan

    # Conservative: add 1 to numerator and denominator
    perm_pval = (np.sum(np.abs(valid_null) >= np.abs(observed_coef)) + 1) / (len(valid_null) + 1)
    return perm_pval


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="GPU-accelerated permutation testing")
    parser.add_argument("--celltype", type=str, default=None,
                        help="Cell type name (e.g., LEC, CD8_T_cells)")
    parser.add_argument("--annotation-level", type=str, default=None,
                        help="Annotation level: HyMy or L1")
    parser.add_argument("--n-perms", type=int, default=N_PERMUTATIONS,
                        help=f"Number of permutations (default: {N_PERMUTATIONS})")
    parser.add_argument("--adata-path", type=str, default=None,
                        help="Path to h5ad file")
    args = parser.parse_args()

    # Resolve annotation level
    annotation_level = (args.annotation_level
                        or os.environ.get("ANNOTATION_LEVEL", "HyMy"))

    # Resolve cell type (from arg, env var, or error)
    celltype_index = int(os.environ.get("CELLTYPE_INDEX", "0"))
    if args.celltype:
        celltype_name = args.celltype
    elif celltype_index > 0:
        celltype_name = CELLTYPE_INDEX_MAP.get(celltype_index)
        if celltype_name is None:
            sys.exit(f"Invalid CELLTYPE_INDEX: {celltype_index}")
    else:
        sys.exit("Must specify --celltype or set CELLTYPE_INDEX env var")

    n_perms = args.n_perms
    if args.adata_path:
        adata_path = args.adata_path
    elif ADATA_PATH is not None:
        adata_path = str(ADATA_PATH)
    else:
        sys.exit("No h5ad file found. Set ADATA_PATH env var or use --adata-path.")

    # Resolve analysis name and model type
    analysis_name = os.environ.get("ANALYSIS_NAME", "hymy_distance_correlation")
    use_poisson = "v2" in analysis_name
    k_neighbors = int(os.environ.get("K_NEIGHBORS", "1"))

    # Configure paths — 3-tier env var resolution (matches utils.R)
    query_celltype_env = os.environ.get("QUERY_CELLTYPE", "")
    celltype_col_env = os.environ.get("CELLTYPE_COLUMN", "")

    if query_celltype_env and celltype_col_env:
        query_celltype = query_celltype_env
        celltype_col = celltype_col_env
        suffix = f"_{query_celltype}"
    elif annotation_level == "L1":
        query_celltype = "IL1B_myeloid"
        celltype_col = "cell_type_assignment_L1"
        suffix = "_L1"
    else:
        query_celltype = "HyMy_GMM"
        celltype_col = "cell_type_with_HyMy_GMM"  # h5ad column name
        suffix = ""

    if OUTPUT_DIR:
        results_base = Path(OUTPUT_DIR) / f"spatial_analysis{suffix}" / analysis_name
    else:
        results_base = PROJECT_ROOT / "results" / f"spatial_analysis{suffix}" / analysis_name

    ct_dir = results_base / "per_celltype" / celltype_name

    print("=" * 70)
    print("GPU-Accelerated Permutation Testing")
    print("=" * 70)
    print(f"Analysis name:    {analysis_name}")
    print(f"Model:            {'Poisson GLM + offset' if use_poisson else 'Logistic regression'}")
    print(f"K neighbors:      {k_neighbors}")
    print(f"Annotation level: {annotation_level}")
    print(f"Query cell type:  {query_celltype}")
    print(f"Target cell type: {celltype_name}")
    print(f"Cell type column: {celltype_col}")
    print(f"N permutations:   {n_perms}")
    print(f"Results dir:      {ct_dir}")

    # Check GPU
    if torch.cuda.is_available():
        device = torch.device("cuda")
        print(f"GPU: {torch.cuda.get_device_name(0)}")
        print(f"GPU memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
    else:
        print("WARNING: No GPU detected, falling back to CPU (will be slow)")
        device = torch.device("cpu")

    # Check meta-analysis results exist (from R script with N_PERMUTATIONS=0)
    meta_file = ct_dir / "meta_analysis_results.csv"
    if not meta_file.exists():
        sys.exit(f"meta_analysis_results.csv not found at {meta_file}\n"
                 f"Run R script with N_PERMUTATIONS=0 first.")

    meta_results = pd.read_csv(meta_file)
    print(f"\nLoaded meta-analysis results: {len(meta_results)} genes")

    # Determine genes for permutation testing
    # Top N by absolute coefficient + priority genes (same logic as R script)
    top_by_effect = (meta_results
                     .sort_values("combined_coef", key=abs, ascending=False)
                     .head(PERM_TOP_N)["gene"].tolist())

    # Load h5ad to get available genes
    print(f"\nLoading h5ad: {adata_path}")
    t0 = time.time()
    import anndata as ad
    adata = ad.read_h5ad(adata_path)
    print(f"  Loaded in {time.time() - t0:.1f}s: {adata.shape[0]} cells x {adata.shape[1]} genes")

    available_genes = set(adata.var_names)
    analyzed_genes = set(meta_results["gene"].tolist())

    priority_in_data = [g for g in PERM_PRIORITY_GENES
                        if g in available_genes and g in analyzed_genes]
    perm_genes = list(dict.fromkeys(top_by_effect + priority_in_data))  # Deduplicate, preserve order

    print(f"  Permutation genes: {len(top_by_effect)} top by effect + "
          f"{len(priority_in_data)} priority = {len(perm_genes)} total")

    # Compute total counts per cell (for Poisson offset) BEFORE filtering
    if use_poisson:
        print("\nComputing total counts per cell for Poisson offset...")
        if sparse.issparse(adata.X):
            total_counts_all = np.array(adata.X.sum(axis=1)).flatten()
        else:
            total_counts_all = np.array(adata.X.sum(axis=1)).flatten()

        # Sanity check: .X should contain integer counts (not normalized)
        if sparse.issparse(adata.X):
            sample_vals = adata.X[:100].toarray().flatten()
        else:
            sample_vals = adata.X[:100].flatten()
        non_zero = sample_vals[sample_vals > 0]
        if len(non_zero) > 0 and not np.allclose(non_zero, np.round(non_zero)):
            print("  WARNING: .X appears to contain non-integer values!")
            print("  Checking for 'counts' layer...")
            if "counts" in adata.layers:
                print("  Using adata.layers['counts'] instead")
                if sparse.issparse(adata.layers["counts"]):
                    total_counts_all = np.array(
                        adata.layers["counts"].sum(axis=1)
                    ).flatten()
                else:
                    total_counts_all = np.array(
                        adata.layers["counts"].sum(axis=1)
                    ).flatten()
            else:
                print("  WARNING: No 'counts' layer found. Proceeding with .X")

        log_total_counts_all = np.log(np.maximum(total_counts_all, 1.0))
        print(f"  Total counts range: {total_counts_all.min():.0f} - {total_counts_all.max():.0f}")
        print(f"  Median total counts: {np.median(total_counts_all):.0f}")

    # Condition filtering
    if CONDITION_COL and CONDITION_VAL:
        print(f"\nFiltering to {CONDITION_COL} == {CONDITION_VAL}...")
        cond_mask = adata.obs[CONDITION_COL] == CONDITION_VAL
    elif CONDITION_COL:
        print(f"\nNo CONDITION_VALUE set; using all cells (CONDITION_COLUMN={CONDITION_COL} present but unfiltered)")
        cond_mask = np.ones(len(adata), dtype=bool)
    elif "group" in adata.obs.columns and not INPUT_PATH:
        # Legacy CeMM behavior: filter to TDLN
        print("\nFiltering to TDLN (legacy)...")
        cond_mask = adata.obs["group"] == "TDLN"
    else:
        print("\nNo condition filtering applied (using all cells)")
        cond_mask = np.ones(len(adata), dtype=bool)

    cond_indices_in_full = np.where(cond_mask)[0]
    adata = adata[cond_mask].copy()
    print(f"  Cells after filtering: {adata.shape[0]}")

    # Subset total counts to filtered cells
    if use_poisson:
        log_total_counts_filtered = log_total_counts_all[cond_indices_in_full]

    # Extract data (using configurable column names)
    # Spatial coordinates: try obsm["spatial"], then obs columns
    x_col_env = os.environ.get("X_COLUMN", "")
    y_col_env = os.environ.get("Y_COLUMN", "")
    if x_col_env and y_col_env:
        # User-specified coordinate columns
        coords = np.column_stack([
            adata.obs[x_col_env].values.astype(np.float32),
            adata.obs[y_col_env].values.astype(np.float32),
        ])
        print(f"  Using coordinate columns: {x_col_env}, {y_col_env}")
    elif "spatial" in adata.obsm:
        coords = np.array(adata.obsm["spatial"], dtype=np.float32)
        print("  Using adata.obsm['spatial'] for coordinates")
    else:
        # Auto-detect from obs columns
        for xc, yc in [("spatial_x", "spatial_y"), ("x", "y"),
                        ("x_centroid", "y_centroid")]:
            if xc in adata.obs.columns and yc in adata.obs.columns:
                coords = np.column_stack([
                    adata.obs[xc].values.astype(np.float32),
                    adata.obs[yc].values.astype(np.float32),
                ])
                print(f"  Auto-detected coordinate columns: {xc}, {yc}")
                break
        else:
            sys.exit("Could not find spatial coordinates. "
                     "Set X_COLUMN and Y_COLUMN env vars.")

    cell_types = adata.obs[celltype_col].values
    sample_ids = adata.obs[SAMPLE_COL].values

    # Identify query cells
    query_mask = cell_types == query_celltype
    n_query = query_mask.sum()
    print(f"  Query cells ({query_celltype}): {n_query}")

    # Query cells per sample (for stratified permutation)
    query_per_sample = {}
    for samp in np.unique(sample_ids):
        n = int((query_mask & (sample_ids == samp)).sum())
        if n > 0:
            query_per_sample[samp] = n
    print(f"  Query per sample: {query_per_sample}")

    # Identify target cells
    target_types = TARGET_CELLTYPES[celltype_name]
    target_mask = np.isin(cell_types, target_types)
    n_target = target_mask.sum()
    print(f"  Target cells ({celltype_name}): {n_target}")

    # Valid samples (enough target cells)
    target_sample_ids = sample_ids[target_mask]
    unique_samples, sample_counts = np.unique(target_sample_ids, return_counts=True)
    valid_samples = set(unique_samples[sample_counts >= MIN_CELLS_PER_SAMPLE])
    valid_samples = valid_samples & set(query_per_sample.keys())
    print(f"  Valid samples: {len(valid_samples)}")

    if len(valid_samples) < 2:
        sys.exit("Need at least 2 valid samples for meta-analysis")

    # Filter target to valid samples
    target_in_valid = target_mask & np.isin(sample_ids, list(valid_samples))
    target_indices = np.where(target_in_valid)[0]
    target_barcodes = adata.obs_names[target_indices]
    target_sample_ids = sample_ids[target_indices]

    # Get expression matrix for target cells
    print("\nExtracting expression matrix for permutation genes...")
    gene_indices = [list(adata.var_names).index(g) for g in perm_genes
                    if g in available_genes]
    perm_genes_available = [g for g in perm_genes if g in available_genes]

    if sparse.issparse(adata.X):
        expr_target = np.array(adata.X[target_indices][:, gene_indices].toarray(),
                               dtype=np.float32)
    else:
        expr_target = np.array(adata.X[target_indices][:, gene_indices],
                               dtype=np.float32)

    print(f"  Expression matrix: {expr_target.shape}")

    # Subset log_total_counts to target cells (for Poisson)
    if use_poisson:
        log_total_counts_target = log_total_counts_filtered[target_indices]
        print(f"  Log total counts for target: {log_total_counts_target.shape}")

    # Move coordinates to GPU
    print(f"\nMoving coordinates to {device}...")
    all_coords_gpu = torch.tensor(coords, dtype=torch.float32, device=device)
    target_coords_gpu = torch.tensor(coords[target_indices], dtype=torch.float32,
                                     device=device)

    # Get observed coefficients from meta-analysis
    meta_dict = dict(zip(meta_results["gene"], meta_results["combined_coef"]))

    # Filter query_per_sample to valid samples only
    query_per_sample_valid = {s: n for s, n in query_per_sample.items()
                              if s in valid_samples}

    # ==========================================================================
    # Run Permutation Tests
    # ==========================================================================

    print(f"\n{'=' * 70}")
    print(f"Running {n_perms} permutations for {len(perm_genes_available)} genes on {device}")
    print(f"{'=' * 70}")

    rng = np.random.default_rng(42)
    results = []
    t_start = time.time()

    for i, gene in enumerate(perm_genes_available):
        t_gene = time.time()

        gene_col_idx = perm_genes_available.index(gene)
        gene_expr = expr_target[:, gene_col_idx]

        # v1: binarize expression; v2: use raw counts
        if use_poisson:
            gene_data = gene_expr.astype(np.float64)
        else:
            gene_data = (gene_expr > 0).astype(np.float64)

        observed = meta_dict.get(gene, np.nan)
        if not np.isfinite(observed):
            results.append({"gene": gene, "perm_pval": np.nan})
            continue

        perm_pval = run_permutation_test_gpu(
            gene_data=gene_data,
            target_coords_gpu=target_coords_gpu,
            all_coords_gpu=all_coords_gpu,
            sample_ids_target=target_sample_ids,
            sample_ids_all=sample_ids,
            query_per_sample=query_per_sample_valid,
            observed_coef=observed,
            n_perms=n_perms,
            device=device,
            rng=rng,
            k_neighbors=k_neighbors,
            use_poisson=use_poisson,
            log_total_counts_target=(log_total_counts_target if use_poisson
                                     else None),
        )

        results.append({"gene": gene, "perm_pval": perm_pval})

        elapsed = time.time() - t_gene
        total_elapsed = time.time() - t_start
        eta = total_elapsed / (i + 1) * (len(perm_genes_available) - i - 1)

        if (i + 1) % 10 == 0 or (i + 1) == len(perm_genes_available):
            print(f"  [{i+1}/{len(perm_genes_available)}] {gene}: "
                  f"perm_pval={perm_pval:.4f} ({elapsed:.1f}s) "
                  f"ETA: {eta/60:.1f}min")

    total_time = time.time() - t_start
    print(f"\nCompleted in {total_time/60:.1f} minutes")

    # ==========================================================================
    # Save Results
    # ==========================================================================

    results_df = pd.DataFrame(results)
    output_path = ct_dir / "permutation_pvals.csv"
    results_df.to_csv(output_path, index=False)
    print(f"\nSaved: {output_path}")
    print(f"  Genes tested: {len(results_df)}")
    print(f"  Significant (p<0.05): {(results_df['perm_pval'] < 0.05).sum()}")

    # Clean up GPU memory
    if torch.cuda.is_available():
        torch.cuda.empty_cache()

    print(f"\n{'=' * 70}")
    print("GPU Permutation Testing Complete!")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
