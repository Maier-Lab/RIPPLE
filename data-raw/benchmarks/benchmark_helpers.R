# ============================================================================
# Shared helpers for RIPPLE benchmarks
# ============================================================================
# Generates synthetic SpatialExperiment datasets with controlled spatial
# structure and optional distance-dependent gene expression gradients.
#
# Used by: bench_null.R, bench_power.R, bench_comparison.R
# ============================================================================

library(SpatialExperiment)
library(S4Vectors)
library(Matrix)

# ---------------------------------------------------------------------------
# generate_benchmark_data()
#
# Args:
#   n_samples       - number of biological replicates
#   n_gradient_neg  - genes with negative gradient (induced near query)
#   n_gradient_pos  - genes with positive gradient (repressed near query)
#   n_background    - genes with no distance dependence
#   beta            - gradient strength (log-rate change per um); only used
#                     if n_gradient_neg > 0 or n_gradient_pos > 0
#   cells_per_sample - named list of cell counts per type per sample
#   field_um        - size of spatial field (um)
#   tumor_radius    - radius of query cell cluster (um)
#   background_rate - baseline Poisson rate for all genes
#   seed            - random seed (set before generation)
#
# Returns: a SpatialExperiment with columns cell_type, sample_id in colData
# ---------------------------------------------------------------------------
generate_benchmark_data <- function(
    n_samples = 3,
    n_gradient_neg = 0,
    n_gradient_pos = 0,
    n_background = 50,
    beta = -0.01,
    cells_per_sample = list(Tumor = 30, T_cell = 120, Fibroblast = 50),
    field_um = 500,
    tumor_radius = 60,
    background_rate = 3,
    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_genes <- n_gradient_neg + n_gradient_pos + n_background
  gene_names <- character(n_genes)
  idx <- 0
  if (n_gradient_neg > 0) {
    gene_names[idx + seq_len(n_gradient_neg)] <-
      paste0("GRAD_NEG_", seq_len(n_gradient_neg))
    idx <- idx + n_gradient_neg
  }
  if (n_gradient_pos > 0) {
    gene_names[idx + seq_len(n_gradient_pos)] <-
      paste0("GRAD_POS_", seq_len(n_gradient_pos))
    idx <- idx + n_gradient_pos
  }
  if (n_background > 0) {
    gene_names[idx + seq_len(n_background)] <-
      paste0("BG_", sprintf("%02d", seq_len(n_background)))
  }

  sample_disc <- function(n, cx, cy, radius) {
    r <- radius * sqrt(runif(n))
    theta <- runif(n, 0, 2 * pi)
    cbind(cx + r * cos(theta), cy + r * sin(theta))
  }

  build_one <- function(sid) {
    cx <- field_um / 2 + runif(1, -30, 30)
    cy <- field_um / 2 + runif(1, -30, 30)

    n_tumor <- cells_per_sample$Tumor
    n_tcell <- cells_per_sample$T_cell
    n_fib <- cells_per_sample$Fibroblast
    n_total <- n_tumor + n_tcell + n_fib

    tumor_xy <- sample_disc(n_tumor, cx, cy, tumor_radius)
    tcell_xy <- cbind(
      runif(n_tcell, 0, field_um),
      runif(n_tcell, 0, field_um)
    )
    fib_r <- runif(n_fib, 100, 200)
    fib_theta <- runif(n_fib, 0, 2 * pi)
    fib_xy <- cbind(
      pmax(0, pmin(field_um, cx + fib_r * cos(fib_theta))),
      pmax(0, pmin(field_um, cy + fib_r * sin(fib_theta)))
    )

    xy <- rbind(tumor_xy, tcell_xy, fib_xy)
    colnames(xy) <- c("x", "y")
    ct <- c(
      rep("Tumor", n_tumor), rep("T_cell", n_tcell),
      rep("Fibroblast", n_fib)
    )

    # Distance to nearest tumor cell (fast via RANN kd-tree)
    dist_to_tumor <- RANN::nn2(tumor_xy, xy, k = 1)$nn.dists[, 1]

    # Cell-level total counts (varies to test offset robustness)
    total_counts <- rpois(n_total, lambda = background_rate * n_genes)
    total_counts <- pmax(total_counts, n_genes)
    size_factor <- total_counts / (background_rate * n_genes)

    # Vectorized Poisson sampling
    # Base rate per (gene, cell): background_rate * size_factor[cell]
    # Dimension: n_genes rows x n_total cols
    rate_mat <- matrix(background_rate, nrow = n_genes, ncol = n_total)
    rate_mat <- rate_mat * rep(size_factor, each = n_genes)

    # Apply gradient modulation only to T_cell target cells for gradient genes
    n_grad <- n_gradient_neg + n_gradient_pos
    if (n_grad > 0) {
      tcell_idx <- which(ct == "T_cell")
      if (length(tcell_idx) > 0) {
        if (n_gradient_neg > 0) {
          # exp(beta * d) per T cell, repeated across the n_gradient_neg genes
          mod <- exp(beta * dist_to_tumor[tcell_idx])
          rate_mat[seq_len(n_gradient_neg), tcell_idx] <-
            rate_mat[seq_len(n_gradient_neg), tcell_idx] *
              rep(mod, each = n_gradient_neg)
        }
        if (n_gradient_pos > 0) {
          pos_rows <- (n_gradient_neg + 1):(n_gradient_neg + n_gradient_pos)
          mod <- exp(-beta * dist_to_tumor[tcell_idx])
          rate_mat[pos_rows, tcell_idx] <-
            rate_mat[pos_rows, tcell_idx] *
              rep(mod, each = n_gradient_pos)
        }
      }
    }

    counts <- matrix(rpois(length(rate_mat), pmax(as.vector(rate_mat), 0.1)),
      nrow = n_genes, ncol = n_total
    )
    rownames(counts) <- gene_names
    cnames <- paste0(sid, "_cell", seq_len(n_total))
    colnames(counts) <- cnames

    list(
      counts = counts,
      meta = data.frame(
        cell_type = ct, sample_id = sid,
        row.names = cnames, stringsAsFactors = FALSE
      ),
      coords = xy
    )
  }

  sids <- paste0("sample_", seq_len(n_samples))
  samples <- lapply(sids, build_one)

  all_counts <- do.call(cbind, lapply(samples, `[[`, "counts"))
  all_meta <- do.call(rbind, lapply(samples, `[[`, "meta"))
  all_coords <- do.call(rbind, lapply(samples, `[[`, "coords"))
  rownames(all_coords) <- rownames(all_meta)

  SpatialExperiment(
    assays        = list(counts = methods::as(all_counts, "CsparseMatrix")),
    colData       = DataFrame(all_meta),
    spatialCoords = as.matrix(all_coords)
  )
}

# ---------------------------------------------------------------------------
# run_ripple_quiet()
#
# Runs run_ripple() in a temp directory with minimal output.
# Returns the results data.table.
# ---------------------------------------------------------------------------
run_ripple_quiet <- function(spe, ...) {
  out_dir <- tempfile("ripple_bench_")
  dir.create(out_dir, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  results <- run_ripple(
    input           = spe,
    query_celltype  = "Tumor",
    celltype_column = "cell_type",
    sample_column   = "sample_id",
    output_dir      = out_dir,
    analysis_name   = "bench",
    verbose         = FALSE,
    n_permutations  = 0,
    ...
  )
  results
}
