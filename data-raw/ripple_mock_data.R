# ============================================================================
# Build the ripple_mock_data dataset
# ============================================================================
# This script generates a small synthetic SpatialExperiment used in examples,
# vignettes, and tests. It has a planted distance-dependent gradient so RIPPLE
# can be verified to recover the known effect.
#
# Layout (per sample):
#   - Tumor cells form a dense central cluster
#   - T cells are scattered, with density biased toward the tumor edge
#   - Fibroblasts form a looser ring at the periphery
#
# Planted effects on T cells:
#   - 5 "induced" genes  (exponential decay with distance; expression high
#     near tumor, drops off at ~40 um)
#   - 5 "repressed" genes (expression low near tumor, rises with distance)
#   - 40 background genes (no distance dependence, Poisson noise only)
#
# Run with: source("data-raw/ripple_mock_data.R")
# ============================================================================

set.seed(42)

library(SpatialExperiment)
library(S4Vectors)
library(Matrix)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
n_samples <- 3
cells_per_sample <- list(Tumor = 30, T_cell = 100, Fibroblast = 70)
n_induced <- 5
n_repressed <- 5
n_background <- 40

field_um <- 500 # size of the spatial field
tumor_radius <- 60 # radius of the tumor cluster
decay_scale_um <- 40 # characteristic decay distance for planted genes

background_rate <- 2 # mean counts per cell for background genes
induced_low <- 1 # floor counts far from tumor
induced_high <- 12 # peak counts adjacent to tumor
repressed_low <- 1 # counts near tumor
repressed_high <- 8 # counts far from tumor

# ---------------------------------------------------------------------------
# Gene catalog
# ---------------------------------------------------------------------------
induced_genes <- paste0("INDUCED_", seq_len(n_induced))
repressed_genes <- paste0("REPRESSED_", seq_len(n_repressed))
background_genes <- paste0("BG_", sprintf("%02d", seq_len(n_background)))
all_genes <- c(induced_genes, repressed_genes, background_genes)
n_genes <- length(all_genes)

# ---------------------------------------------------------------------------
# Helper: sample points inside a 2D disc around a center
# ---------------------------------------------------------------------------
sample_disc <- function(n, cx, cy, radius) {
  r <- radius * sqrt(stats::runif(n))
  theta <- stats::runif(n, 0, 2 * pi)
  cbind(cx + r * cos(theta), cy + r * sin(theta))
}

# ---------------------------------------------------------------------------
# Build one sample at a time, then concatenate
# ---------------------------------------------------------------------------
build_sample <- function(sample_id) {
  # Tumor center varies slightly per sample to avoid identical layouts
  tumor_cx <- field_um / 2 + stats::runif(1, -30, 30)
  tumor_cy <- field_um / 2 + stats::runif(1, -30, 30)

  # Tumor cells: dense disc at center
  n_tumor <- cells_per_sample$Tumor
  tumor_xy <- sample_disc(n_tumor, tumor_cx, tumor_cy, tumor_radius)

  # T cells: uniform over the field
  n_tcell <- cells_per_sample$T_cell
  tcell_xy <- cbind(
    stats::runif(n_tcell, 0, field_um),
    stats::runif(n_tcell, 0, field_um)
  )

  # Fibroblasts: annular ring around tumor (radius 100-200)
  n_fib <- cells_per_sample$Fibroblast
  fib_r <- stats::runif(n_fib, 100, 200)
  fib_theta <- stats::runif(n_fib, 0, 2 * pi)
  fib_xy <- cbind(
    tumor_cx + fib_r * cos(fib_theta),
    tumor_cy + fib_r * sin(fib_theta)
  )
  # Clip to field
  fib_xy[, 1] <- pmax(0, pmin(field_um, fib_xy[, 1]))
  fib_xy[, 2] <- pmax(0, pmin(field_um, fib_xy[, 2]))

  all_xy <- rbind(tumor_xy, tcell_xy, fib_xy)
  colnames(all_xy) <- c("x", "y")

  cell_types <- c(
    rep("Tumor", n_tumor),
    rep("T_cell", n_tcell),
    rep("Fibroblast", n_fib)
  )

  # Distance from each cell to nearest Tumor cell
  # For Tumor cells themselves, distance is 0 (or tiny within the cluster)
  n_total <- length(cell_types)
  dist_to_tumor <- numeric(n_total)
  for (i in seq_len(n_total)) {
    dists <- sqrt(rowSums((tumor_xy - matrix(all_xy[i, ], n_tumor, 2,
      byrow = TRUE
    ))^2))
    dist_to_tumor[i] <- min(dists)
  }

  # Build counts matrix (genes x cells)
  counts <- matrix(0L, nrow = n_genes, ncol = n_total)
  rownames(counts) <- all_genes
  colnames(counts) <- paste0(sample_id, "_cell", seq_len(n_total))

  for (i in seq_len(n_total)) {
    if (cell_types[i] == "T_cell") {
      d <- dist_to_tumor[i]
      decay_factor <- exp(-d / decay_scale_um)

      # Induced: high near tumor, low far
      for (g in induced_genes) {
        rate <- induced_low + (induced_high - induced_low) * decay_factor
        counts[g, i] <- stats::rpois(1, rate)
      }
      # Repressed: low near tumor, high far
      for (g in repressed_genes) {
        rate <- repressed_high - (repressed_high - repressed_low) * decay_factor
        counts[g, i] <- stats::rpois(1, rate)
      }
      # Background: flat rate
      for (g in background_genes) {
        counts[g, i] <- stats::rpois(1, background_rate)
      }
    } else {
      # Non-T cells: all genes at background rate
      for (g in all_genes) {
        counts[g, i] <- stats::rpois(1, background_rate)
      }
    }
  }

  list(
    counts = counts,
    meta = data.frame(
      cell_type = cell_types,
      sample_id = sample_id,
      stringsAsFactors = FALSE,
      row.names = colnames(counts)
    ),
    coords = all_xy
  )
}

# ---------------------------------------------------------------------------
# Assemble all samples
# ---------------------------------------------------------------------------
sample_names <- paste0("sample_", seq_len(n_samples))
samples <- lapply(sample_names, build_sample)

all_counts <- do.call(cbind, lapply(samples, `[[`, "counts"))
all_meta <- do.call(rbind, lapply(samples, `[[`, "meta"))
all_coords <- do.call(rbind, lapply(samples, `[[`, "coords"))
rownames(all_coords) <- rownames(all_meta)

message("Synthetic dataset built:")
message("  Genes: ", nrow(all_counts))
message("  Cells: ", ncol(all_counts))
message("  Samples: ", length(unique(all_meta$sample_id)))
message(
  "  Cell types: ",
  paste(names(table(all_meta$cell_type)), collapse = ", ")
)
print(table(all_meta$cell_type, all_meta$sample_id))

# Convert counts to a sparse matrix to match real data
all_counts_sparse <- methods::as(all_counts, "CsparseMatrix")

# ---------------------------------------------------------------------------
# Build a SpatialExperiment
# ---------------------------------------------------------------------------
ripple_mock_data <- SpatialExperiment::SpatialExperiment(
  assays        = list(counts = all_counts_sparse),
  colData       = S4Vectors::DataFrame(all_meta),
  spatialCoords = as.matrix(all_coords)
)

message("\nSpatialExperiment object:")
print(ripple_mock_data)

# ---------------------------------------------------------------------------
# Save to data/ for the package
# ---------------------------------------------------------------------------
usethis::use_data(ripple_mock_data, overwrite = TRUE, compress = "xz")
message("\nSaved: data/ripple_mock_data.rda")
