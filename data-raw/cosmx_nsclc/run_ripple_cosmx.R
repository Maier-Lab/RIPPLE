# Run RIPPLE on the CosMx NSCLC dataset (He et al., 2022)
#
# Setup:
#   - Merge per-patient tumor labels into a single "tumor" label
#   - Use 'patient' as biological replicate (N = 5)
#   - Extract spatial coordinates from reducedDim('spatial') into colData
#
# Run with: Rscript data-raw/cosmx_nsclc/run_ripple_cosmx.R

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
})

cat("=== Step 1: Load h5ad ===\n")
sce <- readH5AD("data-raw/cosmx_nsclc/cosmx_nsclc.h5ad")
cat("Loaded:", nrow(sce), "genes x", ncol(sce), "cells\n")

cat("\n=== Step 2: Inspect spatial coordinates ===\n")
spatial_coords <- reducedDim(sce, "spatial")
cat("spatial dim:", dim(spatial_coords), "\n")
cat("class:", class(spatial_coords), "\n")
cat("first 3 rows:\n")
print(head(spatial_coords, 3))

# Per-patient coordinate ranges (sanity check that coords are stitched globally)
cd <- colData(sce)
cat("\nCoordinate ranges per patient:\n")
for (p in unique(cd$patient)) {
  idx <- cd$patient == p
  xr <- range(spatial_coords[idx, 1])
  yr <- range(spatial_coords[idx, 2])
  cat(sprintf(
    "  %-8s x = [%8.0f, %8.0f]  y = [%8.0f, %8.0f]\n",
    p, xr[1], xr[2], yr[1], yr[2]
  ))
}

cat("\n=== Step 3: Merge tumor subtypes ===\n")
ct <- as.character(cd$cell_type)
n_tumor_before <- sum(grepl("^tumor", ct))
ct[grepl("^tumor", ct)] <- "tumor"
n_tumor_after <- sum(ct == "tumor")
cat("Tumor cells before merge:", n_tumor_before, "\n")
cat("Tumor cells after merge:", n_tumor_after, "\n")
cat("\nFinal cell type distribution:\n")
print(sort(table(ct), decreasing = TRUE))

cat("\n=== Step 4: Build SpatialExperiment for RIPPLE ===\n")
# Use only assay 'counts' (raw integer counts) and a clean colData
new_cd <- DataFrame(
  cell_type = ct,
  patient = cd$patient,
  sample = cd$sample
)

# Make sure rownames are cell barcodes
rownames(new_cd) <- colnames(sce)

spe <- SpatialExperiment(
  assays        = list(counts = assay(sce, "counts")),
  colData       = new_cd,
  spatialCoords = as.matrix(spatial_coords)
)
cat("Built SPE:", nrow(spe), "genes x", ncol(spe), "cells\n")

# Free memory
rm(sce)
gc(verbose = FALSE)

cat("\n=== Step 5: Run RIPPLE ===\n")
suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
})

out_dir <- "data-raw/cosmx_nsclc/ripple_output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

results <- run_ripple(
  input                = spe,
  query_celltype       = "tumor",
  celltype_column      = "cell_type",
  sample_column        = "patient",
  output_dir           = out_dir,
  analysis_name        = "tumor_ripple",
  k_neighbors          = 1,
  max_distance_um      = 200,
  min_cells_per_sample = 30,
  min_expr_pct         = 0.01,
  min_expr_floor       = 25,
  n_permutations       = 0,
  verbose              = TRUE
)

cat("\n=== Step 6: Top hits ===\n")
cat("Total gene x cell_type rows:", nrow(results), "\n")
cat("Significant (fisher_fdr < 0.05):", sum(results$fisher_fdr < 0.05, na.rm = TRUE), "\n")

cat("\nTop 25 by fisher_fdr:\n")
print(results[order(fisher_fdr)][
  1:25,
  .(gene, cell_type, median_coef, fisher_fdr, sign_consistency, n_samples)
])

cat("\nSignificant hits per cell type:\n")
sig_counts <- results[fisher_fdr < 0.05, .N, by = cell_type][order(-N)]
print(sig_counts)

cat("\n=== Done ===\n")
cat("Results saved in:", out_dir, "\n")
