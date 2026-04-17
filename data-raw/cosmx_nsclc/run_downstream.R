# Run downstream RIPPLE stages on CosMx NSCLC
# Prerequisite: run_ripple_cosmx.R must have completed (Stage 1)
#
# This script runs:
#   - Stage 5: Atlas visualization + fGSEA pathway enrichment
#   - Stage 6: Ligand-receptor integration
#
# Run with: Rscript data-raw/cosmx_nsclc/run_downstream.R

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(Matrix)
})

# Load package
suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
})

results_dir <- "data-raw/cosmx_nsclc/ripple_output/tumor_ripple"
h5ad_path <- "data-raw/cosmx_nsclc/cosmx_nsclc.h5ad"

# =========================================================================
# Stage 5: Atlas + fGSEA
# =========================================================================
cat("\n", strrep("=", 70), "\n")
cat("STAGE 5: Atlas Visualization + fGSEA\n")
cat(strrep("=", 70), "\n")

run_ripple_atlas(
  results_dir = results_dir,
  output_dir  = file.path(results_dir, "atlas"),
  query_label = "Tumor",
  run_fgsea   = TRUE,
  organism    = "human",
  fgsea_seed  = 42,
  verbose     = TRUE
)

# =========================================================================
# Stage 6: L-R Integration
# =========================================================================
cat("\n", strrep("=", 70), "\n")
cat("STAGE 6: Ligand-Receptor Integration\n")
cat(strrep("=", 70), "\n")

# Load the SPE (rebuild from h5ad with merged tumor labels)
cat("Loading data for L-R integration...\n")
sce <- readH5AD(h5ad_path)
cd <- colData(sce)
ct <- as.character(cd$cell_type)
ct[grepl("^tumor", ct)] <- "tumor"

new_cd <- DataFrame(
  cell_type = ct,
  patient = cd$patient,
  sample = cd$sample
)
rownames(new_cd) <- colnames(sce)
spatial_coords <- reducedDim(sce, "spatial")

spe <- SpatialExperiment(
  assays        = list(counts = assay(sce, "counts")),
  colData       = new_cd,
  spatialCoords = as.matrix(spatial_coords)
)
rm(sce)
gc(verbose = FALSE)

cat("Running L-R integration...\n")
lr_results <- run_ripple_lr(
  results_dir = results_dir,
  input = spe,
  query_celltype = "tumor",
  celltype_column = "cell_type",
  sample_column = "patient",
  organism = "human",
  expr_threshold_pct = 5,
  fdr_threshold = 0.05,
  output_dir = file.path(dirname(results_dir), "lr_integration"),
  verbose = TRUE
)

cat("\n", strrep("=", 70), "\n")
cat("L-R RESULTS SUMMARY\n")
cat(strrep("=", 70), "\n")
cat("Total L-R pairs found:", nrow(lr_results), "\n")

if (nrow(lr_results) > 0 && "combined_score" %in% names(lr_results)) {
  cat("\nTop 20 L-R pairs by combined score:\n")
  print(lr_results[order(-combined_score)][
    1:min(20, nrow(lr_results)),
    .(ligand, receptor, cell_type, combined_score, reproducibility)
  ])
}

# Check for the paper's key L-R pairs
cat("\n=== Checking paper's key L-R pairs ===\n")
key_pairs <- data.frame(
  ligand = c("CD274", "CTLA4", "SPP1", "EGFR", "ERBB2", "MIF"),
  receptor = c("PDCD1", "CD86", "CD44", "EGF", "NRG1", "CD74"),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(key_pairs))) {
  lig <- key_pairs$ligand[i]
  rec <- key_pairs$receptor[i]
  if (nrow(lr_results) > 0) {
    match <- lr_results[ligand == lig & receptor == rec]
    if (nrow(match) > 0) {
      cat(sprintf(
        "  %s -> %s: FOUND in %s (score=%.3f)\n",
        lig, rec, match$cell_type[1], match$combined_score[1]
      ))
    } else {
      # Check reverse
      match2 <- lr_results[ligand == rec & receptor == lig]
      if (nrow(match2) > 0) {
        cat(sprintf(
          "  %s -> %s: FOUND (reversed) in %s (score=%.3f)\n",
          rec, lig, match2$cell_type[1], match2$combined_score[1]
        ))
      } else {
        cat(sprintf("  %s -> %s: not found\n", lig, rec))
      }
    }
  }
}

cat("\n=== All stages complete ===\n")
cat("Atlas output:", file.path(results_dir, "atlas"), "\n")
cat("L-R output:", file.path(dirname(results_dir), "lr_integration"), "\n")
