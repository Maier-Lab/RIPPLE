# Run only Stage 6 (L-R integration) on CosMx NSCLC
# Prerequisite: Stage 1 must be complete
# Run with: Rscript data-raw/cosmx_nsclc/run_lr_only.R

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(Matrix)
})
pkgload::load_all(".", quiet = TRUE)

results_dir <- "data-raw/cosmx_nsclc/ripple_output/tumor_ripple"
h5ad_path <- "data-raw/cosmx_nsclc/cosmx_nsclc.h5ad"

cat("Loading h5ad and rebuilding SPE...\n")
sce <- readH5AD(h5ad_path)
cd <- colData(sce)
ct <- as.character(cd$cell_type)
ct[grepl("^tumor", ct)] <- "tumor"

new_cd <- DataFrame(
  cell_type = ct,
  patient = cd$patient
)
rownames(new_cd) <- colnames(sce)

spe <- SpatialExperiment(
  assays = list(counts = assay(sce, "counts")),
  colData = new_cd,
  spatialCoords = as.matrix(reducedDim(sce, "spatial"))
)
rm(sce)
gc(verbose = FALSE)
cat("SPE ready:", ncol(spe), "cells\n")

cat("\nRunning L-R integration...\n")
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
cat("L-R RESULTS\n")
cat(strrep("=", 70), "\n")

if (is.null(lr_results) || nrow(lr_results) == 0) {
  cat("No L-R results produced.\n")
} else {
  cat("Total L-R pairs:", nrow(lr_results), "\n")

  if ("combined_score" %in% names(lr_results)) {
    cat("\nTop 25 by combined score:\n")
    print(lr_results[order(-combined_score)][
      1:min(25, nrow(lr_results)),
      .(ligand, receptor, cell_type, combined_score, reproducibility)
    ])
  }

  # Check for key immune checkpoint pairs
  cat("\n=== Key L-R pairs from He et al. 2022 ===\n")
  check_pairs <- list(
    c("CD274", "PDCD1"),
    c("PDCD1", "CD274"),
    c("MIF", "CD74"),
    c("CD74", "MIF"),
    c("SPP1", "CD44"),
    c("CTLA4", "CD86"),
    c("CD86", "CTLA4"),
    c("TGFB1", "TGFBR1"),
    c("TGFB1", "TGFBR2"),
    c("CCL2", "CCR2"),
    c("CXCL12", "CXCR4")
  )

  for (pair in check_pairs) {
    match <- lr_results[ligand == pair[1] & receptor == pair[2]]
    if (nrow(match) > 0) {
      for (i in seq_len(nrow(match))) {
        cat(sprintf(
          "  %s -> %s in %s: score=%.3f, repro=%.2f\n",
          pair[1], pair[2], match$cell_type[i],
          match$combined_score[i],
          match$reproducibility[i]
        ))
      }
    }
  }
}

cat("\nDone.\n")
