# Test RIPPLE on real data: CosMx NSCLC (He et al., 2022)
# This script:
#   1. Loads the h5ad file via zellkonverter
#   2. Inspects the data structure
#   3. Runs run_ripple() with tumor cells as query
#   4. Reports results
#
# Run with: Rscript data-raw/cosmx_nsclc/test_ripple_cosmx.R

library(zellkonverter)
library(SpatialExperiment)
library(SummarizedExperiment)

cat("=== Step 1: Load h5ad ===\n")
h5ad_path <- "data-raw/cosmx_nsclc/cosmx_nsclc.h5ad"
if (!file.exists(h5ad_path)) {
  stop("h5ad file not found at: ", h5ad_path)
}

sce <- readH5AD(h5ad_path)
cat("Class:", paste(class(sce), collapse = ", "), "\n")
cat("Dims:", nrow(sce), "genes x", ncol(sce), "cells\n")

cat("\n=== Step 2: Inspect metadata ===\n")
cd <- as.data.frame(colData(sce))
cat("colData columns:", paste(names(cd), collapse = ", "), "\n")

# Check for cell type column
celltype_candidates <- grep("cell.?type|cluster|annotation|leiden|celltype",
  names(cd),
  ignore.case = TRUE, value = TRUE
)
cat("Cell type candidates:", paste(celltype_candidates, collapse = ", "), "\n")

# Check for sample/patient column
sample_candidates <- grep("sample|patient|slide|fov|section|donor|batch",
  names(cd),
  ignore.case = TRUE, value = TRUE
)
cat("Sample candidates:", paste(sample_candidates, collapse = ", "), "\n")

# Print first few rows
cat("\nFirst 5 colData rows:\n")
print(head(cd, 5))

# Cell type distribution
if (length(celltype_candidates) > 0) {
  ct_col <- celltype_candidates[1]
  cat("\nCell type distribution (", ct_col, "):\n", sep = "")
  print(sort(table(cd[[ct_col]]), decreasing = TRUE))
}

# Sample distribution
if (length(sample_candidates) > 0) {
  for (s in sample_candidates) {
    cat("\nSample column '", s, "':\n", sep = "")
    print(table(cd[[s]]))
  }
}

# Check for spatial coordinates
cat("\n=== Spatial coordinates ===\n")
coord_candidates <- grep("^x$|^y$|centroid|spatial|coord",
  names(cd),
  ignore.case = TRUE, value = TRUE
)
cat("Coordinate candidates in colData:", paste(coord_candidates, collapse = ", "), "\n")

# Check reducedDims
rd_names <- reducedDimNames(sce)
cat("reducedDimNames:", paste(rd_names, collapse = ", "), "\n")

# Check assay names
cat("assayNames:", paste(assayNames(sce), collapse = ", "), "\n")

cat("\n=== Step 3: Summary ===\n")
cat("Ready to decide on:\n")
cat("  - celltype_column = ?\n")
cat("  - sample_column = ?\n")
cat("  - query_celltype = ?\n")
cat("  - coordinate columns = ?\n")
