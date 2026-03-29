#!/usr/bin/env Rscript
# =============================================================================
# Merge Distance Correlation Results from Array Jobs (Logistic Regression)
# =============================================================================
# Run after all array jobs complete to create combined summary files.
#
# Uses logistic regression: P(expressing) ~ distance to HyMy
# Gradient score = log-odds coefficient (negative = HyMy-induced)
#
# Usage:
#   Rscript merge_distance_correlation_results.R
#   ANNOTATION_LEVEL=L1 Rscript merge_distance_correlation_results.R
#
# This script:
#   1. Reads per-celltype meta_analysis_results.csv files
#   2. Combines them into all_genes_results.csv
#   3. Creates top_gradient_genes.csv (top 50 per cell type by FDR)
#   4. Creates decay_pattern_summary.csv (pattern counts per cell type)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# =============================================================================
# Configuration
# =============================================================================

ANNOTATION_LEVEL <- Sys.getenv("ANNOTATION_LEVEL", unset = "HyMy")
ANALYSIS_NAME <- Sys.getenv("ANALYSIS_NAME", unset = "hymy_distance_correlation_v2")
FDR_THRESHOLD <- 0.05

# Set output base path based on annotation level
if (.Platform$OS.type == "windows") {
  BASE_PATH <- "N:/lab_maier/Projects/mXenium"
} else {
  BASE_PATH <- "/nobackup/lab_maier/Projects/mXenium"
}

if (ANNOTATION_LEVEL == "HyMy") {
  OUTPUT_BASE <- file.path(BASE_PATH, "CMM/results/spatial_analysis", ANALYSIS_NAME)
} else {
  OUTPUT_BASE <- file.path(BASE_PATH, "CMM/results/spatial_analysis_L1", ANALYSIS_NAME)
}

message(strrep("=", 70))
message("Merge Distance Correlation Results")
message(strrep("=", 70))
message("Annotation level: ", ANNOTATION_LEVEL)
message("Output base: ", OUTPUT_BASE)

# =============================================================================
# Find and Load Per-Celltype Results
# =============================================================================

celltype_dir <- file.path(OUTPUT_BASE, "per_celltype")
if (!dir.exists(celltype_dir)) {
  stop("Per-celltype directory not found: ", celltype_dir)
}

celltype_dirs <- list.dirs(celltype_dir, recursive = FALSE)
message("\nFound ", length(celltype_dirs), " cell type directories:")
for (d in celltype_dirs) {
  message("  - ", basename(d))
}

# Load and combine all meta-analysis results
all_results <- rbindlist(lapply(celltype_dirs, function(d) {
  f <- file.path(d, "meta_analysis_results.csv")
  if (file.exists(f)) {
    dt <- fread(f)
    dt[, cell_type := basename(d)]
    message("  Loaded ", basename(d), ": ", nrow(dt), " genes")
    return(dt)
  } else {
    message("  WARNING: ", basename(d), " - meta_analysis_results.csv not found")
    return(NULL)
  }
}), fill = TRUE)

if (nrow(all_results) == 0) {
  stop("No results found to merge!")
}

message("\nTotal genes loaded: ", nrow(all_results))

# =============================================================================
# Create Summary Files
# =============================================================================

summary_dir <- file.path(OUTPUT_BASE, "summary")
dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)

# 1. All results combined
fwrite(all_results, file.path(summary_dir, "all_genes_results.csv"))
message("\nSaved: summary/all_genes_results.csv")

# 2. Top genes per cell type (by FDR)
top_genes <- all_results[fdr < FDR_THRESHOLD][order(fdr), head(.SD, 50), by = cell_type]
fwrite(top_genes, file.path(summary_dir, "top_gradient_genes.csv"))
message("Saved: summary/top_gradient_genes.csv (", nrow(top_genes), " genes)")

# 2b. Top genes by Fisher FDR (if available)
if ("fisher_fdr" %in% names(all_results)) {
  top_fisher <- all_results[fisher_fdr < FDR_THRESHOLD][order(fisher_fdr), head(.SD, 50), by = cell_type]
  fwrite(top_fisher, file.path(summary_dir, "top_gradient_genes_fisher.csv"))
  message("Saved: summary/top_gradient_genes_fisher.csv (", nrow(top_fisher), " genes)")
}

# 3. Decay pattern summary
if ("decay_pattern" %in% names(all_results)) {
  decay_summary <- all_results[fdr < FDR_THRESHOLD, .N, by = .(cell_type, decay_pattern)]
  setorder(decay_summary, cell_type, -N)
  fwrite(decay_summary, file.path(summary_dir, "decay_pattern_summary.csv"))
  message("Saved: summary/decay_pattern_summary.csv")
}

# 4. Permutation results summary (if available)
if ("perm_pval" %in% names(all_results)) {
  perm_summary <- all_results[!is.na(perm_pval), .(
    n_tested = .N,
    n_perm_sig = sum(perm_pval < 0.05),
    pct_perm_sig = round(mean(perm_pval < 0.05) * 100, 1)
  ), by = cell_type]
  fwrite(perm_summary, file.path(summary_dir, "permutation_summary.csv"))
  message("Saved: summary/permutation_summary.csv")
}

# =============================================================================
# Print Summary Report
# =============================================================================

message("\n", strrep("=", 70))
message("Results Summary")
message(strrep("=", 70))

has_fisher <- "fisher_fdr" %in% names(all_results)
summary_stats <- all_results[, .(
  n_genes = .N,
  n_significant = sum(fdr < FDR_THRESHOLD, na.rm = TRUE),
  n_fisher_sig = if (has_fisher) sum(fisher_fdr < FDR_THRESHOLD, na.rm = TRUE) else NA_integer_,
  n_perm_tested = sum(!is.na(perm_pval)),
  n_perm_sig = sum(perm_pval < 0.05, na.rm = TRUE)
), by = cell_type]

print(summary_stats)

message("\n", strrep("=", 70))
message("Done!")
message(strrep("=", 70))
