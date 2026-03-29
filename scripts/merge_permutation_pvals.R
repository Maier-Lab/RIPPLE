#!/usr/bin/env Rscript
# =============================================================================
# RIPPLE Stage 3: Merge GPU Permutation P-values into Meta-Analysis Results
# =============================================================================
#
# Reads the GPU-produced permutation_pvals.csv and merges it into the
# R-produced meta_analysis_results.csv (from N_PERMUTATIONS=0 run).
#
# Usage:
#   Rscript merge_permutation_pvals.R
#   QUERY_CELLTYPE=MyType CELLTYPE_COLUMN=my_col Rscript merge_permutation_pvals.R
#   ANALYSIS_NAME=hymy_distance_correlation_v2 Rscript merge_permutation_pvals.R
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# Source lightweight config (no heavy package loads needed)
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  for (i in seq_len(sys.nframe())) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  }
  return(getwd())
}
source(file.path(get_script_dir(), "config.R"))

results_base <- file.path(OUTPUT_ROOT, ANALYSIS_NAME)

cat("==============================================\n")
cat("Merge GPU Permutation P-values\n")
cat("==============================================\n")
cat("Analysis name:", ANALYSIS_NAME, "\n")
cat("Annotation level:", ANNOTATION_LEVEL, "\n")
cat("Results base:", results_base, "\n")
cat("Date:", format(Sys.time()), "\n")
cat("==============================================\n\n")

# Cell types to process
CELL_TYPES <- c("LEC", "FRC", "BEC", "CD4_T_cells", "CD8_T_cells",
                "gdT_cells", "Macrophages", "Monocyte", "Fibroblasts_mac",
                "cDC1", "cDC2", "mature_migDC", "B_cells", "Plasma_cell")

ct_dirs <- file.path(results_base, "per_celltype", CELL_TYPES)
names(ct_dirs) <- CELL_TYPES

n_merged <- 0
n_skipped <- 0

for (ct in CELL_TYPES) {
  ct_dir <- ct_dirs[ct]
  meta_file <- file.path(ct_dir, "meta_analysis_results.csv")
  perm_file <- file.path(ct_dir, "permutation_pvals.csv")

  # Check meta-analysis results exist
  if (!file.exists(meta_file)) {
    message("  [SKIP] ", ct, ": meta_analysis_results.csv not found")
    n_skipped <- n_skipped + 1
    next
  }

  # Check GPU permutation results exist
  if (!file.exists(perm_file)) {
    message("  [SKIP] ", ct, ": permutation_pvals.csv not found (GPU job not run yet?)")
    n_skipped <- n_skipped + 1
    next
  }

  meta <- fread(meta_file)
  perm <- fread(perm_file)

  # Validate columns
  if (!"gene" %in% names(perm) || !"perm_pval" %in% names(perm)) {
    message("  [ERROR] ", ct, ": permutation_pvals.csv missing gene/perm_pval columns")
    n_skipped <- n_skipped + 1
    next
  }

  # Count how many genes had perm_pval already
  n_had_pval <- sum(!is.na(meta$perm_pval))
  if (n_had_pval > 0) {
    message("  [NOTE] ", ct, ": ", n_had_pval, " genes already had perm_pval; overwriting with GPU results")
  }

  # Remove old perm_pval column and merge new one
  meta[, perm_pval := NULL]
  meta <- merge(meta, perm[, .(gene, perm_pval)], by = "gene", all.x = TRUE)
  setDT(meta)

  # Summary stats
  n_tested <- sum(!is.na(perm$perm_pval))
  n_sig <- sum(perm$perm_pval < 0.05, na.rm = TRUE)
  n_genes <- nrow(meta)

  # Save updated meta-analysis results
  fwrite(meta, meta_file)

  message(sprintf("  [OK] %s: %d/%d genes with perm_pval (%d significant at p<0.05)",
                  ct, n_tested, n_genes, n_sig))
  n_merged <- n_merged + 1
}

cat("\n==============================================\n")
cat("Summary\n")
cat("==============================================\n")
cat("Cell types merged:", n_merged, "\n")
cat("Cell types skipped:", n_skipped, "\n")

if (n_merged > 0) {
  cat("\nNext steps:\n")
  cat("  1. Run merge_distance_correlation_results.R to update summary tables\n")
  cat("  2. Run Stage 2 if desired: sbatch run_hymy_distance_correlation_stage2_array.sh\n")
}

cat("==============================================\n")
