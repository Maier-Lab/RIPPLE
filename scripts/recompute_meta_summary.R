#!/usr/bin/env Rscript
# =============================================================================
# Recompute Meta-Analysis Summary with Fisher's Combined P-value
# =============================================================================
#
# Reads existing coef_per_sample.csv files and adds new columns to
# meta_analysis_results.csv:
#   - median_coef: median of per-sample coefficients (equal mouse weighting)
#   - fisher_stat: Fisher's combined test statistic (-2 * sum(log(p_i)))
#   - fisher_pval: Fisher's combined p-value (with sign consistency gate)
#   - fisher_fdr: BH-adjusted Fisher p-value
#
# This script does NOT re-run any GLMs. It only reads existing per-sample
# results and computes summary statistics.
#
# Fisher's method: X^2 = -2 * sum(log(p_i)), df = 2k
# Sign consistency gate: genes with < 75% directional agreement across mice
#   get fisher_pval = 1 (not significant), because Fisher's method does not
#   account for effect direction.
#
# Usage:
#   Rscript recompute_meta_summary.R
#   ANALYSIS_NAME=hymy_distance_correlation_v2 Rscript recompute_meta_summary.R
#   QUERY_CELLTYPE=MyType CELLTYPE_COLUMN=my_col Rscript recompute_meta_summary.R
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# =============================================================================
# Configuration (sourced from config.R — single source of truth)
# =============================================================================

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

SIGN_CONSISTENCY_THRESHOLD <- 1.0  # All mice must agree on direction (4/4)

OUTPUT_BASE <- file.path(OUTPUT_ROOT, ANALYSIS_NAME)

message(strrep("=", 70))
message("Recompute Meta-Analysis Summary (Fisher's + Sign Consistency)")
message(strrep("=", 70))
message("Annotation level: ", ANNOTATION_LEVEL)
message("Analysis name: ", ANALYSIS_NAME)
message("Output base: ", OUTPUT_BASE)
message("Sign consistency threshold: ", SIGN_CONSISTENCY_THRESHOLD)

# =============================================================================
# Fisher's Combined P-value Function
# =============================================================================

compute_fisher_pval <- function(pvals, coefs, min_samples = 2,
                                 sign_threshold = 0.75) {
  # Filter to valid entries
  valid <- !is.na(pvals) & !is.na(coefs) & pvals > 0
  valid_pvals <- pvals[valid]
  valid_coefs <- coefs[valid]
  n_valid <- length(valid_pvals)

  if (n_valid < min_samples) {
    return(list(fisher_stat = NA_real_, fisher_pval = NA_real_,
                median_coef = NA_real_, n_valid = n_valid))
  }

  # Median coefficient (equal mouse weighting)
  median_coef <- median(valid_coefs)

  # Sign consistency check
  n_pos <- sum(valid_coefs > 0)
  n_neg <- sum(valid_coefs < 0)
  sc <- max(n_pos, n_neg) / n_valid

  if (sc < sign_threshold) {
    # Contradictory directions — not significant
    return(list(fisher_stat = NA_real_, fisher_pval = 1.0,
                median_coef = median_coef, n_valid = n_valid))
  }

  # Fisher's method: X^2 = -2 * sum(log(p_i)), df = 2k
  clamped_pvals <- pmax(valid_pvals, 1e-15)  # avoid log(0)
  fisher_stat <- -2 * sum(log(clamped_pvals))
  fisher_pval <- pchisq(fisher_stat, df = 2 * n_valid, lower.tail = FALSE)

  list(fisher_stat = fisher_stat, fisher_pval = fisher_pval,
       median_coef = median_coef, n_valid = n_valid)
}

# =============================================================================
# Process Each Cell Type
# =============================================================================

celltype_dir <- file.path(OUTPUT_BASE, "per_celltype")
if (!dir.exists(celltype_dir)) {
  stop("Per-celltype directory not found: ", celltype_dir)
}

ct_dirs <- list.dirs(celltype_dir, recursive = FALSE)
message("\nFound ", length(ct_dirs), " cell type directories")

total_updated <- 0L

for (ct_path in ct_dirs) {
  ct_name <- basename(ct_path)
  coef_file <- file.path(ct_path, "coef_per_sample.csv")
  meta_file <- file.path(ct_path, "meta_analysis_results.csv")

  if (!file.exists(coef_file)) {
    message("  ", ct_name, ": coef_per_sample.csv not found, skipping")
    next
  }
  if (!file.exists(meta_file)) {
    message("  ", ct_name, ": meta_analysis_results.csv not found, skipping")
    next
  }

  # Read per-sample coefficients
  coef_data <- fread(coef_file)
  meta_data <- fread(meta_file)

  # Compute Fisher's p-value per gene
  genes <- unique(coef_data$gene)
  fisher_results <- rbindlist(lapply(genes, function(g) {
    gene_data <- coef_data[gene == g]
    result <- compute_fisher_pval(
      pvals = gene_data$pval,
      coefs = gene_data$coef,
      sign_threshold = SIGN_CONSISTENCY_THRESHOLD
    )
    data.table(
      gene = g,
      median_coef = result$median_coef,
      fisher_stat = result$fisher_stat,
      fisher_pval = result$fisher_pval
    )
  }))

  # BH correction within this cell type
  fisher_results[, fisher_fdr := p.adjust(fisher_pval, method = "BH")]

  # Merge into existing meta_analysis_results
  # Remove old Fisher columns if they exist (from previous runs)
  cols_to_remove <- intersect(
    c("median_coef", "fisher_stat", "fisher_pval", "fisher_fdr"),
    names(meta_data)
  )
  if (length(cols_to_remove) > 0) {
    meta_data[, (cols_to_remove) := NULL]
  }

  meta_updated <- merge(meta_data, fisher_results, by = "gene", all.x = TRUE)

  # Overwrite
  fwrite(meta_updated, meta_file)

  # Report
  n_fisher_sig <- sum(meta_updated$fisher_fdr < 0.05, na.rm = TRUE)
  n_meta_sig <- sum(meta_updated$fdr < 0.05, na.rm = TRUE)
  n_direction_blocked <- sum(meta_updated$fisher_pval == 1, na.rm = TRUE)
  message(sprintf("  %s: %d genes | meta FDR<0.05: %d | Fisher FDR<0.05: %d | direction-blocked: %d",
                  ct_name, nrow(meta_updated), n_meta_sig, n_fisher_sig, n_direction_blocked))
  total_updated <- total_updated + nrow(meta_updated)
}

# =============================================================================
# Summary
# =============================================================================

message("\n", strrep("=", 70))
message("Done! Updated ", total_updated, " gene entries across ", length(ct_dirs), " cell types")
message("New columns: median_coef, fisher_stat, fisher_pval, fisher_fdr")
message(strrep("=", 70))
