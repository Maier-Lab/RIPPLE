#!/usr/bin/env Rscript
#' =============================================================================
#' RIPPLE Stage 1: Distance Correlation Analysis (Poisson GLM)
#' =============================================================================
#'
#' Per-sample Poisson GLM: gene expression ~ distance to query cell type,
#' with cell-size offset. Results are combined via Fisher's combined p-value.
#'
#' Model: glm(counts ~ distance + offset(log(total_counts)), family=poisson)
#'
#' Coefficient Interpretation:
#' - Coef < 0: expression rate DECREASES with distance = query-INDUCED
#' - Coef > 0: expression rate INCREASES with distance = query-REPRESSED
#' - Units: log-rate change per um of distance
#'
#' Usage:
#'   Rscript hymy_distance_correlation_v2.R
#'   QUERY_CELLTYPE=MyType CELLTYPE_COLUMN=my_col Rscript hymy_distance_correlation_v2.R
#'   ANNOTATION_LEVEL=L1 Rscript hymy_distance_correlation_v2.R
#'
#' Author: CMM Project
#' =============================================================================

# =============================================================================
# Setup
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(RANN)
  library(parallel)
  library(meta)
  library(viridis)
  library(scales)
  library(pheatmap)
})

# Set seed for reproducibility
set.seed(42)

# Source utilities
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  return(getwd())
}
script_dir <- get_script_dir()
source(file.path(script_dir, "utils.R"))

# =============================================================================
# Configuration
# =============================================================================

# k-NN distance parameter (default k=1 = same as v1; k=5 for sanity check)
K_NEIGHBORS <- as.integer(Sys.getenv("K_NEIGHBORS", unset = "1"))

ANALYSIS_NAME <- paste0("hymy_distance_correlation_v2",
                        if (K_NEIGHBORS > 1) paste0("_k", K_NEIGHBORS) else "")

# Inherited from config.R (via utils.R): QUERY_CELLTYPE, CELLTYPE_COL, OUTPUT_SUFFIX, QUERY_LABEL
OUTPUT_BASE <- file.path(OUTPUT_ROOT, ANALYSIS_NAME)

ensure_dir(OUTPUT_BASE)
ensure_dir(file.path(OUTPUT_BASE, "per_celltype"))
ensure_dir(file.path(OUTPUT_BASE, "summary"))
ensure_dir(file.path(OUTPUT_BASE, "plots"))
ensure_dir(file.path(OUTPUT_BASE, "plots", "forest_plots"))

# Statistical parameters
MIN_CELLS_PER_SAMPLE <- 30       # Minimum cells of target type per sample
MIN_EXPR_PCT <- 0.01             # Minimum fraction of cells expressing per sample (1%)
MIN_EXPR_FLOOR <- 25             # Absolute floor for expressing cells per sample
MIN_EXPR_CELLS <- 5              # Minimum cells for stable GLM fit
FDR_THRESHOLD <- 0.05            # Significance threshold
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
N_PERMUTATIONS <- as.integer(Sys.getenv("N_PERMUTATIONS", unset = "500"))
PERM_TOP_N <- 100                # Run permutations on top N genes by effect size

# Distance parameters
MAX_DISTANCE_UM <- 200           # Maximum distance to consider (µm)

# Priority genes for permutation testing
PERM_PRIORITY_GENES <- c(
  # --- CC Chemokines ---
  "Ccl1", "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl6", "Ccl7", "Ccl8", "Ccl9",
  "Ccl11", "Ccl12", "Ccl17", "Ccl19", "Ccl20", "Ccl21a", "Ccl21b", "Ccl21c",
  "Ccl22", "Ccl24", "Ccl25", "Ccl27a", "Ccl27b", "Ccl28",
  # --- CXC Chemokines ---
  "Cxcl1", "Cxcl2", "Cxcl3", "Cxcl4", "Cxcl5", "Cxcl7", "Cxcl9", "Cxcl10",
  "Cxcl11", "Cxcl12", "Cxcl13", "Cxcl14", "Cxcl15", "Cxcl16", "Cxcl17",
  # --- CX3C / XC Chemokines ---
  "Cx3cl1", "Xcl1",
  # --- CC Chemokine Receptors ---
  "Ccr1", "Ccr2", "Ccr3", "Ccr4", "Ccr5", "Ccr6", "Ccr7", "Ccr8", "Ccr9", "Ccr10",
  # --- CXC Chemokine Receptors ---
  "Cxcr1", "Cxcr2", "Cxcr3", "Cxcr4", "Cxcr5", "Cxcr6",
  # --- CX3C / XC / Atypical Chemokine Receptors ---
  "Cx3cr1", "Xcr1", "Ackr1", "Ackr2", "Ackr3", "Ackr4",
  # --- Interleukins ---
  "Il1a", "Il1b", "Il1rn", "Il2", "Il3", "Il4", "Il5", "Il6", "Il7", "Il9",
  "Il10", "Il11", "Il12a", "Il12b", "Il13", "Il14", "Il15", "Il16",
  "Il17a", "Il17b", "Il17c", "Il17d", "Il17f",
  "Il18", "Il19", "Il20", "Il21", "Il22", "Il23a", "Il24", "Il25",
  "Il27", "Il31", "Il33", "Il34",
  "Il36a", "Il36b", "Il36g", "Il36rn",
  "Il1f5", "Il1f9", "Il1f10",
  # --- Interleukin Receptors ---
  "Il1r1", "Il1r2", "Il1rl1", "Il1rl2", "Il1rap",
  "Il2ra", "Il2rb", "Il2rg",
  "Il3ra", "Il4ra", "Il5ra", "Il6ra", "Il6st", "Il7r", "Il9r",
  "Il10ra", "Il10rb", "Il11ra1", "Il12rb1", "Il12rb2",
  "Il13ra1", "Il13ra2", "Il15ra",
  "Il17ra", "Il17rb", "Il17rc", "Il17rd", "Il17re",
  "Il18r1", "Il18rap", "Il20ra", "Il20rb", "Il21r",
  "Il22ra1", "Il22ra2", "Il23r", "Il27ra", "Il31ra",
  # --- Interferons ---
  "Ifna1", "Ifna2", "Ifna4", "Ifna5", "Ifna6", "Ifna7", "Ifna9",
  "Ifna11", "Ifna12", "Ifna13", "Ifna14", "Ifnab",
  "Ifnb1", "Ifne", "Ifnk", "Ifng", "Ifnl2", "Ifnl3",
  # --- Interferon Receptors ---
  "Ifnar1", "Ifnar2", "Ifngr1", "Ifngr2", "Ifnlr1",
  # --- TNF Superfamily Ligands ---
  "Tnf", "Lta", "Ltb", "Fasl", "Cd40lg",
  "Tnfsf4", "Tnfsf8", "Tnfsf9", "Tnfsf10", "Tnfsf11", "Tnfsf12",
  "Tnfsf13", "Tnfsf13b", "Tnfsf14", "Tnfsf15", "Tnfsf18",
  # --- TNF Receptor Superfamily ---
  "Tnfrsf1a", "Tnfrsf1b", "Tnfrsf4", "Fas", "Cd40",
  "Tnfrsf8", "Tnfrsf9", "Tnfrsf10b", "Tnfrsf11a", "Tnfrsf11b",
  "Tnfrsf12a", "Tnfrsf13b", "Tnfrsf13c", "Tnfrsf14",
  "Tnfrsf17", "Tnfrsf18", "Tnfrsf19", "Tnfrsf21", "Tnfrsf25",
  # --- Colony Stimulating Factors + Receptors ---
  "Csf1", "Csf2", "Csf3", "Csf1r", "Csf2ra", "Csf2rb", "Csf2rb2", "Csf3r",
  # --- TGF-beta Family + Receptors ---
  "Tgfb1", "Tgfb2", "Tgfb3", "Tgfbr1", "Tgfbr2", "Tgfbr3",
  # --- VEGF Family + Receptors ---
  "Vegfa", "Vegfb", "Vegfc", "Vegfd", "Flt1", "Kdr", "Flt4", "Nrp1", "Nrp2",
  # --- gp130 Family (IL-6 related) ---
  "Lif", "Osm", "Cntf", "Lifr", "Osmr", "Cntfr",
  # --- Other Cytokines ---
  "Tslp", "Mif", "Spp1", "Kitl", "Kit", "Epo", "Epor", "Thpo", "Mpl"
)

# Target cell types: user-specified or auto-detect (populated after data load)
TARGET_CELLTYPES_ENV <- Sys.getenv("TARGET_CELLTYPES", unset = "")
if (nchar(TARGET_CELLTYPES_ENV) > 0) {
  target_names <- trimws(strsplit(TARGET_CELLTYPES_ENV, ",")[[1]])
} else {
  target_names <- NULL  # Will be populated after data load
}

# =============================================================================
# Cell Type Selection (for SLURM array job parallelization)
# =============================================================================
# Note: CELLTYPE_INDEX selection is applied AFTER data load and auto-detection.
# See below where TARGET_CELLTYPES is finalized.

CELLTYPE_INDEX <- as.integer(Sys.getenv("CELLTYPE_INDEX", unset = "0"))

# Positive control genes (for validation)
POSITIVE_CONTROLS <- c("Csf3", "Il33", "Cxcl12")

message(strrep("=", 70))
message(paste0(QUERY_LABEL, " Distance Correlation Analysis v2 (Poisson GLM)"))
message(strrep("=", 70))
message("Annotation level: ", ANNOTATION_LEVEL)
message("Query cell type: ", QUERY_CELLTYPE)
message("K neighbors: ", K_NEIGHBORS)
message("Model: Poisson GLM (with cell size offset)")
message("Output directory: ", OUTPUT_BASE)
message("N cores: ", N_CORES)

# =============================================================================
# Helper Functions
# =============================================================================

#' Fit Poisson GLM of expression counts vs distance (per sample)
#'
#' Models counts ~ distance + offset(log(total_counts)) using Poisson regression.
#' The offset accounts for cell size (total transcripts per cell).
#'
#' @param count_vec Vector of raw integer counts for a single gene
#' @param dist_vec Vector of distances to nearest query cell
#' @param total_counts Vector of total counts per cell (for offset)
#' @param min_cells Minimum cells for analysis
#' @return List with coef (log-rate per µm), se, n_cells, pval, dispersion
fit_poisson <- function(count_vec, dist_vec, total_counts, min_cells = MIN_EXPR_CELLS) {

  # Remove NAs and invalid values
  valid_idx <- !is.na(count_vec) & !is.na(dist_vec) &
               is.finite(count_vec) & is.finite(dist_vec) &
               !is.na(total_counts) & total_counts > 0
  count_vec <- count_vec[valid_idx]
  dist_vec <- dist_vec[valid_idx]
  log_total <- log(total_counts[valid_idx])

  na_result <- list(coef = NA_real_, se = NA_real_, n_cells = length(count_vec),
                    pval = NA_real_, dispersion = NA_real_)

  # Need enough cells
  if (length(count_vec) < min_cells) return(na_result)

  # Need some non-zero counts for the model to be meaningful
  if (sum(count_vec > 0) < min_cells) return(na_result)

  # Poisson GLM with offset for cell size
  fit <- tryCatch({
    suppressWarnings(glm(count_vec ~ dist_vec + offset(log_total), family = poisson))
  }, error = function(e) NULL)

  if (is.null(fit) || !fit$converged) return(na_result)

  coef_summary <- summary(fit)$coefficients
  if (!"dist_vec" %in% rownames(coef_summary)) return(na_result)

  # Overdispersion diagnostic: residual deviance / residual df
  # Values >> 1 suggest negative binomial would be more appropriate
  dispersion <- fit$deviance / fit$df.residual

  list(
    coef = coef_summary["dist_vec", "Estimate"],      # log-rate change per µm
    se = coef_summary["dist_vec", "Std. Error"],
    n_cells = length(count_vec),
    pval = coef_summary["dist_vec", "Pr(>|z|)"],      # Wald z-test
    dispersion = dispersion
  )
}


#' Run meta-analysis across samples
#'
#' Uses random-effects meta-analysis to combine log-rate coefficients
#'
#' @param coefs Vector of log-rate coefficient estimates
#' @param ses Vector of standard errors
#' @param sample_ids Vector of sample identifiers
#' @return List with combined effect, se, pval, i2, n_samples
run_meta_analysis <- function(coefs, ses, sample_ids) {
  # Remove NAs
  valid_idx <- !is.na(coefs) & !is.na(ses) & ses > 0

  if (sum(valid_idx) < 2) {
    return(list(
      combined_coef = NA_real_,
      combined_se = NA_real_,
      pval = NA_real_,
      i2 = NA_real_,
      n_samples = sum(valid_idx)
    ))
  }

  coefs <- coefs[valid_idx]
  ses <- ses[valid_idx]
  sample_ids <- sample_ids[valid_idx]

  meta_result <- tryCatch({
    meta::metagen(
      TE = coefs,
      seTE = ses,
      studlab = sample_ids,
      random = TRUE,
      method.tau = "REML"
    )
  }, error = function(e) NULL)

  if (is.null(meta_result)) {
    return(list(
      combined_coef = NA_real_,
      combined_se = NA_real_,
      pval = NA_real_,
      i2 = NA_real_,
      n_samples = length(coefs)
    ))
  }

  list(
    combined_coef = meta_result$TE.random,
    combined_se = meta_result$seTE.random,
    pval = meta_result$pval.random,
    i2 = meta_result$I2,
    n_samples = length(coefs)
  )
}


#' Run permutation test for a single gene (Poisson version)
#'
#' Shuffles query cell labels and recalculates distance-expression coefficient.
#'
#' @param count_vec Count vector for target cells
#' @param coords_target Coordinates of target cells
#' @param coords_all Coordinates of ALL cells
#' @param sample_ids_target Sample IDs for target cells
#' @param n_perms Number of permutations
#' @param observed_coef The observed combined coefficient
#' @param sample_ids_all Sample IDs for all cells
#' @param query_per_sample Named vector with query cell counts per sample
#' @param k_neighbors Number of nearest neighbors for distance
#' @param total_counts_target Total counts per target cell (for offset)
#' @return List with null_coefs vector and empirical pval
run_permutation_test <- function(count_vec, coords_target, coords_all,
                                  sample_ids_target, n_perms,
                                  observed_coef, sample_ids_all,
                                  query_per_sample, k_neighbors,
                                  total_counts_target) {

  null_coefs <- numeric(n_perms)
  unique_samples <- names(query_per_sample)

  for (i in seq_len(n_perms)) {
    # STRATIFIED SAMPLING: Sample pseudo-query cells WITHIN each sample
    pseudo_query_coords_list <- lapply(unique_samples, function(samp) {
      samp_mask <- sample_ids_all == samp
      samp_coords <- coords_all[samp_mask, , drop = FALSE]
      n_to_sample <- query_per_sample[samp]

      if (n_to_sample > 0 && nrow(samp_coords) >= n_to_sample) {
        pseudo_idx <- sample(nrow(samp_coords), n_to_sample)
        samp_coords[pseudo_idx, , drop = FALSE]
      } else {
        matrix(nrow = 0, ncol = 2)
      }
    })

    pseudo_query_coords <- do.call(rbind, pseudo_query_coords_list)

    if (nrow(pseudo_query_coords) < 5) {
      null_coefs[i] <- NA
      next
    }

    # Calculate distances to pseudo-query cells
    effective_k <- min(k_neighbors, nrow(pseudo_query_coords))
    nn_result <- RANN::nn2(pseudo_query_coords, coords_target, k = effective_k)
    if (effective_k == 1) {
      perm_distances <- pmin(as.vector(nn_result$nn.dists), MAX_DISTANCE_UM)
    } else {
      perm_distances <- pmin(rowMeans(nn_result$nn.dists), MAX_DISTANCE_UM)
    }

    # Calculate Poisson coefficients per sample using inverse-variance weighting
    coefs <- numeric(length(unique_samples))
    ses <- numeric(length(unique_samples))

    for (j in seq_along(unique_samples)) {
      samp <- unique_samples[j]
      idx <- which(sample_ids_target == samp)

      if (length(idx) >= MIN_CELLS_PER_SAMPLE) {
        samp_counts <- count_vec[idx]
        samp_dist <- perm_distances[idx]
        samp_log_total <- log(total_counts_target[idx])

        if (length(samp_counts) >= MIN_CELLS_PER_SAMPLE && sum(samp_counts > 0) >= MIN_EXPR_CELLS) {
          fit <- tryCatch({
            suppressWarnings(glm(samp_counts ~ samp_dist + offset(samp_log_total), family = poisson))
          }, error = function(e) NULL)

          if (!is.null(fit) && fit$converged) {
            coef_summary <- summary(fit)$coefficients
            if (nrow(coef_summary) >= 2) {
              coefs[j] <- coef_summary[2, "Estimate"]
              ses[j] <- coef_summary[2, "Std. Error"]
            }
          }
        }
      }
    }

    # Simple inverse-variance weighted mean (faster than full meta-analysis)
    valid <- !is.na(coefs) & !is.na(ses) & ses > 0
    if (sum(valid) >= 2) {
      weights <- 1 / (ses[valid]^2)
      null_coefs[i] <- sum(weights * coefs[valid]) / sum(weights)
    } else {
      null_coefs[i] <- NA
    }
  }

  # Calculate empirical p-value (two-sided)
  null_coefs <- null_coefs[!is.na(null_coefs)]
  if (length(null_coefs) < 10) {
    return(list(null_coefs = null_coefs, perm_pval = NA_real_))
  }

  perm_pval <- (sum(abs(null_coefs) >= abs(observed_coef)) + 1) / (length(null_coefs) + 1)

  list(null_coefs = null_coefs, perm_pval = perm_pval)
}


#' Run permutation tests for multiple genes in parallel (Poisson version)
run_permutation_tests_parallel <- function(genes, count_matrix, target_barcodes,
                                            coords_target, coords_all,
                                            sample_ids_target, sample_ids_all,
                                            query_per_sample, observed_coefs,
                                            n_perms, n_cores, k_neighbors,
                                            total_counts_target) {

  message(sprintf("    Running %d permutations for %d genes using %d cores...",
                  n_perms, length(genes), n_cores))

  if (.Platform$OS.type == "unix") {
    results <- parallel::mclapply(genes, function(g) {
      count_vec <- as.numeric(count_matrix[g, target_barcodes])
      obs_coef <- observed_coefs[g]

      perm_result <- run_permutation_test(
        count_vec = count_vec,
        coords_target = coords_target,
        coords_all = coords_all,
        sample_ids_target = sample_ids_target,
        n_perms = n_perms,
        observed_coef = obs_coef,
        sample_ids_all = sample_ids_all,
        query_per_sample = query_per_sample,
        k_neighbors = k_neighbors,
        total_counts_target = total_counts_target
      )

      data.table(gene = g, perm_pval = perm_result$perm_pval)
    }, mc.cores = n_cores)
  } else {
    message("    (Windows detected - running sequentially)")
    results <- lapply(genes, function(g) {
      count_vec <- as.numeric(count_matrix[g, target_barcodes])
      obs_coef <- observed_coefs[g]

      perm_result <- run_permutation_test(
        count_vec = count_vec,
        coords_target = coords_target,
        coords_all = coords_all,
        sample_ids_target = sample_ids_target,
        n_perms = n_perms,
        observed_coef = obs_coef,
        sample_ids_all = sample_ids_all,
        query_per_sample = query_per_sample,
        k_neighbors = k_neighbors,
        total_counts_target = total_counts_target
      )

      data.table(gene = g, perm_pval = perm_result$perm_pval)
    })
  }

  rbindlist(results)
}


#' Classify decay pattern by fitting multiple Poisson models
#'
#' Uses Poisson regression on raw counts with offset (consistent with main analysis).
#'
#' @param count_vec Vector of raw counts
#' @param dist_vec Vector of distances
#' @param total_counts Vector of total counts per cell (for offset)
#' @return Character string: "linear", "exponential", "step_10um", "step_25um", "step_50um", or "none"
classify_decay_pattern <- function(count_vec, dist_vec, total_counts) {

  valid_idx <- !is.na(count_vec) & !is.na(dist_vec) &
               is.finite(count_vec) & is.finite(dist_vec) &
               !is.na(total_counts) & total_counts > 0
  count_vec <- count_vec[valid_idx]
  dist_vec <- dist_vec[valid_idx]
  log_total <- log(total_counts[valid_idx])

  if (length(count_vec) < 30) return("insufficient_data")
  if (sum(count_vec > 0) < 5) return("no_variation")

  # Fit linear Poisson model with offset
  fit_linear <- tryCatch({
    fit <- suppressWarnings(glm(count_vec ~ dist_vec + offset(log_total), family = poisson))
    if (fit$converged) fit else NULL
  }, error = function(e) NULL)

  # Fit step models at different thresholds with offset
  fit_step_10 <- tryCatch({
    fit <- suppressWarnings(glm(count_vec ~ I(dist_vec < 10) + offset(log_total), family = poisson))
    if (fit$converged) fit else NULL
  }, error = function(e) NULL)

  fit_step_25 <- tryCatch({
    fit <- suppressWarnings(glm(count_vec ~ I(dist_vec < 25) + offset(log_total), family = poisson))
    if (fit$converged) fit else NULL
  }, error = function(e) NULL)

  fit_step_50 <- tryCatch({
    fit <- suppressWarnings(glm(count_vec ~ I(dist_vec < 50) + offset(log_total), family = poisson))
    if (fit$converged) fit else NULL
  }, error = function(e) NULL)

  # Fit exponential decay (approximate via binned mean rates)
  fit_exp <- tryCatch({
    bins <- cut(dist_vec, breaks = seq(0, max(dist_vec) + 10, by = 20), include.lowest = TRUE)
    bin_rates <- tapply(count_vec / exp(log_total), bins, mean)
    bin_mids <- tapply(dist_vec, bins, mean)

    valid_bins <- !is.na(bin_rates) & !is.na(bin_mids) & bin_rates > 0
    if (sum(valid_bins) < 3) return(NULL)

    bin_rates <- bin_rates[valid_bins]
    bin_mids <- bin_mids[valid_bins]
    bin_n <- tapply(count_vec, bins, length)[valid_bins]

    log_rates <- log(bin_rates)

    nls(log_rates ~ a * exp(-b * bin_mids) + c,
        start = list(a = max(log_rates) - min(log_rates), b = 0.02, c = min(log_rates)),
        weights = bin_n,
        control = list(maxiter = 100, warnOnly = TRUE))
  }, error = function(e) NULL, warning = function(w) NULL)

  aics <- c(
    linear = if (!is.null(fit_linear)) AIC(fit_linear) else Inf,
    exponential = if (!is.null(fit_exp)) AIC(fit_exp) + 10 else Inf,
    step_10um = if (!is.null(fit_step_10)) AIC(fit_step_10) else Inf,
    step_25um = if (!is.null(fit_step_25)) AIC(fit_step_25) else Inf,
    step_50um = if (!is.null(fit_step_50)) AIC(fit_step_50) else Inf
  )

  if (all(is.infinite(aics))) return("none")
  names(aics)[which.min(aics)]
}


#' Calculate gradient score from Poisson regression coefficient
calculate_gradient_score <- function(coef) {
  coef
}


#' Create forest plot for meta-analysis
create_forest_plot <- function(coefs, ses, sample_ids, gene, cell_type, output_path) {
  valid_idx <- !is.na(coefs) & !is.na(ses) & ses > 0
  if (sum(valid_idx) < 2) return(invisible(NULL))

  coefs <- coefs[valid_idx]
  ses <- ses[valid_idx]
  sample_ids <- sample_ids[valid_idx]

  plot_data <- data.table(
    sample = sample_ids,
    coef = coefs,
    se = ses,
    lower = coefs - 1.96 * ses,
    upper = coefs + 1.96 * ses
  )

  meta_result <- run_meta_analysis(coefs, ses, sample_ids)

  if (!is.na(meta_result$combined_coef)) {
    meta_row <- data.table(
      sample = "Combined",
      coef = meta_result$combined_coef,
      se = meta_result$combined_se,
      lower = meta_result$combined_coef - 1.96 * meta_result$combined_se,
      upper = meta_result$combined_coef + 1.96 * meta_result$combined_se
    )
    plot_data <- rbind(plot_data, meta_row)
    plot_data[, is_combined := sample == "Combined"]
  } else {
    plot_data[, is_combined := FALSE]
  }

  plot_data[, sample := factor(sample, levels = rev(unique(sample)))]

  y_labels <- levels(plot_data$sample)
  y_faces <- ifelse(y_labels == "Combined", "bold", "plain")

  # Sign consistency annotation
  n_neg <- sum(coefs < 0)
  n_pos <- sum(coefs > 0)
  sign_text <- sprintf("%d/%d mice agree on sign", max(n_neg, n_pos), length(coefs))

  i2_display <- if (is.na(meta_result$i2)) "N/A" else sprintf("%.1f%%", meta_result$i2 * 100)
  pval_display <- if (is.na(meta_result$pval)) "N/A" else sprintf("%.2e", meta_result$pval)

  p <- ggplot(plot_data, aes(x = coef, y = sample)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
    geom_point(aes(shape = is_combined, size = is_combined)) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
    scale_size_manual(values = c("FALSE" = 3, "TRUE" = 5), guide = "none") +
    labs(
      x = "Log-rate coefficient (per \u00b5m)",
      y = NULL,
      title = sprintf("%s in %s (Poisson GLM)", gene, cell_type),
      subtitle = sprintf("p = %s, I\u00b2 = %s | %s", pval_display, i2_display, sign_text)
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(face = y_faces)
    )

  ggsave(output_path, p, width = 6, height = 4)
  invisible(p)
}


#' Create gradient volcano plot
create_gradient_volcano <- function(results, cell_type, output_path) {
  plot_data <- copy(results)

  plot_data[, neg_log10_fdr := -log10(fdr)]
  plot_data[neg_log10_fdr > 50, neg_log10_fdr := 50]

  plot_data[, significant := fdr < FDR_THRESHOLD]

  top_genes <- head(plot_data[significant == TRUE][order(fdr)], 20)

  max_score <- max(abs(plot_data$gradient_score), na.rm = TRUE) * 1.1

  plot_data[, decay_pattern := factor(decay_pattern,
    levels = c("linear", "exponential", "step_10um", "step_25um", "step_50um",
               "none", "no_variation", "insufficient_data", "not_significant", "undetermined"))]

  p <- ggplot(plot_data, aes(x = gradient_score, y = neg_log10_fdr)) +
    geom_point(aes(color = decay_pattern, size = significant), alpha = 0.6) +
    scale_color_brewer(palette = "Set2", name = "Decay Pattern") +
    scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2.5), guide = "none") +
    geom_hline(yintercept = -log10(FDR_THRESHOLD), linetype = "dashed",
               color = "grey40") +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey60") +
    xlim(-max_score, max_score) +
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5
    ) +
    labs(
      title = sprintf("Distance-Expression Analysis (Poisson GLM, k=%d): %s", K_NEIGHBORS, cell_type),
      subtitle = sprintf("%d genes significant (FDR < %.2f)",
                        sum(plot_data$significant, na.rm = TRUE), FDR_THRESHOLD),
      x = paste0("Log-rate coefficient (negative = ", QUERY_LABEL, "-induced)"),
      y = "-log10(FDR)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  ggsave(output_path, p, width = 10, height = 8)
  invisible(p)
}


#' Create decay example plots showing proportion expressing vs distance
#' (Kept as proportion for visual clarity, even though model uses counts)
create_decay_examples <- function(obj, distances, target_mask, genes,
                                   decay_patterns, output_path) {
  target_barcodes <- colnames(obj)[target_mask]
  target_distances <- distances[target_mask]

  available_genes <- intersect(genes, rownames(GetAssayData(obj, layer = "data")))
  if (length(available_genes) == 0) {
    message("    No genes available for decay examples plot")
    return(invisible(NULL))
  }
  if (length(available_genes) < length(genes)) {
    message(sprintf("    Note: %d/%d genes available for decay plot",
                    length(available_genes), length(genes)))
  }
  genes <- available_genes

  expr_data <- GetAssayData(obj, layer = "data")[genes, target_barcodes, drop = FALSE]

  plot_list <- lapply(genes, function(g) {
    expr_vec <- as.numeric(expr_data[g, ])
    pattern <- if (g %in% names(decay_patterns)) decay_patterns[g] else "unknown"

    df <- data.table(
      distance = target_distances,
      expressing = as.integer(expr_vec > 0)
    )

    df[, dist_bin := cut(distance, breaks = seq(0, MAX_DISTANCE_UM, by = 10),
                          include.lowest = TRUE)]

    bin_stats <- df[, .(
      prop_expressing = mean(expressing),
      n_cells = .N,
      se = sqrt(mean(expressing) * (1 - mean(expressing)) / .N)
    ), by = dist_bin]

    bin_stats[, dist_mid := as.numeric(sub("\\(|\\[", "", sub(",.*", "", as.character(dist_bin)))) + 5]
    bin_stats <- bin_stats[!is.na(dist_mid)]

    ggplot(bin_stats, aes(x = dist_mid, y = prop_expressing)) +
      geom_ribbon(aes(ymin = pmax(0, prop_expressing - 1.96 * se),
                      ymax = pmin(1, prop_expressing + 1.96 * se)),
                  fill = "#E74C3C", alpha = 0.2) +
      geom_line(color = "#E74C3C", linewidth = 1) +
      geom_point(aes(size = n_cells), color = "#E74C3C", alpha = 0.7) +
      scale_size_continuous(range = c(1, 4), guide = "none") +
      ylim(0, min(1, max(bin_stats$prop_expressing, na.rm = TRUE) * 1.2)) +
      labs(
        title = g,
        subtitle = sprintf("Pattern: %s", pattern),
        x = paste0("Distance to ", QUERY_LABEL, " (um)"),
        y = "Proportion Expressing"
      ) +
      theme_bw(base_size = 10) +
      theme(plot.title = element_text(face = "bold", size = 11))
  })

  combined <- wrap_plots(plot_list, ncol = 3)
  ggsave(output_path, combined, width = 12, height = 4 * ceiling(length(genes) / 3))
  invisible(combined)
}


#' Create coefficient strip plot for per-sample reproducibility
#'
#' Shows per-sample coefficients as points with SE bars for top significant genes.
#' Quick visual of whether all mice agree on direction.
create_coefficient_strips <- function(coef_results, meta_results, cell_type, output_path,
                                       n_top = 20) {
  sig_genes <- meta_results[fdr < FDR_THRESHOLD][order(fdr)]$gene
  if (length(sig_genes) == 0) {
    message("    No significant genes for coefficient strip plot")
    return(invisible(NULL))
  }
  genes_to_plot <- head(sig_genes, n_top)

  plot_data <- coef_results[gene %in% genes_to_plot & !is.na(coef)]

  if (nrow(plot_data) == 0) return(invisible(NULL))

  # Order genes by combined coefficient
  gene_order <- meta_results[gene %in% genes_to_plot][order(combined_coef)]$gene
  plot_data[, gene := factor(gene, levels = gene_order)]

  # Add sign consistency annotation from meta_results
  consistency <- meta_results[gene %in% genes_to_plot, .(gene, sign_consistency)]
  plot_data <- merge(plot_data, consistency, by = "gene", all.x = TRUE)

  p <- ggplot(plot_data, aes(x = coef, y = gene, color = sample_id)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = coef - 1.96 * se, xmax = coef + 1.96 * se),
                   height = 0.2, alpha = 0.5) +
    geom_point(size = 2.5, alpha = 0.8) +
    scale_color_brewer(palette = "Set1", name = "Sample") +
    labs(
      x = "Log-rate coefficient (per \u00b5m)",
      y = NULL,
      title = sprintf("Per-Sample Coefficients: %s (Poisson GLM)", cell_type),
      subtitle = sprintf("Top %d significant genes | Points = individual mice",
                         length(genes_to_plot))
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggsave(output_path, p, width = 8, height = 0.4 * length(genes_to_plot) + 2)
  invisible(p)
}


# =============================================================================
# Main Analysis Function
# =============================================================================

#' Run Poisson distance correlation analysis for a single cell type
run_celltype_analysis <- function(obj, cell_data, cell_type_name, target_types,
                                   query_celltype, celltype_col, output_dir,
                                   coords_all, sample_ids_all, query_per_sample,
                                   count_matrix, total_counts_all) {

  message("\n", strrep("-", 60))
  message("Analyzing: ", cell_type_name)
  message(strrep("-", 60))

  ensure_dir(output_dir)

  # Identify target cells (filtered data — condition already applied)
  cell_data[, is_target := get(celltype_col) %in% target_types]
  target_data <- cell_data[is_target == TRUE]

  message("  Target cells: ", nrow(target_data))

  if (nrow(target_data) < MIN_CELLS_PER_SAMPLE * 2) {
    message("  Insufficient cells for analysis")
    return(NULL)
  }

  # Get unique samples
  samples <- unique(target_data[[SAMPLE_COL]])
  message("  Samples: ", length(samples))

  # Check each sample has enough cells
  sample_counts <- target_data[, .N, by = c(SAMPLE_COL)]
  valid_samples <- sample_counts[N >= MIN_CELLS_PER_SAMPLE][[SAMPLE_COL]]
  message("  Valid samples (>= ", MIN_CELLS_PER_SAMPLE, " cells): ", length(valid_samples))

  if (length(valid_samples) < 2) {
    message("  Need at least 2 valid samples for meta-analysis")
    return(NULL)
  }
  if (length(valid_samples) == 2) {
    message("  WARNING: Only 2 valid samples - meta-analysis will have low power")
  }

  # Get target cell barcodes
  target_barcodes <- target_data[get(SAMPLE_COL) %in% valid_samples]$barcode

  # Get raw counts for target cells
  target_counts <- count_matrix[, target_barcodes, drop = FALSE]

  # Filter genes with sufficient expressing cells PER SAMPLE (Two-Tier)
  target_valid <- target_data[get(SAMPLE_COL) %in% valid_samples]
  sample_ids_for_filter <- droplevels(as.factor(target_valid[[SAMPLE_COL]]))

  cells_per_sample <- table(sample_ids_for_filter)
  threshold_per_sample <- pmax(ceiling(cells_per_sample * MIN_EXPR_PCT), MIN_EXPR_FLOOR)
  message("  Per-sample expression thresholds (max(1%, ", MIN_EXPR_FLOOR, ")):")
  for (s in names(threshold_per_sample)) {
    message("    ", s, ": ", threshold_per_sample[s], " cells (of ", cells_per_sample[s], ")")
  }

  threshold_floor <- setNames(rep(MIN_EXPR_FLOOR, length(cells_per_sample)),
                              names(cells_per_sample))

  # Count expressing cells per gene per sample (count > 0 is same as binary filter)
  expressing_counts <- sapply(rownames(target_counts), function(g) {
    count_vec <- target_counts[g, ]
    tapply(count_vec > 0, sample_ids_for_filter, sum)
  })

  # --- Tier 1: Strict filter ---
  genes_pass_strict <- apply(expressing_counts, 2, function(counts) {
    all(counts >= threshold_per_sample[names(counts)], na.rm = TRUE)
  })
  genes_strict <- names(genes_pass_strict[genes_pass_strict])

  # --- Tier 2: Lenient filter (priority genes) ---
  priority_in_panel <- intersect(PERM_PRIORITY_GENES, rownames(target_counts))
  genes_pass_lenient <- apply(expressing_counts[, priority_in_panel, drop = FALSE], 2, function(counts) {
    sum(counts >= threshold_floor[names(counts)], na.rm = TRUE) >= 2
  })
  genes_lenient <- names(genes_pass_lenient[genes_pass_lenient])
  genes_lenient_only <- setdiff(genes_lenient, genes_strict)

  genes_to_analyze <- union(genes_strict, genes_lenient_only)
  message("  Genes passing strict filter: ", length(genes_strict))
  if (length(genes_lenient_only) > 0) {
    message("  Priority genes rescued by lenient filter: ", length(genes_lenient_only),
            " (", paste(head(genes_lenient_only, 10), collapse = ", "),
            if (length(genes_lenient_only) > 10) ", ..." else "", ")")
  }
  message("  Total genes to analyze: ", length(genes_to_analyze))

  if (length(genes_to_analyze) < 10) {
    message("  Too few genes passed filtering")
    return(NULL)
  }

  # Step 1: Calculate Poisson regression coefficients per sample for each gene
  message("  Step 1: Calculating per-sample Poisson coefficients...")

  coef_results <- rbindlist(lapply(genes_to_analyze, function(g) {
    gene_counts <- as.numeric(target_counts[g, target_barcodes])

    sample_results <- rbindlist(lapply(valid_samples, function(samp) {
      samp_idx <- which(target_valid[[SAMPLE_COL]] == samp)
      if (length(samp_idx) < MIN_CELLS_PER_SAMPLE) {
        return(data.table(
          gene = g, sample_id = samp,
          coef = NA_real_, se = NA_real_,
          n_cells = length(samp_idx), pval = NA_real_,
          dispersion = NA_real_
        ))
      }

      samp_counts <- gene_counts[samp_idx]
      samp_dist <- target_valid[samp_idx]$dist_to_query
      samp_total <- total_counts_all[target_barcodes[samp_idx]]

      fit_result <- fit_poisson(samp_counts, samp_dist, samp_total)

      data.table(
        gene = g, sample_id = samp,
        coef = fit_result$coef, se = fit_result$se,
        n_cells = fit_result$n_cells, pval = fit_result$pval,
        dispersion = fit_result$dispersion
      )
    }))

    sample_results
  }), fill = TRUE)

  # Save per-sample coefficients
  fwrite(coef_results, file.path(output_dir, "coef_per_sample.csv"))
  message("  Saved: coef_per_sample.csv")

  # Step 2: Meta-analysis across samples
  message("  Step 2: Running meta-analysis...")

  meta_results <- rbindlist(lapply(genes_to_analyze, function(g) {
    gene_data <- coef_results[gene == g]

    meta_result <- run_meta_analysis(
      coefs = gene_data$coef,
      ses = gene_data$se,
      sample_ids = gene_data[["sample_id"]]
    )

    # Calculate expression statistics from counts
    gene_counts <- as.numeric(target_counts[g, target_barcodes])

    # Per-sample sign consistency
    valid_coefs <- gene_data$coef[!is.na(gene_data$coef)]
    n_pos <- sum(valid_coefs > 0)
    n_neg <- sum(valid_coefs < 0)
    n_valid <- n_pos + n_neg
    sign_con <- if (n_valid > 0) max(n_pos, n_neg) / n_valid else NA_real_

    # Median dispersion across samples
    med_disp <- median(gene_data$dispersion, na.rm = TRUE)

    data.table(
      gene = g,
      combined_coef = meta_result$combined_coef,
      combined_se = meta_result$combined_se,
      pval = meta_result$pval,
      i2 = meta_result$i2,
      n_samples = meta_result$n_samples,
      mean_expr = mean(gene_counts, na.rm = TRUE),
      pct_expr = mean(gene_counts > 0, na.rm = TRUE),
      n_positive_samples = n_pos,
      n_negative_samples = n_neg,
      sign_consistency = sign_con,
      median_dispersion = med_disp
    )
  }), fill = TRUE)

  # Calculate FDR
  meta_results[, fdr := p.adjust(pval, method = "BH")]

  # Step 2b: Permutation testing
  message("  Step 2b: Running permutation tests...")

  top_by_effect <- meta_results[order(-abs(combined_coef))][1:min(PERM_TOP_N, .N)]$gene
  priority_in_data <- intersect(PERM_PRIORITY_GENES, genes_to_analyze)
  top_genes_for_perm <- unique(c(top_by_effect, priority_in_data))

  message(sprintf("    Permutation genes: %d top by effect + %d priority = %d total (after dedup)",
                  length(top_by_effect), length(priority_in_data), length(top_genes_for_perm)))

  if (length(top_genes_for_perm) > 0 && N_PERMUTATIONS > 0) {
    target_valid <- target_data[get(SAMPLE_COL) %in% valid_samples]
    coords_target <- as.matrix(target_valid[, ..coord_cols])

    observed_coefs <- setNames(
      meta_results[gene %in% top_genes_for_perm]$combined_coef,
      meta_results[gene %in% top_genes_for_perm]$gene
    )

    query_per_sample_valid <- query_per_sample[names(query_per_sample) %in% valid_samples]

    perm_results <- run_permutation_tests_parallel(
      genes = top_genes_for_perm,
      count_matrix = target_counts,
      target_barcodes = target_barcodes,
      coords_target = coords_target,
      coords_all = coords_all,
      sample_ids_target = target_valid[[SAMPLE_COL]],
      sample_ids_all = sample_ids_all,
      query_per_sample = query_per_sample_valid,
      observed_coefs = observed_coefs,
      n_perms = N_PERMUTATIONS,
      n_cores = N_CORES,
      k_neighbors = K_NEIGHBORS,
      total_counts_target = total_counts_all[target_barcodes]
    )

    meta_results <- merge(meta_results, perm_results, by = "gene", all.x = TRUE)
    setDT(meta_results)

    fwrite(perm_results, file.path(output_dir, "permutation_pvals.csv"))
    message("    Saved: permutation_pvals.csv")
    message("    Genes with perm_pval < 0.05: ", sum(perm_results$perm_pval < 0.05, na.rm = TRUE))
  } else {
    meta_results[, perm_pval := NA_real_]
    message("    Skipping permutation tests (N_PERMUTATIONS = ", N_PERMUTATIONS, ")")
  }

  # Step 3: Calculate gradient scores
  message("  Step 3: Calculating gradient scores...")
  meta_results[, gradient_score := calculate_gradient_score(combined_coef)]

  # Step 4: Classify decay patterns (for significant genes)
  message("  Step 4: Classifying decay patterns (per-sample majority vote)...")

  sig_genes <- meta_results[fdr < FDR_THRESHOLD]$gene
  message("  Significant genes: ", length(sig_genes))

  if (length(sig_genes) > 0) {
    target_valid <- target_data[get(SAMPLE_COL) %in% valid_samples]

    decay_patterns <- sapply(sig_genes, function(gene) {
      gene_counts <- as.numeric(target_counts[gene, target_barcodes])

      per_sample_patterns <- sapply(valid_samples, function(samp) {
        samp_idx <- which(target_valid[[SAMPLE_COL]] == samp)
        if (length(samp_idx) < MIN_CELLS_PER_SAMPLE) return(NA_character_)

        samp_counts <- gene_counts[samp_idx]
        samp_dist <- target_valid[samp_idx]$dist_to_query
        samp_total <- total_counts_all[target_barcodes[samp_idx]]

        tryCatch(
          classify_decay_pattern(samp_counts, samp_dist, samp_total),
          error = function(e) NA_character_
        )
      })

      per_sample_patterns <- per_sample_patterns[!is.na(per_sample_patterns)]
      if (length(per_sample_patterns) == 0) return("undetermined")

      pattern_counts <- table(per_sample_patterns)
      names(pattern_counts)[which.max(pattern_counts)]
    })

    meta_results[, decay_pattern := "not_significant"]
    meta_results[gene %in% sig_genes, decay_pattern := decay_patterns[match(gene, sig_genes)]]
  } else {
    meta_results[, decay_pattern := "not_significant"]
  }

  # Save results
  fwrite(meta_results, file.path(output_dir, "meta_analysis_results.csv"))
  message("  Saved: meta_analysis_results.csv")

  decay_summary <- meta_results[fdr < FDR_THRESHOLD, .N, by = decay_pattern]
  fwrite(decay_summary, file.path(output_dir, "decay_classification.csv"))

  gradient_results <- meta_results[, .(gene, gradient_score, combined_coef, fdr, decay_pattern,
                                        sign_consistency, median_dispersion)]
  gradient_results <- gradient_results[order(gradient_score)]
  fwrite(gradient_results, file.path(output_dir, "gradient_scores.csv"))

  # Step 5: Create visualizations
  message("  Step 5: Creating visualizations...")

  # Gradient volcano plot
  create_gradient_volcano(
    meta_results, cell_type_name,
    file.path(output_dir, "gradient_volcano.pdf")
  )

  # Forest plots for top genes
  top_sig_genes <- head(meta_results[fdr < FDR_THRESHOLD][order(fdr)]$gene, 20)

  if (length(top_sig_genes) > 0) {
    forest_dir <- file.path(output_dir, "forest_plots")
    ensure_dir(forest_dir)

    for (g in top_sig_genes) {
      gene_data <- coef_results[gene == g]
      create_forest_plot(
        coefs = gene_data$coef,
        ses = gene_data$se,
        sample_ids = gene_data[["sample_id"]],
        gene = g,
        cell_type = cell_type_name,
        output_path = file.path(forest_dir, sprintf("%s_forest.pdf", g))
      )
    }

    # Decay example plots
    decay_patterns_named <- setNames(
      meta_results[gene %in% head(top_sig_genes, 10)]$decay_pattern,
      meta_results[gene %in% head(top_sig_genes, 10)]$gene
    )

    create_decay_examples(
      obj = obj,
      distances = cell_data$dist_to_query,
      target_mask = cell_data[, is_target],
      genes = head(top_sig_genes, 10),
      decay_patterns = decay_patterns_named,
      output_path = file.path(output_dir, "decay_examples.pdf")
    )

    # Coefficient strip plot (new in v2)
    create_coefficient_strips(
      coef_results = coef_results,
      meta_results = meta_results,
      cell_type = cell_type_name,
      output_path = file.path(output_dir, "coefficient_strips.pdf")
    )
  }

  # Add cell type to results
  meta_results[, cell_type := cell_type_name]

  message("  Analysis complete for ", cell_type_name)

  return(meta_results)
}


# =============================================================================
# Main Execution
# =============================================================================

message("\n", strrep("=", 70))
message("Loading Data")
message(strrep("=", 70))

# Load Seurat object
obj <- load_seurat()

# Merge HyMy annotations if using legacy HyMy annotation level
if (ANNOTATION_LEVEL == "HyMy") {
  message("Merging HyMy annotations...")
  obj <- merge_hymy_annotations(obj)
}

# Verify counts layer exists
counts_check <- tryCatch({
  test_counts <- GetAssayData(obj, layer = "counts")
  message("Counts layer verified: ", nrow(test_counts), " genes x ", ncol(test_counts), " cells")
  TRUE
}, error = function(e) {
  stop("FATAL: No 'counts' layer found in Seurat object. ",
       "Poisson GLM requires raw integer counts. Error: ", e$message)
})

# Get count matrix
message("Loading count matrix...")
count_matrix_full <- GetAssayData(obj, layer = "counts")

# Extract metadata as data.table
cell_data <- as.data.table(obj@meta.data, keep.rownames = "barcode")

# Resolve condition column
if (nchar(CONDITION_COL) > 0 && CONDITION_COL %in% names(cell_data)) {
  cell_data[, condition := get(CONDITION_COL)]
  message("Using condition column: ", CONDITION_COL)
} else if ("condition" %in% names(cell_data)) {
  message("Using existing 'condition' column")
} else if ("group" %in% names(cell_data)) {
  cell_data[, condition := group]
  message("Aliasing 'group' -> 'condition'")
} else {
  cell_data[, condition := "all"]
  message("No condition column found; analyzing all samples")
}

if (!CELLTYPE_COL %in% names(cell_data)) {
  stop("Cell type column not found: ", CELLTYPE_COL)
}

message("\nData summary (before filtering):")
message("  Total cells: ", nrow(cell_data))
message("  Samples: ", length(unique(cell_data[[SAMPLE_COL]])))
message("  Conditions: ", paste(unique(cell_data$condition), collapse = ", "))

# =============================================================================
# Filter by condition if specified
# =============================================================================

if (nchar(CONDITION_VAL) > 0) {
  message("\nFiltering to condition == '", CONDITION_VAL, "'...")
  cell_data <- cell_data[condition == CONDITION_VAL]
} else {
  message("\nNo condition filter; analyzing all ", uniqueN(cell_data$condition), " conditions")
}

message("Data summary (after filtering):")
message("  Total cells: ", nrow(cell_data))
message("  Samples: ", length(unique(cell_data[[SAMPLE_COL]])))
message("  Cell types: ", length(unique(cell_data[[CELLTYPE_COL]])))

# Subset count matrix to TDLN barcodes
count_matrix_tdln <- count_matrix_full[, cell_data$barcode, drop = FALSE]

# Compute total counts per cell (for Poisson offset)
total_counts_tdln <- colSums(count_matrix_full[, cell_data$barcode, drop = FALSE])
total_counts_tdln <- setNames(as.numeric(total_counts_tdln), cell_data$barcode)
message("Total counts per cell: median=", round(median(total_counts_tdln)),
        ", range=[", round(min(total_counts_tdln)), "-", round(max(total_counts_tdln)), "]")

# Free full count matrix
rm(count_matrix_full)
gc()

# =============================================================================
# Calculate distances to query cells
# =============================================================================

message("\n", strrep("=", 70))
message("Calculating Distances to Query Cells (k=", K_NEIGHBORS, ")")
message(strrep("=", 70))

coord_cols <- get_coord_columns(cell_data)
coords <- as.matrix(cell_data[, ..coord_cols])
sample_ids_all <- cell_data[[SAMPLE_COL]]

query_mask <- cell_data[[CELLTYPE_COL]] == QUERY_CELLTYPE
n_query <- sum(query_mask)
message("Query cells (", QUERY_CELLTYPE, "): ", n_query)

if (n_query < 10) {
  stop("Too few query cells for analysis (", n_query, ")")
}

# Calculate query cells per sample (for stratified permutation)
query_per_sample <- cell_data[query_mask == TRUE, .N, by = c(SAMPLE_COL)]
query_per_sample <- setNames(query_per_sample$N, query_per_sample[[SAMPLE_COL]])
message("Query cells per sample:")
for (samp in names(query_per_sample)) {
  message("  ", samp, ": ", query_per_sample[samp])
}

query_coords <- coords[query_mask, , drop = FALSE]

# Calculate distance from each cell to nearest k query cells
effective_k <- min(K_NEIGHBORS, nrow(query_coords))
message("Computing ", effective_k, "-NN distances...")
nn_result <- nn2(query_coords, coords, k = effective_k)

if (effective_k == 1) {
  cell_data[, dist_to_query := as.vector(nn_result$nn.dists)]
} else {
  cell_data[, dist_to_query := rowMeans(nn_result$nn.dists)]
}

# Cap distances at maximum
cell_data[dist_to_query > MAX_DISTANCE_UM, dist_to_query := MAX_DISTANCE_UM]

message("Distance distribution:")
message("  Min: ", round(min(cell_data$dist_to_query), 1), " \u00b5m")
message("  Median: ", round(median(cell_data$dist_to_query), 1), " \u00b5m")
message("  Max: ", round(max(cell_data$dist_to_query), 1), " \u00b5m")

# =============================================================================
# Finalize Target Cell Types (auto-detect if not specified)
# =============================================================================

if (is.null(target_names)) {
  all_types <- unique(cell_data[[CELLTYPE_COL]])
  target_names <- setdiff(all_types, c(QUERY_CELLTYPE, NA_character_))
  target_names <- sort(target_names)  # deterministic order
  message("Auto-detected ", length(target_names), " target cell types: ",
          paste(target_names, collapse = ", "))
}
# Each cell type is its own target (no aggregation)
TARGET_CELLTYPES <- as.list(setNames(target_names, target_names))

# Apply CELLTYPE_INDEX selection for array jobs
if (CELLTYPE_INDEX > 0) {
  if (CELLTYPE_INDEX > length(target_names)) {
    stop("CELLTYPE_INDEX=", CELLTYPE_INDEX, " exceeds ", length(target_names), " target cell types")
  }
  selected <- target_names[CELLTYPE_INDEX]
  TARGET_CELLTYPES <- TARGET_CELLTYPES[selected]
  message("Selected cell type #", CELLTYPE_INDEX, ": ", selected)
}

# =============================================================================
# QC Diagnostics
# =============================================================================

message("\n", strrep("=", 70))
message("QC Diagnostics")
message(strrep("=", 70))

qc_dir <- file.path(OUTPUT_BASE, "qc")
ensure_dir(qc_dir)

# 1. Distance distribution histogram
message("Creating distance distribution histogram...")
p_dist <- ggplot(cell_data, aes(x = dist_to_query)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.7) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = 55, y = Inf, label = "50\u00b5m threshold", vjust = 2, hjust = 0, color = "red") +
  labs(
    title = sprintf("Distance to Nearest %s (k=%d)", QUERY_CELLTYPE, K_NEIGHBORS),
    x = "Distance (\u00b5m)",
    y = "Number of Cells"
  ) +
  theme_bw(base_size = 12)
ggsave(file.path(qc_dir, "distance_distribution.pdf"), p_dist, width = 8, height = 5)
message("  Saved: qc/distance_distribution.pdf")

# 2. Per-sample summary table
sample_summary <- cell_data[, .(
  n_total = .N,
  n_query = sum(get(CELLTYPE_COL) == QUERY_CELLTYPE),
  median_dist = median(dist_to_query),
  pct_within_50um = mean(dist_to_query <= 50) * 100,
  n_cell_types = uniqueN(get(CELLTYPE_COL))
), by = c(SAMPLE_COL, "condition")]
setorder(sample_summary, condition)
fwrite(sample_summary, file.path(qc_dir, "sample_summary.csv"))
message("  Saved: qc/sample_summary.csv")

message("\nPer-sample summary:")
print(sample_summary)

# 3. Cell type counts (filtered data)
filtered_celltype_counts <- cell_data[, .N, by = c(CELLTYPE_COL)]
setnames(filtered_celltype_counts, CELLTYPE_COL, "cell_type")
setorder(filtered_celltype_counts, -N)
fwrite(filtered_celltype_counts, file.path(qc_dir, "filtered_celltype_counts.csv"))
message("\nCell type counts (top 10):")
print(head(filtered_celltype_counts, 10))

# 4. Count matrix summary
message("\nCount matrix: ", nrow(count_matrix_tdln), " genes x ", ncol(count_matrix_tdln), " cells")

# =============================================================================
# Run Analysis for Each Target Cell Type
# =============================================================================

message("\n", strrep("=", 70))
message("Running Distance Correlation Analysis (Poisson GLM)")
message(strrep("=", 70))

all_results <- list()

for (ct_name in names(TARGET_CELLTYPES)) {
  ct_types <- TARGET_CELLTYPES[[ct_name]]
  ct_output_dir <- file.path(OUTPUT_BASE, "per_celltype", ct_name)

  result <- tryCatch({
    run_celltype_analysis(
      obj = obj,
      cell_data = cell_data,
      cell_type_name = ct_name,
      target_types = ct_types,
      query_celltype = QUERY_CELLTYPE,
      celltype_col = CELLTYPE_COL,
      output_dir = ct_output_dir,
      coords_all = coords,
      sample_ids_all = sample_ids_all,
      query_per_sample = query_per_sample,
      count_matrix = count_matrix_tdln,
      total_counts_all = total_counts_tdln
    )
  }, error = function(e) {
    message("  ERROR: ", e$message)
    NULL
  })

  if (!is.null(result)) {
    all_results[[ct_name]] <- result
  }
}

# =============================================================================
# Summary Statistics
# =============================================================================

message("\n", strrep("=", 70))
message("Creating Summary")
message(strrep("=", 70))

if (length(all_results) > 0) {
  combined_results <- rbindlist(all_results, fill = TRUE)

  fwrite(combined_results, file.path(OUTPUT_BASE, "summary", "all_genes_results.csv"))
  message("Saved: summary/all_genes_results.csv")

  top_genes <- combined_results[fdr < FDR_THRESHOLD][
    order(fdr), head(.SD, 50), by = cell_type
  ]
  fwrite(top_genes, file.path(OUTPUT_BASE, "summary", "top_gradient_genes.csv"))
  message("Saved: summary/top_gradient_genes.csv")

  decay_summary <- combined_results[fdr < FDR_THRESHOLD, .N,
                                     by = .(cell_type, decay_pattern)]
  decay_summary_wide <- dcast(decay_summary, cell_type ~ decay_pattern,
                               value.var = "N", fill = 0)
  fwrite(decay_summary_wide, file.path(OUTPUT_BASE, "summary", "decay_pattern_summary.csv"))
  message("Saved: summary/decay_pattern_summary.csv")

  # Overdispersion summary (new in v2)
  disp_summary <- combined_results[, .(
    median_dispersion = median(median_dispersion, na.rm = TRUE),
    pct_overdispersed = mean(median_dispersion > 2, na.rm = TRUE) * 100
  ), by = cell_type]
  fwrite(disp_summary, file.path(OUTPUT_BASE, "summary", "dispersion_summary.csv"))
  message("Saved: summary/dispersion_summary.csv")
  message("\nOverdispersion summary:")
  print(disp_summary)

  # Sign consistency summary (new in v2)
  sign_summary <- combined_results[fdr < FDR_THRESHOLD, .(
    n_sig = .N,
    pct_all_agree = mean(sign_consistency == 1, na.rm = TRUE) * 100,
    median_sign_consistency = median(sign_consistency, na.rm = TRUE)
  ), by = cell_type]
  fwrite(sign_summary, file.path(OUTPUT_BASE, "summary", "sign_consistency_summary.csv"))
  message("Saved: summary/sign_consistency_summary.csv")
  message("\nSign consistency among significant genes:")
  print(sign_summary)

  # Summary heatmap (only when multiple cell types)
  n_celltypes <- length(unique(combined_results$cell_type))

  if (n_celltypes < 2) {
    message("\nSkipping summary heatmap (only 1 cell type analyzed)")
    message("  Run merge after all array jobs complete for cross-cell-type summary")
  } else {
    message("\nCreating summary heatmap...")

    top_all <- combined_results[fdr < FDR_THRESHOLD][order(fdr)][1:min(100, .N)]

    if (nrow(top_all) > 0) {
      heatmap_data <- dcast(combined_results[gene %in% top_all$gene],
                            gene ~ cell_type, value.var = "gradient_score")

      heatmap_matrix <- as.matrix(heatmap_data[, -1, with = FALSE])
      rownames(heatmap_matrix) <- heatmap_data$gene
      heatmap_matrix[is.na(heatmap_matrix)] <- 0

      row_annot <- data.frame(
        decay = top_all$decay_pattern[match(rownames(heatmap_matrix), top_all$gene)]
      )
      rownames(row_annot) <- rownames(heatmap_matrix)

      pdf(file.path(OUTPUT_BASE, "plots", "heatmap_top_gradient_genes.pdf"),
          width = 10, height = 12)

      max_val <- max(abs(heatmap_matrix), na.rm = TRUE)
      if (is.na(max_val) || max_val == 0) max_val <- 1

      cluster_rows <- nrow(heatmap_matrix) >= 2
      cluster_cols <- ncol(heatmap_matrix) >= 2

      pheatmap(
        heatmap_matrix,
        main = "Top Gradient Genes Across Cell Types (Poisson GLM)",
        color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
        breaks = seq(-max_val, max_val, length.out = 101),
        cluster_rows = cluster_rows,
        cluster_cols = cluster_cols,
        show_rownames = nrow(heatmap_matrix) <= 50,
        annotation_row = row_annot,
        fontsize = 8
      )

      dev.off()
      message("Saved: plots/heatmap_top_gradient_genes.pdf")
    }
  }

  # Decay pattern bar plot
  if (nrow(decay_summary) > 0) {
    p_decay <- ggplot(decay_summary, aes(x = cell_type, y = N, fill = decay_pattern)) +
      geom_col(position = "stack") +
      scale_fill_brewer(palette = "Set2", name = "Decay Pattern") +
      labs(
        title = "Decay Pattern Distribution by Cell Type (Poisson GLM)",
        x = NULL,
        y = "Number of Genes"
      ) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
      )

    ggsave(file.path(OUTPUT_BASE, "plots", "decay_pattern_barplot.pdf"),
           p_decay, width = 10, height = 6)
    message("Saved: plots/decay_pattern_barplot.pdf")
  }

  # Print summary statistics
  message("\n", strrep("=", 70))
  message("Analysis Summary")
  message(strrep("=", 70))

  for (ct in names(all_results)) {
    ct_result <- all_results[[ct]]
    n_sig <- sum(ct_result$fdr < FDR_THRESHOLD, na.rm = TRUE)
    n_negative <- sum(ct_result$fdr < FDR_THRESHOLD &
                       ct_result$gradient_score < 0, na.rm = TRUE)
    n_positive <- sum(ct_result$fdr < FDR_THRESHOLD &
                       ct_result$gradient_score > 0, na.rm = TRUE)

    message(sprintf("  %s: %d significant genes (%d %s-induced, %d %s-repressed)",
                    ct, n_sig, n_negative, QUERY_LABEL, n_positive, QUERY_LABEL))
  }

  # ==========================================================================
  # Positive Control Validation
  # ==========================================================================
  message("\n", strrep("-", 60))
  message("Positive Control Validation")
  message(strrep("-", 60))
  message("Expected: CSF3, IL33, CXCL12 should show NEGATIVE gradient in FRC/LEC")
  message(paste0("(negative = higher expression rate near ", QUERY_LABEL, " = ", QUERY_LABEL, "-induced)"))
  message("")

  for (ct in c("FRC", "LEC")) {
    if (ct %in% names(all_results)) {
      ct_result <- all_results[[ct]]
      for (ctrl_gene in POSITIVE_CONTROLS) {
        ctrl_row <- ct_result[gene == ctrl_gene]
        if (nrow(ctrl_row) > 0) {
          score <- round(ctrl_row$gradient_score, 5)
          fdr_val <- signif(ctrl_row$fdr, 3)
          sig_marker <- if (!is.na(fdr_val) && fdr_val < FDR_THRESHOLD) "*" else ""
          direction <- if (!is.na(score) && score < 0) paste0(QUERY_LABEL, "-induced") else "NOT induced"
          sign_con <- if (!is.na(ctrl_row$sign_consistency)) {
            sprintf("sign: %d/%d agree", ctrl_row$n_negative_samples + ctrl_row$n_positive_samples -
                      min(ctrl_row$n_negative_samples, ctrl_row$n_positive_samples),
                    ctrl_row$n_negative_samples + ctrl_row$n_positive_samples)
          } else "N/A"
          message(sprintf("  %s in %s: coef=%.5f, FDR=%.2e %s [%s] (%s)",
                          ctrl_gene, ct, score, fdr_val, sig_marker, direction, sign_con))
        } else {
          message(sprintf("  %s in %s: not found in results", ctrl_gene, ct))
        }
      }
    }
  }

  validation_results <- rbindlist(lapply(names(all_results), function(ct) {
    ct_result <- all_results[[ct]]
    ct_result[gene %in% POSITIVE_CONTROLS,
              .(cell_type = ct, gene, gradient_score, fdr, decay_pattern,
                sign_consistency, median_dispersion)]
  }), fill = TRUE)

  if (nrow(validation_results) > 0) {
    fwrite(validation_results, file.path(OUTPUT_BASE, "summary", "positive_control_validation.csv"))
    message("\nSaved: summary/positive_control_validation.csv")
  }
}

message("\n", strrep("=", 70))
message("Analysis Complete!")
message("Output directory: ", OUTPUT_BASE)
message("Model: Poisson GLM | K neighbors: ", K_NEIGHBORS)
message(strrep("=", 70))
