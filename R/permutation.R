#' @title Permutation Testing Functions
#'
#' @description Functions for label-permutation testing to validate query cell
#'   specificity of distance-expression gradients, and for merging GPU
#'   permutation results.
#'
#' @name permutation
NULL

#' Single-gene permutation test
#'
#' Validates a distance-expression gradient by shuffling query cell labels
#' within each sample and recalculating the distance-expression coefficient.
#' Uses stratified sampling to preserve per-sample query cell counts.
#'
#' @param counts Integer vector of raw transcript counts for target cells.
#' @param coords_target Numeric matrix (n_target x 2). Coordinates of target
#'   cells.
#' @param coords_all Numeric matrix (n_all x 2). Coordinates of ALL cells
#'   (for drawing pseudo-query cells from).
#' @param sample_ids Character vector (length n_target). Sample IDs for each
#'   target cell.
#' @param n_perms Integer. Number of permutations (default: 500).
#' @param observed_coef Numeric. The observed combined coefficient from the
#'   real analysis.
#' @param sample_ids_all Character vector (length n_all). Sample IDs for all
#'   cells.
#' @param query_per_sample Named integer vector. Number of query cells per
#'   sample (names = sample IDs).
#' @param k_neighbors Integer. Number of nearest neighbors for distance
#'   calculation (default: 1).
#' @param total_counts Numeric vector (length n_target). Total counts per
#'   target cell (for Poisson offset).
#' @param max_distance Numeric. Maximum distance to consider in um
#'   (default: 200).
#' @param min_cells_per_sample Integer. Minimum target cells per sample
#'   (default: 30).
#' @param min_expr_cells Integer. Minimum expressing cells for GLM fit
#'   (default: 5).
#'
#' @return A list with:
#' \describe{
#'   \item{\code{null_coefs}}{Numeric vector of null distribution coefficients.}
#'   \item{\code{perm_pval}}{Numeric. Two-sided empirical p-value.}
#' }
#'
#' @details For each permutation:
#' \enumerate{
#'   \item Randomly sample pseudo-query cells WITHIN each sample (stratified),
#'     preserving the original query cell count per sample.
#'   \item Compute distances from all target cells to the nearest pseudo-query cell.
#'   \item Fit per-sample Poisson GLMs and take the median of per-sample
#'     coefficients -- the same statistic as the observed \code{median_coef},
#'     so the empirical p-value compares like with like.
#' }
#'
#' The empirical p-value is calculated as:
#'   \code{(sum(|null| >= |observed|) + 1) / (n_valid_perms + 1)}
#'
#' @examples
#' \dontrun{
#' result <- run_permutation_test(
#'   counts = rpois(100, 5),
#'   coords_target = matrix(runif(200), ncol = 2),
#'   coords_all = matrix(runif(2000), ncol = 2),
#'   sample_ids = rep(c("s1", "s2"), each = 50),
#'   n_perms = 100,
#'   observed_coef = -0.005,
#'   sample_ids_all = rep(c("s1", "s2"), each = 500),
#'   query_per_sample = c(s1 = 50, s2 = 60),
#'   k_neighbors = 1,
#'   total_counts = rpois(100, 5000)
#' )
#' }
#'
#' @importFrom RANN nn2
#' @importFrom stats glm poisson median
#' @export
run_permutation_test <- function(counts, coords_target, coords_all,
                                 sample_ids, n_perms = 500,
                                 observed_coef, sample_ids_all,
                                 query_per_sample, k_neighbors = 1,
                                 total_counts,
                                 max_distance = 200,
                                 min_cells_per_sample = 30,
                                 min_expr_cells = 5) {
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
      perm_distances <- pmin(as.vector(nn_result$nn.dists), max_distance)
    } else {
      perm_distances <- pmin(rowMeans(nn_result$nn.dists), max_distance)
    }

    # Calculate Poisson coefficients per sample using inverse-variance weighting
    coefs <- numeric(length(unique_samples))
    ses <- numeric(length(unique_samples))

    for (j in seq_along(unique_samples)) {
      samp <- unique_samples[j]
      idx <- which(sample_ids == samp)

      if (length(idx) >= min_cells_per_sample) {
        samp_counts <- counts[idx]
        samp_dist <- perm_distances[idx]
        samp_log_total <- log(total_counts[idx])

        if (length(samp_counts) >= min_cells_per_sample &&
          sum(samp_counts > 0) >= min_expr_cells) {
          fit <- tryCatch(
            {
              suppressWarnings(stats::glm(samp_counts ~ samp_dist + offset(samp_log_total),
                family = stats::poisson()
              ))
            },
            error = function(e) NULL
          )

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

    # Null statistic must match the observed statistic (median_coef, the
    # equal-weight median of per-sample coefficients from compute_fisher_pval).
    # ses > 0 identifies samples whose GLM actually converged.
    valid <- !is.na(coefs) & !is.na(ses) & ses > 0
    if (sum(valid) >= 2) {
      null_coefs[i] <- stats::median(coefs[valid])
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


#' Run permutation tests for multiple genes
#'
#' Batch version of \code{\link{run_permutation_test}} that loops over a vector
#' of genes and returns a \code{data.table} of gene-level permutation p-values.
#' Used internally by \code{\link{run_ripple}} to validate the top significant
#' distance-expression gradients.
#'
#' @param genes Character vector. Gene names to test.
#' @param count_matrix Sparse or dense matrix. Raw count matrix (genes x cells).
#' @param target_barcodes Character vector. Barcodes of target cells.
#' @param coords_target Numeric matrix (n_target x 2). Coordinates of target
#'   cells.
#' @param coords_all Numeric matrix (n_all x 2). Coordinates of ALL cells.
#' @param sample_ids_target Character vector. Sample IDs for target cells.
#' @param sample_ids_all Character vector. Sample IDs for all cells.
#' @param query_per_sample Named integer vector. Number of query cells per
#'   sample.
#' @param observed_coefs Named numeric vector. Observed combined coefficients
#'   per gene (names = gene names).
#' @param n_perms Integer. Number of permutations per gene.
#' @param k_neighbors Integer. Number of nearest neighbors for distance
#'   calculation.
#' @param max_distance_um Numeric. Maximum distance in micrometers.
#' @param min_cells_per_sample Integer. Minimum target cells per sample for
#'   GLM fitting.
#' @param min_expr_cells Integer. Minimum expressing cells for GLM fitting.
#' @param total_counts_target Numeric vector. Total UMI counts per target cell
#'   (for Poisson offset).
#'
#' @return A \code{data.table} with columns \code{gene} and \code{perm_pval}.
#'
#' @examples
#' \dontrun{
#' perm_dt <- run_permutation_tests(
#'   genes = c("Cxcl12", "Ccl21a"),
#'   count_matrix = counts,
#'   target_barcodes = barcodes,
#'   coords_target = target_xy,
#'   coords_all = all_xy,
#'   sample_ids_target = target_samples,
#'   sample_ids_all = all_samples,
#'   query_per_sample = c(s1 = 50, s2 = 60),
#'   observed_coefs = c(Cxcl12 = -0.005, Ccl21a = -0.003),
#'   n_perms = 500,
#'   k_neighbors = 1,
#'   max_distance_um = 200,
#'   min_cells_per_sample = 30,
#'   min_expr_cells = 5,
#'   total_counts_target = total_counts
#' )
#' }
#'
#' @importFrom RANN nn2
#' @importFrom stats glm poisson median
#' @importFrom data.table data.table rbindlist
#' @export
run_permutation_tests <- function(genes, count_matrix, target_barcodes,
                                  coords_target, coords_all,
                                  sample_ids_target, sample_ids_all,
                                  query_per_sample, observed_coefs,
                                  n_perms, k_neighbors, max_distance_um,
                                  min_cells_per_sample, min_expr_cells,
                                  total_counts_target) {
  unique_samples <- names(query_per_sample)

  results <- lapply(genes, function(g) {
    count_vec <- as.numeric(count_matrix[g, target_barcodes])
    obs_coef <- observed_coefs[g]

    null_coefs <- numeric(n_perms)
    for (i in seq_len(n_perms)) {
      # Stratified sampling within each sample
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

      if (is.null(pseudo_query_coords) || nrow(pseudo_query_coords) < 5) {
        null_coefs[i] <- NA
        next
      }

      eff_k <- min(k_neighbors, nrow(pseudo_query_coords))
      nn_res <- RANN::nn2(pseudo_query_coords, coords_target, k = eff_k)
      if (eff_k == 1) {
        perm_distances <- pmin(as.vector(nn_res$nn.dists), max_distance_um)
      } else {
        perm_distances <- pmin(rowMeans(nn_res$nn.dists), max_distance_um)
      }

      coefs <- numeric(length(unique_samples))
      ses <- numeric(length(unique_samples))
      for (j in seq_along(unique_samples)) {
        samp <- unique_samples[j]
        idx <- which(sample_ids_target == samp)
        if (length(idx) >= min_cells_per_sample) {
          samp_counts <- count_vec[idx]
          samp_dist <- perm_distances[idx]
          samp_log_total <- log(total_counts_target[idx])
          if (sum(samp_counts > 0) >= min_expr_cells) {
            fit <- tryCatch(
              {
                suppressWarnings(stats::glm(
                  samp_counts ~ samp_dist + offset(samp_log_total),
                  family = stats::poisson
                ))
              },
              error = function(e) NULL
            )
            if (!is.null(fit) && fit$converged) {
              cs <- summary(fit)$coefficients
              if (nrow(cs) >= 2) {
                coefs[j] <- cs[2, "Estimate"]
                ses[j] <- cs[2, "Std. Error"]
              }
            }
          }
        }
      }

      # Null statistic must match the observed median_coef (see
      # run_permutation_test). ses > 0 marks samples whose GLM converged.
      valid <- !is.na(coefs) & !is.na(ses) & ses > 0
      if (sum(valid) >= 2) {
        null_coefs[i] <- stats::median(coefs[valid])
      } else {
        null_coefs[i] <- NA
      }
    }

    null_coefs <- null_coefs[!is.na(null_coefs)]
    if (length(null_coefs) < 10) {
      perm_pval <- NA_real_
    } else {
      perm_pval <- (sum(abs(null_coefs) >= abs(obs_coef)) + 1) /
        (length(null_coefs) + 1)
    }

    data.table::data.table(gene = g, perm_pval = perm_pval)
  })

  data.table::rbindlist(results)
}


#' Merge GPU permutation results into meta-analysis CSVs
#'
#' Reads GPU-produced \code{permutation_pvals.csv} files from per-celltype
#' result directories and merges them into the corresponding
#' \code{meta_analysis_results.csv} files.
#'
#' @param results_dir Character. Path to the analysis results directory
#'   (e.g., \code{"./results/spatial_analysis_Tumor/hymy_distance_correlation_v2"}).
#'   Must contain a \code{per_celltype/} subdirectory.
#'
#' @return Invisible integer. Number of cell types successfully merged.
#'
#' @details For each cell type directory under \code{results_dir/per_celltype/}:
#' \enumerate{
#'   \item Reads \code{permutation_pvals.csv} (GPU output with gene and perm_pval columns).
#'   \item Reads \code{meta_analysis_results.csv} (existing R output).
#'   \item Replaces any existing \code{perm_pval} column with GPU results.
#'   \item Overwrites \code{meta_analysis_results.csv} with updated data.
#' }
#'
#' Cell types missing either file are skipped with a message.
#'
#' @examples
#' \dontrun{
#' n_merged <- merge_permutation_results(
#'   "./results/spatial_analysis_Tumor/hymy_distance_correlation_v2"
#' )
#' }
#'
#' @importFrom data.table fread fwrite setDT
#' @export
merge_permutation_results <- function(results_dir) {
  ct_base <- file.path(results_dir, "per_celltype")

  if (!dir.exists(ct_base)) {
    stop("Per-celltype directory not found: ", ct_base)
  }

  cell_types <- basename(list.dirs(ct_base, recursive = FALSE))
  message("Found ", length(cell_types), " cell type directories")

  n_merged <- 0L

  for (ct in cell_types) {
    ct_dir <- file.path(ct_base, ct)
    meta_file <- file.path(ct_dir, "meta_analysis_results.csv")
    perm_file <- file.path(ct_dir, "permutation_pvals.csv")

    # Check meta-analysis results exist
    if (!file.exists(meta_file)) {
      message("  [SKIP] ", ct, ": meta_analysis_results.csv not found")
      next
    }

    # Check GPU permutation results exist
    if (!file.exists(perm_file)) {
      message("  [SKIP] ", ct, ": permutation_pvals.csv not found")
      next
    }

    meta <- data.table::fread(meta_file)
    perm <- data.table::fread(perm_file)

    # Validate columns
    if (!"gene" %in% names(perm) || !"perm_pval" %in% names(perm)) {
      message("  [ERROR] ", ct, ": permutation_pvals.csv missing gene/perm_pval columns")
      next
    }

    # Remove old perm_pval column and merge new one
    if ("perm_pval" %in% names(meta)) {
      meta[, perm_pval := NULL]
    }
    meta <- merge(meta, perm[, .(gene, perm_pval)], by = "gene", all.x = TRUE)
    data.table::setDT(meta)

    # Summary stats
    n_tested <- sum(!is.na(perm$perm_pval))
    n_sig <- sum(perm$perm_pval < 0.05, na.rm = TRUE)
    n_genes <- nrow(meta)

    # Save updated meta-analysis results
    data.table::fwrite(meta, meta_file)

    message(sprintf(
      "  [OK] %s: %d/%d genes with perm_pval (%d significant at p<0.05)",
      ct, n_tested, n_genes, n_sig
    ))
    n_merged <- n_merged + 1L
  }

  message(sprintf("\nMerged %d cell types", n_merged))
  invisible(n_merged)
}
