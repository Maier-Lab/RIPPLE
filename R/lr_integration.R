#' @title Ligand-Receptor Integration Functions
#'
#' @description Functions for connecting RIPPLE gradient results to specific
#'   ligand-receptor communication mechanisms using NicheNet databases.
#'   Implements three complementary approaches: direct L-R mapping, NicheNet
#'   ligand activity prediction, and downstream target enrichment.
#'
#' @name lr_integration
NULL


# ============================================================================
# Internal helpers
# ============================================================================

#' Load NicheNet L-R network for a given organism
#'
#' Loads the ligand-receptor network and ligand-target matrix from the
#' nichenetr package, with gene name standardization.
#'
#' @param organism Character: "mouse" or "human"
#'
#' @return List with components: lr_network (data.table), ligand_target_matrix (matrix)
#'
#' @noRd
.load_nichenet_db <- function(organism = "mouse") {
  if (!requireNamespace("nichenetr", quietly = TRUE)) {
    stop("Package 'nichenetr' is required for L-R integration.\n",
         "Install with: devtools::install_github('saeyslab/nichenetr')",
         call. = FALSE)
  }

  message("  Loading NicheNet L-R network...")

  lr_network_all <- nichenetr::lr_network
  if (is.null(lr_network_all) || nrow(lr_network_all) == 0) {
    stop("Failed to load NicheNet L-R network. ",
         "Ensure nichenetr is properly installed with its data files.",
         call. = FALSE)
  }

  # Standardize gene names
  lr_network_all <- as.data.frame(lr_network_all)
  lr_network_all$ligand <- nichenetr::convert_alias_to_symbols(
    lr_network_all$ligand, organism = organism
  )
  lr_network_all$receptor <- nichenetr::convert_alias_to_symbols(
    lr_network_all$receptor, organism = organism
  )
  lr_network_all$ligand <- make.names(lr_network_all$ligand)
  lr_network_all$receptor <- make.names(lr_network_all$receptor)

  lr_network <- data.table::as.data.table(
    unique(lr_network_all[, c("ligand", "receptor")])
  )

  message("  Loading NicheNet ligand-target matrix...")

  ligand_target_matrix <- nichenetr::ligand_target_matrix
  if (is.null(ligand_target_matrix)) {
    stop("Failed to load NicheNet ligand-target matrix.",
         call. = FALSE)
  }

  colnames(ligand_target_matrix) <- make.names(
    nichenetr::convert_alias_to_symbols(
      colnames(ligand_target_matrix), organism = organism
    )
  )
  rownames(ligand_target_matrix) <- make.names(
    nichenetr::convert_alias_to_symbols(
      rownames(ligand_target_matrix), organism = organism
    )
  )

  # Filter L-R network to ligands present in target matrix
  lr_network <- lr_network[ligand %in% colnames(ligand_target_matrix)]
  ligand_target_matrix <- ligand_target_matrix[, unique(lr_network$ligand),
                                                 drop = FALSE]

  message(sprintf("  L-R network: %d pairs | Ligand-target matrix: %d targets x %d ligands",
                  nrow(lr_network), nrow(ligand_target_matrix),
                  ncol(ligand_target_matrix)))

  list(lr_network = lr_network, ligand_target_matrix = ligand_target_matrix)
}


#' Compute expression profile for a set of cells
#'
#' @param expr_matrix Sparse expression matrix (genes x cells)
#' @param barcodes Character vector of cell barcodes
#'
#' @return data.table with columns: gene, pct_expr, mean_expr
#' @noRd
.compute_expr_profile <- function(expr_matrix, barcodes) {
  sub_expr <- expr_matrix[, barcodes, drop = FALSE]
  data.table::data.table(
    gene = rownames(expr_matrix),
    pct_expr = as.numeric(Matrix::rowMeans(sub_expr > 0) * 100),
    mean_expr = as.numeric(Matrix::rowMeans(sub_expr))
  )
}


#' Run Fisher's exact test for downstream target enrichment
#'
#' @param ligand_target_matrix Ligand-target weight matrix
#' @param ligand Character: ligand gene name
#' @param gradient_genes Character vector of significant gradient genes
#' @param universe_genes Character vector of all tested genes
#' @param n_top_targets Integer: number of top predicted targets to use
#'
#' @return data.table with one row of enrichment results, or NULL
#' @noRd
.test_downstream_enrichment <- function(ligand_target_matrix, ligand,
                                         gradient_genes, universe_genes,
                                         n_top_targets = 200) {
  if (!ligand %in% colnames(ligand_target_matrix)) return(NULL)

  target_weights <- ligand_target_matrix[, ligand]
  target_weights <- target_weights[names(target_weights) %in% universe_genes]
  target_weights <- sort(target_weights, decreasing = TRUE)
  predicted_targets <- names(utils::head(target_weights, n_top_targets))

  n_universe <- length(universe_genes)
  n_targets <- length(predicted_targets)
  n_gradient <- length(intersect(gradient_genes, universe_genes))
  overlap <- intersect(predicted_targets, gradient_genes)
  n_overlap <- length(overlap)

  # 2x2 contingency table for Fisher's exact test
  a <- n_overlap
  b <- n_targets - n_overlap
  c_val <- n_gradient - n_overlap
  d <- n_universe - n_targets - n_gradient + n_overlap

  fisher_res <- tryCatch({
    stats::fisher.test(matrix(c(a, c_val, b, d), nrow = 2),
                        alternative = "greater")
  }, error = function(e) {
    list(p.value = 1, estimate = 1)
  })

  jaccard <- if ((n_targets + n_gradient - n_overlap) > 0) {
    n_overlap / (n_targets + n_gradient - n_overlap)
  } else {
    0
  }

  data.table::data.table(
    n_predicted_targets = n_targets,
    n_gradient_genes = n_gradient,
    n_overlap = n_overlap,
    overlap_genes = paste(overlap, collapse = ";"),
    pvalue_fisher = fisher_res$p.value,
    odds_ratio = as.numeric(fisher_res$estimate),
    jaccard_index = jaccard
  )
}


# ============================================================================
# Main L-R Integration Function
# ============================================================================

#' Run ligand-receptor integration for RIPPLE gradient results
#'
#' Matches significant gradient genes to ligand-receptor pairs using NicheNet,
#' computing a combined prioritization score that integrates direct L-R mapping,
#' NicheNet ligand activity prediction, and downstream target enrichment.
#'
#' @param results_dir Character. Path to RIPPLE results directory (Stage 1 or
#'   merged). Must contain \code{summary/all_genes_results.csv} and optionally
#'   \code{per_celltype/<ct>/coef_per_sample.csv} for per-sample scoring.
#' @param input A Seurat, SingleCellExperiment, or SpatialExperiment object,
#'   or a path to an \code{.rds} file containing one of these. Used to extract
#'   the normalized expression matrix for ligand/receptor filtering.
#' @param query_celltype Character. Query cell type label as it appears in
#'   the cell type metadata column.
#' @param celltype_column Character. Name of the metadata column containing
#'   cell type annotations.
#' @param sample_column Character. Name of the metadata column containing
#'   sample identifiers. Default: "sample_id".
#' @param condition_column Character or NULL. Condition column name for
#'   subsetting samples. Default: NULL (use all samples).
#' @param condition_value Character or NULL. Condition value to filter to.
#'   Default: NULL.
#' @param organism Character. "mouse" or "human". Default: "mouse".
#' @param expr_threshold_pct Numeric. Minimum expression percentage in query
#'   cells for a ligand/receptor to be considered expressed. Default: 5.
#' @param fdr_threshold Numeric. FDR cutoff for gradient genes. Default: 0.05.
#' @param coef_col Character. Column name for gradient coefficients.
#'   Default: "median_coef".
#' @param fdr_col Character. Column name for FDR values.
#'   Default: "fisher_fdr".
#' @param output_dir Character or NULL. Output directory. If NULL, defaults to
#'   \code{<parent of results_dir>/gradient_lr_integration}.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A \code{data.table} with combined prioritization results including
#'   L-R pairs, scores, expression data, and per-sample reproducibility
#'   information. Columns include:
#' \describe{
#'   \item{ligand, receptor}{Gene names of the L-R pair.}
#'   \item{cell_type}{Target cell type.}
#'   \item{direct_score}{Score from direct L-R mapping.}
#'   \item{nichenet_activity}{NicheNet ligand activity (AUPR).}
#'   \item{enrichment_pval}{Fisher's exact p-value for downstream enrichment.}
#'   \item{combined_score}{Weighted average of all scoring components.}
#'   \item{n_samples_supporting}{Number of samples supporting this pair.}
#'   \item{reproducibility}{Fraction of samples supporting.}
#' }
#'
#' @details
#' The function implements three complementary approaches:
#' \describe{
#'   \item{Part 1: Direct L-R mapping}{Gradient genes that are known receptors
#'     are matched to query-expressed ligands. Per-sample scoring verifies
#'     reproducibility.}
#'   \item{Part 2: NicheNet ligand activity}{Predicts which query-expressed
#'     ligands best explain the gradient gene signature using NicheNet's
#'     ligand-target model.}
#'   \item{Part 3: Downstream enrichment}{Fisher's exact test validates that
#'     predicted downstream targets of top ligands overlap with gradient genes.}
#' }
#'
#' The combined score is a weighted average:
#' \code{(4 * direct + 3 * activity + 2 * enrichment + 1 * expression) / 10},
#' where each component is rescaled to [0, 1].
#'
#' @examples
#' \dontrun{
#' lr_results <- run_ripple_lr(
#'   results_dir     = "output/distance_correlation/",
#'   input           = my_spe,   # or a Seurat/SCE object, or a path to an .rds
#'   query_celltype  = "Neutrophil",
#'   celltype_column = "cell_type",
#'   organism        = "mouse"
#' )
#' # Filter to clean, high-confidence pairs
#' lr_results[combined_score > 0.3]
#' }
#'
#' @importFrom data.table fread fwrite data.table rbindlist setorder copy as.data.table
#' @importFrom Matrix rowMeans
#' @importFrom scales rescale
#' @export
run_ripple_lr <- function(results_dir,
                           input,
                           query_celltype,
                           celltype_column,
                           sample_column = "sample_id",
                           condition_column = NULL,
                           condition_value = NULL,
                           organism = "mouse",
                           expr_threshold_pct = 5,
                           fdr_threshold = 0.05,
                           coef_col = "median_coef",
                           fdr_col = "fisher_fdr",
                           output_dir = NULL,
                           verbose = TRUE) {

  .msg <- function(...) if (isTRUE(verbose)) message(...)
  .ensure_dir <- function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    path
  }

  set.seed(42)

  # --- Setup output ---
  if (is.null(output_dir)) {
    output_dir <- file.path(dirname(results_dir), "gradient_lr_integration")
  }
  .ensure_dir(output_dir)

  .msg(strrep("=", 70))
  .msg("RIPPLE Ligand-Receptor Integration")
  .msg(strrep("=", 70))
  .msg("Results directory: ", results_dir)
  .msg("Query cell type: ", query_celltype)
  .msg("Output: ", output_dir)

  # --- Load gradient results ---
  .msg("\nLoading gradient results...")
  gradient_file <- file.path(results_dir, "summary", "all_genes_results.csv")
  if (!file.exists(gradient_file)) {
    stop("Gradient results not found: ", gradient_file, call. = FALSE)
  }
  all_gradient <- data.table::fread(gradient_file)
  .msg(sprintf("  Loaded %d gene-celltype rows", nrow(all_gradient)))

  # Validate columns
  if (!coef_col %in% names(all_gradient)) {
    if ("combined_coef" %in% names(all_gradient)) {
      .msg("  Column '", coef_col, "' not found, using 'combined_coef'")
      coef_col <- "combined_coef"
    } else {
      stop("Coefficient column '", coef_col, "' not found.", call. = FALSE)
    }
  }
  if (!fdr_col %in% names(all_gradient)) {
    if ("fdr" %in% names(all_gradient)) {
      .msg("  Column '", fdr_col, "' not found, using 'fdr'")
      fdr_col <- "fdr"
    } else {
      stop("FDR column '", fdr_col, "' not found.", call. = FALSE)
    }
  }

  # --- Load NicheNet databases ---
  .msg("\nLoading NicheNet databases...")
  nn_db <- .load_nichenet_db(organism)
  lr_network <- nn_db$lr_network
  ligand_target_matrix <- nn_db$ligand_target_matrix

  # --- Load data (Seurat, SCE, or SpatialExperiment) ---
  .msg("\nLoading input data...")
  if (missing(input)) {
    stop("'input' is required (Seurat/SCE/SPE object or path to an .rds file)",
         call. = FALSE)
  }
  data <- .resolve_input(input, require_expr = TRUE, verbose = verbose)
  meta <- data$meta
  expr_matrix <- data$expr

  if (!celltype_column %in% names(meta)) {
    stop("Cell type column '", celltype_column, "' not found in metadata. ",
         "Available: ", paste(head(names(meta), 20), collapse = ", "),
         call. = FALSE)
  }
  if (!sample_column %in% names(meta)) {
    stop("Sample column '", sample_column, "' not found in metadata. ",
         "Available: ", paste(head(names(meta), 20), collapse = ", "),
         call. = FALSE)
  }
  if (!query_celltype %in% unique(meta[[celltype_column]])) {
    stop("query_celltype '", query_celltype, "' not found in column '",
         celltype_column, "'.", call. = FALSE)
  }

  # Get query cells
  if (!is.null(condition_column) && !is.null(condition_value)) {
    query_barcodes <- meta[get(celltype_column) == query_celltype &
                            get(condition_column) == condition_value]$barcode
    analysis_samples <- unique(meta[get(condition_column) == condition_value][[sample_column]])
    .msg(sprintf("  Query cells (%s, %s): %d", query_celltype,
                 condition_value, length(query_barcodes)))
  } else {
    query_barcodes <- meta[get(celltype_column) == query_celltype]$barcode
    analysis_samples <- unique(meta[[sample_column]])
    .msg(sprintf("  Query cells (%s, all): %d", query_celltype,
                 length(query_barcodes)))
  }

  if (length(query_barcodes) == 0) {
    stop("No query cells found. Check celltype_column and query_celltype.",
         call. = FALSE)
  }

  # Sanitize gene names for downstream matching
  rownames(expr_matrix) <- make.names(rownames(expr_matrix))

  # Compute query expression profile (pooled)
  query_profile <- .compute_expr_profile(expr_matrix, query_barcodes)

  # Per-sample query expression profiles
  query_per_sample <- data.table::rbindlist(lapply(analysis_samples, function(sid) {
    barcodes <- meta[get(celltype_column) == query_celltype &
                      get(sample_column) == sid]$barcode
    if (length(barcodes) == 0) return(NULL)
    prof <- .compute_expr_profile(expr_matrix, barcodes)
    prof[, sample_id := sid]
    prof[, n_query_cells := length(barcodes)]
    prof
  }))

  # Identify query-expressed ligands and receptors
  available_genes <- rownames(expr_matrix)
  available_ligands <- intersect(lr_network$ligand, available_genes)
  available_receptors <- intersect(lr_network$receptor, available_genes)

  query_expressed_ligands <- query_profile[
    gene %in% available_ligands & pct_expr >= expr_threshold_pct
  ]$gene
  query_expressed_receptors <- query_profile[
    gene %in% available_receptors & pct_expr >= expr_threshold_pct
  ]$gene

  .msg(sprintf("  Query-expressed ligands (>=%d%%): %d",
               expr_threshold_pct, length(query_expressed_ligands)))
  .msg(sprintf("  Query-expressed receptors (>=%d%%): %d",
               expr_threshold_pct, length(query_expressed_receptors)))

  # Free Seurat object memory
  rm(obj)
  gc(verbose = FALSE)

  # --- Process each target cell type ---
  target_celltypes <- unique(all_gradient$cell_type)
  .msg(sprintf("\nProcessing %d target cell types...", length(target_celltypes)))

  all_combined <- data.table::data.table()

  for (ct_name in target_celltypes) {
    .msg(sprintf("\n--- %s ---", ct_name))

    ct_output_dir <- .ensure_dir(file.path(output_dir, "per_celltype", ct_name))

    # Get gradient data for this cell type
    ct_gradient <- all_gradient[cell_type == ct_name]
    if (nrow(ct_gradient) == 0) {
      .msg("  No gradient data. Skipping.")
      next
    }
    ct_gradient[, gene_safe := make.names(gene)]

    sig_genes <- ct_gradient[get(fdr_col) < fdr_threshold]
    induced_genes <- sig_genes[get(coef_col) < 0]$gene_safe
    repressed_genes <- sig_genes[get(coef_col) > 0]$gene_safe
    background_genes <- ct_gradient$gene_safe

    .msg(sprintf("  Significant: %d (%d induced, %d repressed)",
                 nrow(sig_genes), length(induced_genes), length(repressed_genes)))

    # Load per-sample coefficients if available
    coef_per_sample_file <- file.path(
      results_dir, "per_celltype", ct_name, "coef_per_sample.csv"
    )
    coef_per_sample <- if (file.exists(coef_per_sample_file)) {
      dt <- data.table::fread(coef_per_sample_file)
      dt[, gene_safe := make.names(gene)]
      dt
    } else {
      NULL
    }

    # =====================================================================
    # Part 1: Direct L-R Mapping (Query ligand -> Target receptor)
    # =====================================================================
    .msg("  Part 1: Direct L-R mapping...")

    sig_receptors <- sig_genes[gene_safe %in% available_receptors]
    direct_a_results <- data.table::data.table()

    if (nrow(sig_receptors) > 0) {
      for (i in seq_len(nrow(sig_receptors))) {
        receptor_gene <- sig_receptors$gene_safe[i]
        receptor_info <- sig_receptors[i]

        cognate_ligands <- lr_network[
          receptor == receptor_gene & ligand %in% query_expressed_ligands
        ]$ligand

        for (lig in cognate_ligands) {
          lig_info <- query_profile[gene == lig]

          # Per-sample scoring
          per_sample_scores <- numeric(0)
          n_supporting <- 0L

          if (!is.null(coef_per_sample)) {
            for (sid in analysis_samples) {
              lig_sample <- query_per_sample[gene == lig & sample_id == sid]
              rec_sample <- coef_per_sample[gene_safe == receptor_gene &
                                             sample_id == sid]

              if (nrow(lig_sample) > 0 && nrow(rec_sample) > 0) {
                lig_pct <- lig_sample$pct_expr[1]
                rec_coef <- rec_sample$coef[1]
                expected_neg <- (receptor_info[[coef_col]] < 0)
                agrees <- if (expected_neg) rec_coef < 0 else rec_coef > 0

                if (lig_pct >= expr_threshold_pct && !is.na(rec_coef) && agrees) {
                  per_sample_scores <- c(per_sample_scores,
                                          abs(rec_coef) * (lig_pct / 100))
                  n_supporting <- n_supporting + 1L
                }
              }
            }
          }

          n_total_samples <- length(analysis_samples)
          reproducibility <- n_supporting / n_total_samples

          median_score <- if (length(per_sample_scores) >= 2) {
            stats::median(per_sample_scores)
          } else if (length(per_sample_scores) == 1) {
            per_sample_scores[1]
          } else {
            abs(receptor_info[[coef_col]]) * (lig_info$pct_expr / 100)
          }

          pct_target_col <- intersect(c("pct_expr", "pct"), names(receptor_info))
          rpt <- if (length(pct_target_col) > 0) {
            receptor_info[[pct_target_col[1]]]
          } else {
            NA_real_
          }

          direct_a_results <- rbind(direct_a_results, data.table::data.table(
            ligand = lig,
            receptor = receptor_gene,
            cell_type = ct_name,
            ligand_pct_query = lig_info$pct_expr,
            ligand_mean_query = lig_info$mean_expr,
            receptor_gradient_coef = receptor_info[[coef_col]],
            receptor_fdr = receptor_info[[fdr_col]],
            receptor_pct_target = rpt,
            n_samples_supporting = n_supporting,
            n_samples_total = n_total_samples,
            reproducibility = reproducibility,
            median_per_sample_score = median_score
          ), fill = TRUE)
        }
      }
    }

    if (nrow(direct_a_results) > 0) {
      direct_a_results[, direct_score := median_per_sample_score * reproducibility]
      data.table::setorder(direct_a_results, -direct_score)
      .msg(sprintf("  Direction A: %d L-R pairs", nrow(direct_a_results)))
    } else {
      .msg("  Direction A: 0 L-R pairs")
    }

    data.table::fwrite(direct_a_results,
                        file.path(ct_output_dir, "direct_lr_pairs.csv"))

    # =====================================================================
    # Part 2: NicheNet Ligand Activity
    # =====================================================================
    .msg("  Part 2: NicheNet ligand activity...")
    activity_results <- data.table::data.table()

    response_induced <- intersect(induced_genes, rownames(ligand_target_matrix))
    response_all_sig <- intersect(sig_genes$gene_safe, rownames(ligand_target_matrix))
    bg_in_matrix <- intersect(background_genes, rownames(ligand_target_matrix))
    potential_ligands <- intersect(query_expressed_ligands,
                                    colnames(ligand_target_matrix))

    response_in_matrix <- if (length(response_induced) >= 10) {
      response_induced
    } else if (length(response_all_sig) >= 10) {
      .msg("  NOTE: <10 induced genes, using all significant as fallback")
      response_all_sig
    } else {
      response_induced
    }

    if (length(response_in_matrix) >= 10 && length(potential_ligands) >= 5) {
      ligand_activity <- nichenetr::predict_ligand_activities(
        geneset = response_in_matrix,
        background_expressed_genes = bg_in_matrix,
        ligand_target_matrix = ligand_target_matrix,
        potential_ligands = potential_ligands
      )
      activity_results <- data.table::as.data.table(ligand_activity)
      if ("test_ligand" %in% names(activity_results)) {
        data.table::setnames(activity_results, "test_ligand", "ligand")
      }
      if ("aupr_corrected" %in% names(activity_results)) {
        data.table::setnames(activity_results, "aupr_corrected", "activity")
      }
      activity_results[, cell_type := ct_name]
      data.table::setorder(activity_results, -activity)
      .msg(sprintf("  Ligand activities: %d ligands scored", nrow(activity_results)))
    } else {
      .msg("  Insufficient genes for ligand activity. Skipping.")
    }

    data.table::fwrite(activity_results,
                        file.path(ct_output_dir, "nichenet_ligand_activity.csv"))

    # =====================================================================
    # Part 3: Downstream Target Enrichment
    # =====================================================================
    .msg("  Part 3: Downstream target enrichment...")

    top_ligands_direct <- if (nrow(direct_a_results) > 0) {
      utils::head(unique(direct_a_results$ligand), 30)
    } else {
      character(0)
    }
    top_ligands_activity <- if (nrow(activity_results) > 0) {
      utils::head(activity_results$ligand, 20)
    } else {
      character(0)
    }
    test_ligands <- unique(c(top_ligands_direct, top_ligands_activity))
    test_ligands <- intersect(test_ligands, colnames(ligand_target_matrix))
    universe_genes <- intersect(background_genes, rownames(ligand_target_matrix))

    enrichment_results <- data.table::data.table()

    if (length(test_ligands) > 0 && length(sig_genes$gene_safe) > 0) {
      enrichment_results <- data.table::rbindlist(lapply(test_ligands, function(lig) {
        res <- .test_downstream_enrichment(
          ligand_target_matrix, lig,
          sig_genes$gene_safe, universe_genes
        )
        if (!is.null(res)) {
          receptors_for_lig <- lr_network[
            ligand == lig & receptor %in% available_receptors
          ]$receptor
          res[, `:=`(ligand = lig,
                     receptors = paste(receptors_for_lig, collapse = ";"),
                     cell_type = ct_name)]
        }
        res
      }), fill = TRUE)

      if (nrow(enrichment_results) > 0) {
        enrichment_results[, fdr_fisher := stats::p.adjust(pvalue_fisher,
                                                            method = "BH")]
        enrichment_results[, enrichment_score :=
          -log10(pvalue_fisher + 1e-50) * sign(log2(pmax(odds_ratio, 0.01)))]
        data.table::setorder(enrichment_results, pvalue_fisher)
        .msg(sprintf("  Enrichment: %d / %d ligands significant (FDR<0.05)",
                     sum(enrichment_results$fdr_fisher < 0.05),
                     nrow(enrichment_results)))
      }
    }

    data.table::fwrite(enrichment_results,
                        file.path(ct_output_dir, "downstream_enrichment.csv"))

    # =====================================================================
    # Combined Prioritization
    # =====================================================================
    .msg("  Combining results...")

    if (nrow(direct_a_results) > 0) {
      combined <- data.table::copy(direct_a_results)

      # Merge NicheNet activity
      if (nrow(activity_results) > 0) {
        activity_cols <- intersect(
          c("ligand", "activity", "auroc"),
          names(activity_results)
        )
        act_merge <- activity_results[, ..activity_cols]
        data.table::setnames(act_merge,
          old = intersect(c("activity", "auroc"), names(act_merge)),
          new = paste0("nichenet_", intersect(c("activity", "auroc"),
                                               names(act_merge)))
        )
        combined <- merge(combined, act_merge, by = "ligand", all.x = TRUE)
      } else {
        combined[, nichenet_activity := NA_real_]
      }

      # Merge enrichment
      if (nrow(enrichment_results) > 0) {
        enrich_merge <- enrichment_results[, .(
          ligand, enrichment_pval = pvalue_fisher,
          enrichment_fdr = fdr_fisher,
          enrichment_odds_ratio = odds_ratio,
          n_downstream_overlap = n_overlap,
          enrichment_score
        )]
        combined <- merge(combined, enrich_merge, by = "ligand", all.x = TRUE)
      } else {
        combined[, `:=`(enrichment_pval = NA_real_, enrichment_fdr = NA_real_,
                        enrichment_odds_ratio = NA_real_,
                        n_downstream_overlap = NA_integer_,
                        enrichment_score = NA_real_)]
      }

      # Combined score (weighted average of rescaled components)
      w_direct <- 4
      w_activity <- 3
      w_enrichment <- 2
      w_expression <- 1

      combined[, scaled_direct := scales::rescale(
        direct_score, to = c(0, 1), na.rm = TRUE
      )]
      combined[, scaled_activity := data.table::fifelse(
        is.na(nichenet_activity), 0,
        scales::rescale(nichenet_activity, to = c(0, 1), na.rm = TRUE)
      )]
      combined[, scaled_enrichment := data.table::fifelse(
        is.na(enrichment_score), 0,
        scales::rescale(pmax(enrichment_score, 0), to = c(0, 1), na.rm = TRUE)
      )]
      combined[, scaled_expression := scales::rescale(
        ligand_pct_query / 100, to = c(0, 1), na.rm = TRUE
      )]

      combined[, combined_score := (
        w_direct * scaled_direct +
          w_activity * scaled_activity +
          w_enrichment * scaled_enrichment +
          w_expression * scaled_expression
      ) / (w_direct + w_activity + w_enrichment + w_expression)]

      data.table::setorder(combined, -combined_score)
      .msg(sprintf("  Combined: %d L-R pairs scored", nrow(combined)))
    } else {
      combined <- data.table::data.table()
      .msg("  No Direction A pairs to combine.")
    }

    data.table::fwrite(combined,
                        file.path(ct_output_dir, "combined_prioritization.csv"))

    all_combined <- rbind(all_combined, combined, fill = TRUE)
  }

  # --- Save overall summary ---
  if (nrow(all_combined) > 0) {
    summary_dir <- .ensure_dir(file.path(output_dir, "summary"))
    data.table::fwrite(all_combined,
                        file.path(summary_dir, "all_lr_pairs_combined.csv"))
    .msg(sprintf("\n  Total L-R pairs: %d across %d cell types",
                 nrow(all_combined),
                 data.table::uniqueN(all_combined$cell_type)))
  }

  .msg("\n=== L-R Integration Complete ===")
  .msg("  Output: ", output_dir)

  all_combined
}


# ============================================================================
# Artifact Classification
# ============================================================================

#' Classify ligand-receptor artifacts
#'
#' Applies a tiered artifact classification to L-R pairs based on query
#' signature leakage, contamination patterns, and expression thresholds.
#' Rules are applied in order (first match wins), producing a single flag
#' per L-R pair.
#'
#' @param lr_results \code{data.table} from \code{run_ripple_lr} with at
#'   minimum the columns: receptor, cell_type. Also uses
#'   \code{receptor_pct_target} if available.
#' @param query_signature Character vector of query cell marker genes.
#'   Receptors matching these on non-myeloid cells are flagged as "artifact".
#'   Default: empty vector.
#' @param contamination_genes Character vector of genes significant in many
#'   cell types (potential segmentation artifacts). Default: empty vector.
#' @param myeloid_celltypes Character vector of myeloid cell type names that
#'   are exempt from some artifact rules (since they share lineage with many
#'   query cell types). Default: empty vector.
#' @param low_expr_threshold Numeric. Expression fraction (0-1) below which
#'   non-myeloid receptor expression triggers a "low_confidence" flag.
#'   Default: 0.02.
#'
#' @return The input \code{data.table} with an added column
#'   \code{artifact_flag} containing one of: "artifact", "suspect",
#'   "low_confidence", or "clean".
#'
#' @details
#' Classification rules (applied in order, first match wins):
#' \describe{
#'   \item{Rule 1 - "artifact"}{Receptor is in \code{query_signature} AND
#'     cell type is NOT in \code{myeloid_celltypes}.}
#'   \item{Rule 2 - "suspect"}{Receptor is in \code{contamination_genes}
#'     AND cell type is NOT myeloid AND receptor expression < 5\%.}
#'   \item{Rule 3 - "low_confidence"}{Cell type is NOT myeloid AND
#'     receptor expression < \code{low_expr_threshold}.}
#'   \item{Default - "clean"}{All other pairs.}
#' }
#'
#' @examples
#' \dontrun{
#' lr <- run_ripple_lr(...)
#' lr <- classify_lr_artifacts(
#'   lr,
#'   query_signature = c("Csf3r", "Ly6g", "S100a8"),
#'   contamination_genes = c("C1qa", "Igkc", "Jchain"),
#'   myeloid_celltypes = c("Monocyte", "Macrophages", "cDC1")
#' )
#' # View clean pairs only
#' lr[artifact_flag == "clean"]
#' }
#'
#' @importFrom data.table copy fifelse
#' @export
classify_lr_artifacts <- function(lr_results,
                                   query_signature = character(0),
                                   contamination_genes = character(0),
                                   myeloid_celltypes = character(0),
                                   low_expr_threshold = 0.02) {

  if (!inherits(lr_results, "data.table")) {
    lr_results <- data.table::as.data.table(lr_results)
  }

  lr_results <- data.table::copy(lr_results)

  if (nrow(lr_results) == 0) {
    lr_results[, artifact_flag := character(0)]
    return(lr_results)
  }

  if (!"receptor" %in% names(lr_results)) {
    stop("Column 'receptor' not found in lr_results.", call. = FALSE)
  }
  if (!"cell_type" %in% names(lr_results)) {
    stop("Column 'cell_type' not found in lr_results.", call. = FALSE)
  }

  # Standardize gene names for matching
  query_signature <- make.names(query_signature)
  contamination_genes <- make.names(contamination_genes)

  # Initialize all as clean
  lr_results[, artifact_flag := "clean"]

  # Rule 1: Receptor is query signature gene on non-myeloid cell type
  if (length(query_signature) > 0) {
    lr_results[receptor %in% query_signature &
               !(cell_type %in% myeloid_celltypes),
               artifact_flag := "artifact"]
  }

  # Rule 2: Receptor in contamination list on non-myeloid at low expression
  has_pct <- "receptor_pct_target" %in% names(lr_results)
  if (length(contamination_genes) > 0 && has_pct) {
    lr_results[artifact_flag == "clean" &
               receptor %in% contamination_genes &
               !(cell_type %in% myeloid_celltypes) &
               receptor_pct_target < 0.05,
               artifact_flag := "suspect"]
  }

  # Rule 3: Very low receptor expression on non-myeloid cells
  if (has_pct) {
    lr_results[artifact_flag == "clean" &
               !(cell_type %in% myeloid_celltypes) &
               receptor_pct_target < low_expr_threshold,
               artifact_flag := "low_confidence"]
  }

  # Summary
  flag_counts <- lr_results[, .N, by = artifact_flag]
  for (i in seq_len(nrow(flag_counts))) {
    message(sprintf("  %s: %d pairs", flag_counts$artifact_flag[i],
                    flag_counts$N[i]))
  }

  lr_results[]
}
