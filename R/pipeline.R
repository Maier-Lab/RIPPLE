#' =============================================================================
#' RIPPLE Pipeline Entry Points
#' =============================================================================
#'
#' Main user-facing functions that wrap the RIPPLE analysis pipeline.
#' These replace the standalone Rscript commands with callable R functions
#' that take explicit parameters.
#'
#' @name ripple-pipeline
#' @keywords internal
NULL


# ============================================================================
# Internal helpers (not exported)
# ============================================================================

#' Resolve a parameter: explicit argument > package option > default
#' @noRd
.resolve <- function(value, option_name, default = NULL) {
  if (!is.null(value)) return(value)
  opt <- getOption(paste0("ripple.", option_name))
  if (!is.null(opt) && nzchar(as.character(opt))) return(opt)
  default
}

#' Conditional message (respects verbose flag)
#' @noRd
.msg <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) message(...)
}

#' Ensure a directory exists
#' @noRd
.ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  path
}


# ============================================================================
# Priority gene list (chemokines, cytokines, receptors)
# ============================================================================

#' Default priority genes for lenient expression filtering
#' @noRd
.default_priority_genes <- function() {
  c(
    # CC Chemokines
    "Ccl1", "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl6", "Ccl7", "Ccl8", "Ccl9",
    "Ccl11", "Ccl12", "Ccl17", "Ccl19", "Ccl20", "Ccl21a", "Ccl21b", "Ccl21c",
    "Ccl22", "Ccl24", "Ccl25", "Ccl27a", "Ccl27b", "Ccl28",
    # CXC Chemokines
    "Cxcl1", "Cxcl2", "Cxcl3", "Cxcl4", "Cxcl5", "Cxcl7", "Cxcl9", "Cxcl10",
    "Cxcl11", "Cxcl12", "Cxcl13", "Cxcl14", "Cxcl15", "Cxcl16", "Cxcl17",
    # CX3C / XC Chemokines
    "Cx3cl1", "Xcl1",
    # Chemokine Receptors
    "Ccr1", "Ccr2", "Ccr3", "Ccr4", "Ccr5", "Ccr6", "Ccr7", "Ccr8", "Ccr9",
    "Ccr10", "Cxcr1", "Cxcr2", "Cxcr3", "Cxcr4", "Cxcr5", "Cxcr6",
    "Cx3cr1", "Xcr1", "Ackr1", "Ackr2", "Ackr3", "Ackr4",
    # Interleukins
    "Il1a", "Il1b", "Il1rn", "Il2", "Il3", "Il4", "Il5", "Il6", "Il7", "Il9",
    "Il10", "Il11", "Il12a", "Il12b", "Il13", "Il15", "Il16",
    "Il17a", "Il17b", "Il17c", "Il17d", "Il17f",
    "Il18", "Il21", "Il22", "Il23a", "Il25", "Il27", "Il33", "Il34",
    # Interleukin Receptors
    "Il1r1", "Il1r2", "Il1rl1", "Il1rap", "Il2ra", "Il2rb", "Il2rg",
    "Il4ra", "Il6ra", "Il6st", "Il7r", "Il10ra", "Il10rb",
    "Il17ra", "Il17rb", "Il17rc", "Il18r1", "Il18rap", "Il21r", "Il23r",
    # Interferons + Receptors
    "Ifna1", "Ifna2", "Ifnb1", "Ifng", "Ifnar1", "Ifnar2", "Ifngr1", "Ifngr2",
    # TNF Superfamily
    "Tnf", "Lta", "Ltb", "Fasl", "Cd40lg",
    "Tnfsf4", "Tnfsf8", "Tnfsf9", "Tnfsf10", "Tnfsf11", "Tnfsf13b", "Tnfsf14",
    "Tnfrsf1a", "Tnfrsf1b", "Tnfrsf4", "Fas", "Cd40",
    "Tnfrsf9", "Tnfrsf11a", "Tnfrsf11b", "Tnfrsf13b", "Tnfrsf13c", "Tnfrsf14",
    # Colony Stimulating Factors
    "Csf1", "Csf2", "Csf3", "Csf1r", "Csf2ra", "Csf2rb", "Csf3r",
    # TGF-beta Family
    "Tgfb1", "Tgfb2", "Tgfb3", "Tgfbr1", "Tgfbr2", "Tgfbr3",
    # VEGF Family
    "Vegfa", "Vegfb", "Vegfc", "Flt1", "Kdr", "Flt4", "Nrp1", "Nrp2",
    # gp130 Family
    "Lif", "Osm", "Cntf", "Lifr", "Osmr",
    # Other Cytokines
    "Tslp", "Mif", "Spp1", "Kitl", "Kit"
  )
}


# ============================================================================
# run_ripple()
# ============================================================================

#' Run RIPPLE distance correlation analysis (Stage 1)
#'
#' Detects genes whose expression changes as a function of distance from
#' a query cell type. Fits per-sample Poisson GLMs with cell-size offset,
#' then combines across replicates via Fisher's combined p-value with a
#' sign consistency gate.
#'
#' @param input_path Path to Seurat \code{.rds} file with raw counts in the
#'   \code{RNA} assay.
#' @param query_celltype Cell type label for the source population (e.g.,
#'   \code{"Tumor"}).
#' @param celltype_column Metadata column containing cell type annotations.
#' @param output_dir Output directory (default: \code{"."}).
#' @param sample_column Metadata column for sample/replicate IDs
#'   (default: \code{"sample_id"}).
#' @param condition_column Metadata column for condition filtering
#'   (default: \code{NULL} = no filter).
#' @param condition_value Which condition to analyze
#'   (default: \code{NULL} = all).
#' @param x_column X coordinate column (default: \code{NULL} = auto-detect).
#' @param y_column Y coordinate column (default: \code{NULL} = auto-detect).
#' @param target_celltypes Character vector of target cell types to analyze
#'   (default: \code{NULL} = auto-detect all non-query types).
#' @param k_neighbors Number of nearest query cells for distance calculation
#'   (default: \code{1}).
#' @param max_distance_um Maximum distance in micrometers to consider
#'   (default: \code{200}).
#' @param n_permutations Number of label permutations for null distribution
#'   (default: \code{0} = skip permutation testing).
#' @param fdr_threshold FDR cutoff for significance
#'   (default: \code{0.05}).
#' @param min_cells_per_sample Minimum cells of target type per sample
#'   (default: \code{30}).
#' @param min_expr_pct Minimum fraction of cells expressing a gene per sample
#'   for the strict filter tier (default: \code{0.01}).
#' @param min_expr_floor Absolute floor for expressing cells per sample; used
#'   as the threshold for the lenient (priority gene) filter tier and as the
#'   lower bound for the strict tier (default: \code{25}).
#' @param priority_genes Character vector of priority gene names that receive
#'   lenient expression filtering. Set to \code{NULL} to use the built-in
#'   chemokine/cytokine/receptor list, or \code{character(0)} to disable
#'   lenient filtering entirely.
#' @param query_signature_genes Character vector of query cell markers for
#'   contamination checking (default: \code{NULL} = no check).
#' @param query_label Display label for the query cell type in plots
#'   (default: same as \code{query_celltype}).
#' @param analysis_name Subdirectory name for output
#'   (default: \code{"ripple"}).
#' @param sign_consistency Minimum fraction of replicates that must agree on
#'   coefficient direction for Fisher's combined p-value to be computed. Set
#'   to \code{1.0} (the default) to require unanimous agreement.
#' @param verbose Print progress messages (default: \code{TRUE}).
#'
#' @return A \code{data.table} with columns: \code{gene}, \code{cell_type},
#'   \code{median_coef}, \code{fisher_pval}, \code{fisher_fdr},
#'   \code{sign_consistency}, \code{combined_coef}, \code{combined_se},
#'   \code{pval}, \code{fdr}, \code{gradient_score}, \code{decay_pattern},
#'   and per-sample coefficient information. Results are also written as CSV
#'   files to \code{output_dir}.
#'
#' @details
#' The analysis proceeds as follows:
#'
#' \enumerate{
#'   \item Load the Seurat object and extract raw counts + metadata.
#'   \item Resolve spatial coordinate columns (user-specified or auto-detected).
#'   \item Optionally filter to a specific condition.
#'   \item Compute distance from every cell to the nearest \code{k_neighbors}
#'     query cells.
#'   \item Auto-detect target cell types if not specified.
#'   \item For each target cell type:
#'     \enumerate{
#'       \item Filter genes using a two-tier system: strict for regular genes
#'         (>= max(1\% of cells, floor) expressing in ALL valid samples),
#'         lenient for priority genes (>= floor in >= 2 samples).
#'       \item Fit per-sample Poisson GLMs:
#'         \code{glm(counts ~ distance + offset(log(total_counts)), family = poisson)}.
#'       \item Run random-effects meta-analysis across samples.
#'       \item Compute Fisher's combined p-value with sign consistency gate.
#'       \item Optionally run permutation tests on top genes.
#'       \item Classify decay patterns for significant genes.
#'       \item Generate volcano plots, forest plots, and coefficient strips.
#'     }
#'   \item Merge results across cell types and write summary files.
#' }
#'
#' @section Coefficient interpretation:
#' \itemize{
#'   \item Negative coefficient: expression rate \emph{decreases} with distance
#'     = gene is \strong{induced} near the query cell type.
#'   \item Positive coefficient: expression rate \emph{increases} with distance
#'     = gene is \strong{repressed} near the query cell type.
#' }
#'
#' @export
run_ripple <- function(
  input_path,
  query_celltype,
  celltype_column,
  output_dir              = ".",
  sample_column           = "sample_id",
  condition_column        = NULL,
  condition_value         = NULL,
  x_column                = NULL,
  y_column                = NULL,
  target_celltypes        = NULL,
  k_neighbors             = 1,
  max_distance_um         = 200,
  n_permutations          = 0,
  fdr_threshold           = 0.05,
  min_cells_per_sample    = 30,
  min_expr_pct            = 0.01,
  min_expr_floor          = 25,
  priority_genes          = NULL,
  query_signature_genes   = NULL,
  query_label             = NULL,
  analysis_name           = "ripple",
  sign_consistency        = 1.0,
  verbose                 = TRUE
) {

  # --------------------------------------------------------------------------
  # 0. Resolve defaults from package options
  # --------------------------------------------------------------------------
  sample_column        <- .resolve(sample_column,        "sample_column",        "sample_id")
  condition_column     <- .resolve(condition_column,      "condition_column",     NULL)
  condition_value      <- .resolve(condition_value,       "condition_value",      NULL)
  k_neighbors          <- .resolve(k_neighbors,           "k_neighbors",          1L)
  max_distance_um      <- .resolve(max_distance_um,       "max_distance_um",      200)

  fdr_threshold        <- .resolve(fdr_threshold,         "fdr_threshold",        0.05)
  min_cells_per_sample <- .resolve(min_cells_per_sample,  "min_cells_per_sample", 30L)
  min_expr_floor       <- .resolve(min_expr_floor,        "min_expr_cells",       25L)
  min_expr_pct         <- .resolve(min_expr_pct,          "min_expr_pct",         0.01)
  sign_consistency     <- .resolve(sign_consistency,      "sign_consistency",     1.0)
  verbose              <- .resolve(verbose,               "verbose",              TRUE)

  if (is.null(query_label)) query_label <- query_celltype
  if (is.null(priority_genes)) priority_genes <- .default_priority_genes()

  # Minimum cells for stable GLM fit (internal constant)
  min_expr_cells_glm <- 5L

  set.seed(42)

  # --------------------------------------------------------------------------
  # 1. Validate inputs
  # --------------------------------------------------------------------------
  if (missing(input_path) || !file.exists(input_path)) {
    stop("input_path does not exist: ", input_path, call. = FALSE)
  }
  if (missing(query_celltype) || !nzchar(query_celltype)) {
    stop("query_celltype must be a non-empty string.", call. = FALSE)
  }
  if (missing(celltype_column) || !nzchar(celltype_column)) {
    stop("celltype_column must be a non-empty string.", call. = FALSE)
  }

  # --------------------------------------------------------------------------
  # 2. Set up output directories
  # --------------------------------------------------------------------------
  analysis_k_suffix <- if (k_neighbors > 1) paste0("_k", k_neighbors) else ""
  analysis_dir_name <- paste0(analysis_name, analysis_k_suffix)

  output_base <- file.path(output_dir, analysis_dir_name)
  .ensure_dir(output_base)
  .ensure_dir(file.path(output_base, "per_celltype"))
  .ensure_dir(file.path(output_base, "summary"))
  .ensure_dir(file.path(output_base, "plots"))
  .ensure_dir(file.path(output_base, "plots", "forest_plots"))
  .ensure_dir(file.path(output_base, "qc"))

  .msg(strrep("=", 70), verbose = verbose)
  .msg("RIPPLE Distance Correlation Analysis (Poisson GLM)", verbose = verbose)
  .msg(strrep("=", 70), verbose = verbose)
  .msg("Query cell type:  ", query_celltype, verbose = verbose)
  .msg("Cell type column: ", celltype_column, verbose = verbose)
  .msg("K neighbors:      ", k_neighbors, verbose = verbose)
  .msg("Output directory: ", output_base, verbose = verbose)

  # --------------------------------------------------------------------------
  # 3. Load and normalize data
  # --------------------------------------------------------------------------
  .msg("\nLoading data...", verbose = verbose)
  data <- .resolve_input(input_path, require_expr = FALSE, verbose = verbose)
  count_matrix_full <- data$counts
  cell_data <- data$meta
  rm(data)

  .msg("Counts verified: ", nrow(count_matrix_full), " genes x ",
       ncol(count_matrix_full), " cells", verbose = verbose)

  if (!celltype_column %in% names(cell_data)) {
    stop("Cell type column '", celltype_column, "' not found in metadata. ",
         "Available: ", paste(head(names(cell_data), 20), collapse = ", "),
         call. = FALSE)
  }

  # --------------------------------------------------------------------------
  # 5. Resolve coordinate columns
  # --------------------------------------------------------------------------
  coord_cols <- get_coord_columns(cell_data, x_col = x_column,
                                  y_col = y_column)
  .msg("Coordinate columns: ", coord_cols[1], ", ", coord_cols[2],
       verbose = verbose)

  # --------------------------------------------------------------------------
  # 6. Resolve condition column and filter
  # --------------------------------------------------------------------------
  if (!is.null(condition_column) && nzchar(condition_column)) {
    if (!condition_column %in% names(cell_data)) {
      stop("Condition column '", condition_column, "' not found in metadata.",
           call. = FALSE)
    }
    cell_data[, condition := get(condition_column)]
  } else {
    cell_data[, condition := "all"]
  }

  .msg("\nData summary (before filtering):", verbose = verbose)
  .msg("  Total cells: ", nrow(cell_data), verbose = verbose)
  .msg("  Samples: ", data.table::uniqueN(cell_data[[sample_column]]),
       verbose = verbose)
  .msg("  Conditions: ",
       paste(unique(cell_data$condition), collapse = ", "),
       verbose = verbose)

  if (!is.null(condition_value) && nzchar(condition_value)) {
    .msg("\nFiltering to condition == '", condition_value, "'...",
         verbose = verbose)
    cell_data <- cell_data[condition == condition_value]
    if (nrow(cell_data) == 0) {
      stop("No cells remain after filtering to condition '", condition_value,
           "'.", call. = FALSE)
    }
  }

  .msg("Data summary (after filtering):", verbose = verbose)
  .msg("  Total cells: ", nrow(cell_data), verbose = verbose)
  .msg("  Samples: ", data.table::uniqueN(cell_data[[sample_column]]),
       verbose = verbose)

  # Subset count matrix to filtered barcodes
  count_matrix <- count_matrix_full[, cell_data$barcode, drop = FALSE]

  # Total counts per cell (for Poisson offset)
  total_counts <- colSums(count_matrix_full[, cell_data$barcode, drop = FALSE])
  total_counts <- stats::setNames(as.numeric(total_counts), cell_data$barcode)

  .msg("Total counts per cell: median=", round(stats::median(total_counts)),
       ", range=[", round(min(total_counts)), "-",
       round(max(total_counts)), "]", verbose = verbose)

  # Free full count matrix
  rm(count_matrix_full)
  gc(verbose = FALSE)

  # --------------------------------------------------------------------------
  # 7. Calculate distances to query cells
  # --------------------------------------------------------------------------
  .msg("\n", strrep("=", 70), verbose = verbose)
  .msg("Calculating Distances to Query Cells (k=", k_neighbors, ")",
       verbose = verbose)
  .msg(strrep("=", 70), verbose = verbose)

  coords <- as.matrix(cell_data[, ..coord_cols])
  sample_ids_all <- cell_data[[sample_column]]

  query_mask <- cell_data[[celltype_column]] == query_celltype
  n_query <- sum(query_mask)
  .msg("Query cells (", query_celltype, "): ", n_query, verbose = verbose)

  if (n_query < 10) {
    stop("Too few query cells for analysis (", n_query, "). ",
         "Need at least 10.", call. = FALSE)
  }

  # Query cells per sample (for stratified permutation)
  query_per_sample_dt <- cell_data[query_mask == TRUE, .N,
                                   by = c(sample_column)]
  query_per_sample <- stats::setNames(query_per_sample_dt$N,
                                      query_per_sample_dt[[sample_column]])
  .msg("Query cells per sample:", verbose = verbose)
  for (samp in names(query_per_sample)) {
    .msg("  ", samp, ": ", query_per_sample[samp], verbose = verbose)
  }

  query_coords <- coords[query_mask, , drop = FALSE]

  # k-NN distance calculation
  effective_k <- min(k_neighbors, nrow(query_coords))
  .msg("Computing ", effective_k, "-NN distances...", verbose = verbose)
  nn_result <- RANN::nn2(query_coords, coords, k = effective_k)

  if (effective_k == 1) {
    cell_data[, dist_to_query := as.vector(nn_result$nn.dists)]
  } else {
    cell_data[, dist_to_query := rowMeans(nn_result$nn.dists)]
  }

  # Cap distances
  cell_data[dist_to_query > max_distance_um, dist_to_query := max_distance_um]

  .msg("Distance distribution:", verbose = verbose)
  .msg("  Min: ", round(min(cell_data$dist_to_query), 1), " um",
       verbose = verbose)
  .msg("  Median: ", round(stats::median(cell_data$dist_to_query), 1), " um",
       verbose = verbose)
  .msg("  Max: ", round(max(cell_data$dist_to_query), 1), " um",
       verbose = verbose)

  # --------------------------------------------------------------------------
  # 8. Determine target cell types
  # --------------------------------------------------------------------------
  if (is.null(target_celltypes) || length(target_celltypes) == 0) {
    all_types <- unique(cell_data[[celltype_column]])
    target_celltypes <- sort(setdiff(all_types, c(query_celltype, NA_character_)))
    .msg("Auto-detected ", length(target_celltypes), " target cell types: ",
         paste(target_celltypes, collapse = ", "), verbose = verbose)
  }

  # --------------------------------------------------------------------------
  # 9. QC diagnostics
  # --------------------------------------------------------------------------
  .msg("\n", strrep("=", 70), verbose = verbose)
  .msg("QC Diagnostics", verbose = verbose)
  .msg(strrep("=", 70), verbose = verbose)

  qc_dir <- file.path(output_base, "qc")

  # Distance distribution histogram
  p_dist <- ggplot2::ggplot(cell_data, ggplot2::aes(x = dist_to_query)) +
    ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "white",
                            alpha = 0.7) +
    ggplot2::geom_vline(xintercept = 50, linetype = "dashed", color = "red",
                        linewidth = 1) +
    ggplot2::labs(
      title = sprintf("Distance to Nearest %s (k=%d)", query_celltype,
                       k_neighbors),
      x = "Distance (um)",
      y = "Number of Cells"
    ) +
    ggplot2::theme_bw(base_size = 12)
  ggplot2::ggsave(file.path(qc_dir, "distance_distribution.pdf"),
                  p_dist, width = 8, height = 5)
  .msg("  Saved: qc/distance_distribution.pdf", verbose = verbose)

  # Per-sample summary
  sample_summary <- cell_data[, .(
    n_total = .N,
    n_query = sum(get(celltype_column) == query_celltype),
    median_dist = stats::median(dist_to_query),
    pct_within_50um = mean(dist_to_query <= 50) * 100,
    n_cell_types = data.table::uniqueN(get(celltype_column))
  ), by = c(sample_column, "condition")]
  data.table::fwrite(sample_summary, file.path(qc_dir, "sample_summary.csv"))
  .msg("  Saved: qc/sample_summary.csv", verbose = verbose)

  # Cell type counts
  celltype_counts <- cell_data[, .N, by = c(celltype_column)]
  data.table::setnames(celltype_counts, celltype_column, "cell_type")
  data.table::setorder(celltype_counts, -N)
  data.table::fwrite(celltype_counts,
                     file.path(qc_dir, "filtered_celltype_counts.csv"))

  # --------------------------------------------------------------------------
  # 10. Main analysis loop: iterate over target cell types
  # --------------------------------------------------------------------------
  .msg("\n", strrep("=", 70), verbose = verbose)
  .msg("Running Distance Correlation Analysis (Poisson GLM)",
       verbose = verbose)
  .msg(strrep("=", 70), verbose = verbose)

  all_results <- list()

  for (ct_name in target_celltypes) {

    .msg("\n", strrep("-", 60), verbose = verbose)
    .msg("Analyzing: ", ct_name, verbose = verbose)
    .msg(strrep("-", 60), verbose = verbose)

    ct_output_dir <- file.path(output_base, "per_celltype", ct_name)
    .ensure_dir(ct_output_dir)

    # Identify target cells
    cell_data[, is_target := get(celltype_column) == ct_name]
    target_data <- cell_data[is_target == TRUE]

    .msg("  Target cells: ", nrow(target_data), verbose = verbose)

    if (nrow(target_data) < min_cells_per_sample * 2) {
      .msg("  Insufficient cells for analysis", verbose = verbose)
      next
    }

    # Valid samples
    sample_counts_dt <- target_data[, .N, by = c(sample_column)]
    valid_samples <- sample_counts_dt[N >= min_cells_per_sample][[sample_column]]
    .msg("  Valid samples (>= ", min_cells_per_sample, " cells): ",
         length(valid_samples), verbose = verbose)

    if (length(valid_samples) < 2) {
      .msg("  Need at least 2 valid samples for meta-analysis",
           verbose = verbose)
      next
    }
    if (length(valid_samples) == 2) {
      .msg("  WARNING: Only 2 valid samples - meta-analysis will have low power",
           verbose = verbose)
    }

    # Target cell barcodes in valid samples
    target_valid <- target_data[get(sample_column) %in% valid_samples]
    target_barcodes <- target_valid$barcode

    # Raw counts for target cells
    target_counts <- count_matrix[, target_barcodes, drop = FALSE]

    # --------------------------------------------------------------------
    # Two-tier gene expression filtering
    # --------------------------------------------------------------------
    sample_ids_for_filter <- droplevels(as.factor(target_valid[[sample_column]]))
    cells_per_sample <- table(sample_ids_for_filter)
    threshold_per_sample <- pmax(ceiling(cells_per_sample * min_expr_pct),
                                 min_expr_floor)
    threshold_floor <- stats::setNames(
      rep(min_expr_floor, length(cells_per_sample)),
      names(cells_per_sample)
    )

    .msg("  Per-sample expression thresholds (max(", min_expr_pct * 100,
         "%, ", min_expr_floor, ")):", verbose = verbose)
    for (s in names(threshold_per_sample)) {
      .msg("    ", s, ": ", threshold_per_sample[s], " cells (of ",
           cells_per_sample[s], ")", verbose = verbose)
    }

    # Count expressing cells per gene per sample
    expressing_counts <- sapply(rownames(target_counts), function(g) {
      count_vec <- target_counts[g, ]
      tapply(count_vec > 0, sample_ids_for_filter, sum)
    })

    # Tier 1: strict filter (all samples must pass)
    genes_pass_strict <- apply(expressing_counts, 2, function(counts) {
      all(counts >= threshold_per_sample[names(counts)], na.rm = TRUE)
    })
    genes_strict <- names(genes_pass_strict[genes_pass_strict])

    # Tier 2: lenient filter (priority genes only)
    genes_lenient_only <- character(0)
    if (length(priority_genes) > 0) {
      priority_in_panel <- intersect(priority_genes, rownames(target_counts))
      if (length(priority_in_panel) > 0) {
        genes_pass_lenient <- apply(
          expressing_counts[, priority_in_panel, drop = FALSE], 2,
          function(counts) {
            sum(counts >= threshold_floor[names(counts)], na.rm = TRUE) >= 2
          }
        )
        genes_lenient <- names(genes_pass_lenient[genes_pass_lenient])
        genes_lenient_only <- setdiff(genes_lenient, genes_strict)
      }
    }

    genes_to_analyze <- union(genes_strict, genes_lenient_only)
    .msg("  Genes passing strict filter: ", length(genes_strict),
         verbose = verbose)
    if (length(genes_lenient_only) > 0) {
      .msg("  Priority genes rescued by lenient filter: ",
           length(genes_lenient_only),
           " (", paste(head(genes_lenient_only, 10), collapse = ", "),
           if (length(genes_lenient_only) > 10) ", ..." else "", ")",
           verbose = verbose)
    }
    .msg("  Total genes to analyze: ", length(genes_to_analyze),
         verbose = verbose)

    if (length(genes_to_analyze) < 10) {
      .msg("  Too few genes passed filtering", verbose = verbose)
      next
    }

    # --------------------------------------------------------------------
    # Step 1: Per-sample Poisson GLM coefficients
    # --------------------------------------------------------------------
    .msg("  Step 1: Calculating per-sample Poisson coefficients...",
         verbose = verbose)

    coef_results <- data.table::rbindlist(lapply(genes_to_analyze, function(g) {
      gene_counts <- as.numeric(target_counts[g, target_barcodes])

      sample_results <- data.table::rbindlist(lapply(valid_samples, function(samp) {
        samp_idx <- which(target_valid[[sample_column]] == samp)
        if (length(samp_idx) < min_cells_per_sample) {
          return(data.table::data.table(
            gene = g, sample_id = samp,
            coef = NA_real_, se = NA_real_,
            n_cells = length(samp_idx), pval = NA_real_,
            dispersion = NA_real_
          ))
        }

        samp_counts <- gene_counts[samp_idx]
        samp_dist <- target_valid[samp_idx]$dist_to_query
        samp_total <- total_counts[target_barcodes[samp_idx]]

        fit_result <- fit_poisson(samp_counts, samp_dist, samp_total,
                                  min_cells = min_expr_cells_glm)

        if (is.null(fit_result)) {
          data.table::data.table(
            gene = g, sample_id = samp,
            coef = NA_real_, se = NA_real_,
            n_cells = length(samp_idx), pval = NA_real_,
            dispersion = NA_real_
          )
        } else {
          data.table::data.table(
            gene = g, sample_id = samp,
            coef = fit_result$beta, se = fit_result$se,
            n_cells = fit_result$n_cells, pval = fit_result$pval,
            dispersion = fit_result$dispersion
          )
        }
      }))

      sample_results
    }), fill = TRUE)

    # Save per-sample coefficients
    data.table::fwrite(coef_results,
                       file.path(ct_output_dir, "coef_per_sample.csv"))
    .msg("  Saved: coef_per_sample.csv", verbose = verbose)

    # --------------------------------------------------------------------
    # Step 2: Meta-analysis across samples
    # --------------------------------------------------------------------
    .msg("  Step 2: Running meta-analysis...", verbose = verbose)

    meta_results <- data.table::rbindlist(lapply(genes_to_analyze, function(g) {
      gene_data <- coef_results[gene == g]

      meta_result <- run_meta_analysis(
        coefs = gene_data$coef,
        ses = gene_data$se,
        sample_ids = gene_data[["sample_id"]]
      )

      # Expression statistics
      gene_counts <- as.numeric(target_counts[g, target_barcodes])

      # Sign consistency
      valid_coefs <- gene_data$coef[!is.na(gene_data$coef)]
      n_pos <- sum(valid_coefs > 0)
      n_neg <- sum(valid_coefs < 0)
      n_valid <- n_pos + n_neg
      sign_con <- if (n_valid > 0) max(n_pos, n_neg) / n_valid else NA_real_

      # Median dispersion
      med_disp <- stats::median(gene_data$dispersion, na.rm = TRUE)

      data.table::data.table(
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

    # BH-adjusted p-value
    meta_results[, fdr := stats::p.adjust(pval, method = "BH")]

    # --------------------------------------------------------------------
    # Step 2b: Fisher's combined p-value
    # --------------------------------------------------------------------
    .msg("  Step 2b: Computing Fisher's combined p-values...",
         verbose = verbose)

    fisher_results <- data.table::rbindlist(lapply(genes_to_analyze, function(g) {
      gene_data <- coef_results[gene == g]
      fisher <- compute_fisher_pval(
        pvals = gene_data$pval,
        coefs = gene_data$coef,
        min_samples = 2L,
        sign_threshold = sign_consistency
      )
      data.table::data.table(
        gene = g,
        median_coef = fisher$median_coef,
        fisher_stat = fisher$fisher_stat,
        fisher_pval = fisher$fisher_pval
      )
    }), fill = TRUE)

    fisher_results[, fisher_fdr := stats::p.adjust(fisher_pval, method = "BH")]
    meta_results <- merge(meta_results, fisher_results, by = "gene",
                          all.x = TRUE)
    data.table::setDT(meta_results)

    # --------------------------------------------------------------------
    # Step 2c: Permutation testing (optional)
    # --------------------------------------------------------------------
    if (n_permutations > 0) {
      .msg("  Step 2c: Running permutation tests (",
           n_permutations, " permutations)...", verbose = verbose)

      perm_top_n <- 100
      top_by_effect <- meta_results[
        order(-abs(combined_coef))][1:min(perm_top_n, .N)]$gene
      priority_in_data <- intersect(priority_genes, genes_to_analyze)
      top_genes_for_perm <- unique(c(top_by_effect, priority_in_data))

      .msg(sprintf("    Permutation genes: %d total (after dedup)",
                    length(top_genes_for_perm)), verbose = verbose)

      coords_target <- as.matrix(target_valid[, ..coord_cols])
      observed_coefs <- stats::setNames(
        meta_results[gene %in% top_genes_for_perm]$combined_coef,
        meta_results[gene %in% top_genes_for_perm]$gene
      )
      query_per_sample_valid <- query_per_sample[
        names(query_per_sample) %in% valid_samples]

      perm_results <- run_permutation_tests(
        genes = top_genes_for_perm,
        count_matrix = target_counts,
        target_barcodes = target_barcodes,
        coords_target = coords_target,
        coords_all = coords,
        sample_ids_target = target_valid[[sample_column]],
        sample_ids_all = sample_ids_all,
        query_per_sample = query_per_sample_valid,
        observed_coefs = observed_coefs,
        n_perms = n_permutations,
        k_neighbors = k_neighbors,
        max_distance_um = max_distance_um,
        min_cells_per_sample = min_cells_per_sample,
        min_expr_cells = min_expr_cells_glm,
        total_counts_target = total_counts[target_barcodes]
      )

      meta_results <- merge(meta_results, perm_results, by = "gene",
                            all.x = TRUE)
      data.table::setDT(meta_results)

      data.table::fwrite(perm_results,
                         file.path(ct_output_dir, "permutation_pvals.csv"))
      .msg("    Saved: permutation_pvals.csv", verbose = verbose)
    } else {
      meta_results[, perm_pval := NA_real_]
    }

    # --------------------------------------------------------------------
    # Step 3: Gradient scores
    # --------------------------------------------------------------------
    meta_results[, gradient_score := combined_coef]

    # --------------------------------------------------------------------
    # Step 4: Decay pattern classification (significant genes)
    # --------------------------------------------------------------------
    .msg("  Step 4: Classifying decay patterns...", verbose = verbose)

    sig_genes <- meta_results[fisher_fdr < fdr_threshold]$gene
    .msg("  Significant genes (Fisher FDR < ", fdr_threshold, "): ",
         length(sig_genes), verbose = verbose)

    if (length(sig_genes) > 0) {
      decay_patterns <- sapply(sig_genes, function(gene) {
        gene_counts <- as.numeric(target_counts[gene, target_barcodes])

        per_sample_patterns <- sapply(valid_samples, function(samp) {
          samp_idx <- which(target_valid[[sample_column]] == samp)
          if (length(samp_idx) < min_cells_per_sample) return(NA_character_)

          samp_counts <- gene_counts[samp_idx]
          samp_dist <- target_valid[samp_idx]$dist_to_query
          samp_total <- total_counts[target_barcodes[samp_idx]]

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
      meta_results[gene %in% sig_genes,
                   decay_pattern := decay_patterns[match(gene, sig_genes)]]
    } else {
      meta_results[, decay_pattern := "not_significant"]
    }

    # Save per-celltype results
    data.table::fwrite(meta_results,
                       file.path(ct_output_dir, "meta_analysis_results.csv"))
    .msg("  Saved: meta_analysis_results.csv", verbose = verbose)

    # Gradient scores CSV
    gradient_results <- meta_results[, .(gene, gradient_score, combined_coef,
                                          fdr, fisher_fdr, decay_pattern,
                                          sign_consistency,
                                          median_dispersion)]
    gradient_results <- gradient_results[order(gradient_score)]
    data.table::fwrite(gradient_results,
                       file.path(ct_output_dir, "gradient_scores.csv"))

    # Decay classification CSV
    decay_summary <- meta_results[fisher_fdr < fdr_threshold,
                                  .N, by = decay_pattern]
    data.table::fwrite(decay_summary,
                       file.path(ct_output_dir, "decay_classification.csv"))

    # --------------------------------------------------------------------
    # Step 5: Visualizations
    # --------------------------------------------------------------------
    .msg("  Step 5: Creating visualizations...", verbose = verbose)

    # Volcano plot
    create_gradient_volcano(meta_results, ct_name,
                             file.path(ct_output_dir, "gradient_volcano.pdf"),
                             fdr_threshold = fdr_threshold,
                             query_label = query_label,
                             k_neighbors = k_neighbors)

    # Forest plots and coefficient strips for top significant genes
    top_sig_genes <- head(
      meta_results[fisher_fdr < fdr_threshold][order(fisher_fdr)]$gene, 20
    )

    if (length(top_sig_genes) > 0) {
      forest_dir <- file.path(ct_output_dir, "forest_plots")
      .ensure_dir(forest_dir)

      for (g in top_sig_genes) {
        gene_data <- coef_results[gene == g]
        create_forest_plot(
          coefs = gene_data$coef,
          ses = gene_data$se,
          sample_ids = gene_data[["sample_id"]],
          gene = g,
          cell_type = ct_name,
          output_path = file.path(forest_dir,
                                  sprintf("%s_forest.pdf", g)),
          query_label = query_label
        )
      }

      # Coefficient strip plot
      create_coefficient_strips(
        coef_results = coef_results,
        meta_results = meta_results,
        cell_type = ct_name,
        output_path = file.path(ct_output_dir, "coefficient_strips.pdf"),
        fdr_threshold = fdr_threshold
      )
    }

    # Add cell type label
    meta_results[, cell_type := ct_name]

    all_results[[ct_name]] <- meta_results
    .msg("  Analysis complete for ", ct_name, verbose = verbose)
  }

  # --------------------------------------------------------------------------
  # 11. Merge results across cell types
  # --------------------------------------------------------------------------
  .msg("\n", strrep("=", 70), verbose = verbose)
  .msg("Creating Summary", verbose = verbose)
  .msg(strrep("=", 70), verbose = verbose)

  if (length(all_results) == 0) {
    .msg("No cell types had sufficient data for analysis.", verbose = verbose)
    return(invisible(data.table::data.table()))
  }

  combined_results <- data.table::rbindlist(all_results, fill = TRUE)

  data.table::fwrite(combined_results,
                     file.path(output_base, "summary", "all_genes_results.csv"))
  .msg("Saved: summary/all_genes_results.csv", verbose = verbose)

  # Top genes by Fisher FDR
  if ("fisher_fdr" %in% names(combined_results)) {
    top_genes <- combined_results[fisher_fdr < fdr_threshold][
      order(fisher_fdr), head(.SD, 50), by = cell_type
    ]
    data.table::fwrite(top_genes,
                       file.path(output_base, "summary",
                                 "top_gradient_genes.csv"))
    .msg("Saved: summary/top_gradient_genes.csv", verbose = verbose)
  }

  # Decay pattern summary
  decay_summary_all <- combined_results[fisher_fdr < fdr_threshold,
                                        .N, by = .(cell_type, decay_pattern)]
  if (nrow(decay_summary_all) > 0) {
    data.table::fwrite(decay_summary_all,
                       file.path(output_base, "summary",
                                 "decay_pattern_summary.csv"))
    .msg("Saved: summary/decay_pattern_summary.csv", verbose = verbose)
  }

  # Dispersion summary
  disp_summary <- combined_results[, .(
    median_dispersion = stats::median(median_dispersion, na.rm = TRUE),
    pct_overdispersed = mean(median_dispersion > 2, na.rm = TRUE) * 100
  ), by = cell_type]
  data.table::fwrite(disp_summary,
                     file.path(output_base, "summary",
                               "dispersion_summary.csv"))

  # Sign consistency summary
  sign_summary <- combined_results[fisher_fdr < fdr_threshold, .(
    n_sig = .N,
    pct_all_agree = mean(sign_consistency == 1, na.rm = TRUE) * 100,
    median_sign_consistency = stats::median(sign_consistency, na.rm = TRUE)
  ), by = cell_type]
  data.table::fwrite(sign_summary,
                     file.path(output_base, "summary",
                               "sign_consistency_summary.csv"))

  # Print summary
  .msg("\n", strrep("=", 70), verbose = verbose)
  .msg("Analysis Summary", verbose = verbose)
  .msg(strrep("=", 70), verbose = verbose)
  for (ct in names(all_results)) {
    ct_result <- all_results[[ct]]
    n_sig <- sum(ct_result$fisher_fdr < fdr_threshold, na.rm = TRUE)
    n_neg <- sum(ct_result$fisher_fdr < fdr_threshold &
                   ct_result$gradient_score < 0, na.rm = TRUE)
    n_pos <- sum(ct_result$fisher_fdr < fdr_threshold &
                   ct_result$gradient_score > 0, na.rm = TRUE)
    .msg(sprintf("  %s: %d significant genes (%d %s-induced, %d %s-repressed)",
                 ct, n_sig, n_neg, query_label, n_pos, query_label),
         verbose = verbose)
  }

  .msg("\n", strrep("=", 70), verbose = verbose)
  .msg("Analysis Complete!", verbose = verbose)
  .msg("Output directory: ", output_base, verbose = verbose)
  .msg(strrep("=", 70), verbose = verbose)

  return(invisible(combined_results))
}


# ============================================================================
# run_ripple_confounder()
# ============================================================================

#' Run RIPPLE confounder control analysis (Stage 4)
#'
#' Validates Stage 1 results by adding distance to a control cell type as a
#' covariate in a bivariate Poisson GLM. For each significant gene from
#' Stage 1, fits:
#'
#' \code{glm(counts ~ dist_to_query + dist_to_control + offset(log(total_counts)),
#'     family = poisson)}
#'
#' Classifies genes as \code{query_specific}, \code{enhanced},
#' \code{niche_driven}, or \code{underpowered} based on whether the
#' query-distance coefficient persists after controlling for the shared
#' niche effect.
#'
#' @param input_path Path to Seurat \code{.rds} file.
#' @param results_dir Path to Stage 1 results directory (the parent directory
#'   containing \code{summary/all_genes_results.csv}).
#' @param query_celltype Query cell type label.
#' @param celltype_column Cell type column name.
#' @param control_celltype Control cell type for confounder analysis (e.g.,
#'   \code{"Monocyte"}).
#' @param output_dir Output directory for Stage 4 results (default: same
#'   parent as \code{results_dir}, with \code{"_stage2"} suffix appended to
#'   the analysis directory name).
#' @param sample_column Metadata column for sample/replicate IDs
#'   (default: \code{"sample_id"}).
#' @param condition_column Metadata column for condition filtering
#'   (default: \code{NULL}).
#' @param condition_value Which condition to analyze (default: \code{NULL}).
#' @param x_column X coordinate column (default: \code{NULL} = auto-detect).
#' @param y_column Y coordinate column (default: \code{NULL} = auto-detect).
#' @param target_celltypes Character vector of target cell types
#'   (default: \code{NULL} = use cell types present in Stage 1 results).
#' @param fdr_threshold FDR cutoff (default: \code{0.05}).
#' @param min_cells_per_sample Minimum cells per sample (default: \code{30}).
#' @param min_control_cells Minimum control cells per sample for reliable
#'   distance calculation (default: \code{30}).
#' @param max_distance_um Maximum distance in micrometers (default: \code{200}).
#' @param sig_column Which column from Stage 1 results to use for selecting
#'   significant genes. Typically \code{"fisher_fdr"} (default) or
#'   \code{"fdr"}.
#' @param query_label Display label for query cell type (default: same as
#'   \code{query_celltype}).
#' @param sign_consistency Sign consistency threshold for Fisher's p-value
#'   (default: \code{1.0}).
#' @param verbose Print progress messages (default: \code{TRUE}).
#'
#' @return A \code{data.table} with Stage 1 vs Stage 2 comparison and
#'   classification for each gene, including columns: \code{gene},
#'   \code{cell_type}, \code{stage1_coef}, \code{stage2_median_coef},
#'   \code{stage2_fisher_fdr}, \code{classification},
#'   \code{control_celltype}.
#'
#' @section Gene classification:
#' \describe{
#'   \item{query_specific}{Fisher FDR < threshold in Stage 4, same sign as
#'     Stage 1. The gradient persists after controlling for the niche.}
#'   \item{enhanced}{query_specific AND absolute coefficient is > 1.1x
#'     Stage 1 value. The control was suppressing the signal.}
#'   \item{niche_driven}{FDR >= threshold AND coefficient attenuated > 50\%.
#'     The gradient was explained by the shared niche.}
#'   \item{underpowered}{FDR >= threshold AND coefficient preserved >= 50\%.
#'     Likely a power issue from SE inflation, not a true niche effect.}
#' }
#'
#' @export
run_ripple_confounder <- function(
  input_path,
  results_dir,
  query_celltype,
  celltype_column,
  control_celltype,
  output_dir              = NULL,
  sample_column           = "sample_id",
  condition_column        = NULL,
  condition_value         = NULL,
  x_column                = NULL,
  y_column                = NULL,
  target_celltypes        = NULL,
  fdr_threshold           = 0.05,
  min_cells_per_sample    = 30,
  min_control_cells       = 30,
  max_distance_um         = 200,
  sig_column              = "fisher_fdr",
  query_label             = NULL,
  sign_consistency        = 1.0,
  verbose                 = TRUE
) {

  # --------------------------------------------------------------------------
  # 0. Resolve defaults
  # --------------------------------------------------------------------------
  sample_column        <- .resolve(sample_column,        "sample_column",        "sample_id")
  condition_column     <- .resolve(condition_column,      "condition_column",     NULL)
  condition_value      <- .resolve(condition_value,       "condition_value",      NULL)
  fdr_threshold        <- .resolve(fdr_threshold,         "fdr_threshold",        0.05)
  min_cells_per_sample <- .resolve(min_cells_per_sample,  "min_cells_per_sample", 30L)
  sign_consistency     <- .resolve(sign_consistency,      "sign_consistency",     1.0)
  max_distance_um      <- .resolve(max_distance_um,       "max_distance_um",      200)
  verbose              <- .resolve(verbose,               "verbose",              TRUE)

  if (is.null(query_label)) query_label <- query_celltype
  min_expr_cells_glm <- 5L

  set.seed(42)

  # --------------------------------------------------------------------------
  # 1. Validate inputs
  # --------------------------------------------------------------------------
  if (missing(input_path) || !file.exists(input_path)) {
    stop("input_path does not exist: ", input_path, call. = FALSE)
  }
  if (missing(control_celltype) || !nzchar(control_celltype)) {
    stop("control_celltype is required for confounder analysis.",
         call. = FALSE)
  }

  # --------------------------------------------------------------------------
  # 2. Load Stage 1 results
  # --------------------------------------------------------------------------
  stage1_summary_file <- file.path(results_dir, "summary",
                                   "all_genes_results.csv")
  if (!file.exists(stage1_summary_file)) {
    stop("Stage 1 results not found: ", stage1_summary_file,
         "\nRun run_ripple() or merge_ripple_results() first.",
         call. = FALSE)
  }

  stage1_all <- data.table::fread(stage1_summary_file)
  .msg("Loaded Stage 1 results: ", nrow(stage1_all), " gene-celltype entries",
       verbose = verbose)

  # Verify significance column
  if (!sig_column %in% names(stage1_all)) {
    .msg("WARNING: '", sig_column, "' not found. Falling back to 'fdr'.",
         verbose = verbose)
    sig_column <- "fdr"
    if (!sig_column %in% names(stage1_all)) {
      stop("Neither 'fisher_fdr' nor 'fdr' found in Stage 1 results.",
           call. = FALSE)
    }
  }

  stage1_sig <- stage1_all[get(sig_column) < fdr_threshold]
  .msg("Stage 1 significant (", sig_column, " < ", fdr_threshold, "): ",
       nrow(stage1_sig), " entries", verbose = verbose)

  if (nrow(stage1_sig) == 0) {
    .msg("No significant genes from Stage 1. Nothing to validate.",
         verbose = verbose)
    return(invisible(data.table::data.table()))
  }

  # --------------------------------------------------------------------------
  # 3. Set up output directory
  # --------------------------------------------------------------------------
  if (is.null(output_dir)) {
    output_dir <- paste0(results_dir, "_stage2")
  }
  .ensure_dir(output_dir)
  .ensure_dir(file.path(output_dir, "per_celltype"))
  .ensure_dir(file.path(output_dir, "summary"))
  .ensure_dir(file.path(output_dir, "plots"))

  .msg(strrep("=", 70), verbose = verbose)
  .msg("RIPPLE Stage 4: Confounder Control (Bivariate Poisson GLM)",
       verbose = verbose)
  .msg(strrep("=", 70), verbose = verbose)
  .msg("Query cell type:   ", query_celltype, verbose = verbose)
  .msg("Control cell type: ", control_celltype, verbose = verbose)
  .msg("Output directory:  ", output_dir, verbose = verbose)

  # --------------------------------------------------------------------------
  # 4. Load and normalize data
  # --------------------------------------------------------------------------
  .msg("\nLoading data...", verbose = verbose)
  data <- .resolve_input(input_path, require_expr = FALSE, verbose = verbose)
  count_matrix_full <- data$counts
  cell_data <- data$meta
  rm(data)

  if (!celltype_column %in% names(cell_data)) {
    stop("Cell type column '", celltype_column, "' not found in metadata.",
         call. = FALSE)
  }

  # Resolve condition
  if (!is.null(condition_column) && nzchar(condition_column)) {
    if (!condition_column %in% names(cell_data)) {
      stop("Condition column '", condition_column, "' not found.",
           call. = FALSE)
    }
    cell_data[, condition := get(condition_column)]
  } else {
    cell_data[, condition := "all"]
  }

  # Filter by condition
  if (!is.null(condition_value) && nzchar(condition_value)) {
    cell_data <- cell_data[condition == condition_value]
  }

  .msg("Data after filtering: ", nrow(cell_data), " cells",
       verbose = verbose)

  # --------------------------------------------------------------------------
  # 5. Calculate distances
  # --------------------------------------------------------------------------
  coord_cols <- get_coord_columns(cell_data, x_col = x_column,
                                    y_col = y_column)
  coords <- as.matrix(cell_data[, ..coord_cols])

  # Distance to query
  query_mask <- cell_data[[celltype_column]] == query_celltype
  .msg("Query cells (", query_celltype, "): ", sum(query_mask),
       verbose = verbose)
  if (sum(query_mask) < 10) {
    stop("Too few query cells (", sum(query_mask), ").", call. = FALSE)
  }
  query_coords <- coords[query_mask, , drop = FALSE]
  nn_query <- RANN::nn2(query_coords, coords, k = 1)
  cell_data[, dist_to_query := pmin(as.vector(nn_query$nn.dists),
                                     max_distance_um)]

  # Distance to control
  control_mask <- cell_data[[celltype_column]] == control_celltype
  n_control_total <- sum(control_mask)
  .msg("Control cells (", control_celltype, "): ", n_control_total,
       verbose = verbose)

  if (n_control_total < min_control_cells) {
    stop("Too few control cells (", n_control_total, ") for reliable ",
         "distance calculation.", call. = FALSE)
  }

  control_coords <- coords[control_mask, , drop = FALSE]
  nn_control <- RANN::nn2(control_coords, coords, k = 1)
  cell_data[, dist_to_control := pmin(as.vector(nn_control$nn.dists),
                                       max_distance_um)]

  # Per-sample control cell counts
  control_per_sample <- cell_data[control_mask == TRUE, .N,
                                  by = c(sample_column)]
  data.table::setnames(control_per_sample, "N", "n_control")

  # Collinearity check
  .msg("\nCollinearity Diagnostics:", verbose = verbose)
  samples_all <- unique(cell_data[[sample_column]])
  for (samp in samples_all) {
    samp_data <- cell_data[get(sample_column) == samp]
    cor_val <- stats::cor(samp_data$dist_to_query, samp_data$dist_to_control,
                          use = "complete.obs", method = "pearson")
    flag <- if (abs(cor_val) > 0.8) " [WARNING: high collinearity]" else ""
    .msg("  ", samp, ": r = ", round(cor_val, 3), flag, verbose = verbose)
  }

  # --------------------------------------------------------------------------
  # 6. Determine target cell types from Stage 1 results
  # --------------------------------------------------------------------------
  if (is.null(target_celltypes) || length(target_celltypes) == 0) {
    target_celltypes <- unique(stage1_sig$cell_type)
    target_celltypes <- setdiff(target_celltypes, c(NA_character_))
  }

  # --------------------------------------------------------------------------
  # 7. Bivariate analysis per cell type
  # --------------------------------------------------------------------------
  all_stage2_results <- list()
  query_specific_label <- paste0(query_label, "_specific")

  for (ct_name in target_celltypes) {
    .msg("\n", strrep("-", 60), verbose = verbose)
    .msg("Analyzing: ", ct_name, verbose = verbose)

    # If the target IS the control, skip
    if (ct_name == control_celltype) {
      .msg("  SKIPPING: target == control cell type", verbose = verbose)
      next
    }

    ct_stage1 <- stage1_sig[cell_type == ct_name]
    if (nrow(ct_stage1) == 0) {
      .msg("  No Stage 1 significant genes for ", ct_name, verbose = verbose)
      next
    }
    sig_genes <- ct_stage1$gene

    # Target cells
    cell_data[, is_target := get(celltype_column) == ct_name]
    target_data <- cell_data[is_target == TRUE]
    .msg("  Target cells: ", nrow(target_data), verbose = verbose)

    if (nrow(target_data) < min_cells_per_sample * 2) {
      .msg("  Insufficient cells", verbose = verbose)
      next
    }

    # Valid samples (enough target AND control cells)
    samp_counts_dt <- target_data[, .N, by = c(sample_column)]
    valid_samples <- samp_counts_dt[N >= min_cells_per_sample][[sample_column]]
    valid_ctrl_samples <- control_per_sample[
      n_control >= min_control_cells][[sample_column]]
    valid_samples <- intersect(valid_samples, valid_ctrl_samples)

    .msg("  Valid samples: ", length(valid_samples), verbose = verbose)

    if (length(valid_samples) < 2) {
      .msg("  Need at least 2 valid samples", verbose = verbose)
      next
    }

    target_valid <- target_data[get(sample_column) %in% valid_samples]
    target_barcodes <- target_valid$barcode

    count_matrix_ct <- count_matrix_full[, target_barcodes, drop = FALSE]
    total_counts_target <- colSums(count_matrix_ct)

    sig_genes <- intersect(sig_genes, rownames(count_matrix_ct))
    if (length(sig_genes) == 0) {
      .msg("  No genes available", verbose = verbose)
      next
    }

    ct_output_dir <- file.path(output_dir, "per_celltype", ct_name)
    .ensure_dir(ct_output_dir)

    # Bivariate GLM per gene per sample
    .msg("  Step 1: Fitting bivariate Poisson GLM...", verbose = verbose)

    coef_results <- data.table::rbindlist(lapply(sig_genes, function(g) {
      gene_counts <- as.numeric(count_matrix_ct[g, target_barcodes])

      data.table::rbindlist(lapply(valid_samples, function(samp) {
        samp_idx <- which(target_valid[[sample_column]] == samp)
        if (length(samp_idx) < min_cells_per_sample) {
          return(data.table::data.table(
            gene = g, sample_id = samp,
            coef = NA_real_, se = NA_real_,
            n_cells = length(samp_idx), pval = NA_real_,
            dispersion = NA_real_
          ))
        }

        samp_counts <- gene_counts[samp_idx]
        samp_dist_query <- target_valid[samp_idx]$dist_to_query
        samp_dist_control <- target_valid[samp_idx]$dist_to_control
        samp_total <- total_counts_target[target_barcodes[samp_idx]]

        fit_result <- fit_poisson_controlled(
          samp_counts, samp_dist_query, samp_dist_control, samp_total,
          min_cells = min_expr_cells_glm
        )

        if (is.null(fit_result)) {
          data.table::data.table(
            gene = g, sample_id = samp,
            coef = NA_real_, se = NA_real_,
            n_cells = length(samp_idx), pval = NA_real_,
            dispersion = NA_real_
          )
        } else {
          data.table::data.table(
            gene = g, sample_id = samp,
            coef = fit_result$beta, se = fit_result$se,
            n_cells = fit_result$n_cells, pval = fit_result$pval,
            dispersion = fit_result$dispersion
          )
        }
      }))
    }), fill = TRUE)

    data.table::fwrite(coef_results,
                       file.path(ct_output_dir, "coef_per_sample.csv"))

    # Fisher's combined p-value
    .msg("  Step 2: Combining with Fisher's method...", verbose = verbose)

    fisher_results <- data.table::rbindlist(lapply(sig_genes, function(g) {
      gene_data <- coef_results[gene == g]
      result <- compute_fisher_pval(
        pvals = gene_data$pval,
        coefs = gene_data$coef,
        min_samples = 2L,
        sign_threshold = sign_consistency
      )
      data.table::data.table(
        gene = g,
        stage2_median_coef = result$median_coef,
        stage2_fisher_pval = result$fisher_pval,
        stage2_n_samples = result$n_valid
      )
    }), fill = TRUE)

    fisher_results[, stage2_fisher_fdr := stats::p.adjust(
      stage2_fisher_pval, method = "BH")]

    # Sign consistency for stage 2
    stage2_sign <- data.table::rbindlist(lapply(sig_genes, function(g) {
      gene_data <- coef_results[gene == g]
      v <- gene_data$coef[!is.na(gene_data$coef)]
      n_p <- sum(v > 0)
      n_n <- sum(v < 0)
      n_v <- n_p + n_n
      sc <- if (n_v > 0) max(n_p, n_n) / n_v else NA_real_
      data.table::data.table(gene = g, stage2_sign_consistency = sc)
    }))
    fisher_results <- merge(fisher_results, stage2_sign, by = "gene",
                            all.x = TRUE)

    # Merge with Stage 1
    .msg("  Step 3: Comparing Stage 1 vs Stage 2...", verbose = verbose)

    # Build Stage 1 subset with standardized column names
    stage1_coef_col <- if ("median_coef" %in% names(ct_stage1)) {
      "median_coef"
    } else {
      "combined_coef"
    }
    stage1_fdr_col <- if ("fisher_fdr" %in% names(ct_stage1)) {
      "fisher_fdr"
    } else {
      "fdr"
    }

    cols_needed <- c("gene", stage1_coef_col, stage1_fdr_col)
    cols_needed <- intersect(cols_needed, names(ct_stage1))
    ct_stage1_sub <- ct_stage1[, ..cols_needed]
    if (stage1_coef_col != "stage1_coef") {
      data.table::setnames(ct_stage1_sub, stage1_coef_col, "stage1_coef",
                           skip_absent = TRUE)
    }
    if (stage1_fdr_col != "stage1_fdr") {
      data.table::setnames(ct_stage1_sub, stage1_fdr_col, "stage1_fdr",
                           skip_absent = TRUE)
    }

    comparison <- merge(ct_stage1_sub, fisher_results, by = "gene",
                        all.x = TRUE)
    data.table::setDT(comparison)

    # Classify genes
    attenuation_threshold <- 0.5
    comparison[, coef_ratio := abs(stage2_median_coef) / abs(stage1_coef)]

    comparison[, classification := data.table::fcase(
      is.na(stage2_fisher_fdr), "no_stage2_result",
      stage2_fisher_fdr < fdr_threshold &
        sign(stage2_median_coef) == sign(stage1_coef), query_specific_label,
      stage2_fisher_fdr >= fdr_threshold &
        coef_ratio < attenuation_threshold, "niche_driven",
      stage2_fisher_fdr >= fdr_threshold &
        coef_ratio >= attenuation_threshold, "underpowered",
      stage2_fisher_fdr < fdr_threshold &
        sign(stage2_median_coef) != sign(stage1_coef), "reversed",
      default = "unclassified"
    )]

    # Detect enhanced genes
    comparison[classification == query_specific_label &
                 abs(stage2_median_coef) > abs(stage1_coef) * 1.1,
               classification := "enhanced"]

    comparison[, coef_ratio := NULL]
    comparison[, cell_type := ct_name]
    comparison[, control_celltype := control_celltype]

    data.table::fwrite(comparison,
                       file.path(ct_output_dir, "stage2_comparison.csv"))
    .msg("  Saved: stage2_comparison.csv", verbose = verbose)

    # Classification summary
    class_summary <- comparison[, .N, by = classification]
    data.table::setorder(class_summary, -N)
    .msg("\n  Classification summary:", verbose = verbose)
    for (i in seq_len(nrow(class_summary))) {
      .msg("    ", class_summary$classification[i], ": ",
           class_summary$N[i], verbose = verbose)
    }

    # Visualization: Stage 1 vs Stage 2 scatter
    plot_data <- comparison[!is.na(stage2_median_coef)]

    if (nrow(plot_data) > 0) {
      max_range <- max(abs(c(plot_data$stage1_coef,
                             plot_data$stage2_median_coef)),
                       na.rm = TRUE) * 1.1

      class_colors <- stats::setNames(
        c("#E74C3C", "#8E44AD", "#3498DB", "#F1C40F", "#F39C12", "grey70"),
        c(query_specific_label, "enhanced", "niche_driven", "underpowered",
          "reversed", "no_stage2_result")
      )

      label_genes <- head(plot_data[order(stage1_fdr)], 15)

      p_scatter <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = stage1_coef, y = stage2_median_coef,
                     color = classification)
      ) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                             color = "grey50") +
        ggplot2::geom_hline(yintercept = 0, linetype = "dotted",
                            color = "grey70") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dotted",
                            color = "grey70") +
        ggplot2::geom_point(alpha = 0.7, size = 2) +
        ggplot2::scale_color_manual(values = class_colors,
                                    name = "Classification") +
        ggrepel::geom_text_repel(
          data = label_genes,
          ggplot2::aes(label = gene),
          size = 3, max.overlaps = 15, box.padding = 0.4
        ) +
        ggplot2::xlim(-max_range, max_range) +
        ggplot2::ylim(-max_range, max_range) +
        ggplot2::labs(
          title = sprintf("Stage 1 vs Stage 2 Coefficients: %s", ct_name),
          subtitle = sprintf("Control: %s | diagonal = no change",
                             control_celltype),
          x = "Stage 1 coefficient (log-rate per um)",
          y = "Stage 2 coefficient (adjusted, log-rate per um)"
        ) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold"),
          aspect.ratio = 1
        )

      ggplot2::ggsave(
        file.path(ct_output_dir, "stage1_vs_stage2_scatter.pdf"),
        p_scatter, width = 8, height = 7
      )
      .msg("  Saved: stage1_vs_stage2_scatter.pdf", verbose = verbose)
    }

    all_stage2_results[[ct_name]] <- comparison
  }

  # --------------------------------------------------------------------------
  # 8. Combined summary
  # --------------------------------------------------------------------------
  if (length(all_stage2_results) > 0) {
    combined <- data.table::rbindlist(all_stage2_results, fill = TRUE)
    data.table::fwrite(combined,
                       file.path(output_dir, "summary",
                                 "stage2_all_results.csv"))
    .msg("\nSaved: summary/stage2_all_results.csv", verbose = verbose)

    class_summary <- combined[, .N, by = .(cell_type, classification)]
    class_wide <- data.table::dcast(class_summary,
                                    cell_type ~ classification,
                                    value.var = "N", fill = 0)
    data.table::fwrite(class_wide,
                       file.path(output_dir, "summary",
                                 "classification_summary.csv"))
    .msg("Saved: summary/classification_summary.csv", verbose = verbose)

    # Print summary
    .msg("\n", strrep("=", 70), verbose = verbose)
    .msg("Stage 4 Analysis Summary", verbose = verbose)
    .msg(strrep("=", 70), verbose = verbose)
    for (ct in names(all_stage2_results)) {
      ct_result <- all_stage2_results[[ct]]
      n_total <- nrow(ct_result)
      n_specific <- sum(ct_result$classification %in%
                          c(query_specific_label, "enhanced"), na.rm = TRUE)
      n_niche <- sum(ct_result$classification == "niche_driven",
                     na.rm = TRUE)
      pct_specific <- round(n_specific / n_total * 100, 1)
      .msg(sprintf("  %s: %d genes -> %d %s-specific (%.1f%%), %d niche-driven",
                   ct, n_total, n_specific, query_label, pct_specific,
                   n_niche),
           verbose = verbose)
    }

    return(invisible(combined))
  }

  .msg("No cell types had sufficient data for Stage 4 analysis.",
       verbose = verbose)
  return(invisible(data.table::data.table()))
}


# ============================================================================
# merge_ripple_results()
# ============================================================================

#' Merge RIPPLE results across cell types
#'
#' Combines per-celltype results into summary tables. If
#' \code{coef_per_sample.csv} files are available in the per-celltype
#' directories, recomputes Fisher's combined p-values for consistent
#' methodology across all cell types.
#'
#' @param results_dir Path to results directory containing
#'   \code{per_celltype/} subdirectories with
#'   \code{meta_analysis_results.csv} files.
#' @param fdr_threshold FDR threshold for selecting top genes
#'   (default: \code{0.05}).
#' @param sign_threshold Sign consistency threshold. Fraction of samples
#'   that must agree on coefficient direction for Fisher's p-value to be
#'   computed (default: \code{1.0}).
#' @param recompute_fisher If \code{TRUE} (default), recompute Fisher's
#'   combined p-values from \code{coef_per_sample.csv} files. If
#'   \code{FALSE}, use existing values from
#'   \code{meta_analysis_results.csv}.
#' @param verbose Print progress messages (default: \code{TRUE}).
#'
#' @return A \code{data.table} with all genes across all cell types,
#'   with updated Fisher's combined p-values. Results are also written to
#'   \code{results_dir/summary/}.
#'
#' @details
#' Output files written to \code{results_dir/summary/}:
#' \describe{
#'   \item{all_genes_results.csv}{All genes from all cell types combined.}
#'   \item{top_gradient_genes.csv}{Top 50 genes per cell type by FDR.}
#'   \item{top_gradient_genes_fisher.csv}{Top 50 by Fisher FDR (if
#'     available).}
#'   \item{decay_pattern_summary.csv}{Decay pattern counts per cell type.}
#'   \item{permutation_summary.csv}{Permutation testing summary (if
#'     available).}
#' }
#'
#' @export
merge_ripple_results <- function(
  results_dir,
  fdr_threshold    = 0.05,
  sign_threshold   = 1.0,
  recompute_fisher = TRUE,
  verbose          = TRUE
) {

  fdr_threshold  <- .resolve(fdr_threshold,  "fdr_threshold",    0.05)
  sign_threshold <- .resolve(sign_threshold, "sign_consistency", 1.0)
  verbose        <- .resolve(verbose,        "verbose",          TRUE)

  .msg(strrep("=", 70), verbose = verbose)
  .msg("Merge RIPPLE Results", verbose = verbose)
  .msg(strrep("=", 70), verbose = verbose)
  .msg("Results dir: ", results_dir, verbose = verbose)

  # Find per-celltype directories
  celltype_dir <- file.path(results_dir, "per_celltype")
  if (!dir.exists(celltype_dir)) {
    stop("Per-celltype directory not found: ", celltype_dir, call. = FALSE)
  }

  celltype_dirs <- list.dirs(celltype_dir, recursive = FALSE)
  .msg("\nFound ", length(celltype_dirs), " cell type directories",
       verbose = verbose)

  # Load and combine results
  all_results <- data.table::rbindlist(lapply(celltype_dirs, function(d) {
    f <- file.path(d, "meta_analysis_results.csv")
    if (file.exists(f)) {
      dt <- data.table::fread(f)
      dt[, cell_type := basename(d)]
      .msg("  Loaded ", basename(d), ": ", nrow(dt), " genes",
           verbose = verbose)
      return(dt)
    } else {
      .msg("  WARNING: ", basename(d),
           " - meta_analysis_results.csv not found", verbose = verbose)
      return(NULL)
    }
  }), fill = TRUE)

  if (nrow(all_results) == 0) {
    stop("No results found to merge!", call. = FALSE)
  }

  .msg("\nTotal genes loaded: ", nrow(all_results), verbose = verbose)

  # Recompute Fisher's combined p-values from per-sample data
  if (recompute_fisher) {
    .msg("\nRecomputing Fisher's combined p-values...", verbose = verbose)

    for (d in celltype_dirs) {
      ct_name <- basename(d)
      coef_file <- file.path(d, "coef_per_sample.csv")
      if (!file.exists(coef_file)) {
        .msg("  ", ct_name, ": no coef_per_sample.csv, skipping recompute",
             verbose = verbose)
        next
      }

      coef_data <- data.table::fread(coef_file)
      genes <- unique(coef_data$gene)

      fisher_dt <- data.table::rbindlist(lapply(genes, function(g) {
        gd <- coef_data[gene == g]
        result <- compute_fisher_pval(
          pvals = gd$pval,
          coefs = gd$coef,
          min_samples = 2L,
          sign_threshold = sign_threshold
        )
        data.table::data.table(
          gene = g,
          median_coef = result$median_coef,
          fisher_pval = result$fisher_pval
        )
      }))
      fisher_dt[, fisher_fdr := stats::p.adjust(fisher_pval, method = "BH")]

      # Update all_results for this cell type
      ct_idx <- which(all_results$cell_type == ct_name)
      if (length(ct_idx) > 0) {
        matched <- match(all_results$gene[ct_idx], fisher_dt$gene)
        all_results[ct_idx, median_coef := fisher_dt$median_coef[matched]]
        all_results[ct_idx, fisher_pval := fisher_dt$fisher_pval[matched]]
        all_results[ct_idx, fisher_fdr := fisher_dt$fisher_fdr[matched]]
      }

      .msg("  ", ct_name, ": recomputed Fisher for ", nrow(fisher_dt),
           " genes", verbose = verbose)
    }
  }

  # Write summary files
  summary_dir <- file.path(results_dir, "summary")
  .ensure_dir(summary_dir)

  # All results
  data.table::fwrite(all_results,
                     file.path(summary_dir, "all_genes_results.csv"))
  .msg("\nSaved: summary/all_genes_results.csv", verbose = verbose)

  # Top genes by FDR
  top_genes <- all_results[fdr < fdr_threshold][
    order(fdr), head(.SD, 50), by = cell_type]
  data.table::fwrite(top_genes,
                     file.path(summary_dir, "top_gradient_genes.csv"))
  .msg("Saved: summary/top_gradient_genes.csv (", nrow(top_genes), " genes)",
       verbose = verbose)

  # Top genes by Fisher FDR
  if ("fisher_fdr" %in% names(all_results)) {
    top_fisher <- all_results[fisher_fdr < fdr_threshold][
      order(fisher_fdr), head(.SD, 50), by = cell_type]
    data.table::fwrite(top_fisher,
                       file.path(summary_dir, "top_gradient_genes_fisher.csv"))
    .msg("Saved: summary/top_gradient_genes_fisher.csv (",
         nrow(top_fisher), " genes)", verbose = verbose)
  }

  # Decay pattern summary
  if ("decay_pattern" %in% names(all_results)) {
    decay_summary <- all_results[fdr < fdr_threshold, .N,
                                 by = .(cell_type, decay_pattern)]
    data.table::setorder(decay_summary, cell_type, -N)
    data.table::fwrite(decay_summary,
                       file.path(summary_dir, "decay_pattern_summary.csv"))
    .msg("Saved: summary/decay_pattern_summary.csv", verbose = verbose)
  }

  # Permutation summary
  if ("perm_pval" %in% names(all_results)) {
    perm_summary <- all_results[!is.na(perm_pval), .(
      n_tested = .N,
      n_perm_sig = sum(perm_pval < 0.05),
      pct_perm_sig = round(mean(perm_pval < 0.05) * 100, 1)
    ), by = cell_type]
    data.table::fwrite(perm_summary,
                       file.path(summary_dir, "permutation_summary.csv"))
    .msg("Saved: summary/permutation_summary.csv", verbose = verbose)
  }

  # Print summary report
  .msg("\n", strrep("=", 70), verbose = verbose)
  .msg("Results Summary", verbose = verbose)
  .msg(strrep("=", 70), verbose = verbose)

  has_fisher <- "fisher_fdr" %in% names(all_results)
  summary_stats <- all_results[, .(
    n_genes = .N,
    n_significant = sum(fdr < fdr_threshold, na.rm = TRUE),
    n_fisher_sig = if (has_fisher) sum(fisher_fdr < fdr_threshold,
                                       na.rm = TRUE) else NA_integer_,
    n_perm_tested = sum(!is.na(perm_pval)),
    n_perm_sig = sum(perm_pval < 0.05, na.rm = TRUE)
  ), by = cell_type]
  if (verbose) print(summary_stats)

  .msg("\n", strrep("=", 70), verbose = verbose)
  .msg("Done!", verbose = verbose)
  .msg(strrep("=", 70), verbose = verbose)

  return(invisible(all_results))
}


