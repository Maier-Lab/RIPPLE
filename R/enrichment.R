#' @title Pathway Enrichment and Gene Specificity Functions
#'
#' @description Functions for running pathway enrichment analysis (fGSEA),
#'   classifying gene specificity across cell types, and binning expression
#'   data by distance for decay curve visualization.
#'
#' @name enrichment
NULL


# ============================================================================
# Gene Set Loading Helper
# ============================================================================

#' Load gene sets from msigdbr
#'
#' Retrieves gene set collections using msigdbr, mapping collection shorthand
#' names to the appropriate category/subcategory parameters.
#'
#' @param collection Character: "hallmark", "kegg", "reactome", "go_bp"
#' @param organism Character: "mouse" or "human"
#'
#' @return Named list of character vectors (pathway name -> gene symbols)
#'
#' @noRd
.load_msigdbr_sets <- function(collection, organism = "mouse") {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Package 'msigdbr' is required for gene set loading.\n",
         "Install with: BiocManager::install('msigdbr')",
         call. = FALSE)
  }

  species <- switch(
    tolower(organism),
    "mouse" = "Mus musculus",
    "human" = "Homo sapiens",
    stop("Unsupported organism '", organism, "'. Use 'mouse' or 'human'.",
         call. = FALSE)
  )

  # Map shorthand names to msigdbr category/subcategory
  params <- switch(
    tolower(collection),
    "hallmark" = list(category = "H", subcategory = NULL),
    "kegg"     = list(category = "C2", subcategory = "CP:KEGG"),
    "reactome" = list(category = "C2", subcategory = "CP:REACTOME"),
    "go_bp"    = list(category = "C5", subcategory = "GO:BP"),
    stop("Unknown gene set collection '", collection,
         "'. Use one of: hallmark, kegg, reactome, go_bp",
         call. = FALSE)
  )

  if (is.null(params$subcategory)) {
    gs_df <- msigdbr::msigdbr(species = species, category = params$category)
  } else {
    gs_df <- msigdbr::msigdbr(species = species, category = params$category,
                               subcategory = params$subcategory)
  }

  split(gs_df$gene_symbol, gs_df$gs_name)
}


#' Clean pathway names for display
#'
#' Strips common prefixes (e.g., HALLMARK_, KEGG_, REACTOME_) and converts
#' underscores to spaces with title case.
#'
#' @param x Character vector of pathway names
#'
#' @return Character vector of cleaned names
#'
#' @noRd
.clean_pathway_name <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  x <- sub("^KEGG_", "", x)
  x <- sub("^REACTOME_", "", x)
  x <- sub("^GOBP_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}


# ============================================================================
# fGSEA Pathway Enrichment
# ============================================================================

#' Run fGSEA pathway enrichment on RIPPLE gradient results
#'
#' Ranks genes by gradient coefficient and runs fast Gene Set Enrichment
#' Analysis per cell type using the fgsea package.
#'
#' For RIPPLE results, negative coefficients indicate genes induced near the
#' query cell type, so negative NES values correspond to pathways enriched
#' among query-induced genes.
#'
#' @param results \code{data.table} with columns: gene, cell_type, and the
#'   columns specified by \code{coef_col} and \code{fdr_col}.
#' @param gene_sets Character string ("hallmark", "kegg", "reactome", "go_bp")
#'   or a named list of character vectors (pathway name -> gene symbols).
#' @param organism Character: "mouse" or "human" (used when \code{gene_sets}
#'   is a string, for msigdbr species mapping). Default: "mouse".
#' @param coef_col Character. Column name for the ranking statistic.
#'   Default: "median_coef".
#' @param fdr_col Character. Column name for the significance filter.
#'   Default: "fisher_fdr".
#' @param fdr_threshold Numeric. Only include genes with FDR below this value
#'   in the ranking. Default: 1 (i.e., all genes included).
#' @param min_size Integer. Minimum gene set size for fGSEA. Default: 15.
#' @param max_size Integer. Maximum gene set size for fGSEA. Default: 500.
#' @param n_perm Integer. Number of fGSEA permutations. Default: 10000.
#' @param min_genes Integer. Minimum number of ranked genes required to run
#'   fGSEA for a cell type. Cell types with fewer genes are skipped.
#'   Default: 100. Lower for small custom gene panels.
#' @param seed Integer. Random seed for fGSEA permutations. Set to ensure
#'   reproducibility. Default: 42.
#'
#' @return A \code{data.table} with columns: cell_type, pathway, pathway_clean,
#'   pval, padj, ES, NES, size, leadingEdge. The leadingEdge column contains
#'   comma-separated gene names.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Loads gene sets from msigdbr if \code{gene_sets} is a string,
#'     or uses the provided named list directly.
#'   \item For each cell type, creates a ranked gene list sorted by
#'     \code{coef_col} values (named by gene).
#'   \item Runs \code{fgsea::fgsea()} with the specified parameters.
#'   \item Combines results across cell types.
#' }
#'
#' Cell types with fewer than \code{min_genes} ranked genes are skipped.
#'
#' @examples
#' \dontrun{
#' results <- data.table::fread("all_genes_results.csv")
#' gsea_res <- run_ripple_fgsea(results, gene_sets = "hallmark", organism = "mouse")
#' # Top enriched pathways
#' gsea_res[padj < 0.05][order(pval)]
#' }
#'
#' @importFrom data.table data.table rbindlist setorder copy
#' @export
run_ripple_fgsea <- function(results,
                              gene_sets = "hallmark",
                              organism = "mouse",
                              coef_col = "median_coef",
                              fdr_col = "fisher_fdr",
                              fdr_threshold = 1,
                              min_size = 15,
                              max_size = 500,
                              n_perm = 10000,
                              min_genes = 100,
                              seed = 42) {

  # --- Validate inputs ---
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required for pathway enrichment.\n",
         "Install with: BiocManager::install('fgsea')",
         call. = FALSE)
  }

  if (!inherits(results, "data.table")) {
    results <- data.table::as.data.table(results)
  }

  required_cols <- c("gene", "cell_type", coef_col)
  missing <- setdiff(required_cols, names(results))
  if (length(missing) > 0) {
    stop("Missing required columns in results: ",
         paste(missing, collapse = ", "),
         call. = FALSE)
  }

  # --- Load gene sets ---
  if (is.character(gene_sets) && length(gene_sets) == 1) {
    message("Loading ", gene_sets, " gene sets for ", organism, "...")
    pathways <- .load_msigdbr_sets(gene_sets, organism)
  } else if (is.list(gene_sets) && !is.null(names(gene_sets))) {
    pathways <- gene_sets
  } else {
    stop("gene_sets must be a character string (e.g., 'hallmark') or ",
         "a named list of gene vectors.",
         call. = FALSE)
  }
  message("  Gene sets loaded: ", length(pathways), " pathways")

  # --- Optional FDR filter ---
  dt <- data.table::copy(results)
  if (fdr_col %in% names(dt) && fdr_threshold < 1) {
    dt <- dt[!is.na(get(fdr_col)) & get(fdr_col) < fdr_threshold]
    message("  Filtered to FDR < ", fdr_threshold, ": ", nrow(dt), " rows")
  }

  # --- Run fGSEA per cell type ---
  cell_types <- unique(dt$cell_type)
  message("  Running fGSEA across ", length(cell_types), " cell types...")

  all_fgsea <- data.table::rbindlist(lapply(cell_types, function(ct) {
    ct_data <- dt[cell_type == ct & !is.na(get(coef_col))]
    stats <- stats::setNames(ct_data[[coef_col]], ct_data$gene)
    stats <- sort(stats)

    if (length(stats) < min_genes) {
      message("    ", ct, ": too few genes (", length(stats), " < ",
              min_genes, "), skipping")
      return(NULL)
    }

    # Seed immediately before fGSEA to ensure reproducibility across runs
    # (fgsea uses a multilevel random permutation internally)
    set.seed(seed)
    res <- fgsea::fgsea(
      pathways = pathways,
      stats = stats,
      minSize = min_size,
      maxSize = max_size,
      nPermSimple = n_perm
    )

    res <- data.table::as.data.table(res)
    res[, cell_type := ct]
    res[, pathway_clean := .clean_pathway_name(pathway)]

    # Convert leadingEdge list column to comma-separated string
    res[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ",")]

    res
  }), fill = TRUE)

  if (nrow(all_fgsea) == 0) {
    message("  No fGSEA results produced.")
    return(data.table::data.table(
      cell_type = character(), pathway = character(),
      pathway_clean = character(), pval = numeric(),
      padj = numeric(), ES = numeric(), NES = numeric(),
      size = integer(), leadingEdge = character()
    ))
  }

  # --- Select and order output columns ---
  out_cols <- c("cell_type", "pathway", "pathway_clean", "pval", "padj",
                "ES", "NES", "size", "leadingEdge")
  out_cols <- intersect(out_cols, names(all_fgsea))
  all_fgsea <- all_fgsea[, ..out_cols]

  data.table::setorder(all_fgsea, cell_type, pval)

  n_sig <- sum(all_fgsea$padj < 0.05, na.rm = TRUE)
  message("  fGSEA complete: ", nrow(all_fgsea), " pathway-celltype tests, ",
          n_sig, " significant (padj < 0.05)")

  all_fgsea
}


# ============================================================================
# Gene Specificity Classification
# ============================================================================

#' Classify gene specificity across cell types
#'
#' For each significant gene, counts how many cell types it is significant in
#' and classifies it as specific, moderate, ubiquitous, or contamination.
#'
#' Genes significant in many cell types are more likely to be segmentation
#' artifacts (query transcripts leaking into neighboring cells) than genuine
#' paracrine effects. This function provides a simple classification to flag
#' such potential artifacts.
#'
#' @param results \code{data.table} with columns: gene, cell_type, and a
#'   significance column specified by \code{fdr_col}.
#' @param fdr_col Character. Column name for the significance measure.
#'   Default: "fisher_fdr".
#' @param fdr_threshold Numeric. Significance cutoff. Default: 0.05.
#' @param contamination_threshold Integer. Genes significant in at least this
#'   many cell types are flagged as "contamination". Default: 4.
#'
#' @return A \code{data.table} with columns:
#' \describe{
#'   \item{gene}{Character. Gene name.}
#'   \item{n_celltypes}{Integer. Number of cell types where the gene is
#'     significant.}
#'   \item{celltypes}{Character. Comma-separated list of cell type names.}
#'   \item{specificity_class}{Character. One of: "specific" (1 cell type),
#'     "moderate" (2-3 cell types), "ubiquitous" (4+ but below contamination
#'     threshold), or "contamination" (>= contamination_threshold).}
#' }
#'
#' @examples
#' \dontrun{
#' results <- data.table::fread("all_genes_results.csv")
#' spec <- classify_gene_specificity(results, fdr_threshold = 0.05)
#' table(spec$specificity_class)
#' # Flag contamination genes for downstream filtering
#' contam_genes <- spec[specificity_class == "contamination"]$gene
#' }
#'
#' @importFrom data.table data.table uniqueN fifelse
#' @export
classify_gene_specificity <- function(results,
                                       fdr_col = "fisher_fdr",
                                       fdr_threshold = 0.05,
                                       contamination_threshold = 4) {

  if (!inherits(results, "data.table")) {
    results <- data.table::as.data.table(results)
  }

  required_cols <- c("gene", "cell_type", fdr_col)
  missing <- setdiff(required_cols, names(results))
  if (length(missing) > 0) {
    stop("Missing required columns in results: ",
         paste(missing, collapse = ", "),
         call. = FALSE)
  }

  # Filter to significant genes
  sig <- results[!is.na(get(fdr_col)) & get(fdr_col) < fdr_threshold]

  if (nrow(sig) == 0) {
    message("No significant genes at FDR < ", fdr_threshold)
    return(data.table::data.table(
      gene = character(), n_celltypes = integer(),
      celltypes = character(), specificity_class = character()
    ))
  }

  # Count cell types per gene
  gene_counts <- sig[, .(
    n_celltypes = data.table::uniqueN(cell_type),
    celltypes = paste(sort(unique(as.character(cell_type))), collapse = ", ")
  ), by = gene]

  # Classify
  gene_counts[, specificity_class := data.table::fifelse(
    n_celltypes == 1, "specific",
    data.table::fifelse(
      n_celltypes <= 3, "moderate",
      data.table::fifelse(
        n_celltypes >= contamination_threshold, "contamination",
        "ubiquitous"
      )
    )
  )]

  gene_counts[]
}


# ============================================================================
# Distance Binning for Decay Curves
# ============================================================================

#' Bin expression data by distance for decay curves
#'
#' Groups cells into equal-width distance bins and computes per-bin summary
#' statistics (proportion expressing, mean expression, standard errors) for
#' visualization of expression decay curves away from query cells.
#'
#' @param counts Numeric vector of raw counts for one gene (one value per cell).
#' @param distances Numeric vector of distances to the nearest query cell
#'   (same length as \code{counts}).
#' @param n_bins Integer. Number of distance bins. Default: 20.
#' @param max_distance Numeric. Maximum distance to include. Cells beyond this
#'   distance are excluded. Default: 200.
#'
#' @return A \code{data.table} with columns:
#' \describe{
#'   \item{bin_center}{Numeric. Center of the distance bin.}
#'   \item{bin_start}{Numeric. Left edge of the bin.}
#'   \item{bin_end}{Numeric. Right edge of the bin.}
#'   \item{n_cells}{Integer. Number of cells in the bin.}
#'   \item{n_expressing}{Integer. Number of cells with count > 0.}
#'   \item{prop_expressing}{Numeric. Proportion of cells expressing (0-1).}
#'   \item{se_prop}{Numeric. Standard error of the proportion.}
#'   \item{mean_expr}{Numeric. Mean expression value in the bin.}
#'   \item{se_expr}{Numeric. Standard error of the mean expression.}
#' }
#'
#' @details
#' Standard error of the proportion is computed as
#' \code{sqrt(p * (1-p) / n)} (binomial SE). Standard error of the mean
#' expression uses \code{sd / sqrt(n)}.
#'
#' Bins with zero cells are excluded from the output.
#'
#' @examples
#' \dontrun{
#' # For a single gene in a single sample
#' binned <- bin_decay_data(
#'   counts = gene_counts_vector,
#'   distances = dist_to_query_vector,
#'   n_bins = 20,
#'   max_distance = 200
#' )
#' # Plot decay curve
#' library(ggplot2)
#' ggplot(binned, aes(x = bin_center, y = prop_expressing)) +
#'   geom_point() +
#'   geom_errorbar(aes(ymin = prop_expressing - se_prop,
#'                     ymax = prop_expressing + se_prop), width = 2) +
#'   labs(x = "Distance to query", y = "Proportion expressing")
#' }
#'
#' @importFrom data.table data.table
#' @export
bin_decay_data <- function(counts, distances, n_bins = 20, max_distance = 200) {

  if (length(counts) != length(distances)) {
    stop("counts and distances must have the same length.", call. = FALSE)
  }

  # Filter to valid distances
  valid <- !is.na(distances) & !is.na(counts) &
    distances >= 0 & distances <= max_distance
  counts <- counts[valid]
  distances <- distances[valid]

  if (length(counts) == 0) {
    return(data.table::data.table(
      bin_center = numeric(), bin_start = numeric(), bin_end = numeric(),
      n_cells = integer(), n_expressing = integer(),
      prop_expressing = numeric(), se_prop = numeric(),
      mean_expr = numeric(), se_expr = numeric()
    ))
  }

  # Create equal-width bins
  bin_width <- max_distance / n_bins
  bin_idx <- pmin(floor(distances / bin_width) + 1L, n_bins)

  # Compute per-bin statistics
  results <- data.table::rbindlist(lapply(seq_len(n_bins), function(i) {
    mask <- bin_idx == i
    n <- sum(mask)

    if (n == 0) return(NULL)

    bin_counts <- counts[mask]
    n_expr <- sum(bin_counts > 0)
    prop <- n_expr / n
    mean_val <- mean(bin_counts)
    sd_val <- if (n > 1) stats::sd(bin_counts) else 0

    data.table::data.table(
      bin_center = (i - 0.5) * bin_width,
      bin_start = (i - 1) * bin_width,
      bin_end = i * bin_width,
      n_cells = n,
      n_expressing = n_expr,
      prop_expressing = prop,
      se_prop = sqrt(prop * (1 - prop) / n),
      mean_expr = mean_val,
      se_expr = sd_val / sqrt(n)
    )
  }))

  results[]
}
