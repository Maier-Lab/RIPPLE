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
      call. = FALSE
    )
  }

  species <- switch(tolower(organism),
    "mouse" = "Mus musculus",
    "human" = "Homo sapiens",
    stop("Unsupported organism '", organism, "'. Use 'mouse' or 'human'.",
      call. = FALSE
    )
  )

  # Map shorthand names to msigdbr category/subcategory
  params <- switch(tolower(collection),
    "hallmark" = list(category = "H", subcategory = NULL),
    "kegg" = list(category = "C2", subcategory = "CP:KEGG"),
    "reactome" = list(category = "C2", subcategory = "CP:REACTOME"),
    "go_bp" = list(category = "C5", subcategory = "GO:BP"),
    stop("Unknown gene set collection '", collection,
      "'. Use one of: hallmark, kegg, reactome, go_bp",
      call. = FALSE
    )
  )

  if (is.null(params$subcategory)) {
    gs_df <- msigdbr::msigdbr(species = species, category = params$category)
  } else {
    gs_df <- msigdbr::msigdbr(
      species = species, category = params$category,
      subcategory = params$subcategory
    )
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
#' @param exclude_contamination Logical. If \code{TRUE}, drops genes flagged
#'   as cross-cell-type contamination from the ranked list before fGSEA.
#'   When the input does not already carry a \code{specificity_class}
#'   column, \code{\link{classify_gene_specificity}} is called internally
#'   with \code{contamination_threshold} to compute it. Default \code{FALSE}.
#'   See "Limitations" below before turning this on.
#' @param contamination_threshold Integer. Minimum number of cell types in
#'   which a gene must be significant (at
#'   \code{contamination_sig_threshold}) to be classified as
#'   contamination. Used only when \code{exclude_contamination = TRUE}
#'   and \code{specificity_class} is not already present. Default:
#'   \code{4L} -- illustrative; tune to your panel size and annotation
#'   granularity (see "Choosing \code{contamination_threshold}" in
#'   \code{\link{classify_gene_specificity}}).
#' @param contamination_sig_threshold Numeric or character. Stricter
#'   significance bar for the contamination tally. \code{NULL} (default)
#'   uses \code{fdr_col} < 0.05 (any "*" hit counts). Pass \code{"**"}
#'   (FDR < 0.01), \code{"***"} (FDR < 0.001), or a numeric to require a
#'   higher significance bar -- this makes the filter \strong{looser}
#'   (fewer genes flagged). See \code{\link{classify_gene_specificity}}
#'   for the full discussion.
#' @param exclude_specificity_class Optional character vector of
#'   specificity-class labels (e.g. \code{"contamination"}) to drop from the
#'   ranked list. Lower-level alternative to \code{exclude_contamination}
#'   when you want to drop other classes (e.g. \code{c("contamination",
#'   "ubiquitous")}). Requires the input to carry a \code{specificity_class}
#'   column. Default \code{NULL}.
#'
#' @return A \code{data.table} with columns: cell_type, pathway, pathway_clean,
#'   pval, padj, ES, NES, size, leadingEdge. The leadingEdge column contains
#'   comma-separated gene names.
#'
#' @section Limitations of contamination filtering:
#' Setting \code{exclude_contamination = TRUE} (or
#' \code{exclude_specificity_class = "contamination"}) is a useful sanity
#' check, but it is a \strong{heuristic} that can both over- and
#' under-correct. Be aware of the following before reporting filtered
#' enrichments:
#' \itemize{
#'   \item \strong{Loses real biology when contamination + induction co-occur.}
#'     A gene can be both genuinely induced near the query \emph{and}
#'     contaminated by ambient RNA from query cells (e.g. MIF in the CosMx
#'     NSCLC walkthrough). Filtering removes the gene from every cell type's
#'     ranking, including cell types where the induction is real.
#'   \item \strong{The threshold is panel-dependent.} The default
#'     \code{contamination_threshold = 4} assumes ~10-20 cell types in the
#'     panel. With 3 cell types the cutoff is unreachable; with 30 fine
#'     subtypes it is too lenient. Re-tune for your annotation granularity.
#'   \item \strong{Removes genes globally, not per cell type.} A gene
#'     classified as contamination in the dataset as a whole is dropped from
#'     every cell type's ranking, even where its expression is genuinely
#'     specific.
#'   \item \strong{Multi-cell-type biology is not contamination.} Many
#'     cytokines, MHC genes, and stress-response genes are legitimately
#'     expressed across many cell types. Filtering may remove signal that
#'     reflects real shared biology.
#'   \item \strong{Rank-based tests are panel-size sensitive.} Dropping
#'     genes shifts the ranks of all remaining genes; NES values from
#'     filtered and unfiltered runs are not directly comparable. Always
#'     report whether the filter was applied and ideally show the unfiltered
#'     result alongside.
#' }
#' Recommended use: run fGSEA both with and without the filter, compare
#' top pathways, and treat large discrepancies as a flag for further
#' investigation rather than evidence that one is "right".
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Loads gene sets from msigdbr if \code{gene_sets} is a string,
#'     or uses the provided named list directly.
#'   \item Optionally classifies and drops contamination-class genes (see
#'     \code{exclude_contamination}).
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
                             seed = 42,
                             exclude_contamination = FALSE,
                             contamination_threshold = 4L,
                             contamination_sig_threshold = NULL,
                             exclude_specificity_class = NULL) {
  # --- Validate inputs ---
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required for pathway enrichment.\n",
      "Install with: BiocManager::install('fgsea')",
      call. = FALSE
    )
  }

  if (!inherits(results, "data.table")) {
    results <- data.table::as.data.table(results)
  }

  required_cols <- c("gene", "cell_type", coef_col)
  missing <- setdiff(required_cols, names(results))
  if (length(missing) > 0) {
    stop("Missing required columns in results: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
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
      call. = FALSE
    )
  }
  message("  Gene sets loaded: ", length(pathways), " pathways")

  # --- Optional FDR filter ---
  dt <- data.table::copy(results)
  if (fdr_col %in% names(dt) && fdr_threshold < 1) {
    dt <- dt[!is.na(get(fdr_col)) & get(fdr_col) < fdr_threshold]
    message("  Filtered to FDR < ", fdr_threshold, ": ", nrow(dt), " rows")
  }

  # --- Optional contamination filter (high-level convenience) ---
  if (isTRUE(exclude_contamination)) {
    if (!"specificity_class" %in% names(dt)) {
      spec_dt <- classify_gene_specificity(
        dt,
        fdr_col                     = fdr_col,
        fdr_threshold               = 0.05,
        contamination_threshold     = contamination_threshold,
        contamination_sig_threshold = contamination_sig_threshold
      )
      contam_genes <- spec_dt[
        specificity_class == "contamination", unique(gene)
      ]
      n_before <- nrow(dt)
      dt <- dt[!gene %in% contam_genes]
      sig_label <- if (is.null(contamination_sig_threshold)) {
        "FDR < 0.05"
      } else if (is.character(contamination_sig_threshold)) {
        paste0("significance ", contamination_sig_threshold)
      } else {
        paste0("FDR < ", contamination_sig_threshold)
      }
      message(
        "  exclude_contamination: dropped ", length(contam_genes),
        " gene(s) flagged as contamination ",
        "(>= ", contamination_threshold, " cell types at ",
        sig_label, "); ",
        n_before - nrow(dt), " row(s) removed; ",
        nrow(dt), " row(s) remain"
      )
    } else {
      n_before <- nrow(dt)
      dt <- dt[specificity_class != "contamination" |
                 is.na(specificity_class)]
      message(
        "  exclude_contamination: dropped ",
        n_before - nrow(dt),
        " row(s) with specificity_class == 'contamination'; ",
        nrow(dt), " row(s) remain"
      )
    }
  }

  # --- Optional specificity-class exclusion (lower-level) ---
  if (!is.null(exclude_specificity_class)) {
    if (!"specificity_class" %in% names(dt)) {
      stop(
        "exclude_specificity_class is set but `specificity_class` is not a ",
        "column on `results`. Run classify_gene_specificity() or ",
        "run_ripple_atlas() first to add this column.",
        call. = FALSE
      )
    }
    n_before <- nrow(dt)
    dt <- dt[!specificity_class %in% exclude_specificity_class]
    message(
      "  Dropped ", n_before - nrow(dt),
      " row(s) with specificity_class in {",
      paste(exclude_specificity_class, collapse = ", "),
      "}; ", nrow(dt), " rows remain"
    )
  }

  # --- Run fGSEA per cell type ---
  cell_types <- unique(dt$cell_type)
  message("  Running fGSEA across ", length(cell_types), " cell types...")

  all_fgsea <- data.table::rbindlist(lapply(cell_types, function(ct) {
    ct_data <- dt[cell_type == ct & !is.na(get(coef_col))]
    stats <- stats::setNames(ct_data[[coef_col]], ct_data$gene)
    stats <- sort(stats)

    if (length(stats) < min_genes) {
      message(
        "    ", ct, ": too few genes (", length(stats), " < ",
        min_genes, "), skipping"
      )
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
  out_cols <- c(
    "cell_type", "pathway", "pathway_clean", "pval", "padj",
    "ES", "NES", "size", "leadingEdge"
  )
  out_cols <- intersect(out_cols, names(all_fgsea))
  all_fgsea <- all_fgsea[, ..out_cols]

  data.table::setorder(all_fgsea, cell_type, pval)

  n_sig <- sum(all_fgsea$padj < 0.05, na.rm = TRUE)
  message(
    "  fGSEA complete: ", nrow(all_fgsea), " pathway-celltype tests, ",
    n_sig, " significant (padj < 0.05)"
  )

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
#' @param fdr_threshold Numeric. Loose significance cutoff used to decide
#'   what counts as "significant in a cell type" for the
#'   specific/moderate/ubiquitous classes and for the
#'   \code{n_celltypes} count. Default: 0.05.
#' @param contamination_threshold Integer. Genes significant in at least this
#'   many cell types (at \code{contamination_sig_threshold}) are flagged
#'   as "contamination". Default: \code{4}. \strong{This default is
#'   illustrative, not universal -- you should pick a value appropriate
#'   to your dataset.} See "Choosing \code{contamination_threshold}"
#'   below.
#' @param contamination_sig_threshold Numeric or character. Stricter
#'   significance cutoff used \emph{only} for the contamination flag. Lets
#'   you require, e.g., \code{**} significance (FDR < 0.01) before a gene
#'   counts toward the contamination tally, while still allowing \code{*}
#'   (FDR < 0.05) for the broader specific/moderate/ubiquitous classes.
#'   Accepts a numeric FDR (e.g. \code{0.01}) or a star string:
#'   \code{"*"} -> 0.05, \code{"**"} -> 0.01, \code{"***"} -> 0.001,
#'   \code{"****"} -> 1e-4. Default \code{NULL} -> use \code{fdr_threshold}
#'   (back-compat: any "*"-significant gene contributes). Tightening this
#'   makes the contamination filter \strong{looser} (fewer genes flagged)
#'   because each cell-type "hit" must clear a higher bar.
#'
#' @return A \code{data.table} with columns:
#' \describe{
#'   \item{gene}{Character. Gene name.}
#'   \item{n_celltypes}{Integer. Number of cell types where the gene is
#'     significant at \code{fdr_threshold} (the loose cutoff).}
#'   \item{n_celltypes_strict}{Integer. Number of cell types where the gene
#'     is significant at \code{contamination_sig_threshold}. Equal to
#'     \code{n_celltypes} when the two thresholds match (the default).}
#'   \item{celltypes}{Character. Comma-separated list of cell type names
#'     (using the loose \code{fdr_threshold}).}
#'   \item{specificity_class}{Character. One of: "specific" (1 cell type),
#'     "moderate" (2-3 cell types), "ubiquitous" (4+ but below contamination
#'     threshold), or "contamination" (\code{n_celltypes_strict >=
#'     contamination_threshold}).}
#' }
#'
#' @section Choosing \code{contamination_threshold}:
#' There is no universally correct cutoff. The right value depends on
#' three things:
#' \itemize{
#'   \item \strong{How many cell types are in your panel.} The default of
#'     4 is sensible for ~10-20 cell types (the typical scale of a
#'     mid-resolution immune / stromal annotation, including the original
#'     HyMy / TDLN analyses). With ~3-5 cell types the cutoff is barely
#'     reachable -- everything will be classified as "specific" or
#'     "moderate". With ~30+ fine subtypes the cutoff is too lenient --
#'     many genes that are biologically expressed across a few related
#'     subtypes will pass without being flagged.
#'   \item \strong{Annotation granularity.} Coarse types (e.g. "T cell",
#'     "myeloid") share fewer marker genes than fine subtypes (e.g.
#'     "CD4 memory", "CD8 effector", "Treg"). Fine annotations need a
#'     larger cutoff to avoid flagging real lineage-shared biology as
#'     contamination.
#'   \item \strong{Biological prior.} Cytokines, MHC genes, ribosomal
#'     genes, and stress-response programs are legitimately expressed
#'     across many cell types -- the cross-cell-type heuristic cannot
#'     tell them apart from ambient RNA on its own. Inspect the flagged
#'     genes and decide whether the cutoff is doing what you want.
#' }
#' \strong{Recommended workflow:} run with the default first, look at
#' the flagged gene list (\code{spec[specificity_class ==
#' "contamination"]}), and adjust up or down. A useful sanity check is
#' to plot \code{n_celltypes} as a histogram across all significant
#' genes -- the cutoff should sit in the right tail, separating clear
#' multi-cell-type genes from cell-type-specific signal. Whatever value
#' you choose, report it explicitly in the methods.
#'
#' @section Tuning \code{contamination_sig_threshold} (the significance bar):
#' The contamination flag has \emph{two} dials:
#' \itemize{
#'   \item \code{contamination_threshold} -- how many cell types
#'     (default \code{4}).
#'   \item \code{contamination_sig_threshold} -- how strictly significant
#'     each "hit" must be (default \code{NULL} -> uses
#'     \code{fdr_threshold}, i.e. any \code{*} gene counts).
#' }
#' Tightening the significance bar (e.g. \code{"**"} for FDR < 0.01)
#' makes the filter \strong{looser overall} -- a gene only counts as a
#' "hit" in cell types where the signal is unambiguous, so fewer genes
#' reach the cell-type tally and fewer are flagged. Use this when the
#' default flags genes you believe are real biology, or when you want
#' the contamination class to capture only the most blatant ambient-RNA
#' offenders. Conversely, leave it at the default (or use \code{NULL})
#' to be more aggressive about flagging.
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
                                      contamination_threshold = 4,
                                      contamination_sig_threshold = NULL) {
  if (!inherits(results, "data.table")) {
    results <- data.table::as.data.table(results)
  }

  required_cols <- c("gene", "cell_type", fdr_col)
  missing <- setdiff(required_cols, names(results))
  if (length(missing) > 0) {
    stop("Missing required columns in results: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  # Resolve contamination_sig_threshold:
  #   NULL -> use fdr_threshold (back-compat)
  #   "*" / "**" / "***" / "****" -> 0.05 / 0.01 / 0.001 / 1e-4
  #   numeric -> as-is
  if (is.null(contamination_sig_threshold)) {
    cs_thr <- fdr_threshold
  } else if (is.character(contamination_sig_threshold)) {
    star_map <- c(`*` = 0.05, `**` = 0.01, `***` = 0.001, `****` = 1e-4)
    if (!contamination_sig_threshold %in% names(star_map)) {
      stop(
        "contamination_sig_threshold must be NULL, a numeric FDR, or one ",
        "of '*', '**', '***', '****'.", call. = FALSE
      )
    }
    cs_thr <- unname(star_map[contamination_sig_threshold])
  } else if (is.numeric(contamination_sig_threshold) &&
             length(contamination_sig_threshold) == 1) {
    cs_thr <- contamination_sig_threshold
  } else {
    stop(
      "contamination_sig_threshold must be NULL, a single numeric FDR, ",
      "or one of '*', '**', '***', '****'.", call. = FALSE
    )
  }
  if (cs_thr > fdr_threshold) {
    warning(
      "contamination_sig_threshold (", cs_thr, ") is looser than ",
      "fdr_threshold (", fdr_threshold, "); contamination flag will be ",
      "MORE aggressive than the overall significance gate.",
      call. = FALSE
    )
  }

  # Loose set: significant at fdr_threshold (drives n_celltypes,
  # specific/moderate/ubiquitous classes, and the celltypes string).
  sig <- results[!is.na(get(fdr_col)) & get(fdr_col) < fdr_threshold]

  if (nrow(sig) == 0) {
    message("No significant genes at FDR < ", fdr_threshold)
    return(data.table::data.table(
      gene = character(), n_celltypes = integer(),
      n_celltypes_strict = integer(),
      celltypes = character(), specificity_class = character()
    ))
  }

  gene_counts <- sig[, .(
    n_celltypes = data.table::uniqueN(cell_type),
    celltypes   = paste(sort(unique(as.character(cell_type))), collapse = ", ")
  ), by = gene]

  # Strict set: significant at contamination_sig_threshold (drives the
  # contamination flag only).
  strict <- results[!is.na(get(fdr_col)) & get(fdr_col) < cs_thr,
                    .(n_celltypes_strict = data.table::uniqueN(cell_type)),
                    by = gene]
  gene_counts <- merge(gene_counts, strict, by = "gene", all.x = TRUE)
  gene_counts[is.na(n_celltypes_strict), n_celltypes_strict := 0L]

  # Classify: contamination first (uses strict count), then
  # specific/moderate/ubiquitous (uses loose count).
  gene_counts[, specificity_class := data.table::fifelse(
    n_celltypes_strict >= contamination_threshold, "contamination",
    data.table::fifelse(
      n_celltypes == 1, "specific",
      data.table::fifelse(
        n_celltypes <= 3, "moderate",
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
#'   geom_errorbar(aes(
#'     ymin = prop_expressing - se_prop,
#'     ymax = prop_expressing + se_prop
#'   ), width = 2) +
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

    if (n == 0) {
      return(NULL)
    }

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
