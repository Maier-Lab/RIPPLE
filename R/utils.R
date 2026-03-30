#' @title General Utility Functions
#'
#' @description Miscellaneous helper functions for statistical calculations,
#'   gene scoring, neighborhood analysis, and enrichment analysis.
#'
#' @name utils
NULL

#' Shannon entropy
#'
#' Computes Shannon entropy on the log2 scale for a vector of counts.
#'
#' @param counts Numeric vector of counts (non-negative).
#'
#' @return Numeric. Shannon entropy in bits (log2 scale). Returns 0 if all
#'   counts are zero.
#'
#' @details Computes \code{-sum(p * log2(p))} where \code{p = counts / sum(counts)}.
#'   Zero-count categories are excluded from the summation.
#'
#' @examples
#' shannon_entropy(c(10, 10, 10))  # Uniform: log2(3) = 1.585
#' shannon_entropy(c(100, 0, 0))   # Pure: 0
#' shannon_entropy(c(0, 0, 0))     # Empty: 0
#'
#' @export
shannon_entropy <- function(counts) {
  if (sum(counts) == 0) return(0)

  probs <- counts / sum(counts)
  probs <- probs[probs > 0]  # Remove zeros for log

  -sum(probs * log2(probs))
}


#' Calculate enrichment scores
#'
#' Computes enrichment statistics by comparing observed counts to expected
#' counts across categories.
#'
#' @param observed Named numeric vector of observed counts.
#' @param expected Named numeric vector of expected counts. Names must
#'   overlap with \code{observed}.
#'
#' @return A \code{data.table} with columns:
#' \describe{
#'   \item{\code{category}}{Character. Category names.}
#'   \item{\code{observed}}{Numeric. Observed counts.}
#'   \item{\code{expected}}{Numeric. Expected counts.}
#'   \item{\code{enrichment}}{Numeric. Observed / expected ratio.}
#'   \item{\code{log2_enrichment}}{Numeric. Log2 of the enrichment ratio.}
#' }
#'
#' @details Only categories present in both \code{observed} and \code{expected}
#'   are included. Categories with zero expected counts get \code{NA} for
#'   enrichment and log2_enrichment.
#'
#' @examples
#' \dontrun{
#' obs <- c(A = 50, B = 10, C = 40)
#' exp <- c(A = 33, B = 33, C = 33)
#' calculate_enrichment(obs, exp)
#' }
#'
#' @importFrom data.table data.table
#' @export
calculate_enrichment <- function(observed, expected) {
  # Ensure same order
  common_names <- intersect(names(observed), names(expected))
  observed <- observed[common_names]
  expected <- expected[common_names]

  # Avoid division by zero
  expected_safe <- ifelse(expected == 0, NA_real_, expected)

  data.table::data.table(
    category = common_names,
    observed = observed,
    expected = expected,
    enrichment = observed / expected_safe,
    log2_enrichment = log2(observed / expected_safe)
  )
}


#' Score cells for a gene signature
#'
#' Scores each cell in a Seurat object for a gene signature using
#' \code{Seurat::AddModuleScore}.
#'
#' @param obj A Seurat object.
#' @param genes Character vector of gene names for the signature.
#' @param name Character. Name for the score column in metadata.
#' @param ctrl Integer. Number of control genes per bin (default: 100).
#'
#' @return The Seurat object with the score added to metadata under the
#'   column name specified by \code{name}.
#'
#' @details Genes not found in the Seurat object are silently skipped with
#'   a message. If fewer than 2 genes are available, the score is set to
#'   \code{NA} for all cells. The function corrects for the "1" suffix that
#'   \code{AddModuleScore} appends to the name.
#'
#' @examples
#' \dontrun{
#' obj <- score_gene_signature(obj, genes = c("CD3D", "CD3E"), name = "T_cell_score")
#' }
#'
#' @importFrom Seurat AddModuleScore
#' @export
score_gene_signature <- function(obj, genes, name, ctrl = 100) {
  # Filter to available genes
  available <- intersect(genes, rownames(obj))
  missing <- setdiff(genes, rownames(obj))

  if (length(missing) > 0) {
    message(sprintf("Note: %d genes not found: %s",
                    length(missing), paste(head(missing, 3), collapse = ", ")))
  }

  if (length(available) < 2) {
    warning("Fewer than 2 genes available for scoring")
    obj@meta.data[[name]] <- NA_real_
    return(obj)
  }

  obj <- Seurat::AddModuleScore(
    obj,
    features = list(available),
    name = name,
    ctrl = ctrl,
    seed = 42
  )

  # AddModuleScore appends "1" to name
  colnames(obj@meta.data)[colnames(obj@meta.data) == paste0(name, "1")] <- name

  return(obj)
}


#' Score multiple gene modules
#'
#' Scores cells for multiple gene signatures by calling
#' \code{\link{score_gene_signature}} for each module.
#'
#' @param obj A Seurat object.
#' @param modules Named list of character vectors, where each element is a
#'   set of gene names for one module.
#' @param prefix Character. Prefix for score column names (default: "module_").
#'
#' @return The Seurat object with score columns added for each module.
#'
#' @examples
#' \dontrun{
#' modules <- list(
#'   inflammatory = c("IL1B", "TNF", "IL6"),
#'   proliferation = c("MKI67", "TOP2A")
#' )
#' obj <- score_multiple_modules(obj, modules)
#' }
#'
#' @export
score_multiple_modules <- function(obj, modules, prefix = "module_") {
  for (module_name in names(modules)) {
    score_name <- paste0(prefix, module_name)
    message("Scoring module: ", module_name)
    obj <- score_gene_signature(obj, modules[[module_name]], score_name)
  }
  return(obj)
}


#' P-value from permutation null distribution
#'
#' Computes an empirical p-value by comparing an observed statistic to a
#' null distribution from permutations.
#'
#' @param observed Numeric. The observed test statistic.
#' @param null_distribution Numeric vector. Null statistics from permutations.
#' @param alternative Character. Type of test: \code{"two.sided"} (default),
#'   \code{"greater"}, or \code{"less"}.
#'
#' @return Numeric. The empirical p-value, corrected for finite permutations
#'   (minimum: \code{1 / (n_perms + 1)}).
#'
#' @examples
#' \dontrun{
#' null <- rnorm(1000)
#' permutation_pvalue(2.5, null)            # two-sided
#' permutation_pvalue(2.5, null, "greater") # one-sided
#' }
#'
#' @export
permutation_pvalue <- function(observed, null_distribution, alternative = "two.sided") {
  n_perms <- length(null_distribution)

  p <- switch(alternative,
    "two.sided" = sum(abs(null_distribution) >= abs(observed)) / n_perms,
    "greater" = sum(null_distribution >= observed) / n_perms,
    "less" = sum(null_distribution <= observed) / n_perms,
    stop("Unknown alternative: ", alternative)
  )

  # Correct for finite permutations
  p <- max(p, 1 / (n_perms + 1))

  return(p)
}


#' Calculate neighborhood entropy for each cell
#'
#' Computes the Shannon entropy of cell type composition in each cell's
#' neighborhood, based on a precomputed kNN graph.
#'
#' @param cell_types Character vector of cell type labels for all cells.
#' @param knn_result List. Result from \code{\link{build_knn_graph}}, containing
#'   an \code{indices} matrix.
#'
#' @return Numeric vector of entropy values (one per cell). Higher entropy
#'   indicates more diverse neighborhoods; lower entropy indicates more
#'   homogeneous neighborhoods.
#'
#' @details For each cell, looks up its k nearest neighbors in the kNN graph,
#'   counts the cell types among those neighbors, and computes Shannon entropy
#'   (log2 scale). The cell itself is not included in its own neighborhood
#'   (assumed excluded by \code{build_knn_graph}).
#'
#' @examples
#' \dontrun{
#' coords <- matrix(runif(200), ncol = 2)
#' types <- sample(c("A", "B", "C"), 100, replace = TRUE)
#' knn <- build_knn_graph(coords, k = 10)
#' entropies <- calculate_neighborhood_entropy(types, knn)
#' }
#'
#' @export
calculate_neighborhood_entropy <- function(cell_types, knn_result) {
  n <- nrow(knn_result$indices)
  entropies <- numeric(n)

  for (i in seq_len(n)) {
    neighbor_types <- cell_types[knn_result$indices[i, ]]
    counts <- table(neighbor_types)
    entropies[i] <- shannon_entropy(counts)
  }

  return(entropies)
}
