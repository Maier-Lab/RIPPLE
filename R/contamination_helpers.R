#' Family-level ambient-RNA gene patterns
#'
#' Return the regex used to flag broadly-expressed ambient-RNA gene families:
#' immunoglobulin chains (IGH / IGK / IGL, J-chain), ribosomal proteins
#' (RPL, RPS, MRPL, MRPS), and mitochondrial transcripts (MT-/mt-). These
#' families are universal ambient-RNA culprits in spatial transcriptomics
#' and are not cell-type-specific signal, so they should be flagged as
#' contamination even when the cross-cell-type heuristic
#' (\code{classify_gene_specificity()}) doesn't catch them (e.g. when
#' their source cell type was excluded for low replicate counts).
#'
#' Symbols are case-sensitive in most reference annotations: human
#' symbols are uppercase (\code{IGHA1}, \code{RPL10}), mouse are
#' title-case (\code{Igha1}, \code{Rpl10}). \code{species = "both"} (the
#' safe default) matches either convention.
#'
#' @param species One of \code{"human"}, \code{"mouse"}, or \code{"both"}.
#'   Default \code{"both"} (case-insensitive match for the same family
#'   roots).
#' @return A character scalar -- a Perl-compatible regular expression
#'   anchored at the start of the symbol.
#' @seealso \code{\link{find_ambient_family_genes}} to apply the pattern
#'   directly to a vector of gene symbols.
#' @export
#' @examples
#' default_ambient_family_pattern("human")
#' default_ambient_family_pattern("mouse")
#' default_ambient_family_pattern("both")
default_ambient_family_pattern <- function(species = c("both", "human", "mouse")) {
  species <- match.arg(species)
  human <- "IG[HKL]|JCHAIN|RP[LS][0-9]|MRP[LS][0-9]|MT-"
  mouse <- "Ig[hkl]|Jchain|Rp[ls][0-9]|Mrp[ls][0-9]|mt-"
  body <- switch(species,
                 human = human,
                 mouse = mouse,
                 both  = paste(human, mouse, sep = "|"))
  paste0("^(", body, ")")
}

#' Find ambient-family gene symbols in a vector
#'
#' Applies \code{\link{default_ambient_family_pattern}} to a vector of
#' gene symbols and returns the matches. Use this in contamination-filter
#' construction so that the Fig 3 (Xenium / mouse) and Fig 4 (CosMx /
#' human) pipelines cannot drift on the family regex.
#'
#' @param genes Character vector of gene symbols.
#' @param species One of \code{"human"}, \code{"mouse"}, \code{"both"}, or
#'   \code{"auto"}. \code{"auto"} (default) infers species from the
#'   majority case convention of the input (treats >70\% all-uppercase
#'   symbols as human, otherwise mouse / both).
#' @return Character vector -- the subset of \code{genes} that match.
#' @export
#' @examples
#' find_ambient_family_genes(c("IGHA1", "EPCAM", "RPL10", "TP53"))
#' find_ambient_family_genes(c("Igha1", "Epcam", "Rpl10", "Trp53"))
find_ambient_family_genes <- function(genes,
                                       species = c("auto", "both",
                                                   "human", "mouse")) {
  species <- match.arg(species)
  if (species == "auto") {
    upper_frac <- mean(grepl("^[A-Z0-9_.\\-]+$", genes, perl = TRUE),
                       na.rm = TRUE)
    species <- if (isTRUE(upper_frac > 0.7)) "human" else "both"
  }
  pattern <- default_ambient_family_pattern(species)
  grep(pattern, genes, value = TRUE, perl = TRUE)
}
