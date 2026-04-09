#' @title RIPPLE Atlas Visualization
#'
#' @description Generates multi-panel publication-quality figures from merged
#'   RIPPLE results, including volcanos, dotplots, heatmaps, and pathway
#'   enrichment panels.
#'
#' @name atlas
NULL


# ============================================================================
# Internal helpers
# ============================================================================

#' Build a diverging color palette
#' @noRd
.diverging_palette <- function(n = 100) {
  grDevices::colorRampPalette(c("#B2182B", "white", "#2166AC"))(n)
}


#' Build the gene counts stacked bar chart (Panel 1)
#' @noRd
.panel_gene_counts <- function(sig_results, coef_col, fdr_threshold,
                                query_label, contamination_threshold) {
  induced_label <- paste0(query_label, "-induced")
  repressed_label <- paste0(query_label, "-repressed")

  count_by_dir <- sig_results[, .(
    total = .N,
    specific = sum(specificity_class %in% c("specific", "moderate")),
    contamination = sum(specificity_class %in% c("ubiquitous", "contamination"))
  ), by = .(cell_type, direction)]

  count_long <- data.table::melt(
    count_by_dir,
    id.vars = c("cell_type", "direction"),
    measure.vars = c("specific", "contamination"),
    variable.name = "type", value.name = "count"
  )

  direction_colors <- stats::setNames(
    c("#D95F02", "#FDAE6B", "#1B9E77", "#A1D99B"),
    c(paste0("specific.", induced_label),
      paste0("contamination.", induced_label),
      paste0("specific.", repressed_label),
      paste0("contamination.", repressed_label))
  )

  direction_labels <- stats::setNames(
    c("Induced (specific)", "Induced (ubiquitous)",
      "Repressed (specific)", "Repressed (ubiquitous)"),
    names(direction_colors)
  )

  ggplot2::ggplot(count_long,
    ggplot2::aes(x = .data$cell_type, y = .data$count,
                 fill = interaction(.data$type, .data$direction))) +
    ggplot2::geom_col(position = "stack", width = 0.7) +
    ggplot2::scale_fill_manual(values = direction_colors,
                                labels = direction_labels, name = NULL) +
    ggplot2::labs(
      title = "Significant Spatial Gradients per Cell Type",
      subtitle = paste0("FDR < ", fdr_threshold,
                        " | Faded = potential contamination (sig in >=",
                        contamination_threshold, " cell types)"),
      x = NULL, y = "Number of significant genes"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9, color = "grey40"),
      legend.position = "top",
      legend.text = ggplot2::element_text(size = 8)
    )
}


#' Build the dotplot of cell-type-specific genes (Panel 2)
#' @noRd
.panel_dotplot <- function(all_results, sig_results, coef_col, sig_col,
                            fdr_threshold, query_label,
                            top_per_ct = 5, top_n = 50) {

  specific_hits <- sig_results[specificity_class %in% c("specific", "moderate")]
  if (nrow(specific_hits) < 2) return(NULL)

  dotplot_genes <- specific_hits[
    order(-abs(get(coef_col))),
    utils::head(.SD, top_per_ct),
    by = cell_type
  ]$gene
  dotplot_genes <- unique(dotplot_genes)
  if (length(dotplot_genes) > top_n) dotplot_genes <- dotplot_genes[seq_len(top_n)]
  if (length(dotplot_genes) < 2) return(NULL)

  dot_data <- all_results[gene %in% dotplot_genes]
  dot_data[, neg_log10_fdr := -log10(pmax(get(sig_col), 1e-50))]
  dot_data[, neg_log10_fdr_capped := pmin(neg_log10_fdr, 20)]

  coef_limit <- max(abs(dot_data[[coef_col]]), na.rm = TRUE)
  palette <- .diverging_palette(100)

  ggplot2::ggplot(dot_data, ggplot2::aes(x = .data$cell_type, y = .data$gene)) +
    ggplot2::geom_point(ggplot2::aes(
      size = .data$neg_log10_fdr_capped,
      color = .data[[coef_col]]
    )) +
    ggplot2::scale_color_gradientn(
      colors = palette,
      limits = c(-coef_limit, coef_limit),
      name = "Log-rate\ncoefficient"
    ) +
    ggplot2::scale_size_continuous(
      range = c(0.5, 5),
      name = expression(-log[10](FDR)),
      breaks = c(1, 5, 10, 20)
    ) +
    ggplot2::labs(
      title = "Cell-Type-Specific Spatial Gradients",
      subtitle = paste0("Top ", top_per_ct,
                        " specific genes per cell type | ",
                        "Red = induced near ", query_label),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 7),
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9, color = "grey40"),
      panel.grid.major = ggplot2::element_line(color = "grey92")
    )
}


#' Build multi-panel volcano (Panel 3)
#' @noRd
.panel_multi_volcano <- function(all_results, sig_results, coef_col, sig_col,
                                  fdr_threshold, query_label,
                                  contamination_genes) {

  induced_label <- paste0(query_label, "-induced")
  repressed_label <- paste0(query_label, "-repressed")
  direction_colors <- stats::setNames(
    c("#B2182B", "#2166AC"),
    c(induced_label, repressed_label)
  )

  volcano_data <- data.table::copy(all_results)
  volcano_data[, neg_log10_fdr := -log10(pmax(get(sig_col), 1e-50))]
  volcano_data[, is_contamination := gene %in% contamination_genes]

  # Label top specific genes per cell type
  label_genes <- sig_results[!gene %in% contamination_genes][
    order(-abs(get(coef_col))), utils::head(.SD, 20), by = cell_type
  ]
  label_set <- unique(label_genes[, .(gene, cell_type)])
  volcano_data[, show_label := FALSE]
  for (i in seq_len(nrow(label_set))) {
    volcano_data[gene == label_set$gene[i] & cell_type == label_set$cell_type[i],
                 show_label := TRUE]
  }

  x_max <- max(abs(volcano_data[[coef_col]]), na.rm = TRUE) * 1.05
  y_max <- max(volcano_data$neg_log10_fdr, na.rm = TRUE) * 1.05

  p <- ggplot2::ggplot(volcano_data,
    ggplot2::aes(x = .data[[coef_col]], y = .data$neg_log10_fdr)) +
    # Non-significant
    ggplot2::geom_point(
      data = volcano_data[get(sig_col) >= fdr_threshold & !is_contamination],
      color = "grey80", size = 0.3, alpha = 0.3
    ) +
    # Contamination
    ggplot2::geom_point(
      data = volcano_data[is_contamination == TRUE & get(sig_col) < fdr_threshold],
      color = "#FDAE6B", size = 1.2, alpha = 0.6, shape = 4, stroke = 0.6
    ) +
    # Significant
    ggplot2::geom_point(
      data = volcano_data[get(sig_col) < fdr_threshold & !is_contamination],
      ggplot2::aes(color = .data$direction), size = 0.5, alpha = 0.4
    ) +
    ggrepel::geom_text_repel(
      data = volcano_data[show_label == TRUE],
      ggplot2::aes(label = .data$gene),
      size = 2.2, max.overlaps = 20, segment.size = 0.2,
      min.segment.length = 0, fontface = "italic"
    ) +
    ggplot2::geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed",
                         color = "grey50", linewidth = 0.3) +
    ggplot2::scale_color_manual(
      values = c(direction_colors, "ns" = "grey70"), guide = "none"
    ) +
    ggplot2::coord_cartesian(xlim = c(-x_max, x_max), ylim = c(0, y_max)) +
    ggplot2::facet_wrap(~ cell_type, ncol = 4, scales = "free_y") +
    ggplot2::labs(
      title = "Spatial Gradient Volcanos Across Cell Types",
      subtitle = paste0("Labels = top cell-type-specific genes | ",
                        "x = potential contamination (sig in >= many cell types)"),
      x = paste0("Coefficient [", coef_col,
                 "] (negative = ", query_label, "-induced)"),
      y = expression(-log[10](FDR))
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9, color = "grey40"),
      strip.text = ggplot2::element_text(face = "bold", size = 9)
    )

  p
}


#' Build specificity heatmap (Panel 6)
#' @noRd
.panel_specific_heatmap <- function(all_results, sig_results, gene_spec,
                                     coef_col, sig_col, fdr_threshold,
                                     query_label, output_path) {

  specific_sig <- sig_results[specificity_class %in% c("specific", "moderate")]
  top_specific <- specific_sig[order(-abs(get(coef_col))),
                                utils::head(.SD, 8), by = cell_type]$gene
  top_specific <- unique(top_specific)

  if (length(top_specific) < 3) return(invisible(NULL))

  heatmap_data <- data.table::dcast(
    all_results[gene %in% top_specific],
    gene ~ cell_type, value.var = coef_col, fill = 0
  )
  heatmap_mat <- as.matrix(heatmap_data[, -1, with = FALSE])
  rownames(heatmap_mat) <- heatmap_data$gene

  sig_annot <- data.table::dcast(
    all_results[gene %in% top_specific],
    gene ~ cell_type, value.var = sig_col, fill = 1
  )
  sig_mat <- as.matrix(sig_annot[, -1, with = FALSE])
  rownames(sig_mat) <- sig_annot$gene
  display_mat <- ifelse(sig_mat < fdr_threshold, "*", "")

  max_abs <- max(abs(heatmap_mat), na.rm = TRUE)
  breaks <- seq(-max_abs, max_abs, length.out = 101)

  specificity_colors <- c("specific" = "#2166AC", "moderate" = "#92C5DE")
  row_annot <- data.frame(
    specificity = gene_spec[match(top_specific, gene)]$specificity_class,
    row.names = top_specific
  )
  annot_colors <- list(specificity = specificity_colors)

  grDevices::pdf(output_path,
                  width = max(8, ncol(heatmap_mat) * 0.8 + 2),
                  height = max(6, nrow(heatmap_mat) * 0.3 + 2))

  pheatmap::pheatmap(
    heatmap_mat,
    color = .diverging_palette(100),
    breaks = breaks,
    display_numbers = display_mat,
    fontsize_number = 8,
    cluster_rows = nrow(heatmap_mat) >= 2,
    cluster_cols = ncol(heatmap_mat) >= 2,
    annotation_row = row_annot,
    annotation_colors = annot_colors,
    main = paste0("Cell-Type-Specific Gradient Genes (red = induced near ",
                  query_label, ")"),
    fontsize_row = 8, fontsize_col = 10, angle_col = 45
  )

  grDevices::dev.off()
  invisible(output_path)
}


# ============================================================================
# Main Atlas Function
# ============================================================================

#' Generate RIPPLE atlas visualization
#'
#' Creates multi-panel publication-quality figures from merged RIPPLE results,
#' including volcanos, dotplots, heatmaps, and pathway enrichment panels.
#'
#' @param results_dir Character. Path to merged results directory (from
#'   \code{merge_ripple_results}). Must contain a \code{summary/} subdirectory
#'   with \code{all_genes_results.csv}.
#' @param output_dir Character or NULL. Output directory for figures. If NULL,
#'   defaults to \code{results_dir/plots}.
#' @param query_label Character. Display label for the query cell type.
#'   Default: "Query".
#' @param coef_col Character. Column name for the effect size coefficient.
#'   Default: "median_coef".
#' @param sig_col Character. Column name for the significance measure.
#'   Default: "fisher_fdr".
#' @param fdr_threshold Numeric. Significance cutoff. Default: 0.05.
#' @param contamination_threshold Integer. Flag genes significant in at least
#'   this many cell types as potential contamination. Default: 4.
#' @param top_per_celltype Integer. Number of top genes per cell type in the
#'   dotplot. Default: 5.
#' @param run_fgsea Logical. Whether to run pathway enrichment (requires
#'   fgsea and msigdbr packages). Default: TRUE.
#' @param gene_sets Character or named list. Gene sets for fGSEA. See
#'   \code{\link{run_ripple_fgsea}} for options. Default: "hallmark".
#' @param organism Character. "mouse" or "human" for gene set collections.
#'   Default: "mouse".
#' @param stage2_dir Character or NULL. Path to Stage 2 results directory.
#'   If provided, adds classification panels. Default: NULL.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return Invisible list of generated plot file paths.
#'
#' @details
#' The function generates the following panels (saved as individual PDFs):
#' \describe{
#'   \item{P1}{Gene counts bar chart with specificity breakdown.}
#'   \item{P2}{Dotplot of top cell-type-specific genes.}
#'   \item{P3}{Multi-panel volcano per cell type.}
#'   \item{P6}{Heatmap of top specific genes.}
#'   \item{P7}{Contamination candidates table (CSV).}
#'   \item{P8-9}{fGSEA pathway enrichment heatmap and dotplot (if
#'     \code{run_fgsea = TRUE}).}
#'   \item{P10}{Stage 2 classification breakdown (if \code{stage2_dir}
#'     is provided).}
#' }
#'
#' The results directory must contain
#' \code{summary/all_genes_results.csv} with at minimum the columns:
#' gene, cell_type, and the columns specified by \code{coef_col} and
#' \code{sig_col}.
#'
#' @examples
#' \dontrun{
#' run_ripple_atlas(
#'   results_dir = "output/distance_correlation/",
#'   query_label = "Neutrophil",
#'   organism = "mouse"
#' )
#' }
#'
#' @importFrom data.table fread fwrite copy dcast uniqueN fifelse rbindlist
#' @importFrom ggplot2 ggplot aes geom_col geom_point labs theme_bw theme
#'   element_text scale_fill_manual scale_color_manual scale_color_gradientn
#'   scale_size_continuous coord_cartesian facet_wrap geom_hline ggsave
#' @importFrom ggrepel geom_text_repel
#' @importFrom pheatmap pheatmap
#' @export
run_ripple_atlas <- function(results_dir,
                              output_dir = NULL,
                              query_label = "Query",
                              coef_col = "median_coef",
                              sig_col = "fisher_fdr",
                              fdr_threshold = 0.05,
                              contamination_threshold = 4,
                              top_per_celltype = 5,
                              run_fgsea = TRUE,
                              gene_sets = "hallmark",
                              organism = "mouse",
                              fgsea_min_genes = 100,
                              fgsea_seed = 42,
                              stage2_dir = NULL,
                              verbose = TRUE) {

  # --- Setup ---
  .msg <- function(...) if (isTRUE(verbose)) message(...)
  .ensure_dir <- function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    path
  }

  induced_label <- paste0(query_label, "-induced")
  repressed_label <- paste0(query_label, "-repressed")

  if (is.null(output_dir)) {
    output_dir <- file.path(results_dir, "plots")
  }
  .ensure_dir(output_dir)

  saved_plots <- character()

  # --- Load data ---
  .msg("=== RIPPLE Atlas ===")
  .msg("Loading results from: ", results_dir)

  results_file <- file.path(results_dir, "summary", "all_genes_results.csv")
  if (!file.exists(results_file)) {
    stop("Results file not found: ", results_file,
         "\nExpected summary/all_genes_results.csv inside results_dir.",
         call. = FALSE)
  }

  all_results <- data.table::fread(results_file)
  .msg("  Loaded ", nrow(all_results), " gene-celltype combinations")

  # Validate columns
  if (!coef_col %in% names(all_results)) {
    # Fallback: try combined_coef
    if ("combined_coef" %in% names(all_results)) {
      .msg("  Column '", coef_col, "' not found, falling back to 'combined_coef'")
      coef_col <- "combined_coef"
    } else {
      stop("Coefficient column '", coef_col, "' not found in results.",
           call. = FALSE)
    }
  }
  if (!sig_col %in% names(all_results)) {
    if ("fdr" %in% names(all_results)) {
      .msg("  Column '", sig_col, "' not found, falling back to 'fdr'")
      sig_col <- "fdr"
    } else {
      stop("Significance column '", sig_col, "' not found in results.",
           call. = FALSE)
    }
  }

  # --- Classify direction ---
  all_results[, direction := data.table::fifelse(
    get(sig_col) < fdr_threshold & get(coef_col) < 0, induced_label,
    data.table::fifelse(
      get(sig_col) < fdr_threshold & get(coef_col) > 0, repressed_label, "ns"
    )
  )]

  sig_results <- all_results[get(sig_col) < fdr_threshold]
  .msg("  Significant genes (", sig_col, " < ", fdr_threshold, "): ",
       nrow(sig_results))

  # --- Classify gene specificity ---
  .msg("Classifying gene specificity...")
  gene_spec <- classify_gene_specificity(
    all_results, fdr_col = sig_col, fdr_threshold = fdr_threshold,
    contamination_threshold = contamination_threshold
  )

  # Merge specificity into results
  all_results <- merge(
    all_results,
    gene_spec[, .(gene, n_celltypes, specificity_class)],
    by = "gene", all.x = TRUE
  )
  all_results[is.na(n_celltypes), `:=`(n_celltypes = 0L,
                                         specificity_class = "ns")]

  sig_results <- merge(
    sig_results,
    gene_spec[, .(gene, n_celltypes, specificity_class)],
    by = "gene", all.x = TRUE
  )

  .msg("  Specific: ", sum(gene_spec$specificity_class == "specific"),
       " | Moderate: ", sum(gene_spec$specificity_class == "moderate"),
       " | Contamination: ", sum(gene_spec$specificity_class == "contamination"))

  # Cell type ordering by specific gene count
  ct_order <- sig_results[
    specificity_class %in% c("specific", "moderate"),
    .N, by = cell_type
  ][order(-N)]$cell_type
  ct_remaining <- setdiff(unique(all_results$cell_type), ct_order)
  ct_order <- c(ct_order, ct_remaining)
  all_results[, cell_type := factor(cell_type, levels = ct_order)]
  sig_results[, cell_type := factor(cell_type, levels = ct_order)]

  contamination_genes <- gene_spec[specificity_class == "contamination"]$gene

  # --- Panel 1: Gene counts bar chart ---
  .msg("Panel 1: Gene counts bar chart...")
  tryCatch({
    p1 <- .panel_gene_counts(sig_results, coef_col, fdr_threshold,
                              query_label, contamination_threshold)
    f1 <- file.path(output_dir, "gene_counts_by_celltype.pdf")
    ggplot2::ggsave(f1, p1, width = 9, height = 6)
    saved_plots <- c(saved_plots, f1)
    .msg("  Saved: gene_counts_by_celltype.pdf")
  }, error = function(e) .msg("  [ERROR] Panel 1: ", conditionMessage(e)))

  # --- Panel 2: Dotplot ---
  .msg("Panel 2: Dotplot (cell-type-specific genes)...")
  tryCatch({
    p2 <- .panel_dotplot(all_results, sig_results, coef_col, sig_col,
                          fdr_threshold, query_label,
                          top_per_ct = top_per_celltype)
    if (!is.null(p2)) {
      n_genes <- length(unique(
        sig_results[specificity_class %in% c("specific", "moderate")]$gene
      ))
      f2 <- file.path(output_dir, "dotplot_specific_genes.pdf")
      ggplot2::ggsave(f2, p2, width = 10, height = max(6, n_genes * 0.22))
      saved_plots <- c(saved_plots, f2)
      .msg("  Saved: dotplot_specific_genes.pdf")
    } else {
      .msg("  [SKIP] Too few specific genes for dotplot")
    }
  }, error = function(e) .msg("  [ERROR] Panel 2: ", conditionMessage(e)))

  # --- Panel 3: Multi-panel volcano ---
  .msg("Panel 3: Multi-panel volcano...")
  tryCatch({
    n_ct <- data.table::uniqueN(all_results$cell_type)
    p3 <- .panel_multi_volcano(all_results, sig_results, coef_col, sig_col,
                                fdr_threshold, query_label,
                                contamination_genes)
    f3 <- file.path(output_dir, "multi_volcano.pdf")
    ggplot2::ggsave(f3, p3, width = 14, height = 3.5 * ceiling(n_ct / 4))
    saved_plots <- c(saved_plots, f3)
    .msg("  Saved: multi_volcano.pdf")
  }, error = function(e) .msg("  [ERROR] Panel 3: ", conditionMessage(e)))

  # --- Panel 6: Specific genes heatmap ---
  .msg("Panel 6: Specific genes heatmap...")
  tryCatch({
    f6 <- file.path(output_dir, "specific_genes_heatmap.pdf")
    .panel_specific_heatmap(all_results, sig_results, gene_spec,
                             coef_col, sig_col, fdr_threshold,
                             query_label, f6)
    if (file.exists(f6)) {
      saved_plots <- c(saved_plots, f6)
      .msg("  Saved: specific_genes_heatmap.pdf")
    } else {
      .msg("  [SKIP] Too few specific genes for heatmap")
    }
  }, error = function(e) .msg("  [ERROR] Panel 6: ", conditionMessage(e)))

  # --- Panel 7: Contamination candidates ---
  .msg("Panel 7: Contamination candidates...")
  contamination_list <- gene_spec[specificity_class == "contamination"][
    order(-n_celltypes)
  ]
  if (nrow(contamination_list) > 0) {
    f7 <- file.path(output_dir, "contamination_candidates.csv")
    data.table::fwrite(contamination_list, f7)
    saved_plots <- c(saved_plots, f7)
    .msg("  Saved: contamination_candidates.csv (",
         nrow(contamination_list), " genes)")
  } else {
    .msg("  No contamination candidates found")
  }

  # --- Panels 8-9: fGSEA pathway enrichment ---
  if (isTRUE(run_fgsea)) {
    has_fgsea <- requireNamespace("fgsea", quietly = TRUE) &&
                 requireNamespace("msigdbr", quietly = TRUE)

    if (has_fgsea) {
      .msg("Panels 8-9: fGSEA pathway enrichment...")
      tryCatch({
        fgsea_results <- run_ripple_fgsea(
          all_results, gene_sets = gene_sets, organism = organism,
          coef_col = coef_col, fdr_col = sig_col,
          min_genes = fgsea_min_genes, seed = fgsea_seed
        )

        if (nrow(fgsea_results) > 0) {
          # Save combined table
          f_gsea_csv <- file.path(output_dir, "fgsea_all_celltypes.csv")
          data.table::fwrite(fgsea_results, f_gsea_csv)
          saved_plots <- c(saved_plots, f_gsea_csv)
          .msg("  Saved: fgsea_all_celltypes.csv")

          # NES heatmap for significant pathways
          sig_pathways <- fgsea_results[padj < 0.05, unique(pathway)]
          if (length(sig_pathways) >= 2) {
            nes_wide <- data.table::dcast(
              fgsea_results[pathway %in% sig_pathways],
              pathway_clean ~ cell_type, value.var = "NES", fill = 0
            )
            nes_mat <- as.matrix(nes_wide[, -1, with = FALSE])
            rownames(nes_mat) <- nes_wide$pathway_clean

            padj_wide <- data.table::dcast(
              fgsea_results[pathway %in% sig_pathways],
              pathway_clean ~ cell_type, value.var = "padj", fill = 1
            )
            padj_mat <- as.matrix(padj_wide[, -1, with = FALSE])
            rownames(padj_mat) <- padj_wide$pathway_clean
            star_mat <- ifelse(padj_mat < 0.05, "*", "")

            max_nes <- max(abs(nes_mat), na.rm = TRUE)
            breaks <- seq(-max_nes, max_nes, length.out = 101)

            f8 <- file.path(output_dir, "fgsea_heatmap.pdf")
            grDevices::pdf(f8,
              width = max(9, ncol(nes_mat) * 0.8 + 4),
              height = max(6, nrow(nes_mat) * 0.35 + 2)
            )
            pheatmap::pheatmap(
              nes_mat, color = .diverging_palette(100), breaks = breaks,
              display_numbers = star_mat, fontsize_number = 10,
              cluster_rows = TRUE, cluster_cols = TRUE,
              main = paste0("Pathway Enrichment (NES)\nRed = induced near ",
                            query_label, " | Blue = repressed near ",
                            query_label),
              fontsize_row = 9, fontsize_col = 10, angle_col = 45
            )
            grDevices::dev.off()
            saved_plots <- c(saved_plots, f8)
            .msg("  Saved: fgsea_heatmap.pdf (", length(sig_pathways),
                 " pathways)")
          }

          # Pathway dotplot (significant pathways only)
          sig_pw_names <- fgsea_results[padj < 0.05, unique(pathway)]
          dotplot_data <- fgsea_results[pathway %in% sig_pw_names]
          if (nrow(dotplot_data) > 0) {
            dotplot_data[, neg_log10_padj := -log10(pmax(padj, 1e-20))]
            pw_order <- dotplot_data[
              , .(mean_nes = mean(NES, na.rm = TRUE)), by = pathway_clean
            ][order(mean_nes)]$pathway_clean
            dotplot_data[, pathway_clean := factor(pathway_clean,
                                                    levels = pw_order)]
            nes_lim <- max(abs(dotplot_data$NES), na.rm = TRUE)

            p9 <- ggplot2::ggplot(
              dotplot_data[padj < 0.05],
              ggplot2::aes(x = .data$cell_type, y = .data$pathway_clean)
            ) +
              ggplot2::geom_point(ggplot2::aes(
                size = .data$neg_log10_padj, color = .data$NES
              )) +
              ggplot2::scale_color_gradientn(
                colors = .diverging_palette(100),
                limits = c(-nes_lim, nes_lim), name = "NES"
              ) +
              ggplot2::scale_size_continuous(
                range = c(1.5, 6),
                name = expression(-log[10](padj))
              ) +
              ggplot2::labs(
                title = "Pathway Enrichment per Cell Type",
                subtitle = paste0("Pathways significant in >= 1 cell type",
                                  " (padj < 0.05)"),
                x = NULL, y = NULL
              ) +
              ggplot2::theme_bw(base_size = 10) +
              ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                plot.title = ggplot2::element_text(face = "bold", size = 12),
                plot.subtitle = ggplot2::element_text(size = 9,
                                                       color = "grey40")
              )

            n_pw <- length(sig_pw_names)
            f9 <- file.path(output_dir, "fgsea_dotplot.pdf")
            ggplot2::ggsave(f9, p9,
              width = max(7, n_ct * 0.5 + 3),
              height = max(5, n_pw * 0.3 + 2)
            )
            saved_plots <- c(saved_plots, f9)
            .msg("  Saved: fgsea_dotplot.pdf")
          }
        }
      }, error = function(e) .msg("  [ERROR] fGSEA panels: ",
                                   conditionMessage(e)))
    } else {
      .msg("[SKIP] Panels 8-9: fgsea/msigdbr not installed")
    }
  }

  # --- Panel 10: Stage 2 classification (optional) ---
  if (!is.null(stage2_dir)) {
    .msg("Panel 10: Stage 2 classification breakdown...")
    stage2_file <- file.path(stage2_dir, "summary", "stage2_all_results.csv")

    if (file.exists(stage2_file)) {
      tryCatch({
        stage2_results <- data.table::fread(stage2_file)
        .msg("  Stage 2 loaded: ", nrow(stage2_results), " rows")

        query_specific_label <- paste0(query_label, "_specific")
        classification_colors <- stats::setNames(
          c("#1B9E77", "#66A61E", "#E7298A", "#E6AB02", "#7570B3", "grey70"),
          c(query_specific_label, "enhanced", "niche_driven",
            "underpowered", "reversed", "not_tested")
        )

        class_by_ct <- stage2_results[, .N, by = .(cell_type, classification)]
        class_by_ct[, cell_type := factor(cell_type, levels = ct_order)]

        p10 <- ggplot2::ggplot(
          class_by_ct,
          ggplot2::aes(x = .data$cell_type, y = .data$N,
                       fill = .data$classification)
        ) +
          ggplot2::geom_col(position = "stack", width = 0.7) +
          ggplot2::scale_fill_manual(values = classification_colors,
                                     name = "Classification") +
          ggplot2::labs(
            title = paste0("Stage 2: ", query_label,
                           "-Specific vs Niche-Driven Gene Classification"),
            x = NULL, y = "Number of Stage 1 significant genes"
          ) +
          ggplot2::theme_bw(base_size = 11) +
          ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            plot.title = ggplot2::element_text(face = "bold", size = 12),
            legend.position = "top"
          )

        f10 <- file.path(output_dir, "stage2_classification_breakdown.pdf")
        ggplot2::ggsave(f10, p10, width = 10, height = 6)
        saved_plots <- c(saved_plots, f10)
        .msg("  Saved: stage2_classification_breakdown.pdf")
      }, error = function(e) .msg("  [ERROR] Panel 10: ",
                                   conditionMessage(e)))
    } else {
      .msg("  [SKIP] Stage 2 file not found: ", stage2_file)
    }
  }

  # --- Summary ---
  .msg("\n=== Atlas Complete ===")
  .msg("  Output directory: ", output_dir)
  .msg("  Panels saved: ", length(saved_plots))

  invisible(saved_plots)
}
