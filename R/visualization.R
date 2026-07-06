#' @title Visualization Functions
#'
#' @description Plotting functions for spatial data, volcano plots, decay curves,
#'   and other RIPPLE analysis visualizations.
#'
#' @name visualization
NULL

#' RIPPLE default ggplot theme
#'
#' A minimal ggplot theme used across all RIPPLE plotting functions:
#' no panel border (the "black box" around the plot is removed), no minor
#' grid lines, visible axis lines and ticks, and readable axis text sizes.
#' Derived from \code{ggplot2::theme_bw()}.
#'
#' @param base_size Numeric. Base font size in points (default: 12).
#' @param base_family Character. Base font family (default: "").
#'
#' @return A ggplot2 theme object.
#'
#' @examples
#' \dontrun{
#' ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) +
#'   ggplot2::geom_point() +
#'   theme_ripple()
#' }
#'
#' @importFrom ggplot2 theme_bw theme element_line element_blank element_text
#'   element_rect
#' @export
theme_ripple <- function(base_size = 12, base_family = "") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        color = "grey92", linewidth = 0.3
      ),
      axis.line = ggplot2::element_line(color = "grey30", linewidth = 0.4),
      axis.ticks = ggplot2::element_line(color = "grey30", linewidth = 0.4),
      axis.text = ggplot2::element_text(
        size = ggplot2::rel(0.9), color = "grey20"
      ),
      axis.title = ggplot2::element_text(
        size = ggplot2::rel(1.0), color = "grey20"
      ),
      strip.background = ggplot2::element_rect(
        fill = "grey95", color = NA
      ),
      strip.text = ggplot2::element_text(face = "bold")
    )
}


#' Single-sample spatial plot
#'
#' Creates a clean spatial plot for a single sample using \code{theme_void()}
#' and \code{coord_fixed()} for proper aspect ratio.
#'
#' @param data A \code{data.table} or \code{data.frame} with \code{x}, \code{y}
#'   (or \code{spatial_x}, \code{spatial_y}) and the color column.
#' @param color_var Character. Name of the column to color by.
#' @param title Character or NULL. Plot title.
#' @param palette Named character vector of colors (for categorical).
#' @param point_size Numeric. Point size (default: 0.5).
#' @param alpha Numeric. Transparency (default: 0.8).
#' @param continuous Logical. Whether the color variable is continuous
#'   (default: FALSE).
#' @param highlight_value Character or NULL. If set, highlights only cells with
#'   this value and greys out all others.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' plot_spatial_single(cell_data, "cell_type", title = "Sample 1")
#' plot_spatial_single(cell_data, "cell_type", highlight_value = "Tumor")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point coord_fixed theme_void theme
#'   element_text labs scale_color_manual scale_color_viridis_c
#' @importFrom data.table copy fifelse
#' @export
plot_spatial_single <- function(data, color_var, title = NULL,
                                palette = NULL, point_size = 0.5, alpha = 0.8,
                                continuous = FALSE, highlight_value = NULL) {
  # Handle coordinate column names
  if ("spatial_x" %in% names(data) && !"x" %in% names(data)) {
    data <- data.table::copy(data)
    data[, x := spatial_x]
    data[, y := spatial_y]
  }

  # If highlighting specific value, mask others
  if (!is.null(highlight_value)) {
    data <- data.table::copy(data)
    data[, plot_color := data.table::fifelse(
      get(color_var) == highlight_value, highlight_value, "Other"
    )]
    data[, plot_color := factor(plot_color, levels = c(highlight_value, "Other"))]
    color_var <- "plot_color"

    if (is.null(palette)) {
      palette <- c("#E74C3C", "lightgrey")
      names(palette) <- c(highlight_value, "Other")
    } else {
      palette <- c(palette[highlight_value], "lightgrey")
      names(palette) <- c(highlight_value, "Other")
    }
  }

  p <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data$x, y = .data$y,
    color = .data[[color_var]]
  )) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 8)
    ) +
    ggplot2::labs(title = title)

  # Apply color scale
  if (continuous) {
    p <- p + ggplot2::scale_color_viridis_c(option = "magma", name = color_var)
  } else if (!is.null(palette)) {
    p <- p + ggplot2::scale_color_manual(values = palette, name = "Cell Type")
  }

  return(p)
}


#' Multi-sample spatial plots via patchwork
#'
#' Creates individual spatial plots per sample and combines them using
#' \code{patchwork::wrap_plots}. This avoids the known ggplot2 incompatibility
#' between \code{coord_fixed()} and \code{facet_wrap(scales = "free")}.
#'
#' @param data A \code{data.table} with \code{x}, \code{y}, sample column,
#'   and color column.
#' @param color_var Character. Name of the column to color by.
#' @param sample_col Character. Name of the sample column (default: "sample").
#' @param title Character or NULL. Overall plot title.
#' @param palette Named character vector of colors.
#' @param point_size Numeric. Point size (default: 0.3).
#' @param alpha Numeric. Transparency (default: 0.6).
#' @param ncol Integer or NULL. Number of columns in layout (default: auto).
#' @param continuous Logical. Whether the color variable is continuous.
#'
#' @return A \code{patchwork} combined plot.
#'
#' @examples
#' \dontrun{
#' plot_spatial_by_sample(cell_data, "cell_type", sample_col = "sample_id")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point coord_fixed theme_void theme
#'   element_text labs scale_color_manual scale_color_viridis_c
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom data.table copy
#' @export
plot_spatial_by_sample <- function(data, color_var, sample_col = "sample",
                                   title = NULL, palette = NULL,
                                   point_size = 0.3, alpha = 0.6,
                                   ncol = NULL, continuous = FALSE) {
  # Handle coordinate column names
  if ("spatial_x" %in% names(data) && !"x" %in% names(data)) {
    data <- data.table::copy(data)
    data[, x := spatial_x]
    data[, y := spatial_y]
  }

  samples <- unique(data[[sample_col]])
  n_samples <- length(samples)

  if (is.null(ncol)) {
    ncol <- min(4, ceiling(sqrt(n_samples)))
  }

  # Create individual plots
  plot_list <- lapply(samples, function(samp) {
    samp_data <- data[get(sample_col) == samp]

    p <- ggplot2::ggplot(samp_data, ggplot2::aes(
      x = .data$x, y = .data$y,
      color = .data[[color_var]]
    )) +
      ggplot2::geom_point(size = point_size, alpha = alpha) +
      ggplot2::coord_fixed() +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 9, hjust = 0.5),
        legend.position = "none"
      ) +
      ggplot2::labs(title = samp)

    if (continuous) {
      p <- p + ggplot2::scale_color_viridis_c(option = "magma")
    } else if (!is.null(palette)) {
      p <- p + ggplot2::scale_color_manual(values = palette)
    }

    return(p)
  })

  # Combine with patchwork
  combined <- patchwork::wrap_plots(plot_list, ncol = ncol)

  if (!is.null(title)) {
    combined <- combined + patchwork::plot_annotation(title = title)
  }

  return(combined)
}


#' Gradient volcano plot
#'
#' Creates a volcano plot showing distance-expression gradient results with
#' genes colored by decay pattern and labeled if significant.
#'
#' @param results A \code{data.table} with analysis results. Must contain at
#'   minimum a coefficient column and an FDR column.
#' @param coef_col Character. Name of the coefficient column
#'   (default: "median_coef").
#' @param fdr_col Character. Name of the FDR column (default: "fisher_fdr").
#' @param fdr_threshold Numeric. Significance threshold for FDR
#'   (default: 0.05).
#' @param query_label Character. Display label for the query cell type
#'   (default: "Query").
#' @param target_label Character or NULL. Display label for the target cell
#'   type being plotted. When supplied (and \code{title} is NULL), the title
#'   becomes "<target_label> near <query_label>", which self-documents
#'   one-vs-one runs. Default \code{NULL}.
#' @param n_label Integer. Maximum number of genes to label (default: 20).
#'   Ignored when \code{label_genes} is supplied.
#' @param title Character or NULL. Plot title. If NULL, auto-generated
#'   (uses \code{target_label} when provided).
#' @param exclude_specificity_class Optional character vector. Specificity-class
#'   labels to exclude from the *label-selection pool only* (genes still
#'   appear as dots, just without text labels). Typical use:
#'   \code{exclude_specificity_class = "broad"} to skip the broad-expression
#'   class (heuristic proxy for ambient-RNA / query-signature genes flagged
#'   by the cross-cell-type heuristic). If set, the input \code{results}
#'   must carry a \code{specificity_class} column (produced by
#'   \code{classify_gene_specificity()} or \code{run_ripple_atlas()}).
#'   Default \code{NULL} preserves the prior labelling behaviour. Ignored
#'   when \code{label_genes} is supplied (the user has explicitly chosen
#'   which genes to label).
#' @param label_genes Optional character vector. When supplied, these exact
#'   genes are labelled (subject to availability in \code{results}), instead
#'   of the top-N significant by FDR. Useful for manuscript figures where
#'   you want to highlight a curated panel of positive controls or pathway
#'   markers regardless of where they fall in the FDR ranking. Genes that
#'   are missing from the data trigger a one-time message and are silently
#'   dropped. Default \code{NULL} (use top-N).
#' @param color_by_direction Logical. When \code{TRUE}, colour points by
#'   gradient direction (induced near query, repressed near query, not
#'   significant) using \code{direction_palette}, and drop the
#'   size-by-significance aesthetic in favour of a fixed point size. This
#'   produces the manuscript-style volcano with three clearly demarcated
#'   classes. Default \code{FALSE} (size-by-significance, no colour).
#' @param direction_palette Optional named character vector of length 3 with
#'   names \code{"induced"}, \code{"repressed"}, \code{"ns"}. Override the
#'   default palette (\code{induced = "#B2182B"}, \code{repressed = "#2166AC"},
#'   \code{ns = "grey75"}). Ignored when \code{color_by_direction = FALSE}.
#' @param direction_labels Optional named character vector of length 3 with
#'   names \code{"induced"}, \code{"repressed"}, \code{"ns"}. Override the
#'   legend labels. Defaults to
#'   \code{paste("Induced near", query_label)}, \code{paste("Repressed near",
#'   query_label)}, and \code{"Not significant"}. Ignored when
#'   \code{color_by_direction = FALSE}.
#' @param x_axis_label Optional character or expression. Override the default
#'   x-axis label. Pass an \code{expression()} for math/unit formatting (e.g.
#'   \code{expression("Gradient score " * beta ~ "(per " * mu * "m)")}).
#'   Default \code{NULL} uses the generic
#'   \code{"Coefficient (negative = <query_label>-induced)"}.
#' @param point_size Numeric. Fixed point size used when
#'   \code{color_by_direction = TRUE}. Default \code{0.8}. Ignored when
#'   \code{color_by_direction = FALSE} (size aesthetic is used instead).
#' @param base_size Numeric. Base font size for \code{theme_ripple()}.
#'   Default \code{12}. Pass the figure-wide \code{BASE_SIZE} from a
#'   multi-panel render so axis text scales consistently with surrounding
#'   panels (avoids double-applying \code{theme_ripple()} downstream).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' plot_gradient_volcano(
#'   results,
#'   coef_col = "median_coef",
#'   fdr_col = "fisher_fdr",
#'   query_label = "Tumor"
#' )
#'
#' # Skip broad-class genes when picking labels:
#' plot_gradient_volcano(
#'   results,
#'   exclude_specificity_class = "broad"
#' )
#'
#' # Curated label set (e.g. T zone positive controls):
#' plot_gradient_volcano(
#'   results,
#'   label_genes = c("Ccr7", "Sell", "Lef1", "Tcf7", "Cd69", "S1pr1")
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline xlim labs
#'   theme_classic theme element_text scale_size_manual
#' @importFrom ggrepel geom_text_repel
#' @importFrom data.table copy
#' @export
plot_gradient_volcano <- function(results, coef_col = "median_coef",
                                  fdr_col = "fisher_fdr",
                                  fdr_threshold = 0.05,
                                  query_label = "Query",
                                  target_label = NULL,
                                  n_label = 20,
                                  title = NULL,
                                  exclude_specificity_class = NULL,
                                  label_genes = NULL,
                                  color_by_direction = FALSE,
                                  direction_palette = NULL,
                                  direction_labels = NULL,
                                  x_axis_label = NULL,
                                  point_size = 0.8,
                                  base_size = 12) {
  plot_data <- data.table::copy(results)

  plot_data[, neg_log10_fdr := -log10(get(fdr_col))]
  plot_data[neg_log10_fdr > 50, neg_log10_fdr := 50]

  plot_data[, significant := get(fdr_col) < fdr_threshold]

  if (isTRUE(color_by_direction)) {
    default_palette <- c(induced = "#B2182B", repressed = "#2166AC",
                         ns = "grey75")
    default_labels  <- c(
      induced   = paste("Induced near", query_label),
      repressed = paste("Repressed near", query_label),
      ns        = "Not significant"
    )
    pal <- if (is.null(direction_palette)) default_palette else {
      missing_keys <- setdiff(c("induced", "repressed", "ns"),
                              names(direction_palette))
      if (length(missing_keys) > 0) {
        stop("direction_palette must have names: induced, repressed, ns",
             call. = FALSE)
      }
      direction_palette[c("induced", "repressed", "ns")]
    }
    lbl <- if (is.null(direction_labels)) default_labels else {
      missing_keys <- setdiff(c("induced", "repressed", "ns"),
                              names(direction_labels))
      if (length(missing_keys) > 0) {
        stop("direction_labels must have names: induced, repressed, ns",
             call. = FALSE)
      }
      direction_labels[c("induced", "repressed", "ns")]
    }
    plot_data[, direction := data.table::fifelse(
      significant & get(coef_col) < 0, lbl[["induced"]],
      data.table::fifelse(significant & get(coef_col) > 0, lbl[["repressed"]],
                          lbl[["ns"]]))]
    # Map labels back to colours via a named palette keyed by the labels
    # themselves (so ggplot's scale_colour_manual matches on `direction`).
    direction_color_map <- stats::setNames(
      c(pal[["induced"]], pal[["repressed"]], pal[["ns"]]),
      c(lbl[["induced"]], lbl[["repressed"]], lbl[["ns"]])
    )
  }

  if (!is.null(label_genes)) {
    # Curated label set -- use exactly the genes the user asked for, ignoring
    # the top-N FDR ranking and any specificity-class exclusion (they have
    # explicitly chosen these genes).
    requested <- unique(as.character(label_genes))
    available <- intersect(requested, plot_data$gene)
    missing <- setdiff(requested, available)
    if (length(missing) > 0) {
      message(
        sprintf("plot_gradient_volcano(): %d of %d label_genes not present in `results`: %s",
                length(missing), length(requested),
                paste(utils::head(missing, 10), collapse = ", "))
      )
    }
    top_genes <- plot_data[gene %in% available]
  } else {
    label_pool <- plot_data[significant == TRUE]
    if (!is.null(exclude_specificity_class)) {
      if (!"specificity_class" %in% names(label_pool)) {
        stop(
          "exclude_specificity_class is set but `specificity_class` is not a ",
          "column on `results`. Run classify_gene_specificity() or ",
          "run_ripple_atlas() first to add this column.",
          call. = FALSE
        )
      }
      label_pool <- label_pool[
        !specificity_class %in% exclude_specificity_class
      ]
    }
    top_genes <- utils::head(label_pool[order(get(fdr_col))], n_label)
  }

  coef_values <- plot_data[[coef_col]]
  max_score <- max(abs(coef_values), na.rm = TRUE) * 1.1

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data[[coef_col]], y = .data$neg_log10_fdr)
  )

  if (isTRUE(color_by_direction)) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(colour = .data$direction),
                          size = point_size, alpha = 0.6) +
      ggplot2::scale_colour_manual(values = direction_color_map, name = NULL)
  } else {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(size = .data$significant), alpha = 0.6) +
      ggplot2::scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2.5),
                                  guide = "none")
  }

  p <- p +
    ggplot2::geom_hline(
      yintercept = -log10(fdr_threshold), linetype = "dashed",
      color = "grey40"
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "grey60") +
    ggplot2::xlim(-max_score, max_score) +
    ggrepel::geom_text_repel(
      data = top_genes,
      ggplot2::aes(label = .data$gene),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      seed = 1,
      min.segment.length = 0
    ) +
    ggplot2::labs(
      title = if (!is.null(title)) {
        title
      } else if (!is.null(target_label)) {
        sprintf("%s near %s", target_label, query_label)
      } else {
        "Distance-Expression Gradient Volcano"
      },
      subtitle = sprintf(
        "%d genes significant (FDR < %.2f)",
        sum(plot_data$significant, na.rm = TRUE), fdr_threshold
      ),
      x = if (!is.null(x_axis_label)) x_axis_label else
            paste0("Coefficient (negative = ", query_label, "-induced)"),
      y = "-log10(FDR)"
    ) +
    theme_ripple(base_size = base_size) +
    ggplot2::theme(
      legend.position = if (isTRUE(color_by_direction)) "bottom" else "right",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  return(p)
}


#' Gradient curve plot
#'
#' Plots a per-cell expression statistic as a function of distance from
#' query cells. The shape can be monotonic decay, monotonic increase,
#' biphasic, or noisy -- hence "gradient curve" rather than "decay curve".
#' Two modes:
#' \itemize{
#'   \item \strong{Per-sample mode} (recommended): each biological
#'     replicate is drawn as a faint line and the cross-sample mean is
#'     overlaid in bold, with a 95\% CI ribbon (mean +/- 1.96 * SE across
#'     samples per bin, with samples weighted equally). Triggered when
#'     \code{sample_col} is supplied, or auto-detected when
#'     \code{bin_stats} contains a \code{sample_id} column (the default
#'     output of \code{bin_decay_data(sample_ids = ...)}). This matches
#'     the package's \emph{N = samples, not cells} philosophy and is the
#'     mode used in the companion manuscript.
#'   \item \strong{Pooled mode}: one curve from a pre-pooled bin table,
#'     with a 95\% CI ribbon derived from a per-bin SE column
#'     (mean +/- 1.96 * SE). The ribbon here is a binomial Wald CI on the
#'     pooled-cell proportion (denominator = total cells, not number of
#'     samples), which \strong{overstates precision} when replicates
#'     disagree and lets samples with more cells dominate the curve.
#'     Kept for back-compat; prefer per-sample mode for any
#'     manuscript-bound figure.
#' }
#' To switch modes: call \code{bin_decay_data()} with \code{sample_ids}
#' for per-sample, or without it for pooled. \code{plot_gradient_curve()}
#' picks the right mode from the input table.
#' The gene's meta-analysis gradient score and FDR are shown in the
#' subtitle.
#'
#' @param bin_stats A \code{data.table} of binned decay statistics. In
#'   pooled mode, must contain the columns named by \code{x_col},
#'   \code{y_col}, optionally \code{se_col} for the ribbon, and
#'   optionally \code{n_cells} for the point-size encoding (see
#'   "Dot size encoding" below). In per-sample mode, must additionally
#'   contain \code{sample_col}; per-bin means and SEs are computed
#'   internally and \code{n_cells}, if present, is summed across samples
#'   per bin to drive the dot-size encoding.
#' @param gene_name Character. Gene name for the plot title.
#' @param cell_type Character. Cell type name for context and for looking up
#'   the meta-analysis row when \code{results} is passed.
#' @param gradient_score Numeric or NULL. Per-gene meta-analysis gradient
#'   score (a.k.a. \code{median_coef} / \code{combined_coef}). Shown in the
#'   subtitle. If \code{NULL} and \code{results} is supplied, looked up
#'   automatically from \code{results}.
#' @param fdr Numeric or NULL. Per-gene Fisher FDR. Shown in the subtitle.
#'   If \code{NULL} and \code{results} is supplied, looked up automatically.
#' @param results Optional \code{data.table} of RIPPLE results (e.g. output of
#'   \code{merge_ripple_results()}). If supplied and \code{gradient_score} /
#'   \code{fdr} are \code{NULL}, the row matching \code{gene_name} and
#'   \code{cell_type} is used. Must contain \code{gene}, \code{cell_type},
#'   and one of \code{median_coef} / \code{gradient_score} /
#'   \code{combined_coef}, and one of \code{fisher_fdr} / \code{fdr}.
#' @param sample_col Character or NULL. Name of the per-sample identifier
#'   column in \code{bin_stats}. When non-NULL, switches to per-sample
#'   mode: thin lines per sample plus a bold mean line with a 95\% CI
#'   ribbon (mean +/- 1.96 * SE across samples per bin). When \code{NULL}
#'   (default), the function auto-detects per-sample mode if
#'   \code{bin_stats} contains a \code{sample_id} column (the default
#'   output of \code{bin_decay_data(sample_ids = ...)}); otherwise it
#'   falls back to pooled mode.
#' @param x_col Character. Name of the distance / bin-centre column on
#'   \code{bin_stats}. Default \code{"bin_center"} (matches
#'   \code{bin_decay_data()} output). Falls back to \code{"dist_mid"} or
#'   \code{"bin_mid_um"} if the default isn't present.
#' @param y_col Character. Name of the response column (proportion
#'   expressing, mean rate, etc.) on \code{bin_stats}. Default
#'   \code{"prop_expressing"}.
#' @param se_col Character. Name of the per-bin SE column used for the
#'   ribbon in pooled mode. Ignored in per-sample mode (SE is computed
#'   across samples). Default \code{"se_prop"} (matches
#'   \code{bin_decay_data()}).
#' @param min_cells_per_bin Integer. Per-sample mode only. Drops
#'   sample-bin observations supported by fewer than this many cells
#'   before aggregating across samples. Default \code{10L}, matching the
#'   HyMy companion-manuscript script. Requires a \code{n_cells} column
#'   on the input; silently ignored if absent. Set to \code{0L} to
#'   disable.
#' @param min_samples_per_bin Integer. Per-sample mode only. Drops bins
#'   backed by fewer than this many samples after the
#'   \code{min_cells_per_bin} filter. Default \code{2L}: at least two
#'   samples must contribute to a bin for it to appear in the mean curve.
#'   Set to \code{0L} to disable.
#' @param y_lab Character or NULL. Y-axis label. Default \code{NULL} ->
#'   "P(expressing)" in pooled mode, "Mean expression rate" in per-sample
#'   mode.
#' @param query_label Character. Display label for the query cell type
#'   (default: "Query").
#' @param max_distance Numeric. Maximum distance for x-axis (default: 200).
#' @param color Character or \code{NULL}. Line / ribbon / point colour.
#'   Default \code{NULL} picks a sign-aware colour based on
#'   \code{gradient_score}: red (\code{#B2182B}) when negative (induced
#'   near query), blue (\code{#2166AC}) when positive (repressed near
#'   query). Falls back to \code{#E74C3C} if the gradient sign cannot
#'   be determined (no \code{gradient_score} and no \code{results}
#'   lookup). Matches the \code{plot_gradient_volcano()} direction
#'   palette so the volcano and curve panels read consistently.
#' @param sample_alpha Numeric. Alpha for the per-sample lines in
#'   per-sample mode (default \code{0.4}).
#' @param sample_linewidth Numeric. Linewidth for per-sample lines in
#'   per-sample mode (default \code{0.4}).
#'
#' @return A \code{ggplot} object, or NULL if \code{bin_stats} is NULL or
#'   empty.
#'
#' @section CI interpretation -- the two modes are not the same statistic:
#' Both modes draw a \code{mean +/- 1.96 * SE} ribbon, but the SE is
#' constructed differently and the resulting intervals answer different
#' questions. \strong{Read this before reporting confidence intervals.}
#' \describe{
#'   \item{\strong{Pooled mode (binomial Wald CI of a pooled-cell
#'     proportion)}}{The \code{se_col} is taken as-is from
#'     \code{bin_stats}. When the input comes from
#'     \code{\link{bin_decay_data}} (or the original HyMy
#'     distance-correlation analysis), this column is the classical
#'     binomial standard error \code{sqrt(p * (1 - p) / n)} of the
#'     proportion expressing, with \code{n} = total cells in the bin
#'     pooled across all samples. The ribbon is therefore a normal-
#'     approximation 95\% Wald interval for that pooled-cell proportion,
#'     clamped to the 0-to-1 range. This matches the original HyMy plot exactly,
#'     but it \strong{overstates precision when there is between-sample
#'     variability}: the effective N is closer to the number of
#'     replicates than the number of cells.}
#'   \item{\strong{Per-sample mode (cross-sample CI of the mean)}}{The
#'     SE is computed inside the function as
#'     \code{sd_across_samples / sqrt(n_samples)} per bin, using the
#'     per-sample \code{y_col} values supplied in \code{bin_stats}. The
#'     ribbon is a normal-approximation 95\% CI for the cross-sample
#'     \emph{mean} of \code{y_col}. This captures between-replicate
#'     variability that pooled mode hides, and is the more honest
#'     interval when replicate consistency matters (e.g. for the paper).
#'     With small N (3-5 replicates) the normal approximation is rough;
#'     the per-sample lines themselves are the more reliable visual.}
#' }
#' Neither mode corrects for spatial autocorrelation between cells in
#' the same bin; both inherit the independence assumption of the
#' underlying RIPPLE Poisson GLM.
#'
#' @section Dot size encoding:
#' The bold dots on the curve are sized by the number of cells
#' contributing to each bin. The size legend is hidden by design (the
#' encoding is informational, not a primary axis), but the visual
#' meaning is consistent across modes:
#' \describe{
#'   \item{\strong{Pooled mode}}{One dot per bin at \code{(x_col, y_col)};
#'     size = the bin's \code{n_cells} value (all cells, all samples
#'     pooled). Bigger dot = tighter mean estimate.}
#'   \item{\strong{Per-sample mode}}{One dot per bin at
#'     \code{(x_col, mean across samples)}; size = the sum of
#'     \code{n_cells} across samples in that bin (i.e. total cells
#'     contributing to the cross-sample mean). Faint per-sample lines
#'     have no points -- only the bold mean is dotted. Bigger dot =
#'     more cells contributed to that bin's cross-sample mean.}
#' }
#' Why this matters: bins close to the query population are typically
#' supported by many more cells than bins at the distance cap, so dots
#' taper as distance increases. A small dot at the right edge of the
#' plot is a signal that the apparent value there is supported by few
#' observations. If the input \code{bin_stats} table has no
#' \code{n_cells} column, all points render at a single small size and
#' this encoding is silently absent.
#'
#' @examples
#' \dontrun{
#' # Pooled mode
#' plot_gradient_curve(
#'   bin_stats      = binned_data,
#'   gene_name      = "Cxcl12",
#'   cell_type      = "LEC",
#'   gradient_score = -0.005,
#'   fdr            = 0.001
#' )
#'
#' # Per-sample mode: thin lines per replicate + bold mean
#' plot_gradient_curve(
#'   bin_stats   = curve_per_sample,        # has a sample_id column
#'   gene_name   = "Ccr7",
#'   cell_type   = "T_cell",
#'   sample_col  = "sample_id",
#'   x_col       = "bin_mid_um",
#'   y_col       = "mean_rate",
#'   y_lab       = "Mean expression rate (per UMI)",
#'   results     = merged_results
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point
#'   scale_size_continuous scale_x_continuous ylim labs theme_bw theme
#'   element_text
#' @importFrom data.table as.data.table copy setnames
#' @export
plot_gradient_curve <- function(bin_stats, gene_name, cell_type,
                                gradient_score = NULL, fdr = NULL,
                                results = NULL,
                                sample_col = NULL,
                                x_col = "bin_center",
                                y_col = "prop_expressing",
                                se_col = "se_prop",
                                y_lab = NULL,
                                query_label = "Query",
                                max_distance = 200,
                                color = NULL,
                                sample_alpha = 0.4,
                                sample_linewidth = 0.4,
                                min_cells_per_bin = 10L,
                                min_samples_per_bin = 2L) {
  if (is.null(bin_stats) || nrow(bin_stats) == 0) {
    return(NULL)
  }

  # Auto-lookup from results table if stats not provided directly
  if ((is.null(gradient_score) || is.null(fdr)) && !is.null(results)) {
    res <- data.table::as.data.table(results)
    target_ct <- cell_type
    target_gene <- gene_name
    row <- res[gene == target_gene & cell_type == target_ct]
    if (nrow(row) >= 1) {
      if (is.null(gradient_score)) {
        coef_col <- intersect(
          c("median_coef", "gradient_score", "combined_coef"),
          names(row)
        )[1]
        if (!is.na(coef_col)) gradient_score <- row[[coef_col]][1]
      }
      if (is.null(fdr)) {
        fdr_col <- intersect(c("fisher_fdr", "fdr"), names(row))[1]
        if (!is.na(fdr_col)) fdr <- row[[fdr_col]][1]
      }
    }
  }

  # Sign-aware default colour: red = induced near query (negative coef),
  # blue = repressed near query (positive coef). Matches the
  # plot_gradient_volcano() direction palette. If the caller passes
  # `color` explicitly, that wins.
  if (is.null(color)) {
    if (!is.null(gradient_score) && !is.na(gradient_score)) {
      color <- if (gradient_score < 0) "#B2182B" else "#2166AC"
    } else {
      color <- "#E74C3C"  # legacy default when sign is unknown
    }
  }

  sub_text <- ""
  if (!is.null(gradient_score) && !is.null(fdr)) {
    sub_text <- sprintf(
      "Gradient score = %.4f  |  FDR = %.1e", gradient_score, fdr
    )
  }

  dt <- data.table::as.data.table(data.table::copy(bin_stats))

  # Auto-detect per-sample input. If the caller didn't pass sample_col but
  # bin_stats has a sample_id column, switch on per-sample mode automatically.
  # Matches bin_decay_data(sample_ids = ...) output directly.
  if (is.null(sample_col) && "sample_id" %in% names(dt)) {
    sample_col <- "sample_id"
  }

  # Fall back from a default x_col that isn't in the table to one that is.
  if (!(x_col %in% names(dt))) {
    for (cand in c("bin_center", "dist_mid", "bin_mid_um")) {
      if (cand %in% names(dt)) {
        x_col <- cand
        break
      }
    }
  }
  # Same for the pooled-mode SE column.
  if (!(se_col %in% names(dt))) {
    for (cand in c("se_prop", "se")) {
      if (cand %in% names(dt)) {
        se_col <- cand
        break
      }
    }
  }

  required <- c(x_col, y_col)
  if (!is.null(sample_col)) required <- c(required, sample_col)
  miss <- setdiff(required, names(dt))
  if (length(miss) > 0) {
    stop("Missing required columns in bin_stats: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }

  per_sample <- !is.null(sample_col)

  if (per_sample) {
    # Filter sample-bin observations supported by fewer than
    # min_cells_per_bin cells, then aggregate. Matches the HyMy
    # companion-manuscript distance-correlation script: low-cell bins
    # produce unstable per-sample proportions (a bin with 3 cells gives
    # prop_expressing of 0, 1/3, 2/3, or 1) and the cross-sample mean
    # gets dragged around by those.
    has_n_cells <- "n_cells" %in% names(dt)
    if (has_n_cells && min_cells_per_bin > 0L) {
      dt <- dt[get("n_cells") >= min_cells_per_bin]
    }

    # Compute mean curve and across-sample SE per bin
    mean_dt <- if (has_n_cells) {
      dt[, list(
        mean_val   = mean(get(y_col), na.rm = TRUE),
        se_val     = stats::sd(get(y_col), na.rm = TRUE) /
                       sqrt(sum(!is.na(get(y_col)))),
        ncells_sum = sum(get("n_cells"), na.rm = TRUE),
        n_samples  = sum(!is.na(get(y_col)))
      ), by = c(x_col)]
    } else {
      dt[, list(
        mean_val   = mean(get(y_col), na.rm = TRUE),
        se_val     = stats::sd(get(y_col), na.rm = TRUE) /
                       sqrt(sum(!is.na(get(y_col)))),
        n_samples  = sum(!is.na(get(y_col)))
      ), by = c(x_col)]
    }

    # Drop bins backed by too few samples -- the across-sample SD is
    # undefined or degenerate with a single observation, and a single
    # mouse driving the tail of the curve looks more confident than the
    # data warrants.
    if (min_samples_per_bin > 0L) {
      mean_dt <- mean_dt[n_samples >= min_samples_per_bin]
    }

    if (is.null(y_lab)) y_lab <- "Mean expression rate"
    y_max <- max(c(dt[[y_col]], mean_dt$mean_val + 1.96 * mean_dt$se_val),
                 na.rm = TRUE)

    p <- ggplot2::ggplot()
    # Faint per-sample lines (background)
    p <- p + ggplot2::geom_line(
      data = dt,
      ggplot2::aes(
        x = .data[[x_col]], y = .data[[y_col]],
        group = .data[[sample_col]]
      ),
      colour = color, alpha = sample_alpha, linewidth = sample_linewidth
    )
    # Mean +/- 1.96 * SE ribbon (95% CI of the mean across samples)
    p <- p + ggplot2::geom_ribbon(
      data = mean_dt,
      ggplot2::aes(x = .data[[x_col]],
                   ymin = pmax(0, .data$mean_val - 1.96 * .data$se_val),
                   ymax = .data$mean_val + 1.96 * .data$se_val),
      fill = color, alpha = 0.18
    )
    # Bold mean line (foreground)
    p <- p + ggplot2::geom_line(
      data = mean_dt,
      ggplot2::aes(x = .data[[x_col]], y = .data$mean_val),
      colour = color, linewidth = 1.1
    )
    if (has_n_cells) {
      p <- p + ggplot2::geom_point(
        data = mean_dt,
        ggplot2::aes(x = .data[[x_col]], y = .data$mean_val,
                     size = .data$ncells_sum),
        colour = color, alpha = 0.85
      ) +
        ggplot2::scale_size_continuous(range = c(0.8, 3.5), guide = "none")
    }
  } else {
    # Pooled mode (back-compat)
    if (is.null(y_lab)) y_lab <- "P(expressing)"
    has_se <- se_col %in% names(dt)
    y_max  <- max(dt[[y_col]], na.rm = TRUE) * 1.3

    p <- ggplot2::ggplot(dt, ggplot2::aes(
      x = .data[[x_col]], y = .data[[y_col]]
    ))
    if (has_se) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = pmax(0, .data[[y_col]] - 1.96 * .data[[se_col]]),
          ymax = pmin(1, .data[[y_col]] + 1.96 * .data[[se_col]])
        ),
        fill = color, alpha = 0.15
      )
    }
    p <- p +
      ggplot2::geom_line(color = color, linewidth = 0.8) +
      ggplot2::geom_point(
        ggplot2::aes(size = .data$n_cells),
        color = color, alpha = 0.7
      ) +
      ggplot2::scale_size_continuous(range = c(0.8, 3.5), guide = "none")
  }

  p <- p +
    ggplot2::scale_x_continuous(breaks = seq(0, max_distance, by = 50)) +
    ggplot2::ylim(0, min(1, y_max)) +
    ggplot2::labs(
      # Name both the gene and the target cell type; the query is on the x-axis.
      title    = if (!is.null(cell_type) && nzchar(cell_type)) {
        paste0(gene_name, " in ", cell_type)
      } else {
        gene_name
      },
      subtitle = sub_text,
      x        = paste0("Distance to ", query_label, " (um)"),
      y        = y_lab
    ) +
    theme_ripple(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold.italic", size = 11),
      plot.subtitle = ggplot2::element_text(
        size = 9, color = "grey20", face = "bold"
      )
    )

  return(p)
}


#' Decay curve plot (deprecated)
#'
#' Deprecated alias for \code{\link{plot_gradient_curve}}. Renamed because
#' "decay" is misleading -- many gradient curves are not monotonically
#' decaying (positive coefficients show the inverse pattern, and biphasic
#' curves are common). Calls through to \code{plot_gradient_curve()} with
#' the same arguments.
#'
#' @param ... Arguments passed to \code{\link{plot_gradient_curve}}.
#'
#' @return A \code{ggplot} object.
#'
#' @keywords internal
#' @export
plot_decay_curve <- function(...) {
  .Deprecated("plot_gradient_curve")
  plot_gradient_curve(...)
}


#' Proportion-expressing curve plot
#'
#' Specialised wrapper around \code{\link{plot_gradient_curve}} for the
#' canonical RIPPLE / HyMy plot: the **proportion of cells expressing a
#' gene** (counts > 0) as a function of distance from query cells. This is
#' the y-axis used in the original HyMy distance-correlation analysis and
#' is the most direct visual readout of a spatial gradient -- it does not
#' depend on cell-size normalisation, and the binomial scale (0 to 1)
#' makes the effect size easy to read.
#'
#' Use \code{\link{plot_gradient_curve}} when you want to plot something
#' other than proportion expressing (e.g. mean expression rate per UMI,
#' library-size-normalised mean, or a custom statistic).
#'
#' All other behaviour -- pooled vs per-sample mode, the bold mean overlay
#' with a 95\% CI ribbon, the subtitle with gradient score and FDR -- is
#' identical to \code{\link{plot_gradient_curve}}.
#'
#' \strong{CI interpretation differs between modes.} See the
#' "CI interpretation -- the two modes are not the same statistic"
#' section of \code{\link{plot_gradient_curve}} for the full discussion.
#' Briefly: pooled mode uses the supplied per-bin SE (binomial Wald,
#' \code{sqrt(p * (1 - p) / n)} on pooled-cell n, matching the original
#' HyMy plot); per-sample mode uses cross-sample
#' \code{sd / sqrt(n_samples)} (more honest about replicate
#' variability).
#'
#' \strong{Dot size encodes cell count per bin}, not significance. See
#' the "Dot size encoding" section of \code{\link{plot_gradient_curve}}
#' for the full discussion. Briefly: bigger dot = more cells contributed
#' to that bin's mean. The size legend is hidden by design.
#'
#' @inheritParams plot_gradient_curve
#' @param y_col Character. Name of the proportion-expressing column on
#'   \code{bin_stats}. Default \code{"prop_expressing"} (the column written
#'   by \code{\link{bin_decay_data}}).
#' @param y_lab Character or NULL. Y-axis label. Default
#'   \code{"Proportion expressing"}.
#' @param color Character. Line / ribbon / point color. Default
#'   \code{"#2C7FB8"} (a more "neutral" blue than the gradient-curve
#'   default red, since proportion-expressing plots usually depict the
#'   raw signal rather than a directional gradient).
#'
#' @return A \code{ggplot} object, or NULL if \code{bin_stats} is NULL or
#'   empty.
#'
#' @examples
#' \dontrun{
#' # Pooled mode: pre-binned proportion-expressing data
#' plot_prop_curve(
#'   bin_stats      = binned_prop,
#'   gene_name      = "Cxcl12",
#'   cell_type      = "T_cell",
#'   gradient_score = -0.005,
#'   fdr            = 1e-3
#' )
#'
#' # Per-sample mode: thin per-replicate lines + bold mean + 95% CI
#' plot_prop_curve(
#'   bin_stats   = binned_prop_per_sample,   # has a sample_id column
#'   gene_name   = "Ccr7",
#'   cell_type   = "T_cell",
#'   sample_col  = "sample_id",
#'   x_col       = "bin_mid_um",
#'   results     = merged_results
#' )
#' }
#'
#' @seealso \code{\link{plot_gradient_curve}} for the general-purpose
#'   version, and \code{\link{bin_decay_data}} for producing the
#'   \code{prop_expressing} column from raw counts and distances.
#'
#' @export
plot_prop_curve <- function(bin_stats, gene_name, cell_type,
                            gradient_score = NULL, fdr = NULL,
                            results = NULL,
                            sample_col = NULL,
                            x_col = "dist_mid",
                            y_col = "prop_expressing",
                            se_col = "se",
                            y_lab = "Proportion expressing",
                            query_label = "Query",
                            max_distance = 200,
                            color = "#2C7FB8",
                            sample_alpha = 0.4,
                            sample_linewidth = 0.4) {
  plot_gradient_curve(
    bin_stats        = bin_stats,
    gene_name        = gene_name,
    cell_type        = cell_type,
    gradient_score   = gradient_score,
    fdr              = fdr,
    results          = results,
    sample_col       = sample_col,
    x_col            = x_col,
    y_col            = y_col,
    se_col           = se_col,
    y_lab            = y_lab,
    query_label      = query_label,
    max_distance     = max_distance,
    color            = color,
    sample_alpha     = sample_alpha,
    sample_linewidth = sample_linewidth
  )
}


#' fGSEA dot plot
#'
#' Dotplot summary of fGSEA pathway enrichment across cell types. Each
#' point represents a (pathway, cell type) pair: dot size encodes
#' -log10(padj) and dot color encodes the normalized enrichment score (NES)
#' on a diverging red-blue scale. This is the companion visualization to
#' \code{plot_fgsea_heatmap()} / the heatmap in \code{run_ripple_atlas()};
#' the heatmap shows a filled NES grid, while the dotplot additionally
#' conveys significance.
#'
#' @param fgsea_results A \code{data.table} of fGSEA results, typically
#'   \code{fgsea_all_celltypes.csv} written by \code{run_ripple_atlas()}
#'   or the output of \code{run_ripple_fgsea()}. Must contain columns
#'   \code{cell_type}, \code{pathway_clean} (or \code{pathway}), \code{NES},
#'   and \code{padj}.
#' @param padj_threshold Numeric. Only pathways reaching this padj in at
#'   least one cell type are shown (default: 0.05). Set to 1 to show all
#'   pathways.
#' @param pathways Optional character vector to subset to specific pathways
#'   (matched against \code{pathway_clean} or \code{pathway}). If NULL, all
#'   pathways passing \code{padj_threshold} are shown.
#' @param top_n Integer or NULL. If set, limits to the top N pathways by
#'   mean absolute NES across cell types (useful when the heatmap has many
#'   rows).
#' @param title Character or NULL. Plot title (default: auto-generated).
#' @param subtitle Character or NULL. Plot subtitle (default: auto-generated).
#'
#' @return A \code{ggplot} object. Plot height scales with the number of
#'   pathways; a reasonable \code{ggsave} height is
#'   \code{max(5, n_pathways * 0.3 + 2)}.
#'
#' @examples
#' \dontrun{
#' fgsea_dt <- data.table::fread("ripple_atlas/fgsea_all_celltypes.csv")
#' plot_fgsea_dotplot(fgsea_dt)
#' plot_fgsea_dotplot(fgsea_dt, top_n = 20)
#' plot_fgsea_dotplot(fgsea_dt, pathways = c("HALLMARK_HYPOXIA",
#'                                           "HALLMARK_GLYCOLYSIS"))
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradientn
#'   scale_size_continuous labs theme_bw theme element_text
#' @importFrom data.table as.data.table copy
#' @export
plot_fgsea_dotplot <- function(fgsea_results,
                               padj_threshold = 0.05,
                               pathways = NULL,
                               top_n = NULL,
                               title = NULL,
                               subtitle = NULL) {
  dt <- data.table::as.data.table(data.table::copy(fgsea_results))

  if (!"pathway_clean" %in% names(dt)) {
    if ("pathway" %in% names(dt)) {
      dt[, pathway_clean := pathway]
    } else {
      stop("fgsea_results must contain a 'pathway_clean' or 'pathway' column.",
        call. = FALSE
      )
    }
  }

  required <- c("cell_type", "pathway_clean", "NES", "padj")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop("fgsea_results is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.null(pathways)) {
    dt <- dt[pathway_clean %in% pathways | pathway %in% pathways]
  } else {
    sig_pw <- dt[padj < padj_threshold, unique(pathway_clean)]
    dt <- dt[pathway_clean %in% sig_pw]
  }

  if (nrow(dt) == 0) {
    stop(
      "No pathways pass padj_threshold = ", padj_threshold,
      ". Try relaxing the threshold or passing a `pathways` vector.",
      call. = FALSE
    )
  }

  pw_order <- dt[, .(mean_abs_nes = mean(abs(NES), na.rm = TRUE)),
    by = pathway_clean
  ][order(-mean_abs_nes)]

  if (!is.null(top_n) && is.numeric(top_n) && top_n > 0) {
    keep <- pw_order$pathway_clean[seq_len(min(top_n, nrow(pw_order)))]
    dt <- dt[pathway_clean %in% keep]
    pw_order <- pw_order[pathway_clean %in% keep]
  }

  # Order: strongest mean absolute NES at the top
  dt[, pathway_clean := factor(pathway_clean, levels = rev(pw_order$pathway_clean))]

  dt[, neg_log10_padj := -log10(pmax(padj, 1e-20))]
  dt[, is_sig := padj < padj_threshold]

  nes_lim <- max(abs(dt$NES), na.rm = TRUE)
  if (!is.finite(nes_lim) || nes_lim == 0) nes_lim <- 1

  # Negative NES = induced near query (matches red in volcano);
  # positive NES = repressed near query (matches blue in volcano).
  diverging_palette <- grDevices::colorRampPalette(
    c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
      "#F7F7F7",
      "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC")
  )

  if (is.null(title)) title <- "Pathway enrichment per cell type"
  if (is.null(subtitle)) {
    subtitle <- sprintf(
      "Pathways with padj < %s in >= 1 cell type | dot size = -log10(padj) | red NES = induced, blue NES = repressed",
      format(padj_threshold)
    )
  }

  ggplot2::ggplot(
    dt,
    ggplot2::aes(x = .data$cell_type, y = .data$pathway_clean)
  ) +
    ggplot2::geom_point(ggplot2::aes(
      size = .data$neg_log10_padj, color = .data$NES,
      alpha = .data$is_sig
    )) +
    ggplot2::scale_color_gradientn(
      colors = diverging_palette(100),
      limits = c(-nes_lim, nes_lim),
      name = "NES"
    ) +
    ggplot2::scale_size_continuous(
      range = c(1.5, 6),
      name = expression(-log[10](padj))
    ) +
    ggplot2::scale_alpha_manual(
      values = c("TRUE" = 1, "FALSE" = 0.35),
      guide = "none"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = NULL, y = NULL
    ) +
    theme_ripple(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9, color = "grey40")
    )
}


#' Gene category dot plot
#'
#' Faceted dot plot of gradient coefficients for a curated gene panel
#' organized into biological categories. Each category occupies its own
#' horizontal strip (facet). Dot size encodes -log10(FDR), dot fill encodes
#' coefficient sign and magnitude on a diverging scale, and a black border
#' marks genes passing the FDR threshold.
#'
#' This is the visualization used for signature-level summaries
#' (e.g. CD8 exhaustion panel, macrophage polarization panel) in the RIPPLE
#' manuscript and companion analyses.
#'
#' @param results A \code{data.table} or \code{data.frame} of RIPPLE results
#'   containing, at minimum, \code{gene}, a coefficient column, and an FDR
#'   column. Typically the output of \code{merge_ripple_results()} filtered
#'   to a single cell type.
#' @param gene_categories Named list mapping category label to a character
#'   vector of gene symbols. Category names become facet strip labels; their
#'   order in the list becomes the facet order (top to bottom) unless
#'   \code{category_order} is supplied. Genes appearing in multiple categories
#'   are assigned to the first match.
#' @param coef_col Character. Coefficient column (default: \code{"median_coef"}).
#' @param fdr_col Character. FDR column (default: \code{"fisher_fdr"}).
#' @param fdr_threshold Numeric. Significance threshold used for the black
#'   border (default: 0.05).
#' @param query_label Character. Display label for the query cell type,
#'   used in subtitle (default: \code{"Query"}).
#' @param category_order Character vector or NULL. Explicit category order
#'   (top to bottom). If NULL, uses the order of \code{gene_categories}.
#' @param max_neg_log10_fdr Numeric. Cap for the -log10(FDR) size aesthetic
#'   so a single extreme gene does not dominate the size scale
#'   (default: 20).
#' @param title Character or NULL. Plot title. If NULL, auto-generated.
#' @param subtitle Character or NULL. Plot subtitle. If NULL, auto-generated.
#' @param size_breaks Numeric vector. Legend breaks for the dot-size scale
#'   (default: \code{c(2, 5, 10, 20)}).
#'
#' @return A \code{ggplot} object. Height scales with the number of genes;
#'   a reasonable \code{ggsave} height is \code{max(5, n_genes * 0.28 + 2.5)}.
#'
#' @examples
#' \dontrun{
#' exhaustion_panel <- list(
#'   "Inhibitory receptors" = c("Pdcd1", "Lag3", "Havcr2", "Tigit", "Ctla4"),
#'   "Transcription factors" = c("Tox", "Eomes", "Tbx21", "Tcf7"),
#'   "Effector molecules" = c("Gzmb", "Gzma", "Prf1", "Ifng", "Tnf"),
#'   "Tpex markers" = c("Slamf6", "Xcl1", "Il7r", "Cxcr5")
#' )
#' cd8 <- all_results[cell_type == "CD8_T_cells"]
#' plot_gene_category_dotplot(cd8, exhaustion_panel, query_label = "Tumor")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_vline geom_point scale_fill_gradientn
#'   scale_color_manual scale_size_continuous facet_grid labs theme_bw theme
#'   element_text element_rect element_blank unit
#' @importFrom data.table as.data.table copy setorderv
#' @export
plot_gene_category_dotplot <- function(results,
                                       gene_categories,
                                       coef_col = "median_coef",
                                       fdr_col = "fisher_fdr",
                                       fdr_threshold = 0.05,
                                       query_label = "Query",
                                       category_order = NULL,
                                       max_neg_log10_fdr = 20,
                                       title = NULL,
                                       subtitle = NULL,
                                       size_breaks = c(2, 5, 10, 20)) {
  if (!is.list(gene_categories) || is.null(names(gene_categories)) ||
    any(names(gene_categories) == "")) {
    stop("`gene_categories` must be a named list of character vectors.",
      call. = FALSE
    )
  }

  plot_data <- data.table::as.data.table(data.table::copy(results))
  for (col in c("gene", coef_col, fdr_col)) {
    if (!col %in% names(plot_data)) {
      stop("Required column not found in results: ", col, call. = FALSE)
    }
  }

  gene_to_cat <- character()
  for (cat_name in names(gene_categories)) {
    new_genes <- setdiff(gene_categories[[cat_name]], names(gene_to_cat))
    gene_to_cat[new_genes] <- cat_name
  }

  plot_data <- plot_data[plot_data$gene %in% names(gene_to_cat)]
  if (nrow(plot_data) == 0) {
    stop("None of the genes in `gene_categories` were found in `results`.",
      call. = FALSE
    )
  }

  plot_data[, category := gene_to_cat[as.character(gene)]]

  if (is.null(category_order)) {
    category_order <- names(gene_categories)
  }
  category_order <- intersect(category_order, unique(plot_data$category))
  plot_data[, category := factor(category, levels = category_order)]

  plot_data[, neg_log10_fdr := -log10(pmax(get(fdr_col), 1e-300))]
  plot_data[, neg_log10_fdr_capped := pmin(neg_log10_fdr, max_neg_log10_fdr)]
  plot_data[, is_sig := get(fdr_col) < fdr_threshold]

  data.table::setorderv(plot_data, c("category", coef_col))
  plot_data[, gene := factor(gene, levels = unique(gene))]

  coef_lim <- max(abs(plot_data[[coef_col]]), na.rm = TRUE)
  if (!is.finite(coef_lim) || coef_lim == 0) coef_lim <- 1e-4

  diverging_palette <- grDevices::colorRampPalette(
    c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
      "#F7F7F7",
      "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC")
  )

  if (is.null(title)) {
    title <- "Gene category gradient dot plot"
  }
  if (is.null(subtitle)) {
    subtitle <- paste0(
      "Negative gradient score = higher expression rate near ", query_label,
      " | Black border = FDR < ", format(fdr_threshold, nsmall = 0)
    )
  }

  ggplot2::ggplot(plot_data, ggplot2::aes(
    x = .data[[coef_col]], y = .data$gene
  )) +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed", color = "grey50"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = .data$neg_log10_fdr_capped,
        fill = .data[[coef_col]],
        color = .data$is_sig
      ),
      shape = 21, stroke = 0.6
    ) +
    ggplot2::scale_fill_gradientn(
      colors = diverging_palette(100),
      limits = c(-coef_lim, coef_lim),
      name = "Gradient\nscore"
    ) +
    ggplot2::scale_color_manual(
      values = c("TRUE" = "black", "FALSE" = "grey70"),
      guide = "none"
    ) +
    ggplot2::scale_size_continuous(
      range = c(1.5, 6),
      name = expression(-log[10](FDR)),
      breaks = size_breaks
    ) +
    ggplot2::facet_grid(category ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Gradient score (per um)",
      y = NULL
    ) +
    theme_ripple(base_size = 11) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9, face = "italic"),
      strip.text.y = ggplot2::element_text(size = 8, face = "bold", angle = 0),
      strip.background = ggplot2::element_rect(fill = "grey95"),
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 8, color = "grey40"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.key.size = ggplot2::unit(0.4, "cm")
    )
}


#' Ensure directory exists
#'
#' Creates a directory (and all parent directories) if it does not already exist.
#'
#' @param path Character. Directory path to create.
#' @return The path, invisibly.
#' @noRd
ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  return(invisible(path))
}


#' Save data.table to CSV
#'
#' Saves a \code{data.table} to a CSV file, creating the output directory
#' if needed.
#'
#' @param dt A \code{data.table} to save.
#' @param output_dir Character. Directory to save into.
#' @param filename Character. Filename (without extension; \code{.csv} is appended).
#' @return Invisible NULL.
#'
#' @importFrom data.table fwrite
#' @noRd
save_results <- function(dt, output_dir, filename) {
  ensure_dir(output_dir)
  filepath <- file.path(output_dir, paste0(filename, ".csv"))
  data.table::fwrite(dt, filepath)
  message("Saved: ", filepath)
  invisible(NULL)
}


#' Save ggplot to file
#'
#' Saves a ggplot object to one or more file formats, creating the output
#' directory if needed.
#'
#' @param p A \code{ggplot} object.
#' @param output_dir Character. Directory to save into.
#' @param filename Character. Filename without extension.
#' @param width Numeric. Plot width in inches (default: 8).
#' @param height Numeric. Plot height in inches (default: 6).
#' @param formats Character vector of file formats (default: \code{c("png", "pdf")}).
#' @return Invisible NULL.
#'
#' @importFrom ggplot2 ggsave
#' @noRd
save_plot <- function(p, output_dir, filename, width = 8, height = 6,
                      formats = c("png", "pdf")) {
  ensure_dir(output_dir)

  for (fmt in formats) {
    filepath <- file.path(output_dir, paste0(filename, ".", fmt))
    ggplot2::ggsave(filepath, p, width = width, height = height, dpi = 300)
    message("Saved: ", filepath)
  }
  invisible(NULL)
}


#' Create gradient volcano plot for distance-expression analysis
#'
#' Creates a volcano plot showing genes colored by decay pattern and sized by
#' significance. Uses \code{gradient_score} on the x-axis and
#' \code{-log10(FDR)} on the y-axis, with the top significant genes labeled.
#'
#' @param results A \code{data.table} with per-gene meta-analysis results.
#'   Must contain \code{gradient_score}, \code{gene}, and an FDR column
#'   (\code{fisher_fdr} or \code{fdr}). Optionally contains
#'   \code{decay_pattern} for coloring.
#' @param cell_type Character. Cell type name for the plot title.
#' @param output_path Character. File path to save the plot (e.g., a PDF).
#' @param fdr_threshold Numeric. Significance threshold (default: 0.05).
#' @param query_label Character. Display label for the query cell type
#'   (default: "Query").
#' @param k_neighbors Integer. Number of nearest neighbors used in the
#'   analysis, shown in the title (default: 1).
#'
#' @return The \code{ggplot} object, invisibly.
#'
#' @examples
#' \dontrun{
#' create_gradient_volcano(
#'   results = meta_results,
#'   cell_type = "LEC",
#'   output_path = "volcano.pdf",
#'   query_label = "Tumor"
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline xlim labs
#'   theme_classic theme element_text scale_size_manual scale_color_brewer
#'   ggsave
#' @importFrom ggrepel geom_text_repel
#' @importFrom data.table copy
#' @export
create_gradient_volcano <- function(results, cell_type, output_path,
                                    fdr_threshold = 0.05,
                                    query_label = "Query",
                                    k_neighbors = 1) {
  plot_data <- data.table::copy(results)

  # Use fisher_fdr if available, otherwise fall back to fdr
  fdr_col <- if ("fisher_fdr" %in% names(plot_data)) "fisher_fdr" else "fdr"
  plot_data[, neg_log10_fdr := -log10(get(fdr_col))]
  plot_data[neg_log10_fdr > 50, neg_log10_fdr := 50]
  plot_data[, significant := get(fdr_col) < fdr_threshold]

  top_genes <- head(plot_data[significant == TRUE][order(get(fdr_col))], 20)
  max_score <- max(abs(plot_data$gradient_score), na.rm = TRUE) * 1.1

  if ("decay_pattern" %in% names(plot_data)) {
    plot_data[, decay_pattern := factor(decay_pattern,
      levels = c(
        "linear", "exponential", "step_10um", "step_25um",
        "step_50um", "none", "no_variation", "insufficient_data",
        "not_significant", "undetermined"
      )
    )]
  }

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = gradient_score, y = neg_log10_fdr)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        color = if ("decay_pattern" %in% names(plot_data)) {
          decay_pattern
        } else {
          NULL
        },
        size = significant
      ),
      alpha = 0.6
    ) +
    ggplot2::scale_size_manual(
      values = c("FALSE" = 1, "TRUE" = 2.5),
      guide = "none"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(fdr_threshold),
      linetype = "dashed", color = "grey40"
    ) +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "solid",
      color = "grey60"
    ) +
    ggplot2::xlim(-max_score, max_score) +
    ggrepel::geom_text_repel(
      data = top_genes,
      ggplot2::aes(label = gene),
      size = 3, max.overlaps = 20, box.padding = 0.5
    ) +
    ggplot2::labs(
      title = sprintf("%s near %s", cell_type, query_label),
      subtitle = sprintf(
        "Distance-Expression gradients (Poisson GLM, k=%d) | %d genes significant (FDR < %.2f)",
        k_neighbors, sum(plot_data$significant, na.rm = TRUE),
        fdr_threshold
      ),
      x = paste0(
        "Gradient score (negative = ", query_label,
        "-induced)"
      ),
      y = paste0("-log10(", fdr_col, ")")
    ) +
    theme_ripple(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  if ("decay_pattern" %in% names(plot_data)) {
    p <- p + ggplot2::scale_color_brewer(
      palette = "Set2",
      name = "Decay Pattern"
    )
  }

  ggplot2::ggsave(output_path, p, width = 10, height = 8)
  invisible(p)
}


#' Forest plot of per-sample gradient coefficients
#'
#' Shows per-sample Poisson GLM coefficients with 95\% confidence intervals
#' for one or many genes. Useful for assessing cross-sample reproducibility
#' of distance-expression gradients.
#'
#' Two modes:
#' \itemize{
#'   \item \strong{Single-gene mode} (legacy): pass \code{coefs}, \code{ses},
#'     and \code{sample_ids} as parallel vectors plus a scalar \code{gene}.
#'     The plot adds an inverse-variance weighted combined estimate (diamond)
#'     row computed from the supplied per-sample coefficients.
#'   \item \strong{Multi-gene mode}: pass a long-format \code{data.frame} or
#'     \code{data.table} as \code{coefs} containing the columns
#'     \code{gene}, \code{sample_id}, \code{coef}, \code{se}. Each gene
#'     becomes a y-axis row; per-sample coefficients are dodged horizontally
#'     and coloured by sample. No combined diamond is drawn (the panel is
#'     about replicate consistency across a curated gene panel, not summary
#'     statistics for a single gene).
#' }
#'
#' @param coefs In single-gene mode, a numeric vector of per-sample
#'   coefficients. In multi-gene mode, a long-format data.frame /
#'   data.table with columns \code{gene}, \code{sample_id}, \code{coef},
#'   \code{se}.
#' @param ses Numeric vector of per-sample standard errors. Required in
#'   single-gene mode; ignored in multi-gene mode.
#' @param sample_ids Character vector of sample identifiers, same length as
#'   \code{coefs}. Required in single-gene mode; ignored in multi-gene mode.
#' @param gene Character. Single gene name (single-gene mode), or character
#'   vector of genes to display (multi-gene mode). In multi-gene mode the
#'   vector also controls the y-axis order (top-to-bottom = first-to-last).
#'   If omitted in multi-gene mode, all genes in \code{coefs} are used in
#'   the order they first appear.
#' @param cell_type Character. Cell type name for the plot title.
#' @param output_path Optional character. If supplied, save the plot to
#'   this path with \code{ggplot2::ggsave()}. Default \code{NULL} returns
#'   the ggplot without writing to disk.
#' @param query_label Character. Display label for the query cell type
#'   (default: "Query").
#' @param show_combined Logical. Single-gene mode only. If \code{TRUE}
#'   (default), append a "Median" summary row drawn as a larger diamond at
#'   the median of the per-sample coefficients. The subtitle reports the
#'   sign-consistency count across samples. Set \code{FALSE} when the panel
#'   should display only per-sample replicate consistency without a
#'   summary marker. Ignored in multi-gene mode (no summary row is ever
#'   drawn there).
#'
#' @return A \code{ggplot} object (returned visibly when \code{output_path}
#'   is \code{NULL}; invisibly otherwise). Returns invisible \code{NULL}
#'   in single-gene mode if fewer than two valid samples.
#'
#' @examples
#' \dontrun{
#' # Single-gene mode (legacy):
#' create_forest_plot(
#'   coefs = c(-0.005, -0.003, -0.004),
#'   ses = c(0.001, 0.002, 0.001),
#'   sample_ids = c("S1", "S2", "S3"),
#'   gene = "Cxcl12",
#'   cell_type = "LEC",
#'   output_path = "Cxcl12_forest.pdf"
#' )
#'
#' # Single-gene mode without the combined diamond (replicate-consistency
#' # focus, no meta-analytic summary):
#' create_forest_plot(
#'   coefs = c(-0.005, -0.003, -0.004),
#'   ses = c(0.001, 0.002, 0.001),
#'   sample_ids = c("S1", "S2", "S3"),
#'   gene = "IGFBP5", cell_type = "fibroblast",
#'   show_combined = FALSE
#' )
#'
#' # Multi-gene mode (per-sample table):
#' coefs_dt <- data.table::data.table(
#'   gene = rep(c("Ccr7", "Sell", "Lef1"), each = 4),
#'   sample_id = rep(paste0("M", 1:4), 3),
#'   coef = rnorm(12, -0.005, 0.002),
#'   se = runif(12, 0.001, 0.003)
#' )
#' create_forest_plot(coefs_dt,
#'                    gene = c("Ccr7", "Sell", "Lef1"),
#'                    cell_type = "T_cell_all")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_vline geom_errorbar geom_point
#'   scale_shape_manual scale_size_manual scale_colour_brewer
#'   position_dodge labs theme element_text ggsave
#' @importFrom data.table data.table as.data.table
#' @export
create_forest_plot <- function(coefs, ses = NULL, sample_ids = NULL,
                               gene, cell_type,
                               output_path = NULL, query_label = "Query",
                               show_combined = TRUE) {
  if (is.data.frame(coefs)) {
    return(.forest_multi_gene(coefs, gene, cell_type, output_path,
                              query_label))
  }

  valid_idx <- !is.na(coefs) & !is.na(ses) & ses > 0
  if (sum(valid_idx) < 2) {
    return(invisible(NULL))
  }

  coefs <- coefs[valid_idx]
  ses <- ses[valid_idx]
  sample_ids <- sample_ids[valid_idx]

  plot_data <- data.table::data.table(
    sample = sample_ids,
    coef = coefs,
    se = ses,
    lower = coefs - 1.96 * ses,
    upper = coefs + 1.96 * ses
  )

  if (isTRUE(show_combined)) {
    median_beta <- stats::median(coefs)
    meta_row <- data.table::data.table(
      sample = "Median",
      coef = median_beta,
      se = NA_real_,
      lower = NA_real_,
      upper = NA_real_
    )
    plot_data <- rbind(plot_data, meta_row)
    plot_data[, is_combined := sample == "Median"]
  } else {
    plot_data[, is_combined := FALSE]
  }

  plot_data[, sample := factor(sample, levels = rev(unique(sample)))]

  n_neg <- sum(coefs < 0)
  n_pos <- sum(coefs > 0)
  sign_text <- sprintf(
    "%d/%d samples agree on sign",
    max(n_neg, n_pos), length(coefs)
  )
  subtitle_text <- sign_text

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = coef, y = sample)) +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed",
      color = "grey50"
    ) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lower, xmax = upper),
      orientation = "y", width = 0.2
    ) +
    ggplot2::geom_point(ggplot2::aes(
      shape = is_combined,
      size = is_combined
    )) +
    ggplot2::scale_shape_manual(
      values = c("FALSE" = 16, "TRUE" = 18),
      guide = "none"
    ) +
    ggplot2::scale_size_manual(
      values = c("FALSE" = 3, "TRUE" = 5),
      guide = "none"
    ) +
    ggplot2::labs(
      x = "Gradient score (per um)",
      y = NULL,
      title = sprintf("%s in %s (Poisson GLM)", gene, cell_type),
      subtitle = subtitle_text
    ) +
    theme_ripple(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  if (!is.null(output_path)) {
    ggplot2::ggsave(output_path, p, width = 6, height = 4)
    return(invisible(p))
  }
  p
}

# Internal -- multi-gene forest panel.
# Each gene becomes a y-axis row, per-sample coefs dodged and coloured by
# sample. No combined diamond (the panel is about replicate consistency).
.forest_multi_gene <- function(coefs, gene, cell_type, output_path,
                               query_label) {
  dt <- data.table::as.data.table(coefs)
  required <- c("gene", "sample_id", "coef", "se")
  miss <- setdiff(required, names(dt))
  if (length(miss) > 0) {
    stop(
      "Multi-gene mode requires columns: ", paste(required, collapse = ", "),
      ". Missing: ", paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  if (!missing(gene) && !is.null(gene) && length(gene) > 0) {
    requested <- as.character(gene)
    available <- intersect(requested, dt$gene)
    if (length(available) == 0) {
      stop("None of the requested genes are present in `coefs`.",
           call. = FALSE)
    }
    missing_genes <- setdiff(requested, available)
    if (length(missing_genes) > 0) {
      message(
        sprintf("create_forest_plot(): %d of %d genes not present: %s",
                length(missing_genes), length(requested),
                paste(utils::head(missing_genes, 10), collapse = ", "))
      )
    }
    dt <- dt[gene %in% available]
    gene_levels <- rev(available)
  } else {
    gene_levels <- rev(unique(dt$gene))
  }
  dt <- dt[!is.na(coef) & !is.na(se) & se > 0]
  if (nrow(dt) < 2) {
    return(invisible(NULL))
  }
  dt[, ci_lo := coef - 1.96 * se]
  dt[, ci_hi := coef + 1.96 * se]
  dt[, gene := factor(gene, levels = gene_levels)]

  p <- ggplot2::ggplot(dt, ggplot2::aes(x = coef, y = gene,
                                         colour = sample_id)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey50") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
      orientation = "y", width = 0.3, alpha = 0.7,
      position = ggplot2::position_dodge(width = 0.6)
    ) +
    ggplot2::geom_point(
      size = 2.2, alpha = 0.9,
      position = ggplot2::position_dodge(width = 0.6)
    ) +
    ggplot2::scale_colour_brewer(palette = "Set2", name = "Sample") +
    ggplot2::labs(
      x = "Per-sample gradient score (per um)",
      y = NULL,
      title = sprintf("Per-sample replicate consistency in %s", cell_type)
    ) +
    theme_ripple(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  if (!is.null(output_path)) {
    ggplot2::ggsave(output_path, p,
                    width = 7, height = max(3, length(gene_levels) * 0.6))
    return(invisible(p))
  }
  p
}


#' Create coefficient strip plot for per-sample reproducibility
#'
#' Shows per-sample coefficients for the top significant genes, allowing
#' visual assessment of effect consistency across biological replicates.
#' Genes are ordered by combined coefficient magnitude.
#'
#' @param coef_results A \code{data.table} with per-sample coefficient
#'   results. Must contain columns \code{gene}, \code{coef}, \code{se},
#'   and \code{sample_id}.
#' @param meta_results A \code{data.table} with per-gene meta-analysis
#'   results. Must contain \code{gene}, \code{combined_coef}, and an FDR
#'   column (\code{fisher_fdr} or \code{fdr}).
#' @param cell_type Character. Cell type name for the plot title.
#' @param output_path Character. File path to save the plot.
#' @param fdr_threshold Numeric. Significance threshold (default: 0.05).
#' @param n_top Integer. Maximum number of top genes to plot (default: 20).
#'
#' @return The \code{ggplot} object invisibly, or invisible \code{NULL} if
#'   no significant genes.
#'
#' @examples
#' \dontrun{
#' create_coefficient_strips(
#'   coef_results = per_sample_coefs,
#'   meta_results = meta_results,
#'   cell_type = "LEC",
#'   output_path = "coefficient_strips.pdf"
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_vline geom_errorbarh geom_point
#'   scale_color_brewer labs theme_bw theme element_text ggsave
#' @export
create_coefficient_strips <- function(coef_results, meta_results, cell_type,
                                      output_path, fdr_threshold = 0.05,
                                      n_top = 20) {
  fdr_col <- if ("fisher_fdr" %in% names(meta_results)) {
    "fisher_fdr"
  } else {
    "fdr"
  }

  sig_genes <- meta_results[get(fdr_col) < fdr_threshold][
    order(get(fdr_col))
  ]$gene
  if (length(sig_genes) == 0) {
    return(invisible(NULL))
  }

  genes_to_plot <- head(sig_genes, n_top)
  plot_data <- coef_results[gene %in% genes_to_plot & !is.na(coef)]
  if (nrow(plot_data) == 0) {
    return(invisible(NULL))
  }

  order_col <- if ("median_coef" %in% names(meta_results)) {
    "median_coef"
  } else {
    "combined_coef"  # legacy CSV compatibility
  }
  gene_order <- meta_results[gene %in% genes_to_plot][
    order(get(order_col))
  ]$gene
  plot_data[, gene := factor(gene, levels = gene_order)]

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = coef, y = gene,
      color = sample_id
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed",
      color = "grey50"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = coef - 1.96 * se, xmax = coef + 1.96 * se),
      orientation = "y", width = 0.2, alpha = 0.5
    ) +
    ggplot2::geom_point(size = 2.5, alpha = 0.8) +
    ggplot2::scale_color_brewer(palette = "Set1", name = "Sample") +
    ggplot2::labs(
      x = "Gradient score (per um)",
      y = NULL,
      title = sprintf(
        "Per-Sample Gradient Scores: %s (Poisson GLM)",
        cell_type
      ),
      subtitle = sprintf(
        "Top %d significant genes | Points = individual samples",
        length(genes_to_plot)
      )
    ) +
    theme_ripple(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggplot2::ggsave(output_path, p,
    width = 8, height = 0.4 * length(genes_to_plot) + 2
  )
  invisible(p)
}


#' Plot k-selection diagnostics for choosing k_neighbors
#'
#' Generates two diagnostic plots to help users choose an appropriate k value
#' for the \code{k_neighbors} parameter in \code{\link{run_ripple}}:
#' \enumerate{
#'   \item Mean distance to the k nearest query cells vs k, per target cell type
#'   \item Standard deviation of that distance vs k, per target cell type
#' }
#'
#' At low k, the distance to the nearest query cell is noisy (sensitive to
#' one stray query cell). As k increases, the mean distance stabilizes. The
#' "elbow" where the SD flattens suggests a reasonable k choice. Cell types
#' closer to the query population will have lower mean distances at all k.
#'
#' @param input A Seurat, SCE, or SpatialExperiment object, or a path to an
#'   \code{.rds} file.
#' @param query_celltype Character. Query cell type label.
#' @param celltype_column Character. Cell type column name.
#' @param sample_column Character. Sample column (default: \code{"sample_id"}).
#'   Diagnostics are computed per-sample then averaged.
#' @param k_range Integer vector of k values to evaluate (default: \code{1:10}).
#' @param max_distance_um Numeric. Distance cap (default: 200).
#' @param x_column Character or NULL. X coordinate column (default: auto-detect).
#' @param y_column Character or NULL. Y coordinate column (default: auto-detect).
#' @param verbose Logical. Print progress (default: TRUE).
#'
#' @return A \code{patchwork} object with two panels. The underlying summary
#'   data is returned invisibly as a \code{data.table} with columns: k,
#'   cell_type, mean_dist, sd_dist.
#'
#' @examples
#' \dontrun{
#' data(ripple_mock_data)
#' plot_k_diagnostics(
#'   ripple_mock_data,
#'   query_celltype = "Tumor",
#'   celltype_column = "cell_type"
#' )
#' }
#'
#' @importFrom data.table data.table rbindlist
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_bw
#'   scale_x_continuous
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom RANN nn2
#' @export
plot_k_diagnostics <- function(input,
                               query_celltype,
                               celltype_column,
                               sample_column = "sample_id",
                               k_range = 1:10,
                               max_distance_um = 200,
                               x_column = NULL,
                               y_column = NULL,
                               verbose = TRUE) {
  .msg <- function(...) if (isTRUE(verbose)) message(...)

  # Load data
  data <- .resolve_input(input, require_expr = FALSE, verbose = verbose)
  cell_data <- data$meta
  rm(data)

  coord_cols <- get_coord_columns(cell_data, x_col = x_column, y_col = y_column)
  coords <- as.matrix(cell_data[, ..coord_cols])

  if (!celltype_column %in% names(cell_data)) {
    stop("celltype_column '", celltype_column, "' not found.", call. = FALSE)
  }

  celltypes_all <- cell_data[[celltype_column]]
  query_mask <- !is.na(celltypes_all) & celltypes_all == query_celltype
  if (sum(query_mask) < 1) {
    stop("No query cells found for '", query_celltype, "'.", call. = FALSE)
  }

  # Target cell types (everything except query)
  all_types <- unique(celltypes_all)
  target_types <- setdiff(all_types, c(query_celltype, NA_character_))

  max_k <- max(k_range)

  # Distances MUST be measured within each sample's own coordinate frame. A
  # single pooled nn2 can match a target cell to a query cell in a different
  # sample whenever coordinate frames overlap (e.g. per-slide/per-patient
  # coordinates that both start near the origin), producing meaningless
  # distances -- exactly the failure this diagnostic is meant to catch.
  use_samples <- !is.null(sample_column) && sample_column %in% names(cell_data)
  if (!use_samples) {
    warning("sample_column '", sample_column, "' not found; computing k ",
      "diagnostics pooled across all cells. If coordinate frames differ per ",
      "sample, the resulting distances are not meaningful.", call. = FALSE)
    sample_vec <- rep("__all__", nrow(cell_data))
  } else {
    sample_vec <- as.character(cell_data[[sample_column]])
  }
  samples <- unique(sample_vec[!is.na(sample_vec)])

  # Cap k at the smallest per-sample query count so every sample supports
  # every k in the range.
  query_per_sample <- vapply(samples, function(s) {
    sum(query_mask & sample_vec == s)
  }, integer(1))
  samples_with_query <- samples[query_per_sample >= 1]
  if (length(samples_with_query) < length(samples)) {
    warning(length(samples) - length(samples_with_query),
      " sample(s) have no '", query_celltype, "' cells; their cells are ",
      "excluded from the k diagnostics.", call. = FALSE)
  }
  min_query <- min(query_per_sample[query_per_sample >= 1])
  effective_k <- min(max_k, min_query)
  if (effective_k < max_k) {
    .msg("  Smallest per-sample query count is ", min_query,
      "; capping k_range at ", effective_k)
    k_range <- k_range[k_range <= effective_k]
  }

  .msg("Computing per-sample distances for k = 1 to ", effective_k, "...")

  # Fill an (n_cells x effective_k) distance matrix, one nn2 call per sample.
  nn_dists <- matrix(NA_real_, nrow = nrow(cell_data), ncol = effective_k)
  for (s in samples_with_query) {
    cell_idx <- which(sample_vec == s)
    q_idx <- which(query_mask & sample_vec == s)
    nn_s <- RANN::nn2(
      coords[q_idx, , drop = FALSE],
      coords[cell_idx, , drop = FALSE],
      k = effective_k
    )
    nn_dists[cell_idx, ] <- nn_s$nn.dists
  }

  # Compute per-k, per-cell-type stats
  .msg("Computing per-cell-type distance statistics...")
  summary_list <- lapply(k_range, function(k) {
    if (k == 1) {
      dists <- nn_dists[, 1]
    } else {
      dists <- rowMeans(nn_dists[, 1:k, drop = FALSE])
    }
    dists <- pmin(dists, max_distance_um)

    data.table::rbindlist(lapply(target_types, function(ct) {
      ct_mask <- !is.na(celltypes_all) & celltypes_all == ct
      ct_dists <- dists[ct_mask]
      data.table::data.table(
        k = k,
        cell_type = ct,
        mean_dist = mean(ct_dists, na.rm = TRUE),
        sd_dist = stats::sd(ct_dists, na.rm = TRUE),
        n_cells = sum(ct_mask & !is.na(dists))
      )
    }))
  })

  summary_dt <- data.table::rbindlist(summary_list)

  # Plot 1: Mean distance vs k
  p1 <- ggplot2::ggplot(
    summary_dt,
    ggplot2::aes(
      x = .data$k, y = .data$mean_dist,
      color = .data$cell_type
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = k_range) +
    ggplot2::labs(
      title = "Mean distance to query vs k",
      subtitle = "Steep rise = query cells are clustered, use low k; flat = k has little effect",
      x = "k (number of nearest query cells)",
      y = "Mean distance (um)",
      color = "Cell type"
    ) +
    theme_ripple(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")

  # Plot 2: SD of distance vs k
  p2 <- ggplot2::ggplot(
    summary_dt,
    ggplot2::aes(
      x = .data$k, y = .data$sd_dist,
      color = .data$cell_type
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = k_range) +
    ggplot2::labs(
      title = "SD of distance to query vs k",
      subtitle = "Flat = k=1 is sufficient; sharp drop = higher k stabilizes the estimate",
      x = "k (number of nearest query cells)",
      y = "SD of distance (um)",
      color = "Cell type"
    ) +
    theme_ripple(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")

  combined <- patchwork::wrap_plots(p1, p2, ncol = 2) +
    patchwork::plot_annotation(
      title = paste0("k-selection diagnostics (query: ", query_celltype, ")"),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 13)
      )
    )

  print(combined)
  invisible(summary_dt)
}


#' Stacked bar of gene counts per cell type, split by specificity
#'
#' Builds a stacked bar chart of significant gene counts per target cell type,
#' split by direction (induced vs repressed) and by specificity class
#' (specific + moderate = cell-type-restricted signal; broad =
#' broadly-expressed, the heuristic flag for potential ambient RNA / segmentation
#' artefacts).
#'
#' Matches the stacked barplot used in the original HyMy / TDLN analysis and
#' extends the plain induced/repressed bar chart by surfacing the
#' broad-class portion of each bar.
#'
#' @param results \code{data.table} or data.frame with columns \code{gene},
#'   \code{cell_type}, a significance column, and a coefficient column.
#'   Typically the output of \code{\link{run_ripple}} or
#'   \code{\link{merge_ripple_results}} (i.e. \code{all_genes_results.csv}).
#' @param fdr_col Character. Significance column name (default:
#'   \code{"fisher_fdr"}).
#' @param coef_col Character. Coefficient column name (default:
#'   \code{"median_coef"}).
#' @param fdr_threshold Numeric. FDR cutoff (default: \code{0.05}).
#' @param broad_threshold Integer. Number of cell types at or above
#'   which a gene is flagged as "broad" (default: \code{4}).
#'   The default is illustrative -- adapt to your panel size and
#'   annotation granularity. See "Choosing \code{broad_threshold}" in
#'   \code{\link{classify_gene_specificity}} for guidance.
#' @param query_label Character. Display label for the query cell type
#'   (default: \code{"Query"}). Used in legend entries like
#'   "Query-induced (specific)".
#' @param cell_type_order Character vector or \code{NULL}. Optional ordering
#'   of cell types along the x axis; \code{NULL} means order by total
#'   significant genes, most to fewest.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \dontrun{
#' results <- fread("all_genes_results.csv")
#' plot_gene_counts_by_celltype(results, query_label = "Tumor")
#' }
#'
#' @importFrom data.table copy fifelse melt
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual labs theme_bw
#'   theme element_text
#' @export
plot_gene_counts_by_celltype <- function(results,
                                         fdr_col = "fisher_fdr",
                                         coef_col = "median_coef",
                                         fdr_threshold = 0.05,
                                         broad_threshold = 4,
                                         query_label = "Query",
                                         cell_type_order = NULL) {
  if (!inherits(results, "data.table")) {
    results <- data.table::as.data.table(results)
  }

  spec <- classify_gene_specificity(
    results,
    fdr_col         = fdr_col,
    fdr_threshold   = fdr_threshold,
    broad_threshold = broad_threshold
  )
  if (nrow(spec) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = "No significant genes at the chosen threshold") +
        ggplot2::theme_void()
    )
  }

  sig <- results[!is.na(get(fdr_col)) & get(fdr_col) < fdr_threshold]
  sig <- merge(sig, spec[, .(gene, specificity_class)], by = "gene",
               all.x = TRUE)
  sig[, direction := data.table::fifelse(
    get(coef_col) < 0,
    paste0(query_label, "-induced"),
    paste0(query_label, "-repressed")
  )]
  sig[, type := data.table::fifelse(
    specificity_class %in% c("specific", "moderate"),
    "specific", "broad"
  )]

  count_by_dir <- sig[, .N, by = .(cell_type, direction, type)]
  setnames(count_by_dir, "N", "count")

  if (is.null(cell_type_order)) {
    ct_totals <- sig[, .N, by = cell_type][order(-N), cell_type]
    cell_type_order <- as.character(ct_totals)
  }
  count_by_dir[, cell_type := factor(cell_type, levels = cell_type_order)]

  induced_label <- paste0(query_label, "-induced")
  repressed_label <- paste0(query_label, "-repressed")

  fill_keys <- c(
    paste0("specific.", induced_label),
    paste0("broad.",    induced_label),
    paste0("specific.", repressed_label),
    paste0("broad.",    repressed_label)
  )
  fill_colours <- stats::setNames(
    c("#B2182B", "#F4A582", "#2166AC", "#92C5DE"),
    fill_keys
  )
  fill_labels <- stats::setNames(
    c(paste0(induced_label, " (specific)"),
      paste0(induced_label, " (broad)"),
      paste0(repressed_label, " (specific)"),
      paste0(repressed_label, " (broad)")),
    fill_keys
  )

  ggplot2::ggplot(
    count_by_dir,
    ggplot2::aes(x = cell_type, y = count,
                 fill = interaction(type, direction))
  ) +
    ggplot2::geom_col(position = "stack", width = 0.7) +
    ggplot2::scale_fill_manual(values = fill_colours,
                               labels = fill_labels, name = NULL) +
    ggplot2::labs(
      x = NULL, y = "Significant genes",
      title = "Gradient genes per cell type",
      subtitle = sprintf(
        "FDR < %g | faded = broad-expression class (significant in >= %d cell types)",
        fdr_threshold, broad_threshold)
    ) +
    theme_ripple(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "grey40"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 8)
    )
}


#' Summary bar of gene counts by specificity class
#'
#' Produces a compact bar chart of significant genes, one bar per specificity
#' class (specific, moderate, broad). Designed for the
#' "breakdown" figure used in the TDLN analysis: it shows at a glance how much
#' of the significant-gene pool is cell-type-restricted versus broadly
#' expressed across many cell types (the heuristic flag for potential
#' ambient-RNA / segmentation artefacts).
#'
#' @inheritParams plot_gene_counts_by_celltype
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \dontrun{
#' results <- fread("all_genes_results.csv")
#' plot_specificity_breakdown(results, broad_threshold = 5)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_text scale_fill_manual labs
#'   theme_bw theme element_text
#' @export
plot_specificity_breakdown <- function(results,
                                       fdr_col = "fisher_fdr",
                                       fdr_threshold = 0.05,
                                       broad_threshold = 4) {
  if (!inherits(results, "data.table")) {
    results <- data.table::as.data.table(results)
  }

  spec <- classify_gene_specificity(
    results,
    fdr_col         = fdr_col,
    fdr_threshold   = fdr_threshold,
    broad_threshold = broad_threshold
  )
  if (nrow(spec) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = "No significant genes at the chosen threshold") +
        ggplot2::theme_void()
    )
  }

  # Three classes: specific (1 cell type), moderate (2 up to broad_threshold-1),
  # broad (>= broad_threshold). broad_threshold is the single boundary.
  class_levels <- c("specific", "moderate", "broad")
  counts <- spec[, .(n_genes = .N), by = specificity_class]
  counts <- counts[specificity_class %in% class_levels]
  counts[, specificity_class := factor(specificity_class, levels = class_levels)]
  counts <- counts[order(specificity_class)]

  class_colours <- c(specific = "#1B7837", moderate = "#7FBC41",
                     broad = "#B2182B")

  moderate_label <- if (broad_threshold > 2L) {
    sprintf("Moderate (2 to %d cell types)", broad_threshold - 1L)
  } else {
    "Moderate"
  }
  class_labels <- c(
    specific = "Specific (1 cell type)",
    moderate = moderate_label,
    broad    = sprintf("Broad (>= %d cell types)", broad_threshold)
  )

  ggplot2::ggplot(counts,
                  ggplot2::aes(x = specificity_class, y = n_genes,
                               fill = specificity_class)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = n_genes), vjust = -0.4, size = 3.5) +
    ggplot2::scale_fill_manual(values = class_colours,
                               labels = class_labels, name = NULL,
                               drop = FALSE) +
    ggplot2::scale_x_discrete(labels = class_labels, drop = FALSE) +
    ggplot2::labs(
      x = NULL, y = "Number of significant genes",
      title = "Gene specificity breakdown",
      subtitle = sprintf(
        "FDR < %g | 'specific' + 'moderate' = likely real biology, 'broad' >= %d cell types",
        fdr_threshold, broad_threshold)
    ) +
    theme_ripple(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "grey40"),
      axis.text.x = ggplot2::element_text(angle = 20, hjust = 1),
      legend.position = "none"
    )
}


#' RIPPLE QC dashboard
#'
#' Multi-panel QC dashboard composed from the artefacts written by
#' \code{\link{run_ripple}}. Bundles replicate-quality panels (per-sample
#' distance density, cell composition, sign consistency, dispersion) with
#' specificity panels (gene counts split by direction and specificity
#' class, specificity breakdown, top widely-shared genes with optional
#' query-marker bleed-through highlighting -- the latter is the
#' smoking-gun pattern for ambient-RNA contamination).
#'
#' Reads only files written under \code{results_dir}:
#' \itemize{
#'   \item \code{summary/all_genes_results.csv} (required)
#'   \item \code{qc/cell_distances.csv.gz} (optional; enables panels 1+2)
#' }
#'
#' @param results_dir Path to a RIPPLE output directory (the one containing
#'   \code{summary/} and \code{qc/} subdirectories).
#' @param query_signature_genes Optional character vector of query marker
#'   genes. When supplied, genes in the bleed-through panel that are also
#'   query markers are highlighted in red -- this is the smoking-gun pattern
#'   for ambient RNA / segmentation contamination (a query marker showing as
#'   "induced" across many target cell types).
#' @param fdr_threshold Numeric. Significance cutoff (default: 0.05).
#' @param broad_threshold Integer. Min number of cell types where a
#'   gene is significant to be flagged as the broad-expression class
#'   (default: \code{4}). The default is illustrative -- tune for your
#'   panel size and annotation granularity. See
#'   \code{\link{classify_gene_specificity}} for guidance.
#' @param query_label Character. Display label for the query cell type
#'   (default: "Query").
#' @param top_n_bleed Integer. Number of top widely-shared genes to show in
#'   the bleed-through panel (default: 15).
#' @param output_file Character or NULL. If set, also saves the assembled
#'   patchwork to this path via \code{ggsave}.
#' @param width,height Numeric. \code{ggsave} dimensions when
#'   \code{output_file} is set (defaults: 14 x 18 inches).
#'
#' @return A \code{patchwork} object (returned invisibly when
#'   \code{output_file} is set).
#'
#' @examples
#' \dontrun{
#' ripple_plot_qc(
#'   results_dir = "results/spatial_analysis_Tumor/ripple",
#'   query_signature_genes = c("KRT19", "EPCAM", "MKI67"),
#'   query_label = "Tumor",
#'   output_file = "qc_dashboard.pdf"
#' )
#' }
#'
#' @importFrom data.table fread as.data.table
#' @importFrom ggplot2 ggplot aes geom_violin geom_col geom_histogram
#'   geom_boxplot geom_hline geom_vline scale_fill_manual stat_summary
#'   labs theme element_text unit ggsave
#' @importFrom patchwork wrap_plots plot_annotation
#' @export
ripple_plot_qc <- function(results_dir,
                           query_signature_genes = NULL,
                           fdr_threshold = 0.05,
                           broad_threshold = 4L,
                           query_label = "Query",
                           top_n_bleed = 15L,
                           output_file = NULL,
                           width = 14, height = 18) {
  qc_dir <- file.path(results_dir, "qc")
  summary_dir <- file.path(results_dir, "summary")

  results_path <- file.path(summary_dir, "all_genes_results.csv")
  if (!file.exists(results_path)) {
    stop("Missing required file: ", results_path, call. = FALSE)
  }
  results <- data.table::fread(results_path)

  # Issue #14: warn once here if the broad-expression class is unreachable.
  # The specificity/bleed-through panels below re-run classify_gene_specificity
  # and would otherwise repeat the same warning, so those calls are muffled.
  n_ct_qc <- data.table::uniqueN(results$cell_type)
  if (n_ct_qc > 0 && n_ct_qc < broad_threshold) {
    warning(
      "broad_threshold (", broad_threshold, ") exceeds the number of cell ",
      "types in the results (", n_ct_qc, "); the specificity and bleed-through ",
      "panels cannot flag broad-expression / contaminating genes. Run on all ",
      "cell types or lower broad_threshold.", call. = FALSE
    )
  }

  dist_path <- file.path(qc_dir, "cell_distances.csv.gz")
  cell_distances <- if (file.exists(dist_path)) {
    data.table::fread(dist_path)
  } else {
    message(
      "Note: ", dist_path, " not found -- ",
      "distance and composition panels will be skipped."
    )
    NULL
  }

  panels <- list()

  # -- Panel 1: per-sample distance distribution (violin) --
  if (!is.null(cell_distances)) {
    panels$dist <- ggplot2::ggplot(
      cell_distances, ggplot2::aes(x = sample_id, y = dist_to_query)
    ) +
      ggplot2::geom_violin(
        fill = "steelblue", color = "grey30",
        alpha = 0.6, scale = "width"
      ) +
      ggplot2::stat_summary(
        fun = stats::median, geom = "point",
        color = "white", size = 1.6
      ) +
      ggplot2::labs(
        title = "Per-sample distance to query",
        subtitle = "Heavy skew or sample-specific spikes flag tissue-boundary effects",
        x = NULL, y = paste0("Distance to ", query_label, " (um)")
      ) +
      theme_ripple(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 8, color = "grey40"),
        axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
      )

    # -- Panel 2: cells per cell type per sample (stacked bar) --
    cell_comp <- cell_distances[, .N, by = .(sample_id, cell_type)]
    panels$comp <- ggplot2::ggplot(
      cell_comp, ggplot2::aes(x = sample_id, y = N, fill = cell_type)
    ) +
      ggplot2::geom_col(position = "stack", width = 0.7) +
      ggplot2::labs(
        title = "Cell composition per sample",
        subtitle = "Watch for sample dropouts or annotation imbalance",
        x = NULL, y = "Cells", fill = "Cell type"
      ) +
      theme_ripple(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 8, color = "grey40"),
        axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
        legend.position = "right",
        legend.text = ggplot2::element_text(size = 8),
        legend.key.size = ggplot2::unit(0.4, "cm")
      )
  }

  # -- Panel 3: sign-consistency histogram --
  if ("sign_consistency" %in% names(results)) {
    panels$sign <- ggplot2::ggplot(
      results[!is.na(sign_consistency)],
      ggplot2::aes(x = sign_consistency)
    ) +
      ggplot2::geom_histogram(
        bins = 20, fill = "#4393C3", color = "white", alpha = 0.85
      ) +
      ggplot2::geom_vline(
        xintercept = 1.0, linetype = "dashed", color = "grey30"
      ) +
      ggplot2::labs(
        title = "Sign consistency across replicates",
        subtitle = "Pile-up at 1.0 = many genes have replicate-consistent direction",
        x = "Fraction of samples agreeing on sign",
        y = "Genes (across all cell types)"
      ) +
      theme_ripple(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 8, color = "grey40")
      )
  }

  # -- Panel 4: median dispersion per cell type --
  if (all(c("median_dispersion", "cell_type") %in% names(results))) {
    panels$disp <- ggplot2::ggplot(
      results[!is.na(median_dispersion)],
      ggplot2::aes(x = cell_type, y = median_dispersion)
    ) +
      ggplot2::geom_boxplot(
        fill = "#92C5DE", color = "grey30", outlier.size = 0.4
      ) +
      ggplot2::geom_hline(
        yintercept = c(0.3, 2), linetype = "dashed", color = "grey50"
      ) +
      ggplot2::labs(
        title = "Median dispersion per cell type",
        subtitle = "Poisson appropriate when most genes sit between 0.3 and 2",
        x = NULL, y = "Median dispersion"
      ) +
      theme_ripple(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 8, color = "grey40"),
        axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
      )
  }

  # -- Panel 5: significant genes per cell type, split by direction x specificity --
  panels$counts <- plot_gene_counts_by_celltype(
    results,
    fdr_threshold   = fdr_threshold,
    broad_threshold = broad_threshold,
    query_label     = query_label
  )

  # -- Panel 6: specificity breakdown --
  # (unreachable-broad_threshold warning already emitted once above)
  panels$spec <- suppressWarnings(plot_specificity_breakdown(
    results,
    fdr_threshold   = fdr_threshold,
    broad_threshold = broad_threshold
  ))

  # -- Panel 7: bleed-through (top widely-shared genes, optional query-marker flag) --
  spec_dt <- suppressWarnings(classify_gene_specificity(
    results,
    fdr_threshold   = fdr_threshold,
    broad_threshold = broad_threshold
  ))
  if (nrow(spec_dt) > 0) {
    top_bleed <- spec_dt[order(-n_celltypes)][seq_len(min(top_n_bleed, .N))]
    if (is.null(query_signature_genes)) {
      top_bleed[, is_query_marker := FALSE]
    } else {
      top_bleed[, is_query_marker := gene %in% query_signature_genes]
    }
    gene_order <- as.character(top_bleed$gene)
    top_bleed[, gene := factor(gene, levels = rev(gene_order))]

    bleed_subtitle <- if (is.null(query_signature_genes)) {
      "Pass `query_signature_genes` to flag query-marker bleed-through (e.g. KRT19 in non-tumor cells)"
    } else {
      "Red bars = query markers showing up significant in many cell types (likely contamination)"
    }

    panels$bleed <- ggplot2::ggplot(
      top_bleed,
      ggplot2::aes(x = n_celltypes, y = gene, fill = is_query_marker)
    ) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::scale_fill_manual(
        values = c("FALSE" = "#92C5DE", "TRUE" = "#B2182B"),
        labels = c("FALSE" = "Other", "TRUE" = "Query marker"),
        name = NULL,
        drop = FALSE
      ) +
      ggplot2::labs(
        title = sprintf(
          "Top %d most widely-shared significant genes",
          nrow(top_bleed)
        ),
        subtitle = bleed_subtitle,
        x = "# cell types where significant", y = NULL
      ) +
      theme_ripple(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 8, color = "grey40"),
        axis.text.y = ggplot2::element_text(face = "italic", size = 9),
        legend.position = if (is.null(query_signature_genes)) "none" else "bottom"
      )
  }

  combined <- patchwork::wrap_plots(panels, ncol = 2) +
    patchwork::plot_annotation(
      title = sprintf("RIPPLE QC dashboard -- query: %s", query_label),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14)
      )
    )

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, combined, width = width, height = height)
    message("Saved: ", output_file)
    return(invisible(combined))
  }

  combined
}


# =============================================================================
# Stage 4 (confounder) classification plots
# =============================================================================

# Internal: canonical class order, palette, and renaming for Stage 4 outputs.
# The package pipeline emits a dynamic class name "<query_label>_specific"
# and the legacy names niche_driven / no_stage2_result; the new plot
# functions canonicalise to confounder_specific / confounder_driven /
# no_conf_result so a single palette serves any dataset.
#
# Palette (coolors base ffc759-607196-babfd1-e8e9ed, with #DB4C70 for
# enhanced and #643A71 for the rare "reversed" class):
#   confounder_specific = #607196  (steel blue, the canonical "good" class)
#   enhanced            = #DB4C70  (rose, strongest positive signal)
#   confounder_driven   = #FFC759  (warm amber, signal explained by control)
#   underpowered        = #BABFD1  (light lavender-grey, neutral)
#   reversed            = #643A71  (deep purple, sign-flip)
#   no_conf_result      = #E8E9ED  (very light grey, "no info")
.confounder_class_levels <- function() {
  c("confounder_specific", "enhanced", "confounder_driven",
    "underpowered", "reversed", "no_conf_result")
}

.confounder_class_palette <- function() {
  c(
    confounder_specific = "#607196",
    enhanced            = "#DB4C70",
    confounder_driven   = "#FFC759",
    underpowered        = "#BABFD1",
    reversed            = "#643A71",
    no_conf_result      = "#E8E9ED"
  )
}

.canonicalize_confounder_classification <- function(x) {
  x <- as.character(x)
  # Anything ending in "_specific" (the dynamic <query_label>_specific) ->
  # confounder_specific. Already-canonical "confounder_specific" is a no-op.
  x[grepl("_specific$", x)] <- "confounder_specific"
  # Legacy -> canonical
  x[x == "niche_driven"]     <- "confounder_driven"
  x[x == "no_stage2_result"] <- "no_conf_result"
  x
}


#' Stage 4 classification per cell type (stacked bar)
#'
#' Stacked bar chart of how Stage 1-significant genes break down across the
#' Stage 4 (confounder-controlled) classification, one bar per target cell
#' type. Cell types are ordered by total gene count descending; the total
#' is annotated above each bar. Use this for the "did the gradient survive
#' confounder control" check across an entire annotation.
#'
#' Class names are canonicalised so the same palette serves any dataset:
#' the dynamic \code{<query_label>_specific} produced by
#' \code{run_ripple_confounder()} becomes \code{confounder_specific}; legacy
#' \code{niche_driven} and \code{no_stage2_result} become
#' \code{confounder_driven} and \code{no_conf_result}.
#'
#' @param stage4_results A \code{data.table} or \code{data.frame} with one
#'   row per gene-per-cell_type, typically the output of
#'   \code{\link{run_ripple_confounder}} (\code{stage2_multitarget_results.csv}
#'   or \code{stage2_all_results.csv} aggregated across cell types).
#' @param classification_column Character. Column carrying the Stage 4
#'   class. Default \code{"classification"}.
#' @param cell_type_column Character. Column with the target cell type.
#'   Default \code{"cell_type"}.
#' @param query_label,control_label Optional character. Used in the
#'   subtitle when supplied (e.g. "Bivariate model: Tumor query, CAF
#'   control"). \code{NULL} suppresses the subtitle.
#' @param title Character or NULL. Plot title.
#' @param base_size Numeric. \code{theme_ripple} base size (default 13).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' s4 <- data.table::fread(system.file(
#'   "extdata/naive_trc_cached/confounder_cxcl12/stage2_multitarget_results.csv",
#'   package = "ripple"))
#' plot_confounder_bar(
#'   s4,
#'   query_label   = "TRC-Ccl21a",
#'   control_label = "TRC-Cxcl12"
#' )
#' }
#'
#' @importFrom data.table as.data.table copy
#' @importFrom ggplot2 ggplot aes geom_bar geom_text scale_fill_manual
#'   scale_y_continuous expansion labs theme element_text margin unit
#' @export
plot_confounder_bar <- function(stage4_results,
                                classification_column = "classification",
                                cell_type_column      = "cell_type",
                                query_label           = NULL,
                                control_label         = NULL,
                                title                 = NULL,
                                base_size             = 13) {
  dt <- data.table::as.data.table(data.table::copy(stage4_results))
  required <- c(classification_column, cell_type_column)
  miss <- setdiff(required, names(dt))
  if (length(miss) > 0) {
    stop("Missing required columns in stage4_results: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }

  dt[, classification := .canonicalize_confounder_classification(
    get(classification_column)
  )]
  if (cell_type_column != "cell_type") {
    data.table::setnames(dt, cell_type_column, "cell_type")
  }

  class_levels <- .confounder_class_levels()
  class_cols   <- .confounder_class_palette()

  dt[, classification := factor(classification, levels = class_levels)]
  ct_totals <- dt[, .(total = .N), by = cell_type][order(-total)]
  dt[, cell_type := factor(cell_type, levels = ct_totals$cell_type)]

  if (is.null(title)) title <- "Stage 4 classification per target cell type"
  subtitle <- if (!is.null(query_label) && !is.null(control_label)) {
    sprintf("Bivariate model: %s query, %s control", query_label, control_label)
  } else if (!is.null(query_label)) {
    sprintf("Bivariate model: %s query", query_label)
  } else {
    NULL
  }

  ggplot2::ggplot(dt, ggplot2::aes(x = cell_type, fill = classification)) +
    ggplot2::geom_bar() +
    ggplot2::geom_text(
      data = ct_totals,
      ggplot2::aes(x = cell_type, y = total, label = total),
      inherit.aes = FALSE, vjust = -0.4, size = 3.2, colour = "grey25"
    ) +
    ggplot2::scale_fill_manual(values = class_cols, name = NULL,
                               drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::labs(
      x = NULL, y = "Stage-1 significant genes",
      title = title,
      subtitle = subtitle
    ) +
    theme_ripple(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "grey25",
                                            margin = ggplot2::margin(b = 6)),
      axis.text.x   = ggplot2::element_text(angle = 30, hjust = 1, size = 10),
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(0.4, "cm")
    )
}


#' Stage 4 confounder scatter (Stage 1 vs Stage 2 gradient score)
#'
#' Scatter of per-gene Stage 1 gradient score vs Stage 2 (bivariate)
#' gradient score for a single target cell type, coloured by the Stage 4
#' classification. Genes on the dashed \eqn{y = x} diagonal kept their
#' Stage 1 magnitude after controlling for the confounder; genes that
#' fell toward zero on the y-axis are confounder-driven; genes pushed
#' further from zero on the y-axis are enhanced. Same class palette and
#' canonicalisation as \code{\link{plot_confounder_bar}}.
#'
#' @param stage4_results A \code{data.table} or \code{data.frame} with
#'   per-gene Stage 4 results for a single target cell type, typically
#'   \code{stage2_all_results.csv} from
#'   \code{\link{run_ripple_confounder}} filtered to one cell type.
#' @param stage1_coef_column Character. Stage 1 gradient score column.
#'   Default \code{"stage1_coef"}.
#' @param stage2_coef_column Character. Stage 2 (bivariate) gradient
#'   score column. Default \code{"stage2_median_coef"}.
#' @param classification_column Character. Stage 4 class column.
#'   Default \code{"classification"}.
#' @param gene_column Character. Gene name column. Default \code{"gene"}.
#' @param label_genes Optional character vector of genes to label via
#'   \code{ggrepel::geom_text_repel}. \code{NULL} (default) = no labels.
#' @param query_label,control_label Optional character. \code{query_label}
#'   appears parenthetically on the x-axis label and in the subtitle
#'   prefix; \code{control_label} appears in the y-axis label
#'   ("controlled for ...").
#' @param title Character or NULL. Plot title.
#' @param base_size Numeric. \code{theme_ripple} base size (default 13).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' s4 <- data.table::fread(system.file(
#'   "extdata/naive_trc_cached/confounder_cxcl12/stage2_all_results.csv",
#'   package = "ripple"))
#' plot_confounder_scatter(
#'   s4,
#'   label_genes   = c("Ccr7", "Sell", "Lef1", "Tcf7", "Cxcr4", "Gzma"),
#'   query_label   = "TRC-Ccl21a",
#'   control_label = "TRC-Cxcl12"
#' )
#' }
#'
#' @importFrom data.table as.data.table copy
#' @importFrom ggplot2 ggplot aes geom_abline geom_hline geom_vline geom_point
#'   scale_colour_manual coord_cartesian labs theme element_text margin unit
#' @importFrom ggrepel geom_text_repel
#' @export
plot_confounder_scatter <- function(stage4_results,
                                    stage1_coef_column    = "stage1_coef",
                                    stage2_coef_column    = "stage2_median_coef",
                                    classification_column = "classification",
                                    gene_column           = "gene",
                                    label_genes           = NULL,
                                    query_label           = NULL,
                                    control_label         = NULL,
                                    title                 = NULL,
                                    base_size             = 13) {
  dt <- data.table::as.data.table(data.table::copy(stage4_results))
  required <- c(stage1_coef_column, stage2_coef_column,
                classification_column, gene_column)
  miss <- setdiff(required, names(dt))
  if (length(miss) > 0) {
    stop("Missing required columns in stage4_results: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }

  dt[, classification := .canonicalize_confounder_classification(
    get(classification_column)
  )]
  class_levels <- .confounder_class_levels()
  class_cols   <- .confounder_class_palette()
  dt[, classification := factor(classification, levels = class_levels)]

  # Per-class counts for subtitle, in canonical order, only classes present
  class_summary <- dt[, .N, by = classification][order(classification)]
  count_str <- paste(
    apply(class_summary, 1, function(r) sprintf("%s (%s)", r[[1]], r[[2]])),
    collapse = "  |  "
  )
  subtitle <- if (!is.null(query_label)) {
    paste0(query_label, ":  ", count_str)
  } else {
    count_str
  }

  axis_max <- max(
    abs(c(dt[[stage1_coef_column]], dt[[stage2_coef_column]])),
    na.rm = TRUE
  ) * 1.05
  if (!is.finite(axis_max) || axis_max == 0) axis_max <- 1

  x_lab <- if (!is.null(query_label)) {
    bquote("Stage 1 gradient score " * beta ~ "(" * .(query_label) * ")")
  } else {
    expression("Stage 1 gradient score " * beta)
  }
  y_lab <- if (!is.null(control_label)) {
    bquote("Stage 4 gradient score " * beta ~
             "(controlled for " * .(control_label) * ")")
  } else {
    expression("Stage 4 gradient score " * beta)
  }

  p <- ggplot2::ggplot(
    dt, ggplot2::aes(
      x = .data[[stage1_coef_column]],
      y = .data[[stage2_coef_column]],
      colour = .data$classification
    )
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey50") +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey80") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey80") +
    ggplot2::geom_point(alpha = 0.55, size = 1.6) +
    ggplot2::scale_colour_manual(values = class_cols, name = NULL,
                                 drop = FALSE) +
    ggplot2::coord_cartesian(xlim = c(-axis_max, axis_max),
                             ylim = c(-axis_max, axis_max)) +
    ggplot2::labs(x = x_lab, y = y_lab,
                  title = title, subtitle = subtitle) +
    theme_ripple(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "grey25",
                                            margin = ggplot2::margin(b = 6)),
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(0.4, "cm")
    )

  if (!is.null(label_genes)) {
    label_dt <- dt[get(gene_column) %in% label_genes]
    if (nrow(label_dt) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = label_dt,
        ggplot2::aes(label = .data[[gene_column]]),
        size = 3.3, max.overlaps = 20, box.padding = 0.4,
        seed = 1, min.segment.length = 0,
        show.legend = FALSE, inherit.aes = TRUE
      )
    }
  }

  p
}


#' Stage 4 confounder attenuation-ratio plot
#'
#' Sibling of \code{\link{plot_confounder_scatter}} that plots the
#' classifier's underlying statistic on the y-axis directly. Per gene:
#' x = Stage 1 gradient score, y = \code{stage2_median_coef / stage1_coef}
#' (the attenuation ratio used by \code{run_ripple_confounder()} to assign
#' the Stage 4 class). Reference lines at y = 1.1 (enhancement cutoff),
#' y = 1.0 ("no change"), y = 0.5 (attenuation cutoff), and y = 0
#' (sign-flip line) make the class boundaries readable directly off the
#' plot.
#'
#' Easier to interpret than the raw scatter when overplotting on the
#' diagonal makes per-class assignment hard to see, because the y-axis
#' position alone tells the reader which Stage 4 class each gene
#' belongs to.
#'
#' @inheritParams plot_confounder_scatter
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' Genes with \code{stage1_coef == 0} (would produce \code{Inf}) or with
#' a missing \code{stage2_median_coef} (\code{no_conf_result} class) are
#' dropped from the plot with a one-line \code{message()} reporting the
#' count. They still appear in the per-class subtitle counts so the
#' total class breakdown is visible.
#'
#' @examples
#' \dontrun{
#' s4 <- data.table::fread(system.file(
#'   "extdata/naive_trc_cached/confounder_cxcl12/stage2_all_results.csv",
#'   package = "ripple"))
#' plot_confounder_ratio(
#'   s4,
#'   label_genes   = c("Ccr7", "Sell", "Lef1", "Tcf7", "Cxcr4", "Gzma"),
#'   query_label   = "TRC-Ccl21a",
#'   control_label = "TRC-Cxcl12"
#' )
#' }
#'
#' @seealso \code{\link{plot_confounder_scatter}} for the raw stage 1 vs
#'   stage 2 scatter.
#'
#' @importFrom data.table as.data.table copy
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline geom_point
#'   scale_colour_manual labs theme element_text margin unit
#' @importFrom ggrepel geom_text_repel
#' @export
plot_confounder_ratio <- function(stage4_results,
                                  stage1_coef_column    = "stage1_coef",
                                  stage2_coef_column    = "stage2_median_coef",
                                  classification_column = "classification",
                                  gene_column           = "gene",
                                  label_genes           = NULL,
                                  query_label           = NULL,
                                  control_label         = NULL,
                                  title                 = NULL,
                                  base_size             = 13) {
  dt <- data.table::as.data.table(data.table::copy(stage4_results))
  required <- c(stage1_coef_column, stage2_coef_column,
                classification_column, gene_column)
  miss <- setdiff(required, names(dt))
  if (length(miss) > 0) {
    stop("Missing required columns in stage4_results: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }

  dt[, classification := .canonicalize_confounder_classification(
    get(classification_column)
  )]
  class_levels <- .confounder_class_levels()
  class_cols   <- .confounder_class_palette()
  dt[, classification := factor(classification, levels = class_levels)]

  # Per-class counts (computed BEFORE dropping rows so no_conf_result and
  # any divide-by-zero genes still appear in the subtitle).
  class_summary <- dt[, .N, by = classification][order(classification)]
  count_str <- paste(
    apply(class_summary, 1, function(r) sprintf("%s (%s)", r[[1]], r[[2]])),
    collapse = "  |  "
  )
  subtitle <- if (!is.null(query_label)) {
    paste0(query_label, ":  ", count_str)
  } else {
    count_str
  }

  # Compute the ratio and drop unplottable rows
  dt[, ratio := get(stage2_coef_column) / get(stage1_coef_column)]
  n_total <- nrow(dt)
  dt <- dt[is.finite(ratio)]
  n_dropped <- n_total - nrow(dt)
  if (n_dropped > 0) {
    message(
      "  plot_confounder_ratio: dropped ", n_dropped,
      " gene(s) with stage1_coef == 0 or missing stage2_median_coef ",
      "(still counted in the subtitle); ", nrow(dt), " plotted."
    )
  }

  if (nrow(dt) == 0) {
    stop("No genes remain after dropping non-finite ratios.", call. = FALSE)
  }

  x_lab <- if (!is.null(query_label)) {
    bquote("Stage 1 gradient score " * beta ~ "(" * .(query_label) * ")")
  } else {
    expression("Stage 1 gradient score " * beta)
  }
  y_lab <- if (!is.null(control_label)) {
    bquote(beta["confounder-controlled"] / beta[Stage~1] ~
             "(" * .(control_label) * ")")
  } else {
    expression(beta["confounder-controlled"] / beta[Stage~1])
  }

  p <- ggplot2::ggplot(
    dt, ggplot2::aes(
      x = .data[[stage1_coef_column]],
      y = .data$ratio,
      colour = .data$classification
    )
  ) +
    # Reference lines, drawn first so points sit on top
    ggplot2::geom_hline(yintercept = 0,   linewidth = 0.3,
                        colour = "grey80") +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dotted",
                        colour = "grey55") +
    ggplot2::geom_hline(yintercept = 1.0, linetype = "dashed",
                        colour = "grey45") +
    ggplot2::geom_hline(yintercept = 1.1, linetype = "dotted",
                        colour = "grey55") +
    ggplot2::geom_vline(xintercept = 0,   linewidth = 0.3,
                        colour = "grey80") +
    ggplot2::geom_point(alpha = 0.6, size = 1.6) +
    ggplot2::scale_colour_manual(values = class_cols, name = NULL,
                                 drop = FALSE) +
    ggplot2::labs(x = x_lab, y = y_lab,
                  title = title, subtitle = subtitle) +
    theme_ripple(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "grey25",
                                            margin = ggplot2::margin(b = 6)),
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(0.4, "cm")
    )

  if (!is.null(label_genes)) {
    label_dt <- dt[get(gene_column) %in% label_genes]
    if (nrow(label_dt) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = label_dt,
        ggplot2::aes(label = .data[[gene_column]]),
        size = 3.3, max.overlaps = 20, box.padding = 0.4,
        seed = 1, min.segment.length = 0,
        show.legend = FALSE, inherit.aes = TRUE
      )
    }
  }

  p
}


