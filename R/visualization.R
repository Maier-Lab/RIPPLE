#' @title Visualization Functions
#'
#' @description Plotting functions for spatial data, volcano plots, decay curves,
#'   and other RIPPLE analysis visualizations.
#'
#' @name visualization
NULL

#' Basic spatial scatter plot
#'
#' Creates a spatial scatter plot of cells colored by a variable.
#' Supports both categorical (character/factor) and continuous coloring.
#'
#' @param coords Numeric matrix of spatial coordinates (n x 2).
#' @param color_by Vector to color points by (character, factor, or numeric).
#' @param point_size Numeric. Point size (default: 0.5).
#' @param alpha Numeric. Transparency (default: 0.5).
#' @param palette Named character vector of colors for categorical variables,
#'   or NULL for default palettes.
#' @param title Character or NULL. Plot title.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' coords <- matrix(runif(200), ncol = 2)
#' types <- sample(c("A", "B", "C"), 100, replace = TRUE)
#' plot_spatial_scatter(coords, types, title = "Cell Types")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point coord_fixed labs theme_minimal
#'   theme scale_color_manual scale_color_viridis_c
#' @export
plot_spatial_scatter <- function(coords, color_by, point_size = 0.5, alpha = 0.5,
                                 palette = NULL, title = NULL) {
  df <- data.frame(x = coords[, 1], y = coords[, 2], color = color_by)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, color = .data$color)) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = "X (um)", y = "Y (um)", title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right")

  if (is.factor(color_by) || is.character(color_by)) {
    if (!is.null(palette)) {
      p <- p + ggplot2::scale_color_manual(values = palette)
    }
  } else {
    p <- p + ggplot2::scale_color_viridis_c()
  }

  return(p)
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
      get(color_var) == highlight_value, highlight_value, "Other")]
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

  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$x, y = .data$y,
                                           color = .data[[color_var]])) +
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

    p <- ggplot2::ggplot(samp_data, ggplot2::aes(x = .data$x, y = .data$y,
                                                   color = .data[[color_var]])) +
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
#' @param n_label Integer. Maximum number of genes to label (default: 20).
#' @param title Character or NULL. Plot title. If NULL, auto-generated.
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
                                  n_label = 20,
                                  title = NULL) {
  plot_data <- data.table::copy(results)

  plot_data[, neg_log10_fdr := -log10(get(fdr_col))]
  plot_data[neg_log10_fdr > 50, neg_log10_fdr := 50]

  plot_data[, significant := get(fdr_col) < fdr_threshold]

  top_genes <- head(plot_data[significant == TRUE][order(get(fdr_col))], n_label)

  coef_values <- plot_data[[coef_col]]
  max_score <- max(abs(coef_values), na.rm = TRUE) * 1.1

  p <- ggplot2::ggplot(plot_data,
                        ggplot2::aes(x = .data[[coef_col]], y = .data$neg_log10_fdr)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$significant), alpha = 0.6) +
    ggplot2::scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2.5), guide = "none") +
    ggplot2::geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed",
                         color = "grey40") +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "grey60") +
    ggplot2::xlim(-max_score, max_score) +
    ggrepel::geom_text_repel(
      data = top_genes,
      ggplot2::aes(label = .data$gene),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5
    ) +
    ggplot2::labs(
      title = if (!is.null(title)) title else "Distance-Expression Gradient Volcano",
      subtitle = sprintf("%d genes significant (FDR < %.2f)",
                         sum(plot_data$significant, na.rm = TRUE), fdr_threshold),
      x = paste0("Coefficient (negative = ", query_label, "-induced)"),
      y = "-log10(FDR)"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  return(p)
}


#' Decay curve plot
#'
#' Creates a decay curve showing the proportion of expressing cells as a
#' function of distance from query cells. Combines per-sample thin lines
#' with a pooled mean and 95\% CI ribbon.
#'
#' @param bin_stats A \code{data.table} with binned decay statistics. Must
#'   contain columns \code{dist_mid}, \code{prop_expressing}, \code{n_cells},
#'   and \code{se}.
#' @param gene_name Character. Gene name for the plot title.
#' @param cell_type Character. Cell type name for context.
#' @param meta_coef Numeric or NULL. Meta-analysis coefficient for annotation.
#' @param meta_fdr Numeric or NULL. Meta-analysis FDR for annotation.
#' @param query_label Character. Display label for the query cell type
#'   (default: "Query").
#' @param max_distance Numeric. Maximum distance for x-axis (default: 200).
#' @param color Character. Line/ribbon color (default: "#E74C3C").
#'
#' @return A \code{ggplot} object, or NULL if bin_stats is NULL or empty.
#'
#' @examples
#' \dontrun{
#' plot_decay_curve(
#'   bin_stats = binned_data,
#'   gene_name = "Cxcl12",
#'   cell_type = "LEC",
#'   meta_coef = -0.005,
#'   meta_fdr = 0.001
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point
#'   scale_size_continuous scale_x_continuous ylim labs theme_bw theme
#'   element_text
#' @export
plot_decay_curve <- function(bin_stats, gene_name, cell_type,
                             meta_coef = NULL, meta_fdr = NULL,
                             query_label = "Query",
                             max_distance = 200,
                             color = "#E74C3C") {
  if (is.null(bin_stats) || nrow(bin_stats) == 0) return(NULL)

  # Build subtitle from meta-analysis info
  sub_text <- ""
  if (!is.null(meta_coef) && !is.null(meta_fdr)) {
    sub_text <- sprintf("coef = %.4f | FDR = %.1e", meta_coef, meta_fdr)
  }

  p <- ggplot2::ggplot(bin_stats, ggplot2::aes(x = .data$dist_mid,
                                                y = .data$prop_expressing)) +
    ggplot2::geom_ribbon(ggplot2::aes(
      ymin = pmax(0, .data$prop_expressing - 1.96 * .data$se),
      ymax = pmin(1, .data$prop_expressing + 1.96 * .data$se)),
      fill = color, alpha = 0.15) +
    ggplot2::geom_line(color = color, linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(size = .data$n_cells),
                         color = color, alpha = 0.7) +
    ggplot2::scale_size_continuous(range = c(0.8, 3.5), guide = "none") +
    ggplot2::scale_x_continuous(breaks = seq(0, max_distance, by = 50)) +
    ggplot2::ylim(0, min(1, max(bin_stats$prop_expressing, na.rm = TRUE) * 1.3)) +
    ggplot2::labs(
      title = gene_name,
      subtitle = sub_text,
      x = paste0("Distance to ", query_label, " (um)"),
      y = "P(expressing)"
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold.italic", size = 11),
      plot.subtitle = ggplot2::element_text(size = 7, color = "grey40")
    )

  return(p)
}


#' Violin plot helper
#'
#' Creates a violin plot with overlaid box plot for group comparisons.
#'
#' @param data A \code{data.frame} or \code{data.table}.
#' @param x Character. Column name for x-axis (groups).
#' @param y Character. Column name for y-axis (values).
#' @param fill Character or NULL. Column name for fill color. If NULL, uses
#'   the \code{x} column.
#' @param palette Named character vector of colors, or NULL.
#' @param title Character or NULL. Plot title.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' plot_violin(data, x = "condition", y = "entropy", title = "Entropy by Condition")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot labs theme_bw theme
#'   element_text scale_fill_manual
#' @export
plot_violin <- function(data, x, y, fill = NULL, palette = NULL, title = NULL) {
  if (is.null(fill)) fill <- x

  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x]], y = .data[[y]],
                                           fill = .data[[fill]])) +
    ggplot2::geom_violin(scale = "width", trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white") +
    ggplot2::labs(title = title) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (!is.null(palette)) {
    p <- p + ggplot2::scale_fill_manual(values = palette)
  }

  return(p)
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
