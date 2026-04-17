#' @title Visualization Functions
#'
#' @description Plotting functions for spatial data, volcano plots, decay curves,
#'   and other RIPPLE analysis visualizations.
#'
#' @name visualization
NULL

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

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data[[coef_col]], y = .data$neg_log10_fdr)
  ) +
    ggplot2::geom_point(ggplot2::aes(size = .data$significant), alpha = 0.6) +
    ggplot2::scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2.5), guide = "none") +
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
      box.padding = 0.5
    ) +
    ggplot2::labs(
      title = if (!is.null(title)) title else "Distance-Expression Gradient Volcano",
      subtitle = sprintf(
        "%d genes significant (FDR < %.2f)",
        sum(plot_data$significant, na.rm = TRUE), fdr_threshold
      ),
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
  if (is.null(bin_stats) || nrow(bin_stats) == 0) {
    return(NULL)
  }

  # Build subtitle from meta-analysis info
  sub_text <- ""
  if (!is.null(meta_coef) && !is.null(meta_fdr)) {
    sub_text <- sprintf("coef = %.4f | FDR = %.1e", meta_coef, meta_fdr)
  }

  p <- ggplot2::ggplot(bin_stats, ggplot2::aes(
    x = .data$dist_mid,
    y = .data$prop_expressing
  )) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = pmax(0, .data$prop_expressing - 1.96 * .data$se),
        ymax = pmin(1, .data$prop_expressing + 1.96 * .data$se)
      ),
      fill = color, alpha = 0.15
    ) +
    ggplot2::geom_line(color = color, linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(size = .data$n_cells),
      color = color, alpha = 0.7
    ) +
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
      title = sprintf(
        "Distance-Expression Analysis (Poisson GLM, k=%d): %s",
        k_neighbors, cell_type
      ),
      subtitle = sprintf(
        "%d genes significant (FDR < %.2f)",
        sum(plot_data$significant, na.rm = TRUE),
        fdr_threshold
      ),
      x = paste0(
        "Log-rate coefficient (negative = ", query_label,
        "-induced)"
      ),
      y = paste0("-log10(", fdr_col, ")")
    ) +
    ggplot2::theme_classic(base_size = 12) +
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


#' Create forest plot for a single gene
#'
#' Shows per-sample Poisson GLM coefficients with 95\% confidence intervals
#' and the inverse-variance weighted combined estimate (diamond). Useful for
#' assessing cross-sample reproducibility of a distance-expression gradient.
#'
#' @param coefs Numeric vector. Per-sample coefficients.
#' @param ses Numeric vector. Per-sample standard errors.
#' @param sample_ids Character vector. Sample identifiers (same length as
#'   \code{coefs}).
#' @param gene Character. Gene name for the plot title.
#' @param cell_type Character. Cell type name for the plot title.
#' @param output_path Character. File path to save the plot.
#' @param query_label Character. Display label for the query cell type
#'   (default: "Query").
#'
#' @return The \code{ggplot} object invisibly, or invisible \code{NULL} if
#'   fewer than two valid samples.
#'
#' @examples
#' \dontrun{
#' create_forest_plot(
#'   coefs = c(-0.005, -0.003, -0.004),
#'   ses = c(0.001, 0.002, 0.001),
#'   sample_ids = c("S1", "S2", "S3"),
#'   gene = "Cxcl12",
#'   cell_type = "LEC",
#'   output_path = "Cxcl12_forest.pdf"
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_vline geom_errorbarh geom_point
#'   scale_shape_manual scale_size_manual labs theme_bw theme element_text
#'   ggsave
#' @importFrom data.table data.table
#' @export
create_forest_plot <- function(coefs, ses, sample_ids, gene, cell_type,
                               output_path, query_label = "Query") {
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

  meta_result <- run_meta_analysis(coefs, ses, sample_ids)

  if (!is.na(meta_result$combined_coef)) {
    meta_row <- data.table::data.table(
      sample = "Combined",
      coef = meta_result$combined_coef,
      se = meta_result$combined_se,
      lower = meta_result$combined_coef - 1.96 * meta_result$combined_se,
      upper = meta_result$combined_coef + 1.96 * meta_result$combined_se
    )
    plot_data <- rbind(plot_data, meta_row)
    plot_data[, is_combined := sample == "Combined"]
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

  i2_display <- if (is.na(meta_result$i2)) {
    "N/A"
  } else {
    sprintf("%.1f%%", meta_result$i2 * 100)
  }
  pval_display <- if (is.na(meta_result$combined_pval)) {
    "N/A"
  } else {
    sprintf("%.2e", meta_result$combined_pval)
  }

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
      x = "Log-rate coefficient (per um)",
      y = NULL,
      title = sprintf("%s in %s (Poisson GLM)", gene, cell_type),
      subtitle = sprintf(
        "p = %s, I2 = %s | %s",
        pval_display, i2_display, sign_text
      )
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  ggplot2::ggsave(output_path, p, width = 6, height = 4)
  invisible(p)
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

  gene_order <- meta_results[gene %in% genes_to_plot][
    order(combined_coef)
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
      x = "Log-rate coefficient (per um)",
      y = NULL,
      title = sprintf(
        "Per-Sample Coefficients: %s (Poisson GLM)",
        cell_type
      ),
      subtitle = sprintf(
        "Top %d significant genes | Points = individual samples",
        length(genes_to_plot)
      )
    ) +
    ggplot2::theme_bw(base_size = 10) +
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

  query_mask <- cell_data[[celltype_column]] == query_celltype
  if (sum(query_mask) < 1) {
    stop("No query cells found for '", query_celltype, "'.", call. = FALSE)
  }
  query_coords <- coords[query_mask, , drop = FALSE]

  # Target cell types (everything except query)
  all_types <- unique(cell_data[[celltype_column]])
  target_types <- setdiff(all_types, c(query_celltype, NA_character_))

  max_k <- max(k_range)
  .msg("Computing distances for k = 1 to ", max_k, "...")

  # One nn2 call with max k
  effective_k <- min(max_k, nrow(query_coords))
  if (effective_k < max_k) {
    .msg("  Only ", nrow(query_coords), " query cells; capping k_range at ",
      effective_k)
    k_range <- k_range[k_range <= effective_k]
  }

  nn <- RANN::nn2(query_coords, coords, k = effective_k)

  # Compute per-k, per-cell-type stats
  .msg("Computing per-cell-type distance statistics...")
  summary_list <- lapply(k_range, function(k) {
    if (k == 1) {
      dists <- as.vector(nn$nn.dists[, 1])
    } else {
      dists <- rowMeans(nn$nn.dists[, 1:k, drop = FALSE])
    }
    dists <- pmin(dists, max_distance_um)

    data.table::rbindlist(lapply(target_types, function(ct) {
      ct_mask <- cell_data[[celltype_column]] == ct
      ct_dists <- dists[ct_mask]
      data.table::data.table(
        k = k,
        cell_type = ct,
        mean_dist = mean(ct_dists, na.rm = TRUE),
        sd_dist = stats::sd(ct_dists, na.rm = TRUE),
        n_cells = sum(ct_mask)
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
    ggplot2::theme_bw(base_size = 11) +
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
    ggplot2::theme_bw(base_size = 11) +
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
