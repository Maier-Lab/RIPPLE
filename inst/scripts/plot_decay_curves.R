#!/usr/bin/env Rscript
# =============================================================================
# RIPPLE: Decay Curve Plots for Selected Genes
# =============================================================================
#
# Generates proportion-expressing vs distance-to-query decay plots for:
#   1. Top summary figure genes (top 10 per cell type for CD8, CD4, FRC, LEC, BEC)
#   2. Tex signature genes in CD8 T cells
#
# These plots show how the probability of expression changes as a function
# of distance from the nearest query cell.
#
# FROZEN LEGACY: env-var-driven standalone original that predates the package.
# It hand-rolls distance (RANN::nn2) and binning instead of
# calculate_distance_to_type() / bin_decay_data(), and uses the old
# CONTAMINATION_THRESHOLD "contamination" relabel. Use the package functions
# (e.g. plot_gradient_curve, plot_decay_curve) for new work. Kept for reference.
#
# Usage:
#   Rscript plot_decay_curves.R
#   QUERY_CELLTYPE=MyType CELLTYPE_COLUMN=my_col Rscript plot_decay_curves.R
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(Seurat)
  library(RANN)
})

# Source shared utilities
script_dir <- if (interactive()) {
  dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
}
source(file.path(script_dir, "utils.R"))

# =============================================================================
# Configuration
# =============================================================================

# Inherited from config.R (via utils.R): QUERY_CELLTYPE, CELLTYPE_COL, OUTPUT_SUFFIX, QUERY_LABEL
MAX_DISTANCE_UM <- 200
BIN_WIDTH <- 10  # um per bin
FDR_THRESHOLD <- 0.05
TOP_N_PER_CT <- 10
CONTAMINATION_THRESHOLD <- 4

RESULTS_BASE <- file.path(OUTPUT_ROOT, "hymy_distance_correlation")

OUTPUT_DIR <- file.path(RESULTS_BASE, "plots", "decay_curves")
ensure_dir(OUTPUT_DIR)

# Query signature genes to exclude (from utils.R: QUERY_SIGNATURE / HYMY_SIGNATURE)
HYMY_SIGNATURE_GENES <- QUERY_SIGNATURE

# Tex signature genes
TEX_SIGNATURE <- c("Havcr2", "Lag3", "Entpd1", "Cd38", "Cd244a")

# Auto-discover cell types from results directories
available_celltypes <- basename(list.dirs(
  file.path(RESULTS_BASE, "per_celltype"), recursive = FALSE
))
if (length(available_celltypes) > 0) {
  KEY_CELLTYPES <- available_celltypes
  # Build identity mapping from discovered cell types
  CELLTYPE_MAPPING <- setNames(as.list(available_celltypes), available_celltypes)
  message("Auto-discovered ", length(available_celltypes), " cell types from results")
} else if (QUERY_CELLTYPE %in% c("HyMy_GMM", "IL1B_myeloid") && nchar(INPUT_PATH) == 0) {
  # Legacy HyMy cell types and mapping
  KEY_CELLTYPES <- c("CD8_T_cells", "CD4_T_cells", "FRC", "LEC", "BEC")
  CELLTYPE_MAPPING <- list(
    LEC = "LEC",
    FRC = "FRC",
    BEC = "BEC",
    CD4_T_cells = c("Naive_CD4", "Tfh", "Treg"),
    CD8_T_cells = c("Naive_CD8", "Activated_CD8", "Cytotoxic_CD8", "Tpex"),
    gdT_cells = "gdT_cell",
    Macrophages = "Macrophages",
    Monocyte = "Monocyte",
    Fibroblasts_mac = "Fibroblasts_mac",
    cDC1 = "cDC1",
    cDC2 = "cDC2",
    mature_migDC = "mature_migDC"
  )
  message("Using legacy HyMy cell type mapping")
} else {
  KEY_CELLTYPES <- c()
  CELLTYPE_MAPPING <- list()
  message("WARNING: No cell types discovered; will attempt from results files")
}

message(strrep("=", 70))
message("Decay Curve Plots")
message(strrep("=", 70))
message("Annotation level: ", ANNOTATION_LEVEL)
message("Output: ", OUTPUT_DIR)

# =============================================================================
# Determine genes to plot
# =============================================================================

# Load meta-analysis results to identify top genes
celltype_dir <- file.path(RESULTS_BASE, "per_celltype")
all_results <- rbindlist(lapply(list.dirs(celltype_dir, recursive = FALSE), function(d) {
  f <- file.path(d, "meta_analysis_results.csv")
  if (!file.exists(f)) return(NULL)
  dt <- fread(f)
  dt[, cell_type := basename(d)]
  dt
}), fill = TRUE)

sig_results <- all_results[fdr < FDR_THRESHOLD]

# Identify contamination genes (sig in >= CONTAMINATION_THRESHOLD cell types)
gene_ct_counts <- sig_results[, .(n_celltypes = uniqueN(cell_type)), by = gene]
contamination_genes <- gene_ct_counts[n_celltypes >= CONTAMINATION_THRESHOLD]$gene

# Select top genes per cell type (same logic as summary figure Panel G)
top_genes_per_ct <- sig_results[
  cell_type %in% KEY_CELLTYPES &
  !gene %in% contamination_genes &
  !gene %in% HYMY_SIGNATURE_GENES
][order(-abs(combined_coef)), head(.SD, TOP_N_PER_CT), by = cell_type]

# Build the gene-celltype request list
decay_requests <- rbind(
  top_genes_per_ct[, .(gene, cell_type, source = "top_summary")],
  data.table(gene = TEX_SIGNATURE, cell_type = "CD8_T_cells", source = "tex_signature")
)

# Deduplicate (a gene might be in both top_summary and tex_signature for CD8)
decay_requests <- unique(decay_requests, by = c("gene", "cell_type"))

message("\nGenes to plot:")
for (ct in KEY_CELLTYPES) {
  ct_genes <- decay_requests[cell_type == ct]
  if (nrow(ct_genes) > 0) {
    message("  ", ct, ": ", paste(ct_genes$gene, collapse = ", "),
            " (", nrow(ct_genes), " genes)")
  }
}

# =============================================================================
# Load Seurat object and compute distances
# =============================================================================

message("\nLoading Seurat object...")
obj <- load_seurat()

if (USE_HYMY_ANNOTATION) {
  message("Merging HyMy annotations...")
  obj <- merge_hymy_annotations(obj)
}

cell_data <- as.data.table(obj@meta.data, keep.rownames = "barcode")

# Set condition (generalized)
if (nchar(CONDITION_COL) > 0 && CONDITION_COL %in% names(cell_data)) {
  cell_data[, condition := get(CONDITION_COL)]
} else if ("condition" %in% names(cell_data)) {
  # use existing
} else if ("group" %in% names(cell_data)) {
  cell_data[, condition := group]
} else {
  cell_data[, condition := "all"]
}

condition_label <- if (nchar(CONDITION_VAL) > 0) CONDITION_VAL else "all conditions"

# Filter by condition if specified
if (nchar(CONDITION_VAL) > 0) {
  cell_data <- cell_data[condition == CONDITION_VAL]
  message(condition_label, " cells: ", nrow(cell_data))
} else {
  message("All cells (no condition filter): ", nrow(cell_data))
}

# Get coordinates (dynamic)
coord_cols <- get_coord_columns(cell_data)
coords <- as.matrix(cell_data[, ..coord_cols])

# Identify query cells and compute distances
query_mask <- cell_data[[CELLTYPE_COL]] == QUERY_CELLTYPE
query_coords <- coords[query_mask, , drop = FALSE]
message("Query cells (", QUERY_CELLTYPE, "): ", sum(query_mask))

nn_result <- nn2(query_coords, coords, k = 1)
cell_data[, dist_to_query := pmin(as.vector(nn_result$nn.dists), MAX_DISTANCE_UM)]

# Get expression matrix (we'll subset per gene)
expr_matrix <- GetAssayData(obj, layer = "data")
# Subset to filtered barcodes
filtered_barcodes <- cell_data$barcode
expr_matrix <- expr_matrix[, filtered_barcodes, drop = FALSE]

message("Expression matrix: ", nrow(expr_matrix), " genes x ", ncol(expr_matrix), " cells")

# =============================================================================
# Helper: compute binned decay data
# =============================================================================

compute_decay_data <- function(gene, cell_type_name) {
  # Map aggregated name to actual annotation values
  target_types <- if (cell_type_name %in% names(CELLTYPE_MAPPING)) {
    CELLTYPE_MAPPING[[cell_type_name]]
  } else {
    cell_type_name
  }
  target_idx <- cell_data[[CELLTYPE_COL]] %in% target_types
  if (sum(target_idx) < 50) return(NULL)

  target_barcodes <- cell_data[target_idx, barcode]
  target_distances <- cell_data[target_idx, dist_to_query]

  # Check gene exists
  if (!gene %in% rownames(expr_matrix)) return(NULL)

  # Get expression
  expr_vec <- as.numeric(expr_matrix[gene, target_barcodes])

  df <- data.table(
    distance = target_distances,
    expressing = as.integer(expr_vec > 0)
  )

  # Bin by distance
  df[, dist_bin := cut(distance,
                        breaks = seq(0, MAX_DISTANCE_UM, by = BIN_WIDTH),
                        include.lowest = TRUE)]

  # Compute stats per bin
  bin_stats <- df[, .(
    prop_expressing = mean(expressing),
    n_cells = .N,
    se = sqrt(mean(expressing) * (1 - mean(expressing)) / .N)
  ), by = dist_bin]

  # Get bin midpoints
  bin_stats[, dist_mid := as.numeric(sub("\\(|\\[", "", sub(",.*", "", as.character(dist_bin)))) + BIN_WIDTH / 2]
  bin_stats <- bin_stats[!is.na(dist_mid)]
  bin_stats[, gene := gene]
  bin_stats[, cell_type := cell_type_name]

  bin_stats
}

# =============================================================================
# Helper: make one decay plot
# =============================================================================

make_decay_plot <- function(bin_stats, gene_name, ct_name, meta_info = NULL) {
  if (is.null(bin_stats) || nrow(bin_stats) == 0) return(NULL)

  # Build subtitle from meta-analysis info
  sub_text <- ""
  if (!is.null(meta_info)) {
    sub_text <- sprintf("coef = %.4f | FDR = %.1e | pattern: %s",
                        meta_info$combined_coef,
                        meta_info$fdr,
                        ifelse("decay_pattern" %in% names(meta_info) && !is.na(meta_info$decay_pattern),
                               meta_info$decay_pattern, "N/A"))
  }

  ggplot(bin_stats, aes(x = dist_mid, y = prop_expressing)) +
    geom_ribbon(aes(ymin = pmax(0, prop_expressing - 1.96 * se),
                    ymax = pmin(1, prop_expressing + 1.96 * se)),
                fill = "#E74C3C", alpha = 0.15) +
    geom_line(color = "#E74C3C", linewidth = 0.8) +
    geom_point(aes(size = n_cells), color = "#E74C3C", alpha = 0.7) +
    scale_size_continuous(range = c(0.8, 3.5), guide = "none") +
    scale_x_continuous(breaks = seq(0, MAX_DISTANCE_UM, by = 50)) +
    ylim(0, min(1, max(bin_stats$prop_expressing, na.rm = TRUE) * 1.3)) +
    labs(
      title = gene_name,
      subtitle = sub_text,
      x = paste0("Distance to ", QUERY_LABEL, " (um)"),
      y = "P(expressing)"
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold.italic", size = 11),
      plot.subtitle = element_text(size = 7, color = "grey40")
    )
}

# =============================================================================
# Generate decay plots per cell type
# =============================================================================

message("\nGenerating decay plots...")

for (ct in KEY_CELLTYPES) {
  ct_requests <- decay_requests[cell_type == ct]
  if (nrow(ct_requests) == 0) next

  genes_to_plot <- ct_requests$gene
  message("\n  ", ct, ": ", length(genes_to_plot), " genes")

  # Compute decay data for all genes
  plot_list <- lapply(genes_to_plot, function(g) {
    bin_data <- compute_decay_data(g, ct)
    if (is.null(bin_data)) {
      message("    [SKIP] ", g, " — insufficient data")
      return(NULL)
    }

    # Get meta-analysis info for subtitle
    meta_row <- all_results[gene == g & cell_type == ct]
    meta_info <- if (nrow(meta_row) > 0) meta_row[1] else NULL

    make_decay_plot(bin_data, g, ct, meta_info)
  })

  # Remove NULLs
  plot_list <- Filter(Negate(is.null), plot_list)

  if (length(plot_list) == 0) {
    message("    [SKIP] No valid plots for ", ct)
    next
  }

  # Combine into multi-panel figure
  ncols <- min(4, length(plot_list))
  nrows <- ceiling(length(plot_list) / ncols)
  combined <- wrap_plots(plot_list, ncol = ncols)

  out_file <- file.path(OUTPUT_DIR, paste0("decay_", ct, ".pdf"))
  ggsave(out_file, combined, width = 4 * ncols, height = 3.5 * nrows)
  message("    Saved: ", basename(out_file), " (", length(plot_list), " genes)")
}

# =============================================================================
# Tex signature: dedicated panel with per-sample overlay
# =============================================================================
# A more detailed view: show each sample as a separate line for the
# Tex signature genes, to visualize inter-sample consistency.

message("\n--- Tex Signature Decay (per-sample) ---")

tex_genes_available <- intersect(TEX_SIGNATURE, rownames(expr_matrix))
message("Tex genes on panel: ", paste(tex_genes_available, collapse = ", "))

if (length(tex_genes_available) > 0) {
  # Compute per-sample decay for CD8 T cells (aggregated from subtypes)
  cd8_subtypes <- CELLTYPE_MAPPING[["CD8_T_cells"]]
  cd8_idx <- cell_data[[CELLTYPE_COL]] %in% cd8_subtypes
  cd8_data <- cell_data[cd8_idx]
  cd8_barcodes <- cd8_data$barcode

  tex_per_sample <- rbindlist(lapply(tex_genes_available, function(g) {
    expr_vec <- as.numeric(expr_matrix[g, cd8_barcodes])

    df <- data.table(
      distance = cd8_data$dist_to_query,
      expressing = as.integer(expr_vec > 0),
      sample_id = cd8_data[[SAMPLE_COL]]
    )

    df[, dist_bin := cut(distance, breaks = seq(0, MAX_DISTANCE_UM, by = BIN_WIDTH),
                          include.lowest = TRUE)]

    # Per-sample bin stats
    bin_stats <- df[, .(
      prop_expressing = mean(expressing),
      n_cells = .N
    ), by = .(dist_bin, sample_id)]

    bin_stats[, dist_mid := as.numeric(sub("\\(|\\[", "", sub(",.*", "", as.character(dist_bin)))) + BIN_WIDTH / 2]
    bin_stats <- bin_stats[!is.na(dist_mid)]
    bin_stats[, gene := g]
    bin_stats
  }))

  # Also compute pooled (all samples) for the ribbon
  tex_pooled <- rbindlist(lapply(tex_genes_available, function(g) {
    bin_data <- compute_decay_data(g, "CD8_T_cells")
    if (!is.null(bin_data)) bin_data[, gene := g]
    bin_data
  }))

  if (nrow(tex_per_sample) > 0) {
    # Get meta-analysis info for each gene
    tex_meta <- all_results[gene %in% tex_genes_available & cell_type == "CD8_T_cells"]

    tex_per_sample[, gene_label := gene]
    if (nrow(tex_meta) > 0) {
      for (g in tex_genes_available) {
        meta_row <- tex_meta[gene == g]
        if (nrow(meta_row) > 0) {
          tex_per_sample[gene == g, gene_label := sprintf(
            "%s (coef=%.4f, FDR=%.1e)",
            g, meta_row$combined_coef[1], meta_row$fdr[1]
          )]
        }
      }
    }

    # Order genes by name for consistent layout
    gene_labels_ordered <- unique(tex_per_sample[order(gene)]$gene_label)
    tex_per_sample[, gene_label := factor(gene_label, levels = gene_labels_ordered)]

    # Per-sample overlay plot
    p_tex <- ggplot(tex_per_sample, aes(x = dist_mid, y = prop_expressing)) +
      geom_line(aes(color = sample_id, group = sample_id), alpha = 0.6, linewidth = 0.5)

    # Add pooled ribbon if available
    if (!is.null(tex_pooled) && nrow(tex_pooled) > 0) {
      tex_pooled[, gene_label := gene]
      if (nrow(tex_meta) > 0) {
        for (g in tex_genes_available) {
          meta_row <- tex_meta[gene == g]
          if (nrow(meta_row) > 0) {
            tex_pooled[gene == g, gene_label := sprintf(
              "%s (coef=%.4f, FDR=%.1e)",
              g, meta_row$combined_coef[1], meta_row$fdr[1]
            )]
          }
        }
      }
      tex_pooled[, gene_label := factor(gene_label, levels = gene_labels_ordered)]

      p_tex <- p_tex +
        geom_ribbon(data = tex_pooled, inherit.aes = FALSE,
                    aes(x = dist_mid, y = prop_expressing,
                        ymin = pmax(0, prop_expressing - 1.96 * se),
                        ymax = pmin(1, prop_expressing + 1.96 * se)),
                    fill = "grey30", alpha = 0.15) +
        geom_line(data = tex_pooled, inherit.aes = FALSE,
                  aes(x = dist_mid, y = prop_expressing),
                  color = "black", linewidth = 1, alpha = 0.8)
    }

    p_tex <- p_tex +
      facet_wrap(~ gene_label, scales = "free_y", ncol = 3) +
      scale_x_continuous(breaks = seq(0, MAX_DISTANCE_UM, by = 50)) +
      scale_color_brewer(palette = "Set2", name = "Sample") +
      labs(
        title = "Tex Signature Gene Decay in CD8 T Cells",
        subtitle = paste0("Per-sample lines + pooled mean (black) | ",
                          "Distance to nearest ", QUERY_LABEL),
        x = paste0("Distance to ", QUERY_LABEL, " (um)"),
        y = "P(expressing)"
      ) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        strip.text = element_text(face = "bold.italic", size = 9),
        legend.position = "bottom"
      )

    out_tex <- file.path(OUTPUT_DIR, "decay_tex_signature_CD8.pdf")
    ggsave(out_tex, p_tex,
           width = 12, height = 4 * ceiling(length(tex_genes_available) / 3))
    message("  Saved: decay_tex_signature_CD8.pdf")
  }
}

# =============================================================================
# Done
# =============================================================================

message("\n", strrep("=", 70))
message("Decay curve plots complete!")
message("Output: ", OUTPUT_DIR)
message(strrep("=", 70))
