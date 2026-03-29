# =============================================================================
# Hypothesis Visualizations — Decay Plots for Key Biological Themes
# =============================================================================
#
# Creates publication-quality decay plots showing proportion of cells expressing
# a gene as a function of distance to HyMy. Per-sample curves demonstrate
# consistency across biological replicates (N=4 TDLN).
#
# Themes:
#   1. CD8 T cell exhaustion vs stem-like phenotype
#   2. CD4 Foxp3 (Treg enrichment near HyMy)
#   3. LEC stromal remodeling
#   4. FRC inflammatory remodeling (HyMy-specific vs niche-driven)
#   5. Chemokine recruitment (multi-cell-type)
#   6. Inflammation & anti-inflammation balance (multi-cell-type)
#   7. BEC vascular activation
#
# Usage:
#   Rscript hypothesis_visualizations.R
#   ANNOTATION_LEVEL=L1 Rscript hypothesis_visualizations.R
#
# =============================================================================

script_dir <- if (exists("script_dir")) script_dir else {
  tryCatch(dirname(normalizePath(sys.frame(1)$ofile)),
           error = function(e) getwd())
}
source(file.path(script_dir, "utils.R"))

suppressPackageStartupMessages({
  library(ggrepel)
})

# =============================================================================
# Configuration
# =============================================================================

MAX_DISTANCE_UM <- 200
BIN_WIDTH <- 10
MIN_CELLS_PER_BIN <- 10  # Skip bins with very few cells

# Output directory
ANALYSIS_NAME <- Sys.getenv("ANALYSIS_NAME", unset = "hymy_distance_correlation")
VIZ_DIR <- file.path(OUTPUT_ROOT, ANALYSIS_NAME, "hypothesis_figures")
ensure_dir(VIZ_DIR)

# Cell type mapping (aggregated names → actual column values)
# Must match hymy_distance_correlation.R TARGET_CELLTYPES
CELLTYPE_MAP <- list(
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
  mature_migDC = "mature_migDC",
  B_cells = c("B_cell", "Follicular_B"),
  Plasma_cell = "Plasma_cell"
)

# Color palette
SAMPLE_COLORS <- c(
  "m16" = "#E74C3C", "m17" = "#3498DB", "m18" = "#2ECC71", "m19" = "#9B59B6"
)
INDUCED_COLOR <- "#B2182B"   # red — negative coefficient (higher near HyMy)
REPRESSED_COLOR <- "#2166AC" # blue — positive coefficient (lower near HyMy)
NICHE_COLOR <- "#95A5A6"     # gray — niche-driven (control)

# =============================================================================
# Gene Themes
# =============================================================================

themes <- list(
  cd8_exhaustion = list(
    cell_type = "CD8_T_cells",
    title = "CD8 T Cell Exhaustion Near HyMy",
    bin_width = 10,
    genes = list(
      list(gene = "Havcr2", label = "Tim-3", direction = "induced"),
      list(gene = "Pdcd1", label = "PD-1", direction = "induced"),
      list(gene = "Entpd1", label = "CD39", direction = "induced"),
      list(gene = "Tcf7", label = "TCF1 (stem-like)", direction = "repressed"),
      list(gene = "Sell", label = "CD62L (naive)", direction = "repressed"),
      list(gene = "Il7r", label = "CD127 (memory)", direction = "repressed")
    ),
    ncol = 3
  ),

  cd4_foxp3 = list(
    cell_type = "CD4_T_cells",
    title = "CD4 T Cell Foxp3 Induction Near HyMy",
    genes = list(
      list(gene = "Foxp3", label = "Foxp3 (Treg)", direction = "induced"),
      list(gene = "Ctla4", label = "CTLA-4", direction = "induced"),
      list(gene = "Entpd1", label = "CD39", direction = "induced"),
      list(gene = "Tnfrsf9", label = "4-1BB", direction = "induced"),
      list(gene = "Themis", label = "Themis (naive)", direction = "repressed"),
      list(gene = "Trbc2", label = "TCR-beta (naive)", direction = "repressed")
    ),
    ncol = 3
  ),

  lec_remodeling = list(
    cell_type = "LEC",
    title = "LEC Stromal Remodeling Near HyMy",
    genes = list(
      list(gene = "Spon1", label = "Spondin-1 (ECM)", direction = "induced"),
      list(gene = "Hspg2", label = "Perlecan (ECM)", direction = "induced"),
      list(gene = "Il6ra", label = "IL-6Rα", direction = "repressed"),
      list(gene = "Sesn1", label = "Sestrin-1 (stress)", direction = "repressed"),
      list(gene = "Ifngr1", label = "IFNγR1", direction = "repressed"),
      list(gene = "Cd27", label = "CD27", direction = "repressed")
    ),
    ncol = 3
  ),

  frc_remodeling = list(
    cell_type = "FRC",
    title = "FRC Inflammatory Remodeling Near HyMy",
    genes = list(
      list(gene = "Cxcl5", label = "CXCL5 (HyMy-specific)", direction = "induced"),
      list(gene = "Cebpb", label = "C/EBPβ (HyMy-specific)", direction = "induced"),
      list(gene = "Sod2", label = "SOD2 (HyMy-specific)", direction = "induced"),
      list(gene = "Cxcl10", label = "CXCL10 (HyMy-specific)", direction = "induced"),
      list(gene = "Col4a1", label = "Col4a1 (niche-driven)", direction = "niche"),
      list(gene = "Col1a2", label = "Col1a2 (niche-driven)", direction = "niche")
    ),
    ncol = 3
  ),

  # --- Multi-cell-type themes (per-gene cell_type override) ---

  chemokine_recruitment = list(
    cell_type = NULL,  # multi-cell-type — each gene specifies its own
    title = "Chemokine & Recruitment Signaling Near HyMy",
    genes = list(
      list(gene = "Ccr1",   label = "CCR1 (CD8 T)",   direction = "induced", cell_type = "CD8_T_cells"),
      list(gene = "Ccr5",   label = "CCR5 (CD4 T)",   direction = "induced", cell_type = "CD4_T_cells"),
      list(gene = "Cxcl10", label = "CXCL10 (FRC)",   direction = "induced", cell_type = "FRC"),
      list(gene = "Cxcl5",  label = "CXCL5 (FRC)",    direction = "induced", cell_type = "FRC"),
      list(gene = "Ccl2",   label = "CCL2 (cDC2)",    direction = "induced", cell_type = "cDC2"),
      list(gene = "Cxcl9",  label = "CXCL9 (CD8 T)",  direction = "induced", cell_type = "CD8_T_cells")
    ),
    ncol = 3
  ),

  inflammation_balance = list(
    cell_type = NULL,  # multi-cell-type
    title = "Inflammation & Anti-Inflammation Balance Near HyMy",
    genes = list(
      list(gene = "Il1b",     label = "IL-1\u03b2 (CD4 T)",   direction = "induced",  cell_type = "CD4_T_cells"),
      list(gene = "Tnfrsf9",  label = "4-1BB (CD4 T)",        direction = "induced",  cell_type = "CD4_T_cells"),
      list(gene = "Il1rn",    label = "IL-1RA (Plasma)",       direction = "induced",  cell_type = "Plasma_cell"),
      list(gene = "Il10",     label = "IL-10 (Plasma)",        direction = "induced",  cell_type = "Plasma_cell"),
      list(gene = "Tgfbr1",   label = "TGF\u03b2R1 (FRC)",    direction = "induced",  cell_type = "FRC"),
      list(gene = "Tnfrsf1b", label = "TNFR2 (Mono)",         direction = "repressed", cell_type = "Monocyte")
    ),
    ncol = 3
  ),

  bec_vascular_activation = list(
    cell_type = "BEC",
    title = "BEC Vascular Activation Near HyMy",
    genes = list(
      list(gene = "Selp",   label = "P-selectin",       direction = "induced"),
      list(gene = "Sele",   label = "E-selectin",       direction = "induced"),
      list(gene = "Tek",    label = "Tie2",             direction = "induced"),
      list(gene = "Plvap",  label = "PLVAP",            direction = "induced"),
      list(gene = "Ackr1",  label = "ACKR1/DARC",      direction = "induced"),
      list(gene = "Pecam1", label = "CD31",             direction = "repressed")
    ),
    ncol = 3
  ),

  csf3_il33_csf1_stromal = list(
    cell_type = NULL,  # multi-cell-type — FRC and LEC
    title = "Stromal Ligands (CSF3, IL-33, CSF1) Near HyMy",
    genes = list(
      list(gene = "Csf3", label = "CSF3 (FRC)",  direction = "induced", cell_type = "FRC"),
      list(gene = "Csf3", label = "CSF3 (LEC)",  direction = "induced", cell_type = "LEC"),
      list(gene = "Il33", label = "IL-33 (FRC)",  direction = "induced", cell_type = "FRC"),
      list(gene = "Il33", label = "IL-33 (LEC)",  direction = "induced", cell_type = "LEC"),
      list(gene = "Csf1", label = "CSF1 (FRC)",   direction = "induced", cell_type = "FRC"),
      list(gene = "Csf1", label = "CSF1 (LEC)",   direction = "induced", cell_type = "LEC")
    ),
    ncol = 3
  ),

  lec_emt = list(
    cell_type = "LEC",
    title = "Epithelial\u2013Mesenchymal Transition (EndMT) in LEC Near HyMy",
    genes = list(
      list(gene = "Ccn1",     label = "CCN1/CYR61",      direction = "induced"),
      list(gene = "Cxcl12",   label = "CXCL12/SDF-1",    direction = "induced"),
      list(gene = "Serpine1",  label = "Serpine1/PAI-1",   direction = "induced"),
      list(gene = "Thbs1",    label = "Thrombospondin-1", direction = "induced"),
      list(gene = "Vegfc",    label = "VEGF-C",           direction = "induced"),
      list(gene = "Spp1",     label = "Osteopontin",      direction = "induced")
    ),
    ncol = 3
  ),

  lec_il6_jak_stat3 = list(
    cell_type = "LEC",
    title = "IL-6/JAK/STAT3 Signaling in LEC Near HyMy",
    genes = list(
      list(gene = "Socs3",  label = "SOCS3",     direction = "induced"),
      list(gene = "Pim1",   label = "PIM1",      direction = "induced"),
      list(gene = "Osmr",   label = "OSMR",      direction = "induced"),
      list(gene = "Hmox1",  label = "HO-1",      direction = "induced"),
      list(gene = "Il7",    label = "IL-7",       direction = "induced"),
      list(gene = "Tlr2",   label = "TLR2",      direction = "induced")
    ),
    ncol = 3
  ),

  frc_inflammation = list(
    cell_type = "FRC",
    title = "FRC Inflammatory Signaling Near HyMy",
    genes = list(
      list(gene = "Cxcl1",    label = "CXCL1",          direction = "induced"),
      list(gene = "Ccl7",     label = "CCL7",            direction = "induced"),
      list(gene = "Tnfsf13b", label = "BAFF/BLyS",      direction = "induced"),
      list(gene = "Il1rap",   label = "IL-1RAcP",        direction = "induced"),
      list(gene = "Il1rn",    label = "IL-1RA",          direction = "induced"),
      list(gene = "C1qb",     label = "C1qB (Complement)", direction = "induced")
    ),
    ncol = 3
  ),

  cd8_proliferation = list(
    cell_type = "CD8_T_cells",
    title = "CD8 T Cell Proliferation Near HyMy",
    bin_width = 10,
    genes = list(
      list(gene = "Mki67",  label = "Ki-67",       direction = "induced"),
      list(gene = "Top2a",  label = "TOP2A",       direction = "induced"),
      list(gene = "Cdk1",   label = "CDK1",        direction = "induced"),
      list(gene = "Ccna2",  label = "Cyclin A2",   direction = "induced"),
      list(gene = "Foxm1",  label = "FOXM1",       direction = "induced"),
      list(gene = "Birc5",  label = "Survivin",    direction = "induced")
    ),
    ncol = 3
  )
)

# =============================================================================
# Load Data
# =============================================================================

message("\n", strrep("=", 70))
message("Loading Data")
message(strrep("=", 70))

obj <- load_seurat()
if (USE_HYMY_ANNOTATION) {
  obj <- merge_hymy_annotations(obj)
}

CELLTYPE_COL <- CELLTYPE_COLUMN

# Convert metadata to data.table
cell_data <- as.data.table(obj@meta.data, keep.rownames = "barcode")

# Add condition from 'group' column
cell_data[, condition := group]

message("Total cells: ", nrow(cell_data))
message("TDLN cells: ", sum(cell_data$condition == "TDLN"))

# =============================================================================
# Calculate Distances to HyMy (TDLN only)
# =============================================================================

message("\n", strrep("=", 70))
message("Calculating Distances to HyMy (TDLN only)")
message(strrep("=", 70))

# Filter to TDLN
tdln_data <- cell_data[condition == "TDLN"]

# Get spatial coordinates
coords_all <- as.matrix(tdln_data[, .(spatial_x, spatial_y)])

# Get query cell coordinates
query_mask <- tdln_data[[CELLTYPE_COL]] == QUERY_CELLTYPE
query_coords <- coords_all[query_mask, , drop = FALSE]
message("Query cells (", QUERY_CELLTYPE, ") in TDLN: ", sum(query_mask))

# Calculate nearest neighbor distances
nn_result <- nn2(query_coords, coords_all, k = 1)
tdln_data[, dist_to_query := pmin(as.vector(nn_result$nn.dists), MAX_DISTANCE_UM)]

message("Distance distribution (TDLN):")
message("  Median: ", round(median(tdln_data$dist_to_query), 1), " µm")
message("  Mean: ", round(mean(tdln_data$dist_to_query), 1), " µm")

# =============================================================================
# Load Meta-Analysis Results (for annotation)
# =============================================================================

results_path <- file.path(OUTPUT_ROOT, ANALYSIS_NAME, "summary",
                          "all_genes_results.csv")
if (file.exists(results_path)) {
  all_results <- fread(results_path)
  message("Loaded meta-analysis results: ", nrow(all_results), " gene-celltype pairs")
} else {
  message("WARNING: Meta-analysis results not found at ", results_path)
  all_results <- NULL
}

# =============================================================================
# Decay Plot Function
# =============================================================================

#' Create a single decay plot for one gene in one cell type
#'
#' @param tdln_dt data.table of TDLN cells with dist_to_query
#' @param obj Seurat object (for expression data)
#' @param cell_type Cell type to filter (from CELLTYPE_COL)
#' @param gene Gene name
#' @param label Display label for the gene
#' @param direction "induced", "repressed", or "niche"
#' @param results_dt Optional meta-analysis results for annotation
#' @return ggplot object
make_decay_plot <- function(tdln_dt, obj, ct_name, g_name, label, direction,
                            results_dt = NULL, bin_width = BIN_WIDTH) {

  # Map aggregated cell type name to actual column values
  ct_values <- if (ct_name %in% names(CELLTYPE_MAP)) CELLTYPE_MAP[[ct_name]] else ct_name
  target_cells <- copy(tdln_dt[get(CELLTYPE_COL) %in% ct_values])

  if (nrow(target_cells) == 0) {
    message("  No cells for ", ct_name, " — skipping ", g_name)
    return(NULL)
  }

  # Check gene exists in the Seurat object
  if (!g_name %in% rownames(obj)) {
    message("  Gene ", g_name, " not found in expression matrix — skipping")
    return(NULL)
  }

  # Get expression for target cells
  expr_vec <- as.numeric(GetAssayData(obj, layer = "data")[g_name, target_cells$barcode])
  target_cells[, expressing := as.integer(expr_vec > 0)]

  # Create distance bins
  target_cells[, dist_bin := cut(dist_to_query,
                                  breaks = seq(0, MAX_DISTANCE_UM, by = bin_width),
                                  include.lowest = TRUE)]
  target_cells[, dist_mid := as.numeric(sub("\\(|\\[", "",
                                             sub(",.*", "", as.character(dist_bin)))) + bin_width / 2]

  # Per-sample bin statistics
  sample_bin_stats <- target_cells[!is.na(dist_mid), .(
    prop_expressing = mean(expressing),
    n_cells = .N
  ), by = .(sample_id, dist_mid)]

  # Filter bins with too few cells
  sample_bin_stats <- sample_bin_stats[n_cells >= MIN_CELLS_PER_BIN]

  # Pooled bin statistics (across samples)
  pooled_stats <- sample_bin_stats[, .(
    mean_prop = mean(prop_expressing),
    se_prop = sd(prop_expressing) / sqrt(.N),
    n_samples = .N
  ), by = dist_mid]
  pooled_stats <- pooled_stats[n_samples >= 2]  # Need ≥2 samples for SE

  if (nrow(pooled_stats) < 3) {
    message("  Too few valid bins for ", g_name, " in ", ct_name, " — skipping")
    return(NULL)
  }

  # Color by direction
  line_color <- switch(direction,
                       induced = INDUCED_COLOR,
                       repressed = REPRESSED_COLOR,
                       niche = NICHE_COLOR,
                       INDUCED_COLOR)

  # Short sample labels (extract mouse ID)
  sample_bin_stats[, sample_short := sub(".*__(m[0-9]+)$", "\\1", sample_id)]

  # Build plot
  p <- ggplot() +
    # Per-sample lines (thin, semi-transparent)
    geom_line(data = sample_bin_stats,
              aes(x = dist_mid, y = prop_expressing, group = sample_id,
                  color = sample_short),
              linewidth = 0.4, alpha = 0.5) +
    # Pooled mean with SE ribbon
    geom_ribbon(data = pooled_stats,
                aes(x = dist_mid,
                    ymin = pmax(0, mean_prop - 1.96 * se_prop),
                    ymax = pmin(1, mean_prop + 1.96 * se_prop)),
                fill = line_color, alpha = 0.15) +
    geom_line(data = pooled_stats,
              aes(x = dist_mid, y = mean_prop),
              color = line_color, linewidth = 1.2) +
    scale_color_manual(values = SAMPLE_COLORS, name = "Sample") +
    ylim(0, NA) +
    labs(
      title = label,
      x = "Distance to HyMy (µm)",
      y = "P(expressing)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11, color = line_color),
      legend.position = "none"
    )

  # Annotate with summary stats (prefer median_coef + Fisher FDR if available)
  if (!is.null(results_dt)) {
    gene_result <- results_dt[cell_type == ct_name & gene == g_name]
    if (nrow(gene_result) == 1) {
      coef_col <- if ("median_coef" %in% names(gene_result)) "median_coef" else "combined_coef"
      fdr_col <- if ("fisher_fdr" %in% names(gene_result)) "fisher_fdr" else "fdr"
      coef_val <- gene_result[[coef_col]]
      fdr_val <- gene_result[[fdr_col]]
      coef_str <- sprintf("median β = %.4f", coef_val)
      fdr_str <- if (fdr_val < 0.001) sprintf("FDR = %.1e", fdr_val) else sprintf("FDR = %.3f", fdr_val)
      p <- p + annotate("text",
                         x = MAX_DISTANCE_UM * 0.98, y = Inf,
                         label = paste0(coef_str, "\n", fdr_str),
                         hjust = 1, vjust = 1.3, size = 3,
                         fontface = "italic", color = "gray30")
    }
  }

  return(p)
}

# =============================================================================
# Generate All Theme Panels
# =============================================================================

message("\n", strrep("=", 70))
message("Generating Decay Plots")
message(strrep("=", 70))

for (theme_name in names(themes)) {
  theme <- themes[[theme_name]]
  is_multi_ct <- is.null(theme$cell_type)
  message("\n--- Theme: ", theme$title, " ---")
  if (is_multi_ct) {
    message("  Cell type: multi-cell-type (per-gene)")
  } else {
    message("  Cell type: ", theme$cell_type)
  }
  message("  Genes: ", length(theme$genes))

  # Use theme-level bin_width if specified, otherwise global default
  theme_bin_width <- if (!is.null(theme$bin_width)) theme$bin_width else BIN_WIDTH

  plot_list <- list()
  for (g_info in theme$genes) {
    # Per-gene cell_type override: use gene-level if present, else theme-level
    gene_ct <- if (!is.null(g_info$cell_type)) g_info$cell_type else theme$cell_type
    if (is.null(gene_ct)) {
      message("  WARNING: No cell_type for ", g_info$gene, " — skipping")
      next
    }
    message("  Processing: ", g_info$gene, " (", g_info$label, ") in ", gene_ct)
    p <- make_decay_plot(
      tdln_dt = tdln_data,
      obj = obj,
      ct_name = gene_ct,
      g_name = g_info$gene,
      label = g_info$label,
      direction = g_info$direction,
      results_dt = all_results,
      bin_width = theme_bin_width
    )
    if (!is.null(p)) {
      plot_key <- if (is_multi_ct) paste0(g_info$gene, "_", gene_ct) else g_info$gene
      plot_list[[plot_key]] <- p
    }
  }

  if (length(plot_list) == 0) {
    message("  No valid plots generated — skipping theme")
    next
  }

  # Add shared legend from one representative plot
  # Recreate one plot with legend to extract it
  legend_plot <- plot_list[[1]] +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 1))

  # Subtitle: show cell type for single-ct themes, generic for multi-ct
  if (is_multi_ct) {
    subtitle_text <- "Multiple cell types  |  N = 4 TDLN samples  |  Thin lines = per-sample, thick = pooled mean \u00b1 95% CI"
  } else {
    subtitle_text <- sprintf("%s  |  N = 4 TDLN samples  |  Thin lines = per-sample, thick = pooled mean \u00b1 95%% CI",
                             theme$cell_type)
  }

  # Combine plots
  n_col <- theme$ncol
  n_row <- ceiling(length(plot_list) / n_col)
  combined <- wrap_plots(plot_list, ncol = n_col) +
    plot_annotation(
      title = theme$title,
      subtitle = subtitle_text,
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 9, color = "gray40")
      )
    )

  # Save
  width <- n_col * 4
  height <- n_row * 3.5 + 1  # Extra for title
  out_path <- file.path(VIZ_DIR, sprintf("decay_%s.pdf", theme_name))
  ggsave(out_path, combined, width = width, height = height)
  message("  Saved: ", out_path)
}

# =============================================================================
# Bonus: Combined Foxp3-Only Highlight Panel
# =============================================================================

message("\n--- Foxp3 Highlight Panel ---")

foxp3_plot <- make_decay_plot(
  tdln_dt = tdln_data,
  obj = obj,
  ct_name = "CD4_T_cells",
  g_name = "Foxp3",
  label = "Foxp3 — Treg Enrichment Near HyMy",
  direction = "induced",
  results_dt = all_results
)

if (!is.null(foxp3_plot)) {
  # Enhanced single-panel version with legend
  foxp3_enhanced <- foxp3_plot +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 14)
    ) +
    guides(color = guide_legend(nrow = 1, title = "Sample"))

  ggsave(file.path(VIZ_DIR, "foxp3_decay_highlight.pdf"),
         foxp3_enhanced, width = 6, height = 5)
  message("  Saved: foxp3_decay_highlight.pdf")
}

# =============================================================================
# Summary
# =============================================================================

message("\n", strrep("=", 70))
message("Hypothesis Visualizations Complete")
message("Output directory: ", VIZ_DIR)
message(strrep("=", 70))
