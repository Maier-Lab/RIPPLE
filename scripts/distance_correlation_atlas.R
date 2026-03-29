#!/usr/bin/env Rscript
# =============================================================================
# RIPPLE Stage 5: Distance Correlation Atlas — Summary Figures
# =============================================================================
#
# Creates publication-quality summary figures from the merged distance
# correlation results across all cell types.
#
# IMPORTANT: Genes significant in many cell types are MORE likely to be
# segmentation artifacts (query transcripts leaking into neighbors) than
# genuine paracrine effects. This script therefore focuses on CELL-TYPE-
# SPECIFIC genes (significant in <=2 cell types) as the biologically
# interesting hits.
#
# Prerequisites:
#   Rscript merge_permutation_pvals.R
#   Rscript merge_distance_correlation_results.R
#   (Optional) Stage 2 results for classification panels
#
# Usage:
#   Rscript distance_correlation_atlas.R
#   QUERY_CELLTYPE=MyType CELLTYPE_COLUMN=my_col Rscript distance_correlation_atlas.R
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(pheatmap)
})

# Optional: fgsea pathway enrichment (skip panels 8-9 if unavailable)
HAS_FGSEA <- requireNamespace("fgsea", quietly = TRUE) &&
             requireNamespace("msigdbr", quietly = TRUE)
if (HAS_FGSEA) {
  suppressPackageStartupMessages({
    library(fgsea)
    library(msigdbr)
  })
  message("fgsea + msigdbr available: pathway enrichment enabled")
} else {
  message("NOTE: fgsea/msigdbr not installed — pathway panels will be skipped")
  message("  Install with: BiocManager::install(c('fgsea', 'msigdbr'))")
}

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

# Inherited from config.R (via utils.R): QUERY_CELLTYPE, CELLTYPE_COL, OUTPUT_SUFFIX, QUERY_LABEL, ANALYSIS_NAME
CONTROL_CELLTYPE <- Sys.getenv("CONTROL_CELLTYPE", unset = "Monocyte")
FDR_THRESHOLD <- 0.05
CONTAMINATION_THRESHOLD <- 4  # Genes significant in >=N cell types flagged as potential contamination
TOP_PER_CELLTYPE <- 5         # Top genes per cell type in dot plot
TOP_N_DOTPLOT <- 50           # Max genes in dot plot

RESULTS_BASE <- file.path(OUTPUT_ROOT, ANALYSIS_NAME)

OUTPUT_DIR <- file.path(RESULTS_BASE, "plots")
ensure_dir(OUTPUT_DIR)

# Color palettes
induced_label <- paste0(QUERY_LABEL, "-induced")
repressed_label <- paste0(QUERY_LABEL, "-repressed")
DIRECTION_COLORS <- setNames(c("#B2182B", "#2166AC"), c(induced_label, repressed_label))
DIVERGING_PALETTE <- colorRampPalette(c("#B2182B", "white", "#2166AC"))
SPECIFICITY_COLORS <- c("specific" = "#2166AC", "moderate" = "#92C5DE",
                         "ubiquitous" = "#F4A582", "contamination" = "#B2182B")
query_specific_label <- paste0(QUERY_LABEL, "_specific")
CLASSIFICATION_COLORS <- setNames(
  c("#1B9E77", "#66A61E", "#E7298A", "#E6AB02", "#7570B3", "grey70"),
  c(query_specific_label, "enhanced", "niche_driven", "underpowered", "reversed", "not_tested")
)

# Query signature genes for leakage check (E3) — inherited from utils.R
HYMY_SIGNATURE_GENES <- QUERY_SIGNATURE

# Known ligand-receptor pairs for E1 highlighting (optional defaults).
# These are common mouse L-R pairs; they will be silently skipped if not
# present in the user's data (e.g., different organism or gene panel).
LR_PAIRS <- list(
  c("Csf3", "Csf3r"),
  c("Il33", "Il1rl1"),
  c("Cxcl12", "Cxcr4"),
  c("Ccl21a", "Ccr7"),
  c("Il6", "Il6ra"),
  c("Vegfa", "Flt4"),
  c("Il1b", "Il1r1"),
  c("Tnf", "Tnfrsf1a")
)

message(strrep("=", 70))
message("Distance Correlation Atlas")
message(strrep("=", 70))
message("Annotation level: ", ANNOTATION_LEVEL)
message("Results base: ", RESULTS_BASE)
message("Output: ", OUTPUT_DIR)

# =============================================================================
# Load Data
# =============================================================================

# Load all per-cell-type results (full gene set, not just significant)
celltype_dir <- file.path(RESULTS_BASE, "per_celltype")
ct_dirs <- list.dirs(celltype_dir, recursive = FALSE)

all_results <- rbindlist(lapply(ct_dirs, function(d) {
  f <- file.path(d, "meta_analysis_results.csv")
  if (!file.exists(f)) return(NULL)
  dt <- fread(f)
  dt[, cell_type := basename(d)]
  dt
}), fill = TRUE)

message("\nLoaded ", nrow(all_results), " gene-celltype combinations")
message("Cell types: ", paste(unique(all_results$cell_type), collapse = ", "))
n_ct <- uniqueN(all_results$cell_type)

# Load or compute median_coef per gene × cell type.
# Rationale: with small N (e.g. 4 replicates), the meta-analysis TE.random
# weights by precision (cell count), which can let one cell-rich sample
# dominate the ranking. The median treats each sample as an equally weighted
# independent test set.
if ("median_coef" %in% names(all_results)) {
  message("median_coef already in meta_analysis_results.csv — using it directly")
} else {
  # Fallback: compute from coef_per_sample.csv
  per_sample_all <- rbindlist(lapply(ct_dirs, function(d) {
    f <- file.path(d, "coef_per_sample.csv")
    if (!file.exists(f)) return(NULL)
    dt <- fread(f)
    dt[, cell_type := basename(d)]
    dt
  }), fill = TRUE)

  if (nrow(per_sample_all) > 0) {
    median_coefs <- per_sample_all[!is.na(coef), .(
      median_coef = median(coef)
    ), by = .(gene, cell_type)]
    all_results <- merge(all_results, median_coefs, by = c("gene", "cell_type"), all.x = TRUE)
    message("Computed median_coef from per-sample data (",
            nrow(per_sample_all), " sample-gene entries)")
  } else {
    message("WARNING: No per-sample coefficient data found — median_coef unavailable")
    all_results[, median_coef := combined_coef]  # fallback
  }
}

# Select primary significance and coefficient columns
# Prefer Fisher FDR (equal mouse weighting) over meta-analysis FDR if available
SIG_COL <- if ("fisher_fdr" %in% names(all_results)) "fisher_fdr" else "fdr"
COEF_COL <- if ("median_coef" %in% names(all_results)) "median_coef" else "combined_coef"
message("Primary significance column: ", SIG_COL)
message("Primary coefficient column: ", COEF_COL)

# Classify direction for significant genes (using primary sig/coef columns)
all_results[, direction := fifelse(
  get(SIG_COL) < FDR_THRESHOLD & get(COEF_COL) < 0, induced_label,
  fifelse(get(SIG_COL) < FDR_THRESHOLD & get(COEF_COL) > 0, repressed_label, "ns")
)]

sig_results <- all_results[get(SIG_COL) < FDR_THRESHOLD]
message("Significant genes (", SIG_COL, " < ", FDR_THRESHOLD, "): ", nrow(sig_results))

# =============================================================================
# Specificity Classification
# =============================================================================
# Genes significant in many cell types are likely segmentation artifacts
# (query transcripts leaking into neighbors). Flag them accordingly.

gene_ct_counts <- sig_results[, .(
  n_celltypes = uniqueN(cell_type),
  celltypes = paste(sort(unique(as.character(cell_type))), collapse = ", ")
), by = gene]

gene_ct_counts[, specificity := fifelse(
  n_celltypes == 1, "specific",
  fifelse(n_celltypes <= 3, "moderate",
  fifelse(n_celltypes < CONTAMINATION_THRESHOLD, "ubiquitous", "contamination"))
)]

# Add specificity info to all_results
all_results <- merge(all_results,
                     gene_ct_counts[, .(gene, n_celltypes, specificity)],
                     by = "gene", all.x = TRUE)
all_results[is.na(n_celltypes), `:=`(n_celltypes = 0, specificity = "ns")]
sig_results <- merge(sig_results,
                     gene_ct_counts[, .(gene, n_celltypes, specificity)],
                     by = "gene", all.x = TRUE)

message("\nSpecificity breakdown of significant genes:")
message("  Specific (1 cell type):    ", sum(gene_ct_counts$specificity == "specific"))
message("  Moderate (2-3 cell types): ", sum(gene_ct_counts$specificity == "moderate"))
message("  Ubiquitous (4+ cell types):", sum(gene_ct_counts$specificity %in% c("ubiquitous", "contamination")))
message("  Contamination flag (>=", CONTAMINATION_THRESHOLD, "): ",
        sum(gene_ct_counts$specificity == "contamination"))

# Cell type ordering by number of SPECIFIC significant genes (not total)
ct_order <- sig_results[specificity %in% c("specific", "moderate"),
                        .N, by = cell_type][order(-N)]$cell_type
# Add any cell types that only have ubiquitous hits
ct_remaining <- setdiff(unique(all_results$cell_type), ct_order)
ct_order <- c(ct_order, ct_remaining)
all_results[, cell_type := factor(cell_type, levels = ct_order)]
sig_results[, cell_type := factor(cell_type, levels = ct_order)]

# =============================================================================
# Load Stage 2 Results (if available)
# =============================================================================

STAGE2_BASE <- paste0(RESULTS_BASE, "_stage2")
stage2_file <- file.path(STAGE2_BASE, "summary", "stage2_all_results.csv")
HAS_STAGE2 <- file.exists(stage2_file)

if (HAS_STAGE2) {
  message("\nLoading Stage 2 results from: ", stage2_file)
  stage2_results <- fread(stage2_file)
  message("  Stage 2 gene-celltype combinations: ", nrow(stage2_results))

  # Normalize column names: v2 uses stage2_median_coef/stage2_fisher_fdr, v1 uses stage2_coef/stage2_fdr
  if ("stage2_median_coef" %in% names(stage2_results) && !"stage2_coef" %in% names(stage2_results)) {
    setnames(stage2_results, "stage2_median_coef", "stage2_coef")
  }
  if ("stage2_fisher_fdr" %in% names(stage2_results) && !"stage2_fdr" %in% names(stage2_results)) {
    setnames(stage2_results, "stage2_fisher_fdr", "stage2_fdr")
  }

  # Merge Stage 2 classifications into all_results
  stage2_cols <- intersect(c("stage2_coef", "stage2_se", "stage2_fdr"), names(stage2_results))
  stage2_key <- stage2_results[, c("gene", "cell_type", stage2_cols, "classification"), with = FALSE]
  all_results <- merge(all_results, stage2_key, by = c("gene", "cell_type"), all.x = TRUE)
  all_results[is.na(classification), classification := "not_tested"]
  sig_results <- merge(sig_results, stage2_key, by = c("gene", "cell_type"), all.x = TRUE)
  sig_results[is.na(classification), classification := "not_tested"]

  # Re-apply factor ordering (merge() drops factor levels from merge keys)
  all_results[, cell_type := factor(cell_type, levels = ct_order)]
  sig_results[, cell_type := factor(cell_type, levels = ct_order)]

  # Stage 2 summary
  class_summary <- stage2_results[, .N, by = classification]
  message("\nStage 2 classification summary:")
  for (i in seq_len(nrow(class_summary))) {
    message("  ", class_summary$classification[i], ": ", class_summary$N[i])
  }
} else {
  message("\n[NOTE] Stage 2 results not found — Stage 2 panels will be skipped")
  message("  Expected: ", stage2_file)
  all_results[, classification := "not_tested"]
  sig_results[, classification := "not_tested"]
}

# =============================================================================
# Panel 1: Significant Gene Counts Bar Chart (with specificity breakdown)
# =============================================================================

message("\nPanel 1: Gene counts bar chart (with specificity)...")

count_data <- sig_results[, .N, by = .(cell_type, direction, specificity)]

# Stacked bar: direction on x facet, specificity as fill
count_by_dir <- sig_results[, .(
  total = .N,
  specific = sum(specificity %in% c("specific", "moderate")),
  contamination = sum(specificity %in% c("ubiquitous", "contamination"))
), by = .(cell_type, direction)]

count_long <- melt(count_by_dir,
                   id.vars = c("cell_type", "direction"),
                   measure.vars = c("specific", "contamination"),
                   variable.name = "type", value.name = "count")

p1 <- ggplot(count_long, aes(x = cell_type, y = count, fill = interaction(type, direction))) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(
    values = setNames(
      c("#D95F02", "#FDAE6B", "#1B9E77", "#A1D99B"),
      c(paste0("specific.", induced_label),
        paste0("contamination.", induced_label),
        paste0("specific.", repressed_label),
        paste0("contamination.", repressed_label))
    ),
    labels = setNames(
      c("Induced (specific)", "Induced (ubiquitous)",
        "Repressed (specific)", "Repressed (ubiquitous)"),
      c(paste0("specific.", induced_label),
        paste0("contamination.", induced_label),
        paste0("specific.", repressed_label),
        paste0("contamination.", repressed_label))
    ),
    name = NULL
  ) +
  labs(
    title = "Significant Spatial Gradients per Cell Type",
    subtitle = paste0("FDR < ", FDR_THRESHOLD,
                      " | Faded = potential contamination (sig in \u2265",
                      CONTAMINATION_THRESHOLD, " cell types)"),
    x = NULL,
    y = "Number of significant genes"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "top",
    legend.text = element_text(size = 8)
  )

ggsave(file.path(OUTPUT_DIR, "gene_counts_by_celltype.pdf"), p1,
       width = 9, height = 6)
message("  Saved: gene_counts_by_celltype.pdf")

# =============================================================================
# Panel 2: Dot Plot — Top Cell-Type-SPECIFIC Genes
# =============================================================================

message("Panel 2: Dot plot (cell-type-specific genes)...")

# For each cell type, pick top genes by abs(coef) that are NOT contamination
specific_hits <- sig_results[specificity %in% c("specific", "moderate")]
dotplot_genes <- specific_hits[order(-abs(get(COEF_COL))),
                               head(.SD, TOP_PER_CELLTYPE),
                               by = cell_type]$gene
dotplot_genes <- unique(dotplot_genes)

# If too few, also add top specific genes regardless of cell type
if (length(dotplot_genes) < 10) {
  extra <- specific_hits[order(-abs(get(COEF_COL)))]$gene
  dotplot_genes <- unique(c(dotplot_genes, extra))[1:min(TOP_N_DOTPLOT, length(unique(c(dotplot_genes, extra))))]
}
if (length(dotplot_genes) > TOP_N_DOTPLOT) {
  dotplot_genes <- dotplot_genes[1:TOP_N_DOTPLOT]
}

if (length(dotplot_genes) >= 2) {
  # Get data for these genes across all cell types
  dot_data <- all_results[gene %in% dotplot_genes]
  dot_data[, neg_log10_fdr := -log10(pmax(get(SIG_COL), 1e-50))]
  dot_data[, neg_log10_fdr_capped := pmin(neg_log10_fdr, 20)]

  # Order genes by the cell type they're most significant in, then by coefficient
  gene_primary_ct <- dot_data[get(SIG_COL) < FDR_THRESHOLD][
    order(get(SIG_COL)), head(.SD, 1), by = gene][, .(gene, primary_ct = cell_type)]
  gene_order_dt <- merge(
    dot_data[, .(mean_coef = mean(get(COEF_COL), na.rm = TRUE)), by = gene],
    gene_primary_ct, by = "gene", all.x = TRUE
  )
  gene_order_dt[, primary_ct := factor(primary_ct, levels = ct_order)]
  gene_order <- gene_order_dt[order(primary_ct, mean_coef)]$gene
  dot_data[, gene := factor(gene, levels = gene_order)]

  # Limit for diverging scale
  coef_limit <- max(abs(dot_data[[COEF_COL]]), na.rm = TRUE)

  p2 <- ggplot(dot_data, aes(x = cell_type, y = gene)) +
    geom_point(aes(
      size = neg_log10_fdr_capped,
      color = .data[[COEF_COL]]
    )) +
    scale_color_gradientn(
      colors = DIVERGING_PALETTE(100),
      limits = c(-coef_limit, coef_limit),
      name = "Log-rate\ncoefficient"
    ) +
    scale_size_continuous(
      range = c(0.5, 5),
      name = expression(-log[10](FDR)),
      breaks = c(1, 5, 10, 20)
    ) +
    labs(
      title = "Cell-Type-Specific Spatial Gradients",
      subtitle = paste0("Top ", TOP_PER_CELLTYPE,
                        " specific genes per cell type (sig in \u22643 cell types) | ",
                        "Red = induced near ", QUERY_LABEL),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      panel.grid.major = element_line(color = "grey92")
    )

  ggsave(file.path(OUTPUT_DIR, "dotplot_specific_genes.pdf"), p2,
         width = 10, height = max(6, length(dotplot_genes) * 0.22))
  message("  Saved: dotplot_specific_genes.pdf")
} else {
  message("  [SKIP] Too few specific genes for dot plot")
  p2 <- NULL
}

# =============================================================================
# Panel 3: Multi-Panel Volcano (highlight specific genes, flag contamination)
# =============================================================================

message("Panel 3: Multi-panel volcano (with permutation validation)...")

contamination_genes <- gene_ct_counts[specificity == "contamination"]$gene

volcano_data <- copy(all_results)
volcano_data[, neg_log10_fdr := -log10(pmax(get(SIG_COL), 1e-50))]
volcano_data[, is_contamination := gene %in% contamination_genes]

# Permutation validation layer: classify significant genes as perm-validated or FDR-only
HAS_PERM <- "perm_pval" %in% names(volcano_data) && sum(!is.na(volcano_data$perm_pval)) > 0
if (HAS_PERM) {
  volcano_data[, perm_status := fifelse(
    get(SIG_COL) >= FDR_THRESHOLD, "ns",
    fifelse(is_contamination, "contamination",
    fifelse(!is.na(perm_pval) & perm_pval < 0.05, "perm_validated",
    fifelse(!is.na(perm_pval) & perm_pval >= 0.05, "fdr_only", "no_perm")))
  )]
  n_validated <- sum(volcano_data$perm_status == "perm_validated", na.rm = TRUE)
  n_fdr_only <- sum(volcano_data$perm_status == "fdr_only", na.rm = TRUE)
  message("  Perm-validated: ", n_validated, " | FDR-only (not perm-confirmed): ", n_fdr_only)
} else {
  volcano_data[, perm_status := "no_perm"]
  message("  No permutation p-values available — all points shown as filled")
}

# For labeling: top 20 significant non-contamination genes per cell type (by effect size)
label_genes <- sig_results[!gene %in% contamination_genes][
  order(-abs(get(COEF_COL))), head(.SD, 20), by = cell_type]
label_set <- unique(label_genes[, .(gene, cell_type)])

volcano_data[, show_label := FALSE]
for (i in seq_len(nrow(label_set))) {
  volcano_data[gene == label_set$gene[i] & cell_type == label_set$cell_type[i],
               show_label := TRUE]
}

# Symmetric x-axis — use full range with 5% margin so no genes are cropped
x_max <- max(abs(volcano_data[[COEF_COL]]), na.rm = TRUE) * 1.05
y_max <- max(volcano_data$neg_log10_fdr, na.rm = TRUE) * 1.05

# Build the volcano with permutation-validation encoding:
#   Filled circle (16) = perm-validated (FDR < 0.05 AND perm_pval < 0.05)
#   Open circle  (1)  = FDR-only (FDR < 0.05 but perm_pval >= 0.05)
#   Cross        (4)  = contamination candidate
#   Dot          (16) = non-significant (grey, small)
p3 <- ggplot(volcano_data, aes(x = .data[[COEF_COL]], y = neg_log10_fdr)) +
  # Layer 1: Non-significant genes (small grey background)
  geom_point(
    data = volcano_data[get(SIG_COL) >= FDR_THRESHOLD & !is_contamination],
    color = "grey80", size = 0.3, alpha = 0.3
  ) +
  # Layer 2: Contamination candidates (cross marker)
  geom_point(
    data = volcano_data[is_contamination == TRUE & get(SIG_COL) < FDR_THRESHOLD],
    color = "#FDAE6B", size = 1.2, alpha = 0.6, shape = 4, stroke = 0.6
  )

if (HAS_PERM) {
  p3 <- p3 +
    # Layer 3: FDR-only genes (open circles — not perm-validated)
    geom_point(
      data = volcano_data[perm_status == "fdr_only"],
      aes(color = direction), size = 1.5, alpha = 0.7, shape = 1, stroke = 0.5
    ) +
    # Layer 4: Perm-validated genes (filled circles — high confidence)
    geom_point(
      data = volcano_data[perm_status == "perm_validated"],
      aes(color = direction), size = 1, alpha = 0.6, shape = 16
    )
} else {
  p3 <- p3 +
    # No perm data: all significant genes as filled
    geom_point(
      data = volcano_data[get(SIG_COL) < FDR_THRESHOLD & !is_contamination],
      aes(color = direction), size = 0.5, alpha = 0.4
    )
}

# B1: Sign consistency overlay — flag genes with inconsistent direction across samples
HAS_SIGN_CONSISTENCY <- "sign_consistency" %in% names(volcano_data)
if (HAS_SIGN_CONSISTENCY) {
  inconsistent_data <- volcano_data[get(SIG_COL) < FDR_THRESHOLD & !is.na(sign_consistency) &
                                     sign_consistency < 0.75]
  n_inconsistent <- nrow(inconsistent_data)
  message("  Sign consistency: ", n_inconsistent, " significant genes with sign_consistency < 0.75")
}

p3 <- p3 +
  # Labels for top specific genes
  geom_text_repel(
    data = volcano_data[show_label == TRUE],
    aes(label = gene),
    size = 2.2, max.overlaps = 20, segment.size = 0.2,
    min.segment.length = 0, fontface = "italic"
  ) +
  geom_hline(yintercept = -log10(FDR_THRESHOLD), linetype = "dashed",
             color = "grey50", linewidth = 0.3)

# Add grey ring overlay for inconsistent genes (B1)
if (HAS_SIGN_CONSISTENCY && n_inconsistent > 0) {
  p3 <- p3 +
    geom_point(
      data = inconsistent_data,
      color = "grey50", size = 2.5, shape = 1, stroke = 0.8, alpha = 0.7
    )
}

p3 <- p3 +
  scale_color_manual(
    values = c(DIRECTION_COLORS, "ns" = "grey70"),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(-x_max, x_max), ylim = c(0, y_max)) +
  facet_wrap(~ cell_type, ncol = 4, scales = "free_y") +
  labs(
    title = "Spatial Gradient Volcanos Across Cell Types",
    subtitle = if (HAS_PERM) {
      paste0("Filled = perm-validated (", n_validated, ") | ",
             "Open = FDR-only (", n_fdr_only, ") | ",
             "\u00d7 = contamination (\u2265", CONTAMINATION_THRESHOLD, " cell types)",
             if (HAS_SIGN_CONSISTENCY && n_inconsistent > 0)
               paste0(" | \u25cb = inconsistent sign (", n_inconsistent, ")") else "")
    } else {
      paste0("Labels = top cell-type-specific genes | ",
             "\u00d7 = potential contamination (sig in \u2265",
             CONTAMINATION_THRESHOLD, " cell types)")
    },
    x = paste0("Coefficient [", COEF_COL, "] (negative = ", QUERY_LABEL, "-induced)"),
    y = expression(-log[10](FDR))
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    strip.text = element_text(face = "bold", size = 9)
  )

ggsave(file.path(OUTPUT_DIR, "multi_volcano.pdf"), p3,
       width = 14, height = 3.5 * ceiling(n_ct / 4))
message("  Saved: multi_volcano.pdf")

# =============================================================================
# Panel 4: Permutation Validation Scatter
# =============================================================================

message("Panel 4: Permutation validation scatter...")

# Note: intentionally using meta-analysis fdr (not Fisher) — permutation validates meta-analysis
perm_data <- all_results[!is.na(perm_pval) & !is.na(fdr)]
perm_data[, neg_log10_fdr := -log10(pmax(fdr, 1e-50))]
perm_data[, neg_log10_perm := -log10(pmax(perm_pval, 1e-10))]

if (nrow(perm_data) > 0) {
  # Correlation
  cor_val <- cor(perm_data$neg_log10_fdr, perm_data$neg_log10_perm,
                 use = "complete.obs", method = "spearman")

  p4 <- ggplot(perm_data, aes(x = neg_log10_fdr, y = neg_log10_perm)) +
    geom_point(aes(color = cell_type), size = 1, alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red", alpha = 0.5) +
    geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "red", alpha = 0.5) +
    annotate("text", x = Inf, y = -Inf,
             label = paste0("Spearman rho = ", round(cor_val, 3)),
             hjust = 1.1, vjust = -0.5, size = 3.5, color = "grey30") +
    scale_color_brewer(palette = "Set3", name = "Cell type") +
    labs(
      title = "Meta-analysis FDR vs Permutation P-value",
      subtitle = "Agreement between parametric and non-parametric significance",
      x = expression(-log[10](FDR)),
      y = expression(-log[10](perm~p-value))
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      legend.position = "right"
    )

  ggsave(file.path(OUTPUT_DIR, "permutation_validation.pdf"), p4,
         width = 8, height = 6)
  message("  Saved: permutation_validation.pdf")
} else {
  message("  [SKIP] No permutation p-values available")
  p4 <- NULL
}

# =============================================================================
# Panel 5: Specificity vs Contamination Overview
# =============================================================================

message("Panel 5: Specificity overview...")

spec_summary <- gene_ct_counts[, .N, by = specificity]
spec_summary[, specificity := factor(specificity,
                                      levels = c("specific", "moderate", "ubiquitous", "contamination"))]
spec_summary <- spec_summary[order(specificity)]

p5 <- ggplot(spec_summary, aes(x = specificity, y = N, fill = specificity)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = N), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = SPECIFICITY_COLORS, guide = "none") +
  scale_x_discrete(labels = c(
    "specific" = "Specific\n(1 cell type)",
    "moderate" = "Moderate\n(2-3 cell types)",
    "ubiquitous" = paste0("Ubiquitous\n(4-", CONTAMINATION_THRESHOLD - 1, " cell types)"),
    "contamination" = paste0("Contamination\n(\u2265", CONTAMINATION_THRESHOLD, " cell types)")
  )) +
  labs(
    title = "Gene Specificity Distribution",
    subtitle = paste0("Genes sig in many cell types likely reflect segmentation artifacts, ",
                      "not paracrine signaling"),
    x = NULL,
    y = "Number of significant genes"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(OUTPUT_DIR, "specificity_distribution.pdf"), p5,
       width = 7, height = 5)
message("  Saved: specificity_distribution.pdf")

# =============================================================================
# Panel 6: Top Specific Genes Heatmap
# =============================================================================

message("Panel 6: Specific genes heatmap...")

# Top specific genes per cell type (not contamination)
specific_sig <- sig_results[specificity %in% c("specific", "moderate")]
top_specific <- specific_sig[order(-abs(get(COEF_COL))),
                              head(.SD, 8), by = cell_type]$gene
top_specific <- unique(top_specific)

if (length(top_specific) >= 3) {
  # Build matrix: genes x cell types
  heatmap_data <- dcast(
    all_results[gene %in% top_specific],
    gene ~ cell_type,
    value.var = COEF_COL,
    fill = 0
  )
  heatmap_mat <- as.matrix(heatmap_data[, -1])
  rownames(heatmap_mat) <- heatmap_data$gene

  # Significance annotation (star if FDR < threshold)
  sig_annot <- dcast(
    all_results[gene %in% top_specific],
    gene ~ cell_type,
    value.var = SIG_COL,
    fill = 1
  )
  sig_mat <- as.matrix(sig_annot[, -1])
  rownames(sig_mat) <- sig_annot$gene
  display_mat <- ifelse(sig_mat < FDR_THRESHOLD, "*", "")

  # Row annotation: specificity category
  row_annot <- data.frame(
    specificity = gene_ct_counts[match(top_specific, gene)]$specificity,
    row.names = top_specific
  )

  # Color limits
  max_abs <- max(abs(heatmap_mat), na.rm = TRUE)
  breaks <- seq(-max_abs, max_abs, length.out = 101)

  annot_colors <- list(specificity = SPECIFICITY_COLORS[c("specific", "moderate")])

  pdf(file.path(OUTPUT_DIR, "specific_genes_heatmap.pdf"),
      width = max(8, ncol(heatmap_mat) * 0.8 + 2),
      height = max(6, nrow(heatmap_mat) * 0.3 + 2))

  pheatmap(
    heatmap_mat,
    color = DIVERGING_PALETTE(100),
    breaks = breaks,
    display_numbers = display_mat,
    fontsize_number = 8,
    cluster_rows = nrow(heatmap_mat) >= 2,
    cluster_cols = ncol(heatmap_mat) >= 2,
    annotation_row = row_annot,
    annotation_colors = annot_colors,
    main = paste0("Cell-Type-Specific Gradient Genes (red = induced near",
                  QUERY_LABEL, ")"),
    fontsize_row = 8,
    fontsize_col = 10,
    angle_col = 45
  )

  dev.off()
  message("  Saved: specific_genes_heatmap.pdf (",
          length(top_specific), " genes)")
} else {
  message("  [SKIP] Too few specific genes for heatmap")
}

# =============================================================================
# Panel 7: Contamination Candidates List
# =============================================================================

message("Panel 7: Contamination candidates...")

contamination_list <- gene_ct_counts[specificity == "contamination"][order(-n_celltypes)]
if (nrow(contamination_list) > 0) {
  fwrite(contamination_list, file.path(OUTPUT_DIR, "contamination_candidates.csv"))
  message("  Saved: contamination_candidates.csv (",
          nrow(contamination_list), " genes in >=", CONTAMINATION_THRESHOLD, " cell types)")

  # Show top contamination genes
  message("  Top contamination candidates:")
  for (i in seq_len(min(10, nrow(contamination_list)))) {
    r <- contamination_list[i]
    message("    ", r$gene, " (", r$n_celltypes, " cell types: ", r$celltypes, ")")
  }
} else {
  message("  No contamination candidates found")
}

# =============================================================================
# Panel D1: I² Heterogeneity Distribution
# =============================================================================
# Shows how consistent the spatial gradients are across samples.
# High I² indicates large between-sample variability in effects.

message("\nPanel D1: I² heterogeneity distribution...")

if ("i2" %in% names(all_results)) {
  i2_data <- all_results[get(SIG_COL) < FDR_THRESHOLD & !is.na(i2)]

  if (nrow(i2_data) > 0) {
    # Classify I² levels
    i2_data[, heterogeneity := fifelse(
      i2 < 25, "Low (<25%)",
      fifelse(i2 < 75, "Moderate (25-75%)", "High (>75%)")
    )]
    i2_data[, heterogeneity := factor(heterogeneity,
                                       levels = c("Low (<25%)", "Moderate (25-75%)", "High (>75%)"))]

    p_d1 <- ggplot(i2_data, aes(x = i2, fill = heterogeneity)) +
      geom_histogram(bins = 30, alpha = 0.8, color = "white", linewidth = 0.2) +
      scale_fill_manual(values = c("Low (<25%)" = "#2166AC",
                                    "Moderate (25-75%)" = "#F4A582",
                                    "High (>75%)" = "#B2182B"),
                        name = "Heterogeneity") +
      geom_vline(xintercept = c(25, 75), linetype = "dashed", color = "grey40") +
      facet_wrap(~ cell_type, scales = "free_y", ncol = 4) +
      labs(
        title = expression(I^2 ~ "Heterogeneity of Significant Spatial Gradients"),
        subtitle = paste0("Higher I\u00B2 = inconsistent effects across samples | ",
                          "Only FDR < ", FDR_THRESHOLD, " genes shown"),
        x = expression(I^2 ~ "(% variation due to between-sample heterogeneity)"),
        y = "Number of genes"
      ) +
      theme_bw(base_size = 9) +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        strip.text = element_text(face = "bold"),
        legend.position = "top"
      )

    ggsave(file.path(OUTPUT_DIR, "heterogeneity_I2_distribution.pdf"), p_d1,
           width = 14, height = 3.5 * ceiling(n_ct / 4))
    message("  Saved: heterogeneity_I2_distribution.pdf")

    # Summary stats
    i2_summary <- i2_data[, .(
      n_genes = .N,
      median_i2 = median(i2, na.rm = TRUE),
      pct_low = sum(i2 < 25) / .N * 100,
      pct_high = sum(i2 > 75) / .N * 100
    ), by = cell_type][order(-median_i2)]
    fwrite(i2_summary, file.path(OUTPUT_DIR, "heterogeneity_I2_summary.csv"))
    message("  Saved: heterogeneity_I2_summary.csv")
  } else {
    message("  [SKIP] No significant genes with I² values")
  }
} else {
  message("  [SKIP] I² column not found in results")
}

# =============================================================================
# Panel D2: Sample Contribution Heatmap
# =============================================================================
# Shows which samples contribute most to each cell type's gradient.
# Identifies if effects are driven by specific samples.

message("Panel D2: Sample contribution heatmap...")

# Load per-sample coefficients from one cell type as example structure
sample_coef_file <- file.path(celltype_dir, ct_order[1], "coef_per_sample.csv")
if (file.exists(sample_coef_file)) {
  # Load all per-sample coefficients
  all_sample_coefs <- rbindlist(lapply(as.character(ct_order), function(ct) {
    f <- file.path(celltype_dir, ct, "coef_per_sample.csv")
    if (!file.exists(f)) return(NULL)
    dt <- fread(f)
    dt[, cell_type := ct]
    dt
  }), fill = TRUE)

  # Detect sample column in coef_per_sample.csv (may be sample_id or SAMPLE_COL name)
  SAMPLE_COL_CSV <- if (SAMPLE_COL %in% names(all_sample_coefs)) SAMPLE_COL else if ("sample_id" %in% names(all_sample_coefs)) "sample_id" else NULL

  if (nrow(all_sample_coefs) > 0 && !is.null(SAMPLE_COL_CSV)) {
    # For top specific genes, show coefficient pattern across samples
    top_genes_d2 <- specific_sig[order(-abs(get(COEF_COL)))][1:min(30, .N)]$gene

    d2_data <- all_sample_coefs[gene %in% top_genes_d2 & !is.na(coef)]

    if (nrow(d2_data) > 10) {
      # Build matrix: genes x samples (for first cell type with most genes)
      d2_wide <- dcast(d2_data, as.formula(paste("gene + cell_type ~", SAMPLE_COL_CSV)), value.var = "coef", fill = NA)

      # Heatmap for each cell type separately
      for (ct in unique(d2_wide$cell_type)[1:min(3, uniqueN(d2_wide$cell_type))]) {
        ct_data <- d2_wide[cell_type == ct]
        if (nrow(ct_data) < 3) next

        ct_mat <- as.matrix(ct_data[, -c(1,2), with = FALSE])
        rownames(ct_mat) <- ct_data$gene

        # Remove all-NA columns
        ct_mat <- ct_mat[, colSums(!is.na(ct_mat)) > 0, drop = FALSE]
        if (ncol(ct_mat) < 2) next

        max_abs <- max(abs(ct_mat), na.rm = TRUE)
        if (is.na(max_abs) || max_abs == 0) next
        breaks <- seq(-max_abs, max_abs, length.out = 101)

        pdf(file.path(OUTPUT_DIR, paste0("sample_contribution_heatmap_", ct, ".pdf")),
            width = max(6, ncol(ct_mat) * 1.2 + 2),
            height = max(5, nrow(ct_mat) * 0.25 + 2))

        pheatmap(
          ct_mat,
          color = DIVERGING_PALETTE(100),
          breaks = breaks,
          na_col = "grey90",
          cluster_rows = nrow(ct_mat) >= 2,
          cluster_cols = ncol(ct_mat) >= 2,
          main = paste0(ct, ": Per-Sample Coefficients (red = ", QUERY_LABEL, "-induced)"),
          fontsize_row = 8,
          fontsize_col = 10,
          angle_col = 45
        )

        dev.off()
        message("  Saved: sample_contribution_heatmap_", ct, ".pdf")
      }
    } else {
      message("  [SKIP] Too few data points for sample contribution heatmap")
    }
  } else {
    message("  [SKIP] Sample coefficient data missing sample ID column")
  }
} else {
  message("  [SKIP] Per-sample coefficient files not found")
}

# =============================================================================
# Panel E1: Ligand-Receptor Pair Highlighting
# =============================================================================
# Highlights known L-R pairs in the results to validate spatial communication.

message("Panel E1: Ligand-receptor highlighting...")

lr_results <- rbindlist(lapply(LR_PAIRS, function(pair) {
  ligand <- pair[1]
  receptor <- pair[2]

  cols_e1 <- c("gene", "cell_type", COEF_COL, SIG_COL, "specificity")
  cols_e1 <- intersect(cols_e1, names(all_results))
  ligand_data <- all_results[gene == ligand, ..cols_e1]
  receptor_data <- all_results[gene == receptor, ..cols_e1]

  if (nrow(ligand_data) == 0 && nrow(receptor_data) == 0) return(NULL)

  rbind(
    if (nrow(ligand_data) > 0) cbind(ligand_data, role = "ligand", pair = paste0(ligand, "-", receptor)) else NULL,
    if (nrow(receptor_data) > 0) cbind(receptor_data, role = "receptor", pair = paste0(ligand, "-", receptor)) else NULL
  )
}), fill = TRUE)

if (nrow(lr_results) > 0) {
  # Show only significant pairs
  lr_sig <- lr_results[get(SIG_COL) < FDR_THRESHOLD]

  if (nrow(lr_sig) > 0) {
    lr_sig[, neg_log10_fdr := -log10(pmax(get(SIG_COL), 1e-50))]
    coef_lim <- max(abs(lr_sig[[COEF_COL]]), na.rm = TRUE) * 1.1

    p_e1 <- ggplot(lr_sig, aes(x = cell_type, y = gene)) +
      geom_point(aes(size = neg_log10_fdr, color = .data[[COEF_COL]], shape = role)) +
      scale_color_gradientn(
        colors = DIVERGING_PALETTE(100),
        limits = c(-coef_lim, coef_lim),
        name = "Log-rate\ncoefficient"
      ) +
      scale_size_continuous(range = c(2, 6), name = expression(-log[10](FDR))) +
      scale_shape_manual(values = c("ligand" = 16, "receptor" = 17), name = "Role") +
      facet_wrap(~ pair, scales = "free_y", ncol = 2) +
      labs(
        title = "Known Ligand-Receptor Pairs in Spatial Gradients",
        subtitle = paste0("Only FDR < ", FDR_THRESHOLD, " shown | ",
                          "Red = induced near ", QUERY_LABEL),
        x = NULL, y = NULL
      ) +
      theme_bw(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        strip.text = element_text(face = "bold"),
        legend.position = "right"
      )

    ggsave(file.path(OUTPUT_DIR, "ligand_receptor_pairs.pdf"), p_e1,
           width = 12, height = max(6, length(unique(lr_sig$pair)) * 2))
    message("  Saved: ligand_receptor_pairs.pdf (",
            uniqueN(lr_sig$pair), " pairs with significant hits)")
  } else {
    message("  [SKIP] No L-R pairs reached significance")
  }

  # Save full L-R table
  fwrite(lr_results[order(pair, cell_type)],
         file.path(OUTPUT_DIR, "ligand_receptor_full_table.csv"))
  message("  Saved: ligand_receptor_full_table.csv")
} else {
  message("  [SKIP] No L-R pairs found in results")
}

# =============================================================================
# Panel E3: Query Signature Gene Leakage Check
# =============================================================================
# If query signature genes show gradients in OTHER cell types, this suggests
# segmentation artifacts (transcript leakage) rather than true expression.

message("Panel E3: Query signature leakage check...")

hymy_sig_results <- all_results[gene %in% HYMY_SIGNATURE_GENES]

if (nrow(hymy_sig_results) > 0) {
  # Exclude query cell type itself if present
  # Exclude the query cell type itself (and common myeloid aliases) from leakage check
  exclude_pattern <- paste(c(QUERY_CELLTYPE, CONTROL_CELLTYPE), collapse = "|")
  hymy_sig_results <- hymy_sig_results[!grepl(exclude_pattern, cell_type, ignore.case = TRUE)]

  if (nrow(hymy_sig_results) > 0) {
    # Heatmap of query signature genes in non-query cell types
    e3_wide <- dcast(hymy_sig_results, gene ~ cell_type, value.var = COEF_COL, fill = 0)
    e3_mat <- as.matrix(e3_wide[, -1])
    rownames(e3_mat) <- e3_wide$gene

    # FDR matrix for stars
    e3_fdr <- dcast(hymy_sig_results, gene ~ cell_type, value.var = SIG_COL, fill = 1)
    e3_fdr_mat <- as.matrix(e3_fdr[, -1])
    rownames(e3_fdr_mat) <- e3_fdr$gene
    e3_stars <- ifelse(e3_fdr_mat < FDR_THRESHOLD, "*", "")

    max_abs <- max(abs(e3_mat), na.rm = TRUE)
    if (max_abs > 0) {
      breaks <- seq(-max_abs, max_abs, length.out = 101)

      pdf(file.path(OUTPUT_DIR, "hymy_signature_leakage_check.pdf"),
          width = max(8, ncol(e3_mat) * 0.7 + 2),
          height = max(5, nrow(e3_mat) * 0.4 + 2))

      pheatmap(
        e3_mat,
        color = DIVERGING_PALETTE(100),
        breaks = breaks,
        display_numbers = e3_stars,
        fontsize_number = 10,
        cluster_rows = TRUE, cluster_cols = TRUE,
        main = paste0(QUERY_LABEL, " Signature Genes in OTHER Cell Types\n",
                      "(Red = detected near ", QUERY_LABEL, " = potential leakage)"),
        fontsize_row = 10,
        fontsize_col = 10,
        angle_col = 45
      )

      dev.off()
      message("  Saved: hymy_signature_leakage_check.pdf")

      # Count how many signature genes show leakage
      sig_leakage <- hymy_sig_results[get(SIG_COL) < FDR_THRESHOLD & get(COEF_COL) < 0]
      if (nrow(sig_leakage) > 0) {
        message("  WARNING: ", uniqueN(sig_leakage$gene), " query signature genes ",
                "show ", QUERY_LABEL, "-induced gradients in other cell types:")
        for (g in unique(sig_leakage$gene)) {
          cts <- paste(sig_leakage[gene == g]$cell_type, collapse = ", ")
          message("    ", g, ": ", cts)
        }
      } else {
        message("  OK: No significant query signature leakage detected")
      }
    }
  } else {
    message("  [SKIP] No non-myeloid cell types to check for leakage")
  }
} else {
  message("  [SKIP] No query signature genes found in results")
}

# =============================================================================
# Stage 2 Panels (if available)
# =============================================================================

if (HAS_STAGE2) {
  message("\n--- Stage 2 Visualizations ---")

  # -------------------------------------------------------------------------
  # Panel P10: Stage 2 Classification Breakdown
  # -------------------------------------------------------------------------
  message("Panel P10: Stage 2 classification breakdown...")

  class_by_ct <- stage2_results[, .N, by = .(cell_type, classification)]
  class_by_ct[, cell_type := factor(cell_type, levels = ct_order)]

  p10 <- ggplot(class_by_ct, aes(x = cell_type, y = N, fill = classification)) +
    geom_col(position = "stack", width = 0.7) +
    scale_fill_manual(values = CLASSIFICATION_COLORS, name = "Classification") +
    labs(
      title = paste0("Stage 2: ", QUERY_LABEL, "-Specific vs Niche-Driven Gene Classification"),
      subtitle = paste0(query_specific_label, " = gradient persists after controlling for ", CONTROL_CELLTYPE, " distance\n",
                        "niche_driven = gradient explained by tissue compartment, not ", QUERY_LABEL, " specifically"),
      x = NULL,
      y = "Number of Stage 1 significant genes"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      legend.position = "top"
    )

  ggsave(file.path(OUTPUT_DIR, "stage2_classification_breakdown.pdf"), p10,
         width = 10, height = 6)
  message("  Saved: stage2_classification_breakdown.pdf")

  # -------------------------------------------------------------------------
  # Panel P11: Stage 1 vs Stage 2 Coefficient Scatter
  # -------------------------------------------------------------------------
  message("Panel P11: Stage 1 vs Stage 2 coefficient scatter...")

  scatter_data <- stage2_results[!is.na(stage1_coef) & !is.na(stage2_coef)]

  if (nrow(scatter_data) > 0) {
    coef_lim <- max(abs(c(scatter_data$stage1_coef, scatter_data$stage2_coef)), na.rm = TRUE) * 1.1

    p11 <- ggplot(scatter_data, aes(x = stage1_coef, y = stage2_coef)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
      geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
      geom_point(aes(color = classification), alpha = 0.5, size = 1) +
      scale_color_manual(values = CLASSIFICATION_COLORS, name = "Classification") +
      coord_fixed(xlim = c(-coef_lim, coef_lim), ylim = c(-coef_lim, coef_lim)) +
      facet_wrap(~ cell_type, ncol = 4) +
      labs(
        title = "Stage 1 vs Stage 2 Coefficients",
        subtitle = paste0("Diagonal = no change after ", CONTROL_CELLTYPE, " control | ",
                          "Points below diagonal = ", QUERY_LABEL, " effect attenuated"),
        x = paste0("Stage 1 coefficient (univariate: dist_to_", QUERY_LABEL, " only)"),
        y = paste0("Stage 2 coefficient (bivariate: dist_to_", QUERY_LABEL, " + dist_to_", CONTROL_CELLTYPE, ")")
      ) +
      theme_bw(base_size = 9) +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        strip.text = element_text(face = "bold"),
        legend.position = "top"
      )

    ggsave(file.path(OUTPUT_DIR, "stage1_vs_stage2_scatter.pdf"), p11,
           width = 14, height = 3.5 * ceiling(n_ct / 4))
    message("  Saved: stage1_vs_stage2_scatter.pdf")
  }

  # -------------------------------------------------------------------------
  # Stage 2: Query-Specific Genes Summary Table
  # -------------------------------------------------------------------------
  hymy_specific_genes <- stage2_results[classification == query_specific_label][
    order(cell_type, stage2_fdr)]
  if (nrow(hymy_specific_genes) > 0) {
    fwrite(hymy_specific_genes, file.path(OUTPUT_DIR, "stage2_hymy_specific_genes.csv"))
    message("  Saved: stage2_hymy_specific_genes.csv (",
            nrow(hymy_specific_genes), " gene-celltype pairs)")
  }

  niche_driven_genes <- stage2_results[classification == "niche_driven"][
    order(cell_type, stage1_fdr)]
  if (nrow(niche_driven_genes) > 0) {
    fwrite(niche_driven_genes, file.path(OUTPUT_DIR, "stage2_niche_driven_genes.csv"))
    message("  Saved: stage2_niche_driven_genes.csv (",
            nrow(niche_driven_genes), " gene-celltype pairs)")
  }
}

# =============================================================================
# Panel P12: Sign Consistency Distribution (v2 diagnostic)
# =============================================================================
# Stacked proportional bar per cell type: fraction of significant genes by
# sign agreement across samples

if ("sign_consistency" %in% names(all_results)) {
  message("\nPanel P12: Sign consistency distribution...")

  sign_data <- sig_results[!is.na(sign_consistency)]

  if (nrow(sign_data) > 0) {
    sign_data[, consistency_bin := fifelse(
      sign_consistency == 1.0, "All agree (1.0)",
      fifelse(sign_consistency >= 0.75, "Majority (0.75-1.0)",
              "Split (<0.75)")
    )]
    sign_data[, consistency_bin := factor(consistency_bin,
      levels = c("All agree (1.0)", "Majority (0.75-1.0)", "Split (<0.75)"))]

    sign_summary <- sign_data[, .N, by = .(cell_type, consistency_bin)]
    sign_summary[, total := sum(N), by = cell_type]
    sign_summary[, fraction := N / total]

    p12 <- ggplot(sign_summary, aes(x = cell_type, y = fraction, fill = consistency_bin)) +
      geom_col(position = "stack", width = 0.7) +
      scale_fill_manual(values = c(
        "All agree (1.0)" = "#2166AC",
        "Majority (0.75-1.0)" = "#92C5DE",
        "Split (<0.75)" = "#F4A582"
      ), name = "Sign consistency") +
      scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
      labs(
        title = "Sign Consistency of Significant Genes Across Samples",
        subtitle = "Fraction of samples where coefficient sign matches the meta-analysis direction",
        x = NULL, y = "Fraction of significant genes"
      ) +
      theme_bw(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"),
        legend.position = "top"
      )

    ggsave(file.path(OUTPUT_DIR, "sign_consistency_distribution.pdf"), p12,
           width = 10, height = 5)
    message("  Saved: sign_consistency_distribution.pdf")
    message("  Genes with split sign (<0.75): ",
            nrow(sign_data[consistency_bin == "Split (<0.75)"]), " / ", nrow(sign_data))
  } else {
    message("  [SKIP] No significant genes with sign_consistency data")
  }
} else {
  message("\nPanel P12: [SKIP] sign_consistency column not found (v1 results)")
}


# =============================================================================
# Panel P13: Overdispersion QC (v2 diagnostic)
# =============================================================================
# Boxplot of median_dispersion per cell type for significant genes.
# Poisson expectation = 1.0; values >> 1 indicate overdispersion.

if ("median_dispersion" %in% names(all_results)) {
  message("\nPanel P13: Overdispersion QC...")

  disp_data <- sig_results[!is.na(median_dispersion)]

  if (nrow(disp_data) > 0) {
    p13 <- ggplot(disp_data, aes(x = cell_type, y = median_dispersion)) +
      geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3, fill = "#D1E5F0") +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "firebrick", linewidth = 0.5) +
      scale_y_log10() +
      annotation_logticks(sides = "l") +
      labs(
        title = "Overdispersion of Significant Genes (Poisson Residuals)",
        subtitle = "Dashed red line = Poisson expectation (dispersion = 1.0)",
        x = NULL, y = "Median dispersion (log scale)"
      ) +
      theme_bw(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold")
      )

    ggsave(file.path(OUTPUT_DIR, "overdispersion_qc.pdf"), p13,
           width = 10, height = 5)
    message("  Saved: overdispersion_qc.pdf")
    message("  Median dispersion across all sig genes: ",
            round(median(disp_data$median_dispersion, na.rm = TRUE), 2))
  } else {
    message("  [SKIP] No significant genes with median_dispersion data")
  }
} else {
  message("\nPanel P13: [SKIP] median_dispersion column not found (v1 results)")
}


# =============================================================================
# Panel B2: Induced vs Repressed Pathway Split
# =============================================================================
# Separate fgsea for genes with negative vs positive coefficients

message("\nPanel B2: Induced vs repressed pathway enrichment...")

if (HAS_FGSEA) {
  # Split genes by direction
  induced_genes <- all_results[get(SIG_COL) < FDR_THRESHOLD & get(COEF_COL) < 0, unique(gene)]
  repressed_genes <- all_results[get(SIG_COL) < FDR_THRESHOLD & get(COEF_COL) > 0, unique(gene)]

  # Load Hallmark
  hallmark_df <- msigdbr(species = "Mus musculus", category = "H")
  pathways <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)

  clean_pathway_name <- function(x) {
    x <- sub("^HALLMARK_", "", x)
    x <- gsub("_", " ", x)
    tools::toTitleCase(tolower(x))
  }

  # Hypergeometric test for enrichment
  run_enrichment <- function(genes, universe, pathways, direction_label) {
    rbindlist(lapply(names(pathways), function(pw) {
      pw_genes <- intersect(pathways[[pw]], universe)
      if (length(pw_genes) < 5) return(NULL)

      overlap <- length(intersect(genes, pw_genes))
      if (overlap < 2) return(NULL)

      # Fisher's exact test
      a <- overlap
      b <- length(genes) - overlap
      c <- length(pw_genes) - overlap
      d <- length(universe) - length(genes) - c

      if (any(c(a, b, c, d) < 0)) return(NULL)

      test <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")

      data.table(
        pathway = pw,
        pathway_clean = clean_pathway_name(pw),
        direction = direction_label,
        overlap = overlap,
        pw_size = length(pw_genes),
        pval = test$p.value,
        odds_ratio = test$estimate
      )
    }))
  }

  universe <- unique(all_results$gene)
  enrich_induced <- run_enrichment(induced_genes, universe, pathways, induced_label)
  enrich_repressed <- run_enrichment(repressed_genes, universe, pathways, repressed_label)

  enrich_both <- rbind(enrich_induced, enrich_repressed)

  if (nrow(enrich_both) > 0) {
    enrich_both[, padj := p.adjust(pval, method = "BH")]
    enrich_both <- enrich_both[order(pval)]

    fwrite(enrich_both, file.path(OUTPUT_DIR, "pathway_enrichment_by_direction.csv"))
    message("  Saved: pathway_enrichment_by_direction.csv")

    # Plot top pathways
    plot_data <- enrich_both[padj < 0.1 | pval < 0.05][
      order(direction, pval), head(.SD, 15), by = direction]

    if (nrow(plot_data) > 0) {
      plot_data[, neg_log10_p := -log10(pmax(pval, 1e-20))]

      p_b2 <- ggplot(plot_data, aes(x = neg_log10_p, y = reorder(pathway_clean, neg_log10_p))) +
        geom_col(aes(fill = direction), width = 0.7) +
        geom_text(aes(label = overlap), hjust = -0.2, size = 3) +
        scale_fill_manual(values = DIRECTION_COLORS, name = NULL) +
        facet_wrap(~ direction, scales = "free", ncol = 2) +
        labs(
          title = paste0("Pathway Enrichment: ", QUERY_LABEL, "-Induced vs ", QUERY_LABEL, "-Repressed Genes"),
          subtitle = "Fisher's exact test | Numbers = overlapping genes",
          x = expression(-log[10](p-value)),
          y = NULL
        ) +
        theme_bw(base_size = 10) +
        theme(
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          strip.text = element_text(face = "bold"),
          legend.position = "none"
        )

      ggsave(file.path(OUTPUT_DIR, "pathway_enrichment_by_direction.pdf"), p_b2,
             width = 12, height = max(6, nrow(plot_data) * 0.3))
      message("  Saved: pathway_enrichment_by_direction.pdf")
    }
  }
} else {
  message("  [SKIP] fgsea not available for B2")
}

# =============================================================================
# Panels 8-9: fgsea Pathway Enrichment (Hallmark)
# =============================================================================
# Ranking metric: median_coef (median of per-sample Poisson coefficients).
# This treats each mouse as an equally weighted independent test set, avoiding
# the precision-weighting of meta-analysis TE.random which can let one cell-rich
# mouse dominate the ranking.
#
# Both rankings (median_coef and combined_coef) are run for comparison.
# Primary plots and tables use median_coef; comparison table saved separately.
#
# NES < 0 → pathway genes enriched among query-induced genes (higher near query)
# NES > 0 → pathway genes enriched among query-repressed genes (lower near query)
#
# NOTE: The Xenium panel covers ~5K genes, not the full transcriptome. GSEA
# results should be interpreted cautiously — pathway coverage varies.

# Helper: run fgsea across cell types with a given ranking column
run_fgsea_per_celltype <- function(data, gene_sets, rank_col, min_genes = 100,
                                   min_size = 10, max_size = 500,
                                   clean_name_fn = identity, label = rank_col) {
  rbindlist(lapply(levels(droplevels(data$cell_type)), function(ct) {
    ct_data <- data[cell_type == ct & !is.na(get(rank_col))]
    stats <- setNames(ct_data[[rank_col]], ct_data$gene)
    stats <- sort(stats)
    if (length(stats) < min_genes) {
      message("    ", ct, ": too few genes (", length(stats), "), skipping [", label, "]")
      return(NULL)
    }
    res <- fgsea(pathways = gene_sets, stats = stats,
                 minSize = min_size, maxSize = max_size)
    res[, cell_type := ct]
    if (!is.null(clean_name_fn)) {
      res[, pathway_clean := clean_name_fn(pathway)]
    }
    res
  }), fill = TRUE)
}

p8 <- NULL
if (HAS_FGSEA) {
  message("\nPanel 8-9: fgsea pathway enrichment (Hallmark)...")

  # Load Hallmark gene sets for mouse
  hallmark_df <- msigdbr(species = "Mus musculus", category = "H")
  pathways <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)

  # Clean pathway names for display
  clean_pathway_name <- function(x) {
    x <- sub("^HALLMARK_", "", x)
    x <- gsub("_", " ", x)
    x <- tools::toTitleCase(tolower(x))
    x
  }

  # Primary: rank by median_coef (equal mouse weighting)
  set.seed(42)
  message("  Running fgsea with median_coef ranking...")
  fgsea_all <- run_fgsea_per_celltype(
    all_results, pathways, "median_coef",
    clean_name_fn = clean_pathway_name, label = "median_coef"
  )

  # Comparison: rank by combined_coef (meta-analysis TE.random)
  message("  Running fgsea with combined_coef ranking (comparison)...")
  fgsea_meta <- run_fgsea_per_celltype(
    all_results, pathways, "combined_coef",
    clean_name_fn = clean_pathway_name, label = "combined_coef"
  )

  # Save ranking comparison table
  if (nrow(fgsea_all) > 0 && nrow(fgsea_meta) > 0) {
    comparison <- merge(
      fgsea_all[, .(cell_type, pathway, pathway_clean,
                     NES_median = NES, padj_median = padj)],
      fgsea_meta[, .(cell_type, pathway,
                      NES_meta = NES, padj_meta = padj)],
      by = c("cell_type", "pathway"), all = TRUE
    )
    fwrite(comparison, file.path(OUTPUT_DIR, "fgsea_hallmark_ranking_comparison.csv"))
    # Summarize agreement
    both_sig <- comparison[padj_median < 0.05 & padj_meta < 0.05]
    median_only <- comparison[padj_median < 0.05 & (padj_meta >= 0.05 | is.na(padj_meta))]
    meta_only <- comparison[padj_meta < 0.05 & (padj_median >= 0.05 | is.na(padj_median))]
    message(sprintf("  Ranking comparison (Hallmark): %d both sig, %d median-only, %d meta-only",
                    nrow(both_sig), nrow(median_only), nrow(meta_only)))
  }

  if (nrow(fgsea_all) > 0) {
    # Save per-cell-type and combined tables
    # (leadingEdge is a list column — convert to string for CSV export)
    fgsea_export <- copy(fgsea_all)
    fgsea_export[, leadingEdge := sapply(leadingEdge, paste, collapse = ",")]

    for (ct in unique(fgsea_export$cell_type)) {
      ct_res <- fgsea_export[cell_type == ct][order(pval)]
      fwrite(ct_res[, .(pathway, pathway_clean, pval, padj, ES, NES, size, leadingEdge)],
             file.path(OUTPUT_DIR, paste0("fgsea_hallmark_", ct, ".csv")))
    }

    fwrite(fgsea_export[, .(cell_type, pathway, pathway_clean, pval, padj, ES, NES, size, leadingEdge)],
           file.path(OUTPUT_DIR, "fgsea_hallmark_all_celltypes.csv"))
    message("  Saved: fgsea_hallmark_all_celltypes.csv + per-cell-type tables")

    # --- Panel 8: NES Heatmap ---
    sig_pathways <- fgsea_all[padj < 0.05, unique(pathway)]

    if (length(sig_pathways) >= 2) {
      # Build NES matrix
      nes_wide <- dcast(fgsea_all[pathway %in% sig_pathways],
                        pathway_clean ~ cell_type, value.var = "NES", fill = 0)
      nes_mat <- as.matrix(nes_wide[, -1])
      rownames(nes_mat) <- nes_wide$pathway_clean

      # Significance stars
      padj_wide <- dcast(fgsea_all[pathway %in% sig_pathways],
                         pathway_clean ~ cell_type, value.var = "padj", fill = 1)
      padj_mat <- as.matrix(padj_wide[, -1])
      rownames(padj_mat) <- padj_wide$pathway_clean
      star_mat <- ifelse(padj_mat < 0.05, "*", "")

      # Diverging color limits
      max_nes <- max(abs(nes_mat), na.rm = TRUE)
      breaks <- seq(-max_nes, max_nes, length.out = 101)

      pdf(file.path(OUTPUT_DIR, "fgsea_hallmark_heatmap.pdf"),
          width = max(9, ncol(nes_mat) * 0.8 + 4),
          height = max(6, nrow(nes_mat) * 0.35 + 2))

      pheatmap(
        nes_mat,
        color = DIVERGING_PALETTE(100),
        breaks = breaks,
        display_numbers = star_mat,
        fontsize_number = 10,
        cluster_rows = TRUE, cluster_cols = TRUE,
        main = paste0("Hallmark Pathway Enrichment (NES)\n",
                      "Red = induced near ", QUERY_LABEL,
                      " | Blue = repressed near ", QUERY_LABEL),
        fontsize_row = 9,
        fontsize_col = 10,
        angle_col = 45
      )

      dev.off()
      message("  Saved: fgsea_hallmark_heatmap.pdf (", length(sig_pathways), " pathways)")
    } else {
      message("  [SKIP] Fewer than 2 significant pathways for heatmap")
    }

    # --- Panel 9: Pathway Dot Plot (significant pathways only) ---
    # Show pathways significant (padj < 0.05) in at least one cell type.

    sig_pw_names <- fgsea_all[padj < 0.05, unique(pathway)]
    dotplot_data <- fgsea_all[pathway %in% sig_pw_names]
    dotplot_data[, neg_log10_padj := -log10(pmax(padj, 1e-20))]
    dotplot_data[, is_sig := padj < 0.05]

    # Order pathways by mean NES (induced-near-query at bottom)
    pw_order <- dotplot_data[, .(mean_nes = mean(NES, na.rm = TRUE)),
                              by = pathway_clean][order(mean_nes)]$pathway_clean
    dotplot_data[, pathway_clean := factor(pathway_clean, levels = pw_order)]

    # NES limits for symmetric color scale
    nes_lim <- max(abs(dotplot_data$NES), na.rm = TRUE)

    if (nrow(dotplot_data[is_sig == TRUE]) > 0) {
      p8 <- ggplot(dotplot_data[is_sig == TRUE], aes(x = cell_type, y = pathway_clean)) +
        geom_point(aes(size = neg_log10_padj, color = NES)) +
        scale_color_gradientn(
          colors = DIVERGING_PALETTE(100),
          limits = c(-nes_lim, nes_lim),
          name = "NES"
        ) +
        scale_size_continuous(
          range = c(1.5, 6),
          name = expression(-log[10](padj)),
          breaks = c(2, 5, 10, 20)
        ) +
        labs(
          title = "Hallmark Pathway Enrichment per Cell Type",
          subtitle = paste0("Pathways significant in \u22651 cell type (padj < 0.05)\n",
                            "Red = induced near ", QUERY_LABEL,
                            " | Blue = repressed near ", QUERY_LABEL),
          x = NULL, y = NULL
        ) +
        theme_bw(base_size = 10) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y = element_text(size = 9),
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          panel.grid.major = element_line(color = "grey92"),
          legend.key.size = unit(0.4, "cm")
        )

      n_pw <- length(sig_pw_names)
      ggsave(file.path(OUTPUT_DIR, "fgsea_hallmark_dotplot.pdf"), p8,
             width = max(7, n_ct * 0.5 + 3),
             height = max(5, n_pw * 0.3 + 2))
      message("  Saved: fgsea_hallmark_dotplot.pdf (", n_pw, " sig pathways x ", n_ct, " cell types)")
    } else {
      message("  [SKIP] No significant pathways for dot plot")
    }

    # --- Summary ---
    message("\n  Pathway enrichment summary:")
    for (ct in levels(droplevels(all_results$cell_type))) {
      if (ct %in% fgsea_all$cell_type) {
        n_sig <- sum(fgsea_all[cell_type == ct]$padj < 0.05, na.rm = TRUE)
        top_pw <- fgsea_all[cell_type == ct & padj < 0.05][order(pval)]$pathway_clean[1]
        if (!is.na(top_pw)) {
          message(sprintf("    %-18s  %2d sig pathways (top: %s)", ct, n_sig, top_pw))
        } else {
          message(sprintf("    %-18s  %2d sig pathways", ct, n_sig))
        }
      }
    }
  }
} else {
  message("\n[SKIP] Panels 8-9: fgsea not available")
}

# =============================================================================
# Panels 8b-9b: Custom Gene Module Enrichment
# =============================================================================
# Curated gene modules relevant to spatial gradient biology. Unlike Hallmark pathways,
# these are small (5-18 genes), focused, and have known panel coverage.
# Uses the same fgsea machinery but with custom gene sets.

if (HAS_FGSEA) {
  message("\nPanel 8b-9b: Custom gene module enrichment...")

  CUSTOM_MODULES <- list(
    # --- Signaling ---
    inflammatory_cytokines = c("Il1b","Il1a","Il6","Tnf","Il12a","Il12b","Il18","Il33","Il1rn","Ifng"),
    anti_inflammatory = c("Il10","Tgfb1","Tgfb2","Tgfb3","Il1rn","Socs1","Socs3","Arg1","Ido1"),
    nfkb_targets = c("Nfkbia","Nfkbiz","Tnfaip3","Bcl2l1","Birc3","Ccl2","Ccl5","Cxcl1","Cxcl2","Cxcl10","Icam1","Vcam1","Il6","Tnf","Ptgs2","Sod2"),
    stat3_targets = c("Socs3","Socs1","Bcl2","Bcl2l1","Mcl1","Myc","Ccnd1","Vegfa","Mmp2","Mmp9"),
    interferon_response = c("Irf1","Irf3","Irf7","Irf9","Stat1","Stat2","Isg15","Isg20","Oasl2","Ifit1","Ifit3","Ifitm1","Mx1","Mx2","Oas1a","Oas2","Oasl1","Ifit2","Ifitm3","Bst2"),

    # --- Chemokines ---
    cc_chemokines = c("Ccl2","Ccl3","Ccl4","Ccl5","Ccl6","Ccl7","Ccl8","Ccl9","Ccl11","Ccl12","Ccl17","Ccl19","Ccl20","Ccl21a","Ccl22","Ccl24","Ccl25","Ccl28"),
    cxc_chemokines = c("Cxcl1","Cxcl2","Cxcl3","Cxcl4","Cxcl5","Cxcl9","Cxcl10","Cxcl11","Cxcl12","Cxcl13","Cxcl14","Cxcl16"),
    chemokine_receptors = c("Ccr1","Ccr2","Ccr3","Ccr5","Ccr6","Ccr7","Ccr9","Cxcr2","Cxcr3","Cxcr4","Cxcr5","Cxcr6","Cx3cr1","Xcr1","Ackr1","Ackr2","Ackr3","Ackr4"),

    # --- Growth factors / Lymphangiogenesis ---
    csf_signaling = c("Csf1","Csf2","Csf3","Csf1r","Csf2ra","Csf2rb","Csf2rb2","Csf3r"),
    vegf_lymphangiogenesis = c("Vegfa","Vegfb","Vegfc","Vegfd","Flt1","Kdr","Flt4","Nrp1","Nrp2","Prox1","Lyve1","Pdpn"),

    # --- T cell states ---
    tcell_exhaustion = c("Pdcd1","Lag3","Havcr2","Tigit","Ctla4","Tox","Entpd1","Cd38","Cd244a","Cd48"),
    tcell_naive_memory = c("Sell","Ccr7","Tcf7","Lef1","Il7r","Bach2","Id3"),

    # --- Antigen presentation / Immune checkpoint ---
    antigen_presentation = c("H2-Aa","H2-Ab1","H2-Eb1","Cd74","Cd80","Cd86","Cd40","Tap1","Tap2","B2m","Ciita"),
    immune_checkpoint = c("Cd274","Pdcd1lg2","Cd80","Cd86","Cd276","Lgals9","Ido1","Arg1"),

    # --- Myeloid polarization ---
    myeloid_proinflam = c("Nos2","Il1b","Il6","Tnf","Cxcl9","Cxcl10","Cd86","Irf1","Irf5"),
    myeloid_antiinflam = c("Arg1","Mrc1","Cd163","Retnla","Chil3","Il10","Tgfb1"),
    phagocytosis = c("Mertk","Axl","Tyro3","Cd36","Msr1","Cd68","Fcgr1","Fcgr3","Marco"),
    complement = c("C1qa","C1qb","C1qc","C3","C4b","C3ar1","C5ar1","Cfb","Cfh","Cfd","Serping1"),

    # --- Tissue remodeling ---
    ecm_remodeling = c("Col1a1","Col1a2","Col3a1","Col4a1","Col4a2","Fn1","Lama4","Lamb1","Mmp2","Mmp9","Mmp14","Timp1","Timp2","Acta2","Fap"),

    # --- Cell fate ---
    proliferation = c("Mki67","Top2a","Pcna","Bub1b","Ccnb1","Ccnb2","Ccna2","Cdk1","Cdc20","Aurka","Aurkb","Birc5"),
    apoptosis_pro = c("Bax","Bak1","Bad","Bbc3","Casp3","Casp7","Casp8","Casp9","Cycs"),
    apoptosis_anti = c("Bcl2","Bcl2l1","Mcl1","Birc5","Xiap"),

    # --- TNF / IL receptor families ---
    tnf_superfamily = c("Tnf","Lta","Ltb","Fasl","Cd40lg","Tnfsf4","Tnfsf9","Tnfsf10","Tnfsf11","Tnfsf13b","Tnfsf14"),
    tnf_receptors = c("Tnfrsf1a","Tnfrsf1b","Fas","Cd40","Tnfrsf4","Tnfrsf9","Tnfrsf11a","Tnfrsf11b","Tnfrsf13b","Tnfrsf13c","Tnfrsf14","Tnfrsf17"),
    il_receptors = c("Il1r1","Il1r2","Il1rl1","Il2ra","Il2rb","Il2rg","Il4ra","Il6ra","Il6st","Il7r","Il10ra","Il10rb","Il15ra","Il17ra","Il18r1","Il21r","Il27ra","Il31ra"),

    # --- Metabolic ---
    glycolysis = c("Hk1","Hk2","Pfkfb3","Pkm","Ldha","Slc2a1","Eno1","Gapdh","Pgk1","Aldoa")

    # NOTE: query_signature module deliberately excluded -- it reflects segmentation
    # leakage (query transcripts in neighbors), not biology. See Panel E3 instead.
  )

  # Filter modules to only include genes present in the results
  available_genes <- unique(all_results$gene)
  custom_modules_filtered <- lapply(CUSTOM_MODULES, function(genes) {
    genes[genes %in% available_genes]
  })

  # Report coverage
  message("  Module coverage on Xenium panel:")
  for (name in names(CUSTOM_MODULES)) {
    n_total <- length(CUSTOM_MODULES[[name]])
    n_avail <- length(custom_modules_filtered[[name]])
    pct <- round(100 * n_avail / n_total)
    if (n_avail < 3) {
      message(sprintf("    %-26s  %2d / %2d  (%3d%%)  [DROPPED: <3 genes]", name, n_avail, n_total, pct))
    } else {
      message(sprintf("    %-26s  %2d / %2d  (%3d%%)", name, n_avail, n_total, pct))
    }
  }

  # Drop modules with fewer than 3 genes on panel
  custom_modules_filtered <- custom_modules_filtered[sapply(custom_modules_filtered, length) >= 3]
  message("  Modules passing filter (>=3 genes): ", length(custom_modules_filtered))

  # Clean module names helper
  clean_module_name <- function(x) {
    tools::toTitleCase(gsub("_", " ", x))
  }

  # Primary: rank by median_coef
  message("  Running custom module fgsea with median_coef ranking...")
  custom_fgsea_all <- run_fgsea_per_celltype(
    all_results, custom_modules_filtered, "median_coef",
    min_genes = 50, min_size = 3,
    clean_name_fn = clean_module_name, label = "median_coef"
  )
  if (nrow(custom_fgsea_all) > 0) {
    custom_fgsea_all[, module_clean := clean_module_name(pathway)]
  }

  # Comparison: rank by combined_coef
  message("  Running custom module fgsea with combined_coef ranking (comparison)...")
  custom_fgsea_meta <- run_fgsea_per_celltype(
    all_results, custom_modules_filtered, "combined_coef",
    min_genes = 50, min_size = 3,
    clean_name_fn = clean_module_name, label = "combined_coef"
  )

  # Save ranking comparison table
  if (nrow(custom_fgsea_all) > 0 && nrow(custom_fgsea_meta) > 0) {
    custom_comparison <- merge(
      custom_fgsea_all[, .(cell_type, pathway,
                           NES_median = NES, padj_median = padj)],
      custom_fgsea_meta[, .(cell_type, pathway,
                            NES_meta = NES, padj_meta = padj)],
      by = c("cell_type", "pathway"), all = TRUE
    )
    fwrite(custom_comparison, file.path(OUTPUT_DIR, "fgsea_custom_modules_ranking_comparison.csv"))
    both_sig <- custom_comparison[padj_median < 0.05 & padj_meta < 0.05]
    median_only <- custom_comparison[padj_median < 0.05 & (padj_meta >= 0.05 | is.na(padj_meta))]
    meta_only <- custom_comparison[padj_meta < 0.05 & (padj_median >= 0.05 | is.na(padj_median))]
    message(sprintf("  Ranking comparison (custom): %d both sig, %d median-only, %d meta-only",
                    nrow(both_sig), nrow(median_only), nrow(meta_only)))
  }

  if (nrow(custom_fgsea_all) > 0) {
    # Save results
    custom_export <- copy(custom_fgsea_all)
    custom_export[, leadingEdge := sapply(leadingEdge, paste, collapse = ",")]
    fwrite(custom_export[, .(cell_type, pathway, module_clean, pval, padj, ES, NES, size, leadingEdge)],
           file.path(OUTPUT_DIR, "fgsea_custom_modules_all_celltypes.csv"))
    message("  Saved: fgsea_custom_modules_all_celltypes.csv")

    # --- Panel 8b: Custom Module NES Heatmap ---
    # Show ALL modules (not just significant) — use padj stars for significance
    modules_tested <- unique(custom_fgsea_all$pathway)

    if (length(modules_tested) >= 2) {
      nes_wide <- dcast(custom_fgsea_all, module_clean ~ cell_type,
                        value.var = "NES", fill = 0)
      nes_mat <- as.matrix(nes_wide[, -1])
      rownames(nes_mat) <- nes_wide$module_clean

      padj_wide <- dcast(custom_fgsea_all, module_clean ~ cell_type,
                         value.var = "padj", fill = 1)
      padj_mat <- as.matrix(padj_wide[, -1])
      rownames(padj_mat) <- padj_wide$module_clean
      star_mat <- ifelse(padj_mat < 0.001, "***",
                  ifelse(padj_mat < 0.01, "**",
                  ifelse(padj_mat < 0.05, "*", "")))

      max_nes <- max(abs(nes_mat), na.rm = TRUE)
      breaks <- seq(-max_nes, max_nes, length.out = 101)

      pdf(file.path(OUTPUT_DIR, "fgsea_custom_modules_heatmap.pdf"),
          width = max(10, ncol(nes_mat) * 0.8 + 6),
          height = max(8, nrow(nes_mat) * 0.35 + 3))

      pheatmap(
        nes_mat,
        color = DIVERGING_PALETTE(100),
        breaks = breaks,
        display_numbers = star_mat,
        fontsize_number = 8,
        cluster_rows = TRUE, cluster_cols = TRUE,
        main = paste0("Custom Gene Module Enrichment (NES)\n",
                      "Red = induced near ", QUERY_LABEL,
                      " | Blue = repressed near ", QUERY_LABEL,
                      "\n* p<0.05  ** p<0.01  *** p<0.001"),
        fontsize_row = 8,
        fontsize_col = 9,
        angle_col = 45
      )

      dev.off()
      message("  Saved: fgsea_custom_modules_heatmap.pdf (",
              length(modules_tested), " modules x ", n_ct, " cell types)")
    }

    # --- Panel 9b: Custom Module Dot Plot (significant only) ---
    sig_modules <- custom_fgsea_all[padj < 0.05, unique(pathway)]
    dotplot_custom <- custom_fgsea_all[pathway %in% sig_modules]
    dotplot_custom[, neg_log10_padj := -log10(pmax(padj, 1e-20))]
    dotplot_custom[, is_sig := padj < 0.05]

    if (nrow(dotplot_custom[is_sig == TRUE]) > 0) {
      # Order modules by mean NES
      mod_order <- dotplot_custom[, .(mean_nes = mean(NES, na.rm = TRUE)),
                                   by = module_clean][order(mean_nes)]$module_clean
      dotplot_custom[, module_clean := factor(module_clean, levels = mod_order)]

      nes_lim <- max(abs(dotplot_custom$NES), na.rm = TRUE)

      # Order cell types alphabetically
      all_ct <- sort(unique(as.character(dotplot_custom$cell_type)))
      ct_order <- all_ct
      dotplot_custom[, cell_type := factor(cell_type, levels = ct_order)]

      p_custom_dot <- ggplot(dotplot_custom[is_sig == TRUE],
                             aes(x = cell_type, y = module_clean)) +
        geom_point(aes(size = neg_log10_padj, color = NES)) +
        scale_color_gradientn(
          colors = DIVERGING_PALETTE(100),
          limits = c(-nes_lim, nes_lim),
          name = "NES"
        ) +
        scale_size_continuous(
          range = c(1.5, 5),
          name = expression(-log[10](padj)),
          breaks = c(2, 5, 10)
        ) +
        labs(
          title = "Custom Gene Module Enrichment per Cell Type",
          subtitle = paste0("Modules significant in >=1 cell type (padj < 0.05)\n",
                            "Red = induced near ", QUERY_LABEL,
                            " | Blue = repressed near ", QUERY_LABEL),
          x = NULL, y = NULL
        ) +
        theme_bw(base_size = 10) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y = element_text(size = 9),
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          panel.grid.major = element_line(color = "grey92"),
          legend.key.size = unit(0.4, "cm")
        )

      n_mod <- length(sig_modules)
      ggsave(file.path(OUTPUT_DIR, "fgsea_custom_modules_dotplot.pdf"), p_custom_dot,
             width = max(7, n_ct * 0.5 + 3),
             height = max(5, n_mod * 0.35 + 2.5))
      message("  Saved: fgsea_custom_modules_dotplot.pdf (",
              n_mod, " sig modules x ", n_ct, " cell types)")
    } else {
      message("  [SKIP] No significant custom modules for dot plot")
    }

    # --- Panel 10b: Leading Edge Heatmap for Key Cell Types ---
    # For key cell types: show leading edge genes
    # of significant modules as a heatmap of gradient scores.
    # Use top 5 cell types by number of significant modules, or all if fewer
    le_candidates <- custom_fgsea_all[padj < 0.05, .(n_sig = .N), by = cell_type][order(-n_sig)]
    LEADING_EDGE_CELLTYPES <- head(le_candidates$cell_type, 5)
    if (length(LEADING_EDGE_CELLTYPES) == 0) LEADING_EDGE_CELLTYPES <- unique(custom_fgsea_all$cell_type)[1:min(5, uniqueN(custom_fgsea_all$cell_type))]
    message("\nPanel 10b: Leading edge heatmap for key cell types...")

    for (ct in LEADING_EDGE_CELLTYPES) {
      ct_sig <- custom_fgsea_all[cell_type == ct & padj < 0.05]
      if (nrow(ct_sig) == 0) {
        message("    ", ct, ": no significant modules, skipping")
        next
      }

      # Extract leading edge genes per module
      le_data <- rbindlist(lapply(seq_len(nrow(ct_sig)), function(i) {
        mod <- ct_sig$pathway[i]
        mod_clean <- ct_sig$module_clean[i]
        nes <- ct_sig$NES[i]
        padj_val <- ct_sig$padj[i]
        le_genes <- if (is.list(ct_sig$leadingEdge)) {
          unlist(ct_sig$leadingEdge[[i]])
        } else {
          trimws(unlist(strsplit(as.character(ct_sig$leadingEdge[i]), ",")))
        }
        data.table(module = mod, module_clean = mod_clean,
                   gene = le_genes, module_nes = nes, module_padj = padj_val)
      }))

      if (nrow(le_data) == 0) next

      # Merge with gradient scores for this cell type
      ct_results <- all_results[cell_type == ct, .(gene, median_coef = get(COEF_COL),
                                                      sig_fdr = get(SIG_COL))]
      le_data <- merge(le_data, ct_results, by = "gene", all.x = TRUE)

      # Order modules by NES (induced at top, repressed at bottom)
      mod_order <- ct_sig[order(NES)]$module_clean
      le_data[, module_clean := factor(module_clean, levels = mod_order)]

      # Order genes within each module by coefficient
      le_data[, gene := factor(gene, levels = unique(le_data[order(module_clean, median_coef)]$gene))]

      # Significance markers
      le_data[, sig_label := fifelse(!is.na(sig_fdr) & sig_fdr < 0.05, "*", "")]

      # Coefficient limits for symmetric color scale
      coef_lim <- max(abs(le_data$median_coef), na.rm = TRUE)

      p_le <- ggplot(le_data, aes(x = gene, y = module_clean)) +
        geom_tile(aes(fill = median_coef), color = "white", linewidth = 0.3) +
        geom_text(aes(label = sig_label), size = 4, vjust = 0.75) +
        scale_fill_gradientn(
          colors = DIVERGING_PALETTE(100),
          limits = c(-coef_lim, coef_lim),
          name = "Gradient\nscore",
          na.value = "grey90"
        ) +
        labs(
          title = paste0(ct, ": Leading Edge Genes of Significant Modules"),
          subtitle = paste0("Tile color = log-rate gradient score (red = induced near ",
                            QUERY_LABEL, ", blue = repressed)\n",
                            "* = gene individually significant (FDR < 0.05)"),
          x = NULL, y = NULL
        ) +
        theme_bw(base_size = 10) +
        theme(
          axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
          axis.text.y = element_text(size = 9),
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 8, color = "grey40"),
          panel.grid = element_blank(),
          legend.position = "right"
        )

      n_genes <- uniqueN(le_data$gene)
      n_mods <- uniqueN(le_data$module_clean)
      ggsave(file.path(OUTPUT_DIR, paste0("leading_edge_", ct, ".pdf")), p_le,
             width = max(8, n_genes * 0.35 + 4),
             height = max(4, n_mods * 0.5 + 2.5))
      message("    ", ct, ": ", n_mods, " modules, ", n_genes, " leading edge genes")
    }

    # --- Combined leading edge across all 5 cell types ---
    # Tile plot: rows = "CellType: Module", columns = leading edge genes
    message("  Creating combined leading edge figure...")

    combined_le <- rbindlist(lapply(LEADING_EDGE_CELLTYPES, function(ct) {
      ct_sig <- custom_fgsea_all[cell_type == ct & padj < 0.05]
      if (nrow(ct_sig) == 0) return(NULL)

      rbindlist(lapply(seq_len(nrow(ct_sig)), function(i) {
        le_genes <- if (is.list(ct_sig$leadingEdge)) {
          unlist(ct_sig$leadingEdge[[i]])
        } else {
          trimws(unlist(strsplit(as.character(ct_sig$leadingEdge[i]), ",")))
        }
        ct_results <- all_results[cell_type == ct, .(gene, median_coef = get(COEF_COL),
                                                        sig_fdr = get(SIG_COL))]
        le_dt <- data.table(
          cell_type = ct,
          module = ct_sig$pathway[i],
          module_clean = ct_sig$module_clean[i],
          module_nes = ct_sig$NES[i],
          module_padj = ct_sig$padj[i],
          gene = le_genes
        )
        merge(le_dt, ct_results, by = "gene", all.x = TRUE)
      }))
    }))

    if (nrow(combined_le) > 0) {
      # Create row labels: "CellType: Module"
      combined_le[, row_label := paste0(cell_type, ": ", module_clean)]

      # Order rows by cell type then NES
      row_order <- combined_le[, .(nes = module_nes[1]),
                                by = .(cell_type, row_label)][order(cell_type, nes)]$row_label
      combined_le[, row_label := factor(row_label, levels = row_order)]

      # Identify genes that appear in multiple cell types (shared leading edge)
      gene_ct_count <- combined_le[, uniqueN(cell_type), by = gene]
      shared_genes <- gene_ct_count[V1 >= 2]$gene

      # Order genes: shared first (alphabetical), then unique (alphabetical)
      gene_order <- c(sort(shared_genes),
                      sort(setdiff(unique(combined_le$gene), shared_genes)))
      combined_le[, gene := factor(gene, levels = gene_order)]

      combined_le[, sig_label := fifelse(!is.na(sig_fdr) & sig_fdr < 0.05, "*", "")]
      combined_le[, is_shared := gene %in% shared_genes]

      coef_lim <- max(abs(combined_le$median_coef), na.rm = TRUE)

      p_combined_le <- ggplot(combined_le, aes(x = gene, y = row_label)) +
        geom_tile(aes(fill = median_coef), color = "white", linewidth = 0.3) +
        geom_text(aes(label = sig_label), size = 3, vjust = 0.75) +
        scale_fill_gradientn(
          colors = DIVERGING_PALETTE(100),
          limits = c(-coef_lim, coef_lim),
          name = "Gradient\nscore",
          na.value = "grey90"
        ) +
        labs(
          title = "Leading Edge Genes Across Key Cell Types",
          subtitle = paste0("Significant custom modules (padj < 0.05) in top cell types\n",
                            "Red = induced near ", QUERY_LABEL,
                            " | * = gene individually FDR < 0.05"),
          x = NULL, y = NULL
        ) +
        theme_bw(base_size = 9) +
        theme(
          axis.text.x = element_text(angle = 60, hjust = 1, size = 7,
                                     face = ifelse(gene_order %in% shared_genes, "bold", "plain")),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 8, color = "grey40"),
          panel.grid = element_blank(),
          legend.position = "right"
        )

      n_rows <- uniqueN(combined_le$row_label)
      n_cols <- uniqueN(combined_le$gene)
      ggsave(file.path(OUTPUT_DIR, "leading_edge_combined.pdf"), p_combined_le,
             width = max(12, n_cols * 0.3 + 5),
             height = max(6, n_rows * 0.35 + 3))
      message("  Saved: leading_edge_combined.pdf (", n_rows, " rows x ", n_cols, " genes)")

      # Save table
      fwrite(combined_le[, .(cell_type, module_clean, gene, median_coef, sig_fdr,
                             module_nes, module_padj, is_shared)][order(cell_type, module_clean, gene)],
             file.path(OUTPUT_DIR, "leading_edge_genes_table.csv"))
      message("  Saved: leading_edge_genes_table.csv")
    }

    message("  Saved: leading_edge_*.pdf (per cell type)")

    # --- Summary ---
    message("\n  Custom module enrichment summary:")
    for (ct in levels(droplevels(all_results$cell_type))) {
      if (ct %in% custom_fgsea_all$cell_type) {
        n_sig <- sum(custom_fgsea_all[cell_type == ct]$padj < 0.05, na.rm = TRUE)
        top_mod <- custom_fgsea_all[cell_type == ct & padj < 0.05][order(pval)]$module_clean[1]
        if (!is.na(top_mod)) {
          message(sprintf("    %-18s  %2d sig modules (top: %s)", ct, n_sig, top_mod))
        } else {
          message(sprintf("    %-18s  %2d sig modules", ct, n_sig))
        }
      }
    }
  }
} else {
  message("\n[SKIP] Panels 8b-9b: fgsea not available")
}

# =============================================================================
# Summary Panel: Top Genes Dot Plot (key cell types)
# =============================================================================
# Shows the top non-contamination, non-query-signature genes for
# the cell types with the most specific significant genes.

message("\nSummary panel: Top genes dot plot for key cell types...")

# Use top cell types by number of specific significant genes, or all available
SUMMARY_CELLTYPES_CANDIDATES <- sig_results[
  specificity %in% c("specific", "moderate"),
  .(n_sig = .N), by = cell_type
][order(-n_sig)]$cell_type
SUMMARY_CELLTYPES <- head(SUMMARY_CELLTYPES_CANDIDATES, 5)
if (length(SUMMARY_CELLTYPES) == 0) SUMMARY_CELLTYPES <- head(as.character(ct_order), 5)
TOP_N_SUMMARY <- 10  # genes per cell type

# Exclude contamination and query signature genes
top_genes_summary <- sig_results[
  cell_type %in% SUMMARY_CELLTYPES &
  !gene %in% contamination_genes &
  !gene %in% HYMY_SIGNATURE_GENES
][order(-abs(get(COEF_COL))), head(.SD, TOP_N_SUMMARY), by = cell_type]

p_topgenes <- NULL
if (nrow(top_genes_summary) > 0) {
  # Get full data for these genes across the 5 cell types
  topg_genes <- unique(top_genes_summary$gene)
  topg_data <- all_results[gene %in% topg_genes & cell_type %in% SUMMARY_CELLTYPES]
  topg_data[, neg_log10_fdr := -log10(pmax(get(SIG_COL), 1e-50))]
  topg_data[, neg_log10_fdr_capped := pmin(neg_log10_fdr, 20)]

  # Order genes: group by their primary cell type, then by coefficient
  gene_primary <- top_genes_summary[, .(primary_ct = cell_type[which.max(abs(get(COEF_COL)))]),
                                     by = gene]
  gene_primary[, primary_ct := factor(primary_ct, levels = SUMMARY_CELLTYPES)]
  gene_order_tg <- merge(
    topg_data[, .(max_coef = max(abs(get(COEF_COL)), na.rm = TRUE)), by = gene],
    gene_primary, by = "gene"
  )[order(primary_ct, -max_coef)]$gene
  topg_data[, gene := factor(gene, levels = rev(gene_order_tg))]
  topg_data[, cell_type := factor(cell_type, levels = SUMMARY_CELLTYPES)]

  coef_lim_tg <- max(abs(topg_data[[COEF_COL]]), na.rm = TRUE)

  p_topgenes <- ggplot(topg_data, aes(x = cell_type, y = gene)) +
    geom_point(aes(size = neg_log10_fdr_capped, color = .data[[COEF_COL]])) +
    scale_color_gradientn(
      colors = DIVERGING_PALETTE(100),
      limits = c(-coef_lim_tg, coef_lim_tg),
      name = "Log-rate\ncoefficient"
    ) +
    scale_size_continuous(
      range = c(0.5, 5),
      name = expression(-log[10](FDR)),
      breaks = c(2, 5, 10, 20)
    ) +
    labs(
      title = paste0("Top Spatial Gradient Genes (non-", QUERY_LABEL, " signature)"),
      subtitle = paste0("Top ", TOP_N_SUMMARY, " genes per cell type by effect size | ",
                        "Red/blue = induced/repressed near ", QUERY_LABEL),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 9, face = "italic"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      panel.grid.major = element_line(color = "grey92"),
      legend.key.size = unit(0.4, "cm")
    )

  ggsave(file.path(OUTPUT_DIR, "top_genes_key_celltypes.pdf"), p_topgenes,
         width = 6, height = max(5, length(topg_genes) * 0.22 + 2))
  message("  Saved: top_genes_key_celltypes.pdf (",
          length(topg_genes), " genes x ", length(SUMMARY_CELLTYPES), " cell types)")

  # Panel G standalone: narrower + taller for publication
  p_topgenes_G <- p_topgenes +
    labs(tag = "G") +
    theme(
      axis.text.y = element_text(size = 8, face = "italic"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      legend.position = "right",
      legend.key.size = unit(0.35, "cm"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8)
    )
  ggsave(file.path(OUTPUT_DIR, "panel_G_top_genes.pdf"), p_topgenes_G,
         width = 4.5, height = max(8, length(topg_genes) * 0.3 + 2))
  message("  Saved: panel_G_top_genes.pdf (standalone)")
} else {
  message("  [SKIP] No qualifying genes for top genes panel")
}

# =============================================================================
# CD8 T Cell Exhaustion Dot Plot (Standalone)
# =============================================================================
# Focused dot plot of exhaustion/effector genes in CD8_T_cells.
# Combines exhaustion markers, effector molecules, and Tpex/Tex signatures
# to show how query cell proximity shapes the CD8 T cell state.

message("\nCD8 exhaustion dot plot...")

EXHAUSTION_GENES <- unique(c(
  # Exhaustion / inhibitory receptors
  "Pdcd1", "Lag3", "Havcr2", "Tigit", "Ctla4", "Tox", "Entpd1", "Cd38", "Cd244a", "Cd48",
  # Effector molecules
  "Gzmb", "Gzma", "Gzmk", "Prf1", "Ifng", "Tnf", "Nkg7", "Fasl",
  # Transcription factors
  "Eomes", "Tbx21", "Tcf7", "Id3", "Bach2",
  # Tpex markers
  "Slamf6", "Xcl1", "Il7r", "Cxcr5",
  # Activation / costimulatory
  "Cd69", "Icos", "Tnfrsf9", "Tnfrsf4"
))

cd8_exhaust_data <- all_results[cell_type == "CD8_T_cells" & gene %in% EXHAUSTION_GENES]

if (nrow(cd8_exhaust_data) > 0) {
  cd8_exhaust_data[, neg_log10_fdr := -log10(pmax(get(SIG_COL), 1e-50))]
  cd8_exhaust_data[, neg_log10_fdr_capped := pmin(neg_log10_fdr, 20)]
  cd8_exhaust_data[, is_sig := get(SIG_COL) < 0.05]

  # Annotate functional category for grouping on y-axis
  cd8_exhaust_data[, category := fcase(
    gene %in% c("Pdcd1","Lag3","Havcr2","Tigit","Ctla4","Entpd1","Cd38","Cd244a","Cd48"), "Inhibitory receptors",
    gene %in% c("Tox","Eomes","Tbx21","Tcf7","Id3","Bach2"), "Transcription factors",
    gene %in% c("Gzmb","Gzma","Gzmk","Prf1","Ifng","Tnf","Nkg7","Fasl"), "Effector molecules",
    gene %in% c("Slamf6","Xcl1","Il7r","Cxcr5"), "Tpex markers",
    default = "Activation"
  )]
  cat_order <- c("Inhibitory receptors", "Transcription factors",
                 "Effector molecules", "Tpex markers", "Activation")
  cd8_exhaust_data[, category := factor(category, levels = cat_order)]

  # Order genes: by category then by coefficient
  cd8_exhaust_data <- cd8_exhaust_data[order(category, get(COEF_COL))]
  cd8_exhaust_data[, gene := factor(gene, levels = unique(gene))]

  coef_lim_ex <- max(abs(cd8_exhaust_data[[COEF_COL]]), na.rm = TRUE)

  p_cd8_exhaust <- ggplot(cd8_exhaust_data, aes(x = .data[[COEF_COL]], y = gene)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(size = neg_log10_fdr_capped,
                   fill = .data[[COEF_COL]],
                   color = is_sig),
               shape = 21, stroke = 0.6) +
    scale_fill_gradientn(
      colors = DIVERGING_PALETTE(100),
      limits = c(-coef_lim_ex, coef_lim_ex),
      name = "Log-rate\ncoefficient"
    ) +
    scale_color_manual(
      values = c("TRUE" = "black", "FALSE" = "grey70"),
      guide = "none"
    ) +
    scale_size_continuous(
      range = c(1.5, 6),
      name = expression(-log[10](FDR)),
      breaks = c(2, 5, 10, 20)
    ) +
    facet_grid(category ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = paste0("CD8 T Cell Exhaustion & Effector Genes — ",
                     "Spatial Gradient near ", QUERY_LABEL),
      subtitle = paste0("Negative coefficient = higher expression rate near ",
                        QUERY_LABEL, " | Black border = FDR < 0.05"),
      x = "Log-rate coefficient (per µm)", y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 9, face = "italic"),
      strip.text.y = element_text(size = 8, face = "bold", angle = 0),
      strip.background = element_rect(fill = "grey95"),
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8, color = "grey40"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.4, "cm")
    )

  n_genes_ex <- nrow(cd8_exhaust_data)
  ggsave(file.path(OUTPUT_DIR, "cd8_exhaustion_dotplot.pdf"), p_cd8_exhaust,
         width = 5.5, height = max(5, n_genes_ex * 0.28 + 2.5))
  message("  Saved: cd8_exhaustion_dotplot.pdf (", n_genes_ex, " genes)")
} else {
  message("  [SKIP] No exhaustion genes found for CD8_T_cells")
}

# =============================================================================
# Combined Summary Figure (Multi-Panel, Publication-Ready)
# =============================================================================
# Layout:
#   Row 1 (Discovery):  A = gene counts bar chart, B = Stage 2 classification
#   Row 2 (Validation): C = permutation validation,  D = specificity breakdown
#   Row 3 (Hallmark):   E = Hallmark fgsea dot plot (full width)
#   Row 4 (Modules):    F = custom module dot plot (full width)
#   Row 5 (Top Genes):  G = top genes dot plot for key cell types (full width)
#
# Design rationale: Tells the story top-to-bottom:
#   "What did we find?" -> "How confident are we?" -> "What pathways?" ->
#   "What modules?" -> "Which genes?"

message("\nAssembling combined summary figure...")

# --- Row 1: Discovery ---
row1_panels <- list()
if (exists("p1") && !is.null(p1)) {
  row1_panels$A <- p1 + labs(tag = "A")
}
if (exists("p10") && !is.null(p10)) {
  row1_panels$B <- p10 + labs(tag = "B")
}

# --- Row 2: Validation ---
row2_panels <- list()
if (exists("p4") && !is.null(p4)) {
  row2_panels$C <- p4 + labs(tag = "C")
}
if (exists("p5") && !is.null(p5)) {
  row2_panels$D <- p5 + labs(tag = "D")
}

# --- Row 3: Hallmark fgsea ---
row3_panel <- NULL
if (exists("p8") && !is.null(p8)) {
  row3_panel <- p8 + labs(tag = "E")
}

# --- Row 4: Custom modules ---
row4_panel <- NULL
if (exists("p_custom_dot") && !is.null(p_custom_dot)) {
  row4_panel <- p_custom_dot + labs(tag = "F")
}

# --- Row 5: Top genes ---
row5_panel <- NULL
if (exists("p_topgenes") && !is.null(p_topgenes)) {
  row5_panel <- p_topgenes + labs(tag = "G")
}

# Build paired rows
build_row <- function(panels, ncol = 2) {
  if (length(panels) == 2) wrap_plots(panels, ncol = ncol)
  else if (length(panels) == 1) panels[[1]]
  else NULL
}

row1 <- build_row(row1_panels)
row2 <- build_row(row2_panels)

# Assemble all rows
all_rows <- Filter(Negate(is.null), list(row1, row2, row3_panel, row4_panel, row5_panel))

if (length(all_rows) > 0) {
  # Heights: paired rows = 1, full-width pathway/module rows = 1, gene row = 1.2
  row_heights <- sapply(all_rows, function(r) {
    if (identical(r, row1) || identical(r, row2)) 1
    else if (identical(r, row5_panel)) 1.2
    else 1
  })

  combined <- wrap_plots(all_rows, ncol = 1, heights = row_heights) +
    plot_annotation(
      title = paste0(QUERY_LABEL, " Distance Correlation Analysis (",
                     ANNOTATION_LEVEL, " annotation)"),
      subtitle = paste0(nrow(sig_results), " significant gene-celltype pairs across ",
                        n_ct, " cell types | ",
                        sum(gene_ct_counts$specificity %in% c("specific", "moderate")),
                        " specific genes"),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 11, color = "grey30")
      )
    )

  n_rows <- length(all_rows)
  ggsave(file.path(OUTPUT_DIR, "summary_figure.pdf"), combined,
         width = 10, height = 5 * n_rows)
  message("Saved: summary_figure.pdf (", n_rows, " rows)")
} else {
  message("  [SKIP] No panels available for summary figure")
}

# =============================================================================
# Output Tables
# =============================================================================

# Save specific gene list (the biologically interesting ones)
# Include Stage 2 classification if available
cols_to_select <- c("gene", "cell_type", COEF_COL, SIG_COL, "direction",
                     "gradient_score", "decay_pattern", "perm_pval")
if (HAS_STAGE2) {
  s2_out_cols <- intersect(c("classification", "stage2_coef", "stage2_fdr"), names(all_results))
  cols_to_select <- c(cols_to_select, s2_out_cols)
}
# B4: Include v2 per-sample diagnostic columns if present
diag_cols <- c("sign_consistency", "median_dispersion",
               "n_positive_samples", "n_negative_samples")
cols_to_select <- c(cols_to_select, diag_cols)
cols_to_select <- intersect(cols_to_select, names(sig_results))

specific_gene_table <- merge(
  sig_results[specificity %in% c("specific", "moderate"), ..cols_to_select],
  gene_ct_counts[, .(gene, n_celltypes, specificity)],
  by = "gene"
)[order(specificity, cell_type, get(SIG_COL))]

fwrite(specific_gene_table, file.path(OUTPUT_DIR, "specific_significant_genes.csv"))
message("Saved: specific_significant_genes.csv (",
        nrow(specific_gene_table), " gene-celltype pairs)")

# Stage 2 summary in statistics
if (HAS_STAGE2) {
  message("\nStage 2 Classification Summary:")
  class_counts <- sig_results[specificity %in% c("specific", "moderate"),
                               .N, by = classification][order(-N)]
  for (i in seq_len(nrow(class_counts))) {
    message("  ", class_counts$classification[i], ": ", class_counts$N[i])
  }
}

# =============================================================================
# Summary Statistics
# =============================================================================

message("\n", strrep("=", 70))
message("Summary Statistics")
message(strrep("=", 70))

# Per cell type summary
ct_summary <- all_results[, .(
  total_genes = .N,
  sig_genes = sum(get(SIG_COL) < FDR_THRESHOLD),
  sig_specific = sum(get(SIG_COL) < FDR_THRESHOLD & specificity %in% c("specific", "moderate"),
                     na.rm = TRUE),
  sig_contamination = sum(get(SIG_COL) < FDR_THRESHOLD & specificity %in% c("ubiquitous", "contamination"),
                          na.rm = TRUE),
  induced = sum(get(SIG_COL) < FDR_THRESHOLD & get(COEF_COL) < 0),
  repressed = sum(get(SIG_COL) < FDR_THRESHOLD & get(COEF_COL) > 0),
  perm_tested = sum(!is.na(perm_pval)),
  perm_sig = sum(perm_pval < 0.05, na.rm = TRUE)
), by = cell_type][order(-sig_specific)]

fwrite(ct_summary, file.path(OUTPUT_DIR, "celltype_summary.csv"))

for (i in seq_len(nrow(ct_summary))) {
  r <- ct_summary[i]
  message(sprintf("  %-18s  sig: %4d (specific: %3d, contamination: %3d)  perm: %d/%d",
                  r$cell_type, r$sig_genes, r$sig_specific, r$sig_contamination,
                  r$perm_sig, r$perm_tested))
}

message("\nTotal unique significant genes: ", uniqueN(sig_results$gene))
message("  Specific (1-3 cell types): ", sum(gene_ct_counts$specificity %in% c("specific", "moderate")))
message("  Potential contamination (>=", CONTAMINATION_THRESHOLD, "): ",
        sum(gene_ct_counts$specificity == "contamination"))

message("\n", strrep("=", 70))
message("Atlas complete! Output: ", OUTPUT_DIR)
message(strrep("=", 70))
