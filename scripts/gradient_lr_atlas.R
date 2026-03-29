#!/usr/bin/env Rscript
#' =============================================================================
#' RIPPLE Stage 6: Gradient-to-LR Integration Atlas & Summary
#' =============================================================================
#'
#' Merges per-celltype results from gradient_lr_integration.R into summary
#' tables and publication figures.
#'
#' Run AFTER all array jobs from gradient_lr_integration.R are complete.
#'
#' Usage:
#'   Rscript gradient_lr_atlas.R
#'   QUERY_CELLTYPE=MyType CELLTYPE_COLUMN=my_col Rscript gradient_lr_atlas.R
#'
#' Author: CMM Project
#' =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(pheatmap)
  library(dplyr)
  library(tidyr)
})

set.seed(42)

# Source utilities
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  return(getwd())
}
script_dir <- get_script_dir()
source(file.path(script_dir, "utils.R"))

# =============================================================================
# Configuration
# =============================================================================

ANALYSIS_NAME <- "gradient_lr_integration"

# Gradient source routing (matches gradient_lr_integration.R)
GRADIENT_SOURCE <- Sys.getenv("GRADIENT_SOURCE", unset = "hymy_distance_correlation")
gradient_suffix <- sub("^hymy_distance_correlation", "", GRADIENT_SOURCE)
output_name <- if (nchar(gradient_suffix) > 0) {
  paste0(ANALYSIS_NAME, gradient_suffix)
} else {
  ANALYSIS_NAME
}

# Inherited from config.R (via utils.R): QUERY_CELLTYPE, CELLTYPE_COL, OUTPUT_SUFFIX, QUERY_LABEL
OUTPUT_BASE <- file.path(OUTPUT_ROOT, output_name)

# NicheNet cross-reference dirs (legacy paths, only meaningful for HyMy/L1)
if (ANNOTATION_LEVEL == "L1") {
  NICHENET_LEC_DIR <- file.path(PROJECT_ROOT, "results", "spatial_nichenet_L1")
  NICHENET_FRC_DIR <- file.path(PROJECT_ROOT, "results", "spatial_nichenet_L1_FRC")
} else {
  NICHENET_LEC_DIR <- file.path(PROJECT_ROOT, "results", "spatial_nichenet")
  NICHENET_FRC_DIR <- file.path(PROJECT_ROOT, "results", "spatial_nichenet_FRC")
}

PERCELLTYPE_DIR <- file.path(OUTPUT_BASE, "per_celltype")
SUMMARY_DIR <- file.path(OUTPUT_BASE, "summary")
PLOT_DIR <- file.path(OUTPUT_BASE, "plots")

ensure_dir(SUMMARY_DIR)
ensure_dir(PLOT_DIR)

message(strrep("=", 70))
message("Gradient-to-LR Integration Atlas")
message(strrep("=", 70))
message("Annotation level: ", ANNOTATION_LEVEL)
message("Input: ", PERCELLTYPE_DIR)

# =============================================================================
# Load Per-Celltype Results
# =============================================================================

message("\n", strrep("-", 70))
message("Loading per-celltype results...")
message(strrep("-", 70), "\n")

ct_dirs <- list.dirs(PERCELLTYPE_DIR, recursive = FALSE, full.names = TRUE)

if (length(ct_dirs) == 0) {
  stop("No per-celltype results found in: ", PERCELLTYPE_DIR,
       "\nRun gradient_lr_integration.R first.")
}

# Helper: read CSV if it exists and has rows
safe_fread <- function(path) {
  if (file.exists(path)) {
    dt <- fread(path)
    if (nrow(dt) > 0) return(dt)
  }
  return(NULL)
}

# Merge all per-celltype files
all_combined <- rbindlist(lapply(ct_dirs, function(d) {
  safe_fread(file.path(d, "combined_prioritization.csv"))
}), fill = TRUE)

all_direct_a <- rbindlist(lapply(ct_dirs, function(d) {
  safe_fread(file.path(d, "direct_lr_pairs_hymy_to_target.csv"))
}), fill = TRUE)

all_direct_b <- rbindlist(lapply(ct_dirs, function(d) {
  safe_fread(file.path(d, "direct_lr_pairs_target_to_hymy.csv"))
}), fill = TRUE)

all_activity <- rbindlist(lapply(ct_dirs, function(d) {
  safe_fread(file.path(d, "nichenet_ligand_activity.csv"))
}), fill = TRUE)

all_enrichment <- rbindlist(lapply(ct_dirs, function(d) {
  safe_fread(file.path(d, "downstream_target_enrichment.csv"))
}), fill = TRUE)

message(sprintf("  Combined prioritization: %d rows across %d cell types",
                nrow(all_combined),
                length(unique(all_combined$cell_type))))
message(sprintf("  Direction A (%s->Target): %d L-R pairs", QUERY_LABEL, nrow(all_direct_a)))
message(sprintf("  Direction B (Target->%s): %d L-R pairs", QUERY_LABEL, nrow(all_direct_b)))
message(sprintf("  Ligand activity scores: %d", nrow(all_activity)))
message(sprintf("  Enrichment tests: %d", nrow(all_enrichment)))

# =============================================================================
# Artifact Classification
# =============================================================================

message("\n", strrep("-", 70))
message("Classifying potential segmentation artifacts...")
message(strrep("-", 70), "\n")

if (nrow(all_combined) > 0) {
  # Query signature genes — receptors matching these on non-myeloid cells are artifacts
  HYMY_SIGNATURE <- make.names(QUERY_SIGNATURE)

  # Genes significant in >=4 cell types (likely leakage from distance correlation)
  CONTAMINATION_CANDIDATES <- make.names(c("Bub1b", "C1qa", "C1qb", "Cd52", "Cebpb",
                                            "Col6a2", "Csf3r", "Epcam", "Grn", "Hist2h2bb",
                                            "Igkc", "Il1b", "Irf4", "Jchain", "Map1lc3a",
                                            "Mzb1", "Pdcd4", "Pon3", "Rnaset2a", "Sdc1",
                                            "Slpi", "Tent5c"))

  # Lineage-mismatched markers
  MYELOID_MARKERS <- make.names(c("Sirpa", "Csf1r"))
  EPITHELIAL_MARKERS <- make.names(c("Epcam"))
  PDC_MARKERS <- make.names(c("Siglech"))

  MYELOID_CELLTYPES <- c("Monocyte", "Macrophages", "cDC1", "cDC2", "mature_migDC")

  # Initialize all as clean
  all_combined[, artifact_flag := "clean"]

  # Rule 1: Receptor is query signature gene on non-myeloid cell type
  all_combined[receptor %in% HYMY_SIGNATURE &
               !(cell_type %in% MYELOID_CELLTYPES),
               artifact_flag := "artifact"]

  # Rule 2: Receptor in contamination list on non-myeloid at low expression (<5%)
  # NOTE: receptor_pct_target is stored as fraction (0-1), not percentage
  all_combined[artifact_flag == "clean" &
               receptor %in% CONTAMINATION_CANDIDATES &
               !(cell_type %in% MYELOID_CELLTYPES) &
               receptor_pct_target < 0.05,
               artifact_flag := "suspect"]

  # Rule 3: Lineage-mismatched markers at low expression
  all_combined[artifact_flag == "clean" &
               receptor %in% MYELOID_MARKERS &
               !(cell_type %in% MYELOID_CELLTYPES) &
               receptor_pct_target < 0.05,
               artifact_flag := "suspect"]

  all_combined[artifact_flag == "clean" &
               receptor %in% EPITHELIAL_MARKERS,
               artifact_flag := "suspect"]

  all_combined[artifact_flag == "clean" &
               receptor %in% PDC_MARKERS,
               artifact_flag := "suspect"]

  # Rule 4: Very low receptor expression (<2%) on non-myeloid cells
  all_combined[artifact_flag == "clean" &
               !(cell_type %in% MYELOID_CELLTYPES) &
               receptor_pct_target < 0.02,
               artifact_flag := "low_confidence"]

  # Summary
  artifact_summary <- all_combined[, .N, by = artifact_flag]
  message("  Artifact classification:")
  for (i in seq_len(nrow(artifact_summary))) {
    message(sprintf("    %s: %d pairs", artifact_summary$artifact_flag[i],
                    artifact_summary$N[i]))
  }

  n_clean <- sum(all_combined$artifact_flag == "clean")
  message(sprintf("\n  Clean pairs for figures: %d / %d (%.0f%%)",
                  n_clean, nrow(all_combined),
                  100 * n_clean / nrow(all_combined)))
}

# =============================================================================
# Summary Tables
# =============================================================================

message("\n", strrep("-", 70))
message("Generating summary tables...")
message(strrep("-", 70), "\n")

# 1. All L-R pairs combined (Direction A)
fwrite(all_combined, file.path(SUMMARY_DIR, "all_lr_pairs_combined.csv"))

# 2. Top 20 per cell type
if (nrow(all_combined) > 0) {
  top_per_ct <- all_combined[artifact_flag == "clean", head(.SD, 20), by = cell_type]
  fwrite(top_per_ct, file.path(SUMMARY_DIR, "top_lr_pairs_per_celltype.csv"))
  message(sprintf("  Top 20 per cell type: %d rows", nrow(top_per_ct)))
}

# 3. Direction B pairs
fwrite(all_direct_b, file.path(SUMMARY_DIR, "all_target_to_hymy_pairs.csv"))

# 4. Coverage stats: how many gradient genes mapped to L-R pairs
coverage_stats <- rbindlist(lapply(ct_dirs, function(d) {
  ct <- basename(d)
  combined <- safe_fread(file.path(d, "combined_prioritization.csv"))
  activity <- safe_fread(file.path(d, "nichenet_ligand_activity.csv"))

  n_combined <- if (!is.null(combined)) nrow(combined) else 0
  n_unique_receptors <- if (!is.null(combined)) length(unique(combined$receptor)) else 0
  n_unique_ligands <- if (!is.null(combined)) length(unique(combined$ligand)) else 0
  n_activity_ligands <- if (!is.null(activity)) nrow(activity) else 0
  top_activity <- if (!is.null(activity) && nrow(activity) > 0) {
    activity$activity[1]
  } else {
    NA_real_
  }

  data.table(
    cell_type = ct,
    n_lr_pairs = n_combined,
    n_unique_ligands = n_unique_ligands,
    n_unique_receptors = n_unique_receptors,
    n_activity_ligands = n_activity_ligands,
    top_activity_score = top_activity
  )
}))

fwrite(coverage_stats, file.path(SUMMARY_DIR, "coverage_stats.csv"))
message(sprintf("  Coverage stats: %d cell types", nrow(coverage_stats)))

# =============================================================================
# Cross-Reference with Existing NicheNet Results
# =============================================================================

message("\n", strrep("-", 70))
message("Cross-referencing with existing NicheNet results...")
message(strrep("-", 70), "\n")

crossref_nichenet <- function(new_results, nichenet_dir, target_type) {
  nichenet_file <- file.path(nichenet_dir, "all_prioritized_lr_pairs.csv")
  if (!file.exists(nichenet_file)) {
    message(sprintf("  NicheNet %s results not found. Skipping.", target_type))
    return(NULL)
  }

  nichenet <- fread(nichenet_file)
  message(sprintf("  NicheNet %s: %d L-R pairs loaded", target_type, nrow(nichenet)))

  # Filter new results to matching cell type
  new_ct <- new_results[cell_type == target_type]
  if (nrow(new_ct) == 0) {
    message(sprintf("  No gradient L-R pairs for %s. Skipping.", target_type))
    return(NULL)
  }

  # Merge on ligand-receptor pair
  crossref <- merge(
    new_ct[, .(ligand, receptor, combined_score, direct_score,
               nichenet_activity, enrichment_pval)],
    nichenet[, .(ligand, receptor,
                 nichenet_prioritization = prioritization_score)],
    by = c("ligand", "receptor"),
    all = TRUE
  )

  crossref[, source := fifelse(
    !is.na(combined_score) & !is.na(nichenet_prioritization), "both",
    fifelse(!is.na(combined_score), "gradient_only", "nichenet_only")
  )]

  # Agreement stats
  n_both <- sum(crossref$source == "both")
  n_gradient_only <- sum(crossref$source == "gradient_only")
  n_nichenet_only <- sum(crossref$source == "nichenet_only")

  message(sprintf("  %s: %d shared, %d gradient-only, %d nichenet-only",
                  target_type, n_both, n_gradient_only, n_nichenet_only))

  # Rank correlation for shared pairs
  if (n_both >= 5) {
    shared <- crossref[source == "both" &
                         !is.na(combined_score) & !is.na(nichenet_prioritization)]
    if (nrow(shared) >= 5) {
      cor_test <- cor.test(shared$combined_score, shared$nichenet_prioritization,
                           method = "spearman")
      message(sprintf("  Spearman rank correlation: %.3f (p=%.3e)",
                      cor_test$estimate, cor_test$p.value))
    }
  }

  return(crossref)
}

crossref_lec <- crossref_nichenet(all_combined, NICHENET_LEC_DIR, "LEC")
if (!is.null(crossref_lec)) {
  fwrite(crossref_lec, file.path(SUMMARY_DIR, "cross_reference_nichenet_lec.csv"))
}

crossref_frc <- crossref_nichenet(all_combined, NICHENET_FRC_DIR, "FRC")
if (!is.null(crossref_frc)) {
  fwrite(crossref_frc, file.path(SUMMARY_DIR, "cross_reference_nichenet_frc.csv"))
}

# =============================================================================
# Visualization: Bubble Plot (Top L-R Pairs x Cell Types)
# =============================================================================

message("\n", strrep("-", 70))
message("Generating figures...")
message(strrep("-", 70), "\n")

if (nrow(all_combined) > 0) {

  # Get top 5 CLEAN L-R pairs per cell type (by combined score)
  top_pairs <- all_combined[artifact_flag == "clean", head(.SD, 5), by = cell_type]
  top_pair_labels <- unique(paste0(top_pairs$ligand, "-", top_pairs$receptor))

  # Filter to these pairs across all cell types
  all_combined[, pair_label := paste0(ligand, "-", receptor)]
  bubble_data <- all_combined[pair_label %in% top_pair_labels & artifact_flag == "clean"]

  if (nrow(bubble_data) > 0) {
    # Order cell types by number of significant pairs
    ct_order <- coverage_stats[order(-n_lr_pairs)]$cell_type
    bubble_data[, cell_type := factor(cell_type, levels = ct_order)]

    # Order pairs by max combined score
    pair_order <- bubble_data[, .(max_score = max(combined_score, na.rm = TRUE)),
                              by = pair_label]
    setorder(pair_order, -max_score)
    bubble_data[, pair_label := factor(pair_label, levels = pair_order$pair_label)]

    p_bubble <- ggplot(bubble_data, aes(x = cell_type, y = pair_label)) +
      geom_point(aes(size = combined_score,
                     color = receptor_gradient_coef)) +
      scale_size_continuous(range = c(1, 8), name = "Combined\nScore") +
      scale_color_gradient2(low = "blue", mid = "white", high = "red",
                            midpoint = 0,
                            name = "Receptor\nGradient\nCoef",
                            limits = c(-0.01, 0.01),
                            oob = squish) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major = element_line(color = "grey90")) +
      labs(title = paste0("Top L-R Pairs Explaining ", QUERY_LABEL, " Gradients (Clean Only)"),
           subtitle = paste0("Direction A: ", QUERY_LABEL, " ligand -> Target receptor (artifacts excluded)"),
           x = "Target Cell Type", y = "Ligand-Receptor Pair")

    ggsave(file.path(PLOT_DIR, "bubble_plot_top_pairs.pdf"),
           p_bubble, width = 12, height = 10)
    message("  Saved: bubble_plot_top_pairs.pdf")
  }

  # ===========================================================================
  # Heatmap: Combined Score (L-R Pairs x Cell Types)
  # ===========================================================================

  # Get top 30 unique CLEAN pairs across all cell types
  top30_pairs <- all_combined[artifact_flag == "clean",
                              .(max_score = max(combined_score, na.rm = TRUE)),
                              by = pair_label]
  setorder(top30_pairs, -max_score)
  top30_labels <- head(top30_pairs$pair_label, 30)

  heatmap_data <- all_combined[pair_label %in% top30_labels & artifact_flag == "clean",
                               .(pair_label, cell_type, combined_score)]

  if (nrow(heatmap_data) > 0) {
    heatmap_wide <- dcast(heatmap_data, pair_label ~ cell_type,
                          value.var = "combined_score", fill = 0)
    heatmap_mat <- as.matrix(heatmap_wide[, -1, with = FALSE])
    rownames(heatmap_mat) <- heatmap_wide$pair_label

    pdf(file.path(PLOT_DIR, "heatmap_lr_by_celltype.pdf"), width = 12, height = 10)
    pheatmap(heatmap_mat,
             color = colorRampPalette(c("white", "orange", "red"))(100),
             cluster_rows = TRUE, cluster_cols = TRUE,
             main = "Combined Prioritization Score: L-R Pairs x Cell Types (Clean)",
             fontsize_row = 8, fontsize_col = 9)
    dev.off()
    message("  Saved: heatmap_lr_by_celltype.pdf")
  }

  # ===========================================================================
  # Bar Plot: Coverage by Cell Type
  # ===========================================================================

  p_coverage <- ggplot(coverage_stats, aes(x = reorder(cell_type, -n_lr_pairs),
                                           y = n_lr_pairs)) +
    geom_col(fill = "#FC8D62") +
    geom_text(aes(label = n_lr_pairs), vjust = -0.3, size = 3) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0("L-R Pairs per Cell Type (Direction A: ", QUERY_LABEL, " -> Target)"),
         x = "Cell Type", y = "Number of L-R Pairs")

  ggsave(file.path(PLOT_DIR, "coverage_by_celltype.pdf"),
         p_coverage, width = 10, height = 6)
  message("  Saved: coverage_by_celltype.pdf")

  # ===========================================================================
  # Stacked Bar: Artifact Classification by Cell Type
  # ===========================================================================

  artifact_by_ct <- all_combined[, .N, by = .(cell_type, artifact_flag)]
  artifact_by_ct[, artifact_flag := factor(artifact_flag,
    levels = c("clean", "low_confidence", "suspect", "artifact"))]

  artifact_colors <- c("clean" = "#66C2A5", "low_confidence" = "#FFD92F",
                        "suspect" = "#FC8D62", "artifact" = "#E78AC3")

  # Order cell types by total count
  ct_totals <- artifact_by_ct[, .(total = sum(N)), by = cell_type]
  setorder(ct_totals, -total)

  p_artifact <- ggplot(artifact_by_ct,
    aes(x = factor(cell_type, levels = ct_totals$cell_type),
        y = N, fill = artifact_flag)) +
    geom_col() +
    scale_fill_manual(values = artifact_colors, name = "Classification",
                      drop = FALSE) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "L-R Pair Artifact Classification by Cell Type",
         subtitle = paste0("Based on receptor lineage, ", QUERY_LABEL, " signature, and expression thresholds"),
         x = "Cell Type", y = "Number of L-R Pairs")

  ggsave(file.path(PLOT_DIR, "artifact_classification.pdf"),
         p_artifact, width = 10, height = 6)
  message("  Saved: artifact_classification.pdf")

  # ===========================================================================
  # Reproducibility Bar: Per-Mouse Support for Clean Pairs
  # ===========================================================================

  if ("n_samples_supporting" %in% names(all_combined)) {
    clean_reprod <- all_combined[artifact_flag == "clean" &
                                  !is.na(n_samples_supporting)]

    if (nrow(clean_reprod) > 0) {
      reprod_summary <- clean_reprod[, .N, by = .(cell_type, n_samples_supporting)]

      p_reprod <- ggplot(reprod_summary,
        aes(x = factor(cell_type, levels = ct_totals$cell_type),
            y = N, fill = factor(n_samples_supporting))) +
        geom_col(position = "stack") +
        scale_fill_brewer(palette = "YlGnBu", name = "Samples\nSupporting",
                          direction = 1) +
        theme_bw(base_size = 11) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Per-Mouse Reproducibility of Clean L-R Pairs",
             subtitle = "How many TDLN samples independently support each pair",
             x = "Cell Type", y = "Number of L-R Pairs")

      ggsave(file.path(PLOT_DIR, "reproducibility_by_celltype.pdf"),
             p_reprod, width = 10, height = 6)
      message("  Saved: reproducibility_by_celltype.pdf")
    }
  }

  # ===========================================================================
  # Bar Plot: Downstream Enrichment for Top Pairs
  # ===========================================================================

  if (nrow(all_enrichment) > 0) {
    # Top 20 enriched ligands (lowest Fisher p)
    top_enriched <- all_enrichment[fdr_fisher < 0.05]
    setorder(top_enriched, pvalue_fisher)
    top_enriched <- head(top_enriched, 30)

    if (nrow(top_enriched) > 0) {
      top_enriched[, label := paste0(ligand, " (", cell_type, ")")]

      p_enrichment <- ggplot(top_enriched, aes(x = reorder(label, -log10(pvalue_fisher)),
                                               y = -log10(pvalue_fisher))) +
        geom_col(aes(fill = n_overlap)) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
        scale_fill_viridis_c(name = "Overlap\nGenes") +
        coord_flip() +
        theme_bw(base_size = 10) +
        labs(title = "Downstream Target Enrichment",
             subtitle = "Fisher's exact test: predicted targets vs gradient genes",
             x = "", y = "-log10(p-value)")

      ggsave(file.path(PLOT_DIR, "enrichment_validation.pdf"),
             p_enrichment, width = 10, height = 8)
      message("  Saved: enrichment_validation.pdf")
    }
  }

  # ===========================================================================
  # Scatter: NicheNet Activity vs Direct Score (where both exist)
  # ===========================================================================

  scatter_data <- all_combined[artifact_flag == "clean" &
                               !is.na(nichenet_activity) & !is.na(direct_score)]

  if (nrow(scatter_data) > 0) {
    p_scatter <- ggplot(scatter_data, aes(x = direct_score, y = nichenet_activity)) +
      geom_point(aes(color = cell_type), alpha = 0.7, size = 2) +
      geom_text_repel(
        data = scatter_data[, .SD[which.max(nichenet_activity)], by = cell_type],
        aes(label = paste0(ligand, "-", receptor)),
        size = 2.5, max.overlaps = 15
      ) +
      theme_bw(base_size = 11) +
      labs(title = "Direct L-R Score vs NicheNet Ligand Activity",
           x = "Direct Score (gradient × expression)",
           y = "NicheNet Activity (AUPR)",
           color = "Cell Type")

    ggsave(file.path(PLOT_DIR, "method_comparison_scatter.pdf"),
           p_scatter, width = 10, height = 8)
    message("  Saved: method_comparison_scatter.pdf")
  }

  # ===========================================================================
  # Cross-reference Agreement Plot (LEC)
  # ===========================================================================

  if (!is.null(crossref_lec)) {
    shared_lec <- crossref_lec[source == "both" &
                                 !is.na(combined_score) &
                                 !is.na(nichenet_prioritization)]

    if (nrow(shared_lec) >= 3) {
      p_crossref <- ggplot(shared_lec, aes(x = combined_score,
                                           y = nichenet_prioritization)) +
        geom_point(size = 3, alpha = 0.7, color = "#66C2A5") +
        geom_text_repel(aes(label = paste0(ligand, "-", receptor)),
                        size = 2.5, max.overlaps = 15) +
        geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "dashed") +
        theme_bw(base_size = 11) +
        labs(title = "LEC: Gradient-based vs NicheNet Prioritization",
             subtitle = sprintf("n=%d shared pairs", nrow(shared_lec)),
             x = "Gradient Combined Score (this analysis)",
             y = "NicheNet Prioritization Score (existing)")

      ggsave(file.path(PLOT_DIR, "nichenet_crossref_lec.pdf"),
             p_crossref, width = 8, height = 7)
      message("  Saved: nichenet_crossref_lec.pdf")
    }
  }

} else {
  message("  No combined results to plot.")
}

# =============================================================================
# Final Summary
# =============================================================================

message("\n", strrep("=", 70))
message("Atlas Complete")
message(strrep("=", 70))
message(sprintf("Summary tables: %s", SUMMARY_DIR))
message(sprintf("Figures: %s", PLOT_DIR))

# Print key numbers
if (nrow(all_combined) > 0) {
  message("\n--- Key Results ---")
  message(sprintf("Total L-R pairs (Direction A): %d", nrow(all_combined)))
  n_clean <- sum(all_combined$artifact_flag == "clean", na.rm = TRUE)
  n_suspect <- sum(all_combined$artifact_flag == "suspect", na.rm = TRUE)
  n_artifact <- sum(all_combined$artifact_flag == "artifact", na.rm = TRUE)
  n_lowconf <- sum(all_combined$artifact_flag == "low_confidence", na.rm = TRUE)
  message(sprintf("  Clean: %d, Low-confidence: %d, Suspect: %d, Artifact: %d",
                  n_clean, n_lowconf, n_suspect, n_artifact))
  message(sprintf("Cell types with L-R pairs: %d",
                  length(unique(all_combined$cell_type))))
  message(sprintf("Unique ligands (clean): %d",
                  length(unique(all_combined[artifact_flag == "clean"]$ligand))))
  message(sprintf("Unique receptors (clean): %d",
                  length(unique(all_combined[artifact_flag == "clean"]$receptor))))

  if ("n_samples_supporting" %in% names(all_combined)) {
    n_reprod_3plus <- sum(all_combined$artifact_flag == "clean" &
                          all_combined$n_samples_supporting >= 3, na.rm = TRUE)
    message(sprintf("Clean pairs reproduced in ≥3/4 mice: %d / %d",
                    n_reprod_3plus, n_clean))
  }

  if (nrow(all_enrichment) > 0) {
    n_enriched <- sum(all_enrichment$fdr_fisher < 0.05, na.rm = TRUE)
    message(sprintf("Enriched ligands (FDR<0.05): %d / %d",
                    n_enriched, nrow(all_enrichment)))
  }
}
if (nrow(all_direct_b) > 0) {
  message(sprintf("Total L-R pairs (Direction B: Target->%s): %d",
                  QUERY_LABEL, nrow(all_direct_b)))
}

message(sprintf("\nTimestamp: %s", Sys.time()))
