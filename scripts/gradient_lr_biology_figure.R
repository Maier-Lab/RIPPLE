#!/usr/bin/env Rscript
# ==============================================================================
# Curated Biological L-R Pairs Figure — Direction A + B Side-by-Side
# ==============================================================================
# Creates a publication-quality figure highlighting biologically informative
# ligand-receptor pairs, grouped by functional theme.
#
# Direction A: HyMy → Target (HyMy ligands acting on gradient receptors)
# Direction B: Target → HyMy (target ligands induced near HyMy, acting on HyMy)
#
# Author: CMM Project
# Date: 2026-03
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# --- Paths ---
script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg)))
    } else {
      getwd()
    }
  }
)
source(file.path(script_dir, "utils.R"))

ANNOTATION_LEVEL <- Sys.getenv("ANNOTATION_LEVEL", "HyMy")
GRADIENT_SOURCE <- Sys.getenv("GRADIENT_SOURCE", unset = "hymy_distance_correlation")
gradient_suffix <- sub("^hymy_distance_correlation", "", GRADIENT_SOURCE)
output_name <- if (nchar(gradient_suffix) > 0) {
  paste0("gradient_lr_integration", gradient_suffix)
} else {
  "gradient_lr_integration"
}

if (ANNOTATION_LEVEL == "L1") {
  base_dir <- file.path(PROJECT_ROOT, "results", "spatial_analysis_L1", output_name)
} else {
  base_dir <- file.path(PROJECT_ROOT, "results", "spatial_analysis", output_name)
}

out_dir <- file.path(base_dir, "plots")
ensure_dir(out_dir)

# ==============================================================================
# Load data — both directions
# ==============================================================================

cat("Loading combined L-R results...\n")
dt_a <- fread(file.path(base_dir, "summary", "all_lr_pairs_combined.csv"))
clean_a <- dt_a[artifact_flag == "clean"]
cat("  Direction A (HyMy->Target) clean pairs:", nrow(clean_a), "\n")

dt_b <- fread(file.path(base_dir, "summary", "all_target_to_hymy_pairs.csv"))
# Direction B has no artifact_flag — filter by negative ligand_gradient_coef (induced)
clean_b <- dt_b[ligand_gradient_coef < 0]
cat("  Direction B (Target->HyMy) induced pairs:", nrow(clean_b), "\n")

# ==============================================================================
# Define curated biologically informative pairs — Direction A (20 pairs)
# ==============================================================================

curated_a <- rbindlist(list(
  # --- Immune Checkpoint (2) ---
  data.table(ligand = "Cd274",  receptor = "Pdcd1",  cell_type = "CD8_T_cells",
             theme = "Immune checkpoint"),
  data.table(ligand = "Lgals3", receptor = "Anxa2",  cell_type = "CD8_T_cells",
             theme = "Immune checkpoint"),

  # --- Phagocytic / Don't-Eat-Me (3) ---
  data.table(ligand = "Sirpb1a", receptor = "Cd47",  cell_type = "CD8_T_cells",
             theme = "Phagocytic signaling"),
  data.table(ligand = "Cd47",    receptor = "Sirpa", cell_type = "CD4_T_cells",
             theme = "Phagocytic signaling"),
  data.table(ligand = "Csf1",    receptor = "Sirpa", cell_type = "CD4_T_cells",
             theme = "Phagocytic signaling"),

  # --- Chemokine Recruitment (4) ---
  data.table(ligand = "Ccl4",   receptor = "Ccr5",  cell_type = "CD4_T_cells",
             theme = "Chemokine recruitment"),
  data.table(ligand = "Ccl6",   receptor = "Ccr1",  cell_type = "Monocyte",
             theme = "Chemokine recruitment"),
  data.table(ligand = "Ccl8",   receptor = "Ccr1",  cell_type = "Monocyte",
             theme = "Chemokine recruitment"),
  data.table(ligand = "Cxcl16", receptor = "Cxcr6", cell_type = "CD8_T_cells",
             theme = "Chemokine recruitment"),

  # --- Vascular Trafficking (2) ---
  data.table(ligand = "Selplg", receptor = "Selp",  cell_type = "BEC",
             theme = "Vascular trafficking"),
  data.table(ligand = "Selplg", receptor = "Sele",  cell_type = "BEC",
             theme = "Vascular trafficking"),

  # --- IL-1 / Alarmin Signaling (2) ---
  data.table(ligand = "Il1b",   receptor = "Il1r2",  cell_type = "mature_migDC",
             theme = "IL-1 / alarmin signaling"),
  data.table(ligand = "Il1b",   receptor = "Il1r1",  cell_type = "CD8_T_cells",
             theme = "IL-1 / alarmin signaling"),

  # --- TGF-beta / Immunosuppression (3) ---
  data.table(ligand = "Tgfb1",  receptor = "Tgfbr1", cell_type = "FRC",
             theme = "TGF-\u03B2 immunosuppression"),
  data.table(ligand = "Tgfbi",  receptor = "Itgb1",  cell_type = "CD8_T_cells",
             theme = "TGF-\u03B2 immunosuppression"),
  data.table(ligand = "Nampt",  receptor = "Itgb1",  cell_type = "CD8_T_cells",
             theme = "TGF-\u03B2 immunosuppression"),

  # --- ECM / Tissue Remodeling (2) ---
  data.table(ligand = "Mmp9",   receptor = "Cd44",   cell_type = "CD8_T_cells",
             theme = "ECM remodeling"),
  data.table(ligand = "Mmp14",  receptor = "Sdc1",   cell_type = "cDC1",
             theme = "ECM remodeling"),

  # --- Stromal Niche (2) ---
  data.table(ligand = "Ly86",   receptor = "Cd180",  cell_type = "mature_migDC",
             theme = "Stromal niche"),
  data.table(ligand = "Sell",   receptor = "Podxl",  cell_type = "LEC",
             theme = "Stromal niche")
))

# ==============================================================================
# Define curated biologically informative pairs — Direction B (9 pairs)
# ==============================================================================

curated_b <- rbindlist(list(
  # --- Immune Evasion (1) ---
  data.table(ligand = "Gzmb",   receptor = "Mcl1",   cell_type = "CD8_T_cells",
             theme = "Immune evasion"),

  # --- Inflammatory Amplification (1) ---
  data.table(ligand = "Il17a",  receptor = "Il17ra", cell_type = "gdT_cells",
             theme = "Inflammatory amplification"),

  # --- Anti-inflammatory Feedback (1) ---
  data.table(ligand = "Il10",   receptor = "Il10ra", cell_type = "Plasma_cell",
             theme = "Anti-inflammatory feedback"),

  # --- Lymphatic Homing (1) ---
  data.table(ligand = "Cxcl12", receptor = "Cxcr4",  cell_type = "LEC",
             theme = "Lymphatic homing"),

  # --- Survival Signaling (1) ---
  data.table(ligand = "Csf3",   receptor = "Csf3r",  cell_type = "FRC",
             theme = "Survival signaling"),

  # --- Vascular Tethering (1) ---
  data.table(ligand = "Selp",   receptor = "Selplg", cell_type = "BEC",
             theme = "Vascular tethering"),

  # --- Myeloid Cross-recruitment (1) ---
  data.table(ligand = "Ccl2",   receptor = "Ccr1",   cell_type = "cDC2",
             theme = "Myeloid recruitment"),

  # --- IL-1 Decoy Feedback (1) ---
  data.table(ligand = "Il1b",   receptor = "Il1r2",  cell_type = "FRC",
             theme = "IL-1 decoy feedback"),

  # --- ECM Tethering (1) ---
  data.table(ligand = "Vcan",   receptor = "Selplg", cell_type = "FRC",
             theme = "ECM tethering")
))

# ==============================================================================
# Shared display settings
# ==============================================================================

ct_order <- c("CD8_T_cells", "CD4_T_cells", "Monocyte", "cDC1", "cDC2",
              "mature_migDC", "FRC", "LEC", "BEC", "B_cells",
              "gdT_cells", "Macrophages", "Fibroblasts_mac", "Plasma_cell")

ct_labels <- c(
  "CD8_T_cells" = "CD8 T", "CD4_T_cells" = "CD4 T", "Monocyte" = "Mono",
  "cDC1" = "cDC1", "cDC2" = "cDC2", "mature_migDC" = "migDC",
  "FRC" = "FRC", "LEC" = "LEC", "BEC" = "BEC", "B_cells" = "B cell",
  "gdT_cells" = "\u03B3\u03B4T", "Macrophages" = "Mac",
  "Fibroblasts_mac" = "Fibro", "Plasma_cell" = "Plasma"
)

theme_colors_a <- c(
  "Immune checkpoint"             = "#E63946",
  "Phagocytic signaling"          = "#F4845F",
  "Chemokine recruitment"         = "#457B9D",
  "Vascular trafficking"          = "#2A9D8F",
  "IL-1 / alarmin signaling"      = "#E9C46A",
  "TGF-\u03B2 immunosuppression"  = "#264653",
  "ECM remodeling"                = "#A8DADC",
  "Stromal niche"                 = "#B5838D"
)

theme_colors_b <- c(
  "Immune evasion"              = "#9B2226",
  "Inflammatory amplification"  = "#CA6702",
  "Anti-inflammatory feedback"  = "#94D2BD",
  "Lymphatic homing"            = "#0A9396",
  "Survival signaling"          = "#EE9B00",
  "Vascular tethering"          = "#005F73",
  "Myeloid recruitment"         = "#AE2012",
  "IL-1 decoy feedback"         = "#BB3E03",
  "ECM tethering"               = "#E9D8A6"
)

# ==============================================================================
# Merge curated pairs with data
# ==============================================================================

merge_curated <- function(curated, data_dt, direction_label) {
  curated[, pair_key := paste(ligand, receptor, cell_type, sep = "_")]
  data_dt[, pair_key := paste(ligand, receptor, cell_type, sep = "_")]
  merged <- merge(curated, data_dt, by = "pair_key", suffixes = c("", ".data"))

  # Keep curated theme, use data columns for everything else
  drop_cols <- grep("\\.data$", names(merged), value = TRUE)
  if (length(drop_cols) > 0) merged[, (drop_cols) := NULL]

  cat(direction_label, "- matched:", nrow(merged), "of", nrow(curated), "\n")
  if (nrow(merged) < nrow(curated)) {
    missing <- curated[!pair_key %in% merged$pair_key]
    cat("  Missing:", paste(missing$pair_key, collapse = ", "), "\n")
  }
  merged
}

cat("\nMerging curated pairs...\n")
merged_a <- merge_curated(curated_a, clean_a, "Direction A")
merged_b <- merge_curated(curated_b, clean_b, "Direction B")

# ==============================================================================
# Prepare Direction A for plotting
# ==============================================================================

merged_a[, pair_label := paste0(ligand, " \u2192 ", receptor)]

theme_order_a <- c(
  "Immune checkpoint",
  "Phagocytic signaling",
  "Chemokine recruitment",
  "Vascular trafficking",
  "IL-1 / alarmin signaling",
  "TGF-\u03B2 immunosuppression",
  "ECM remodeling",
  "Stromal niche"
)
merged_a[, theme := factor(theme, levels = rev(theme_order_a))]
merged_a <- merged_a[order(theme, -combined_score)]
merged_a[, pair_label := factor(pair_label, levels = rev(unique(pair_label)))]
merged_a[, cell_type := factor(cell_type, levels = ct_order)]

# ==============================================================================
# Prepare Direction B for plotting
# ==============================================================================

merged_b[, pair_label := paste0(ligand, " \u2192 ", receptor)]

theme_order_b <- c(
  "Immune evasion",
  "Inflammatory amplification",
  "Anti-inflammatory feedback",
  "Lymphatic homing",
  "Survival signaling",
  "Vascular tethering",
  "Myeloid recruitment",
  "IL-1 decoy feedback",
  "ECM tethering"
)
merged_b[, theme := factor(theme, levels = rev(theme_order_b))]
merged_b <- merged_b[order(theme, -direct_score)]
merged_b[, pair_label := factor(pair_label, levels = rev(unique(pair_label)))]
merged_b[, cell_type := factor(cell_type, levels = ct_order)]

# ==============================================================================
# Figure 1: Direction A dotplot (standalone)
# ==============================================================================

cat("\nGenerating curated biology figures...\n")

p_dot_a <- ggplot(merged_a, aes(x = cell_type, y = pair_label)) +
  geom_point(aes(size = combined_score, fill = theme),
             shape = 21, color = "gray30", stroke = 0.4) +
  scale_size_continuous(range = c(2, 10), name = "Combined\nscore",
                        breaks = c(0.2, 0.4, 0.6)) +
  scale_fill_manual(values = theme_colors_a, name = "Theme") +
  scale_x_discrete(labels = ct_labels) +
  labs(x = NULL, y = NULL,
       title = "A: HyMy \u2192 Target",
       subtitle = "Bubble size = combined score") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    panel.grid.major = element_line(color = "gray92"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray40")
  ) +
  facet_grid(theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9,
                                         face = "bold.italic"))

ggsave(file.path(out_dir, "curated_biology_dotplot.pdf"), p_dot_a,
       width = 11, height = 10, device = cairo_pdf)
ggsave(file.path(out_dir, "curated_biology_dotplot.png"), p_dot_a,
       width = 11, height = 10, dpi = 200)
cat("  Saved: curated_biology_dotplot.pdf/png (Direction A only)\n")

# ==============================================================================
# Figure 2: Combined side-by-side — Direction A + B with evidence panels
# ==============================================================================

# --- Helper: create dotplot panel ---
make_dotplot <- function(merged, theme_colors, score_col, score_label, title) {
  ggplot(merged, aes(x = cell_type, y = pair_label)) +
    geom_point(aes(size = .data[[score_col]], fill = theme),
               shape = 21, color = "gray30", stroke = 0.4) +
    scale_size_continuous(range = c(2, 9), name = score_label,
                          breaks = pretty(range(merged[[score_col]], na.rm = TRUE), 3)) +
    scale_fill_manual(values = theme_colors, name = "Theme") +
    scale_x_discrete(labels = ct_labels) +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid.major = element_line(color = "gray92"),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.box = "vertical",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      plot.title = element_text(size = 13, face = "bold")
    ) +
    facet_grid(theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8,
                                           face = "bold.italic"))
}

# --- Helper: create gradient lollipop panel ---
make_gradient_panel <- function(merged, gradient_col, theme_colors, title) {
  merged[, grad_display := get(gradient_col) * 1000]
  ggplot(merged, aes(x = grad_display, y = pair_label)) +
    geom_segment(aes(xend = 0, yend = pair_label), linewidth = 0.6,
                 color = "#457B9D") +
    geom_point(size = 2.5, color = "#457B9D") +
    geom_vline(xintercept = 0, color = "gray50", linewidth = 0.3) +
    labs(x = "Coef (\u00D710\u207B\u00B3)", y = NULL, title = title) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(size = 9),
      axis.text.y = element_blank(),
      axis.title.x = element_text(size = 9),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 11, face = "bold")
    )
}

# --- Helper: create reproducibility panel ---
make_repro_panel <- function(merged, theme_colors) {
  ggplot(merged, aes(x = n_samples_supporting, y = pair_label, fill = theme)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    scale_fill_manual(values = theme_colors) +
    scale_x_continuous(breaks = 1:4, limits = c(0, 4.2)) +
    labs(x = "Mice (of 4)", y = NULL, title = "Repro.") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(size = 9),
      axis.text.y = element_blank(),
      axis.title.x = element_text(size = 9),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 11, face = "bold")
    )
}

# --- Build panels ---
p_a_dot  <- make_dotplot(merged_a, theme_colors_a, "combined_score",
                          "Combined\nscore", "A: HyMy \u2192 Target")
p_a_grad <- make_gradient_panel(merged_a, "receptor_gradient_coef",
                                 theme_colors_a, "Receptor\ngradient")
p_a_rep  <- make_repro_panel(merged_a, theme_colors_a)

p_b_dot  <- make_dotplot(merged_b, theme_colors_b, "direct_score",
                          "Direct\nscore", "B: Target \u2192 HyMy")
p_b_grad <- make_gradient_panel(merged_b, "ligand_gradient_coef",
                                 theme_colors_b, "Ligand\ngradient")
p_b_rep  <- make_repro_panel(merged_b, theme_colors_b)

# --- Combine: top row = Dir A, bottom row = Dir B ---
row_a <- p_a_dot + p_a_grad + p_a_rep + plot_layout(widths = c(5, 1.5, 0.8))
row_b <- p_b_dot + p_b_grad + p_b_rep + plot_layout(widths = c(5, 1.5, 0.8))

n_a <- nrow(merged_a)
n_b <- nrow(merged_b)
n_mice_a <- sum(merged_a$n_samples_supporting == 4)
n_mice_b <- sum(merged_b$n_samples_supporting == 4)

p_combined <- (row_a / row_b) +
  plot_layout(heights = c(n_a, n_b)) +
  plot_annotation(
    title = "Gradient-to-LR Integration: Curated Biological Pairs",
    subtitle = paste0(
      "Direction A: ", n_a, " pairs (", n_mice_a, "/", n_a, " in all 4 mice) | ",
      "Direction B: ", n_b, " pairs (", n_mice_b, "/", n_b, " in all 4 mice)"
    ),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray40")
    )
  )

ggsave(file.path(out_dir, "curated_biology_combined.pdf"), p_combined,
       width = 18, height = 18, device = cairo_pdf)
ggsave(file.path(out_dir, "curated_biology_combined.png"), p_combined,
       width = 18, height = 18, dpi = 200)
cat("  Saved: curated_biology_combined.pdf/png\n")

# ==============================================================================
# Figure 3: Bidirectional network
# ==============================================================================

# Collect all cell types from both directions
all_cts <- unique(c(as.character(merged_a$cell_type),
                     as.character(merged_b$cell_type)))

# Best pair per cell type per direction
best_a <- merged_a[, .SD[which.max(combined_score)], by = cell_type]
best_a[, dir := "A"]
best_b <- merged_b[, .SD[which.max(direct_score)], by = cell_type]
best_b[, dir := "B"]

# Angular positions for unique cell types
n_ct <- length(all_cts)
ct_positions <- data.table(
  cell_type = all_cts,
  angle = seq(0, 2 * pi * (1 - 1/n_ct), length.out = n_ct)
)
ct_positions[, x := cos(angle) * 3]
ct_positions[, y := sin(angle) * 3]
ct_positions[, lx := cos(angle) * 3.8]
ct_positions[, ly := sin(angle) * 3.8]

# Merge positions
best_a <- merge(best_a, ct_positions, by = "cell_type")
best_b <- merge(best_b, ct_positions, by = "cell_type")

p_network <- ggplot() +
  # Direction A: solid arrows outward (HyMy → Target)
  geom_segment(data = best_a,
               aes(x = 0, y = 0, xend = x * 0.7, yend = y * 0.7),
               arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
               linewidth = 0.7, color = "#457B9D") +
  # Direction B: dashed arrows inward (Target → HyMy)
  geom_segment(data = best_b,
               aes(x = x * 0.7, y = y * 0.7, xend = 0, yend = 0),
               arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
               linewidth = 0.6, color = "#D62828", linetype = "dashed") +
  # Cell type nodes
  geom_point(data = ct_positions, aes(x = x, y = y),
             size = 12, shape = 21, fill = "white", color = "gray30", stroke = 1) +
  geom_text(data = ct_positions, aes(x = x, y = y,
                                      label = ct_labels[cell_type]),
            size = 2.8, fontface = "bold") +
  # Direction A labels (outside, blue)
  geom_label(data = best_a,
             aes(x = x * 0.42, y = y * 0.42,
                 label = paste0(ligand, "\u2192", receptor)),
             size = 2.2, fill = "#EBF1F5", label.size = 0.15,
             label.padding = unit(0.12, "lines"), color = "#457B9D") +
  # Direction B labels (outside, red)
  geom_label(data = best_b,
             aes(x = x * 0.88, y = y * 0.88,
                 label = paste0(ligand, "\u2192", receptor)),
             size = 2.2, fill = "#FCEAEA", label.size = 0.15,
             label.padding = unit(0.12, "lines"), color = "#D62828") +
  # HyMy center
  annotate("point", x = 0, y = 0, size = 18, shape = 21,
           fill = "#FC8D62", color = "gray30", stroke = 1.5) +
  annotate("text", x = 0, y = 0, label = "HyMy", size = 4.5, fontface = "bold") +
  coord_fixed(xlim = c(-5, 5), ylim = c(-5, 5)) +
  labs(title = "HyMy Bidirectional Communication Network",
       subtitle = paste0("Solid blue = HyMy\u2192Target (Dir A) | ",
                         "Dashed red = Target\u2192HyMy (Dir B)")) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9, color = "gray40", hjust = 0.5)
  )

ggsave(file.path(out_dir, "curated_biology_network.pdf"), p_network,
       width = 10, height = 10, device = cairo_pdf)
ggsave(file.path(out_dir, "curated_biology_network.png"), p_network,
       width = 10, height = 10, dpi = 200)
cat("  Saved: curated_biology_network.pdf/png\n")

# ==============================================================================
# Summary tables
# ==============================================================================

# Direction A table
table_a <- merged_a[order(theme, -combined_score), .(
  Direction = "A: HyMy->Target",
  Theme = theme,
  Ligand = ligand,
  Receptor = receptor,
  `Cell type` = cell_type,
  `Combined score` = round(combined_score, 3),
  `Gradient coef` = round(receptor_gradient_coef, 5),
  `NicheNet AUPR` = round(nichenet_activity, 4),
  `Ligand % HyMy` = round(ligand_pct_hymy, 1),
  `Receptor % target` = round(receptor_pct_target * 100, 1),
  `Mice (of 4)` = n_samples_supporting,
  `Stage 2` = stage2_classification,
  `Enrichment OR` = round(enrichment_odds_ratio, 2)
)]

# Direction B table
table_b <- merged_b[order(theme, -direct_score), .(
  Direction = "B: Target->HyMy",
  Theme = theme,
  Ligand = ligand,
  Receptor = receptor,
  `Cell type` = cell_type,
  `Direct score` = round(direct_score, 5),
  `Gradient coef` = round(ligand_gradient_coef, 5),
  `Ligand % target` = round(ligand_pct_target * 100, 1),
  `Receptor % HyMy` = round(receptor_pct_hymy, 1),
  `Mice (of 4)` = n_samples_supporting
)]

fwrite(table_a, file.path(out_dir, "curated_biology_table_dirA.csv"))
fwrite(table_b, file.path(out_dir, "curated_biology_table_dirB.csv"))

# Combined table (with fill=TRUE for mismatched columns)
table_combined <- rbindlist(list(table_a, table_b), fill = TRUE)
fwrite(table_combined, file.path(out_dir, "curated_biology_table.csv"))
cat("  Saved: curated_biology_table*.csv\n")

# ==============================================================================
# Summary
# ==============================================================================

cat("\n=== CURATED PAIRS SUMMARY ===\n")
cat("Direction A (HyMy->Target):\n")
cat("  Total pairs:", nrow(merged_a), "\n")
cat("  Themes:", length(unique(merged_a$theme)), "\n")
cat("  Cell types:", length(unique(merged_a$cell_type)), "\n")
cat("  All 4 mice:", n_mice_a, "/", n_a, "\n")

cat("Direction B (Target->HyMy):\n")
cat("  Total pairs:", nrow(merged_b), "\n")
cat("  Themes:", length(unique(merged_b$theme)), "\n")
cat("  Cell types:", length(unique(merged_b$cell_type)), "\n")
cat("  All 4 mice:", n_mice_b, "/", n_b, "\n")

cat("\nDone.\n")
