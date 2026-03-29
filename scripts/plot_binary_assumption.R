#!/usr/bin/env Rscript
# Quick QC plot: test the binary expression assumption for Xenium 5K
# Plots nFeature vs nCount per sample — if slope ≈ 1, most genes have 1 count

library(Seurat)
library(ggplot2)
library(patchwork)

# --- Load data ---
seurat_path <- "N:/lab_maier/Projects/mXenium/results/cell_type_assignment/HyMy_annotation/seurat_xenium_filtered.rds"
out_dir <- "N:/lab_maier/Projects/mXenium/CMM/results/spatial_analysis/hymy_distance_correlation/qc"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Loading Seurat object...")
obj <- readRDS(seurat_path)
message("Loaded: ", ncol(obj), " cells x ", nrow(obj), " genes")

# Extract metadata
meta <- obj@meta.data
meta$sample <- meta$sample_id  # adjust if different column name

# Check column names
if (!"sample_id" %in% colnames(meta)) {
  # Try common alternatives
  sample_col <- intersect(c("sample_id", "sample", "orig.ident", "Sample"), colnames(meta))[1]
  message("Using '", sample_col, "' as sample column")
  meta$sample <- meta[[sample_col]]
} else {
  meta$sample <- meta$sample_id
}

message("Samples: ", paste(unique(meta$sample), collapse = ", "))

# Identify nCount/nFeature columns (may be _Xenium or _RNA etc.)
count_col <- grep("^nCount_", colnames(meta), value = TRUE)[1]
feat_col <- grep("^nFeature_", colnames(meta), value = TRUE)[1]
message("Using: ", count_col, ", ", feat_col)

# --- Plot 1: nCount vs nFeature scatter per sample ---
p1 <- ggplot(meta, aes(x = .data[[feat_col]], y = .data[[count_col]])) +
  geom_point(alpha = 0.02, size = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~sample, scales = "free") +
  labs(
    title = "nCount vs nFeature per cell (Xenium 5K)",
    subtitle = "Red dashed line = y=x (perfect binary: every gene has exactly 1 count)",
    x = "nFeature (genes detected)",
    y = "nCount (total counts)"
  ) +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 9))

# --- Plot 2: Distribution of counts-per-gene (excluding zeros) ---
# Sample a subset of cells to keep it manageable
set.seed(42)
n_sample <- min(50000, ncol(obj))
cell_idx <- sample(ncol(obj), n_sample)
expr_mat <- GetAssayData(obj, layer = "counts")[, cell_idx]

# Get all non-zero values
nonzero_vals <- expr_mat@x  # sparse matrix: only stored non-zero entries
message("Non-zero entries: ", length(nonzero_vals))
message("Counts == 1: ", sum(nonzero_vals == 1), " (",
        round(100 * sum(nonzero_vals == 1) / length(nonzero_vals), 1), "%)")
message("Counts == 2: ", sum(nonzero_vals == 2), " (",
        round(100 * sum(nonzero_vals == 2) / length(nonzero_vals), 1), "%)")
message("Counts >= 3: ", sum(nonzero_vals >= 3), " (",
        round(100 * sum(nonzero_vals >= 3) / length(nonzero_vals), 1), "%)")

count_df <- data.frame(counts = nonzero_vals)
# Cap at 10 for visualization
count_df$counts_capped <- pmin(count_df$counts, 10)

p2 <- ggplot(count_df, aes(x = factor(counts_capped))) +
  geom_bar(fill = "steelblue", color = "white") +
  scale_x_discrete(labels = c(1:9, "10+")) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Distribution of non-zero expression counts (50K cell subsample)",
    subtitle = paste0("Count=1: ", round(100 * sum(nonzero_vals == 1) / length(nonzero_vals), 1),
                      "% | Count=2: ", round(100 * sum(nonzero_vals == 2) / length(nonzero_vals), 1),
                      "% | Count>=3: ", round(100 * sum(nonzero_vals >= 3) / length(nonzero_vals), 1), "%"),
    x = "Counts per gene per cell",
    y = "Number of observations"
  ) +
  theme_bw(base_size = 11)

# --- Plot 3: Ratio nCount/nFeature per sample (should be ~1 if binary) ---
meta$counts_per_gene <- meta[[count_col]] / meta[[feat_col]]

p3 <- ggplot(meta, aes(x = counts_per_gene, fill = sample)) +
  geom_histogram(bins = 100, alpha = 0.7) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
  facet_wrap(~sample, scales = "free_y") +
  labs(
    title = "Average counts per gene per cell",
    subtitle = "Red line = 1.0 (perfect binary)",
    x = "nCount / nFeature",
    y = "Number of cells"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", strip.text = element_text(size = 9))

# --- Save as PNG (PDF too large with ~1M points) ---
png(file.path(out_dir, "binary_assumption_ncount_vs_nfeat.png"), width = 14, height = 10, units = "in", res = 150)
print(p1)
dev.off()

png(file.path(out_dir, "binary_assumption_count_distribution.png"), width = 10, height = 6, units = "in", res = 150)
print(p2)
dev.off()

png(file.path(out_dir, "binary_assumption_counts_per_gene.png"), width = 14, height = 10, units = "in", res = 150)
print(p3)
dev.off()

message("\nSaved 3 PNGs to: ", out_dir)
message("Done!")
