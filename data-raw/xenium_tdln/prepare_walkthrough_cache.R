# ============================================================================
# Prepare cached data for the naive TRC T zone walkthrough figure
# ============================================================================
# Produces two small CSVs used by the paper's Fig. 3 walkthrough panels:
#
#   1. k_diagnostics_naive_trc.csv
#        Per-target-cell-type k-diagnostic summary for the naive TRC-Ccl21a
#        query. Feeds the "is TRC-Ccl21a well-distributed?" panel.
#
#   2. decay_curves_naive_trc_tcell.csv
#        Long-format per-gene per-mouse binned decay data (distance bin,
#        mean rate, n_cells) for six marker genes in T_cell_all:
#        Ccr7, Sell, Tcf7, Cxcr4, Gzma (positive controls) and
#        Ccl21a (ambient control). Feeds the decay-curve panel.
#
# Both are written to inst/extdata/naive_trc_cached/. The raw Seurat object
# stays in data-raw/ (gitignored) and is not bundled.
#
# Run with: Rscript data-raw/xenium_tdln/prepare_walkthrough_cache.R
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(RANN)
  devtools::load_all(quiet = TRUE)
})

L1_COL     <- "cell_type_assignment_L1"
L2_COL     <- "cell_type_assignment_L2"
SAMPLE_COL <- "mouse_id"
QUERY_CT   <- "TRC-Ccl21a"
TARGET_CT  <- "T_cell_all"
T_CELL_L1  <- c("Activated_CD8", "Cytotoxic_CD8", "gdT_cell",
                "Naive_CD4", "Naive_CD8", "Tfh", "Tpex", "Treg")
DECAY_GENES <- c("Ccr7", "Sell", "Tcf7", "Cxcr4", "Gzma", "Ccl21a")

out_dir <- "inst/extdata/naive_trc_cached"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load + apply L2 + subset naive -----------------------------------
obj_path <- "data-raw/xenium_tdln/seurat_xenium_filtered.rds"
ann_path <- "data-raw/xenium_tdln/frc_subcluster_annotation.csv"

cat("Loading Seurat object...\n")
obj <- readRDS(obj_path)

ann <- fread(ann_path)
setnames(ann, c("index", "frc_subcluster"))
ann[, frc_subcluster_norm := gsub(" ", "_", frc_subcluster)]
sub_lookup <- setNames(ann$frc_subcluster_norm, ann$index)

md <- obj@meta.data
cell_keys <- rownames(md)
overwrite_idx <- which(cell_keys %in% ann$index)
final <- as.character(md[[L1_COL]])
final[overwrite_idx] <- sub_lookup[cell_keys[overwrite_idx]]
final[final %in% T_CELL_L1] <- TARGET_CT
obj@meta.data[[L2_COL]] <- final

naive_cells <- rownames(obj@meta.data)[obj@meta.data$group == "naive"]
obj <- subset(obj, cells = naive_cells)
cat("Naive subset:", ncol(obj), "cells\n")

# --- 2. Extract coordinates and cell type/mouse info ---------------------
md <- obj@meta.data
coord_cols <- grep("^spatial_[xy]$", colnames(md), value = TRUE)
stopifnot(length(coord_cols) == 2)

cell_info <- data.table(
  barcode   = rownames(md),
  cell_type = as.character(md[[L2_COL]]),
  mouse_id  = as.character(md[[SAMPLE_COL]]),
  x         = md[[coord_cols[1]]],
  y         = md[[coord_cols[2]]]
)

# --- 3. k-diagnostics per target cell type --------------------------------
cat("\nRunning k-diagnostics for TRC-Ccl21a...\n")

# Target populations we care about for the walkthrough figure
targets_of_interest <- c(TARGET_CT, "B_cell", "cDC1", "Macrophages", "LEC")
kdiag_rows <- list()

for (samp in unique(cell_info$mouse_id)) {
  samp_info <- cell_info[mouse_id == samp]
  q_xy <- as.matrix(samp_info[cell_type == QUERY_CT, .(x, y)])
  if (nrow(q_xy) < 10) next

  for (tgt in targets_of_interest) {
    t_xy <- as.matrix(samp_info[cell_type == tgt, .(x, y)])
    if (nrow(t_xy) < 30) next

    kmax <- min(20, nrow(q_xy))
    nn <- RANN::nn2(q_xy, t_xy, k = kmax)
    # For each k, per-cell distance = mean of nn.dists columns 1:k
    for (k in 1:kmax) {
      d <- if (k == 1) nn$nn.dists[, 1] else rowMeans(nn$nn.dists[, 1:k, drop = FALSE])
      kdiag_rows[[paste(samp, tgt, k, sep = "_")]] <- data.table(
        mouse_id  = samp,
        cell_type = tgt,
        k         = k,
        mean_dist = mean(d),
        sd_dist   = sd(d),
        n_cells   = length(d)
      )
    }
  }
  cat("  ", samp, "done\n")
}
kdiag_dt <- rbindlist(kdiag_rows)
fwrite(kdiag_dt, file.path(out_dir, "k_diagnostics_naive_trc.csv"))
cat("Saved", file.path(out_dir, "k_diagnostics_naive_trc.csv"),
    sprintf("(%d rows)\n", nrow(kdiag_dt)))

# --- 4. Binned decay data for selected genes in T_cell_all ---------------
cat("\nComputing binned decay curves for",
    paste(DECAY_GENES, collapse = ", "), "...\n")

counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
all_genes <- rownames(counts_mat)
available <- DECAY_GENES[DECAY_GENES %in% all_genes]
missing   <- setdiff(DECAY_GENES, all_genes)
if (length(missing) > 0) {
  cat("  MISSING:", paste(missing, collapse = ", "), "\n")
}

# Restrict to T_cell_all cells; compute dist to TRC-Ccl21a per mouse
tcell_rows <- list()
n_bins <- 20
max_dist <- 200
bins <- seq(0, max_dist, length.out = n_bins + 1)
bin_mid <- (head(bins, -1) + tail(bins, -1)) / 2

for (samp in unique(cell_info$mouse_id)) {
  s_info <- cell_info[mouse_id == samp]
  q_xy <- as.matrix(s_info[cell_type == QUERY_CT, .(x, y)])
  t_info <- s_info[cell_type == TARGET_CT]
  if (nrow(q_xy) < 10 || nrow(t_info) < 30) next

  t_xy <- as.matrix(t_info[, .(x, y)])
  nn <- RANN::nn2(q_xy, t_xy, k = 1)
  dist_q <- as.vector(nn$nn.dists)

  keep <- dist_q <= max_dist
  t_bc   <- t_info$barcode[keep]
  dist_q <- dist_q[keep]
  total_counts <- Matrix::colSums(counts_mat[, t_bc])
  total_counts[total_counts == 0] <- 1

  bin_idx <- findInterval(dist_q, bins, rightmost.closed = TRUE)
  bin_idx <- pmax(pmin(bin_idx, n_bins), 1)

  for (g in available) {
    y <- as.numeric(counts_mat[g, t_bc])
    rate <- y / total_counts  # per-UMI-normalised rate
    for (b in seq_len(n_bins)) {
      mask <- bin_idx == b
      if (sum(mask) < 5) next
      tcell_rows[[paste(samp, g, b, sep = "_")]] <- data.table(
        mouse_id   = samp,
        gene       = g,
        bin        = b,
        bin_mid_um = bin_mid[b],
        mean_rate  = mean(rate[mask]),
        se_rate    = sd(rate[mask]) / sqrt(sum(mask)),
        n_cells    = sum(mask)
      )
    }
  }
  cat("  ", samp, "done\n")
}

decay_dt <- rbindlist(tcell_rows)
fwrite(decay_dt, file.path(out_dir, "decay_curves_naive_trc_tcell.csv"))
cat("Saved", file.path(out_dir, "decay_curves_naive_trc_tcell.csv"),
    sprintf("(%d rows, %d genes)\n",
            nrow(decay_dt), length(unique(decay_dt$gene))))
cat("\nAll done.\n")
