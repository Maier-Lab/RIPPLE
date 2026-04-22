# ============================================================================
# Naive LN validation v2 - FDC query with stromal subcluster annotation
# ============================================================================
# Uses frc_subcluster_annotation.csv to OVERWRITE the L1 cell type column
# for cells that were originally labeled FRC (or were in the stromal
# reclustering pool). This exposes FDC as a separate identity, enabling the
# canonical CXCL13/CXCR5 positive-control test.
#
# Additional rule: all B cell subtypes (L1 B_cell, Follicular_B_cell,
# GC_B_cell, AND any "B cell" cells coming out of the stromal reclustering)
# are collapsed into a single B_cell_all category to avoid circularity with
# follicular-marker-based subsetting.
#
# Query: FDC
# Expected biology in B_cell_all:
#   Induced near FDC (negative coef):   Cxcr5, Bcl6, Fcer2a (CD23),
#                                       Cr2 (CD21), Cd83
#   Repressed near FDC (positive coef): Gpr183 (EBI2, opposes Cxcr5)
#   Ambient-RNA check (produced by FDC, not B cell): Cxcl13
#
# Run with:
#   Rscript data-raw/xenium_tdln/run_ripple_naive_fdc.R
#
# Output:
#   data-raw/xenium_tdln/naive_fdc_ripple/
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  devtools::load_all(quiet = TRUE)
})

L1_COL            <- "cell_type_assignment_L1"
FINAL_COL         <- "cell_type_assignment_L2"
SAMPLE_COL        <- "mouse_id"
QUERY_CT          <- "FDC"
B_CELL_L1         <- c("B_cell", "Follicular_B_cell", "GC_B_cell")
# Any of these (from stromal subclustering) should also collapse into B_cell_all
B_CELL_SUBCLUSTER <- c("B cell")

obj_path <- "data-raw/xenium_tdln/seurat_xenium_filtered.rds"
ann_path <- "data-raw/xenium_tdln/frc_subcluster_annotation.csv"
out_dir  <- "data-raw/xenium_tdln/naive_fdc_ripple"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load Seurat + annotation CSV --------------------------------------
cat("Loading Seurat object...\n")
obj <- readRDS(obj_path)
cat("Full object:", ncol(obj), "cells x", nrow(obj), "genes\n")

ann <- fread(ann_path)
setnames(ann, c("index", "frc_subcluster"))
cat("Subcluster annotation rows:", nrow(ann), "\n")
print(sort(table(ann$frc_subcluster), decreasing = TRUE))

md <- obj@meta.data

# --- 2. Detect which metadata column matches the CSV's index column ------
candidate_cols <- list(
  rownames  = rownames(md),
  cell_id   = if ("cell_id" %in% names(md)) as.character(md$cell_id) else NULL,
  match_key = if ("match_key" %in% names(md)) as.character(md$match_key) else NULL
)
match_counts <- vapply(candidate_cols, function(vals) {
  if (is.null(vals)) 0L else sum(ann$index %in% vals)
}, integer(1))
cat("\nAnnotation index match counts:\n"); print(match_counts)

best_col <- names(match_counts)[which.max(match_counts)]
best_pct <- 100 * max(match_counts) / nrow(ann)
if (max(match_counts) < nrow(ann) * 0.80) {
  stop(sprintf("No Seurat column matches >=80%% of annotation rows. Max: %d/%d (%.1f%%) in '%s'",
               max(match_counts), nrow(ann), best_pct, best_col))
}
cat(sprintf("Using '%s' to merge annotation (%d/%d, %.1f%% match).\n",
            best_col, max(match_counts), nrow(ann), best_pct))
# Diagnostic: how many Seurat cells will be affected by the overwrite
cat(sprintf("Seurat cells receiving L2 label: %d/%d (%.1f%%)\n",
            sum(candidate_cols[[best_col]] %in% ann$index),
            length(candidate_cols[[best_col]]),
            100 * sum(candidate_cols[[best_col]] %in% ann$index) /
                   length(candidate_cols[[best_col]])))

# Key vector from Seurat cells (same length as ncol(obj))
cell_keys <- candidate_cols[[best_col]]

# --- 3. Build FINAL_COL: subcluster labels overwrite L1 where available ----
# Normalize "B cell" (space) -> "B_cell" so downstream collapsing is consistent
ann[, frc_subcluster_norm := gsub(" ", "_", frc_subcluster)]

# Map: keep L1 by default, overwrite with subcluster when in annotation
sub_lookup <- setNames(ann$frc_subcluster_norm, ann$index)
final <- as.character(md[[L1_COL]])
overwrite_idx <- which(cell_keys %in% ann$index)
final[overwrite_idx] <- sub_lookup[cell_keys[overwrite_idx]]

# Collapse B cell identities
b_cell_labels <- c(B_CELL_L1, B_CELL_SUBCLUSTER,
                   gsub(" ", "_", B_CELL_SUBCLUSTER))
final[final %in% b_cell_labels] <- "B_cell_all"

obj@meta.data[[FINAL_COL]] <- final

cat("\n=== Cell type composition after overwrite + B cell collapse ===\n")
print(sort(table(final), decreasing = TRUE))

# Persist the L2 annotation as a CSV for future re-use without re-running
# the merge/collapse logic. Includes both the original L1 label and the new
# L2 label per cell so downstream scripts can pick whichever they need.
l2_out <- data.table::data.table(
  cell_key = cell_keys,
  mouse_id = as.character(md[[SAMPLE_COL]]),
  group    = as.character(md$group),
  cell_type_assignment_L1 = as.character(md[[L1_COL]]),
  cell_type_assignment_L2 = final
)
l2_csv <- "data-raw/xenium_tdln/cell_type_assignment_L2.csv"
data.table::fwrite(l2_out, l2_csv)
cat("\nSaved L2 annotation:", l2_csv, "(", nrow(l2_out), "cells)\n")

# --- 4. Subset to naive ---------------------------------------------------
cat("\nSubsetting to naive...\n")
naive_cells <- rownames(obj@meta.data)[obj@meta.data$group == "naive"]
obj <- subset(obj, cells = naive_cells)
cat("Naive subset:", ncol(obj), "cells across",
    length(unique(obj@meta.data[[SAMPLE_COL]])), "mice\n")

md_n <- obj@meta.data

cat("\n=== Query (FDC) counts per naive mouse ===\n")
fdc_per_sample <- table(md_n[[SAMPLE_COL]][md_n[[FINAL_COL]] == QUERY_CT])
print(fdc_per_sample)
if (any(fdc_per_sample < 30)) {
  cat("\nWARNING: some samples have <30 FDCs, they may be excluded.\n")
}

cat("\n=== B_cell_all counts per naive mouse ===\n")
print(table(md_n[[SAMPLE_COL]][md_n[[FINAL_COL]] == "B_cell_all"]))

# --- 5. Marker gene availability ------------------------------------------
all_genes <- rownames(obj)
markers <- list(
  induced_expected   = c("Cxcr5", "Bcl6", "Fcer2a", "Cr2", "Cd83"),
  repressed_expected = c("Gpr183"),
  ambient_check      = c("Cxcl13")
)
cat("\n=== Marker gene presence ===\n")
for (group in names(markers)) {
  cat(sprintf("[%s]\n", group))
  for (g in markers[[group]]) {
    hits <- grep(paste0("^", g, "$"), all_genes, ignore.case = TRUE,
                 value = TRUE)
    cat(sprintf("  %-10s -> %s\n", g,
                if (length(hits)) paste(hits, collapse = ", ") else "MISSING"))
  }
}

# --- 6. Run RIPPLE --------------------------------------------------------
cat("\nRunning RIPPLE with FDC as query...\n")
t0 <- proc.time()
results <- run_ripple(
  input            = obj,
  query_celltype   = QUERY_CT,
  celltype_column  = FINAL_COL,
  sample_column    = SAMPLE_COL,
  output_dir       = out_dir,
  analysis_name    = "ripple",
  k_neighbors      = 1,
  max_distance_um  = 200,
  n_permutations   = 0,
  verbose          = TRUE
)
cat(sprintf("\nRIPPLE finished in %.1f minutes\n",
            (proc.time() - t0)[3] / 60))

# --- 7. Positive-control readouts -----------------------------------------
report_block <- function(label, genes, dt_subset) {
  cat(sprintf("\n=== %s ===\n", label))
  hits <- dt_subset[toupper(gene) %in% toupper(genes)]
  if (nrow(hits)) {
    print(hits[order(fisher_fdr),
               .(gene, median_coef, fisher_fdr, sign_consistency,
                 n_samples, pct_expr)])
  } else {
    cat("(none tested; likely expression-filtered)\n")
  }
}

bcell_res <- results[cell_type == "B_cell_all"]
cat(sprintf("\nB_cell_all genes tested: %d  |  significant (FDR<0.05): %d\n",
            nrow(bcell_res), sum(bcell_res$fisher_fdr < 0.05, na.rm = TRUE)))

report_block("B_cell_all: expected INDUCED near FDC",
             markers$induced_expected, bcell_res)
report_block("B_cell_all: expected REPRESSED near FDC (EBI2)",
             markers$repressed_expected, bcell_res)
report_block("B_cell_all: Cxcl13 ambient-RNA check",
             markers$ambient_check, bcell_res)

cat("\n=== Top 20 B_cell_all hits by Fisher FDR ===\n")
print(head(bcell_res[order(fisher_fdr),
           .(gene, median_coef, fisher_fdr, sign_consistency, n_samples)], 20))

cat("\nAll done.\n")
