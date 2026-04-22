# ============================================================================
# Naive LN validation v3 - TRC-Ccl21a query with collapsed T cell targets
# ============================================================================
# Rationale (vs the FDC/B cell test): the follicular/GC biology tested earlier
# is dominated by density effects (CXCR5 is saturated on follicular B cells,
# so there is no continuous per-cell gradient within the follicle) and
# requires active GCs that are rare in naive LN. The T zone around CCL21-
# producing TRCs is a better positive control for RIPPLE because CCR7/CCL21
# signaling produces a genuine continuous transcriptional gradient: naive
# T cells retained by chemokine exposure maintain high Ccr7/Sell/Lef1/Tcf7,
# while cells drifting away toward the T-B border or medulla downregulate
# these maintenance markers and begin expressing effector programs.
#
# Design (asymmetric resolution):
#   Query  = TRC-Ccl21a (fine stromal subtype, from L2 annotation CSV)
#   Target = T_cell_all (all 8 L1 T cell subtypes collapsed:
#            Activated_CD8, Cytotoxic_CD8, gdT_cell, Naive_CD4, Naive_CD8,
#            Tfh, Tpex, Treg)
#
# Expected biology in T_cell_all in naive LN:
#   Induced near TRC-Ccl21a (negative coef):
#     Ccr7, Sell, Lef1, Tcf7, Klf2, Il7r  (naive/retention program)
#   Repressed near TRC-Ccl21a (positive coef):
#     Gzma, Gzmb, Prf1, Nkg7, Cx3cr1      (effector program, T-B border /
#                                          medulla enrichment)
#   Ambient-RNA check: Ccl21a (produced by query, not by T cells)
#
# Run with:
#   Rscript data-raw/xenium_tdln/run_ripple_naive_trc.R
#
# Output:
#   data-raw/xenium_tdln/naive_trc_ripple/
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  devtools::load_all(quiet = TRUE)
})

L1_COL     <- "cell_type_assignment_L1"
FINAL_COL  <- "cell_type_assignment_L2"
SAMPLE_COL <- "mouse_id"
QUERY_CT   <- "TRC-Ccl21a"
T_CELL_L1  <- c("Activated_CD8", "Cytotoxic_CD8", "gdT_cell",
                "Naive_CD4", "Naive_CD8", "Tfh", "Tpex", "Treg")

obj_path <- "data-raw/xenium_tdln/seurat_xenium_filtered.rds"
ann_path <- "data-raw/xenium_tdln/frc_subcluster_annotation.csv"
out_dir  <- "data-raw/xenium_tdln/naive_trc_ripple"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load Seurat + annotation CSV --------------------------------------
cat("Loading Seurat object...\n")
obj <- readRDS(obj_path)
cat("Full object:", ncol(obj), "cells x", nrow(obj), "genes\n")

ann <- fread(ann_path)
setnames(ann, c("index", "frc_subcluster"))
cat("Subcluster annotation rows:", nrow(ann), "\n")

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

best_col <- names(match_counts)[which.max(match_counts)]
best_pct <- 100 * max(match_counts) / nrow(ann)
if (max(match_counts) < nrow(ann) * 0.80) {
  stop(sprintf("No Seurat column matches >=80%% of annotation rows. Max: %d/%d (%.1f%%) in '%s'",
               max(match_counts), nrow(ann), best_pct, best_col))
}
cat(sprintf("Using '%s' to merge annotation (%d/%d, %.1f%% match).\n",
            best_col, max(match_counts), nrow(ann), best_pct))

cell_keys <- candidate_cols[[best_col]]

# --- 3. Build L2 (subcluster overwrites L1) + collapse T cells ------------
ann[, frc_subcluster_norm := gsub(" ", "_", frc_subcluster)]
sub_lookup <- setNames(ann$frc_subcluster_norm, ann$index)

final <- as.character(md[[L1_COL]])
overwrite_idx <- which(cell_keys %in% ann$index)
final[overwrite_idx] <- sub_lookup[cell_keys[overwrite_idx]]

# Collapse all L1 T cell subtypes into T_cell_all
final[final %in% T_CELL_L1] <- "T_cell_all"

obj@meta.data[[FINAL_COL]] <- final

cat("\n=== Composition after L2 overwrite + T cell collapse ===\n")
print(sort(table(final), decreasing = TRUE))

# --- 4. Subset to naive ---------------------------------------------------
cat("\nSubsetting to naive...\n")
naive_cells <- rownames(obj@meta.data)[obj@meta.data$group == "naive"]
obj <- subset(obj, cells = naive_cells)
cat("Naive subset:", ncol(obj), "cells across",
    length(unique(obj@meta.data[[SAMPLE_COL]])), "mice\n")

md_n <- obj@meta.data

cat("\n=== Query (", QUERY_CT, ") counts per naive mouse ===\n", sep = "")
q_per_sample <- table(md_n[[SAMPLE_COL]][md_n[[FINAL_COL]] == QUERY_CT])
print(q_per_sample)
if (any(q_per_sample < 30)) {
  cat("\nWARNING: some samples have <30 query cells; they will be excluded.\n")
}

cat("\n=== T_cell_all counts per naive mouse ===\n")
print(table(md_n[[SAMPLE_COL]][md_n[[FINAL_COL]] == "T_cell_all"]))

# --- 5. Marker gene availability ------------------------------------------
all_genes <- rownames(obj)
markers <- list(
  induced_expected   = c("Ccr7", "Sell", "Lef1", "Tcf7", "Klf2", "Il7r"),
  repressed_expected = c("Gzma", "Gzmb", "Prf1", "Nkg7", "Cx3cr1"),
  ambient_check      = c("Ccl21a", "Ccl21b", "Ccl19")
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
cat("\nRunning RIPPLE with TRC-Ccl21a as query...\n")
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

tcell_res <- results[cell_type == "T_cell_all"]
cat(sprintf("\nT_cell_all genes tested: %d  |  significant (FDR<0.05): %d\n",
            nrow(tcell_res), sum(tcell_res$fisher_fdr < 0.05, na.rm = TRUE)))

report_block("T_cell_all: expected INDUCED near TRC-Ccl21a (naive program)",
             markers$induced_expected, tcell_res)
report_block("T_cell_all: expected REPRESSED near TRC-Ccl21a (effector program)",
             markers$repressed_expected, tcell_res)
report_block("T_cell_all: Ccl21a/b, Ccl19 ambient-RNA check",
             markers$ambient_check, tcell_res)

cat("\n=== Top 20 T_cell_all hits by Fisher FDR ===\n")
print(head(tcell_res[order(fisher_fdr),
           .(gene, median_coef, fisher_fdr, sign_consistency, n_samples)], 20))

cat("\nAll done.\n")
