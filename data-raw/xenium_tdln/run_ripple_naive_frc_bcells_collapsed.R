# ============================================================================
# Naive LN validation - collapsed B cell variant (CXCR5/FDC positive control)
# ============================================================================
# Variant of run_ripple_naive_frc.R where all B cell subtypes in L1
# (B_cell, Follicular_B_cell, GC_B_cell) are collapsed into a single
# "B_cell_all" category. This avoids the circularity of using an annotation
# that was itself derived (in part) from the genes we want to test as
# positive controls.
#
# Biology under test:
#   FRC includes / subsumes FDCs, which produce Cxcl13.
#   Cxcr5+ B cells concentrate near FDCs in the B cell follicle.
#   At the pooled B cell level, mean Cxcr5 per cell should rise with
#   proximity to FRC. Secondary follicular/GC program genes
#   (Bcl6, Fcer2a, Cr2) should behave similarly. Gpr183 (EBI2) opposes
#   Cxcr5 (drives cells out of follicles) and should show the opposite
#   direction.
#
# Run with:
#   Rscript data-raw/xenium_tdln/run_ripple_naive_frc_bcells_collapsed.R
#
# Output:
#   data-raw/xenium_tdln/naive_frc_bcells_collapsed/
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  devtools::load_all(quiet = TRUE)
})

ORIGINAL_CELLTYPE_COL <- "cell_type_assignment_L1"
COLLAPSED_COL         <- "cell_type_L1_bcells_collapsed"
SAMPLE_COL            <- "mouse_id"
QUERY_CT              <- "FRC"

# B cell subtypes in L1 that should be merged
B_CELL_SUBTYPES <- c("B_cell", "Follicular_B_cell", "GC_B_cell")

# Genes to inspect post-hoc for the FDC -> B cell follicle positive control.
# Mouse gene symbols (Xenium panel uses Mouse 5K Pan-Tissue + lymph node add-on).
FOLLICULAR_MARKERS <- list(
  induced_expected   = c("Cxcr5", "Bcl6", "Fcer2a", "Cr2", "Cd83", "Bach2"),
  repressed_expected = c("Gpr183"),  # EBI2, opposes Cxcr5 in follicle positioning
  ambient_check      = c("Cxcl13")    # produced by FRC/FDC, not B cell
                                      # -> shouldn't show per-B-cell gradient
)

obj_path <- "data-raw/xenium_tdln/seurat_xenium_filtered.rds"
out_dir  <- "data-raw/xenium_tdln/naive_frc_bcells_collapsed"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load and subset to naive ------------------------------------------
cat("Loading Seurat object...\n")
obj <- readRDS(obj_path)
cat("Full object:", ncol(obj), "cells x", nrow(obj), "genes\n")

md <- obj@meta.data
stopifnot(ORIGINAL_CELLTYPE_COL %in% names(md),
          "group" %in% names(md),
          SAMPLE_COL %in% names(md))

cat("\nSubsetting to naive samples...\n")
naive_cells <- rownames(md)[md$group == "naive"]
obj <- subset(obj, cells = naive_cells)
cat("Naive subset:", ncol(obj), "cells across",
    length(unique(obj@meta.data[[SAMPLE_COL]])), "mice\n\n")

# --- 2. Build collapsed cell type column ----------------------------------
md_n <- obj@meta.data
original <- as.character(md_n[[ORIGINAL_CELLTYPE_COL]])

cat("Collapsing B cell subtypes:", paste(B_CELL_SUBTYPES, collapse = ", "),
    "-> B_cell_all\n")

collapsed <- original
collapsed[original %in% B_CELL_SUBTYPES] <- "B_cell_all"

obj@meta.data[[COLLAPSED_COL]] <- collapsed
md_n <- obj@meta.data

cat("\n=== Before collapse (L1) ===\n")
b_original_counts <- table(original[original %in% B_CELL_SUBTYPES])
print(b_original_counts)

cat("\n=== After collapse ===\n")
cat("B_cell_all total:", sum(collapsed == "B_cell_all"), "cells\n")
cat("B_cell_all per mouse:\n")
print(table(md_n[[SAMPLE_COL]][collapsed == "B_cell_all"]))

cat("\n=== Query (FRC) counts per mouse ===\n")
print(table(md_n[[SAMPLE_COL]][collapsed == "FRC"]))

# --- 3. Marker availability check -----------------------------------------
all_genes <- rownames(obj)
cat("\n=== Marker gene presence ===\n")
report_marker <- function(label, genes) {
  cat(sprintf("[%s]\n", label))
  for (g in genes) {
    hits <- grep(paste0("^", g, "$"), all_genes,
                 ignore.case = TRUE, value = TRUE)
    cat(sprintf("  %-10s -> %s\n", g,
                if (length(hits)) paste(hits, collapse = ", ") else "MISSING"))
  }
}
report_marker("Expected induced near FRC/FDC", FOLLICULAR_MARKERS$induced_expected)
report_marker("Expected repressed near FRC/FDC", FOLLICULAR_MARKERS$repressed_expected)
report_marker("Ambient-RNA check (produced by FRC)", FOLLICULAR_MARKERS$ambient_check)

# --- 4. Run RIPPLE --------------------------------------------------------
cat("\nRunning RIPPLE with collapsed B cell labels...\n")
t0 <- proc.time()
results <- run_ripple(
  input            = obj,
  query_celltype   = QUERY_CT,
  celltype_column  = COLLAPSED_COL,
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

# --- 5. Positive-control readouts -----------------------------------------
cat("\n=== B_cell_all results for expected INDUCED markers ===\n")
ind_markers <- FOLLICULAR_MARKERS$induced_expected
ind_hits <- results[cell_type == "B_cell_all" &
                    toupper(gene) %in% toupper(ind_markers)]
if (nrow(ind_hits)) {
  print(ind_hits[order(fisher_fdr),
                 .(gene, median_coef, fisher_fdr, sign_consistency,
                   n_samples, pct_expr)])
} else {
  cat("No expected-induced markers were tested (likely expression-filtered).\n")
}

cat("\n=== B_cell_all results for expected REPRESSED markers ===\n")
rep_markers <- FOLLICULAR_MARKERS$repressed_expected
rep_hits <- results[cell_type == "B_cell_all" &
                    toupper(gene) %in% toupper(rep_markers)]
if (nrow(rep_hits)) {
  print(rep_hits[order(fisher_fdr),
                 .(gene, median_coef, fisher_fdr, sign_consistency,
                   n_samples, pct_expr)])
} else {
  cat("No expected-repressed markers were tested.\n")
}

cat("\n=== B_cell_all Cxcl13 (ambient check: should be weak/null) ===\n")
ambient_hits <- results[cell_type == "B_cell_all" &
                        toupper(gene) %in% toupper(FOLLICULAR_MARKERS$ambient_check)]
if (nrow(ambient_hits)) {
  print(ambient_hits[, .(gene, median_coef, fisher_fdr, sign_consistency,
                         n_samples, pct_expr)])
} else {
  cat("Cxcl13 was not tested in B_cell_all (likely expression-filtered,",
      "consistent with it being primarily a stromal gene).\n")
}

cat("\n=== Top 20 B_cell_all hits by Fisher FDR ===\n")
print(head(results[cell_type == "B_cell_all"][order(fisher_fdr),
           .(gene, median_coef, fisher_fdr, sign_consistency, n_samples)], 20))

cat("\nAll done.\n")
