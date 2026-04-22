# ============================================================================
# Validate RIPPLE on naive lymph node Xenium data (TDLN_remodelling project)
# ============================================================================
# Query: FRC (cell_type_assignment_L1)
# Targets: all other L1 cell types (auto-detect)
# Samples: naive only (group == "naive"; expected: m3, m5, m7, m8)
#
# Biological positive controls expected in T cells (near FRCs):
#   Induced (negative beta): CCR7, SELL, LEF1, TCF7, S1PR1
#     - these are T-zone homing / naive markers receiving CCL19/21 from FRCs
#
# Run with:
#   Rscript data-raw/xenium_tdln/run_ripple_naive_frc.R
#
# Output:
#   data-raw/xenium_tdln/naive_frc_ripple/ (RIPPLE per-celltype + summary)
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  devtools::load_all(quiet = TRUE)
})

CELLTYPE_COL <- "cell_type_assignment_L1"
SAMPLE_COL   <- "mouse_id"
QUERY_CT     <- "FRC"
POSITIVE_CTRL_GENES <- c("Ccr7", "Sell", "Lef1", "Tcf7", "S1pr1")

obj_path <- "data-raw/xenium_tdln/seurat_xenium_filtered.rds"
out_dir  <- "data-raw/xenium_tdln/naive_frc_ripple"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load and subset to naive ------------------------------------------
cat("Loading Seurat object...\n")
obj <- readRDS(obj_path)
cat("Full object:", ncol(obj), "cells x", nrow(obj), "genes\n")

md <- obj@meta.data
stopifnot(CELLTYPE_COL %in% names(md), "group" %in% names(md), SAMPLE_COL %in% names(md))

cat("\nGroup breakdown:\n"); print(table(md$group))
cat("\nMouse x group:\n"); print(table(md[[SAMPLE_COL]], md$group))

cat("\nSubsetting to naive...\n")
naive_cells <- rownames(md)[md$group == "naive"]
obj <- subset(obj, cells = naive_cells)
cat("Naive subset:", ncol(obj), "cells across",
    length(unique(obj@meta.data[[SAMPLE_COL]])), "mice\n\n")

# --- 2. QC: L1 counts per sample + marker presence ------------------------
md_n <- obj@meta.data
cat("=== L1 counts per naive mouse ===\n")
print(table(md_n[[CELLTYPE_COL]], md_n[[SAMPLE_COL]]))

cat("\n=== Query cell (FRC) availability ===\n")
frc_per_sample <- table(md_n[[SAMPLE_COL]][md_n[[CELLTYPE_COL]] == QUERY_CT])
print(frc_per_sample)

all_genes <- rownames(obj)
cat("\n=== Positive-control marker availability (case-insensitive) ===\n")
for (g in POSITIVE_CTRL_GENES) {
  hits <- grep(paste0("^", g, "$"), all_genes, ignore.case = TRUE, value = TRUE)
  cat(sprintf("  %-8s -> %s\n", g,
              if (length(hits)) paste(hits, collapse = ", ") else "MISSING"))
}

# --- 3. Run RIPPLE --------------------------------------------------------
cat("\nRunning RIPPLE...\n")
cat("  query_celltype  =", QUERY_CT, "\n")
cat("  celltype_column =", CELLTYPE_COL, "\n")
cat("  sample_column   =", SAMPLE_COL, "\n")
cat("  output_dir      =", out_dir, "\n")

t0 <- proc.time()
results <- run_ripple(
  input            = obj,
  query_celltype   = QUERY_CT,
  celltype_column  = CELLTYPE_COL,
  sample_column    = SAMPLE_COL,
  output_dir       = out_dir,
  analysis_name    = "ripple",
  k_neighbors      = 1,
  max_distance_um  = 200,
  n_permutations   = 0,
  verbose          = TRUE
)
elapsed <- (proc.time() - t0)[3]
cat(sprintf("\nRIPPLE finished in %.1f minutes\n", elapsed / 60))

# --- 4. Quick positive-control check --------------------------------------
cat("\n=== Positive-control marker results (T cell targets) ===\n")
t_targets <- c("Naive_CD4", "Naive_CD8", "Activated_CD8", "Cytotoxic_CD8",
               "Treg", "Tfh", "Tpex", "gdT_cell")
pc_hits <- results[cell_type %in% t_targets &
                   toupper(gene) %in% toupper(POSITIVE_CTRL_GENES)]
if (nrow(pc_hits)) {
  print(pc_hits[, .(gene, cell_type, median_coef, fisher_fdr, sign_consistency,
                    n_samples)][order(cell_type, gene)])
} else {
  cat("No positive-control markers were tested (likely filtered out by expression).\n")
}

cat("\n=== Top 20 most significant genes overall ===\n")
print(head(results[order(fisher_fdr),
                   .(gene, cell_type, median_coef, fisher_fdr, sign_consistency)],
           20))

cat("\nAll done.\n")
