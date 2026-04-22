# ============================================================================
# Naive LN TRC -> T_cell_all k-neighbors sweep
# ============================================================================
# Sensitivity analysis of the positive-control gradient to k_neighbors.
# Re-uses the same query (TRC-Ccl21a) and collapsed target (T_cell_all)
# as run_ripple_naive_trc.R, but restricts to the single target so each
# k value finishes in ~15 min rather than ~3 hours.
#
# k = 1 was already run and lives at
#   data-raw/xenium_tdln/naive_trc_ripple/ripple/per_celltype/T_cell_all/
# (no _k suffix because k=1). This script fills in k = 3, 5, 10. The
# analysis_name suffix ("_k{n}") is added automatically by run_ripple
# when k_neighbors > 1, so each k lands in its own directory.
#
# Run with:
#   Rscript data-raw/xenium_tdln/run_ripple_naive_trc_ksweep.R
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
TARGET_CT  <- "T_cell_all"
T_CELL_L1  <- c("Activated_CD8", "Cytotoxic_CD8", "gdT_cell",
                "Naive_CD4", "Naive_CD8", "Tfh", "Tpex", "Treg")
K_VALUES   <- c(3, 5, 10)

obj_path <- "data-raw/xenium_tdln/seurat_xenium_filtered.rds"
ann_path <- "data-raw/xenium_tdln/frc_subcluster_annotation.csv"
out_dir  <- "data-raw/xenium_tdln/naive_trc_ripple"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load + apply L2 + subset naive (same as the main TRC script) ------
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
obj@meta.data[[FINAL_COL]] <- final

naive_cells <- rownames(obj@meta.data)[obj@meta.data$group == "naive"]
obj <- subset(obj, cells = naive_cells)
cat(sprintf("Naive subset: %d cells across %d mice\n",
            ncol(obj), length(unique(obj@meta.data[[SAMPLE_COL]]))))

# --- 2. Loop over k values ------------------------------------------------
for (k in K_VALUES) {
  cat(sprintf("\n=== Running RIPPLE at k_neighbors = %d ===\n", k))
  t0 <- proc.time()
  run_ripple(
    input            = obj,
    query_celltype   = QUERY_CT,
    celltype_column  = FINAL_COL,
    sample_column    = SAMPLE_COL,
    target_celltypes = TARGET_CT,
    output_dir       = out_dir,
    analysis_name    = "ripple",
    k_neighbors      = k,
    max_distance_um  = 200,
    n_permutations   = 0,
    verbose          = TRUE
  )
  cat(sprintf("  k=%d done in %.1f min\n", k, (proc.time() - t0)[3] / 60))
}

# --- 3. Assemble k-sweep comparison table ---------------------------------
positive_ctrl <- c("Ccr7", "Sell", "Lef1", "Tcf7", "Cd69", "Cxcr4",
                   "Gzma", "Il7r", "S1pr1", "Ets1", "Foxp1")

k_layouts <- list(
  "1"  = file.path(out_dir, "ripple", "per_celltype", TARGET_CT,
                   "meta_analysis_results.csv")
)
for (k in K_VALUES) {
  k_layouts[[as.character(k)]] <- file.path(
    out_dir, paste0("ripple_k", k), "per_celltype", TARGET_CT,
    "meta_analysis_results.csv"
  )
}

all_rows <- list()
for (kv in names(k_layouts)) {
  p <- k_layouts[[kv]]
  if (!file.exists(p)) {
    cat(sprintf("  MISSING: k=%s results file at %s\n", kv, p))
    next
  }
  dt <- fread(p)
  hits <- dt[tolower(gene) %in% tolower(positive_ctrl),
             .(gene, median_coef, fisher_fdr, sign_consistency,
               n_samples, pct_expr)]
  hits[, k := as.integer(kv)]
  all_rows[[kv]] <- hits
}
sweep_dt <- rbindlist(all_rows, fill = TRUE)

cat("\n=== k-sweep summary: positive-control genes in T_cell_all ===\n")
print(sweep_dt[order(gene, k)])

# Wide-format view: rows = gene, cols = k
library(data.table)
wide_fdr <- dcast(sweep_dt, gene ~ k, value.var = "fisher_fdr")
wide_coef <- dcast(sweep_dt, gene ~ k, value.var = "median_coef")
wide_sign <- dcast(sweep_dt, gene ~ k, value.var = "sign_consistency")

cat("\n=== Fisher FDR across k ===\n"); print(wide_fdr)
cat("\n=== median_coef across k ===\n"); print(wide_coef)
cat("\n=== sign_consistency across k ===\n"); print(wide_sign)

out_path <- file.path(out_dir, "ksweep_T_cell_all_results.rds")
saveRDS(list(per_run = sweep_dt,
             fdr_wide = wide_fdr,
             coef_wide = wide_coef,
             sign_wide = wide_sign),
        file = out_path)
cat(sprintf("\nSaved k-sweep summary: %s\n", out_path))
cat("\nAll done.\n")
