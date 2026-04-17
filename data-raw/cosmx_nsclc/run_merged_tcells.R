# Test RIPPLE with merged T cell subtypes
# CD4: T CD4 naive + T CD4 memory → CD4_T (Treg stays separate)
# CD8: T CD8 naive + T CD8 memory → CD8_T

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(Matrix)
})
pkgload::load_all(".", quiet = TRUE)

cat("Loading h5ad...\n")
sce <- readH5AD("data-raw/cosmx_nsclc/cosmx_nsclc.h5ad")
cd <- colData(sce)
ct <- as.character(cd$cell_type)

# Merge tumor subtypes
ct[grepl("^tumor", ct)] <- "tumor"

# Merge T cell subtypes (keep Treg separate)
ct[ct %in% c("T CD4 naive", "T CD4 memory")] <- "CD4_T"
ct[ct %in% c("T CD8 naive", "T CD8 memory")] <- "CD8_T"

cat("\nMerged cell type distribution:\n")
print(sort(table(ct), decreasing = TRUE))

new_cd <- DataFrame(cell_type = ct, patient = cd$patient)
rownames(new_cd) <- colnames(sce)

spe <- SpatialExperiment(
  assays = list(counts = assay(sce, "counts")),
  colData = new_cd,
  spatialCoords = as.matrix(reducedDim(sce, "spatial"))
)
rm(sce)
gc(verbose = FALSE)

cat("\nRunning RIPPLE on merged T cell types...\n")
results <- run_ripple(
  input                = spe,
  query_celltype       = "tumor",
  celltype_column      = "cell_type",
  sample_column        = "patient",
  output_dir           = "data-raw/cosmx_nsclc/ripple_output_merged",
  analysis_name        = "tumor_merged_tcells",
  target_celltypes     = c("CD4_T", "CD8_T", "Treg"),
  min_cells_per_sample = 30,
  min_expr_pct         = 0.01,
  min_expr_floor       = 25,
  n_permutations       = 0,
  verbose              = TRUE
)

cat("\n=== Checkpoint genes after T cell merge ===\n")

cat("\nPDCD1 (PD-1):\n")
print(results[gene == "PDCD1", .(gene, cell_type, median_coef, fisher_fdr, sign_consistency, pct_expr)])

cat("\nCD274 (PD-L1):\n")
print(results[gene == "CD274", .(gene, cell_type, median_coef, fisher_fdr, sign_consistency, pct_expr)])

cat("\nCTLA4:\n")
print(results[gene == "CTLA4", .(gene, cell_type, median_coef, fisher_fdr, sign_consistency, pct_expr)])

cat("\nLAG3:\n")
print(results[gene == "LAG3", .(gene, cell_type, median_coef, fisher_fdr, sign_consistency, pct_expr)])

cat("\nHAVCR2 (TIM-3):\n")
print(results[gene == "HAVCR2", .(gene, cell_type, median_coef, fisher_fdr, sign_consistency, pct_expr)])

cat("\nTIGIT:\n")
print(results[gene == "TIGIT", .(gene, cell_type, median_coef, fisher_fdr, sign_consistency, pct_expr)])

cat("\nTop 15 in CD8_T:\n")
print(results[cell_type == "CD8_T"][order(fisher_fdr)][
  1:15,
  .(gene, median_coef, fisher_fdr, sign_consistency)
])

cat("\nTop 15 in CD4_T:\n")
print(results[cell_type == "CD4_T"][order(fisher_fdr)][
  1:15,
  .(gene, median_coef, fisher_fdr, sign_consistency)
])

cat("\nTop 15 in Treg:\n")
print(results[cell_type == "Treg"][order(fisher_fdr)][
  1:15,
  .(gene, median_coef, fisher_fdr, sign_consistency)
])
