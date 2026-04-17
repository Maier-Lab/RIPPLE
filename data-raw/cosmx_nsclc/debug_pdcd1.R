# Debug: why doesn't PDCD1 show a gradient in CD8 T cells?
library(data.table)

# Load per-sample coefficients
coef_file <- "data-raw/cosmx_nsclc/ripple_output_merged/tumor_merged_tcells/per_celltype/CD8_T/coef_per_sample.csv"
coef_data <- fread(coef_file)

cat("=== PDCD1 per-sample coefficients ===\n")
print(coef_data[gene == "PDCD1"])

cat("\n=== Sign direction per patient ===\n")
pdcd1 <- coef_data[gene == "PDCD1"]
cat("Positive coef (repressed near tumor):", sum(pdcd1$coef > 0, na.rm = TRUE), "\n")
cat("Negative coef (induced near tumor):", sum(pdcd1$coef < 0, na.rm = TRUE), "\n")
cat("NA:", sum(is.na(pdcd1$coef)), "\n")

cat("\n=== For comparison: top gradient gene CXCR4 ===\n")
print(coef_data[gene == "CXCR4"])

cat("\n=== CTLA4 per-sample ===\n")
print(coef_data[gene == "CTLA4"])

cat("\n=== LAG3 per-sample ===\n")
print(coef_data[gene == "LAG3"])

cat("\n=== TIGIT per-sample ===\n")
print(coef_data[gene == "TIGIT"])

# Also check expression levels
cat("\n=== Expression percentages from meta results ===\n")
meta_file <- "data-raw/cosmx_nsclc/ripple_output_merged/tumor_merged_tcells/per_celltype/CD8_T/meta_analysis_results.csv"
meta <- fread(meta_file)
checkpoints <- c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "CXCR4", "CD274")
print(meta[gene %in% checkpoints, .(gene, combined_coef, pval, fisher_fdr, pct_expr, sign_consistency)])
