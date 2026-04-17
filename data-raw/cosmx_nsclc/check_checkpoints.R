d <- data.table::fread("data-raw/cosmx_nsclc/ripple_output/tumor_ripple/summary/all_genes_results.csv")

cat("=== PDCD1 (PD-1 receptor) ===\n")
print(d[gene == "PDCD1", .(gene, cell_type, median_coef, fisher_fdr, sign_consistency, pct_expr)])

cat("\n=== CD274 (PD-L1 ligand) ===\n")
print(d[gene == "CD274", .(gene, cell_type, median_coef, fisher_fdr, sign_consistency, pct_expr)])

cat("\n=== CTLA4 ===\n")
print(d[gene == "CTLA4", .(gene, cell_type, median_coef, fisher_fdr, sign_consistency, pct_expr)])

cat("\n=== CD86 (CTLA4 ligand) ===\n")
print(d[gene == "CD86", .(gene, cell_type, median_coef, fisher_fdr, sign_consistency, pct_expr)])

cat("\n=== Were these genes even tested? ===\n")
cat("Total genes in results:", length(unique(d$gene)), "\n")
cat("PDCD1 in results:", "PDCD1" %in% d$gene, "\n")
cat("CD274 in results:", "CD274" %in% d$gene, "\n")
cat("CTLA4 in results:", "CTLA4" %in% d$gene, "\n")

cat("\n=== L-R results check ===\n")
lr <- data.table::fread("data-raw/cosmx_nsclc/ripple_output/lr_integration/all_lr_pairs.csv")
cat("CD274 as ligand:\n")
print(lr[ligand == "CD274"])
cat("PDCD1 as receptor:\n")
print(lr[receptor == "PDCD1"])
