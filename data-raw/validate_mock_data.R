# Validate that RIPPLE recovers the planted gradient in ripple_mock_data.
# Run with: Rscript data-raw/validate_mock_data.R
# This is not part of the package build; it's a one-off sanity check.

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
})

data(ripple_mock_data)
out_dir <- tempfile("ripple_mock_check_")
dir.create(out_dir)

results <- run_ripple(
  input                = ripple_mock_data,
  query_celltype       = "Tumor",
  celltype_column      = "cell_type",
  sample_column        = "sample_id",
  output_dir           = out_dir,
  min_cells_per_sample = 30,
  min_expr_pct         = 0,
  min_expr_floor       = 10,
  verbose              = FALSE
)

cat("\n=== Overall ===\n")
cat("Total gene x cell_type rows:", nrow(results), "\n")
cat("Columns:", paste(names(results), collapse = ", "), "\n")

cat("\n=== Planted INDUCED genes in T_cell ===\n")
cat("(expected: median_coef < 0, fisher_fdr small)\n")
print(results[
  grepl("^INDUCED", gene) & cell_type == "T_cell",
  .(gene, median_coef, fisher_fdr, sign_consistency)
])

cat("\n=== Planted REPRESSED genes in T_cell ===\n")
cat("(expected: median_coef > 0, fisher_fdr small)\n")
print(results[
  grepl("^REPRESSED", gene) & cell_type == "T_cell",
  .(gene, median_coef, fisher_fdr, sign_consistency)
])

cat("\n=== Background genes in T_cell ===\n")
cat("(expected: median_coef near 0, few significant)\n")
bg <- results[grepl("^BG_", gene) & cell_type == "T_cell"]
cat("  mean(median_coef) =", round(mean(bg$median_coef, na.rm = TRUE), 6), "\n")
cat("  sd(median_coef)   =", round(sd(bg$median_coef, na.rm = TRUE), 6), "\n")
cat(
  "  Significant (fisher_fdr < 0.05):",
  sum(bg$fisher_fdr < 0.05, na.rm = TRUE), "/", nrow(bg), "\n"
)

cat("\n=== Top 15 by fisher_fdr ===\n")
print(results[order(fisher_fdr)][
  1:15, .(gene, cell_type, median_coef, fisher_fdr, sign_consistency)
])

cat("\nMock data output directory:", out_dir, "\n")
