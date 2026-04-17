# Quick look at tumor type distribution across patients
library(zellkonverter)
library(SummarizedExperiment)

sce <- readH5AD("data-raw/cosmx_nsclc/cosmx_nsclc.h5ad")
cd <- as.data.frame(colData(sce))

# Cross-tabulate tumor types vs patients
tumor_cells <- cd[grepl("tumor", cd$cell_type, ignore.case = TRUE), ]
cat("Tumor subtypes x patients:\n")
print(table(tumor_cells$cell_type, tumor_cells$patient))

cat("\nAll cell types x patients:\n")
print(table(cd$cell_type, cd$patient))

cat("\nPatient metadata (histology):\n")
cat("Sample names:\n")
print(unique(cd$sample))
