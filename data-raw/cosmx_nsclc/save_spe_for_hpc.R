# Save the merged SPE as .rds for transfer to HPC
suppressPackageStartupMessages({
  library(zellkonverter)
  library(SpatialExperiment)
  library(SummarizedExperiment)
})

cat("Loading h5ad...\n")
sce <- readH5AD("data-raw/cosmx_nsclc/cosmx_nsclc.h5ad")
cd <- colData(sce)
ct <- as.character(cd$cell_type)
ct[grepl("^tumor", ct)] <- "tumor"

new_cd <- DataFrame(cell_type = ct, patient = cd$patient)
rownames(new_cd) <- colnames(sce)

spe <- SpatialExperiment(
  assays = list(counts = assay(sce, "counts")),
  colData = new_cd,
  spatialCoords = as.matrix(reducedDim(sce, "spatial"))
)
rm(sce)
gc(verbose = FALSE)

cat("Saving SPE:", ncol(spe), "cells x", nrow(spe), "genes\n")
saveRDS(spe, "data-raw/cosmx_nsclc/cosmx_nsclc_merged.rds")
cat("Saved: data-raw/cosmx_nsclc/cosmx_nsclc_merged.rds\n")
cat("Size:", round(file.info("data-raw/cosmx_nsclc/cosmx_nsclc_merged.rds")$size / 1e6), "MB\n")
