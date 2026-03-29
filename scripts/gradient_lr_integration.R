#!/usr/bin/env Rscript
#' =============================================================================
#' RIPPLE Stage 6: Gradient-to-Ligand-Receptor Integration Analysis
#' =============================================================================
#'
#' Connects distance correlation gradient results to specific L-R mechanisms.
#' Answers: "Which L-R pairs between query and target cells explain the
#' observed distance-dependent expression gradients?"
#'
#' Three complementary approaches (Direction A: Query -> Target):
#'   Part 1: Direct L-R mapping -- gradient genes that are known receptors,
#'           matched to query-expressed ligands (+ lightweight reverse for B)
#'   Part 2: NicheNet ligand activity -- which query ligands best predict each
#'           cell type's gradient pattern
#'   Part 3: Downstream target enrichment -- Fisher's exact test validating that
#'           predicted L-R downstream targets overlap with gradient genes
#'
#' Runs as SLURM array job (CELLTYPE_INDEX 1-14).
#'
#' Usage:
#'   Rscript gradient_lr_integration.R
#'   QUERY_CELLTYPE=MyType CELLTYPE_COLUMN=my_col Rscript gradient_lr_integration.R
#'
#' Author: CMM Project
#' =============================================================================

# =============================================================================
# Setup
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(nichenetr)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(scales)
})

set.seed(42)

# Source utilities
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  return(getwd())
}
script_dir <- get_script_dir()
source(file.path(script_dir, "utils.R"))

# =============================================================================
# Configuration
# =============================================================================

ANALYSIS_NAME <- "gradient_lr_integration"
organism <- "mouse"

# Gradient source: which distance correlation results to use
# Default: "hymy_distance_correlation" (v1, logistic)
# For v2 (Poisson GLM): GRADIENT_SOURCE=hymy_distance_correlation_v2
GRADIENT_SOURCE <- Sys.getenv("GRADIENT_SOURCE", unset = "hymy_distance_correlation")

# Derive output suffix from gradient source (e.g., "_v2" for v2 results)
gradient_suffix <- sub("^hymy_distance_correlation", "", GRADIENT_SOURCE)
output_name <- if (nchar(gradient_suffix) > 0) {
  paste0(ANALYSIS_NAME, gradient_suffix)
} else {
  ANALYSIS_NAME
}

# Inherited from config.R (via utils.R): QUERY_CELLTYPE, CELLTYPE_COL, OUTPUT_SUFFIX, QUERY_LABEL
GRADIENT_DIR <- file.path(OUTPUT_ROOT, GRADIENT_SOURCE)
STAGE2_DIR <- file.path(OUTPUT_ROOT, paste0(GRADIENT_SOURCE, "_stage2"))
OUTPUT_BASE <- file.path(OUTPUT_ROOT, output_name)

ensure_dir(OUTPUT_BASE)
ensure_dir(file.path(OUTPUT_BASE, "per_celltype"))
ensure_dir(file.path(OUTPUT_BASE, "summary"))
ensure_dir(file.path(OUTPUT_BASE, "plots"))

# Expression threshold for query cell ligand/receptor filtering
EXPR_THRESHOLD_PCT <- 5  # At least 5% of query cells must express

# NicheNet database cache
NICHENET_CACHE <- file.path(PROJECT_ROOT, "resources", "nichenet")
ensure_dir(NICHENET_CACHE)

# Target cell types: auto-detect from per_celltype directories, or from env var
TARGET_CELLTYPES_ENV <- Sys.getenv("TARGET_CELLTYPES", unset = "")
if (nchar(TARGET_CELLTYPES_ENV) > 0) {
  target_names <- trimws(strsplit(TARGET_CELLTYPES_ENV, ",")[[1]])
  TARGET_CELLTYPES <- as.list(setNames(target_names, target_names))
} else {
  # Auto-detect from gradient results per_celltype directories
  ct_base <- file.path(GRADIENT_DIR, "per_celltype")
  if (dir.exists(ct_base)) {
    target_names <- basename(list.dirs(ct_base, recursive = FALSE))
    message("Auto-detected ", length(target_names), " target cell types from: ", ct_base)
    TARGET_CELLTYPES <- as.list(setNames(target_names, target_names))
  } else {
    # Legacy fallback (backward compatible)
    TARGET_CELLTYPES <- list(
      LEC = "LEC",
      FRC = "FRC",
      BEC = "BEC",
      CD4_T_cells = c("Naive_CD4", "Tfh", "Treg"),
      CD8_T_cells = c("Naive_CD8", "Activated_CD8", "Cytotoxic_CD8", "Tpex"),
      gdT_cells = "gdT_cell",
      Macrophages = "Macrophages",
      Monocyte = "Monocyte",
      Fibroblasts_mac = "Fibroblasts_mac",
      cDC1 = "cDC1",
      cDC2 = "cDC2",
      mature_migDC = "mature_migDC",
      B_cells = c("B_cell", "Follicular_B"),
      Plasma_cell = "Plasma_cell"
    )
  }
}

# =============================================================================
# Cell Type Selection (SLURM array job)
# =============================================================================

CELLTYPE_INDEX <- as.integer(Sys.getenv("CELLTYPE_INDEX", unset = "0"))

if (CELLTYPE_INDEX > 0) {
  celltype_names <- names(TARGET_CELLTYPES)
  if (CELLTYPE_INDEX > length(celltype_names)) {
    stop("CELLTYPE_INDEX ", CELLTYPE_INDEX, " exceeds number of cell types (",
         length(celltype_names), ")")
  }
  selected_ct <- celltype_names[CELLTYPE_INDEX]
  message("\n>>> Array job mode: Processing only ", selected_ct,
          " (index ", CELLTYPE_INDEX, ")")
  TARGET_CELLTYPES <- TARGET_CELLTYPES[selected_ct]
}

message(strrep("=", 70))
message("Gradient-to-Ligand-Receptor Integration Analysis")
message(strrep("=", 70))
message("Annotation level: ", ANNOTATION_LEVEL)
message("Gradient source:  ", GRADIENT_SOURCE)
message("Query cell type:  ", QUERY_CELLTYPE)
message("Output directory: ", OUTPUT_BASE)

# =============================================================================
# Load Gradient Results
# =============================================================================

message("\n", strrep("-", 70))
message("Loading gradient results...")
message(strrep("-", 70), "\n")

# Stage 1: all genes results
gradient_file <- file.path(GRADIENT_DIR, "summary", "all_genes_results.csv")
if (!file.exists(gradient_file)) {
  stop("Stage 1 gradient results not found: ", gradient_file,
       "\nRun hymy_distance_correlation.R and merge_distance_correlation_results.R first.")
}
all_gradient <- fread(gradient_file)
message(sprintf("  Stage 1: %d gene-celltype rows loaded", nrow(all_gradient)))

# Select coefficient and FDR columns: prefer Fisher/median if available (from recompute_meta_summary)
COEF_COL_GRAD <- if ("median_coef" %in% names(all_gradient)) "median_coef" else "combined_coef"
FDR_COL_GRAD <- if ("fisher_fdr" %in% names(all_gradient)) "fisher_fdr" else "fdr"
message(sprintf("  Using coefficient column: %s", COEF_COL_GRAD))
message(sprintf("  Using FDR column: %s", FDR_COL_GRAD))

# Stage 2: classification (optional but recommended)
stage2_file <- file.path(STAGE2_DIR, "summary", "stage2_all_results.csv")
has_stage2 <- file.exists(stage2_file)
if (has_stage2) {
  stage2_results <- fread(stage2_file)
  message(sprintf("  Stage 2: %d classified rows loaded", nrow(stage2_results)))
} else {
  message("  Stage 2: Not found (classification will be NA)")
  stage2_results <- NULL
}

# =============================================================================
# Load NicheNet Databases (with caching)
# =============================================================================

message("\n", strrep("-", 70))
message("Loading NicheNet databases...")
message(strrep("-", 70), "\n")

options(timeout = 300)

# L-R network
lr_cache_path <- file.path(NICHENET_CACHE, "lr_network_mouse_allInfo_30112033.rds")
if (file.exists(lr_cache_path)) {
  message("  Loading L-R network from cache...")
  lr_network_all <- readRDS(lr_cache_path)
} else {
  message("  Downloading L-R network from Zenodo...")
  lr_network_all <- readRDS(url(
    "https://zenodo.org/record/10229222/files/lr_network_mouse_allInfo_30112033.rds"
  ))
  saveRDS(lr_network_all, lr_cache_path)
  message("  Cached to: ", lr_cache_path)
}

# Process gene names (same as run_differential_nichenet_spatial.R:219-224)
lr_network_all <- lr_network_all %>%
  mutate(
    ligand = convert_alias_to_symbols(ligand, organism = organism),
    receptor = convert_alias_to_symbols(receptor, organism = organism)
  ) %>%
  mutate(ligand = make.names(ligand), receptor = make.names(receptor))

lr_network <- as.data.table(lr_network_all %>% distinct(ligand, receptor))

# Remove erroneous CSF1-CSF3R pair (CLAUDE.md convention)
lr_network <- lr_network[!(ligand == "Csf1" & receptor == "Csf3r")]
message(sprintf("  L-R network: %d unique pairs (after removing CSF1-CSF3R)",
                nrow(lr_network)))

# Ligand-target matrix
lt_cache_path <- file.path(NICHENET_CACHE, "ligand_target_matrix_nsga2r_final_mouse.rds")
if (file.exists(lt_cache_path)) {
  message("  Loading ligand-target matrix from cache...")
  ligand_target_matrix <- readRDS(lt_cache_path)
} else {
  message("  Downloading ligand-target matrix from Zenodo...")
  ligand_target_matrix <- readRDS(url(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"
  ))
  saveRDS(ligand_target_matrix, lt_cache_path)
  message("  Cached to: ", lt_cache_path)
}

# Make names syntactically valid
colnames(ligand_target_matrix) <- colnames(ligand_target_matrix) %>%
  convert_alias_to_symbols(organism = organism) %>% make.names()
rownames(ligand_target_matrix) <- rownames(ligand_target_matrix) %>%
  convert_alias_to_symbols(organism = organism) %>% make.names()

# Filter L-R network to ligands present in target matrix
lr_network <- lr_network[ligand %in% colnames(ligand_target_matrix)]
ligand_target_matrix <- ligand_target_matrix[, unique(lr_network$ligand)]

message(sprintf("  Filtered to %d L-R pairs with target predictions", nrow(lr_network)))
message(sprintf("  Ligand-target matrix: %d targets x %d ligands",
                nrow(ligand_target_matrix), ncol(ligand_target_matrix)))

# =============================================================================
# Load Seurat Object & Compute Query Cell Expression Profile (Per-Sample + Pooled)
# =============================================================================

message("\n", strrep("-", 70))
message(paste0("Loading Seurat object and computing ", QUERY_LABEL, " expression profiles..."))
message(strrep("-", 70), "\n")

obj <- load_seurat()
if (ANNOTATION_LEVEL == "HyMy") {
  obj <- merge_hymy_annotations(obj)
}

meta <- as.data.table(obj@meta.data, keep.rownames = "barcode")

# Condition column resolution (standard pattern)
if (nchar(CONDITION_COL) > 0 && CONDITION_COL %in% names(meta)) {
  meta[, condition := get(CONDITION_COL)]
} else if ("condition" %in% names(meta)) {
  # already exists
} else if ("group" %in% names(meta)) {
  meta[, condition := group]
} else {
  meta[, condition := "all"]
}

# Get query cells (with optional condition filtering)
if (nchar(CONDITION_VAL) > 0) {
  query_barcodes <- meta[get(CELLTYPE_COL) == QUERY_CELLTYPE & condition == CONDITION_VAL]$barcode
  message(sprintf("  Query cells (%s, %s): %d", QUERY_CELLTYPE, CONDITION_VAL, length(query_barcodes)))
} else {
  query_barcodes <- meta[get(CELLTYPE_COL) == QUERY_CELLTYPE]$barcode
  message(sprintf("  Query cells (%s, all conditions): %d", QUERY_CELLTYPE, length(query_barcodes)))
}

if (length(query_barcodes) == 0) {
  stop("No query cells found! Check CELLTYPE_COL and QUERY_CELLTYPE.")
}

# Compute query cell expression statistics (sparse -- don't as.matrix the whole thing)
expr_matrix <- GetAssayData(obj, layer = "data")
rownames(expr_matrix) <- make.names(rownames(expr_matrix))

# --- Per-sample query cell expression profiles ---
if (nchar(CONDITION_VAL) > 0) {
  analysis_samples <- unique(meta[condition == CONDITION_VAL][[SAMPLE_COL]])
} else {
  analysis_samples <- unique(meta[[SAMPLE_COL]])
}
message(sprintf("  Analysis samples: %s", paste(analysis_samples, collapse = ", ")))

hymy_per_sample <- rbindlist(lapply(analysis_samples, function(sid) {
  barcodes <- meta[get(CELLTYPE_COL) == QUERY_CELLTYPE & get(SAMPLE_COL) == sid]$barcode
  if (length(barcodes) == 0) return(NULL)
  sample_expr <- expr_matrix[, barcodes, drop = FALSE]
  data.table(
    gene = rownames(expr_matrix),
    pct_hymy = as.numeric(rowMeans(sample_expr > 0) * 100),
    mean_hymy = as.numeric(rowMeans(sample_expr)),
    sample_id = sid,
    n_hymy_cells = length(barcodes)
  )
}))

message(sprintf("  Per-sample %s profiles: %d samples x %d genes", QUERY_LABEL,
                length(unique(hymy_per_sample$sample_id)),
                length(unique(hymy_per_sample$gene))))

# --- Pooled query cell expression (used for NicheNet ligand identification) ---
query_expr <- expr_matrix[, query_barcodes, drop = FALSE]
hymy_pct <- rowMeans(query_expr > 0) * 100
hymy_mean <- rowMeans(query_expr)

hymy_profile <- data.table(
  gene = names(hymy_pct),
  pct_hymy = as.numeric(hymy_pct),
  mean_hymy = as.numeric(hymy_mean)
)

# Query-expressed ligands and receptors (pooled -- for NicheNet input)
available_genes <- rownames(expr_matrix)
available_ligands <- intersect(lr_network$ligand, available_genes)
available_receptors <- intersect(lr_network$receptor, available_genes)

hymy_expressed_ligands <- hymy_profile[gene %in% available_ligands &
                                        pct_hymy >= EXPR_THRESHOLD_PCT]$gene
hymy_expressed_receptors <- hymy_profile[gene %in% available_receptors &
                                          pct_hymy >= EXPR_THRESHOLD_PCT]$gene

message(sprintf("  %s-expressed ligands (pooled, >=%d%%): %d",
                QUERY_LABEL, EXPR_THRESHOLD_PCT, length(hymy_expressed_ligands)))
message(sprintf("  %s-expressed receptors (pooled, >=%d%%): %d",
                QUERY_LABEL, EXPR_THRESHOLD_PCT, length(hymy_expressed_receptors)))

# =============================================================================
# Per Cell Type Analysis
# =============================================================================

for (ct_name in names(TARGET_CELLTYPES)) {

  message("\n", strrep("=", 70))
  message(sprintf("Processing: %s", ct_name))
  message(strrep("=", 70))

  ct_output_dir <- file.path(OUTPUT_BASE, "per_celltype", ct_name)
  ensure_dir(ct_output_dir)

  # -------------------------------------------------------------------------
  # Get gradient data for this cell type
  # -------------------------------------------------------------------------
  ct_gradient <- all_gradient[cell_type == ct_name]
  if (nrow(ct_gradient) == 0) {
    message("  No gradient data for ", ct_name, ". Skipping.")
    next
  }

  # Make gene names syntactically valid for NicheNet matching
  ct_gradient[, gene_safe := make.names(gene)]

  # Merge Stage 2 classification BEFORE extracting sig_genes
  if (has_stage2) {
    s2_ct <- stage2_results[cell_type == ct_name, .(gene, classification)]
    s2_ct[, gene_safe := make.names(gene)]
    ct_gradient <- merge(ct_gradient, s2_ct[, .(gene_safe, classification)],
                         by = "gene_safe", all.x = TRUE)
  } else {
    ct_gradient[, classification := NA_character_]
  }

  # Significant gradient genes (now includes classification column)
  sig_genes <- ct_gradient[get(FDR_COL_GRAD) < 0.05]
  induced_genes <- sig_genes[get(COEF_COL_GRAD) < 0]$gene_safe   # Higher near query
  repressed_genes <- sig_genes[get(COEF_COL_GRAD) > 0]$gene_safe  # Lower near query
  background_genes <- ct_gradient$gene_safe

  message(sprintf("  Total tested genes: %d", nrow(ct_gradient)))
  message(sprintf("  Significant (FDR<0.05): %d (%d induced, %d repressed)",
                  nrow(sig_genes), length(induced_genes), length(repressed_genes)))

  # =========================================================================
  # Part 1: Direct L-R Mapping (Per-Mouse Matched)
  # =========================================================================

  message("\n  --- Part 1: Direct L-R Mapping (per-mouse matched) ---")

  # Load per-sample gradient coefficients for this cell type
  coef_per_sample_file <- file.path(GRADIENT_DIR, "per_celltype", ct_name, "coef_per_sample.csv")
  if (file.exists(coef_per_sample_file)) {
    coef_per_sample <- fread(coef_per_sample_file)
    coef_per_sample[, gene_safe := make.names(gene)]
    message(sprintf("  Per-sample gradient coefficients: %d gene-sample rows loaded",
                    nrow(coef_per_sample)))
  } else {
    message("  WARNING: Per-sample coefficients not found at: ", coef_per_sample_file)
    message("  Falling back to pooled-only scoring.")
    coef_per_sample <- NULL
  }

  # ---- Direction A: Query ligand -> Target receptor ----

  # Find significant gradient genes that are known receptors
  sig_receptors <- sig_genes[gene_safe %in% available_receptors]
  message(sprintf("  Direction A: %d gradient genes are known receptors",
                  nrow(sig_receptors)))

  direct_a_results <- data.table()

  if (nrow(sig_receptors) > 0) {
    for (i in seq_len(nrow(sig_receptors))) {
      receptor_gene <- sig_receptors$gene_safe[i]
      receptor_info <- sig_receptors[i]

      # Find cognate ligands expressed by query cells (using pooled threshold for discovery)
      cognate_ligands <- lr_network[receptor == receptor_gene &
                                      ligand %in% hymy_expressed_ligands]$ligand

      if (length(cognate_ligands) > 0) {
        for (lig in cognate_ligands) {
          lig_info_pooled <- hymy_profile[gene == lig]

          # --- Per-sample matching ---
          per_sample_scores <- numeric(0)
          per_sample_details <- character(0)
          n_supporting <- 0L

          if (!is.null(coef_per_sample)) {
            for (sid in analysis_samples) {
              # Per-sample ligand expression on query cells
              lig_sample <- hymy_per_sample[gene == lig & sample_id == sid]
              # Per-sample receptor gradient coefficient
              rec_sample <- coef_per_sample[gene_safe == receptor_gene & sample_id == sid]

              if (nrow(lig_sample) > 0 && nrow(rec_sample) > 0) {
                lig_pct <- lig_sample$pct_hymy[1]
                rec_coef <- rec_sample$coef[1]

                # Check that per-sample coef direction agrees with meta-level direction
                expected_negative <- (receptor_info[[COEF_COL_GRAD]] < 0)
                direction_agrees <- if (expected_negative) rec_coef < 0 else rec_coef > 0

                # Pair is "supported" if ligand expressed AND gradient present AND direction agrees
                if (lig_pct >= EXPR_THRESHOLD_PCT && !is.na(rec_coef) && direction_agrees) {
                  sample_score <- abs(rec_coef) * (lig_pct / 100)
                  per_sample_scores <- c(per_sample_scores, sample_score)
                  per_sample_details <- c(per_sample_details,
                    sprintf("%s:%.4f", sid, sample_score))
                  n_supporting <- n_supporting + 1L
                }
              }
            }
          }

          n_total_samples <- length(analysis_samples)
          reproducibility <- n_supporting / n_total_samples

          # Compute per-mouse-matched direct score
          if (length(per_sample_scores) >= 2) {
            median_score <- median(per_sample_scores)
          } else if (length(per_sample_scores) == 1) {
            median_score <- per_sample_scores[1]
          } else {
            # Fallback: pooled score (coef_per_sample not available)
            median_score <- abs(receptor_info[[COEF_COL_GRAD]]) *
                            (lig_info_pooled$pct_hymy / 100)
          }

          direct_a_results <- rbind(direct_a_results, data.table(
            ligand = lig,
            receptor = receptor_gene,
            cell_type = ct_name,
            direction = paste0(QUERY_LABEL, "_to_Target"),
            ligand_pct_hymy = lig_info_pooled$pct_hymy,
            ligand_mean_hymy = lig_info_pooled$mean_hymy,
            receptor_gradient_coef = receptor_info[[COEF_COL_GRAD]],
            receptor_fdr = receptor_info[[FDR_COL_GRAD]],
            coef_column_used = COEF_COL_GRAD,
            fdr_column_used = FDR_COL_GRAD,
            receptor_pct_target = receptor_info$pct_expr,
            receptor_mean_target = receptor_info$mean_expr,
            receptor_decay_pattern = receptor_info$decay_pattern,
            stage2_classification = receptor_info$classification,
            n_samples_supporting = n_supporting,
            n_samples_total = n_total_samples,
            reproducibility = reproducibility,
            median_per_sample_score = median_score,
            per_sample_scores = paste(per_sample_details, collapse = ";")
          ), fill = TRUE)
        }
      }
    }
  }

  # Score Direction A pairs: per-mouse median × reproducibility
  # (all receptors already passed FDR < 0.05 gate — no FDR multiplier to avoid double-dipping)
  if (nrow(direct_a_results) > 0) {
    direct_a_results[, direct_score :=
      median_per_sample_score * reproducibility]
    setorder(direct_a_results, -direct_score)
    n_reprod <- sum(direct_a_results$n_samples_supporting >= 3)
    message(sprintf("  Direction A: %d L-R pairs found (%d reproduced in ≥3/4 samples)",
                    nrow(direct_a_results), n_reprod))
  } else {
    message("  Direction A: 0 L-R pairs found")
  }

  fwrite(direct_a_results, file.path(ct_output_dir, "direct_lr_pairs_hymy_to_target.csv"))

  # ---- Direction B: Target ligand -> Query receptor (per-mouse matched) ----

  sig_ligands <- sig_genes[gene_safe %in% available_ligands]
  message(sprintf("  Direction B: %d gradient genes are known ligands",
                  nrow(sig_ligands)))

  direct_b_results <- data.table()

  if (nrow(sig_ligands) > 0) {
    for (i in seq_len(nrow(sig_ligands))) {
      ligand_gene <- sig_ligands$gene_safe[i]
      ligand_info <- sig_ligands[i]

      # Find cognate receptors expressed by query cells
      cognate_receptors <- lr_network[ligand == ligand_gene &
                                        receptor %in% hymy_expressed_receptors]$receptor

      if (length(cognate_receptors) > 0) {
        for (rec in cognate_receptors) {
          rec_info_pooled <- hymy_profile[gene == rec]

          # --- Per-sample matching (Direction B) ---
          per_sample_scores_b <- numeric(0)
          per_sample_details_b <- character(0)
          n_supporting_b <- 0L

          if (!is.null(coef_per_sample)) {
            for (sid in analysis_samples) {
              # Per-sample receptor expression on query cells
              rec_sample <- hymy_per_sample[gene == rec & sample_id == sid]
              # Per-sample ligand gradient coefficient
              lig_sample <- coef_per_sample[gene_safe == ligand_gene & sample_id == sid]

              if (nrow(rec_sample) > 0 && nrow(lig_sample) > 0) {
                rec_pct <- rec_sample$pct_hymy[1]
                lig_coef <- lig_sample$coef[1]

                # Check that per-sample coef direction agrees with meta-level direction
                expected_negative_b <- (ligand_info[[COEF_COL_GRAD]] < 0)
                direction_agrees_b <- if (expected_negative_b) lig_coef < 0 else lig_coef > 0

                if (rec_pct >= EXPR_THRESHOLD_PCT && !is.na(lig_coef) && direction_agrees_b) {
                  sample_score_b <- abs(lig_coef) * (rec_pct / 100)
                  per_sample_scores_b <- c(per_sample_scores_b, sample_score_b)
                  per_sample_details_b <- c(per_sample_details_b,
                    sprintf("%s:%.4f", sid, sample_score_b))
                  n_supporting_b <- n_supporting_b + 1L
                }
              }
            }
          }

          n_total_b <- length(analysis_samples)
          reproducibility_b <- n_supporting_b / n_total_b

          if (length(per_sample_scores_b) >= 2) {
            median_score_b <- median(per_sample_scores_b)
          } else if (length(per_sample_scores_b) == 1) {
            median_score_b <- per_sample_scores_b[1]
          } else {
            median_score_b <- abs(ligand_info[[COEF_COL_GRAD]]) *
                              (rec_info_pooled$pct_hymy / 100)
          }

          direct_b_results <- rbind(direct_b_results, data.table(
            ligand = ligand_gene,
            receptor = rec,
            cell_type = ct_name,
            direction = paste0("Target_to_", QUERY_LABEL),
            ligand_gradient_coef = ligand_info[[COEF_COL_GRAD]],
            ligand_fdr = ligand_info[[FDR_COL_GRAD]],
            coef_column_used = COEF_COL_GRAD,
            fdr_column_used = FDR_COL_GRAD,
            ligand_pct_target = ligand_info$pct_expr,
            ligand_mean_target = ligand_info$mean_expr,
            ligand_decay_pattern = ligand_info$decay_pattern,
            receptor_pct_hymy = rec_info_pooled$pct_hymy,
            receptor_mean_hymy = rec_info_pooled$mean_hymy,
            stage2_classification = ligand_info$classification,
            n_samples_supporting = n_supporting_b,
            n_samples_total = n_total_b,
            reproducibility = reproducibility_b,
            median_per_sample_score = median_score_b,
            per_sample_scores = paste(per_sample_details_b, collapse = ";")
          ), fill = TRUE)
        }
      }
    }
  }

  if (nrow(direct_b_results) > 0) {
    direct_b_results[, direct_score :=
      median_per_sample_score * reproducibility]
    setorder(direct_b_results, -direct_score)
    n_reprod_b <- sum(direct_b_results$n_samples_supporting >= 3)
    message(sprintf("  Direction B: %d L-R pairs found (%d reproduced in ≥3/4 samples)",
                    nrow(direct_b_results), n_reprod_b))
  } else {
    message("  Direction B: 0 L-R pairs found")
  }

  fwrite(direct_b_results, file.path(ct_output_dir, "direct_lr_pairs_target_to_hymy.csv"))

  # =========================================================================
  # Part 2: NicheNet Ligand Activity (Direction A only)
  # =========================================================================

  message(paste0("\n  --- Part 2: NicheNet Ligand Activity (", QUERY_LABEL, " -> Target) ---"))

  activity_results <- data.table()

  # Response = gradient genes INDUCED near query (negative coef)
  # Fallback: if <10 induced, use ALL significant gradient genes
  response_induced <- intersect(induced_genes, rownames(ligand_target_matrix))
  response_all_sig <- intersect(sig_genes$gene_safe, rownames(ligand_target_matrix))
  bg_in_matrix <- intersect(background_genes, rownames(ligand_target_matrix))
  potential_ligands <- intersect(hymy_expressed_ligands, colnames(ligand_target_matrix))

  # Decide which response set to use
  if (length(response_induced) >= 10) {
    response_in_matrix <- response_induced
    response_type <- "induced"
  } else if (length(response_all_sig) >= 10) {
    response_in_matrix <- response_all_sig
    response_type <- "all_significant"
    message(sprintf("  NOTE: Only %d induced genes. Using ALL %d significant gradient genes as fallback.",
                    length(response_induced), length(response_all_sig)))
  } else {
    response_in_matrix <- response_induced
    response_type <- "induced"
  }

  message(sprintf("  Response genes (%s, in matrix): %d",
                  response_type, length(response_in_matrix)))
  message(sprintf("  Background genes (in matrix): %d / %d",
                  length(bg_in_matrix), length(background_genes)))
  message(sprintf("  Potential ligands (%s-expressed, in matrix): %d", QUERY_LABEL,
                  length(potential_ligands)))

  if (length(response_in_matrix) >= 10 && length(potential_ligands) >= 5) {
    message(sprintf("  Running predict_ligand_activities() with %s response set...",
                    response_type))

    ligand_activity <- predict_ligand_activities(
      geneset = response_in_matrix,
      background_expressed_genes = bg_in_matrix,
      ligand_target_matrix = ligand_target_matrix,
      potential_ligands = potential_ligands
    )

    ligand_activity <- as.data.table(ligand_activity)
    setnames(ligand_activity, "test_ligand", "ligand")
    setnames(ligand_activity, "aupr_corrected", "activity")

    # Add query cell expression info
    ligand_activity <- merge(ligand_activity, hymy_profile,
                             by.x = "ligand", by.y = "gene", all.x = TRUE)

    # Add cell type and metadata
    ligand_activity[, cell_type := ct_name]
    ligand_activity[, response_type := response_type]
    ligand_activity[, n_response_genes := length(response_in_matrix)]
    ligand_activity[, n_background_genes := length(bg_in_matrix)]

    # Rank
    setorder(ligand_activity, -activity)
    ligand_activity[, rank := .I]

    activity_results <- ligand_activity
    message(sprintf("  Calculated activities for %d ligands", nrow(activity_results)))
    message(sprintf("  Top 5: %s",
                    paste(head(activity_results$ligand, 5), collapse = ", ")))
  } else {
    message("  Insufficient genes for ligand activity scoring. Skipping.")
  }

  fwrite(activity_results, file.path(ct_output_dir, "nichenet_ligand_activity.csv"))

  # =========================================================================
  # Part 3: Downstream Target Enrichment
  # =========================================================================

  message("\n  --- Part 3: Downstream Target Enrichment ---")

  # Collect top ligands to test (union of Part 1 Direction A + Part 2)
  top_ligands_direct <- if (nrow(direct_a_results) > 0) {
    head(unique(direct_a_results$ligand), 30)
  } else {
    character(0)
  }
  top_ligands_activity <- if (nrow(activity_results) > 0) {
    head(activity_results$ligand, 20)
  } else {
    character(0)
  }
  test_ligands <- unique(c(top_ligands_direct, top_ligands_activity))
  # Only test ligands present in the matrix
  test_ligands <- intersect(test_ligands, colnames(ligand_target_matrix))

  message(sprintf("  Testing %d unique ligands for downstream enrichment", length(test_ligands)))

  # Gradient gene sets for this cell type
  gradient_genes_ct <- sig_genes$gene_safe
  gradient_induced_ct <- induced_genes
  gradient_repressed_ct <- repressed_genes
  universe_genes <- intersect(background_genes, rownames(ligand_target_matrix))

  enrichment_results <- data.table()

  if (length(test_ligands) > 0 && length(gradient_genes_ct) > 0) {
    for (lig in test_ligands) {
      # Get predicted downstream targets from ligand-target matrix
      target_weights <- ligand_target_matrix[, lig]
      # Filter to genes in our universe and take top 200
      target_weights <- target_weights[names(target_weights) %in% universe_genes]
      target_weights <- sort(target_weights, decreasing = TRUE)
      predicted_targets <- names(head(target_weights, 200))

      # Overlap with gradient genes
      overlap_all <- intersect(predicted_targets, gradient_genes_ct)
      overlap_induced <- intersect(predicted_targets, gradient_induced_ct)
      overlap_repressed <- intersect(predicted_targets, gradient_repressed_ct)

      # Fisher's exact test (one-sided: enrichment)
      n_universe <- length(universe_genes)
      n_targets <- length(predicted_targets)
      n_gradient <- length(intersect(gradient_genes_ct, universe_genes))
      n_overlap <- length(overlap_all)

      # 2x2 contingency table
      a <- n_overlap
      b <- n_targets - n_overlap
      c <- n_gradient - n_overlap
      d <- n_universe - n_targets - n_gradient + n_overlap

      fisher_res <- tryCatch({
        fisher.test(matrix(c(a, c, b, d), nrow = 2), alternative = "greater")
      }, error = function(e) {
        list(p.value = 1, estimate = 1)
      })

      # Jaccard index
      jaccard <- if ((n_targets + n_gradient - n_overlap) > 0) {
        n_overlap / (n_targets + n_gradient - n_overlap)
      } else {
        0
      }

      # Find corresponding receptor(s) for this ligand
      receptors_for_lig <- lr_network[ligand == lig & receptor %in% available_receptors]$receptor
      receptor_str <- paste(receptors_for_lig, collapse = ";")

      enrichment_results <- rbind(enrichment_results, data.table(
        ligand = lig,
        receptors = receptor_str,
        cell_type = ct_name,
        n_predicted_targets = n_targets,
        n_gradient_genes = n_gradient,
        n_overlap = n_overlap,
        overlap_genes = paste(overlap_all, collapse = ";"),
        n_overlap_induced = length(overlap_induced),
        n_overlap_repressed = length(overlap_repressed),
        pvalue_fisher = fisher_res$p.value,
        odds_ratio = as.numeric(fisher_res$estimate),
        jaccard_index = jaccard
      ))
    }

    # FDR correction
    if (nrow(enrichment_results) > 0) {
      enrichment_results[, fdr_fisher := p.adjust(pvalue_fisher, method = "BH")]
      enrichment_results[, enrichment_score := -log10(pvalue_fisher + 1e-50) *
                           sign(log2(pmax(odds_ratio, 0.01)))]
      setorder(enrichment_results, pvalue_fisher)

      n_sig_enriched <- sum(enrichment_results$fdr_fisher < 0.05)
      message(sprintf("  Enrichment: %d / %d ligands significant (FDR<0.05)",
                      n_sig_enriched, nrow(enrichment_results)))
    }
  } else {
    message("  No ligands to test or no gradient genes. Skipping enrichment.")
  }

  fwrite(enrichment_results, file.path(ct_output_dir, "downstream_target_enrichment.csv"))

  # =========================================================================
  # Combined Prioritization (Direction A: Query -> Target)
  # =========================================================================

  message("\n  --- Combining results ---")

  # Build combined table from Direction A direct pairs
  if (nrow(direct_a_results) > 0) {
    combined <- copy(direct_a_results)

    # Merge NicheNet activity (by ligand)
    if (nrow(activity_results) > 0) {
      combined <- merge(
        combined,
        activity_results[, .(ligand, nichenet_activity = activity, nichenet_auroc = auroc,
                             nichenet_rank = rank)],
        by = "ligand", all.x = TRUE
      )
    } else {
      combined[, `:=`(nichenet_activity = NA_real_, nichenet_auroc = NA_real_,
                      nichenet_rank = NA_integer_)]
    }

    # Merge enrichment (by ligand)
    if (nrow(enrichment_results) > 0) {
      combined <- merge(
        combined,
        enrichment_results[, .(ligand, enrichment_pval = pvalue_fisher,
                               enrichment_fdr = fdr_fisher,
                               enrichment_odds_ratio = odds_ratio,
                               n_downstream_overlap = n_overlap,
                               overlap_genes, enrichment_score)],
        by = "ligand", all.x = TRUE
      )
    } else {
      combined[, `:=`(enrichment_pval = NA_real_, enrichment_fdr = NA_real_,
                      enrichment_odds_ratio = NA_real_,
                      n_downstream_overlap = NA_integer_,
                      overlap_genes = NA_character_,
                      enrichment_score = NA_real_)]
    }

    # Compute combined score (weighted average of rescaled components)
    w_direct <- 4
    w_activity <- 3
    w_enrichment <- 2
    w_expression <- 1

    combined[, scaled_direct := rescale(direct_score, to = c(0, 1), na.rm = TRUE)]
    combined[, scaled_activity := fifelse(
      is.na(nichenet_activity), 0,
      rescale(nichenet_activity, to = c(0, 1), na.rm = TRUE)
    )]
    combined[, scaled_enrichment := fifelse(
      is.na(enrichment_score), 0,
      rescale(pmax(enrichment_score, 0), to = c(0, 1), na.rm = TRUE)
    )]
    combined[, scaled_expression := rescale(ligand_pct_hymy / 100, to = c(0, 1), na.rm = TRUE)]

    combined[, combined_score := (
      w_direct * scaled_direct +
        w_activity * scaled_activity +
        w_enrichment * scaled_enrichment +
        w_expression * scaled_expression
    ) / (w_direct + w_activity + w_enrichment + w_expression)]

    setorder(combined, -combined_score)
    message(sprintf("  Combined: %d L-R pairs scored", nrow(combined)))
    message(sprintf("  Top 5: %s",
                    paste(head(paste0(combined$ligand, "-", combined$receptor), 5),
                          collapse = ", ")))
  } else {
    combined <- data.table()
    message("  No Direction A pairs to combine.")
  }

  fwrite(combined, file.path(ct_output_dir, "combined_prioritization.csv"))

  message(sprintf("\n  Outputs saved to: %s", ct_output_dir))
}

# =============================================================================
# Summary Statistics
# =============================================================================

message("\n", strrep("=", 70))
message("Analysis Complete")
message(strrep("=", 70))
message(sprintf("Output: %s", OUTPUT_BASE))
message(sprintf("Cell types processed: %d", length(TARGET_CELLTYPES)))
message(sprintf("Timestamp: %s", Sys.time()))
