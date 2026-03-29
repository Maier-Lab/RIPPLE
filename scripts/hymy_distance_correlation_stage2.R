#!/usr/bin/env Rscript
#' =============================================================================
#' RIPPLE Stage 4: Confounder Control (Bivariate Poisson GLM)
#' =============================================================================
#'
#' Validates Stage 1 hits by adding distance-to-nearest-control-cell-type as a
#' covariate, isolating query-specific effects from general niche effects.
#'
#' For each significant gene from Stage 1, fits a bivariate Poisson GLM
#' per sample:
#'   glm(counts ~ dist_to_query + dist_to_control + offset(log(total_counts)),
#'       family = poisson)
#'
#' Gene Classification:
#' - query_specific: significant in both Stage 1 and Stage 2 (gradient persists)
#' - niche_driven: significant in Stage 1, NOT in Stage 2 (explained by niche)
#' - enhanced: larger absolute effect in Stage 2 than Stage 1
#'
#' Usage:
#'   Rscript hymy_distance_correlation_stage2.R
#'   QUERY_CELLTYPE=MyType CELLTYPE_COLUMN=my_col Rscript hymy_distance_correlation_stage2.R
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
  library(patchwork)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(RANN)
  library(parallel)
  library(viridis)
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

# Stage 1 analysis name (set via ANALYSIS_NAME env var)
# Default: hymy_distance_correlation_v2 (Poisson GLM)
STAGE1_ANALYSIS_NAME <- Sys.getenv("ANALYSIS_NAME", unset = "hymy_distance_correlation_v2")
ANALYSIS_NAME <- paste0(STAGE1_ANALYSIS_NAME, "_stage2")

# Significance column in Stage 1 results (fisher_fdr for v2, fdr for v1)
SIG_COL <- "fisher_fdr"

# Inherited from config.R (via utils.R): QUERY_CELLTYPE, CELLTYPE_COL, OUTPUT_SUFFIX, QUERY_LABEL
OUTPUT_BASE <- file.path(OUTPUT_ROOT, ANALYSIS_NAME)
STAGE1_BASE <- file.path(OUTPUT_ROOT, STAGE1_ANALYSIS_NAME)

ensure_dir(OUTPUT_BASE)
ensure_dir(file.path(OUTPUT_BASE, "per_celltype"))
ensure_dir(file.path(OUTPUT_BASE, "summary"))
ensure_dir(file.path(OUTPUT_BASE, "plots"))

# Statistical parameters (match Stage 1)
MIN_CELLS_PER_SAMPLE <- 30
MIN_EXPR_CELLS <- 5
FDR_THRESHOLD <- 0.05
MIN_CONTROL_CELLS <- 30     # Minimum control cells per sample for reliable nn2()

# Distance parameters (match Stage 1)
MAX_DISTANCE_UM <- 200

# Control cell type for distance covariate
# For Monocyte target analysis, use Macrophages as alternative control
CONTROL_CELLTYPE <- "Monocyte"

# Target cell types (same as Stage 1)
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

# =============================================================================
# Cell Type Selection (for SLURM array job parallelization)
# =============================================================================

CELLTYPE_INDEX <- as.integer(Sys.getenv("CELLTYPE_INDEX", unset = "0"))

if (CELLTYPE_INDEX > 0) {
  celltype_names <- names(TARGET_CELLTYPES)
  if (CELLTYPE_INDEX > length(celltype_names)) {
    stop("CELLTYPE_INDEX ", CELLTYPE_INDEX, " exceeds number of cell types (", length(celltype_names), ")")
  }
  selected_ct <- celltype_names[CELLTYPE_INDEX]
  message("\n>>> Array job mode: Processing only ", selected_ct, " (index ", CELLTYPE_INDEX, ")")
  TARGET_CELLTYPES <- TARGET_CELLTYPES[selected_ct]
}

message(strrep("=", 70))
message(paste0(QUERY_LABEL, " Distance Correlation Stage 2: Monocyte Distance Control (Poisson)"))
message(strrep("=", 70))
message("Annotation level: ", ANNOTATION_LEVEL)
message("Query cell type: ", QUERY_CELLTYPE)
message("Control cell type: ", CONTROL_CELLTYPE)
message("Stage 1 analysis: ", STAGE1_ANALYSIS_NAME)
message("Stage 1 results: ", STAGE1_BASE)
message("Significance column: ", SIG_COL)
message("Output directory: ", OUTPUT_BASE)

# =============================================================================
# Helper Functions
# =============================================================================

#' Fit bivariate Poisson GLM with distance control and cell-size offset
#'
#' Models counts ~ dist_hymy + dist_control + offset(log(total_counts))
#' using Poisson regression. Extracts the coefficient for dist_hymy only
#' (partial effect, in log-rate/µm units matching Stage 1 v2).
#'
#' @param count_vec Vector of raw integer counts for a single gene
#' @param dist_hymy Vector of distances to nearest HyMy cell
#' @param dist_control Vector of distances to nearest control cell (monocyte)
#' @param total_counts Vector of total counts per cell (for offset)
#' @param min_cells Minimum cells for analysis
#' @return List with coef, se, n_cells, pval, dispersion for dist_hymy term
fit_poisson_controlled <- function(count_vec, dist_hymy, dist_control, total_counts,
                                    min_cells = MIN_EXPR_CELLS) {
  na_result <- list(coef = NA_real_, se = NA_real_, n_cells = 0L,
                    pval = NA_real_, dispersion = NA_real_)

  # Remove NAs and invalid values
  valid_idx <- !is.na(count_vec) & !is.na(dist_hymy) & !is.na(dist_control) &
               is.finite(count_vec) & is.finite(dist_hymy) & is.finite(dist_control) &
               !is.na(total_counts) & total_counts > 0
  count_vec <- count_vec[valid_idx]
  dist_hymy <- dist_hymy[valid_idx]
  dist_control <- dist_control[valid_idx]
  log_total <- log(total_counts[valid_idx])

  if (length(count_vec) < min_cells) return(na_result)

  # Need some non-zero counts for the model to be meaningful
  if (sum(count_vec > 0) < min_cells) return(na_result)

  # Fit bivariate Poisson GLM with cell-size offset
  fit <- tryCatch({
    suppressWarnings(glm(count_vec ~ dist_hymy + dist_control + offset(log_total),
                         family = poisson))
  }, error = function(e) NULL)

  if (is.null(fit) || !fit$converged) return(na_result)

  coef_summary <- summary(fit)$coefficients
  if (!"dist_hymy" %in% rownames(coef_summary)) return(na_result)

  # Overdispersion diagnostic
  dispersion <- fit$deviance / fit$df.residual

  list(
    coef = coef_summary["dist_hymy", "Estimate"],      # log-rate change per µm
    se = coef_summary["dist_hymy", "Std. Error"],
    n_cells = length(count_vec),
    pval = coef_summary["dist_hymy", "Pr(>|z|)"],      # Wald z-test
    dispersion = dispersion
  )
}


#' Combine per-sample results using Fisher's combined p-value
#'
#' Equal-weight combination across mice with sign consistency gate.
#' Matches Stage 1 v2 methodology (recompute_meta_summary.R).
#'
#' @param coefs Vector of per-sample coefficient estimates
#' @param pvals Vector of per-sample p-values
#' @param sample_ids Vector of sample identifiers
#' @return List with median_coef, fisher_pval, sign_consistency, n_samples
combine_with_fisher <- function(coefs, pvals, sample_ids) {
  valid <- !is.na(coefs) & !is.na(pvals) & pvals > 0

  if (sum(valid) < 2) {
    return(list(
      median_coef = NA_real_, fisher_pval = NA_real_,
      sign_consistency = NA_real_, n_samples = sum(valid)
    ))
  }

  v_coefs <- coefs[valid]
  v_pvals <- pvals[valid]
  n_valid <- length(v_coefs)

  # Sign consistency gate
  n_pos <- sum(v_coefs > 0)
  n_neg <- sum(v_coefs < 0)
  sign_consistency <- max(n_pos, n_neg) / n_valid

  # Fisher's combined p-value (only if all mice agree on direction)
  if (sign_consistency < 1.0) {
    fisher_pval <- 1.0
  } else {
    clamped <- pmax(v_pvals, 1e-15)
    fisher_stat <- -2 * sum(log(clamped))
    fisher_pval <- pchisq(fisher_stat, df = 2 * n_valid, lower.tail = FALSE)
  }

  list(
    median_coef = median(v_coefs),
    fisher_pval = fisher_pval,
    sign_consistency = sign_consistency,
    n_samples = n_valid
  )
}


# =============================================================================
# Load Data
# =============================================================================

message("\n", strrep("=", 70))
message("Loading Data")
message(strrep("=", 70))

# Load Stage 1 merged results
stage1_summary_file <- file.path(STAGE1_BASE, "summary", "all_genes_results.csv")
if (!file.exists(stage1_summary_file)) {
  stop("Stage 1 merged results not found: ", stage1_summary_file,
       "\nRun merge_distance_correlation_results.R first.")
}

stage1_all <- fread(stage1_summary_file)
message("Loaded Stage 1 results: ", nrow(stage1_all), " gene-celltype entries")

# Verify significance column exists
if (!SIG_COL %in% names(stage1_all)) {
  message("WARNING: '", SIG_COL, "' not found in Stage 1 results. Available: ",
          paste(grep("fdr|pval", names(stage1_all), value = TRUE), collapse = ", "))
  # Fallback to 'fdr' if fisher_fdr not available
  SIG_COL <- "fdr"
  message("  Falling back to '", SIG_COL, "'")
}

# Filter to significant genes
stage1_sig <- stage1_all[get(SIG_COL) < FDR_THRESHOLD]
message("Stage 1 significant (", SIG_COL, " < ", FDR_THRESHOLD, "): ", nrow(stage1_sig), " entries")
message("  Cell types: ", paste(unique(stage1_sig$cell_type), collapse = ", "))

if (nrow(stage1_sig) == 0) {
  message("No significant genes from Stage 1 — nothing to validate.")
  quit(save = "no", status = 0)
}

# Load Seurat object
obj <- load_seurat()

if (ANNOTATION_LEVEL == "HyMy") {
  message("Merging HyMy annotations...")
  obj <- merge_hymy_annotations(obj)
}

# Extract metadata
cell_data <- as.data.table(obj@meta.data, keep.rownames = "barcode")

if ("group" %in% names(cell_data)) {
  cell_data[, condition := group]
  message("Using 'group' column for condition")
} else if (!"condition" %in% names(cell_data)) {
  cell_data[, condition := sapply(sample_id, get_condition)]
  message("WARNING: Inferring condition from sample_id")
}

if (!CELLTYPE_COL %in% names(cell_data)) {
  stop("Cell type column not found: ", CELLTYPE_COL)
}

# Filter to TDLN (same as Stage 1)
message("\nFiltering to TDLN samples only...")
cell_data <- cell_data[condition == "TDLN"]

message("Data summary (TDLN only):")
message("  Total cells: ", nrow(cell_data))
message("  Samples: ", length(unique(cell_data$sample_id)))

# =============================================================================
# Calculate Distances to Query Cells
# =============================================================================

message("\n", strrep("=", 70))
message("Calculating Distances")
message(strrep("=", 70))

coords <- as.matrix(cell_data[, .(spatial_x, spatial_y)])

# Distance to query cell type
query_mask <- cell_data[[CELLTYPE_COL]] == QUERY_CELLTYPE
message("Query cells (", QUERY_CELLTYPE, "): ", sum(query_mask))
query_coords <- coords[query_mask, , drop = FALSE]
nn_hymy <- nn2(query_coords, coords, k = 1)
cell_data[, dist_to_hymy := pmin(as.vector(nn_hymy$nn.dists), MAX_DISTANCE_UM)]

# =============================================================================
# Calculate Distances to Control Cells (Monocyte or Macrophages)
# =============================================================================

# Determine control cell type per target
# Default is Monocyte, but for Monocyte target we use Macrophages
message("\nCalculating distance to control cells...")

control_mask <- cell_data[[CELLTYPE_COL]] == CONTROL_CELLTYPE
n_control_total <- sum(control_mask)
message("Control cells (", CONTROL_CELLTYPE, "): ", n_control_total)

if (n_control_total < MIN_CONTROL_CELLS) {
  stop("Too few control cells (", n_control_total, ") for reliable distance calculation")
}

control_coords <- coords[control_mask, , drop = FALSE]
nn_control <- nn2(control_coords, coords, k = 1)
cell_data[, dist_to_control := pmin(as.vector(nn_control$nn.dists), MAX_DISTANCE_UM)]

# For Monocyte target: use Macrophages as alternative control
alt_control <- "Macrophages"
alt_control_mask <- cell_data[[CELLTYPE_COL]] == alt_control
n_alt_control <- sum(alt_control_mask)
message("Alternative control cells (", alt_control, "): ", n_alt_control)

if (n_alt_control >= MIN_CONTROL_CELLS) {
  alt_control_coords <- coords[alt_control_mask, , drop = FALSE]
  nn_alt <- nn2(alt_control_coords, coords, k = 1)
  cell_data[, dist_to_alt_control := pmin(as.vector(nn_alt$nn.dists), MAX_DISTANCE_UM)]
} else {
  message("  WARNING: Too few alternative control cells for Monocyte target analysis")
  cell_data[, dist_to_alt_control := NA_real_]
}

# Per-sample control cell counts (for flagging unreliable samples)
control_per_sample <- cell_data[control_mask == TRUE, .N, by = sample_id]
setnames(control_per_sample, "N", "n_control")
message("\nControl cells per sample:")
for (i in seq_len(nrow(control_per_sample))) {
  flag <- if (control_per_sample$n_control[i] < MIN_CONTROL_CELLS) " [FLAGGED: too few]" else ""
  message("  ", control_per_sample$sample_id[i], ": ", control_per_sample$n_control[i], flag)
}

# =============================================================================
# Collinearity Check
# =============================================================================

message("\n", strrep("=", 70))
message("Collinearity Diagnostics")
message(strrep("=", 70))

# Check correlation between dist_to_hymy and dist_to_control per sample
samples <- unique(cell_data$sample_id)
collinearity_report <- rbindlist(lapply(samples, function(samp) {
  samp_data <- cell_data[sample_id == samp]
  cor_val <- cor(samp_data$dist_to_hymy, samp_data$dist_to_control,
                 use = "complete.obs", method = "pearson")
  data.table(
    sample_id = samp,
    cor_hymy_control = round(cor_val, 3)
  )
}))

message("Pearson correlation (dist_to_hymy vs dist_to_control):")
for (i in seq_len(nrow(collinearity_report))) {
  flag <- if (abs(collinearity_report$cor_hymy_control[i]) > 0.8) " [WARNING: high collinearity]" else ""
  message("  ", collinearity_report$sample_id[i], ": r = ",
          collinearity_report$cor_hymy_control[i], flag)
}
message("  Mean: r = ", round(mean(collinearity_report$cor_hymy_control), 3))

# Save collinearity report
fwrite(collinearity_report, file.path(OUTPUT_BASE, "summary", "collinearity_report.csv"))

high_collinearity <- any(abs(collinearity_report$cor_hymy_control) > 0.8)
if (high_collinearity) {
  message("\n  >>> WARNING: High collinearity detected in some samples.")
  message("  >>> Bivariate model SEs will be inflated. Interpret with caution.")
}

# =============================================================================
# Run Stage 2 Analysis per Cell Type
# =============================================================================

message("\n", strrep("=", 70))
message("Running Stage 2 Bivariate Analysis")
message(strrep("=", 70))

all_stage2_results <- list()

for (ct_name in names(TARGET_CELLTYPES)) {
  ct_types <- TARGET_CELLTYPES[[ct_name]]

  message("\n", strrep("-", 60))
  message("Analyzing: ", ct_name)
  message(strrep("-", 60))

  # Determine which control distance to use
  if (ct_name == "Monocyte") {
    if (all(is.na(cell_data$dist_to_alt_control))) {
      message("  SKIPPING: Monocyte target cannot use Monocyte as control,")
      message("  and alternative control (Macrophages) has too few cells.")
      next
    }
    dist_control_col <- "dist_to_alt_control"
    control_label <- alt_control
    message("  Using alternative control: ", alt_control, " (cannot use Monocyte as both target and control)")
  } else {
    dist_control_col <- "dist_to_control"
    control_label <- CONTROL_CELLTYPE
  }

  # Get Stage 1 significant genes for this cell type
  ct_stage1 <- stage1_sig[cell_type == ct_name]
  if (nrow(ct_stage1) == 0) {
    message("  No Stage 1 significant genes for ", ct_name, " — skipping")
    next
  }
  sig_genes <- ct_stage1$gene
  message("  Stage 1 significant genes: ", length(sig_genes))

  # Identify target cells (TDLN only)
  cell_data[, is_target := get(CELLTYPE_COL) %in% ct_types & condition == "TDLN"]
  target_data <- cell_data[is_target == TRUE]
  message("  Target cells: ", nrow(target_data))

  if (nrow(target_data) < MIN_CELLS_PER_SAMPLE * 2) {
    message("  Insufficient cells — skipping")
    next
  }

  # Valid samples
  sample_counts <- target_data[, .N, by = sample_id]
  valid_samples <- sample_counts[N >= MIN_CELLS_PER_SAMPLE]$sample_id

  # Also require sufficient control cells per sample
  valid_control_samples <- control_per_sample[n_control >= MIN_CONTROL_CELLS]$sample_id
  if (ct_name == "Monocyte") {
    # For monocyte target, check alt control counts
    alt_per_sample <- cell_data[alt_control_mask == TRUE, .N, by = sample_id]
    setnames(alt_per_sample, "N", "n_alt_control")
    valid_control_samples <- alt_per_sample[n_alt_control >= MIN_CONTROL_CELLS]$sample_id
  }
  valid_samples <- intersect(valid_samples, valid_control_samples)

  message("  Valid samples (>= ", MIN_CELLS_PER_SAMPLE, " target cells + >= ",
          MIN_CONTROL_CELLS, " control cells): ", length(valid_samples))

  if (length(valid_samples) < 2) {
    message("  Need at least 2 valid samples — skipping")
    next
  }

  # Get target cell barcodes
  target_valid <- target_data[sample_id %in% valid_samples]
  target_barcodes <- target_valid$barcode

  # Get raw count matrix (Poisson GLM uses counts, not normalized data)
  count_matrix <- GetAssayData(obj, layer = "counts")[, target_barcodes, drop = FALSE]

  # Compute total counts per cell (for Poisson offset)
  total_counts_target <- colSums(count_matrix)
  message("  Total counts per cell: median=", round(median(total_counts_target)),
          ", range=[", round(min(total_counts_target)), "-", round(max(total_counts_target)), "]")

  # Filter to genes that exist in the count matrix
  sig_genes <- intersect(sig_genes, rownames(count_matrix))
  message("  Genes in count matrix: ", length(sig_genes))

  if (length(sig_genes) == 0) {
    message("  No genes available — skipping")
    next
  }

  ct_output_dir <- file.path(OUTPUT_BASE, "per_celltype", ct_name)
  ensure_dir(ct_output_dir)

  # Step 1: Bivariate Poisson GLM per gene per sample
  message("  Step 1: Fitting bivariate Poisson GLM...")

  coef_results <- rbindlist(lapply(sig_genes, function(g) {
    gene_counts <- as.numeric(count_matrix[g, target_barcodes])

    sample_results <- rbindlist(lapply(valid_samples, function(samp) {
      samp_idx <- which(target_valid$sample_id == samp)
      if (length(samp_idx) < MIN_CELLS_PER_SAMPLE) {
        return(data.table(
          gene = g, sample_id = samp,
          coef = NA_real_, se = NA_real_, n_cells = length(samp_idx),
          pval = NA_real_, dispersion = NA_real_
        ))
      }

      samp_counts <- gene_counts[samp_idx]
      samp_dist_hymy <- target_valid[samp_idx]$dist_to_hymy
      samp_dist_control <- target_valid[[dist_control_col]][samp_idx]
      samp_total <- total_counts_target[target_barcodes[samp_idx]]

      fit_result <- fit_poisson_controlled(samp_counts, samp_dist_hymy,
                                            samp_dist_control, samp_total)

      data.table(
        gene = g, sample_id = samp,
        coef = fit_result$coef, se = fit_result$se,
        n_cells = fit_result$n_cells, pval = fit_result$pval,
        dispersion = fit_result$dispersion
      )
    }))

    sample_results
  }), fill = TRUE)

  fwrite(coef_results, file.path(ct_output_dir, "coef_per_sample.csv"))
  message("  Saved: coef_per_sample.csv")

  # Step 2: Fisher's combined p-value across samples
  message("  Step 2: Combining with Fisher's method...")

  fisher_results <- rbindlist(lapply(sig_genes, function(g) {
    gene_data <- coef_results[gene == g]

    result <- combine_with_fisher(
      coefs = gene_data$coef,
      pvals = gene_data$pval,
      sample_ids = gene_data$sample_id
    )

    data.table(
      gene = g,
      stage2_median_coef = result$median_coef,
      stage2_fisher_pval = result$fisher_pval,
      stage2_sign_consistency = result$sign_consistency,
      stage2_n_samples = result$n_samples
    )
  }), fill = TRUE)

  fisher_results[, stage2_fisher_fdr := p.adjust(stage2_fisher_pval, method = "BH")]

  # Step 3: Merge with Stage 1 results
  message("  Step 3: Comparing Stage 1 vs Stage 2...")

  # Build Stage 1 columns dynamically based on what's available
  # Maps: new_name -> old_name_in_stage1
  col_map <- c(gene = "gene")
  if ("median_coef" %in% names(ct_stage1)) {
    col_map["stage1_coef"] <- "median_coef"
  } else {
    col_map["stage1_coef"] <- "combined_coef"
  }
  if ("combined_se" %in% names(ct_stage1)) col_map["stage1_se"] <- "combined_se"
  if ("fisher_pval" %in% names(ct_stage1)) {
    col_map["stage1_pval"] <- "fisher_pval"
  } else {
    col_map["stage1_pval"] <- "pval"
  }
  if ("fisher_fdr" %in% names(ct_stage1)) {
    col_map["stage1_fdr"] <- "fisher_fdr"
  } else {
    col_map["stage1_fdr"] <- "fdr"
  }
  if ("i2" %in% names(ct_stage1)) col_map["stage1_i2"] <- "i2"

  # Build the stage1 subset with renamed columns
  src_cols <- unname(col_map)
  ct_stage1_sub <- ct_stage1[, ..src_cols]
  setnames(ct_stage1_sub, src_cols, names(col_map))

  comparison <- merge(
    ct_stage1_sub,
    fisher_results,
    by = "gene",
    all.x = TRUE
  )
  setDT(comparison)

  # Convenience alias for classification code
  comparison[, stage2_coef := stage2_median_coef]

  # Classify genes
  # Key distinction: "niche_driven" requires coefficient ATTENUATION, not just loss of significance
  # If coefficient is similar but SE increased (power loss), classify as "underpowered"
  ATTENUATION_THRESHOLD <- 0.5  # Coefficient must drop to <50% of Stage 1 to be "niche_driven"

  comparison[, coef_ratio := abs(stage2_coef) / abs(stage1_coef)]

  query_specific_label <- paste0(QUERY_LABEL, "_specific")
  comparison[, classification := fcase(
    is.na(stage2_fisher_fdr), "no_stage2_result",
    # query_specific: significant in Stage 2, same direction
    stage2_fisher_fdr < FDR_THRESHOLD & sign(stage2_coef) == sign(stage1_coef), query_specific_label,
    # Niche-driven: lost significance AND coefficient attenuated (effect explained by control)
    stage2_fisher_fdr >= FDR_THRESHOLD & coef_ratio < ATTENUATION_THRESHOLD, "niche_driven",
    # Underpowered: lost significance but coefficient NOT attenuated (likely SE inflation)
    stage2_fisher_fdr >= FDR_THRESHOLD & coef_ratio >= ATTENUATION_THRESHOLD, "underpowered",
    # Reversed: significant but opposite direction (rare edge case)
    stage2_fisher_fdr < FDR_THRESHOLD & sign(stage2_coef) != sign(stage1_coef), "reversed",
    # Fallback
    default = "unclassified"
  )]

  # Detect "enhanced" genes: absolute effect LARGER in Stage 2 than Stage 1
  # (control proximity was suppressing the query-specific signal)
  comparison[classification == query_specific_label &
               abs(stage2_coef) > abs(stage1_coef) * 1.1,
             classification := "enhanced"]

  # Clean up temporary columns
  comparison[, coef_ratio := NULL]
  comparison[, stage2_coef := NULL]

  # Add cell type and control label
  comparison[, cell_type := ct_name]
  comparison[, control_celltype := control_label]

  # Save per-celltype results
  fwrite(comparison, file.path(ct_output_dir, "stage2_comparison.csv"))
  message("  Saved: stage2_comparison.csv")

  # Print classification summary
  class_summary <- comparison[, .N, by = classification]
  setorder(class_summary, -N)
  message("\n  Classification summary:")
  for (i in seq_len(nrow(class_summary))) {
    message("    ", class_summary$classification[i], ": ", class_summary$N[i])
  }

  # Step 4: Visualizations
  message("  Step 4: Creating visualizations...")

  # 4a. Scatter plot: Stage 1 coef vs Stage 2 coef
  plot_data <- comparison[!is.na(stage2_median_coef)]

  if (nrow(plot_data) > 0) {
    # Determine max range for symmetric axes
    max_range <- max(abs(c(plot_data$stage1_coef, plot_data$stage2_median_coef)), na.rm = TRUE) * 1.1

    # Color palette for classification
    class_colors <- setNames(
      c("#E74C3C", "#8E44AD", "#3498DB", "#F1C40F", "#F39C12", "grey70"),
      c(query_specific_label, "enhanced", "niche_driven", "underpowered", "reversed", "no_stage2_result")
    )

    # Identify genes to label (top 15 by Stage 1 significance)
    label_genes <- head(plot_data[order(stage1_fdr)], 15)

    p_scatter <- ggplot(plot_data, aes(x = stage1_coef, y = stage2_median_coef, color = classification)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
      geom_point(alpha = 0.7, size = 2) +
      scale_color_manual(values = class_colors, name = "Classification") +
      geom_text_repel(
        data = label_genes,
        aes(label = gene),
        size = 3, max.overlaps = 15, box.padding = 0.4
      ) +
      xlim(-max_range, max_range) +
      ylim(-max_range, max_range) +
      labs(
        title = sprintf("Stage 1 vs Stage 2 Coefficients: %s", ct_name),
        subtitle = sprintf("Control: %s | diagonal = no change", control_label),
        x = "Stage 1 coefficient (log-rate per µm)",
        y = "Stage 2 coefficient (adjusted, log-rate per µm)"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        aspect.ratio = 1
      )

    ggsave(file.path(ct_output_dir, "stage1_vs_stage2_scatter.pdf"),
           p_scatter, width = 8, height = 7)
    message("  Saved: stage1_vs_stage2_scatter.pdf")

    # 4b. Volcano plot of Stage 2 adjusted results
    volcano_data <- copy(plot_data)
    volcano_data[, neg_log10_fdr := -log10(stage2_fisher_fdr)]
    volcano_data[neg_log10_fdr > 50, neg_log10_fdr := 50]
    volcano_data[, significant := stage2_fisher_fdr < FDR_THRESHOLD]

    top_volcano <- head(volcano_data[significant == TRUE][order(stage2_fisher_fdr)], 20)
    max_coef <- max(abs(volcano_data$stage2_median_coef), na.rm = TRUE) * 1.1

    p_volcano <- ggplot(volcano_data, aes(x = stage2_median_coef, y = neg_log10_fdr)) +
      geom_point(aes(color = classification, size = significant), alpha = 0.6) +
      scale_color_manual(values = class_colors, name = "Classification") +
      scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2.5), guide = "none") +
      geom_hline(yintercept = -log10(FDR_THRESHOLD), linetype = "dashed", color = "grey40") +
      geom_vline(xintercept = 0, linetype = "solid", color = "grey60") +
      xlim(-max_coef, max_coef) +
      geom_text_repel(
        data = top_volcano,
        aes(label = gene), size = 3, max.overlaps = 20, box.padding = 0.5
      ) +
      labs(
        title = sprintf("Stage 2 Controlled Analysis: %s", ct_name),
        subtitle = sprintf("Control: %s | %d genes significant (Fisher FDR < %.2f)",
                          control_label,
                          sum(volcano_data$significant, na.rm = TRUE),
                          FDR_THRESHOLD),
        x = paste0("Adjusted log-rate coefficient (negative = ", QUERY_LABEL, "-induced)"),
        y = "-log10(Fisher FDR)"
      ) +
      theme_classic(base_size = 12) +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold")
      )

    ggsave(file.path(ct_output_dir, "stage2_volcano.pdf"),
           p_volcano, width = 10, height = 8)
    message("  Saved: stage2_volcano.pdf")
  }

  all_stage2_results[[ct_name]] <- comparison
  message("  Analysis complete for ", ct_name)
}

# =============================================================================
# Combined Summary (skip in array mode — use merge script instead)
# =============================================================================

# In array mode, each job only processes one cell type, so combined summary
# would only contain that one cell type. Skip it here and use a separate
# merge step after all array jobs complete.
SKIP_COMBINED <- CELLTYPE_INDEX > 0

if (SKIP_COMBINED) {
  message("\n[NOTE] Array mode: Skipping combined summary (run merge after all jobs complete)")
  message("  Per-celltype results saved to: ", ct_output_dir)
} else {
  message("\n", strrep("=", 70))
  message("Creating Combined Summary")
  message(strrep("=", 70))
}

if (length(all_stage2_results) > 0 && !SKIP_COMBINED) {
  combined <- rbindlist(all_stage2_results, fill = TRUE)
  fwrite(combined, file.path(OUTPUT_BASE, "summary", "stage2_all_results.csv"))
  message("Saved: summary/stage2_all_results.csv")

  # Classification summary across all cell types
  class_summary <- combined[, .N, by = .(cell_type, classification)]
  class_wide <- dcast(class_summary, cell_type ~ classification, value.var = "N", fill = 0)
  fwrite(class_wide, file.path(OUTPUT_BASE, "summary", "classification_summary.csv"))
  message("Saved: summary/classification_summary.csv")

  # Combined scatter plot (all cell types)
  plot_all <- combined[!is.na(stage2_median_coef)]

  if (nrow(plot_all) > 0) {
    class_colors <- setNames(
      c("#E74C3C", "#8E44AD", "#3498DB", "#F39C12", "grey70"),
      c(query_specific_label, "enhanced", "niche_driven", "reversed", "no_stage2_result")
    )

    max_range <- max(abs(c(plot_all$stage1_coef, plot_all$stage2_median_coef)), na.rm = TRUE) * 1.1

    p_combined <- ggplot(plot_all, aes(x = stage1_coef, y = stage2_median_coef, color = classification)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
      geom_point(alpha = 0.5, size = 1.5) +
      scale_color_manual(values = class_colors, name = "Classification") +
      xlim(-max_range, max_range) +
      ylim(-max_range, max_range) +
      facet_wrap(~ cell_type, scales = "free") +
      labs(
        title = "Stage 1 vs Stage 2: All Cell Types",
        subtitle = "diagonal = no change after monocyte control",
        x = "Stage 1 coefficient",
        y = "Stage 2 coefficient (controlled)"
      ) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold")
      )

    n_facets <- length(unique(plot_all$cell_type))
    ncols <- min(3, n_facets)
    nrows <- ceiling(n_facets / ncols)

    ggsave(file.path(OUTPUT_BASE, "plots", "stage1_vs_stage2_all_celltypes.pdf"),
           p_combined, width = 5 * ncols, height = 4.5 * nrows)
    message("Saved: plots/stage1_vs_stage2_all_celltypes.pdf")

    # Classification bar plot
    class_long <- combined[, .N, by = .(cell_type, classification)]

    p_bar <- ggplot(class_long, aes(x = cell_type, y = N, fill = classification)) +
      geom_col(position = "stack") +
      scale_fill_manual(values = class_colors, name = "Classification") +
      labs(
        title = "Gene Classification: Stage 2 Validation",
        x = NULL,
        y = "Number of Genes"
      ) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold")
      )

    ggsave(file.path(OUTPUT_BASE, "plots", "classification_barplot.pdf"),
           p_bar, width = 10, height = 6)
    message("Saved: plots/classification_barplot.pdf")
  }

  # ==========================================================================
  # Positive Control Validation
  # ==========================================================================
  message("\n", strrep("-", 60))
  message("Positive Control Validation")
  message(strrep("-", 60))
  message(paste0("Expected: CSF3, IL33, CXCL12 should remain ", QUERY_LABEL, "-specific (not niche-driven)"))

  POSITIVE_CONTROLS <- c("Csf3", "Il33", "Cxcl12")

  for (ct in c("FRC", "LEC")) {
    if (ct %in% names(all_stage2_results)) {
      ct_result <- all_stage2_results[[ct]]
      for (ctrl_gene in POSITIVE_CONTROLS) {
        ctrl_row <- ct_result[gene == ctrl_gene]
        if (nrow(ctrl_row) > 0) {
          message(sprintf("  %s in %s: stage1=%.4f, stage2=%.4f, class=%s (Fisher_FDR=%.2e)",
                          ctrl_gene, ct,
                          ctrl_row$stage1_coef, ctrl_row$stage2_median_coef,
                          ctrl_row$classification, ctrl_row$stage2_fisher_fdr))
        } else {
          message(sprintf("  %s in %s: not in Stage 1 significant set", ctrl_gene, ct))
        }
      }
    }
  }

  # Print overall summary
  message("\n", strrep("=", 70))
  message("Analysis Summary")
  message(strrep("=", 70))

  for (ct in names(all_stage2_results)) {
    ct_result <- all_stage2_results[[ct]]
    n_total <- nrow(ct_result)
    n_specific <- sum(ct_result$classification %in% c(query_specific_label, "enhanced"), na.rm = TRUE)
    n_niche <- sum(ct_result$classification == "niche_driven", na.rm = TRUE)
    pct_specific <- round(n_specific / n_total * 100, 1)

    message(sprintf("  %s: %d genes → %d %s-specific (%.1f%%), %d niche-driven",
                    ct, n_total, n_specific, QUERY_LABEL, pct_specific, n_niche))
  }
}

message("\n", strrep("=", 70))
message("Stage 2 Analysis Complete!")
message("Output directory: ", OUTPUT_BASE)
message(strrep("=", 70))
