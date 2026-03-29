#!/bin/bash
#SBATCH --job-name=hymy_dist_corr
#SBATCH --output=/nobackup/lab_maier/Projects/mXenium/CMM/results/logs/hymy_dist_corr_%A_%a.out
#SBATCH --error=/nobackup/lab_maier/Projects/mXenium/CMM/results/logs/hymy_dist_corr_%A_%a.err
#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00
#SBATCH --mem=32g
#SBATCH --array=1-12
#SBATCH --mail-type=END
#SBATCH --mail-user=cmangana@cemm.at

# =============================================================================
# HyMy Distance Correlation (Logistic Regression) - SLURM Array Job
# =============================================================================
# Uses logistic regression: P(expressing) ~ distance to HyMy
# Gradient score = log-odds coefficient (negative = HyMy-induced)
#
# Runs each cell type as a separate job for parallelization.
#
# Cell type mapping:
#   1 = LEC, 2 = FRC, 3 = BEC, 4 = CD4_T_cells, 5 = CD8_T_cells,
#   6 = gdT_cells, 7 = Macrophages, 8 = Monocyte, 9 = Fibroblasts_mac,
#   10 = cDC1, 11 = cDC2, 12 = mature_migDC, 13 = B_cells, 14 = Plasma_cell
#
# Usage:
#   # Run all 9 cell types in parallel
#   sbatch run_hymy_distance_correlation_array.sh
#
#   # Run with L1 annotation (IL1B_myeloid as query)
#   ANNOTATION_LEVEL=L1 sbatch run_hymy_distance_correlation_array.sh
#
#   # Run a single cell type (e.g., FRC only)
#   sbatch --array=2 run_hymy_distance_correlation_array.sh
#
# After completion, merge results:
#   Rscript merge_distance_correlation_results.R
#
# Expected Runtime: ~4-8 hours per cell type (permutation testing is slow)
# =============================================================================

set -e

# Query cell type configuration (pass through to R/Python)
export QUERY_CELLTYPE=${QUERY_CELLTYPE:-}
export CELLTYPE_COLUMN=${CELLTYPE_COLUMN:-}
export QUERY_LABEL=${QUERY_LABEL:-}
export ANNOTATION_LEVEL=${ANNOTATION_LEVEL:-HyMy}

CELLTYPE_INDEX=${SLURM_ARRAY_TASK_ID}

# Map index to name for logging
CELLTYPE_NAMES=("" "LEC" "FRC" "BEC" "CD4_T_cells" "CD8_T_cells" "gdT_cells" "Macrophages" "Monocyte" "Fibroblasts_mac" "cDC1" "cDC2" "mature_migDC" "B_cells" "Plasma_cell")
CELLTYPE_NAME=${CELLTYPE_NAMES[$CELLTYPE_INDEX]}

echo "=============================================="
echo "HyMy Distance Correlation - Array Job"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Cell Type: ${CELLTYPE_NAME} (index ${CELLTYPE_INDEX})"
echo "Annotation Level: ${ANNOTATION_LEVEL}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: 128G"
echo "Date: $(date)"
echo "=============================================="

# Change to script directory
cd /nobackup/lab_maier/Projects/mXenium/CMM/scripts/workflow/scripts/one_off/spatial_analysis

# Create logs directory if needed
mkdir -p /nobackup/lab_maier/Projects/mXenium/CMM/results/logs

# Load conda environment
source /home/cmangana/miniconda3/etc/profile.d/conda.sh
conda activate R_IMC_2024

# Export environment variables for R script
export ANNOTATION_LEVEL
export CELLTYPE_INDEX
export SLURM_CPUS_PER_TASK
export N_PERMUTATIONS=${N_PERMUTATIONS:-0}

# Run analysis for this cell type
echo ""
echo "Running hymy_distance_correlation.R for ${CELLTYPE_NAME}..."
Rscript hymy_distance_correlation.R

echo ""
echo "=============================================="
echo "Cell type ${CELLTYPE_NAME} complete: $(date)"
echo "=============================================="
