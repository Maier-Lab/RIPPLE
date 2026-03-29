#!/bin/bash
#SBATCH --job-name=hymy_dist_s2
#SBATCH --output=/nobackup/lab_maier/Projects/mXenium/CMM/results/logs/hymy_dist_s2_%A_%a.out
#SBATCH --error=/nobackup/lab_maier/Projects/mXenium/CMM/results/logs/hymy_dist_s2_%A_%a.err
#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mem=64g
#SBATCH --array=1-14
#SBATCH --mail-type=END
#SBATCH --mail-user=cmangana@cemm.at

# =============================================================================
# HyMy Distance Correlation Stage 2 - SLURM Array Job
# =============================================================================
# Validates Stage 1 hits by adding monocyte distance as a covariate.
# No permutation testing → faster than Stage 1.
#
# Requires: Stage 1 results merged (run merge_distance_correlation_results.R)
#
# Cell type mapping:
#   1 = LEC, 2 = FRC, 3 = BEC, 4 = CD4_T_cells, 5 = CD8_T_cells,
#   6 = gdT_cells, 7 = Macrophages, 8 = Monocyte*, 9 = Fibroblasts_mac,
#   10 = cDC1, 11 = cDC2, 12 = mature_migDC, 13 = B_cells, 14 = Plasma_cell
#
#   * Monocyte uses Macrophages as alternative control (cannot control for self)
#
# Usage:
#   sbatch run_hymy_distance_correlation_stage2_array.sh                          # v2 (default)
#   ANALYSIS_NAME=hymy_distance_correlation sbatch run_hymy_distance_correlation_stage2_array.sh  # v1
#   ANNOTATION_LEVEL=L1 sbatch run_hymy_distance_correlation_stage2_array.sh      # L1 annotation
#
# =============================================================================

set -e

# Query cell type configuration (pass through to R/Python)
export QUERY_CELLTYPE=${QUERY_CELLTYPE:-}
export CELLTYPE_COLUMN=${CELLTYPE_COLUMN:-}
export QUERY_LABEL=${QUERY_LABEL:-}
export ANNOTATION_LEVEL=${ANNOTATION_LEVEL:-HyMy}
ANALYSIS_NAME=${ANALYSIS_NAME:-hymy_distance_correlation_v2}
CELLTYPE_INDEX=${SLURM_ARRAY_TASK_ID}

CELLTYPE_NAMES=("" "LEC" "FRC" "BEC" "CD4_T_cells" "CD8_T_cells" "gdT_cells" "Macrophages" "Monocyte" "Fibroblasts_mac" "cDC1" "cDC2" "mature_migDC" "B_cells" "Plasma_cell")
CELLTYPE_NAME=${CELLTYPE_NAMES[$CELLTYPE_INDEX]}

echo "=============================================="
echo "HyMy Distance Correlation Stage 2 - Array Job"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Cell Type: ${CELLTYPE_NAME} (index ${CELLTYPE_INDEX})"
echo "Annotation Level: ${ANNOTATION_LEVEL}"
echo "Stage 1 Analysis: ${ANALYSIS_NAME}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: 64G"
echo "Date: $(date)"
echo "=============================================="

cd /nobackup/lab_maier/Projects/mXenium/CMM/scripts/workflow/scripts/one_off/spatial_analysis

mkdir -p /nobackup/lab_maier/Projects/mXenium/CMM/results/logs

source /home/cmangana/miniconda3/etc/profile.d/conda.sh
conda activate R_IMC_2024

export ANNOTATION_LEVEL
export ANALYSIS_NAME
export CELLTYPE_INDEX

echo ""
echo "Running hymy_distance_correlation_stage2.R for ${CELLTYPE_NAME}..."
Rscript hymy_distance_correlation_stage2.R

echo ""
echo "=============================================="
echo "Cell type ${CELLTYPE_NAME} complete: $(date)"
echo "=============================================="
