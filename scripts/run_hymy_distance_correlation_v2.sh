#!/bin/bash
#SBATCH --job-name=hymy_dist_v2
#SBATCH --output=/nobackup/lab_maier/Projects/mXenium/CMM/results/logs/hymy_dist_v2_%A_%a.out
#SBATCH --error=/nobackup/lab_maier/Projects/mXenium/CMM/results/logs/hymy_dist_v2_%A_%a.err
#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00
#SBATCH --mem=64g
#SBATCH --array=1-14
#SBATCH --mail-type=END
#SBATCH --mail-user=cmangana@cemm.at

# =============================================================================
# HyMy Distance Correlation v2 (Poisson GLM) - SLURM Array Job
# =============================================================================
# Uses Poisson GLM: counts ~ distance + offset(log(total_counts))
# Coefficient = log-rate change per um (negative = HyMy-induced)
#
# Cell type mapping:
#   1 = LEC, 2 = FRC, 3 = BEC, 4 = CD4_T_cells, 5 = CD8_T_cells,
#   6 = gdT_cells, 7 = Macrophages, 8 = Monocyte, 9 = Fibroblasts_mac,
#   10 = cDC1, 11 = cDC2, 12 = mature_migDC, 13 = B_cells, 14 = Plasma_cell
#
# Usage:
#   sbatch run_hymy_distance_correlation_v2.sh
#   ANNOTATION_LEVEL=L1 sbatch run_hymy_distance_correlation_v2.sh
#   sbatch --array=1-2 run_hymy_distance_correlation_v2.sh  # LEC+FRC only
# =============================================================================

set -e

ANNOTATION_LEVEL=${ANNOTATION_LEVEL:-HyMy}
CELLTYPE_INDEX=${SLURM_ARRAY_TASK_ID}

CELLTYPE_NAMES=("" "LEC" "FRC" "BEC" "CD4_T_cells" "CD8_T_cells" "gdT_cells" "Macrophages" "Monocyte" "Fibroblasts_mac" "cDC1" "cDC2" "mature_migDC" "B_cells" "Plasma_cell")
CELLTYPE_NAME=${CELLTYPE_NAMES[$CELLTYPE_INDEX]}

echo "=============================================="
echo "HyMy Distance Correlation v2 (Poisson GLM)"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Cell Type: ${CELLTYPE_NAME} (index ${CELLTYPE_INDEX})"
echo "Annotation Level: ${ANNOTATION_LEVEL}"
echo "K Neighbors: 1 (default)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Date: $(date)"
echo "=============================================="

cd /nobackup/lab_maier/Projects/mXenium/CMM/scripts/workflow/scripts/one_off/spatial_analysis

mkdir -p /nobackup/lab_maier/Projects/mXenium/CMM/results/logs

source /home/cmangana/miniconda3/etc/profile.d/conda.sh
conda activate R_IMC_2024

export ANNOTATION_LEVEL
export CELLTYPE_INDEX
export SLURM_CPUS_PER_TASK
export N_PERMUTATIONS=${N_PERMUTATIONS:-0}
export K_NEIGHBORS=1

echo ""
echo "Running hymy_distance_correlation_v2.R for ${CELLTYPE_NAME}..."
Rscript hymy_distance_correlation_v2.R

echo ""
echo "=============================================="
echo "Cell type ${CELLTYPE_NAME} complete: $(date)"
echo "=============================================="
