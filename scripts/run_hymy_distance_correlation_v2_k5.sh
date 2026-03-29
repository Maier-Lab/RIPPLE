#!/bin/bash
#SBATCH --job-name=hymy_dist_v2_k5
#SBATCH --output=/nobackup/lab_maier/Projects/mXenium/CMM/results/logs/hymy_dist_v2_k5_%A_%a.out
#SBATCH --error=/nobackup/lab_maier/Projects/mXenium/CMM/results/logs/hymy_dist_v2_k5_%A_%a.err
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
# HyMy Distance Correlation v2 - k=5 Sanity Check
# =============================================================================
# Same Poisson GLM as main v2, but averages distance to 5 nearest query cells
# instead of single nearest neighbor. This tests robustness of results to the
# distance metric — high overlap with k=1 results confirms stability.
#
# Usage:
#   sbatch run_hymy_distance_correlation_v2_k5.sh
#   ANNOTATION_LEVEL=L1 sbatch run_hymy_distance_correlation_v2_k5.sh
# =============================================================================

set -e

ANNOTATION_LEVEL=${ANNOTATION_LEVEL:-HyMy}
CELLTYPE_INDEX=${SLURM_ARRAY_TASK_ID}

CELLTYPE_NAMES=("" "LEC" "FRC" "BEC" "CD4_T_cells" "CD8_T_cells" "gdT_cells" "Macrophages" "Monocyte" "Fibroblasts_mac" "cDC1" "cDC2" "mature_migDC" "B_cells" "Plasma_cell")
CELLTYPE_NAME=${CELLTYPE_NAMES[$CELLTYPE_INDEX]}

echo "=============================================="
echo "HyMy Distance Correlation v2 - k=5 Sanity Check"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Cell Type: ${CELLTYPE_NAME} (index ${CELLTYPE_INDEX})"
echo "Annotation Level: ${ANNOTATION_LEVEL}"
echo "K Neighbors: 5 (sanity check)"
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
export K_NEIGHBORS=5

echo ""
echo "Running hymy_distance_correlation_v2.R for ${CELLTYPE_NAME} (k=5)..."
Rscript hymy_distance_correlation_v2.R

echo ""
echo "=============================================="
echo "Cell type ${CELLTYPE_NAME} (k=5) complete: $(date)"
echo "=============================================="
