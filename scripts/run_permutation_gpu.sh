#!/bin/bash
#SBATCH --job-name=hymy_perm_gpu
#SBATCH --output=/nobackup/lab_maier/Projects/mXenium/CMM/results/logs/hymy_perm_gpu_%A_%a.out
#SBATCH --error=/nobackup/lab_maier/Projects/mXenium/CMM/results/logs/hymy_perm_gpu_%A_%a.err
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:l4_gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=64g
#SBATCH --array=1-14
#SBATCH --mail-type=END
#SBATCH --mail-user=cmangana@cemm.at

# =============================================================================
# GPU-Accelerated Permutation Testing - SLURM Array Job
# =============================================================================
# Runs the GPU permutation script for each cell type in parallel.
#
# Prerequisites:
#   1. Run R script with N_PERMUTATIONS=0 first (produces meta_analysis_results.csv)
#      N_PERMUTATIONS=0 sbatch run_hymy_distance_correlation_array.sh
#   2. Wait for all R jobs to complete
#
# Cell type mapping (same as R scripts):
#   1 = LEC, 2 = FRC, 3 = BEC, 4 = CD4_T_cells, 5 = CD8_T_cells,
#   6 = gdT_cells, 7 = Macrophages, 8 = Monocyte, 9 = Fibroblasts_mac,
#   10 = cDC1, 11 = cDC2, 12 = mature_migDC, 13 = B_cells, 14 = Plasma_cell
#
# Usage:
#   # Run all 14 cell types in parallel (HyMy annotation, default v1)
#   sbatch run_permutation_gpu.sh
#
#   # Run for v2 Poisson GLM (k=1)
#   ANALYSIS_NAME=hymy_distance_correlation_v2 sbatch run_permutation_gpu.sh
#
#   # Run for v2 Poisson GLM (k=5)
#   ANALYSIS_NAME=hymy_distance_correlation_v2_k5 K_NEIGHBORS=5 sbatch run_permutation_gpu.sh
#
#   # Run with L1 annotation
#   ANNOTATION_LEVEL=L1 sbatch run_permutation_gpu.sh
#
#   # Run a single cell type (e.g., LEC only)
#   sbatch --array=1 run_permutation_gpu.sh
#
#   # Custom permutation count (default: 500)
#   N_PERMUTATIONS=1000 sbatch run_permutation_gpu.sh
#
# After completion, merge results:
#   Rscript merge_permutation_pvals.R
#
# Expected Runtime: ~10-30 minutes per cell type on L4 GPU
# =============================================================================

set -e

# Query cell type configuration (pass through to R/Python)
export QUERY_CELLTYPE=${QUERY_CELLTYPE:-}
export CELLTYPE_COLUMN=${CELLTYPE_COLUMN:-}
export QUERY_LABEL=${QUERY_LABEL:-}
export ANNOTATION_LEVEL=${ANNOTATION_LEVEL:-HyMy}
N_PERMUTATIONS=${N_PERMUTATIONS:-500}
ANALYSIS_NAME=${ANALYSIS_NAME:-hymy_distance_correlation}
K_NEIGHBORS=${K_NEIGHBORS:-1}
CELLTYPE_INDEX=${SLURM_ARRAY_TASK_ID}

# Map index to name for logging
CELLTYPE_NAMES=("" "LEC" "FRC" "BEC" "CD4_T_cells" "CD8_T_cells" "gdT_cells" "Macrophages" "Monocyte" "Fibroblasts_mac" "cDC1" "cDC2" "mature_migDC" "B_cells" "Plasma_cell")
CELLTYPE_NAME=${CELLTYPE_NAMES[$CELLTYPE_INDEX]}

echo "=============================================="
echo "GPU Permutation Testing - Array Job"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Cell Type: ${CELLTYPE_NAME} (index ${CELLTYPE_INDEX})"
echo "Analysis Name: ${ANALYSIS_NAME}"
echo "Annotation Level: ${ANNOTATION_LEVEL}"
echo "K Neighbors: ${K_NEIGHBORS}"
echo "N Permutations: ${N_PERMUTATIONS}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: 64G"
echo "GPU: L4 (24GB)"
echo "Date: $(date)"
echo "=============================================="

# Change to script directory
cd /nobackup/lab_maier/Projects/mXenium/CMM/scripts/workflow/scripts/one_off/spatial_analysis

# Create logs directory if needed
mkdir -p /nobackup/lab_maier/Projects/mXenium/CMM/results/logs

# Load CUDA toolkit
module purge
module load cuda12.2/toolkit/12.2.2

# Load conda environment with PyTorch
source /home/cmangana/miniconda3/etc/profile.d/conda.sh
conda activate "/nobackup/lab_maier/envs/harpy"

# Verify GPU availability
echo ""
echo "GPU Check:"
python -c "import torch; print(f'  PyTorch {torch.__version__}'); print(f'  CUDA available: {torch.cuda.is_available()}'); print(f'  Device: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else \"N/A\"}')"

# Export environment variables for Python script
export ANNOTATION_LEVEL
export CELLTYPE_INDEX
export ANALYSIS_NAME
export K_NEIGHBORS

echo ""
echo "Running run_permutation_gpu.py for ${CELLTYPE_NAME}..."
python run_permutation_gpu.py --n-perms ${N_PERMUTATIONS}

echo ""
echo "=============================================="
echo "Cell type ${CELLTYPE_NAME} complete: $(date)"
echo "=============================================="
