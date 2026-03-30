#!/bin/bash
#SBATCH --job-name=ripple_perm_gpu
#SBATCH --output=logs/ripple_perm_gpu_%A_%a.out
#SBATCH --error=logs/ripple_perm_gpu_%A_%a.err
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=64g
## Set --array via sbatch: sbatch --array=1-N run_permutation_gpu.sh
## Partition/QOS: uncomment and adjust for your cluster
##SBATCH --partition=gpu
##SBATCH --qos=gpu
## Mail notifications (uncomment and set your email)
##SBATCH --mail-type=END
##SBATCH --mail-user=your@email.com

# =============================================================================
# RIPPLE Stage 2: GPU-Accelerated Permutation Testing - SLURM Array Job
# =============================================================================
# Runs the GPU permutation script for each cell type in parallel.
#
# Prerequisites:
#   1. Run R script with N_PERMUTATIONS=0 first (produces meta_analysis_results.csv)
#   2. Wait for all R jobs to complete
#
# Usage:
#   sbatch --array=1-N run_permutation_gpu.sh
#   N_PERMUTATIONS=1000 sbatch --array=1-N run_permutation_gpu.sh
#   sbatch --array=1 run_permutation_gpu.sh  # single cell type
#
# After completion:
#   Rscript merge_permutation_pvals.R
#
# Expected Runtime: ~10-30 minutes per cell type on GPU
# =============================================================================

set -e

# Auto-detect script directory
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Create logs directory
mkdir -p logs

# Pass through RIPPLE env vars
export INPUT_PATH="${INPUT_PATH:-}"
export OUTPUT_DIR="${OUTPUT_DIR:-}"
export QUERY_CELLTYPE="${QUERY_CELLTYPE:-}"
export CELLTYPE_COLUMN="${CELLTYPE_COLUMN:-}"
export SAMPLE_COLUMN="${SAMPLE_COLUMN:-}"
export CONDITION_COLUMN="${CONDITION_COLUMN:-}"
export CONDITION_VALUE="${CONDITION_VALUE:-}"
export X_COLUMN="${X_COLUMN:-}"
export Y_COLUMN="${Y_COLUMN:-}"
export TARGET_CELLTYPES="${TARGET_CELLTYPES:-}"
export CONTROL_CELLTYPE="${CONTROL_CELLTYPE:-}"
export QUERY_LABEL="${QUERY_LABEL:-}"
export ANNOTATION_LEVEL="${ANNOTATION_LEVEL:-}"
export ANALYSIS_NAME="${ANALYSIS_NAME:-}"
export QUERY_SIGNATURE_GENES="${QUERY_SIGNATURE_GENES:-}"
export ADATA_PATH="${ADATA_PATH:-}"

export CELLTYPE_INDEX=${SLURM_ARRAY_TASK_ID}
export N_PERMUTATIONS=${N_PERMUTATIONS:-500}
export K_NEIGHBORS=${K_NEIGHBORS:-1}

echo "=============================================="
echo "RIPPLE Stage 2: GPU Permutation Testing - Array Job"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Cell Type Index: ${CELLTYPE_INDEX}"
echo "Analysis Name: ${ANALYSIS_NAME:-hymy_distance_correlation}"
echo "K Neighbors: ${K_NEIGHBORS}"
echo "N Permutations: ${N_PERMUTATIONS}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Date: $(date)"
echo "=============================================="

# Load CUDA module (conditional)
if [ -n "${CUDA_MODULE:-}" ]; then module load "$CUDA_MODULE"; fi

# Conda environment setup (conditional)
if [ -n "${CONDA_SETUP:-}" ]; then source "$CONDA_SETUP"; fi
if [ -n "${CONDA_ENV:-}" ]; then conda activate "$CONDA_ENV"; fi

# Verify GPU availability
echo ""
echo "GPU Check:"
python -c "import torch; print(f'  PyTorch {torch.__version__}'); print(f'  CUDA available: {torch.cuda.is_available()}'); print(f'  Device: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else \"N/A\"}')"

echo ""
echo "Running run_permutation_gpu.py..."
python run_permutation_gpu.py --n-perms ${N_PERMUTATIONS}

echo ""
echo "=============================================="
echo "Complete: $(date)"
echo "=============================================="
