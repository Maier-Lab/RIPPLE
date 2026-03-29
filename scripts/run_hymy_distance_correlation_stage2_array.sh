#!/bin/bash
#SBATCH --job-name=ripple_stage2
#SBATCH --output=logs/ripple_stage2_%A_%a.out
#SBATCH --error=logs/ripple_stage2_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mem=64g
## Set --array via sbatch: sbatch --array=1-N run_hymy_distance_correlation_stage2_array.sh
## Partition/QOS: uncomment and adjust for your cluster
##SBATCH --partition=shortq
##SBATCH --qos=shortq
## Mail notifications (uncomment and set your email)
##SBATCH --mail-type=END
##SBATCH --mail-user=your@email.com

# =============================================================================
# RIPPLE Stage 4: Distance Correlation Stage 2 - SLURM Array Job
# =============================================================================
# Validates Stage 1 hits by adding control cell type distance as a covariate.
# No permutation testing -- faster than Stage 1.
#
# Requires: Stage 1 results merged (run merge_distance_correlation_results.R)
#
# Usage:
#   sbatch --array=1-N run_hymy_distance_correlation_stage2_array.sh
#   ANALYSIS_NAME=hymy_distance_correlation sbatch --array=1-N ...  # v1 source
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

echo "=============================================="
echo "RIPPLE Stage 4: Distance Correlation Stage 2 - Array Job"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Cell Type Index: ${CELLTYPE_INDEX}"
echo "Analysis Name: ${ANALYSIS_NAME:-hymy_distance_correlation_v2}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Date: $(date)"
echo "=============================================="

# Conda environment setup (conditional)
if [ -n "${CONDA_SETUP:-}" ]; then source "$CONDA_SETUP"; fi
if [ -n "${CONDA_ENV:-}" ]; then conda activate "$CONDA_ENV"; fi

echo ""
echo "Running hymy_distance_correlation_stage2.R..."
Rscript hymy_distance_correlation_stage2.R

echo ""
echo "=============================================="
echo "Complete: $(date)"
echo "=============================================="
