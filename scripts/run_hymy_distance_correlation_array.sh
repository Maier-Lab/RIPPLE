#!/bin/bash
#SBATCH --job-name=ripple_dist_corr
#SBATCH --output=logs/ripple_dist_corr_%A_%a.out
#SBATCH --error=logs/ripple_dist_corr_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00
#SBATCH --mem=32g
## Set --array via sbatch: sbatch --array=1-N run_hymy_distance_correlation_array.sh
## Partition/QOS: uncomment and adjust for your cluster
##SBATCH --partition=tinyq
##SBATCH --qos=tinyq
## Mail notifications (uncomment and set your email)
##SBATCH --mail-type=END
##SBATCH --mail-user=your@email.com

# =============================================================================
# RIPPLE: Distance Correlation (Logistic Regression) - SLURM Array Job
# =============================================================================
# Uses logistic regression: P(expressing) ~ distance to query
# Gradient score = log-odds coefficient (negative = query-induced)
#
# Runs each cell type as a separate job for parallelization.
#
# Usage:
#   sbatch --array=1-N run_hymy_distance_correlation_array.sh
#   sbatch --array=2 run_hymy_distance_correlation_array.sh  # single cell type
#
# Required env vars: INPUT_PATH, QUERY_CELLTYPE, CELLTYPE_COLUMN
# See README.md for full env var reference.
#
# After completion, merge results:
#   Rscript merge_distance_correlation_results.R
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
export N_PERMUTATIONS=${N_PERMUTATIONS:-0}

echo "=============================================="
echo "RIPPLE: Distance Correlation - Array Job"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Cell Type Index: ${CELLTYPE_INDEX}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Date: $(date)"
echo "=============================================="

# Conda environment setup (conditional)
if [ -n "${CONDA_SETUP:-}" ]; then source "$CONDA_SETUP"; fi
if [ -n "${CONDA_ENV:-}" ]; then conda activate "$CONDA_ENV"; fi

echo ""
echo "Running hymy_distance_correlation.R..."
Rscript hymy_distance_correlation.R

echo ""
echo "=============================================="
echo "Complete: $(date)"
echo "=============================================="
