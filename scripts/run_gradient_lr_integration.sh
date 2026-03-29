#!/bin/bash
#SBATCH --job-name=ripple_lr
#SBATCH --output=logs/ripple_lr_%A_%a.out
#SBATCH --error=logs/ripple_lr_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --mem=64g
## Set --array via sbatch: sbatch --array=1-N run_gradient_lr_integration.sh
## Partition/QOS: uncomment and adjust for your cluster
##SBATCH --partition=tinyq
##SBATCH --qos=tinyq
## Mail notifications (uncomment and set your email)
##SBATCH --mail-type=END
##SBATCH --mail-user=your@email.com

# ==============================================================================
# RIPPLE Stage 6: Gradient-to-Ligand-Receptor Integration Analysis
# ==============================================================================
#
# SLURM array job: one job per cell type.
#
# Prerequisites:
#   - Stage 1 distance correlation complete (all_genes_results.csv)
#   - Stage 2 classification complete (stage2_all_results.csv) [recommended]
#
# Usage:
#   sbatch --array=1-N run_gradient_lr_integration.sh
#   GRADIENT_SOURCE=hymy_distance_correlation_v2 sbatch --array=1-N ...
#   sbatch --array=1 run_gradient_lr_integration.sh  # single cell type
# ==============================================================================

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
export GRADIENT_SOURCE=${GRADIENT_SOURCE:-hymy_distance_correlation}

echo "=============================================="
echo "RIPPLE Stage 6: Gradient-to-LR Integration"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Cell Type Index: ${CELLTYPE_INDEX}"
echo "Gradient Source: ${GRADIENT_SOURCE}"
echo "Date: $(date)"
echo "=============================================="

# Conda environment setup (conditional)
if [ -n "${CONDA_SETUP:-}" ]; then source "$CONDA_SETUP"; fi
if [ -n "${CONDA_ENV:-}" ]; then conda activate "$CONDA_ENV"; fi

Rscript gradient_lr_integration.R

echo ""
echo "=============================================="
echo "Complete: $(date)"
echo "=============================================="
