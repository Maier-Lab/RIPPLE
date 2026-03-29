#!/bin/bash
#SBATCH --job-name=ripple_viz
#SBATCH --output=logs/ripple_viz_%j.out
#SBATCH --error=logs/ripple_viz_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=32g
## Partition/QOS: uncomment and adjust for your cluster
##SBATCH --partition=tinyq
##SBATCH --qos=tinyq
## Mail notifications (uncomment and set your email)
##SBATCH --mail-type=END
##SBATCH --mail-user=your@email.com

# =============================================================================
# RIPPLE Stage 5: Hypothesis Visualizations -- Decay Plots
# =============================================================================
# Generates publication-quality decay curves for biological themes.
#
# Usage:
#   sbatch run_hypothesis_viz.sh
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

echo "=============================================="
echo "RIPPLE Stage 5: Hypothesis Visualizations"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Date: $(date)"
echo "=============================================="

# Conda environment setup (conditional)
if [ -n "${CONDA_SETUP:-}" ]; then source "$CONDA_SETUP"; fi
if [ -n "${CONDA_ENV:-}" ]; then conda activate "$CONDA_ENV"; fi

Rscript hypothesis_visualizations.R

echo ""
echo "=============================================="
echo "Complete: $(date)"
echo "=============================================="
