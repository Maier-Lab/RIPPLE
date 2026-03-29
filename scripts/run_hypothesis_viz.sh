#!/bin/bash
#SBATCH --job-name=hymy_viz
#SBATCH --output=/nobackup/lab_maier/Projects/mXenium/CMM/results/logs/hymy_viz_%j.out
#SBATCH --error=/nobackup/lab_maier/Projects/mXenium/CMM/results/logs/hymy_viz_%j.err
#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=32g
#SBATCH --mail-type=END
#SBATCH --mail-user=cmangana@cemm.at

# =============================================================================
# Hypothesis Visualizations — Decay Plots
# =============================================================================
# Generates publication-quality decay curves for key biological themes:
#   1. CD8 exhaustion vs stem-like
#   2. CD4 Foxp3 (Treg enrichment)
#   3. LEC remodeling
#   4. FRC inflammatory remodeling
#
# Usage:
#   sbatch run_hypothesis_viz.sh
#   ANNOTATION_LEVEL=L1 sbatch run_hypothesis_viz.sh
# =============================================================================

set -e

# Query cell type configuration (pass through to R/Python)
export QUERY_CELLTYPE=${QUERY_CELLTYPE:-}
export CELLTYPE_COLUMN=${CELLTYPE_COLUMN:-}
export QUERY_LABEL=${QUERY_LABEL:-}
export ANNOTATION_LEVEL=${ANNOTATION_LEVEL:-HyMy}

echo "=============================================="
echo "Hypothesis Visualizations"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Annotation Level: ${ANNOTATION_LEVEL}"
echo "Date: $(date)"
echo "=============================================="

cd /nobackup/lab_maier/Projects/mXenium/CMM/scripts/workflow/scripts/one_off/spatial_analysis

mkdir -p /nobackup/lab_maier/Projects/mXenium/CMM/results/logs

source /home/cmangana/miniconda3/etc/profile.d/conda.sh
conda activate R_IMC_2024

export ANNOTATION_LEVEL

Rscript hypothesis_visualizations.R

echo ""
echo "=============================================="
echo "Complete: $(date)"
echo "=============================================="
