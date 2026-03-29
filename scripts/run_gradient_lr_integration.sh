#!/bin/bash
#SBATCH --job-name=grad_lr
#SBATCH --output=/nobackup/lab_maier/Projects/mXenium/CMM/results/logs/grad_lr_%A_%a.out
#SBATCH --error=/nobackup/lab_maier/Projects/mXenium/CMM/results/logs/grad_lr_%A_%a.err
#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --mem=64g
#SBATCH --array=1-14
#SBATCH --mail-type=END
#SBATCH --mail-user=cmangana@cemm.at

# ==============================================================================
# Gradient-to-Ligand-Receptor Integration Analysis
# ==============================================================================
#
# SLURM array job: one job per cell type (14 total).
#
# Cell type mapping (same as hymy_distance_correlation.R):
#   1=LEC, 2=FRC, 3=BEC, 4=CD4_T_cells, 5=CD8_T_cells, 6=gdT_cells,
#   7=Macrophages, 8=Monocyte, 9=Fibroblasts_mac, 10=cDC1, 11=cDC2,
#   12=mature_migDC, 13=B_cells, 14=Plasma_cell
#
# Prerequisites:
#   - Stage 1 distance correlation complete (all_genes_results.csv)
#   - Stage 2 classification complete (stage2_all_results.csv) [recommended]
#
# Usage:
#   sbatch run_gradient_lr_integration.sh                  # v1 (logistic, default)
#   GRADIENT_SOURCE=hymy_distance_correlation_v2 sbatch run_gradient_lr_integration.sh  # v2 (Poisson)
#   ANNOTATION_LEVEL=L1 sbatch run_gradient_lr_integration.sh
#
# Test single cell type (e.g., LEC for cross-validation with NicheNet):
#   sbatch --array=1 run_gradient_lr_integration.sh
#
# Author: CMM Project
# Date: 2026-02
# ==============================================================================

set -e

echo "=============================================="
echo "Gradient-to-LR Integration Analysis"
echo "=============================================="
echo "Start time: $(date)"
echo "Job ID: ${SLURM_JOB_ID}, Array Task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo ""

# Environment setup
source /home/cmangana/miniconda3/etc/profile.d/conda.sh
conda activate /nobackup/lab_maier/envs/nichenet_env

# Configuration
# Query cell type configuration (pass through to R/Python)
export QUERY_CELLTYPE=${QUERY_CELLTYPE:-}
export CELLTYPE_COLUMN=${CELLTYPE_COLUMN:-}
export QUERY_LABEL=${QUERY_LABEL:-}
export ANNOTATION_LEVEL=${ANNOTATION_LEVEL:-HyMy}
export CELLTYPE_INDEX=${SLURM_ARRAY_TASK_ID}
export GRADIENT_SOURCE=${GRADIENT_SOURCE:-hymy_distance_correlation}

echo "Configuration:"
echo "  ANNOTATION_LEVEL: ${ANNOTATION_LEVEL}"
echo "  GRADIENT_SOURCE: ${GRADIENT_SOURCE}"
echo "  CELLTYPE_INDEX: ${CELLTYPE_INDEX}"
echo ""

# Run analysis
SCRIPT_DIR="/nobackup/lab_maier/Projects/mXenium/CMM/scripts/workflow/scripts/one_off/spatial_analysis"
cd "${SCRIPT_DIR}"

Rscript gradient_lr_integration.R

echo ""
echo "Completed at: $(date)"
echo "=============================================="
