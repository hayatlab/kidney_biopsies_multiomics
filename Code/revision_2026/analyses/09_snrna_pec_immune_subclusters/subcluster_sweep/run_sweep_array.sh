#!/bin/bash
#SBATCH --job-name=l2_subcluster
#SBATCH --array=0-15%5
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH --output=logs/sweep_%A_%a.out
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

# Subcluster each of the 16 level2 cell types with the PEC-style sweep (crescent machinery removed,
# per-lineage canonical panel added). One array task per cell type. Env: pathfire (classic
# harmonypy 0.0.10 -> scanpy.external harmony works; igraph leiden; umap).

TYPES=(PT TAL EC aPT PC DCT Immune CNT dPT IC-A FIB POD PEC IC-B VSM-P MES)
CT="${TYPES[$SLURM_ARRAY_TASK_ID]}"

DIR=${KB_ANALYSIS}
mkdir -p "$DIR/subcluster_level2_20260728/logs"
cd "$DIR"

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate pathfire
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

set -e   # fail loudly on python error -- placed AFTER conda activate (its scripts are NOT set -u safe)
echo "[$(date)] task $SLURM_ARRAY_TASK_ID -> celltype=$CT  python=$(which python)"
python -u ${KB_REPO}/analyses/09_snrna_pec_immune_subclusters/subcluster_sweep/subcluster_sweep.py --celltype "$CT"
echo "[$(date)] done $CT"
