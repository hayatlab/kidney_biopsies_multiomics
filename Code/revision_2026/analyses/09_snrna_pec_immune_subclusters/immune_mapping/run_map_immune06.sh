#!/bin/bash
### scANVI (scArches-style) transfer of snRNA immune leiden_0.6 subclusters -> Xenium immune cells
#SBATCH --job-name=map_immune06
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --output=logs/map_immune06.%J.out
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate kpmp_repurposing
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd ${KB_ANALYSIS}/immune_mapping_leiden06
set -e   # placed AFTER conda activate (its scripts are not set -u safe)
echo "[$(date)] start on $(hostname); python=$(which python)"
python -u ${KB_REPO}/analyses/09_snrna_pec_immune_subclusters/immune_mapping/map_immune_leiden06_xenium.py
echo "[$(date)] done"
