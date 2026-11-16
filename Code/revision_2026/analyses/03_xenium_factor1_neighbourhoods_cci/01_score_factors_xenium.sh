#!/bin/bash
#SBATCH --job-name=vt01_factors
#SBATCH --output=logs/vt01_factors.%j.out
#SBATCH --error=logs/vt01_factors.%j.out
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

# NB: deliberately no `set -u` -- it breaks conda's activate scripts
#     (ADDR2LINE: unbound variable). `set -e` only after activation.
cd ${KB_ANALYSIS}/vivien_takeover
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate nichecompass
set -e
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
python ${KB_REPO}/analyses/03_xenium_factor1_neighbourhoods_cci/01_score_factors_xenium.py
