#!/bin/bash
#SBATCH --job-name=vt09_ccipost
#SBATCH --output=logs/vt09_ccipost.%j.out
#SBATCH --error=logs/vt09_ccipost.%j.out
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

# no `set -u` -- breaks conda activate scripts (ADDR2LINE: unbound variable)
cd ${KB_ANALYSIS}/vivien_takeover
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate cci_env
set -e
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
python ${KB_REPO}/analyses/03_xenium_factor1_neighbourhoods_cci/09_cci_postprocess.py "$@"
