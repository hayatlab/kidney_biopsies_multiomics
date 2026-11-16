#!/bin/bash
#SBATCH --job-name=vt26_nfkb
#SBATCH --output=logs/vt26_nfkb.%j.out
#SBATCH --error=logs/vt26_nfkb.%j.out
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

# no `set -u` -- breaks conda activate scripts (ADDR2LINE: unbound variable)
# Networks must already be cached: run `python scripts/26_score_nfkb_xenium.py --fetch-nets`
# once on the login node (compute nodes may have no internet).
cd ${KB_ANALYSIS}/vivien_takeover
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate nichecompass
set -e
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
python ${KB_REPO}/analyses/03_xenium_factor1_neighbourhoods_cci/26_score_nfkb_xenium.py "$@"
