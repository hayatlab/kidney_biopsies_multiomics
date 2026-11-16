#!/bin/bash
#SBATCH --job-name=na01_egfr
#SBATCH --output=logs/na01_egfr.%j.out
#SBATCH --error=logs/na01_egfr.%j.out
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

# no `set -u` -- breaks conda activate scripts (ADDR2LINE: unbound variable)
cd ${KB_ANALYSIS}/nell1_app_egfr
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate nichecompass
set -e
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
python ${KB_REPO}/analyses/05_nell1_snrna/01_nell1_app_egfr.py "$@"
