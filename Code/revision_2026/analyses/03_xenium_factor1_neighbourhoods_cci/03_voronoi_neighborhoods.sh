#!/bin/bash
#SBATCH --job-name=vt03_nh
#SBATCH --output=logs/vt03_nh.%j.out
#SBATCH --error=logs/vt03_nh.%j.out
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

# no `set -u` -- breaks conda activate scripts (ADDR2LINE: unbound variable)
cd ${KB_ANALYSIS}/vivien_takeover
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate spatialdata
set -e
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
python ${KB_REPO}/analyses/03_xenium_factor1_neighbourhoods_cci/03_voronoi_neighborhoods.py
