#!/bin/bash
#SBATCH --job-name=leiden_sweep
#SBATCH --output=logs/leiden_sweep.%j.out
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=4:00:00
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"
set -euo pipefail
cd ${KB_PROJECT}
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate kpmp_repurposing
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16
cd review_analyses/celltype_concordance
python 03_leiden_sweep.py
