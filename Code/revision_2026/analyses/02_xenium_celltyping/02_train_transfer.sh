#!/bin/bash
#SBATCH --job-name=scanvi_transfer
#SBATCH --output=logs/scanvi_transfer.%j.out
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=12:00:00
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"
set -euo pipefail
cd ${KB_PROJECT}
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate kpmp_repurposing
export OMP_NUM_THREADS=32 MKL_NUM_THREADS=32
cd review_analyses/celltype_concordance
python 02_train_transfer.py
