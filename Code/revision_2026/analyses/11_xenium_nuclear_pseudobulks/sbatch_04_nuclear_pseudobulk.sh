#!/bin/bash
#SBATCH --job-name=xz_04_nuclear_pseudobulk --partition=defq --cpus-per-task=8 --mem=96G
#SBATCH --output=logs/04_nuclear_pseudobulk.%j.out
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"
set -euo pipefail
cd ${KB_ANALYSIS}/xenium_zenodo_20261027
unset ${!CONDA_BACKUP_@} 2>/dev/null || true
source ${CONDA_ROOT}/etc/profile.d/conda.sh; conda activate nichecompass
export OMP_NUM_THREADS=8
echo "env hash nichecompass: $(conda list --explicit | md5sum | cut -c1-12)  $(date)"
python ${KB_REPO}/analyses/11_xenium_nuclear_pseudobulks/04_nuclear_pseudobulk.py
