#!/bin/bash
#SBATCH --job-name=n20val_pt
#SBATCH --output=niche20_validation_pseudotime.%j.out
#SBATCH --error=niche20_validation_pseudotime.%j.out
#SBATCH --time=06:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=8
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

set -euo pipefail
cd ${KB_PROJECT}
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate spatial_env

python -u ${KB_REPO}/analyses/08_glomerular_niche20_crescent/niche20_validation_pseudotime.py
