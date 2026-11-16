#!/bin/bash

### Job Parameters
#SBATCH --time=12:00:00
#SBATCH --job-name=niche_5reg
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1
#SBATCH --mem=256GB

#SBATCH --output=logs/niche_5regions.%J.out

source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate nichecompass

echo "Environment activated: $(which python)"
python --version
echo "Started $(date)"

cd ${KB_PROJECT}

python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/nichecompass_5regions.py \
    --leiden-res 0.4 \
    --n-epochs 400

echo "Completed $(date)"
