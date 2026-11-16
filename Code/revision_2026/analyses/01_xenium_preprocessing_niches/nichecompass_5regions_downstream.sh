#!/bin/bash

### Job Parameters
#SBATCH --time=08:00:00
#SBATCH --job-name=niche_down
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --mem=128GB

#SBATCH --output=logs/niche_downstream.%J.out

source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate nichecompass

echo "Environment activated: $(which python)"
python --version
echo "Started $(date)"

cd ${KB_PROJECT}

python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/nichecompass_5regions_downstream.py \
    --run-dir ${KB_PROJECT}/niche_compass_results/artifacts/single_sample/22042026_221229

echo "Completed $(date)"
