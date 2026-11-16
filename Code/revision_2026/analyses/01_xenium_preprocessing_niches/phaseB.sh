#!/bin/bash

### Job Parameters
#SBATCH --time=02:00:00
#SBATCH --job-name=xenium_phaseB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --mem=96GB

#SBATCH --output=logs/phaseB.%J.out

source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate omicverse

echo "Environment activated: $(which python)"
python --version
echo "Started $(date)"

cd ${KB_PROJECT}

python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/phaseB_label_transfer_merge.py \
    --n-pcs 50 \
    --knn-k 15

echo "Completed $(date)"
