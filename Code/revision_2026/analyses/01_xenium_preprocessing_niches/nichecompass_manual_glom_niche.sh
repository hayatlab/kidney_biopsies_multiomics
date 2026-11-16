#!/bin/bash

### Job Parameters
#SBATCH --time=04:00:00
#SBATCH --job-name=niche_manual_glom
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --mem=96GB

#SBATCH --output=logs/niche_manual_glom.%J.out

source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate nichecompass

echo "Environment activated: $(which python)"
python --version
echo "Started $(date)"

cd ${KB_PROJECT}

# Step 1: refresh in_glom annotations for the two new Region_1 samples
# (idempotent — overwrites the parquet)
python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/annotate_in_glom_new_samples.py

# Step 2: NicheCompass training on cells with in_glom == True
python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/nichecompass_manual_glom_niche.py \
    --n-epochs 400

echo "Completed $(date)"
