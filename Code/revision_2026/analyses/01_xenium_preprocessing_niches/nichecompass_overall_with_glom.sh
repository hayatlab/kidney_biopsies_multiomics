#!/bin/bash

### Job Parameters
#SBATCH --time=08:00:00
#SBATCH --job-name=niche_overall_glom
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --mem=128GB

#SBATCH --output=logs/niche_overall_glom.%J.out

source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate nichecompass

echo "Environment activated: $(which python)"
python --version
echo "Started $(date)"

cd ${KB_PROJECT}

# Ensure the in_glom annotations parquet for the 2 new Region_1 samples is fresh
python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/annotate_in_glom_new_samples.py

# Re-run the overall (5-region) NicheCompass; output goes to a NEW dir
# under niche_compass_results/artifacts/overall_with_glom/<ts>/ — does not
# overwrite the existing single_sample/22042026_221229 run.
python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/nichecompass_overall_with_glom.py \
    --leiden-res 0.4 \
    --n-epochs 400

echo "Completed $(date)"
