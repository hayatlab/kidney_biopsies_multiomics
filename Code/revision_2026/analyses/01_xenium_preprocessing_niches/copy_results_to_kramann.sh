#!/bin/bash

### Job Parameters
#SBATCH --time=01:00:00
#SBATCH --job-name=copy_glom_results
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --mem=32GB

#SBATCH --output=logs/copy_glom_results.%J.out

source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

# Copies the latest manual_glom_niche + overall_with_glom adata.h5ad files
# to $KB_RESULTS/xenium/leon/new_data/ and
# writes a README.md describing each.
#
# Submit with afterok deps on the two training jobs:
#   sbatch --dependency=afterok:<TRAIN_MANUAL_JOBID>:<TRAIN_OVERALL_JOBID> \
#       review_analyses/copy_results_to_kramann.sh

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate nichecompass

echo "Environment activated: $(which python)"
python --version
echo "Started $(date)"

cd ${KB_PROJECT}

python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/copy_results_to_kramann.py

echo "Completed $(date)"
