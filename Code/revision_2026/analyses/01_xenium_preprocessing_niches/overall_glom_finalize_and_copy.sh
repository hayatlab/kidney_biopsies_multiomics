#!/bin/bash

### Job Parameters
#SBATCH --time=04:00:00
#SBATCH --job-name=glom_finalize_copy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --mem=128GB

#SBATCH --output=logs/glom_finalize_copy.%J.out

source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

# Resume the post-training steps from job 321201 that hit the 8h SLURM
# limit, then run glom_postprocess and copy h5ads to the Kramann_lab share.

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate nichecompass

echo "Environment activated: $(which python)"
python --version
echo "Started $(date)"

cd ${KB_PROJECT}

set -euo pipefail

# 1. Recover overall_with_glom: neighbors + UMAP + leiden + standalone adata write
python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/overall_with_glom_finalize.py --leiden-res 0.4

# 2. Post-processing on the now-complete overall run (heatmap + dotplot + niche-vs-in_glom)
OVR_RUN=$(ls -1dt ${KB_PROJECT}/niche_compass_results/artifacts/overall_with_glom/*/ | head -1 | sed 's:/$::')
echo "[post] overall run_dir=$OVR_RUN"
python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/glom_postprocess.py --run-dir "$OVR_RUN" --niche-key latent_leiden_0.4

# 3. Copy both h5ads + write README to Kramann_lab share
python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/copy_results_to_kramann.py

echo "Completed $(date)"
