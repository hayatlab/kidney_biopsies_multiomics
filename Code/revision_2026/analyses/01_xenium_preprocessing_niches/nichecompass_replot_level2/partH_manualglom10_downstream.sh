#!/bin/bash
### PART G: additional MANUAL-glom NicheCompass run at leiden resolution 0.8, by re-clustering
### the existing trained manual_glom model latent (NO retrain), full downstream battery,
### cell type = level2_refined (post PEC-recovery). Goal: resolve a crescent-lesion niche
### (PEC + immune) analogous to niche 10 of the in_glom_niche run.
#SBATCH --time=10:00:00
#SBATCH --job-name=nc_mg10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --mem=200GB
#SBATCH --output=logs/partH_manualglom10_%J.out
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

set -euo pipefail
MODEL_DIR=${KB_PROJECT}/niche_compass_results/artifacts/manual_glom_niche/14052026_211003
OUT=${KB_RESULTS}/xenium/leon/new_data/manual_glom_niche_based_on_04_res10_level2

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate nichecompass
echo "Environment: $(which python)"; python --version
echo "Started $(date)  | RES=1.0  OUT=$OUT"

cd ${KB_ANALYSIS}
python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/nichecompass_5regions_downstream.py \
    --run-dir "$MODEL_DIR" \
    --leiden-res 1.0 \
    --cell-type-key level2_refined \
    --labeled-h5ad ${KB_ANALYSIS}/celltype_concordance/outputs/xenium_labeled.h5ad \
    --out-dir "$OUT" \
    --de-dir-name niche_marker_results \
    --out-h5ad-name xenium_with_niches.h5ad

echo "Completed $(date)"
