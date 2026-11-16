#!/bin/bash
### PART E (PEC recovery): two-panel cell-boundary per-core plots (niche | level2_refined)
### for overall(0.4) + in_glom with the updated level2_refined. in_glom uses --full-context
### (whole core, non-niche cells grey, glom cells coloured on top).
#SBATCH --time=08:00:00
#SBATCH --job-name=nc_pecrecB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=128GB
#SBATCH --output=logs/partE_pecrec_boundaries.%J.out
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate spatial_env
echo "Environment: $(which python)"; python --version; echo "Started $(date)"

cd ${KB_ANALYSIS}
python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/nichecompass_replot_level2/replot_boundaries_level2.py --version overall
python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/nichecompass_replot_level2/replot_boundaries_level2.py --version in_glom --full-context

echo "Completed $(date)"
