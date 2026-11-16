#!/bin/bash
### PART E (PEC recovery): re-plot cell-type-dependent level2 figures (composition /
### niche x celltype heatmap / niche annotation / per-niche DE) for overall(0.4) + in_glom
### using the updated level2_refined (182 aPT->PEC recovered in step 17). No model rerun.
#SBATCH --time=06:00:00
#SBATCH --job-name=nc_pecrecA
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --mem=256GB
#SBATCH --output=logs/partE_pecrec_replot.%J.out
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate nichecompass
echo "Environment: $(which python)"; python --version; echo "Started $(date)"

cd ${KB_ANALYSIS}
python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/nichecompass_replot_level2/replot_partA.py --version overall
python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/nichecompass_replot_level2/replot_partA.py --version in_glom

echo "Completed $(date)"
