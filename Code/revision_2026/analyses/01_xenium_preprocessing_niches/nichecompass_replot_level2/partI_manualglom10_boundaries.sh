#!/bin/bash
### PART I: cell-boundary per-core plots (niche | level2_refined) for the manual_glom res-1.0 run,
### --full-context (whole core, non-glom cells grey, glom cells coloured) to showcase the crescent
### niche 20. Current level2_refined (post PEC-recovery).
#SBATCH --time=06:00:00
#SBATCH --job-name=nc_mg10b
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=128GB
#SBATCH --output=logs/partI_manualglom10_boundaries.%J.out
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate spatial_env
echo "Environment: $(which python)"; python --version; echo "Started $(date)"

cd ${KB_ANALYSIS}
python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/nichecompass_replot_level2/replot_boundaries_level2.py --version manual_glom_10 --full-context

echo "Completed $(date)"
