#!/bin/bash
#SBATCH --job-name=snrna_figs_v2
#SBATCH --output=logs/snrna_figs_v2.%j.out
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"
# NB: no `set -u` -- it breaks conda's activate scripts.
cd ${KB_ANALYSIS}/snrna_overview_20260831
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate agent_env
set -e
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
python ${KB_REPO}/analyses/04_snrna_overview_composition/03_overview_figures_v2.py
