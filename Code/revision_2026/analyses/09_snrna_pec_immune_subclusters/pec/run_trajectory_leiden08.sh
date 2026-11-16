#!/bin/bash
#SBATCH --job-name=pec_traj_l08
#SBATCH --output=logs/traj_leiden08.%j.out
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"
set -euo pipefail
cd ${KB_PROJECT}
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate ViaEnv
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16
cd review_analyses/pec_resubcluster_20260717
python ${KB_REPO}/analyses/09_snrna_pec_immune_subclusters/pec/trajectory_leiden08.py "$@"
