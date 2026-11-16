#!/bin/bash
#SBATCH --job-name=pec_leiden04_xenium
#SBATCH --output=logs/map_leiden04.%j.out
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=04:00:00
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"
set -euo pipefail
cd ${KB_PROJECT}
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate kpmp_repurposing
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=32 MKL_NUM_THREADS=32
cd review_analyses/pec_resubcluster_20260717/xenium_pec_mapping
python ${KB_REPO}/analyses/09_snrna_pec_immune_subclusters/pec/xenium_pec_mapping/01_map_leiden04_to_xenium.py "$@"      # e.g. sbatch run_map.sh --label leiden_0.1
