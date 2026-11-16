#!/bin/bash
#SBATCH --job-name=vt08_cci
#SBATCH --output=logs/vt08_cci.%A_%a.out
#SBATCH --error=logs/vt08_cci.%A_%a.out
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --array=0-15
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

# no `set -u` -- breaks conda activate scripts (ADDR2LINE: unbound variable)
cd ${KB_ANALYSIS}/vivien_takeover
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate cci_env
set -e
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Explicit, non-overlapping chunk bounds derived here (not inside python), so each output
# file is named by the exact root range it contains and a failed task can be resubmitted alone.
N_ROOTS=258375
N_TASKS=16
CHUNK=$(( (N_ROOTS + N_TASKS - 1) / N_TASKS ))
START=$(( SLURM_ARRAY_TASK_ID * CHUNK ))
END=$(( START + CHUNK ))
if [ "$END" -gt "$N_ROOTS" ]; then END=$N_ROOTS; fi
echo "array task $SLURM_ARRAY_TASK_ID -> roots [$START, $END)"

python ${KB_REPO}/analyses/03_xenium_factor1_neighbourhoods_cci/08_run_local_cci.py --start "$START" --end "$END"
