#!/bin/bash
#SBATCH --job-name=atac_43_egrn_paper
#SBATCH --partition=defq
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=logs/atac_43_egrn_paper.%j.out
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"
cd ${KB_ANALYSIS}/aki_c5_exclusion_20261022/atac_rerun/ATAC
for v in $(env | grep -oE '^CONDA[A-Z_]*'); do unset "$v" || true; done
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate scmega2
set -eo pipefail
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK SMOKE=0
echo "[$(date)] start $SLURM_JOB_NAME; env hash $(conda list --explicit 2>/dev/null | md5sum | cut -c1-12)"
# Sept-24 'paper' variant of the eGRN tutorial: loose peak-to-gene FDR (1e-4) so that the 36-TF pseudotime-curve set
# (incl. JUN, HNF4A) is reproduced on the AKI-free, pt_C5-free embedding (glue_v3). Overwrites the rerun's
# archr_dm_paper folder, which the 26 Oct chain had produced with default knobs (duplicate of archr_dm).
TRAJ=archr_dm OUT_SUFFIX=_paper P2G_FDR=1e-4 Rscript ${KB_REPO}/analyses/10_atac_pt_apt_egrn/ATAC/scripts/43_scmega_egrn_tutorial.R
echo "[$(date)] done $SLURM_JOB_NAME"
