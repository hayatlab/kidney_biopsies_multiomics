#!/bin/bash
### GSEA re-run on the non-PEC-marker-filtered cl9-vs-{0,2,4,11} MAST DE (no MAST re-run)
#SBATCH --time=02:00:00
#SBATCH --job-name=gsea_pecfilt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=32GB
#SBATCH --output=logs/gsea_pecfiltered.%J.out
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate de_env
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8

cd ${KB_ANALYSIS}/pec_resubcluster_20260717/de_9vs0241
echo "Started $(date) on $(hostname)"
Rscript ${KB_REPO}/analyses/09_snrna_pec_immune_subclusters/pec/de_9vs0241/gsea_pecfiltered.R
echo "Completed $(date)"
