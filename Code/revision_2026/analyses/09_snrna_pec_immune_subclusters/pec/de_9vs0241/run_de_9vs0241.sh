#!/bin/bash
### MAST DE (leiden_0.8 cl9 vs cl0/2/4/11) + GSEA on the new PEC data
#SBATCH --time=1-00:00:00
#SBATCH --job-name=de_9vs0241
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=128GB
#SBATCH --output=logs/de_9vs0241.%J.out
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate de_env
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8

cd ${KB_ANALYSIS}/pec_resubcluster_20260717/de_9vs0241
echo "Started $(date) on $(hostname)"
Rscript ${KB_REPO}/analyses/09_snrna_pec_immune_subclusters/pec/de_9vs0241/de_analysis_9vs0241.R
echo "Completed $(date)"
