#!/bin/bash
#SBATCH --job-name=n20volgsea
#SBATCH --output=niche20_volcano_gsea.%j.out
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

cd ${KB_PROJECT}
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate de_env

Rscript ${KB_REPO}/analyses/08_glomerular_niche20_crescent/niche20_volcano_gsea.R
