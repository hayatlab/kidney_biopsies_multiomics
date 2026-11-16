#!/bin/bash
#SBATCH --job-name=vt34b_pairs
#SBATCH --output=logs/vt34b_pairs.%j.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --time=06:00:00
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"
# sbatch inherits the submitting shell's conda state; its deactivate hooks read CONDA_BACKUP_*
# which does not survive the export, and `set -u` then kills the job (job 330122).
set -eo pipefail
for v in $(env | grep -oE '^CONDA[A-Z_]*'); do unset "$v" || true; done
unset PYTHONPATH R_LIBS R_LIBS_USER || true
cd ${KB_ANALYSIS}/vivien_takeover
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate spatialdata
set -u
# The 11 pairs Leon asked for: the top 5 by mean aPT->Immune score, plus the pairs that lead
# the other three rankings (significant-root fraction, Factor-1 correlation, specificity),
# plus CD47-SIRPA from task 6.
# plus the scored pairs whose aPT-side molecule has >= 3 of 5 lines of evidence for being
# downstream of RELA (apt_immune_lr/outputs/rela_downstream_ligands.csv): JAG1 (4/5, incl.
# causal CRISPRa) and SEMA4D / CXCL12 (3/5).
PAIRS="APP_CD74,APP_SORL1,APLP2_PIGR,APP_PLXNA4,APLP2_PLXNA4,SEMA4D_PTPRC,CCL8_CCR2,EFNA5_EPHA4,TGFB1_TGFBR3,PECAM1_CD38,CD47_SIRPA,JAG1_NOTCH1,JAG1_NOTCH2,JAG1_NOTCH3,JAG1_CD46,SEMA4D_PLXNB1,SEMA4D_PLXNB2,CXCL12_DPP4"
python ${KB_REPO}/analyses/03_xenium_factor1_neighbourhoods_cci/34_cci_panel.py --pairs "${1:-$PAIRS}" --outsub cci_lr_panels
