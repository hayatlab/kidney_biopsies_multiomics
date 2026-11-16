#!/bin/bash
### 2026-07-22 fixes, applied to every niche version:
###  (1) eGFR groups: 4 KDIGO-style bins (<30, 30-60, 60-90, >90) in clinical order
###      -- was 3 bins (<30, 30-60, >=60) displayed alphabetically ('30-60','<30','>=60').
###      IFTA bins unchanged but now also ordered (<5, 5-25, 26-50, >50); niche rows
###      now sorted numerically instead of lexically ('2' before '10').
###  (2) hallmark dedup heatmap/dotplot: ULM now scored on log-normalised counts
###      (was raw counts), so the dedup view matches the "normal" matrixplot ranking
###      -- restores TNFA_SIGNALING_VIA_NFKB / INFLAMMATORY_RESPONSE / EMT / ANGIOGENESIS.
###
### Writes to a STAGING dir under the user's own quota, because the Kramann_lab group
### quota is at its 60 TB limit and direct writes die mid-savefig with
### OSError(122) Disk quota exceeded, leaving truncated PNG/PDF behind (job 328171).
### Push into the results tree afterwards with sync_staged_figures.py.
#SBATCH --time=08:00:00
#SBATCH --job-name=nc_egfr_hallmark
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=256GB
#SBATCH --output=logs/partJ_egfr_hallmark_fix.%J.out
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate nichecompass
echo "Environment: $(which python)"; python --version; echo "Started $(date)"

cd ${KB_ANALYSIS}
STAGE=${KB_ANALYSIS}/figures_regen_20260722
mkdir -p "$STAGE"

echo "=== (1) metadata figures: eGFR 4 bins + ordering ==="
python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/nichecompass_replot_level2/regenerate_metadata_figures.py --config all --dest-root "$STAGE"

echo "=== (2) hallmark dedup figures: log-normalised ULM ==="
python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/nichecompass_replot_level2/regenerate_hallmark_dedup.py --config all --dest-root "$STAGE"

echo "=== (3) sync staged figures into the results tree (verifies every file) ==="
python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/nichecompass_replot_level2/sync_staged_figures.py --stage "$STAGE"
echo "sync exit: $?"

echo "Completed $(date)"
