#!/bin/bash
### 2026-07-22 CLINICAL METADATA JOIN FIX -- regenerate the affected figures for the two
### runs Leon asked for: overall (leiden 0.4) and manual_glom res 1.0 (the niche-20 run).
###
### The per-core clinical table was rebuilt by celltype_concordance/15_update_metadata_slope.py
### with an identifier-based join (K-number / Apollo well / H-T-number, 119/119 unambiguous,
### validated against each core's own `condition`), replacing a token matcher that had
### mis-assigned eGFR for 44% of cores. Regenerates every metadata-dependent figure
### (egfr_group_*, ifta_group_*, slope_class_*, EGFR/IFTA/eGFRslope continuous) plus
### condition_* so niche-row ordering matches across panels.
###
### Deliberately NOT re-running the full downstream: that re-clusters the latent and
### rewrites the h5ads, which could renumber niches -- and the crescent is pinned to
### `latent_leiden_1.0 == 20`. Everything else in those folders (composition, GP/pathway,
### DE, boundaries) is independent of the clinical metadata and is unchanged by this fix.
###
### Writes to staging under the ep924610 quota, then syncs with per-file verification
### (Kramann_lab is at its 60 TB limit).
#SBATCH --time=04:00:00
#SBATCH --job-name=nc_joinfix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=256GB
#SBATCH --output=logs/partK_joinfix_regen.%J.out
source "${KB_REPO:?run 'source config/env.sh' in the repository before submitting}/config/env.sh"

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate nichecompass
echo "Environment: $(which python)"; python --version; echo "Started $(date)"

cd ${KB_ANALYSIS}
STAGE=${KB_ANALYSIS}/figures_regen_20260722
mkdir -p "$STAGE"

for CFG in overall04 manual_glom10; do
  echo "=== regenerate metadata figures: $CFG ==="
  python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/nichecompass_replot_level2/regenerate_metadata_figures.py --config "$CFG" --dest-root "$STAGE"
  echo "=== sync: $CFG ==="
  python -u ${KB_REPO}/analyses/01_xenium_preprocessing_niches/nichecompass_replot_level2/sync_staged_figures.py --stage "$STAGE" --config "$CFG"
  echo "sync exit: $?"
done

echo "Completed $(date)"
