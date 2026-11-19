#!/bin/bash
# Slim SLURM driver for the PT -> aPT snRNA/snATAC trajectory and scMEGA eGRN analysis behind Extended Data Fig. 3k-n
# and Supplementary Fig. 8a-d (the AKI / pt_C5 exclusion rerun on the glue_v3 co-embedding, 26 Oct 2026).
#
# It submits only the steps those panels need, with the conda environments, resources and commands of the original
# driver (aki_c5_exclusion_20261022/atac_rerun/run_atac_stage.sh, which chained ~60 scripts incl. alternative
# trajectories, reports and the cis-regulatory reviewer analyses). Every job runs the scripts SHIPPED in this repository
# ($KB_REPO/analyses/10_atac_pt_apt_egrn/...) with the ATAC workspace as working directory, exactly as the original
# jobs ran the workspace copies. Nothing runs on the login node.
#
# Steps and their real inputs (-> step that produces them; "ext" = input from outside this analysis):
#   02  02_attach_glue.py          ArchR gene-score matrix (ext) + glue_v3/ and glue/ glue_for_pseudotime.csv (ext;
#                                  glue/ = the same GLUE fit with ATAC cluster pt_C5 kept, compared in 02/03 only)
#   03  03_embeddings_factor1.py   objects/gsm_with_glue{,_v3}.h5ad (02), MCFA loadings (ext)
#   09  09_rna_trajectory.py       snRNA object adata_20260717.h5ad (ext), glue_v3 csv, loadings
#   12  12_joint_pseudotime.py     tables/embeddings_glue_v3.csv (03), tables/rna_trajectory_percell.csv (09)
#   15  15_dpt_distal_check.py     snRNA object (ext)                    -> tables/dPT_percell_scores.csv
#   19  19_via_joint.py            HANDOFF_joint_pseudotime.csv (12), dPT_percell_scores.csv (15)
#   25  25_scmega_export.py        HANDOFF (12), via_joint_pseudotime.csv (19), snRNA object + peak_mat2.h5ad (ext)
#   26  26_panels_a_b.py           HANDOFF (12), sample_metadata_milo.csv (ext) -> tables/trajectory_pt_traj.csv
#   27  27_scmega_run.R            scmega/inputs (25), trajectory_pt_traj.csv (26) -> scmega/outputs/obj_pair.rds
#   31  31_archr_trajectory.R      obj_pair.rds (27), glue_latent.csv (25), HANDOFF (12) -> tables/archr_trajectories.csv
#   34  34_panels_a_b_archr.py     HANDOFF (12), archr_trajectories.csv (31), sample_metadata_milo.csv (ext)
#                                  [ED 3k-l; it does not need 28 or 33, which the original 34_figs job also ran]
#   43  43_scmega_egrn_tutorial.R  obj_pair.rds + rna_gene_detection.csv (27), archr_trajectories.csv (31); run twice:
#                                  archr_dm (P2G FDR 1e-50: network, S8a) and archr_dm_paper (P2G FDR 1e-4: ED 3m insets)
#   f1  make_factor1_pt_list.py    MCFA loadings (ext) -> the Factor-1 top-positive PT gene list read by 45, 47, 48
#                                  (replaces the cis-regulatory chain c01-c07 for the two columns used)
#   45  45 -> 46 -> 47 -> 48       egrn_tutorial/archr_dm (43), gene list (f1); 47 reads 45's edge table, 48 reads 47's
# Not needed for the panels and therefore not submitted: 13, 14, 16-18, 20-22, 07*_via_gsea, 28, 29, 30, 32, 33, 34b,
# 35-42, 44, the 34 variants without ONLY_USED / for archr_glue, 43 TRAJ=dc1rank and the cis-regulatory chain.
#
# Usage (from the repository root):
#   source config/env.sh
#   bash analyses/10_atac_pt_apt_egrn/run_pipeline.sh            # whole chain
#   SMOKE=1 bash analyses/10_atac_pt_apt_egrn/run_pipeline.sh    # scripts' own smoke mode (12, 19, 25, 27, 31, 43)
#   START=43 bash analyses/10_atac_pt_apt_egrn/run_pipeline.sh   # resubmit from a step on (earlier steps assumed done)
#   DEP=<jobid[:jobid]> ...                                      # extra afterok dependency for steps whose inputs come
#                                                                # from outside this submission
#   LOGDIR=<dir> ...                                             # job scripts + logs (default <workspace>/logs/run_pipeline)
# Steps depend on the steps that produce their inputs (afterok), so independent steps (09, 15) run in parallel.
# The workspace ($KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC) must contain glue_v3/ and glue/
# glue_for_pseudotime.csv (the GLUE co-embeddings);
# objects/, tables/ and scmega/ are created by the steps (on the Aachen HPC they are symlinks into $KB_DATA).
set -eo pipefail
: "${KB_REPO:?run 'source config/env.sh' in the repository first}"
S=$KB_REPO/analyses/10_atac_pt_apt_egrn/ATAC/scripts
A=$KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC
L=${LOGDIR:-$A/logs/run_pipeline}
for f in "$A/glue_v3/glue_for_pseudotime.csv" "$A/glue/glue_for_pseudotime.csv" "$KB_DATA/New_ATAC_data/gene_score_matrix.h5ad" \
         "$KB_DATA/New_ATAC_data/peak_mat2.h5ad" "$KB_DATA/charlotte_new_data/adata_20260717.h5ad" \
         "$KB_ANALYSIS/loadings~hvg~5.csv" "$KB_ANALYSIS/vivien_takeover/figures/milo_donor_composition/sample_metadata_milo.csv"; do
  [ -e "$f" ] || { echo "missing input: $f"; exit 1; }
done
mkdir -p "$L"
SMOKE=${SMOKE:-0}; START=${START:-}; DEP=${DEP:-}; ACTIVE=1; [ -n "$START" ] && ACTIVE=0
declare -A JOB
step(){ # name "producer steps" env cpus mem time cmd...   (' && ' separates commands; each becomes its own line)
  local name=$1 deps=$2 env=$3 cpus=$4 mem=$5 tm=$6; shift 6
  if [ $ACTIVE -eq 0 ]; then [ "$name" = "$START" ] && ACTIVE=1 || return 0; fi
  local f=$L/job_$name.sh after="" d cmds
  cmds=$(printf '%s' "$*" | sed 's/ && /\n/g')
  cat > "$f" <<EOF
#!/bin/bash
#SBATCH --job-name=atac_$name
#SBATCH --partition=defq
#SBATCH --cpus-per-task=$cpus
#SBATCH --mem=$mem
#SBATCH --time=$tm
#SBATCH --output=$L/%x.%j.out
source "$KB_REPO/config/env.sh"
cd "$A"
for v in \$(env | grep -oE '^CONDA[A-Z_]*'); do unset "\$v" || true; done
source "\$CONDA_ROOT/etc/profile.d/conda.sh"
conda activate $env
set -eo pipefail
export LD_LIBRARY_PATH=\$CONDA_PREFIX/lib:\$LD_LIBRARY_PATH OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK SMOKE=$SMOKE
echo "[\$(date)] start \$SLURM_JOB_NAME | env $env \$(conda list --explicit 2>/dev/null | md5sum | cut -c1-12)"
$cmds
echo "[\$(date)] done \$SLURM_JOB_NAME"
EOF
  for d in $deps; do if [ -n "${JOB[$d]:-}" ]; then after="$after:${JOB[$d]}"; fi; done
  if [ -z "$after" ] && [ -n "$DEP" ]; then after=":$DEP"; fi
  JOB[$name]=$(sbatch --parsable ${after:+--dependency=afterok$after} "$f")
  echo "$name ${JOB[$name]} after${after:-: none}" | tee -a "$L/jobids.txt"
}
if [ -z "$START" ]; then : > "$L/jobids.txt"; fi
# name  producers  env        cpus mem   time      command(s)
step 02 ""         agent_env  8  200G 03:00:00 "python $S/02_attach_glue.py"
step 03 "02"       ViaEnv     16 250G 06:00:00 "python $S/03_embeddings_factor1.py"
step 09 ""         ViaEnv     16 500G 12:00:00 "python $S/09_rna_trajectory.py"
step 12 "03 09"    ViaEnv     16 300G 06:00:00 "python $S/12_joint_pseudotime.py"
step 15 ""         ViaEnv     8  300G 04:00:00 "python $S/15_dpt_distal_check.py"
step 19 "12 15"    ViaEnv     16 300G 08:00:00 "python $S/19_via_joint.py"
step 25 "12 19"    ViaEnv     8  300G 04:00:00 "python $S/25_scmega_export.py"
step 26 "12"       agent_env  4  64G  02:00:00 "python $S/26_panels_a_b.py"
step 27 "25 26"    scmega2    16 400G 24:00:00 "Rscript $S/27_scmega_run.R"
step 31 "12 25 27" scmega2    4  64G  02:00:00 "Rscript $S/31_archr_trajectory.R"
step 34 "12 31"    agent_env  4  64G  02:00:00 "TRAJ=archr_dm ONLY_USED=1 FIGDIR=scmega/figures_paper python $S/34_panels_a_b_archr.py"
step 43 "27 31"    scmega2    8  200G 06:00:00 "TRAJ=archr_dm Rscript $S/43_scmega_egrn_tutorial.R && TRAJ=archr_dm OUT_SUFFIX=_paper P2G_FDR=1e-4 Rscript $S/43_scmega_egrn_tutorial.R"
step f1 ""         agent_env  1  8G   00:15:00 "python $KB_REPO/analyses/10_atac_pt_apt_egrn/make_factor1_pt_list.py"
step 45 "43 f1"    agent_env  4  32G  02:00:00 "python $S/45_egrn_paper_views.py --run archr_dm && python $S/46_egrn_module_timing.py && python $S/47_egrn_factor1_tf_side.py --highlight RELA,REL,NFKB1,NFKB2,RELB,FOS,JUN,JUNB,JDP2,FOSL2,HNF4A,RXRA,RXRB,RXRG,PPARA,PPARG && python $S/48_egrn_tf_centric_panels.py"
echo "submitted; job ids in $L/jobids.txt"
