# 10 — PT→aPT snRNA/snATAC trajectory and scMEGA enhancer GRN

## Purpose
Orders proximal-tubule (PT) and adaptive PT (aPT) nuclei of the snRNA-seq / snATAC-seq GLUE co-embedding along one
PT→aPT pseudotime, pairs RNA and ATAC nuclei with scMEGA, fits an ArchR trajectory on the joint diffusion map, infers
the enhancer gene-regulatory network (eGRN) along it and draws the network views of Extended Data Fig. 3j–m and
Supplementary Fig. 3. Everything here is the rerun of 26 Oct 2026 on the **glue_v3** co-embedding, from which the AKI
controls (K95, K96, K98; filtered by cell-id prefix) and the snATAC cluster pt_C5 were removed.

- **Embedding and pseudotime** — `02` attaches the 30-dim GLUE latent to the ArchR gene-score matrix; `03` draws
  UMAP / PHATE / diffusion maps with Factor 1 scores (QC views); `09` RNA-only trajectory and RNA-on-RNA Factor 1 scores;
  `12` one joint pseudotime for RNA and ATAC nuclei (rank of DC1 of the joint kNN graph, oriented by Factor 1) and the
  hand-off table; `15` checks the distal-nephron signal of dPT (its per-cell scores feed `19`); `19` VIA pseudotime as an
  independent check of the DC1 ordering.
- **scMEGA** — `25` exports RNA PT + aPT and ATAC pt_C3 + pt_C6 nuclei; `26` re-ranks the pseudotime within that set
  (`tables/trajectory_pt_traj.csv`, read by `27`); `27` (+ `27_downstream.R`) pairs RNA and ATAC nuclei on the GLUE latent
  and builds the paired object with chromVAR motif activities; `31` fits the ArchR trajectory (`scMEGA::AddTrajectory`)
  in the joint diffusion map (`archr_dm`, the published one) and in the GLUE latent.
- **eGRN** — `43` runs Mina Shaigan's scMEGA eGRN workflow (SelectTFs r > 0.7 / adj. p < 0.01; peak-to-gene links
  r > 0.7, FDR < 1e-50; TF–gene edges r > 0.7), and a second time with peak-to-gene FDR < 1e-4 (`OUT_SUFFIX=_paper`) for
  the per-TF pseudotime curves; `make_factor1_pt_list.py` writes the Factor 1 top-positive PT gene list; `45` network
  views (two disjoint modules), `46` module timing along the pseudotime, `47` TFs upstream of the Factor 1 genes, `48`
  TF-centred regulon rings.

## Analysis workspace
`$KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC` — the working directory of every job. It was created by
`atac_rerun/build_atac_rerun_tree.py` (not shipped), which copied the September scripts, rewrote their root paths,
renamed `glue_v2` to `glue_v3` and injected the AKI filters; the scripts shipped here are those patched copies (list of
injections: `atac_rerun/PATCHES_atac.md` in the workspace). `objects/`, `tables/`, `scmega/` and `figures/` are symlinks
into `$KB_DATA/aki_c5_exclusion_20261022/atac_rerun/ATAC/`. The Factor 1 gene list is written to
`$KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/cis_regulatory_pt_apt/outputs/task3/`.

## Run order
`source config/env.sh`, then `bash analyses/10_atac_pt_apt_egrn/run_pipeline.sh` submits the steps below as SLURM jobs,
each depending (afterok) on the steps that produce its inputs (`SMOKE=1`, `START=<step>`, `DEP=`, `LOGDIR=`: see the
script header). It replaces the original driver `atac_rerun/run_atac_stage.sh` (not shipped), which ran the same
commands in one linear chain together with ~45 scripts that no published panel needs (alternative trajectories,
reports, the cis-regulatory reviewer analyses). Resources and environments are those of the original jobs.

| step | script | conda env | CPUs / memory / time | original job | main outputs (workspace-relative) |
|---|---|---|---|---|---|
| 02 | `ATAC/scripts/02_attach_glue.py` | agent_env | 8 / 200G / 3 h | 339453 | `objects/gsm_with_glue{,_v3}.h5ad` |
| 03 | `ATAC/scripts/03_embeddings_factor1.py` | ViaEnv | 16 / 250G / 6 h | 339454 | `tables/embeddings_glue{,_v3}.csv`, `figures/glue_embeddings_factor1.*` |
| 09 | `ATAC/scripts/09_rna_trajectory.py` | ViaEnv | 16 / 500G / 12 h | 339455 | `tables/rna_trajectory_percell.csv`, `rna_gene_correlations.csv` |
| 12 | `ATAC/scripts/12_joint_pseudotime.py` | ViaEnv | 16 / 300G / 6 h | 339456 | `tables/HANDOFF_joint_pseudotime.csv`, `HANDOFF_atac_cells_pseudotime.csv`, `HANDOFF_params.json` |
| 15 | `ATAC/scripts/15_dpt_distal_check.py` | ViaEnv | 8 / 300G / 4 h | 339457 | `tables/dPT_percell_scores.csv`, `dPT_distal_*.csv` |
| 19 | `ATAC/scripts/19_via_joint.py` | ViaEnv | 16 / 300G / 8 h | 339458 | `tables/via_joint_pseudotime.csv`, `via_joint_summary.json` |
| 25 | `ATAC/scripts/25_scmega_export.py` | ViaEnv | 8 / 300G / 4 h | 339459 | `scmega/inputs/{rna_counts.h5, atac_peaks.h5, glue_latent.csv, metadata.csv}` |
| 26 | `ATAC/scripts/26_panels_a_b.py` | agent_env | 4 / 64G / 2 h | 339460 | `tables/trajectory_pt_traj.csv` (+ DC1-rank panels, `scmega/figures_dc1rank/`) |
| 27 | `ATAC/scripts/27_scmega_run.R` (sources `27_downstream.R`) | scmega2 | 16 / 400G / 24 h | 339461 | `scmega/outputs/obj_pair.rds`, `pairs.csv`, `rna_gene_detection.csv`, DC1-rank GRN tables |
| 31 | `ATAC/scripts/31_archr_trajectory.R` | scmega2 | 4 / 64G / 2 h | 339462 | `tables/archr_trajectories.csv` |
| 34 | `ATAC/scripts/34_panels_a_b_archr.py` with `TRAJ=archr_dm ONLY_USED=1 FIGDIR=scmega/figures_paper` | agent_env | 4 / 64G / 2 h | 339466 | `scmega/figures_paper/panel_a_b_archr_dm.*`, `scmega/outputs/panel_b_donor_densities_archr_dm.csv` |
| 43 | `ATAC/scripts/43_scmega_egrn_tutorial.R` with `TRAJ=archr_dm` | scmega2 | 8 / 200G / 6 h | 339468 | `scmega/egrn_tutorial/archr_dm/` (selected TFs, links, `grn_edges*.csv`, `tf_centrality.csv`, `pseudotime_curves.csv.gz`, figures) |
| 43 | same with `TRAJ=archr_dm OUT_SUFFIX=_paper P2G_FDR=1e-4` | scmega2 | 8 / 64G (wrapper `job_43_egrn_paper_p2g1e-4.sh`) | 340411 (9 Nov 2026) | `scmega/egrn_tutorial/archr_dm_paper/` incl. `figures/per_tf/scMEGA_TF_<TF>.pdf` |
| f1 | `make_factor1_pt_list.py` | agent_env | 1 / 8G | was part C of `cis_regulatory_pt_apt/scripts/07_task3_nfkb.py` (job 339553) | `cis_regulatory_pt_apt/outputs/task3/C4_factor1_top200_positive_genes.csv` |
| 45 | `45_egrn_paper_views.py --run archr_dm`, `46_egrn_module_timing.py`, `47_egrn_factor1_tf_side.py --highlight RELA,REL,NFKB1,NFKB2,RELB,FOS,JUN,JUNB,JDP2,FOSL2,HNF4A,RXRA,RXRB,RXRG,PPARA,PPARG`, `48_egrn_tf_centric_panels.py` | agent_env | 4 / 32G / 2 h | 339554 | `scmega/figures_paper/` (views, tables, meta json) |

`34` needs neither `28` nor `33` (the original 34_figs job ran them first only for other figures). `27` needs `26`
(`tables/trajectory_pt_traj.csv`). `09` and `15` need no output of `02`/`03` and run in parallel with them.
`make_factor1_pt_list.py` replaces the cis-regulatory chain (pseudobulks, limma, peak annotation, motif matching,
peak–gene links, 07) that wrote the gene list in the original run: 45, 47 and 48 read only its columns `variable` and
`Loading`, which the script rebuilds with the same selection code (output equal to the original on these columns, checked
12 Nov 2026). If the list already exists it is compared, not overwritten (`--force`, `--out`).

Environment hashes (md5 of `conda list --explicit`, first 12 characters, 26 Oct 2026): agent_env `f69893914162`,
scmega2 `2d4d4eec9756`, ViaEnv `c8024ca105ab`. scmega2: R 4.3.3, scMEGA 1.1.0, ArchR 1.0.3, Signac 1.14.0,
Seurat 5.3.0, JASPAR2020, BSgenome.Hsapiens.UCSC.hg38, EnsDb.Hsapiens.v86. Human hg38, HGNC symbols, seed 42.

## Inputs
| input | producer |
|---|---|
| `glue_v3/glue_for_pseudotime.csv` (148,571 nuclei: 128,655 RNA PT / aPT / dPT + 19,916 ATAC; columns `coembedding_id, cellid, modality, level1, level2, cluster, 0..29`) and `glue/glue_for_pseudotime.csv` (same GLUE fit with pt_C5 kept; compared in 02/03 only) | GLUE co-embedding delivered by the collaborator on 26 Oct 2026 (external) |
| `$KB_DATA/New_ATAC_data/gene_score_matrix.h5ad`, `$KB_DATA/New_ATAC_data/peak_mat2.h5ad` | ArchR project of the snATAC-seq data (external) |
| `$KB_DATA/charlotte_new_data/adata_20260717.h5ad` | snRNA-seq object with the level2 annotation (external to this analysis) |
| `$KB_ANALYSIS/loadings~hvg~5.csv` | MCFA model, 5 factors on highly variable genes (upstream MCFA workflow, external) |
| `$KB_ANALYSIS/vivien_takeover/figures/milo_donor_composition/sample_metadata_milo.csv` | `vivien_takeover/scripts/22_sample_metadata_table.py` (not shipped). It has a `patient_ID` column; 26 and 34 read only `sample`, `disease`, `egfr_group` |
| motif / genome annotation | R packages JASPAR2020, BSgenome.Hsapiens.UCSC.hg38, EnsDb.Hsapiens.v86 (scmega2) |

## Outputs and the panels that use them
| panel | entry script | output used |
|---|---|---|
| ED Fig. 3j–k | `figures/ED3/ED3j-k_archr_trajectory_egfr_density.py` | 34: `panel_a_b_archr_dm.pdf` (j = panel a, k = panel b′) |
| ED Fig. 3l, network | `figures/ED3/ED3l_egrn_two_modules.py` | 45: `egrn_archr_dm_two_modules.pdf` |
| ED Fig. 3l, JUN / HNF4A insets | `figures/ED3/ED3l_insets_scMEGA_TF_JUN_HNF4A.R` | 43 `_paper`: re-drawn from `pseudotime_curves.csv.gz` + `tf_centrality.csv` (the plotting code of 43, verbatim) |
| ED Fig. 3m | `figures/ED3/ED3m_egrn_module_timing.py` | 46: `egrn_module_timing.pdf`, panels A–C (bottom panel printed from September, see caveats) |
| Supp. Fig. 3a | `figures/S3/S3a_tf_activity_expression_heatmap.R` | 43 `archr_dm`: `tf_heatmap_cor70.pdf`, re-drawn from the saved curves |
| Supp. Fig. 3b | `figures/S3/S3b_egrn_factor1_tf_ranking.py` | 47: `egrn_factor1_tf_ranking.pdf` |
| Supp. Fig. 3c | `figures/S3/S3c_egrn_factor1_targets_overview.py` | 45: `egrn_factor1_targets_overview.pdf` |
| Supp. Fig. 3d | `figures/S3/S3d_egrn_tf_centric_FOSL2.py` | 48: `egrn_tf_centric/FOSL2.pdf` |

Every entry script reproduces its original file pixel-identically (12 Nov 2026; `build/verification/atac_snrna.csv`).
Release patches (output folders `OUT` / `TABLEDIR`, defaults = original folders): `build/patch_rules/fig_atac_snrna.py`.

## Caveats (read before reusing)
- **glue_v3 is a new GLUE fit.** AKI controls (K95, K96, K98) and pt_C5 are excluded; coordinates and pseudotimes are not
  comparable number by number with the September glue_v2 results. The scripts' docstrings still quote September numbers
  (17,637 pairs, 86 TFs, 3,686 edges, 388 targets, modules of 72 / 14 TFs, 46 Factor 1 targets). The rerun has 16,843 of
  18,715 pairs with an ArchR value (ArchR leaves the 10 % of pairs farthest from the curve without one), 80 TFs, 2,963
  edges, 337 targets, modules of 64 / 16 TFs (272 / 136 nodes) and 43 of the 80 Factor 1 genes as targets.
- **ED 3m mixes two runs.** Its top (TF activity) and middle (TF expression) panels are the rerun; the printed bottom panel
  (Target expression) is panel C of the September glue_v2 figure (`$KB_ANALYSIS/ATAC/scmega/figures_paper/egrn_module_timing.pdf`:
  module crossing at bin 60, transition window 58–72; rerun: bin 59, window 56–73). In the rerun the module II lead in
  chromatin is no longer significant (steepest change, permutation p = 0.079; September 0.005) while the lead in target
  expression holds (p = 0.0002). 46's suptitle and its caveat line ("p = 0.52") are September wording; they are not in
  the published panel.
- **Insets and network use different link thresholds.** The JUN / HNF4A curves come from the peak-to-gene FDR < 1e-4 run
  (JUN 517 targets, HNF4A 843); the network beside them from FDR < 1e-50 (JUN 51, HNF4A 57). TF selection (81 TFs) and the
  TF activity / expression curves do not depend on that threshold (identical in both runs).
- **ED 3l inset legend.** The legend added in Illustrator uses ggplot's default colours (red "Targets expression", green
  "TF activity", blue "TF expression"), but the curves keep 43's colours (orange = TF activity, green = targets, blue = TF
  expression): TF activity and Targets expression are swapped in the printed legend.
- **ED 3k misses K68 and counts samples as donors.** K68 (FSGS, eGFR 87 in the corrected metadata) has no disease or eGFR
  group in `sample_metadata_milo.csv`, so its 447 RNA nuclei of pairs are excluded. The "donors" of panels b / b′ are
  samples: K94 and K97 are two biopsies of one person and enter as two curves.
- **Stale labels in 47.** Its figure legend and meta json say "BH-FDR < 0.05 over 86 TFs" and "population = 388 targets,
  46 Factor-1 targets"; the test runs on the rerun network (80 TFs, 337 targets, 43 Factor 1 genes; 9 TFs at FDR < 0.05).
  The published Supp. Fig. 3b does not show that legend.
- **Supp. Fig. 3a is drawn from the `_paper` table.** The `archr_dm` run saves curves only for its 80 network TFs (NR1H4,
  selected but without an edge at FDR < 1e-50, is missing); the `_paper` table holds all 81 selected TFs and its curves
  for the 80 shared TFs are bit-identical to `archr_dm`'s.
- **Factor 1 gene list.** Top 200 genes by |loading| of the PT view of MCFA Factor 1, of which the 80 with a positive
  loading (120 of the 200 are negative); no tie at the 200th gene.
- The pseudotime is an ordering of nuclei, not elapsed time; module "leads" are leads in cell ordering.
