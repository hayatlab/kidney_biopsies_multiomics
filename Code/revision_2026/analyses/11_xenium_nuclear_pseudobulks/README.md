# 11 — Xenium cell type × sample pseudobulks (all transcripts and nuclear transcripts only)

## Purpose
Builds the cell type × sample count tables of the Xenium data that the multicellular factor analysis (MCFA) uses to
compute factor 1 on an imaging platform, which has no droplet-based ambient RNA (Supplementary Fig. 27; reviewer
question on ambient RNA). Two versions with identical cells, labels, groupings and sample ids:

1. all transcripts of a cell (`01_xenium_pseudobulk.py`, 26 Oct 2026);
2. only transcripts that overlap the nucleus (`04_nuclear_pseudobulk.py`, 27 Oct 2026) — the version named in the legend
   of Supplementary Fig. 27. Transcript rule: from each run's `transcripts.parquet`, gene transcripts (`is_gene`) with
   Q-score >= 20 assigned to a cell of the analysed object and with `overlaps_nucleus == 1`. With the Q >= 20 filter the
   per-cell totals reproduce the counts layer of the analysed object exactly (asserted in the script for every run).

Both sum raw counts with `decoupler.pp.pseudobulk(mode="sum")` (decoupler 2.1.4) per sample and per cell type, once for
the ten MCFA views (`mcfa_view`: PT = PT + aPT + dPT; Stromal = FIB + MES + VSM-P; CNT_PC = CNT + PC; IC = IC-A + IC-B;
EC, Immune, TAL, DCT, POD, PEC) and once for the refined labels (`level2_refined`). No normalisation and no filtering in
the main files; `psbulk_cells` and `psbulk_counts` are kept so that the user of the tables can filter. A sample is one
TMA core (one biopsy). Human, HGNC symbols, 468-gene Xenium panel. In the objects the label `dPT` is the population
called aPT2 in the paper.

The factor scores, the comparison with snRNA-seq and the panels of Supplementary Fig. 27 were made from these tables
by the MCFA authors; that code is not in this repository.

## Analysis workspaces
`$KB_ANALYSIS/mcfa_xenium_pseudobulk_20261026/` (step 1) and `$KB_ANALYSIS/xenium_zenodo_20261027/` (step 2; its
`outputs/` is `$KB_DATA/zenodo_xenium_20261027/`).

## Run order (source `config/env.sh` first)
| step | script | conda env | wrapper | output |
|---|---|---|---|---|
| 1 | `01_xenium_pseudobulk.py` | nichecompass | none (login node, about 2 min) | `mcfa_xenium_pseudobulk_20261026/outputs/pseudobulk_{mcfa_view,level2_refined}_{unfiltered,min10cells_1000counts}.h5ad` / `_counts.csv` / `_obs.csv`, `cells_per_sample_level2_refined.csv`, `xenium_panel_genes.csv`, `README.md`, and the private sample key `outputs/private/sample_key_PRIVATE.csv` |
| 2 | `04_nuclear_pseudobulk.py` (imports `common.py`) | nichecompass | `sbatch_04_nuclear_pseudobulk.sh` (8 CPUs, 96 GB; 8 min, 5 GB peak) | `$KB_DATA/zenodo_xenium_20261027/for_mcfa_nuclear_20261027/pseudobulk_nuclear_{mcfa_view,level2_refined}_unfiltered.h5ad` / `_counts.csv` / `_obs.csv`, `per_run_transcript_stats.json`, `README.md`, zip of these files |

Published runs: step 1 on 26 Oct 2026; step 2 SLURM job 339634, 27 Oct 2026. Result of step 2: 258,375 cells, 119
cores; 1,190 pseudobulks for `mcfa_view` (median 121.5 cells, 15,164.5 nuclear transcripts) and 1,904 for
`level2_refined`; nuclear fraction of the Q >= 20 gene transcripts 0.464 overall (0.436–0.487 per run); the total-count
pseudobulks computed from the transcript files equal the tables of step 1 (maximum absolute difference 0).

`common.py` is the helper module of the Xenium data preparation for the public repository (paths, run list, core key,
identifier scan). Only the parts used by step 2 are needed here; the other preparation steps (restricting the raw
bundles to the analysed cores, image crops, de-identified object) are not part of this folder.

Set `KB_NUCLEAR_PSEUDOBULK_OUT` to write the output of step 2 to another folder (default: the folder above).

## Inputs
| input | producer |
|---|---|
| `$KB_ANALYSIS/celltype_concordance/outputs/xenium_labeled.h5ad` (258,375 cells × 468 genes, `layers["counts"]`, `level2_refined`) | `analyses/02_xenium_celltyping` |
| `$KB_XENIUM_RAW/output-<run>/transcripts.parquet` of the five runs listed in `common.py` (step 2) | Xenium Onboard Analysis; external |
| `$KB_ANALYSIS/vivien_takeover/outputs/01_factor_scoring_meta.json` (refined label → MCFA view) | `analyses/03_xenium_factor1_neighbourhoods_cci/01_score_factors_xenium.py` |
| `$KB_ANALYSIS/celltype_concordance/outputs/per_sample_metadata_updated.csv` (K number, eGFR, IF/TA per core) | `analyses/02_xenium_celltyping/15_update_metadata_slope.py` |
| `$KB_ANALYSIS/vivien_takeover/outputs/authoritative_core_egfr.csv` (step 1) | `analyses/03_xenium_factor1_neighbourhoods_cci/10_build_authoritative_egfr.py` |
| `mcfa_xenium_pseudobulk_20261026/outputs/pseudobulk_mcfa_view_unfiltered{.h5ad,_obs.csv}` and `outputs/private/sample_key_PRIVATE.csv` (step 2) | step 1 |
| internal core table `$KB_ANALYSIS/INTERNAL_Xenium_old_vs_new_labels_DO_NOT_SHARE.csv` (internal core label → public `xenium_core_id` and sample name; step 2 via `common.py`) | private, never part of this repository |

## Outputs and the panels that use them
| output | panel |
|---|---|
| `pseudobulk_nuclear_mcfa_view_unfiltered.*` (step 2) | starting tables of Supplementary Fig. 27a–c (scores and panels by the MCFA authors) |
| tables of step 1 | no panel; reference for the equality check in step 2 and for the nuclear fraction (`total_counts_same_cells`) |

## Sample ids
`sample_id` is the snRNA-seq K number of the same biopsy where one exists (54 cores), otherwise `xen_core_###`. K
numbers come from the per-sample metadata table, which overrides the object for two cores of the newer run whose
K number in the object does not fit their diagnosis. A K number that labels two cores is kept as two samples
(`<K>_core1`, `<K>_core2`). The nuclear tables add `zenodo_sample_id`, the sample name of the public Xenium data.

## Re-run check (19 Nov 2026)
The copy of `04_nuclear_pseudobulk.py` in this folder was run with `KB_NUCLEAR_PSEUDOBULK_OUT` pointing to a separate
folder (SLURM job 342291, 7 min 45 s, 5 GB peak) and compared with the tables delivered on 27 Oct 2026: the four CSV
tables and `per_run_transcript_stats.json` are byte-identical, the two h5ad files have identical counts, layers and
obs, and the README differs only in its build date. `01_xenium_pseudobulk.py` was not re-run, because it rewrites the
private sample key in place; its path rewrite was checked for equivalence like every other file (identical syntax tree
with the default paths), and step 2 asserts that its totals equal the tables of step 1.

## Caveats (read before reusing)
- **Step 2 cannot be run from the public data as it stands**: it reads the internal core table and the private sample
  key, which contain pathology numbers. The public Xenium data hold the same transcripts (`overlaps_nucleus` included)
  and name samples by `zenodo_sample_id`, so the counting rule can be applied there; that was not run.
- Nuclear overlap is the per-transcript flag of the instrument software against the DAPI-based nucleus segmentation. It
  reduces misassignment between neighbouring cells but does not remove it (nuclei can touch).
- Cell labels are the scANVI transfer from the snRNA-seq reference that was trained with the two AKI samples present;
  the reference was not retrained after these samples were excluded from the snRNA-seq analyses.
- Step 1 also writes eGFR, IF/TA and progression columns per core; they are the values of 26 Oct 2026 and were not
  updated to the final clinical table (Supplementary Data 1).
- `outputs/private/sample_key_PRIVATE.csv` (written by step 1) maps the sample ids to the internal core labels and must
  stay private. Both scripts end with a scan of the shared files for identifier patterns and stop if one is found.
