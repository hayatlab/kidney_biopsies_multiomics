"""Supplementary Fig. 2a-b — snRNA-seq overview and QC: UMAPs coloured by sex, condition and sample (a); per-sample
transcripts, genes and mitochondrial fraction per nucleus (violins), nuclei per sample and level-2 cell-type
composition per sample (b). 363,656 nuclei, 63 samples.
Panel        : All Supp Figures.pdf, Supplementary Figure 2, panels a and b
Reads        : $KB_ANALYSIS/snrna_overview_20260831/umap_and_obs.parquet
                 (analyses/04_snrna_overview_composition/01b_umap_from_csv.py)
               $KB_ANALYSIS/vivien_takeover/outputs/metadata/sample_metadata_corrected.csv (only the columns sample,
                 disease and sex_clinical are used; producer vivien_takeover/scripts/24_apply_corrected_metadata.py, not
                 part of this repository; the table also holds a patient_ID column, which is never written out)
               figures/_lib/snrna_palette_20260917.json (frozen disease / level1 / level2 palette, see below)
Writes       : $KB_FIG_OUT/S2/S2a-b_snrna_overview_qc/{pdf,png,svg}/ — panel a: umap_sex, umap_disease, umap_sample;
               panel b: violin_nCount_RNA, violin_nFeature_RNA, violin_pct_counts_mt, bar_nuclei_per_sample,
               stacked_level2; plus the script's other figures (umap_overview, umap_level2, stacked_level1), tables/
               (numbers behind each figure) and run_meta_v2.json
Plotting code: analyses/04_snrna_overview_composition/03_overview_figures_v2.py::main (imported unchanged; module-level
               FIG redirected to $KB_FIG_OUT and LABELS to the frozen palette). The original read the palette from
               kidney-biopsy-explorer/public/data/meta/labels.json, which was rewritten on 2026-10-02, after these
               figures (2026-09-17); the frozen copy holds the 2026-09-14 entries, identical to the current ones for
               every key the script uses, so the colours do not change
Illustrator  : a: umap_sex / umap_disease / umap_sample in one row, titles retyped (Sex -> sex, Disease -> condition,
               Sample (63 biopsies) -> sample id), sample legend removed, UMAP axis arrows kept on the first UMAP only;
               b: the five per-sample figures stacked in one column, x-axis titles ("Sample (grouped by disease)")
               removed, disease legends kept; panel letters. The level-2 label "dPT" is printed as in the data (other
               figures call this population aPT2)
Revision     : 2026-11-13, dPT is printed as aPT2 (co-authors' rename; same colour). KB_ORIGINAL_LABELS=1 draws "dPT"
Environment  : conda env agent_env
Verified     : 2026-11-12, PIXEL-IDENTICAL for all 11 figures (PNG byte-identical; all six tables identical) vs
               $KB_ANALYSIS/snrna_overview_20260831/figures/pdf/<figure>.pdf (2026-09-17)
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2] / "config"))
from kbpaths import fig_out  # noqa: E402

SCRIPT = HERE.parents[2] / "analyses" / "04_snrna_overview_composition" / "03_overview_figures_v2.py"
spec = importlib.util.spec_from_file_location("kb_snrna_overview_v2", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)                     # module level: paths, palette constants, matplotlib rcParams only

mod.FIG = Path(fig_out("S2", "S2a-b_snrna_overview_qc"))
mod.LABELS = HERE.parents[1] / "_lib" / "snrna_palette_20260917.json"
import os  # noqa: E402
mod.LEVEL2_DISPLAY = {} if os.environ.get("KB_ORIGINAL_LABELS") == "1" else {"dPT": "aPT2"}
mod.main()
print(f"[S2a-b] panel files in {mod.FIG}/pdf: a = umap_sex, umap_disease, umap_sample; b = violin_nCount_RNA, "
      "violin_nFeature_RNA, violin_pct_counts_mt, bar_nuclei_per_sample, stacked_level2", flush=True)
