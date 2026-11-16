#!/usr/bin/env python
"""Fig. 6g (left) — snRNA-seq PEC subclusters (Leiden 0.8) transferred onto Xenium PECs: glomerulus montages of
the niche-20 / crescent core K50 (IgAN) and of a TIN core with a "healthy" glomerulus (K18).
Panel        : Main Figures.pdf, Figure 6, panel g (left half, "predicted from snRNA-seq using PEC subclusters")
Reads        : $KB_ANALYSIS/pec_resubcluster_20260717/xenium_pec_mapping/pec_leiden08_predictions.csv
                 (analyses/09_snrna_pec_immune_subclusters/pec/xenium_pec_mapping/01_map_leiden04_to_xenium.py
                  --label leiden_0.8, wrapper run_map.sh)
               $KB_ANALYSIS/celltype_concordance/outputs/xenium_labeled.h5ad, obs only (analyses/02_xenium_celltyping/)
               $KB_XENIUM_RAW/output-<region>/{cell,nucleus}_boundaries.parquet (external: Xenium Onboard Analysis)
Writes       : $KB_FIG_OUT/Fig6/Fig6g_pec_subclusters_K50_XEN093.png, $KB_FIG_OUT/Fig6/Fig6g_pec_subclusters_K18_XEN116.png
Plotting code: analyses/09_snrna_pec_immune_subclusters/pec/xenium_pec_mapping/02_plot_pec_boundaries.py::plot_core
               (imported; patch rule build/patch_rules/fig_niches.py adds display/stem/overview keywords whose
               defaults keep the original behaviour); per-region loading extracted verbatim from ::main
               (original run: python 02_plot_pec_boundaries.py --label leiden_0.8 -> boundary_plots_leiden08/per_glomerulus/)
Illustrator  : one glomerulus per core cut out with a circular mask; legend relabelled ("PEC cl N" -> "N PEC",
               "POD (glom tuft)" -> "POD"); titles and per-glomerulus captions removed
Environment  : conda env kpmp_repurposing (the original plots were drawn in conda base; same pixels)
Verified     : 2026-11-12, DIFFERENT vs $KB_ANALYSIS/pec_resubcluster_20260717/xenium_pec_mapping/boundary_plots_leiden08/
               per_glomerulus/<core label>.png (K50 and K18) — only the montage title differs (it named the core by
               its internal label; here by its public id); with the title band masked the images are pixel-identical
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from niches_snrna import ANALYSES_09, CORES, core_display, load_script  # noqa: E402
from kbpaths import fig_out  # noqa: E402
from cores import core_label  # noqa: E402

LABEL = "leiden_0.8"
MIN_PEC = 5

pb = load_script(ANALYSES_09 / "pec" / "xenium_pec_mapping" / "02_plot_pec_boundaries.py")


def main():
    # ---- extracted from 02_plot_pec_boundaries.py::main (label leiden_0.8) ----
    tag = "leiden" + LABEL.split("_", 1)[1].replace(".", "")
    pred_col = f"pred_{tag}"
    pred_df = pd.read_csv(pb.HERE / f"pec_{tag}_predictions.csv")
    pred_df["cell_id"] = pred_df["cell_id"].astype(str)
    pred_df["sample"] = pred_df["sample"].astype(str)
    pred_map = pred_df.set_index(["sample", "cell_id"])[pred_col]

    obs = sc.read_h5ad(pb.XENIUM_LABELED, backed="r").obs[
        ["cell_id", pb.SAMPLE_KEY, pb.TMA_KEY, pb.CT_KEY]].copy()
    for c in obs.columns:
        obs[c] = obs[c].astype(str)

    out_dir = Path(fig_out("Fig6", "Fig6g_pec_subclusters_placeholder")).parent
    cache = {}
    for xen in CORES:
        tma = core_label(xen)
        regions = obs.loc[obs[pb.TMA_KEY] == tma, pb.SAMPLE_KEY].unique()
        if len(regions) != 1:
            raise SystemExit(f"{xen}: found in {len(regions)} Xenium regions (expected 1)")
        s = regions[0]
        if s not in cache:
            folder = pb.XENIUM_BASE / f"output-{s}"
            cb, nb = folder / "cell_boundaries.parquet", folder / "nucleus_boundaries.parquet"
            cell_ids, polys, cents = pb.load_polys(cb)
            nuc_polys = nuc_cents = None
            if nb.exists():
                _, nuc_polys, nuc_cents = pb.load_polys(nb)
            obs_s = obs[obs[pb.SAMPLE_KEY] == s].set_index("cell_id")
            lookup = pd.DataFrame(index=cell_ids)
            lookup[pb.TMA_KEY] = obs_s[pb.TMA_KEY].reindex(cell_ids).to_numpy()
            lookup[pb.CT_KEY] = obs_s[pb.CT_KEY].reindex(cell_ids).to_numpy()
            pmap_s = pred_map.loc[s] if s in pred_map.index.get_level_values(0) else pd.Series(dtype=object)
            lookup[pred_col] = pmap_s.reindex(cell_ids).to_numpy()
            cache[s] = dict(cell_ids=cell_ids, polys=polys, cents=cents, pred=lookup[pred_col].to_numpy(),
                            ct=lookup[pb.CT_KEY].astype(str).to_numpy(), nuc_polys=nuc_polys, nuc_cents=nuc_cents,
                            tma_arr=lookup[pb.TMA_KEY].astype(str).to_numpy())
        r = cache[s]
        n_core_pec = int(np.sum([pd.notna(p) for p in r["pred"][r["tma_arr"] == tma]]))
        if n_core_pec < MIN_PEC:
            raise SystemExit(f"{xen}: only {n_core_pec} mapped PEC cells (the original run skipped cores < {MIN_PEC})")
        reg = dict(cell_ids=r["cell_ids"], polys=r["polys"], cents=r["cents"], pred=r["pred"], ct=r["ct"],
                   nuc_polys=r["nuc_polys"], nuc_cents=r["nuc_cents"], in_core=(r["tma_arr"] == tma))
        display, stem = core_display(xen)
        n_pec, n_glom, n_scat = pb.plot_core(tma, reg, out_dir, out_dir, display=display,
                                             stem=f"Fig6g_pec_subclusters_{stem}", overview=False)
        print(f"[Fig6g PEC] {display}: {n_pec} mapped PEC cells, {n_glom} glomeruli, {n_scat} scattered "
              f"-> {out_dir}/Fig6g_pec_subclusters_{stem}.png", flush=True)


if __name__ == "__main__":
    main()
