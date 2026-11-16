"""STEP 2: regenerate the metadata-dependent NicheCompass figures with the updated
per-core metadata (recovered cores + refreshed slopes), for every niche version.

For each version we merge the coalesced per-core metadata onto obs (canonical names,
no doubling), re-derive egfr_group / ifta_group / EGFR / IFTA_pct, and regenerate:
  - egfr_group_* , ifta_group_*        (enrichment: counts/props/binomial/fisher/
                                        signedlog10fdr heat+clustermap/bubble/stacked/
                                        effectsize/lollipop)  -- now with the ~64 cores
                                        that previously had no eGFR/IFTA
  - EGFR_ecdf/violin , IFTA_ecdf/violin (continuous)
  - NEW slope figures (this update):
      slope_class_*                    (eGFR-slope classification enrichment per niche)
      eGFRslope_ecdf/violin            (continuous total eGFR slope per niche)
Condition figures are unchanged (condition never had gaps) and are left as-is.

Writes into each version's figures/ dir (predic + *_level2 + 0.5/0.6). No model needed.
"""
from __future__ import annotations
from kbpaths import KB_RESULTS  # repository path config (config/kbpaths.py)
import argparse
import sys
import traceback
from pathlib import Path
import numpy as np
import scanpy as sc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))
from metadata_common import apply_core_metadata, load_percore   # noqa: E402
import nichecompass_5regions_downstream as ds                    # noqa: E402

NEW = Path(f"{KB_RESULTS}/xenium/leon/new_data")
ARCH = Path(f"{KB_RESULTS}/xenium/leon/archived")
# Paths follow the 2026-07-03 reorg (predic folders moved to archived/) and the
# 2026-07-22 move of the overall object next to its own figures.
CONFIGS = {
    "overall04": dict(
        adata=NEW / "overall_with_glom_04_level2/xenium_with_niches_overall_with_glom.h5ad",
        niche_key="latent_leiden_0.4",
        dirs=[NEW / "overall_with_glom_04_level2/figures",
              ARCH / "overall_with_glom_04/figures"]),
    "manual_glom": dict(
        adata=ARCH / "manual_glom_niche_based_on_04/xenium_with_niches_manual_glom_niche.h5ad",
        niche_key="latent_leiden_0.6",
        dirs=[NEW / "manual_glom_niche_based_on_04_level2/figures",
              ARCH / "manual_glom_niche_based_on_04/figures",
              ARCH / "manual_glom_niche_based_on_04_level2/figures"]),
    "manual_glom08": dict(
        adata=NEW / "manual_glom_niche_based_on_04_res08_level2/xenium_with_niches.h5ad",
        niche_key="latent_leiden_0.8",
        dirs=[NEW / "manual_glom_niche_based_on_04_res08_level2/figures"]),
    "manual_glom10": dict(
        adata=NEW / "manual_glom_niche_based_on_04_res10_level2/xenium_with_niches.h5ad",
        niche_key="latent_leiden_1.0",
        dirs=[NEW / "manual_glom_niche_based_on_04_res10_level2/figures"]),
    "in_glom": dict(
        adata=NEW / "in_glom_niche_based_on_04/xenium_with_niches_in_glom_niche.h5ad",
        niche_key="latent_leiden_0.6",
        dirs=[NEW / "in_glom_niche_based_on_04/figures",
              NEW / "in_glom_niche_based_on_04_level2/figures",
              ARCH / "in_glom_niche_based_on_04/figures"]),
}
# NOTE (2026-07-22): the `overall_with_glom_05` / `_06` resolution-comparison runs were
# DELETED from the results tree (folders + h5ads) while freeing Kramann_lab quota, so their
# configs are gone. Both are reproducible without retraining — re-cluster the latent
# embedding of `niche_compass_results/artifacts/overall_with_glom/14052026_211003` at leiden
# 0.5 / 0.6 (see partB_downstream.sh) — and their last corrected figure sets are in
# `review_analyses/figures_regen_20260722/{overall05,overall06}/figures` until that
# staging dir is cleaned up.
CATS = [("egfr_group", "egfr_group"), ("ifta_group", "ifta_group"),
        ("slope_class", "slope_class"),
        # `condition` is metadata-independent, but regenerating it here keeps the
        # niche-row ordering consistent with the other enrichment panels
        ("condition", "condition")]
CONT = [("EGFR", "EGFR"), ("IFTA_pct", "IFTA"), ("egfr_slope", "eGFRslope")]


def _num_slope(s):
    s = str(s).strip().replace(",", ".")
    try:
        return float(s)
    except ValueError:
        return np.nan


def regen_dir(adata, fig_dir: Path):
    fig_dir.mkdir(parents=True, exist_ok=True)
    for cat, pref in CATS:
        if cat not in adata.obs or adata.obs[cat].notna().sum() == 0:
            print(f"   [skip] {cat}: no values", flush=True); continue
        try:
            fish = ds.run_categorical_enrichment(adata, cat, fig_dir, pref, min_niche_cells=20)
            if fish is not None and not fish.empty:
                ds.plot_bubble_enrichment(fish, cat, fig_dir / f"{pref}_bubble",
                                          f"{pref} enrichment per niche (bubble)", min_n=20)
                ds.plot_stacked_proportions(adata, cat, fig_dir / f"{pref}_stacked_bars", min_total=50)
                ds.plot_effectsize_heatmap(adata, cat, fig_dir / f"{pref}_effectsize_heatmap")
                pick = fish.sort_values("pval_adj").iloc[0]
                ds.plot_lollipop(adata, cat, pick["niche"],
                                 fig_dir / f"{pref}_lollipop_{pick['niche']}", 20)
            print(f"   [ok] {pref}", flush=True)
        except Exception as exc:
            print(f"   [FAIL] {pref}: {exc!r}", flush=True); traceback.print_exc()
    for vk, tag in CONT:
        if vk not in adata.obs or adata.obs[vk].notna().sum() == 0:
            print(f"   [skip] {vk}: no values", flush=True); continue
        try:
            ds.plot_continuous_ecdf(adata, vk, fig_dir / f"{tag}_ecdf_samplelevel", True)
            ds.plot_violin_by_niche(adata, vk, fig_dir / f"{tag}_ridgeline_violin")
            print(f"   [ok] {tag} continuous", flush=True)
        except Exception as exc:
            print(f"   [FAIL] {vk}: {exc!r}", flush=True); traceback.print_exc()


def load_config(cfg):
    """The config's object with the per-core metadata applied and the derived groups added, exactly as
    run_config uses it (also called by the figure entry scripts)."""
    adata = sc.read_h5ad(cfg["adata"])
    ds.NICHE_KEY = cfg["niche_key"]
    percore = load_percore()
    apply_core_metadata(adata, percore)
    ds.derive_ifta_egfr_groups(adata)                       # egfr_group/ifta_group/EGFR/IFTA_pct
    adata.obs["egfr_slope"] = adata.obs["supp_total_slope"].map(_num_slope)
    adata.obs["slope_class"] = (adata.obs["supp_classification"].astype(str).str.strip()
                                .str.capitalize().replace({"Nan": np.nan, "": np.nan}))
    return adata


def run_config(name, cfg, dest_root=None, max_dirs=None):
    print(f"\n{'='*70}\n[{name}] {cfg['adata'].name}  niche={cfg['niche_key']}\n{'='*70}", flush=True)
    adata = load_config(cfg)
    for k in ("egfr_group", "ifta_group", "slope_class"):
        n = adata.obs[k].notna().sum()
        print(f"  {k}: {n} cells, cats={sorted(set(adata.obs[k].dropna().astype(str)))}", flush=True)
    if dest_root:
        # staging mode: compute once per config into <dest_root>/<config>/figures.
        # (Writing per real dir would recompute identical figures and, worse, two of
        # a config's dirs can share a basename and collide.) Use
        # sync_staged_figures.py to push a staged config into its real figures/ dirs.
        dirs = [Path(dest_root) / name / "figures"]
    else:
        dirs = cfg["dirs"][:max_dirs] if max_dirs else cfg["dirs"]
    for d in dirs:
        print(f" -> {d}", flush=True)
        regen_dir(adata, d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", choices=list(CONFIGS) + ["all"], default="all")
    ap.add_argument("--dest-root", default="", help="scratch base for smoke")
    ap.add_argument("--max-dirs", type=int, default=0)
    args = ap.parse_args()
    todo = list(CONFIGS) if args.config == "all" else [args.config]
    for n in todo:
        run_config(n, CONFIGS[n], dest_root=args.dest_root or None,
                   max_dirs=args.max_dirs or None)
    print("\n[ALL DONE]", ", ".join(todo), flush=True)


if __name__ == "__main__":
    main()
