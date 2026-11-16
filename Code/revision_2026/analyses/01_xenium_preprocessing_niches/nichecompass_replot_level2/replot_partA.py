"""PART A: re-plot the three existing (0.4-based) NicheCompass versions using the
refined `level2_refined` cell typing, WITHOUT re-running the niche model (niches
already stored in each object).

For each version we build a self-contained `*_level2/` sibling folder that holds:
  - regenerated cell-type-dependent figures (level2_refined):
      * res_<res>_niche_composition_level2_refined.*
      * niche_vs_celltype_heatmap[_noannot].*   (level2)
      * niche_annotation_summary.csv            (level2)
      * niche_marker_results/                   (per-niche DE on level2)
  - copies of every annotation-INDEPENDENT figure from the original folder
      (condition / EGFR / IFTA / pathway / GP / roi_flags / niches_per_TMA / …)
Cell-boundary plots are produced separately by replot_boundaries_level2.py
(needs geopandas / spatial_env) and land in <dest>/cell_boundary_plots_per_core/.

The predic figures in the original folders are left untouched.
"""
from __future__ import annotations
from kbpaths import KB_RESULTS  # repository path config (config/kbpaths.py)
import argparse
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))          # review_analyses (for the downstream module)
from level2_common import merge_level2, level2_colour_dict, ORDER   # noqa: E402
import nichecompass_5regions_downstream as ds                       # noqa: E402

NEW = Path(f"{KB_RESULTS}/xenium/leon/new_data")
# the bulky predic overall source object was relocated to archived/ (level2 dest stays in new_data)
ARCHIVED = Path(f"{KB_RESULTS}/xenium/leon/archived")
CT_KEY = "level2_refined"

VERSIONS = {
    "overall": dict(
        adata=ARCHIVED / "overall_with_glom_04/xenium_with_niches_overall_with_glom.h5ad",
        niche_key="latent_leiden_0.4", res=0.4,
        src=ARCHIVED / "overall_with_glom_04", dest=NEW / "overall_with_glom_04_level2"),
    "manual_glom": dict(
        adata=ARCHIVED / "manual_glom_niche_based_on_04/xenium_with_niches_manual_glom_niche.h5ad",
        niche_key="latent_leiden_0.6", res=0.6,
        src=ARCHIVED / "manual_glom_niche_based_on_04",
        dest=NEW / "manual_glom_niche_based_on_04_level2"),
    "in_glom": dict(
        adata=NEW / "in_glom_niche_based_on_04/xenium_with_niches_in_glom_niche.h5ad",
        niche_key="latent_leiden_0.6", res=0.6,
        src=NEW / "in_glom_niche_based_on_04",
        dest=NEW / "in_glom_niche_based_on_04_level2"),
}


def is_celltype_fig(name: str) -> bool:
    """Figures that DEPEND on cell type -> regenerated, so not copied over."""
    return ("_niche_composition_predic" in name
            or name.startswith("niche_vs_celltype_heatmap")
            or name == "niche_annotation_summary.csv")


def heatmap_noannot(adata, niche_key: str, ct_key: str, out_base: Path) -> None:
    """`niche_vs_celltype_heatmap_noannot.*` (level2), mirroring the annotated one."""
    ct = pd.crosstab(adata.obs[niche_key].astype(str),
                     adata.obs[ct_key].astype(str), normalize="index")
    ct = ct[[c for c in ORDER if c in ct.columns] +
            [c for c in ct.columns if c not in ORDER]]
    fig, ax = plt.subplots(figsize=(max(10, 0.5 * ct.shape[1]), 0.6 * ct.shape[0] + 2))
    sns.heatmap(ct, cmap="hot_r", annot=False, cbar_kws={"label": "row prop"},
                linewidths=0.3, linecolor="white", ax=ax)
    ax.set_xlabel(ct_key); ax.set_ylabel("niche")
    ax.set_title("Niche x cell-type (row-normalized)")
    fig.tight_layout()
    fig.savefig(f"{out_base}.pdf", bbox_inches="tight")
    fig.savefig(f"{out_base}.png", bbox_inches="tight", dpi=160)
    plt.close(fig)
    ct.round(4).to_csv(f"{out_base}.csv")
    print(f"  wrote {out_base}.{{pdf,png,csv}} (noannot)", flush=True)


def copy_independent(src: Path, dest: Path) -> None:
    """Copy annotation-independent figures + roi_flags from src to dest."""
    (dest / "figures").mkdir(parents=True, exist_ok=True)
    n_copy = n_skip = 0
    for f in sorted((src / "figures").glob("*")):
        if f.is_file() and not is_celltype_fig(f.name):
            shutil.copy2(f, dest / "figures" / f.name)
            n_copy += 1
        elif f.is_file():
            n_skip += 1
    if (src / "roi_flags").is_dir():
        shutil.copytree(src / "roi_flags", dest / "roi_flags", dirs_exist_ok=True)
    # niche-only per-core spatial dirs (no cell-type panel) are annotation-independent
    extra = [d for d in src.iterdir()
             if d.is_dir() and d.name.startswith("glom_niche_per_core")]
    for d in extra:
        shutil.copytree(d, dest / d.name, dirs_exist_ok=True)
    print(f"  copied {n_copy} annotation-independent figures "
          f"(skipped {n_skip} cell-type figures -> regenerated); roi_flags"
          + (f" + {len(extra)} niche-only per-core dir(s)" if extra else "") + " copied",
          flush=True)


def run_version(cfg: dict, skip_de: bool = False) -> None:
    adata_path, niche_key, res = cfg["adata"], cfg["niche_key"], cfg["res"]
    src, dest = cfg["src"], cfg["dest"]
    fig_dir = dest / "figures"
    de_dir = dest / "niche_marker_results"
    fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*70}\n[version] {src.name} -> {dest.name}\n{'='*70}", flush=True)
    print(f"[load] {adata_path}", flush=True)
    adata = sc.read_h5ad(adata_path)
    adata, n_nan = merge_level2(adata, key=CT_KEY)
    if n_nan:
        print(f"[WARN] {n_nan} cells unmatched for level2_refined (kept as NaN)", flush=True)

    # drive the reused downstream functions via module globals
    ds.NICHE_KEY = niche_key
    ds.LEIDEN_RES = res
    ds.CELL_TYPE_KEY = CT_KEY
    palette = level2_colour_dict(adata, CT_KEY)

    # 1. composition by level2_refined (fixed palette)
    ds.plot_niche_composition(
        adata, niche_key, CT_KEY,
        fig_dir / f"res_{res}_niche_composition_{CT_KEY}",
        title_suffix=f"  (n={adata.n_obs:,} cells)", colour_override=palette)
    # 2. niche x cell-type heatmaps (annotated + noannot)
    ds.heatmap_niche_vs_celltype(adata, fig_dir)
    heatmap_noannot(adata, niche_key, CT_KEY, fig_dir / "niche_vs_celltype_heatmap_noannot")
    # 3. niche annotation summary (dominant level2 per niche)
    ds.write_niche_annotation(adata, fig_dir)
    # 4. per-niche within-niche DE on level2_refined
    if skip_de:
        print("  [skip-de] within-niche DE skipped", flush=True)
    else:
        ds.per_niche_celltype_markers(adata, de_dir)

    # a lightweight annotation table so dest is usable without the 0.5 GB object
    adata.obs[["sample", "cell_id", "TMA_ID", niche_key, CT_KEY]].to_csv(
        dest / "cell_level2_and_niche.csv", index=False)

    # copy the annotation-independent figures + roi_flags
    copy_independent(src, dest)
    print(f"[done] {dest}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--version", choices=list(VERSIONS) + ["all"], default="all")
    ap.add_argument("--skip-de", action="store_true")
    ap.add_argument("--dest-root", type=str, default=str(NEW),
                    help="base dir for the *_level2 output folders (scratch for smoke)")
    args = ap.parse_args()
    root = Path(args.dest_root)
    todo = list(VERSIONS) if args.version == "all" else [args.version]
    for v in todo:
        cfg = dict(VERSIONS[v])
        cfg["dest"] = root / cfg["dest"].name
        run_version(cfg, skip_de=args.skip_de)
    print("\n[ALL DONE]", ", ".join(todo), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
