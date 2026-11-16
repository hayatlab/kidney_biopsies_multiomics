"""Copy the three NicheCompass review-round outputs (their h5ads, full
figure folders, and per-niche marker results) to the shared Kramann_lab
results directory in a per-run layout, and write a README.md describing
each subfolder.

Three runs are handled:

  - `manual_glom_niche`  — NicheCompass on cells with manually-annotated
    `in_glom == True` (polygon-derived for the 2 new Region_1 samples).
  - `overall_with_glom`  — full 5-region NicheCompass with `in_glom`
    merged into obs (re-run of the parent pipeline with the new
    annotations).
  - `in_glom_niche`      — NicheCompass on cells assigned to the parent
    "glom niche" (latent_leiden_0.4 == 5; the "niche 9" in earlier
    numbering). Pre-existing run at
    `niche_compass_results/artifacts/in_glom_niche/27042026_060952/`.

Destination layout:

  $KB_RESULTS/xenium/leon/new_data/
    README.md
    <run>/
      xenium_with_niches_<run>.h5ad
      figures/                 # everything from <run_dir>/figures/ + glom_postprocess/
      niche_marker_results/    # per-niche DE dotplots/heatmaps when present
      cell_boundary_plots_per_core/    # PDF+PNG per TMA core when present
      glom_niche_per_core/             # per-core glom subniche plots when present

Idempotent. Re-running overwrites file contents. Pre-existing files in
the parent `…/leon/` directory are not touched.
"""
from __future__ import annotations
from kbpaths import KB_ANALYSIS, KB_PROJECT, KB_RESULTS  # repository path config (config/kbpaths.py)

import shutil
import sys
from datetime import datetime
from pathlib import Path
from textwrap import dedent

import anndata as ad
import pandas as pd


DEST = Path(f"{KB_RESULTS}/xenium/leon/new_data")

RUNS = {
    "manual_glom_niche": {
        "base": Path(
            f"{KB_PROJECT}/niche_compass_results/artifacts/"
            "manual_glom_niche"
        ),
        "fixed_run_dir": None,  # auto-pick latest
        "adata_name": "xenium_with_niches_manual_glom_niche.h5ad",
        "title": "NicheCompass on cells with manually-annotated `in_glom == True`",
        "body": (
            "Filter: `obs['in_glom'] == True` carried from the parent run for "
            "the 3 original ROI samples + polygon-derived annotations for the "
            "2 newly-added Region_1 samples (Combi 1734 R1, FSGS 42348 R1). "
            "Same TMA blacklist as the parent run plus the Combi D2 IgAN QC-fail "
            "core. NicheCompass hyperparameters identical to the previous "
            "`in_glom_niche` run (`nichecompass_glom_subniche.py`); leiden "
            "resolution **0.6**."
        ),
    },
    "overall_with_glom": {
        "base": Path(
            f"{KB_PROJECT}/niche_compass_results/artifacts/"
            "overall_with_glom"
        ),
        "fixed_run_dir": None,
        "adata_name": "xenium_with_niches_overall_with_glom.h5ad",
        "title": "Full 5-region NicheCompass with `in_glom` overlaid in `obs`",
        "body": (
            "Re-runs the same overall pipeline as `xenium_with_niches.h5ad` "
            "(`nichecompass_5regions.py`, leiden **0.4**) but with `in_glom` "
            "merged into `obs` before training: True values come from the "
            "existing parent adata for the 3 original ROI samples and from "
            "the polygon CSVs in `review_analyses/in_glom_for_new_samples/` "
            "for the 2 new Region_1 samples. Same TMA blacklist + Combi D2 IgAN "
            "exclusion; identical hyperparameters. Use the included "
            "`figures/niche_vs_in_glom.csv` / `figures/niche9_match_res_0.4.csv` to "
            "see which leiden cluster is the analog of the original niche 9."
        ),
    },
    "in_glom_niche": {
        "base": Path(
            f"{KB_PROJECT}/niche_compass_results/artifacts/"
            "in_glom_niche"
        ),
        # pin the existing run (no re-training needed)
        "fixed_run_dir": Path(
            f"{KB_PROJECT}/niche_compass_results/artifacts/"
            "in_glom_niche/27042026_060952"
        ),
        "adata_name": "xenium_with_niches_in_glom_niche.h5ad",
        "title": "NicheCompass on the parent **glom niche** (latent_leiden_0.4 == 5)",
        "body": (
            "Defines the glom subset by the parent overall NicheCompass model's "
            "own clustering — cells assigned to niche 5 (= the niche enriched "
            "for `in_glom == True`; called \"niche 9\" in earlier numbering). "
            "Trained with `nichecompass_glom_subniche.py` (leiden **0.6**), "
            "same hyperparameters as the other glom runs. This is the *model-"
            "driven* glom subset, complementary to the *annotation-driven* one "
            "in `manual_glom_niche/`."
        ),
    },
}

# Patterns to mirror from each run dir into <dest>/<run>/...
COPY_GLOBS = [
    # all files in <run_dir>/figures/ (top level only)
    {"src_subdir": "figures",           "dst_subdir": "figures",
     "patterns": ["*.pdf", "*.png", "*.svg", "*.csv"]},
    # all files in <run_dir>/glom_postprocess/ (top level only)
    {"src_subdir": "glom_postprocess",  "dst_subdir": "figures",
     "patterns": ["*.pdf", "*.png", "*.svg", "*.csv"]},
]
# Whole-directory mirrors (only copied if the source dir exists)
COPY_DIRS = [
    # per-niche DE markers — both naming conventions
    ("niche_marker_results",               "niche_marker_results"),
    ("niche_marker_results_compass_niche", "niche_marker_results"),
    # per-core cell-boundary plots (from plot_niche_cell_boundaries_*.py)
    ("cell_boundary_plots/per_core",       "cell_boundary_plots_per_core"),
    # per-core glom-niche plots (from plot_glom_niches_per_core.py); the
    # in_glom_niche source stores them under figures/, so the in_glom_niche
    # entry below picks them up via figures/ glob as well.
    ("figures/glom_niche_per_core",        "glom_niche_per_core"),
    ("figures/glom_niche_per_core_niche9", "glom_niche_per_core_niche9"),
    # roi_flags from downstream
    ("figures/roi_flags",                  "roi_flags"),
]


def latest_run_dir(base: Path) -> Path | None:
    if not base.exists():
        return None
    runs = [
        p for p in base.iterdir()
        if p.is_dir() and (
            (p / "adata.h5ad").exists() or (p / "xenium_with_niches.h5ad").exists()
        )
    ]
    if not runs:
        return None
    runs.sort(key=lambda p: p.stat().st_mtime)
    return runs[-1]


def pick_source_adata(run_dir: Path, adata_name: str) -> Path:
    """Prefer the DE-enriched adata written by
    `nichecompass_5regions_downstream.py` (named either
    `xenium_with_niches.h5ad` for the new runs or `<adata_name>` for the
    in_glom_niche run) over the bare training-output `adata.h5ad` if it
    exists — it carries the same data plus the DE results in uns.
    """
    candidates = [
        run_dir / adata_name,            # in_glom_niche layout (out-h5ad-name was passed)
        run_dir / "xenium_with_niches.h5ad",   # new-runs default downstream output
        run_dir / "adata.h5ad",          # plain training output
    ]
    for c in candidates:
        if c.is_file():
            return c
    raise FileNotFoundError(f"no adata candidates exist under {run_dir}")


def adata_summary(path: Path) -> dict:
    a = ad.read_h5ad(path)
    info: dict = {"n_obs": int(a.n_obs), "n_vars": int(a.n_vars)}
    info["samples"] = sorted(a.obs["sample"].astype(str).unique().tolist()) \
        if "sample" in a.obs.columns else []
    info["tma_count"] = int(a.obs["TMA_ID"].astype(str).nunique()) \
        if "TMA_ID" in a.obs.columns else None
    info["conditions"] = (
        a.obs["condition"].astype(str).value_counts().to_dict()
        if "condition" in a.obs.columns else {}
    )
    niche_keys = [c for c in a.obs.columns if c.startswith("latent_leiden_")]
    info["niche_keys"] = niche_keys
    info["niche_counts"] = {
        k: a.obs[k].astype(str).value_counts().to_dict() for k in niche_keys
    }
    if "in_glom" in a.obs.columns:
        ig = a.obs["in_glom"].astype(bool)
        info["in_glom_true"] = int(ig.sum())
        info["in_glom_pct"] = float(100.0 * ig.mean())
    else:
        info["in_glom_true"] = None
        info["in_glom_pct"] = None
    info["niche_vs_in_glom"] = None
    if niche_keys and "in_glom" in a.obs.columns:
        nk = niche_keys[0]
        ct = (
            a.obs.groupby(nk, observed=True)["in_glom"]
            .agg(n_total="count", n_in_glom="sum")
        )
        ct["pct_in_glom"] = 100.0 * ct["n_in_glom"] / ct["n_total"]
        total = int(ct["n_in_glom"].sum())
        ct["share_of_in_glom"] = 100.0 * ct["n_in_glom"] / max(total, 1)
        ct = ct.sort_values("pct_in_glom", ascending=False)
        info["niche_vs_in_glom"] = {
            "key": nk,
            "table": ct.head(8).round(2).to_dict(orient="index"),
            "top_niche": str(ct.index[0]),
            "top_pct": float(ct.iloc[0]["pct_in_glom"]),
            "top_share": float(ct.iloc[0]["share_of_in_glom"]),
        }
    return info


def copy_dir_contents(src: Path, dst: Path, patterns: list[str]) -> int:
    """Copy files matching patterns from src (top-level only) into dst."""
    if not src.is_dir():
        return 0
    dst.mkdir(parents=True, exist_ok=True)
    n = 0
    for pat in patterns:
        for f in src.glob(pat):
            if f.is_file():
                shutil.copy2(f, dst / f.name)
                n += 1
    return n


def mirror_dir(src: Path, dst: Path) -> int:
    """Mirror src → dst (top-level files only). No-op + no dst created
    if src is empty or missing."""
    if not src.is_dir():
        return 0
    files = [f for f in sorted(src.iterdir()) if f.is_file()]
    if not files:
        return 0
    dst.mkdir(parents=True, exist_ok=True)
    for f in files:
        shutil.copy2(f, dst / f.name)
    return len(files)


def copy_run(tag: str, cfg: dict, run_dir: Path, dst_root: Path) -> dict:
    """Copy adata + figures + DE + per-core plots for one run."""
    dst = dst_root / tag
    dst.mkdir(parents=True, exist_ok=True)
    inventory: dict = {
        "adata": None, "figures": {}, "extras": {},
    }

    # adata — prefer DE-enriched output of downstream
    try:
        src_adata = pick_source_adata(run_dir, cfg["adata_name"])
    except FileNotFoundError as e:
        print(f"  [{tag}] WARN: {e}", flush=True)
        src_adata = None
    if src_adata is None:
        pass
    else:
        dst_adata = dst / cfg["adata_name"]
        shutil.copy2(src_adata, dst_adata)
        inventory["adata"] = {
            "path": dst_adata, "size_mb": dst_adata.stat().st_size / 1e6,
            "src_name": src_adata.name,
        }
        print(f"  [{tag}] adata ({src_adata.name}) -> {dst_adata.name} "
              f"({inventory['adata']['size_mb']:.1f} MB)", flush=True)

    # figures + glom_postprocess (merge under figures/)
    figures_dst = dst / "figures"
    n_files = 0
    for spec in COPY_GLOBS:
        src = run_dir / spec["src_subdir"]
        target = dst / spec["dst_subdir"]
        n = copy_dir_contents(src, target, spec["patterns"])
        n_files += n
        if n:
            print(f"  [{tag}] {spec['src_subdir']}/ -> {spec['dst_subdir']}/ "
                  f"({n} files)", flush=True)
    if figures_dst.is_dir():
        pdfs = len(list(figures_dst.glob("*.pdf")))
        pngs = len(list(figures_dst.glob("*.png")))
        csvs = len(list(figures_dst.glob("*.csv")))
        svgs = len(list(figures_dst.glob("*.svg")))
        inventory["figures"] = {"total": n_files, "pdf": pdfs, "png": pngs,
                                "csv": csvs, "svg": svgs}

    # extras: directory mirrors
    for src_rel, dst_name in COPY_DIRS:
        src = run_dir / src_rel
        if src.is_dir():
            n = mirror_dir(src, dst / dst_name)
            inventory["extras"][dst_name] = n
            print(f"  [{tag}] {src_rel}/ -> {dst_name}/ ({n} files)", flush=True)

    return inventory


def cleanup_old_toplevel_layout() -> None:
    """Remove leftovers from earlier layouts so the destination matches
    the current per-run structure exactly.

    - Any bare .h5ad at the top level of DEST (= DEST/<name>.h5ad)
    - Any non-adata files at the top level of each DEST/<tag>/ subfolder
      (earlier copy versions wrote heatmaps/dotplots there directly; the
      current layout puts them under DEST/<tag>/figures/).
    """
    # bare .h5ads at the top
    for cfg in RUNS.values():
        stale = DEST / cfg["adata_name"]
        if stale.is_file():
            print(f"[cleanup] removing stale top-level {stale}", flush=True)
            stale.unlink()

    # stale figure files at the top of each per-run subdir
    for tag, cfg in RUNS.items():
        sub = DEST / tag
        if not sub.is_dir():
            continue
        keep_name = cfg["adata_name"]
        for f in sub.iterdir():
            if f.is_file() and f.name != keep_name:
                print(f"[cleanup] removing stale {f}", flush=True)
                f.unlink()


def render_readme(meta: dict) -> str:
    today = datetime.now().strftime("%Y-%m-%d")
    parts: list[str] = []
    parts.append(dedent(f"""
        # NicheCompass review-round outputs

        _Last updated: {today}  ·  copied by `review_analyses/copy_results_to_kramann.py`_

        This directory holds **three** NicheCompass result subfolders produced
        during the review round, alongside the pre-existing
        `xenium_with_niches*.h5ad` files in the parent directory
        (`{KB_RESULTS}/xenium/leon/`).

        ## The three runs

        | folder | niche definition | leiden | cells |
        |---|---|---|---|
        | `overall_with_glom/` | full 5-region NicheCompass (no glom filtering); `in_glom` merged into obs | 0.4 | 258,375 |
        | `manual_glom_niche/` | annotation-driven glom subset (`obs['in_glom'] == True`) | 0.6 | 20,081 |
        | `in_glom_niche/`     | model-driven glom subset (parent niche 5 = the niche enriched for in_glom; "niche 9" in earlier numbering) | 0.6 | ~4,000 |

        ## Layout

        ```
        new_data/
        ├── README.md                                (this file)
        ├── overall_with_glom/
        │   ├── xenium_with_niches_overall_with_glom.h5ad
        │   ├── figures/                             (PDF+PNG figures + CSV tables)
        │   ├── niche_marker_results/                (per-niche DE dotplots/heatmaps when run)
        │   └── (cell-boundary plots if generated)
        ├── manual_glom_niche/
        │   ├── xenium_with_niches_manual_glom_niche.h5ad
        │   ├── figures/
        │   ├── niche_marker_results/
        │   └── cell_boundary_plots_per_core/        (76 TMA cores × {{pdf,png}})
        └── in_glom_niche/
            ├── xenium_with_niches_in_glom_niche.h5ad
            ├── figures/
            ├── niche_marker_results/
            └── glom_niche_per_core/                 (per-core glom niche plots)
        ```

        ## Per-run summaries
    """).strip())
    parts.append("")

    for tag in ("overall_with_glom", "manual_glom_niche", "in_glom_niche"):
        if tag not in meta:
            continue
        m = meta[tag]
        info = m["info"]
        inv = m["inventory"]
        cfg = RUNS[tag]
        parts.append(f"### `{tag}/` — {cfg['title']}")
        parts.append("")
        parts.append(cfg["body"])
        parts.append("")
        primary_nk = info["niche_keys"][0] if info["niche_keys"] else None
        n_clusters = len(info["niche_counts"].get(primary_nk, {})) if primary_nk else 0
        in_glom_note = "  (ALL cells)" if (info["in_glom_pct"] or 0) > 99.9 else ""
        parts.append(dedent(f"""
            **adata** — `{cfg['adata_name']}`
            - Source run dir: `{m['source']}`
            - File size: {inv['adata']['size_mb']:.1f} MB
            - Shape: **{info['n_obs']:,} cells × {info['n_vars']} genes**
            - Samples ({len(info['samples'])}): {', '.join(s.replace('XETG00229__', '') for s in info['samples'])}
            - TMA cores: {info['tma_count']}
            - Condition counts: {info['conditions']}
            - `in_glom == True`: {info['in_glom_true']:,} cells ({info['in_glom_pct']:.1f}%){in_glom_note}
            - Niche key: `{primary_nk or 'n/a'}` with {n_clusters} clusters
        """).strip())
        parts.append("")
        fig_inv = inv.get("figures", {})
        if fig_inv:
            parts.append(
                f"**`figures/`** — {fig_inv.get('total', 0)} files "
                f"({fig_inv.get('pdf', 0)} PDF, {fig_inv.get('png', 0)} PNG, "
                f"{fig_inv.get('svg', 0)} SVG, {fig_inv.get('csv', 0)} CSV)"
            )
            parts.append("")
        ex = inv.get("extras", {})
        if ex:
            parts.append("**Extras:**")
            for name, n in ex.items():
                parts.append(f"- `{name}/` — {n} files")
            parts.append("")

        nvi = info.get("niche_vs_in_glom")
        if nvi is not None and tag in ("overall_with_glom",):
            parts.append("**Niche-9 recovery (this run):**")
            parts.append("")
            parts.append("| niche | n_total | n_in_glom | pct_in_glom | share_of_in_glom |")
            parts.append("|---|---:|---:|---:|---:|")
            for niche, row in nvi["table"].items():
                parts.append(
                    f"| {niche} | {int(row['n_total']):,} | {int(row['n_in_glom']):,} "
                    f"| {row['pct_in_glom']:.2f} | {row['share_of_in_glom']:.2f} |"
                )
            parts.append("")
            parts.append(
                f"Top niche by `pct_in_glom`: `{nvi['top_niche']}` "
                f"({nvi['top_pct']:.1f}% of its cells in-glom, "
                f"{nvi['top_share']:.1f}% of all in_glom cells)."
            )
            parts.append("")

    parts.append("## Reproduction")
    parts.append("")
    parts.append(dedent(f"""
        Scripts (all under `{KB_ANALYSIS}/`):
        - `annotate_in_glom_new_samples.py` — polygon-in-point for the 2 new
          Region_1 samples, writes `in_glom_new_samples.parquet`.
        - `nichecompass_manual_glom_niche.py` (+ `.sh`) — trains the
          `manual_glom_niche` h5ad.
        - `nichecompass_overall_with_glom.py` (+ `.sh`) — trains the
          `overall_with_glom` h5ad.
        - `overall_with_glom_finalize.py` — recovery script for the
          post-training leiden + niche-9 diagnostic if the training job hits
          its SLURM time limit (used for job 321201 / 325550).
        - `glom_postprocess.py` (+ `.sh` / `_auto.sh`) — niche × cell-type
          heatmap (no in-cell numbers) + hallmark pathway dotplot + dedup
          clustermap + niche-vs-in_glom enrichment CSV.
        - `nichecompass_5regions_downstream.py` — standard downstream battery
          (composition / condition / EGFR / IFTA enrichment + GP heatmap +
          hallmark / reactome + per-niche DE).
        - `plot_niche_cell_boundaries_manual_glom.py` (+ `.sh`) — per-core
          two-panel PDFs (niches + cell types, grey nucleus boundaries).
        - `plot_glom_niches_per_core.py` (+ `.sh`) — per-core glom subniche
          spatial plots.
        - `copy_results_to_kramann.py` (+ `.sh`) — this script. Idempotent —
          re-running picks up whatever fresh figures the downstream / per-core
          jobs produced and refreshes this README.

        Run order: annotate → train (both new runs) → glom_postprocess (all 3) →
        nichecompass_5regions_downstream (all 3) → per-core plotters → copy.
    """).strip())
    parts.append("")
    return "\n".join(parts)


def main() -> int:
    DEST.mkdir(parents=True, exist_ok=True)
    print(f"[dest] {DEST}", flush=True)
    cleanup_old_toplevel_layout()

    meta: dict = {}
    for tag, cfg in RUNS.items():
        run_dir = cfg["fixed_run_dir"] or latest_run_dir(cfg["base"])
        if run_dir is None or not (run_dir / "adata.h5ad").exists():
            print(f"[warn] skipping {tag}: no adata.h5ad under {cfg['base']}",
                  flush=True)
            continue
        print(f"\n[{tag}] run_dir={run_dir}", flush=True)
        info = adata_summary(run_dir / "adata.h5ad")
        print(f"  shape={info['n_obs']:,} x {info['n_vars']}  "
              f"niches={info['niche_keys']}  in_glom_true={info['in_glom_true']}",
              flush=True)
        inv = copy_run(tag, cfg, run_dir, DEST)
        meta[tag] = {"source": str(run_dir), "info": info, "inventory": inv}

    if not meta:
        print("[error] nothing copied — no run dirs ready", flush=True)
        return 1

    readme = DEST / "README.md"
    readme.write_text(render_readme(meta))
    print(f"\n[readme] wrote {readme}", flush=True)
    print(f"\n=== DONE ===", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
