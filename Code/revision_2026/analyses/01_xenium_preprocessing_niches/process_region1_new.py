"""Phase A: annotate the two new Xenium Region 1 h5ads with TMA core
identity + whitelist metadata, produce QC artifacts for review.

No label transfer, no merging with the existing aggregate — that is
Phase B, gated on Leon's approval of the outputs from this script.
"""
from __future__ import annotations
from kbpaths import KB_ANALYSIS, KB_XENIUM_RAW  # repository path config (config/kbpaths.py)

import argparse
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from shapely.geometry import Polygon
from shapely.strtree import STRtree

REVIEW_DIR = Path(f"{KB_ANALYSIS}")
DATA_DIR = Path(f"{KB_XENIUM_RAW}/h5ad_files_only_transcript")
WHITELIST_CSV = REVIEW_DIR / "samples_to_use.csv"

METADATA_COLS = [
    "Disease",
    "H.number",
    "patient_id",
    "TMA Block",
    "source of sample",
    "IFTA",
    "eGFR_biopsy",
    "Gender",
    "age at biopsy",
    "Etnicity",
]

REGIONS = [
    dict(
        label="41734_R1",
        slide_id=41734,
        region=1,
        h5ad=DATA_DIR / "output-XETG00229__0041734__Region_1__20250416__104148.h5ad",
        coords_csv=REVIEW_DIR / "coordinates Combi 0041734 region 1.csv",
        sample_name="XETG00229__0041734__Region_1__20250416__104148",
        out_h5ad=REVIEW_DIR / "annotated_0041734_Region_1.h5ad",
    ),
    dict(
        label="42348_R1",
        slide_id=42348,
        region=1,
        h5ad=DATA_DIR / "output-XETG00229__0042348__Region_1__20250416__104148.h5ad",
        coords_csv=REVIEW_DIR / "coordinates FSGS 0042348 region 1.csv",
        sample_name="XETG00229__0042348__Region_1__20250416__104148",
        out_h5ad=REVIEW_DIR / "annotated_0042348_Region_1.h5ad",
    ),
]


def _norm(s: str) -> str:
    return " ".join(str(s).replace("\xa0", " ").split()).strip()


def parse_polygons(coords_csv: Path) -> dict[str, Polygon]:
    df = pd.read_csv(coords_csv, comment="#")
    polys: dict[str, Polygon] = {}
    for sel, grp in df.groupby("Selection", sort=False):
        pts = list(zip(grp["X"].to_numpy(), grp["Y"].to_numpy()))
        if len(pts) < 3:
            print(f"  [warn] selection {sel!r} has <3 vertices, skipping")
            continue
        polys[_norm(sel)] = Polygon(pts)
    return polys


def assign_cells_to_cores(
    xy: np.ndarray, polys: dict[str, Polygon]
) -> np.ndarray:
    """Return an array of core names (str) with NaN-equivalent '' for
    unassigned cells. Uses STRtree for O((N+M) log M) point-in-polygon."""
    names = list(polys.keys())
    geoms = [polys[n] for n in names]
    tree = STRtree(geoms)

    # shapely 2 STRtree.query(points, predicate="within") returns
    # a 2xK int array (input_idx, tree_idx) of matches.
    from shapely.geometry import Point

    pts = [Point(x, y) for x, y in xy]
    hits = tree.query(pts, predicate="within")
    assigned = np.full(len(pts), "", dtype=object)
    # hits shape: (2, K) — row 0 = point index, row 1 = polygon index
    for pi, gi in zip(hits[0], hits[1]):
        # keep first hit (polygons are non-overlapping)
        if assigned[pi] == "":
            assigned[pi] = names[gi]
    return assigned


def load_whitelist(slide_id: int, region: int) -> pd.DataFrame:
    df = pd.read_csv(WHITELIST_CSV, sep=";", engine="python", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    mask = (df["Xenium slide ID"].astype(str) == str(slide_id)) & (
        df["Region"].astype(str) == str(region)
    )
    sub = df.loc[mask].copy()
    # Normalise whitespace: the Excel-originated CSV uses non-breaking
    # spaces in "Sample name Xenium", while the Xenium-Explorer polygon
    # export uses regular spaces. Matching requires normalisation on both
    # sides (polygons are already normalised in parse_polygons).
    sub["Sample name Xenium"] = sub["Sample name Xenium"].map(_norm)
    # Deduplicate rows on Sample name Xenium (CSV can have multi-line
    # quoted descriptions that inflate row count but same sample key)
    sub = sub.drop_duplicates(subset=["Sample name Xenium"], keep="first")
    return sub


def plot_qc(adata: ad.AnnData, label: str, out_dir: Path) -> None:
    xy = np.asarray(adata.obsm["spatial"])
    # 1) Spatial scatter coloured by core
    fig, ax = plt.subplots(figsize=(10, 10))
    cores = adata.obs["TMA_ID"].astype(str).to_numpy()
    unique = sorted(set(cores))
    # consistent colours; "Unassigned" in grey
    palette = plt.cm.tab20(np.linspace(0, 1, max(len(unique), 1)))
    colour_map = {c: palette[i % len(palette)] for i, c in enumerate(unique)}
    if "Unassigned" in colour_map:
        colour_map["Unassigned"] = (0.7, 0.7, 0.7, 0.4)
    colours = np.array([colour_map[c] for c in cores])
    ax.scatter(xy[:, 0], xy[:, 1], c=colours, s=0.3, linewidths=0)
    ax.set_aspect("equal")
    ax.set_title(f"{label}: cells coloured by TMA core")
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    # legend with small swatches
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=colour_map[c], markersize=6, label=c)
        for c in unique
    ]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=6, frameon=False)
    plt.tight_layout()
    fig.savefig(out_dir / f"qc_{label}_spatial_by_core.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # 2) Per-core cell counts
    counts = adata.obs["TMA_ID"].astype(str).value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(max(8, 0.3 * len(counts)), 5))
    ax.bar(range(len(counts)), counts.to_numpy())
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(counts.index, rotation=90, fontsize=7)
    ax.axhline(1000, color="green", lw=0.5, ls="--")
    ax.axhline(500, color="orange", lw=0.5, ls="--")
    ax.set_ylabel("cells")
    ax.set_title(f"{label}: cells per core (dashed: 500 / 1000 thresholds)")
    plt.tight_layout()
    fig.savefig(out_dir / f"qc_{label}_cells_per_core.png", dpi=200)
    plt.close(fig)


def process_region(cfg: dict, subsample: int | None = None) -> dict:
    label = cfg["label"]
    print(f"\n=== processing {label} ===", flush=True)

    adata = ad.read_h5ad(cfg["h5ad"])
    if subsample and adata.n_obs > subsample:
        rng = np.random.default_rng(0)
        idx = rng.choice(adata.n_obs, size=subsample, replace=False)
        adata = adata[idx].copy()
        print(f"  subsampled to {adata.n_obs} cells for smoke test", flush=True)
    else:
        print(f"  loaded {adata.n_obs} cells × {adata.n_vars} genes", flush=True)

    # 1) Polygon parse + core assignment on all cells
    polys = parse_polygons(cfg["coords_csv"])
    print(f"  parsed {len(polys)} core polygons from {cfg['coords_csv'].name}", flush=True)
    xy = np.asarray(adata.obsm["spatial"])
    core_assignment = assign_cells_to_cores(xy, polys)
    n_assigned = int((core_assignment != "").sum())
    print(f"  assigned {n_assigned}/{len(core_assignment)} cells to a core "
          f"({100*n_assigned/len(core_assignment):.1f}%)", flush=True)

    adata.obs["core_raw"] = pd.Categorical(
        np.where(core_assignment == "", "Unassigned", core_assignment)
    )

    # 2) Whitelist filter + metadata join
    whitelist = load_whitelist(cfg["slide_id"], cfg["region"])
    wl_names = set(whitelist["Sample name Xenium"].astype(str))
    all_core_names = set(polys.keys())
    print(f"  whitelist: {len(wl_names)} samples; "
          f"polygons present: {len(all_core_names)}; "
          f"whitelist∩polygons: {len(wl_names & all_core_names)}", flush=True)
    missing_polys = wl_names - all_core_names
    if missing_polys:
        print(f"  [warn] whitelisted samples without polygon: {sorted(missing_polys)}", flush=True)
    extra_polys = all_core_names - wl_names
    if extra_polys:
        print(f"  [info] polygons NOT in whitelist (will be dropped): {sorted(extra_polys)}", flush=True)

    keep_mask = adata.obs["core_raw"].astype(str).isin(wl_names).to_numpy()
    n_drop_unassigned = int((core_assignment == "").sum())
    n_drop_notwhitelist = int((~keep_mask).sum() - n_drop_unassigned)
    print(f"  dropping {n_drop_unassigned} unassigned + {n_drop_notwhitelist} not-in-whitelist "
          f"→ keeping {int(keep_mask.sum())} cells", flush=True)

    kept = adata[keep_mask].copy()

    # Join metadata on core (TMA_ID convention: per-core sample name)
    meta = whitelist[["Sample name Xenium", *METADATA_COLS]].copy()
    meta = meta.rename(columns={"Sample name Xenium": "TMA_ID"})
    kept.obs["TMA_ID"] = kept.obs["core_raw"].astype(str).values
    del kept.obs["core_raw"]
    kept.obs = kept.obs.merge(meta, on="TMA_ID", how="left").set_index(kept.obs.index)

    # Bookkeeping to match reference schema
    kept.obs["sample"] = cfg["sample_name"]
    # condition mapping (per Leon's decision 2026-04-21): collapse
    # iFSGS/MCD into FSGS; other Disease values flow through unchanged.
    # Non-whitelisted categories (Alport / TBMD / C3GN / unknown) are
    # already dropped at the whitelist step, so we only see
    # DN / IgAN / FSGS / iFSGS/MCD on the kept cells.
    condition_map = {"iFSGS/MCD": "FSGS"}
    kept.obs["condition"] = (
        kept.obs["Disease"].astype(str).replace(condition_map)
    )

    # h5ad write requires each column to have a consistent dtype; the
    # joined metadata columns contain NaN + strings from the CSV. Use
    # plain object dtype with empty-string NA fill (anndata <0.11 does
    # not serialise pandas nullable StringArray by default).
    for col in [*METADATA_COLS, "TMA_ID", "sample", "condition"]:
        if col in kept.obs.columns:
            s = kept.obs[col].astype(object)
            kept.obs[col] = s.where(s.notna(), "").astype(str)

    # Re-index to match the existing pattern: "<sample>-<cellid>"
    kept.obs_names = [f"{cfg['sample_name']}-{cid}" for cid in kept.obs["cell_id"].astype(str)]

    # 3) QC + save
    cfg["out_h5ad"].parent.mkdir(parents=True, exist_ok=True)
    plot_qc(kept, label, cfg["out_h5ad"].parent)
    kept.write_h5ad(cfg["out_h5ad"], compression="gzip")
    print(f"  wrote {cfg['out_h5ad']} ({kept.n_obs} cells × {kept.n_vars} genes)", flush=True)

    # Summary row
    summary = dict(
        label=label,
        slide_id=cfg["slide_id"],
        region=cfg["region"],
        total_cells_in=adata.n_obs,
        cells_unassigned=n_drop_unassigned,
        cells_not_whitelisted=n_drop_notwhitelist,
        cells_kept=int(kept.n_obs),
        polygons=len(polys),
        whitelist_samples=len(wl_names),
        whitelist_samples_with_polygon=len(wl_names & all_core_names),
        whitelist_samples_missing_polygon=sorted(missing_polys),
    )
    # Per-core counts + disease
    per_core = (
        kept.obs.groupby("TMA_ID", observed=True)
        .agg(cells=("TMA_ID", "size"), disease=("Disease", "first"))
        .reset_index()
        .sort_values("cells", ascending=False)
    )
    per_core["label"] = label
    return dict(summary=summary, per_core=per_core)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true",
                    help="subsample each region to 500 cells for a quick end-to-end test")
    args = ap.parse_args()

    REVIEW_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    per_core_all = []
    for cfg in REGIONS:
        res = process_region(cfg, subsample=500 if args.smoke else None)
        rows.append(res["summary"])
        per_core_all.append(res["per_core"])

    tag = "_smoke" if args.smoke else ""
    summary_df = pd.DataFrame(rows)
    summary_path = REVIEW_DIR / f"phaseA_summary{tag}.csv"
    summary_df.to_csv(summary_path, index=False)
    print("\n=== summary ===")
    print(summary_df.to_string(index=False))

    per_core_df = pd.concat(per_core_all, ignore_index=True)
    per_core_path = REVIEW_DIR / f"phaseA_per_core_counts{tag}.csv"
    per_core_df.to_csv(per_core_path, index=False)

    low = per_core_df.query("cells < 500")
    if len(low):
        print("\n=== cores with <500 cells (flag for Leon) ===")
        print(low.to_string(index=False))
    else:
        print("\nAll kept cores have >=500 cells.")

    print(f"\nwrote {summary_path} and {per_core_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
