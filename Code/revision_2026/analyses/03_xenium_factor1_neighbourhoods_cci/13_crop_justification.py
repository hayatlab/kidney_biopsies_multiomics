#!/usr/bin/env python
"""Cropping justification for the ribbon panels (Fig 3h, ED 3f, and the NF-kB ED 3f variant).

Vivien HARD-CODED the ribbon crop to Factor 2's range ([-2, 5] for Fig 3h, [-4, 8] for
ED 3f) with no stated criterion. Factor 1 is on a different scale, so a criterion is needed
rather than a copied number. This script draws the evidence Leon asked for:

  A  how many 2-hop neighbourhoods sit at each score value (log y)
  B  the same on a linear axis, so the bulk is readable
  C  how many distinct TMA cores contribute at each score value

C is the part that actually settles it: a spline fitted into a tail that only two cores
reach is describing those two cores, not the cohort. Candidate crops are the largest
contiguous run of bins holding at least N neighbourhoods, evaluated at several thresholds
and drawn together with the 1st-99th percentile crop.

The defaults reproduce the Factor-1 run exactly (same file names). `--score` points the
same rule at another per-cell score (2026-09-10: NF-kB). Bin width matters for a
"neighbourhoods per bin" rule, so `--bins-per-1-99` sets it from the score's own spread:
Factor 1 at 0.25 gives 46 bins across its 1st-99th percentile, and 46 gives another score
the same resolution.

Nothing is cropped here -- this script only produces the evidence and the table.
Human GRCh38, HGNC symbols. Seed 42.
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import argparse, json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

HERE = Path(f"{KB_ANALYSIS}/vivien_takeover")
OUT, FIG = HERE / "outputs", HERE / "figures"
SEED = 42
THRESHOLDS = [100, 250, 500, 1000, 2000]
USED_THRESHOLD = 500          # the rule 05_xenium_figures.py applies (Leon, 2026-09-07)


def longest_run(counts, edges, thr):
    """Largest contiguous stretch of bins each holding >= thr neighbourhoods."""
    ok = counts >= thr
    best, i = (0, None), 0
    while i < len(ok):
        if ok[i]:
            j = i
            while j + 1 < len(ok) and ok[j + 1]:
                j += 1
            if j - i + 1 > best[0]:
                best = (j - i + 1, (i, j))
            i = j + 1
        else:
            i += 1
    if best[1] is None:
        return None
    i, j = best[1]
    return float(edges[i]), float(edges[j + 1]), float(counts[i:j + 1].sum() / counts.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--score", default="Factor1_ulm")
    ap.add_argument("--scores-parquet", default=str(OUT / "xenium_factor_scores.parquet"))
    ap.add_argument("--stem", default="factor1",
                    help="figure = crop_justification_<stem>.pdf; tables/JSON get _<stem> "
                         "unless stem is factor1 (keeps the original file names)")
    ap.add_argument("--bin-width", type=float, default=0.25)
    ap.add_argument("--bins-per-1-99", type=float, default=None,
                    help="if set, bin width = (99th - 1st percentile) / this")
    args = ap.parse_args()
    np.random.seed(SEED)
    sns.set_style("white")
    outdir = FIG / "xenium"; outdir.mkdir(parents=True, exist_ok=True)
    score = args.score
    sfx = "" if args.stem == "factor1" else f"_{args.stem}"

    obs = pd.read_parquet(args.scores_parquet,
                          columns=["sample_cell_id", score, "TMA_ID"]).set_index("sample_cell_id")
    compo = pd.read_csv(OUT / "neighborhoods" / "nh_2hop_compo_10.csv", index_col=0, usecols=[0, 1])
    df = obs.reindex(compo.index)
    n_nan = int(df[score].isna().sum())
    print(f"[load] {len(df)} 2-hop neighbourhoods, {n_nan} ({n_nan/len(df)*100:.3f}%) "
          f"without a {score} value -> dropped", flush=True)
    df = df.dropna(subset=[score])

    q_lo, q_hi = df[score].quantile([0.01, 0.99])
    w = (float(q_hi - q_lo) / args.bins_per_1_99) if args.bins_per_1_99 else args.bin_width
    edges = np.arange(np.floor(df[score].min()), np.ceil(df[score].max()) + w, w)
    counts, _ = np.histogram(df[score], bins=edges)
    binid = np.clip(np.digitize(df[score], edges) - 1, 0, len(edges) - 2)
    ncores = (pd.DataFrame({"b": binid, "c": df["TMA_ID"].to_numpy()})
                .groupby("b")["c"].nunique().reindex(range(len(edges) - 1), fill_value=0).to_numpy())
    centres = (edges[:-1] + edges[1:]) / 2
    print(f"[bins] width {w:.4f}; 1st-99th pct = [{q_lo:.2f}, {q_hi:.2f}], covers "
          f"{float(((df[score] >= q_lo) & (df[score] <= q_hi)).mean())*100:.2f}% "
          f"of neighbourhoods", flush=True)

    cands = {}
    for thr in THRESHOLDS:
        r = longest_run(counts, edges, thr)
        if r is None:
            continue
        lo, hi, cov = r
        sel = (df[score] >= lo) & (df[score] <= hi)
        nc = int(df.loc[sel, "TMA_ID"].nunique())
        # the weakest bin inside the candidate range, in cores
        inside = (centres >= lo) & (centres <= hi)
        cands[thr] = {"lo": lo, "hi": hi, "coverage": cov, "n_cores_inside": nc,
                      "min_cores_per_bin_inside": int(ncores[inside].min())}
        print(f"[cand] >= {thr:5d} nhoods/bin  -> [{lo:+.2f}, {hi:+.2f}]  "
              f"covers {cov*100:6.2f}%  |  {nc} cores, weakest bin has "
              f"{int(ncores[inside].min())} cores", flush=True)

    tab = pd.DataFrame({"bin_lo": edges[:-1], "bin_hi": edges[1:], "centre": centres,
                        "n_neighbourhoods": counts, "n_cores": ncores})
    tab["cum_frac"] = tab["n_neighbourhoods"].cumsum() / tab["n_neighbourhoods"].sum()
    tab.to_csv(outdir / f"crop_justification_bins{sfx}.csv", index=False)

    # ------------------------------------------------------------------ figure
    cols = plt.get_cmap("plasma")(np.linspace(0.05, 0.8, len(cands)))
    fig, axes = plt.subplots(3, 1, figsize=(12, 12), sharex=True)

    for ax, logy in zip(axes[:2], (True, False)):
        ax.bar(centres, counts, width=w * 0.94, color="#4c72b0", linewidth=0)
        if logy:
            ax.set_yscale("log"); ax.set_ylabel("neighbourhoods per bin (log)", fontsize=12)
            for c, (thr, d) in zip(cols, cands.items()):
                ax.axhline(thr, color=c, ls=":", lw=1.0)
                ax.text(edges[-1], thr, f"{thr} ", color=c, fontsize=8,
                        va="bottom", ha="right")
        else:
            ax.set_ylabel("neighbourhoods per bin", fontsize=12)
        for c, (thr, d) in zip(cols, cands.items()):
            ax.axvline(d["lo"], color=c, lw=1.2, alpha=0.9)
            ax.axvline(d["hi"], color=c, lw=1.2, alpha=0.9)
        ax.axvline(q_lo, color="black", ls="--", lw=1.4)
        ax.axvline(q_hi, color="black", ls="--", lw=1.4)
        sns.despine(ax=ax)

    axes[0].set_title(f"How many 2-hop neighbourhoods sit at each {score} score\n"
                      f"n = {len(df):,} neighbourhoods, {df['TMA_ID'].nunique()} cores, "
                      f"bin width {w:.3g}", fontsize=13)
    handles = [plt.Line2D([], [], color="black", ls="--", lw=1.4,
                          label=f"1st-99th pct  [{q_lo:.2f}, {q_hi:.2f}]")]
    handles += [plt.Line2D([], [], color=c, lw=1.2,
                           label=f"$\\geq${thr} nhoods/bin  [{d['lo']:.2f}, {d['hi']:.2f}]  "
                                 f"({d['coverage']*100:.2f}%)"
                                 + ("  (used)" if thr == USED_THRESHOLD else ""))
                for c, (thr, d) in zip(cols, cands.items())]
    axes[0].legend(handles=handles, fontsize=9, frameon=False, loc="lower center",
                   ncol=2, bbox_to_anchor=(0.5, -0.02))

    ax = axes[2]
    ax.bar(centres, ncores, width=w * 0.94, color="#dd8452", linewidth=0)
    ax.axhline(df["TMA_ID"].nunique(), color="grey", ls=":", lw=1)
    ax.text(edges[0], df["TMA_ID"].nunique(), f" all {df['TMA_ID'].nunique()} cores",
            fontsize=8, va="bottom", color="grey")
    for c, (thr, d) in zip(cols, cands.items()):
        ax.axvline(d["lo"], color=c, lw=1.2, alpha=0.9); ax.axvline(d["hi"], color=c, lw=1.2, alpha=0.9)
    ax.axvline(q_lo, color="black", ls="--", lw=1.4); ax.axvline(q_hi, color="black", ls="--", lw=1.4)
    ax.set_ylabel("distinct TMA cores per bin", fontsize=12)
    ax.set_xlabel(f"{score}", fontsize=13)
    ax.set_title("A spline fitted where only a handful of cores reach is describing those "
                 "cores, not the cohort", fontsize=11)
    sns.despine(ax=ax)

    fig.tight_layout()
    for e in ("pdf", "png"):
        fig.savefig(outdir / f"crop_justification_{args.stem}.{e}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    (outdir / f"13_crop_justification_meta{sfx}.json").write_text(json.dumps({
        "score": score, "scores_parquet": args.scores_parquet, "bin_width": w,
        "bins_per_1_99": args.bins_per_1_99, "n_neighbourhoods": int(len(df)),
        "n_cores": int(df["TMA_ID"].nunique()),
        "quantile_crop_1_99": [float(q_lo), float(q_hi)], "used_threshold": USED_THRESHOLD,
        "candidates": cands, "seed": SEED}, indent=2))
    print(f"[done] {outdir}/crop_justification_{args.stem}.pdf", flush=True)


if __name__ == "__main__":
    main()
