#!/usr/bin/env python
"""Stage E.2 -- concatenate the local-CCI chunks and draw Fig 3i (Vivien takeover).

Replaces Code/spatial_cci/2_postprocess_cci_res.py (== the `concat_cci_results.py` of
Vivien's email) and the `vs_factor2` panel of 3_plot_cci_score.py (== 132_plot_local_cci.py).

`cci_score_df.parquet` -- the input to every Fig 3i panel -- was never saved and is not on
disk anywhere, so it is rebuilt here from the 16 chunk pairs produced by Stage E.1.

The full table is 112,447,209 score rows over 254 cell-type pairs. Fig 3i needs only the
three cell-type pairs of interest, so chunks are filtered on load rather than concatenated
whole. A compact all-pairs summary (mean score and mean p per source|target|LR) is written
alongside for the heatmap panels.

Cell types are `level2_refined`, so the original's ('aPT','Immune'), ('Immune','Fib/MyoFib'),
('aPT','Fib/MyoFib') become ('aPT','Immune'), ('Immune','FIB'), ('aPT','FIB').
Factor 1 replaces Factor 2 (Leon, 2026-09-01). LOESS frac=0.4, as in the original.

Human GRCh38, HGNC symbols. Seed 42.
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import argparse, glob, json, time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from statsmodels.nonparametric.smoothers_lowess import lowess

HERE = Path(f"{KB_ANALYSIS}/vivien_takeover")
CCI = HERE / "outputs" / "local_cci_results"
FACTOR = "Factor1_ulm"
SEED = 42
LOESS_FRAC = 0.4
# Same Factor-1 crop as Fig 3h / ED 3f so the Factor-1 panels share an x-axis
# (scripts/13_crop_justification.py). Set to None for the 1st-99th percentile of this table.
CROP = (-7.25, 5.25)
CT_OF_INTEREST = [("aPT", "Immune"), ("Immune", "FIB"), ("aPT", "FIB")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top-n", type=int, default=12, help="LR pairs to label per panel")
    ap.add_argument("--reuse", action="store_true",
                    help="skip the 68-minute chunk concatenation and read the saved "
                         "cci_score_df.parquet (identical input, plotting only)")
    args = ap.parse_args()
    t0 = time.time()
    np.random.seed(SEED)
    sns.set_style("white")

    outdir = HERE / "figures" / "cci"; outdir.mkdir(parents=True, exist_ok=True)

    wanted = set()
    for a, b in CT_OF_INTEREST:
        wanted |= {(a, b), (b, a)}

    if args.reuse:
        df = pd.read_parquet(CCI / "cci_score_df.parquet")
        print(f"[reuse] cci_score_df.parquet {df.shape} "
              f"(built by the full run on 2026-09-03)", flush=True)
        return plot(df, args)

    score_files = sorted(glob.glob(str(CCI / "interaction_scores_*.parquet")))
    pval_files = sorted(glob.glob(str(CCI / "interaction_pvalues_*.parquet")))
    assert len(score_files) == len(pval_files) == 16, (len(score_files), len(pval_files))

    keep = ["interacting_pair", "gene_a", "gene_b", "root", "source", "target", "score"]
    parts, summ = [], []
    for i, (sf, pf) in enumerate(zip(score_files, pval_files), 1):
        s = pd.read_parquet(sf, columns=keep)
        p = pd.read_parquet(pf, columns=["interacting_pair", "root", "source", "target", "pvalue"])
        # compact all-pairs summary before filtering (mean over roots)
        summ.append(s.groupby(["source", "target", "interacting_pair"], observed=True)["score"]
                     .agg(["mean", "size"]))
        m = pd.Series(list(zip(s["source"], s["target"]))).isin(wanted).to_numpy()
        s = s[m]
        p = p[pd.Series(list(zip(p["source"], p["target"]))).isin(wanted).to_numpy()]
        parts.append(s.merge(p, on=["interacting_pair", "root", "source", "target"], how="left"))
        print(f"[load] {i}/16 {Path(sf).name}: kept {len(s):,} of the chunk "
              f"({time.time()-t0:.0f}s)", flush=True)

    df = pd.concat(parts, ignore_index=True)
    print(f"[concat] cci_score_df for the 3 cell-type pairs: {df.shape}", flush=True)

    full = (pd.concat(summ).groupby(level=[0, 1, 2])
              .apply(lambda g: pd.Series({"mean_score": np.average(g["mean"], weights=g["size"]),
                                          "n": g["size"].sum()})))
    full.to_parquet(CCI / "cci_summary_all_pairs.parquet")
    print(f"[save] all-pairs summary {full.shape} -> cci_summary_all_pairs.parquet", flush=True)

    # --- attach Factor 1 per root ------------------------------------------
    obs = pd.read_parquet(HERE / "outputs" / "xenium_factor_scores.parquet",
                          columns=["sample_cell_id", FACTOR, "TMA_ID", "condition"])
    obs = obs.set_index("sample_cell_id")
    df[FACTOR] = df["root"].map(obs[FACTOR])
    n_nan = int(df[FACTOR].isna().sum())
    print(f"[factor] {FACTOR} attached; {n_nan} rows ({n_nan/len(df)*100:.3f}%) without a value",
          flush=True)
    df = df.dropna(subset=[FACTOR])
    df.to_parquet(CCI / "cci_score_df.parquet")
    print(f"[save] cci_score_df.parquet {df.shape}", flush=True)
    return plot(df, args)


def plot(df, args, outdir=None, pairs=None, stem="fig3i_{a}_to_{b}"):
    t0 = time.time()
    outdir = HERE / "figures" / "cci" if outdir is None else Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    pairs = CT_OF_INTEREST if pairs is None else pairs

    # --- Fig 3i: LOESS of interaction score vs Factor 1, per LR pair --------
    q_lo, q_hi = df[FACTOR].quantile([0.01, 0.99])
    lo, hi = CROP if CROP is not None else (q_lo, q_hi)
    n_lo = int((df[FACTOR] < lo).sum()); n_hi = int((df[FACTOR] > hi).sum())
    inside = int(((df[FACTOR] >= lo) & (df[FACTOR] <= hi)).sum())
    print(f"[crop] [{lo:.2f}, {hi:.2f}] holds {inside:,} of {len(df):,} score rows "
          f"({inside/len(df)*100:.2f}%); {n_lo:,} below, {n_hi:,} above. "
          f"1st-99th pct of this table would be [{q_lo:.2f}, {q_hi:.2f}]", flush=True)
    x_seq = np.linspace(lo, hi, 100)
    meta = {"n_rows": int(len(df)), "crop": [float(lo), float(hi)],
            "crop_source": ("shared with Fig 3h / ED 3f -- scripts/13_crop_justification.py"
                            if CROP is not None else "1st-99th percentile of this table"),
            "quantile_crop_1_99": [float(q_lo), float(q_hi)],
            "rows_inside_crop": inside,
            "loess_frac": LOESS_FRAC, "cell_type_pairs": [list(p) for p in pairs]}

    for a, b in pairs:
        sub = df[(df["source"] == a) & (df["target"] == b)]
        if sub.empty:
            print(f"[fig3i] {a}->{b}: no rows, skipped", flush=True); continue
        fits = {}
        for lr, g in sub.groupby("interacting_pair", observed=True):
            g = g.dropna(subset=["score"])
            if len(g) < 50:
                continue
            sm = lowess(g["score"].to_numpy(float), g[FACTOR].to_numpy(float),
                        frac=LOESS_FRAC, xvals=x_seq, return_sorted=False)
            fits[lr] = sm
        fit = pd.DataFrame(fits, index=x_seq)
        rank = fit.max().sort_values(ascending=False)
        fit.to_csv(outdir / f"{stem.format(a=a, b=b)}_loess.csv")

        # LR pairs are named next to their own curve rather than in a legend
        # (Leon, 2026-09-07), matching the published panel. Labels sit at the right-hand
        # end of each curve and are pushed apart vertically so none overlap; the x-axis
        # is extended to hold them, and a short leader connects a moved label to its line.
        fig, ax = plt.subplots(figsize=(11, 7))
        pal = sns.color_palette("tab20", len(rank))
        for k, lr in enumerate(rank.index):
            ax.plot(x_seq, fit[lr], lw=1.8 if k < args.top_n else 0.8,
                    color=pal[k], alpha=1.0 if k < args.top_n else 0.45, zorder=3)

        span = float(np.nanmax(fit.to_numpy()) - np.nanmin(fit.to_numpy()))
        gap = 0.032 * span                       # minimum vertical spacing between labels
        lab = sorted(((float(fit[lr].iloc[-1]), lr, k) for k, lr in enumerate(rank.index)
                      if k < args.top_n), reverse=True)
        ys = [v[0] for v in lab]
        for i in range(1, len(ys)):              # push down, then relax upward
            ys[i] = min(ys[i], ys[i - 1] - gap)
        for i in range(len(ys) - 2, -1, -1):
            ys[i] = max(ys[i], ys[i + 1] + gap)
        x_lab = hi + 0.035 * (hi - lo)
        for (y0, lr, k), y in zip(lab, ys):
            # tab20 alternates light/dark; the light shades are unreadable as text, so the
            # label is drawn at 70 % of the line colour while the line keeps its own shade
            tcol = tuple(0.70 * c for c in pal[k][:3])
            ax.annotate(lr.replace("_", " - "), xy=(hi, y0), xytext=(x_lab, y),
                        color=tcol, fontsize=9, va="center", ha="left",
                        annotation_clip=False, zorder=5,
                        arrowprops=dict(arrowstyle="-", color=tcol, lw=0.6,
                                        shrinkA=0, shrinkB=0, alpha=0.7)
                        if abs(y - y0) > gap * 0.4 else None)
        ax.set_xlim(lo, hi)
        # make room so no label is pushed outside the axes
        y0a, y1a = ax.get_ylim()
        ax.set_ylim(min(y0a, min(ys) - gap), max(y1a, max(ys) + gap))
        ax.set_xlabel(f"{FACTOR} (cropped to [{lo:.1f}, {hi:.1f}])", fontsize=12)
        ax.set_ylabel("interaction score (loess smoothing)", fontsize=12)
        ax.set_title(f"CCI from {a} to {b}", fontsize=14)
        sns.despine(ax=ax)
        fig.subplots_adjust(right=0.72)
        for e in ("pdf", "png"):
            fig.savefig(outdir / f"{stem.format(a=a, b=b)}.{e}", dpi=300, bbox_inches="tight")
        plt.close(fig)
        meta[f"{a}->{b}"] = {"n_rows": int(len(sub)), "n_lr_pairs": int(fit.shape[1]),
                             "top_lr": rank.head(args.top_n).round(3).to_dict()}
        print(f"[fig3i] {a}->{b}: {fit.shape[1]} LR pairs, top = "
              f"{list(rank.head(5).index)}", flush=True)

    (outdir / "09_cci_meta.json").write_text(json.dumps(meta, indent=2, default=str))
    print(f"[done] {time.time()-t0:.0f}s -> {outdir}", flush=True)


if __name__ == "__main__":
    main()
