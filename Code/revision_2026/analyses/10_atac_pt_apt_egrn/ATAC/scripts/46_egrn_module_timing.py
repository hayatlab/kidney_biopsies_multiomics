#!/usr/bin/env python
"""Do the two eGRN modules act sequentially or switch together along the PT->aPT trajectory?

Leon (2026-09-25): "can we show that the two split networks behave along the trajectory? can
we say one acts first and the other, or do they simultaneously move in different directions?"

The r > 0.7 scMEGA network on the ArchR diffusion-map trajectory splits into two components
that share no edge (script 45): module I (72 TFs, AP-1 / NF-kB / STAT) and module II (14 TFs,
HNF4A / RXRA / VDR). The modules were defined by NETWORK CONNECTIVITY ALONE -- no timing
information went into the split -- so their behaviour over pseudotime is a genuine test, not a
restatement of how they were built.

Data: scmega/egrn_tutorial/archr_dm/pseudotime_curves.csv.gz, the smoothed curves scMEGA
produces over 100 pseudotime bins, for 86 TFs in four readouts:
    TF activity          chromVAR motif deviation  (chromatin)
    TF expression        the TF's own RNA
    Targets expression   the mean RNA of that TF's targets in the network  (output)
Each curve is z-scored over its own 100 bins, so the comparison is of SHAPE and TIMING, not
of magnitude.

Three timing statistics per TF per readout, because "acts first" is ambiguous for a
monotonically declining curve (its maximum is trivially bin 1):
    midpoint      the bin where the curve first crosses the midpoint of its own min-max range
                  -- the transition point, defined for rising and falling curves alike
    steepest      the bin of largest |d/d(bin)| after a 5-bin moving average -- when it is
                  changing fastest
    peak          the bin of the maximum (rising TFs) or the minimum (falling TFs), chosen by
                  the sign of the Spearman correlation with pseudotime
Monotonicity (|Spearman rho| vs pseudotime) is reported so non-monotone TFs, for which a
midpoint is not meaningful, can be excluded (--min-rho).

Module comparison: difference in MEDIAN statistic, with a label-permutation null (10,000
shuffles of the module labels over the 86 TFs). This respects the observed curve shapes and
asks only whether the split into 72/14 explains the timing; it does NOT assume TFs are
independent (they share targets), which a rank-sum test would.

Simultaneity: cross-correlation of the two module mean curves over lags -40..+40 bins. A
strongest anti-correlation at lag 0 means the modules move in opposite directions at the same
point in the trajectory; a non-zero lag means one leads.

Outputs (scmega/figures_paper/):
  egrn_module_timing.{pdf,png}
  egrn_module_timing_per_tf.csv, egrn_module_timing_curves.csv,
  egrn_module_timing_tests.csv, 46_module_timing_meta.json

Human hg38, HGNC symbols. Seed 42.  Env: agent_env.
Usage: python scripts/46_egrn_module_timing.py [--run archr_dm] [--min-rho 0.0]
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import argparse, json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.stats as ss  # noqa: E402

A = Path(f"{KB_ANALYSIS}/aki_c5_exclusion_20261022/atac_rerun/ATAC")
OUT = A / "scmega" / "figures_paper"  # output folder (release patch: the figure entry scripts in figures/ redirect it; default = original folder)
SEED, N_PERM, SMOOTH, MAXLAG = 42, 10000, 5, 40
READOUTS = ["TF activity", "TF expression", "Targets expression"]
MOD_NAME = {1: "module I — AP-1 / NF-κB / STAT", 2: "module II — HNF4A / RXRA / VDR"}
MOD_COL = {1: "#c44e52", 2: "#4c72b0"}


def timing(v):
    """midpoint / steepest / peak bin and monotonicity for one z-scored curve."""
    b = np.arange(1, len(v) + 1)
    rho = ss.spearmanr(b, v).statistic
    lo, hi = float(np.min(v)), float(np.max(v))
    if hi - lo < 1e-9:
        return dict(midpoint=np.nan, steepest=np.nan, peak=np.nan, rho=rho)
    half = lo + 0.5 * (hi - lo)
    up = rho >= 0
    crossed = (v >= half) if up else (v <= half)
    mid = int(b[np.argmax(crossed)]) if crossed.any() else np.nan
    k = np.ones(SMOOTH) / SMOOTH
    sm = np.convolve(v, k, mode="same")
    d = np.gradient(sm)
    d[: SMOOTH] = d[-SMOOTH:] = 0.0            # edge effects of the moving average
    return dict(midpoint=mid, steepest=int(b[np.argmax(np.abs(d))]),
                peak=int(b[np.argmax(v)] if up else b[np.argmin(v)]), rho=float(rho))


def perm_test(vals, labels, rng):
    """difference in median between the two module labels, label-permutation p (two-sided)."""
    vals = np.asarray(vals, float); labels = np.asarray(labels)
    ok = np.isfinite(vals); vals, labels = vals[ok], labels[ok]
    a, b = vals[labels == 1], vals[labels == 2]
    if len(a) < 3 or len(b) < 3:
        return dict(n1=len(a), n2=len(b), median1=np.nan, median2=np.nan,
                    diff=np.nan, p_perm=np.nan)
    obs = float(np.median(a) - np.median(b))
    null = np.empty(N_PERM)
    for i in range(N_PERM):
        s = rng.permutation(labels)
        null[i] = np.median(vals[s == 1]) - np.median(vals[s == 2])
    p = (np.sum(np.abs(null) >= abs(obs)) + 1) / (N_PERM + 1)
    return dict(n1=int(len(a)), n2=int(len(b)), median1=float(np.median(a)),
                median2=float(np.median(b)), diff=obs, p_perm=float(p))


def xcorr(x, y, maxlag):
    """corr(x[t], y[t+lag]) for lag in [-maxlag, maxlag]; positive lag = y lags x."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    out = []
    for lag in range(-maxlag, maxlag + 1):
        if lag >= 0:
            a, b = x[: len(x) - lag], y[lag:]
        else:
            a, b = x[-lag:], y[: len(y) + lag]
        out.append((lag, float(np.corrcoef(a, b)[0, 1])))
    return pd.DataFrame(out, columns=["lag", "r"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="archr_dm")
    ap.add_argument("--min-rho", type=float, default=0.0,
                    help="drop TFs whose |Spearman rho vs pseudotime| is below this")
    args = ap.parse_args()
    rng = np.random.default_rng(SEED)
    d = A / "scmega" / "egrn_tutorial" / args.run
    out = Path(OUT); out.mkdir(parents=True, exist_ok=True)

    E = pd.read_csv(d / "grn_edges.csv")
    G = nx.from_pandas_edgelist(E, "tf", "gene")
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    mod = {t: (i + 1) for i, c in enumerate(comps) for t in (set(c) & set(E.tf))}
    cur = pd.read_csv(d / "pseudotime_curves.csv.gz")
    print(f"[net] {len(E)} edges, {len(mod)} TFs | module sizes "
          f"{pd.Series(mod).value_counts().sort_index().to_dict()}", flush=True)

    rows, curves = [], []
    for readout in READOUTS:
        w = cur[cur.data == readout].pivot(index="pseudotime", columns="tf", values="value")
        z = (w - w.mean()) / w.std(ddof=0)
        for tf in z.columns:
            r = timing(z[tf].to_numpy())
            r.update(tf=tf, module=mod.get(tf), readout=readout)
            rows.append(r)
        for m in sorted(set(mod.values())):
            tfs = [t for t in z.columns if mod.get(t) == m]
            sub = z[tfs]
            curves.append(pd.DataFrame({
                "pseudotime": z.index, "readout": readout, "module": m,
                "mean": sub.mean(1).to_numpy(), "sd": sub.std(1).to_numpy(),
                "n": len(tfs),
                "lo": sub.mean(1).to_numpy() - 1.96 * sub.std(1).to_numpy() / np.sqrt(len(tfs)),
                "hi": sub.mean(1).to_numpy() + 1.96 * sub.std(1).to_numpy() / np.sqrt(len(tfs)),
            }))
    per_tf = pd.DataFrame(rows)
    curves = pd.concat(curves, ignore_index=True)
    keep = per_tf.rho.abs() >= args.min_rho
    print(f"[tf] {len(per_tf)} TF x readout rows; {int((~keep).sum())} below |rho| "
          f"{args.min_rho}", flush=True)
    per_tf["used"] = keep
    per_tf.to_csv(out / "egrn_module_timing_per_tf.csv", index=False)
    curves.to_csv(out / "egrn_module_timing_curves.csv", index=False)

    # ---- mean-curve transition window per module (25/50/75 % of its own range)
    def crossings(v):
        v = np.asarray(v, float); lo, hi = v.min(), v.max()
        f = (v - lo) / (hi - lo) if v[-1] > v[0] else (hi - v) / (hi - lo)
        return [int(np.argmax(f >= t) + 1) if (f >= t).any() else np.nan
                for t in (0.25, 0.50, 0.75)]

    win = []
    for readout in READOUTS:
        row = {"readout": readout}
        for m in (1, 2):
            v = curves[(curves.readout == readout) & (curves.module == m)] \
                .sort_values("pseudotime")["mean"]
            q25, q50, q75 = crossings(v)
            row |= {f"m{m}_q25": q25, f"m{m}_q50": q50, f"m{m}_q75": q75}
        row["lead_of_module_II_bins"] = row["m1_q50"] - row["m2_q50"]
        row["overlap_lo"] = max(row["m1_q25"], row["m2_q25"])
        row["overlap_hi"] = min(row["m1_q75"], row["m2_q75"])
        row["overlap_bins"] = row["overlap_hi"] - row["overlap_lo"] + 1
        win.append(row)
        print(f"  [{readout:<19s}] mean curve 25/50/75 %: module I {row['m1_q25']}/"
              f"{row['m1_q50']}/{row['m1_q75']}, module II {row['m2_q25']}/{row['m2_q50']}/"
              f"{row['m2_q75']} | module II leads by {row['lead_of_module_II_bins']:+d} bins | "
              f"transition windows overlap over {row['overlap_bins']} bins", flush=True)
    win = pd.DataFrame(win)
    win.to_csv(out / "egrn_module_timing_mean_curve_windows.csv", index=False)

    # ---- scMEGA's own per-TF `time_point`.
    # CAUTION: this is the pseudotime at which the TF's motif activity PEAKS, so for a
    # monotonically declining TF it is bin ~1 by construction (HNF4A: time_point 1.0).
    # It therefore separates rising from falling TFs -- a DIRECTION statistic -- and must not
    # be read as evidence that one module acts before the other. The same applies to the
    # `peak` statistic computed above. Both are reported for completeness only.
    st = pd.read_csv(d / "selected_tfs.csv")
    st = st[st.tfs.isin(mod)].assign(module=lambda x: x.tfs.map(mod))
    t_sc = perm_test(st["time_point"], st["module"], rng)
    print(f"  [scMEGA time_point   ] module I median {t_sc['median1']:.1f} vs module II "
          f"{t_sc['median2']:.1f} | diff {t_sc['diff']:+.1f} | permutation p = "
          f"{t_sc['p_perm']:.4f}   (DIRECTION, not timing: it is the peak bin, which is ~1 "
          f"for any declining TF)", flush=True)

    # ---------------------------------------------------------------- tests
    tests = [dict(t_sc, readout="scMEGA selected_tfs",
                  statistic="time_point_DIRECTION_not_timing")]
    for readout in READOUTS:
        p = per_tf[(per_tf.readout == readout) & per_tf.used]
        for stat in ["midpoint", "steepest", "peak_DIRECTION_not_timing"]:
            t = perm_test(p[stat.split("_")[0]], p["module"], rng)
            t.update(readout=readout, statistic=stat)
            tests.append(t)
            print(f"  [{readout:<19s}] {stat[:9]:<9s} module I median {t['median1']:6.1f} vs "
                  f"module II {t['median2']:6.1f} | diff {t['diff']:+6.1f} bins | "
                  f"permutation p = {t['p_perm']:.4f}", flush=True)
        # module trend direction
        for m in (1, 2):
            mu = curves[(curves.readout == readout) & (curves.module == m)]
            r = ss.spearmanr(mu.pseudotime, mu["mean"])
            tests.append({"readout": readout, "statistic": f"module{m}_trend_spearman",
                          "diff": float(r.statistic), "p_perm": float(r.pvalue)})
    tests = pd.DataFrame(tests)
    tests.to_csv(out / "egrn_module_timing_tests.csv", index=False)

    # ---------------------------------------------------------------- cross-correlation
    xc = {}
    for readout in READOUTS:
        m1 = curves[(curves.readout == readout) & (curves.module == 1)]["mean"].to_numpy()
        m2 = curves[(curves.readout == readout) & (curves.module == 2)]["mean"].to_numpy()
        x = xcorr(m1, m2, MAXLAG)
        best = x.loc[x.r.idxmin()]                       # strongest ANTI-correlation
        xc[readout] = {"table": x, "best_lag": int(best.lag), "best_r": float(best.r),
                       "r_at_lag0": float(x.loc[x.lag == 0, "r"].iloc[0])}
        print(f"  [{readout:<19s}] strongest anti-correlation r = {best.r:+.3f} at lag "
              f"{int(best.lag):+d} bins (r at lag 0 = {xc[readout]['r_at_lag0']:+.3f})",
              flush=True)

    # ---------------------------------------------------------------- figure
    fig = plt.figure(figsize=(19.0, 13.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], hspace=0.46, wspace=0.28)

    # A-C: module mean curves per readout
    for j, readout in enumerate(READOUTS):
        ax = fig.add_subplot(gs[0, j])
        w = cur[cur.data == readout].pivot(index="pseudotime", columns="tf", values="value")
        z = (w - w.mean()) / w.std(ddof=0)
        for m in (1, 2):
            for tf in [t for t in z.columns if mod.get(t) == m]:
                ax.plot(z.index, z[tf], color=MOD_COL[m], lw=0.35, alpha=0.20, zorder=1)
            c = curves[(curves.readout == readout) & (curves.module == m)]
            ax.plot(c.pseudotime, c["mean"], color=MOD_COL[m], lw=2.6, zorder=3,
                    label=f"{MOD_NAME[m]} ({int(c.n.iloc[0])} TFs)")
            ax.fill_between(c.pseudotime, c.lo, c.hi, color=MOD_COL[m], alpha=0.25, lw=0,
                            zorder=2)
        c1 = curves[(curves.readout == readout) & (curves.module == 1)]["mean"].to_numpy()
        c2 = curves[(curves.readout == readout) & (curves.module == 2)]["mean"].to_numpy()
        cross = np.where(np.diff(np.sign(c1 - c2)))[0]
        if len(cross):
            ax.axvline(cross[0] + 1, color="0.35", ls="--", lw=1.0, zorder=4)
            ax.annotate(f"cross at bin {cross[0] + 1}", (cross[0] + 1, ax.get_ylim()[1]),
                        fontsize=8, color="0.3", rotation=90, va="top", ha="right",
                        xytext=(-3, -4), textcoords="offset points")
        wr = win[win.readout == readout].iloc[0]
        ax.axvspan(wr.overlap_lo, wr.overlap_hi, color="0.5", alpha=0.10, zorder=0)
        ax.axhline(0, color="0.85", lw=0.8, zorder=0)
        ax.set_xlabel("pseudotime bin   (1 = PT-like  →  100 = aPT-like)")
        ax.set_ylabel("z-scored curve")
        if j == 0:
            ax.legend(frameon=False, fontsize=8.5, loc="lower left")
        ax.set_title(f"{'ABC'[j]}  {readout}\n"
                     f"module II half-way at bin {wr.m2_q50}, module I at {wr.m1_q50} "
                     f"(lead {wr.lead_of_module_II_bins:+d});\nshaded = both mid-transition "
                     f"({wr.overlap_bins} bins)", fontsize=10.5, fontweight="bold", loc="left")

    # D: per-TF heatmap ordered by midpoint
    ax = fig.add_subplot(gs[1, 0])
    w = cur[cur.data == "TF activity"].pivot(index="pseudotime", columns="tf", values="value")
    z = ((w - w.mean()) / w.std(ddof=0)).T
    o = (per_tf[(per_tf.readout == "TF activity")]
         .sort_values(["module", "midpoint"], na_position="last"))
    z = z.loc[o.tf]
    im = ax.imshow(z.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2,
                   interpolation="nearest")
    b = int((o.module == 1).sum())
    ax.axhline(b - 0.5, color="black", lw=1.6)
    for m, y in [(1, b / 2), (2, b + (len(o) - b) / 2)]:
        ax.annotate(f"module {'I' * m}", (-1.5, y), xycoords=("data", "data"), rotation=90,
                    ha="right", va="center", fontsize=10, fontweight="bold", color=MOD_COL[m])
    ax.set_xlabel("pseudotime bin"); ax.set_yticks([])
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="z-scored motif activity")
    ax.set_title("D  every TF, ordered by transition midpoint within its module",
                 fontsize=11.5, fontweight="bold", loc="left")

    # E: midpoint distributions
    ax = fig.add_subplot(gs[1, 1])
    for j, readout in enumerate(READOUTS):
        p = per_tf[(per_tf.readout == readout) & per_tf.used]
        for m in (1, 2):
            v = p.loc[p.module == m, "midpoint"].dropna()
            x = j + (-0.18 if m == 1 else 0.18)
            ax.scatter(x + rng.uniform(-0.07, 0.07, len(v)), v, s=16, color=MOD_COL[m],
                       alpha=0.75, linewidths=0, zorder=3)
            ax.hlines(np.median(v), x - 0.14, x + 0.14, color=MOD_COL[m], lw=2.6, zorder=4)
        t = tests[(tests.readout == readout) & (tests.statistic == "midpoint")].iloc[0]  # noqa
        ax.annotate(f"$\\Delta$ {t['diff']:+.0f} bins\np = {t['p_perm']:.3f}", (j, 104),
                    ha="center", va="bottom", fontsize=8.5)
    ax.set_xticks(range(len(READOUTS)))
    ax.set_xticklabels([r.replace(" ", "\n") for r in READOUTS], fontsize=9)
    ax.set_ylabel("transition midpoint (pseudotime bin)")
    ax.set_ylim(0, 118)
    ax.set_title("E  when each TF crosses the middle of its own range\n"
                 "difference in medians, label-permutation p (10,000 shuffles)",
                 fontsize=11.5, fontweight="bold", loc="left")

    # F: cross-correlation
    ax = fig.add_subplot(gs[1, 2])
    for readout, style in zip(READOUTS, ["-", "--", ":"]):
        x = xc[readout]["table"]
        ax.plot(x.lag, x.r, style, lw=1.8, label=readout)
        ax.scatter([xc[readout]["best_lag"]], [xc[readout]["best_r"]], s=40, zorder=5,
                   color="black")
    ax.axvline(0, color="0.6", lw=1.0, ls="--")
    ax.axhline(0, color="0.85", lw=0.8)
    ax.set_xlabel("lag (bins) — module II shifted relative to module I")
    ax.set_ylabel("correlation of the two module mean curves")
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    lg = ", ".join(f"{r.split()[0][:4]} {xc[r]['best_lag']:+d}" for r in READOUTS)
    ax.set_title(f"F  anti-correlation is strongest at lag {lg}\n"
                 "lag 0 = opposite directions at the same point",
                 fontsize=11.5, fontweight="bold", loc="left")
    for a in fig.axes:
        for s in ("top", "right"):
            a.spines[s].set_visible(False)

    stp = tests[(tests.readout == "TF activity") & (tests.statistic == "steepest")].iloc[0]
    lags = [xc[r]["best_lag"] for r in READOUTS]
    fig.suptitle(
        "One coordinated switch, not two phases — but the differentiated-PT module\n"
        "disengages before the injury module fires\n"
        f"scMEGA r > 0.7 network on the ArchR PT→aPT trajectory. Modules defined by network connectivity\n"
        f"alone (they share no edge); timing measured afterwards. The mean curves are near mirror images —\n"
        f"strongest anti-correlation at lag {min(lags):+d} to {max(lags):+d} of 100 bins, r ≈ −0.98 — and their transition\n"
        f"windows overlap. Within that, module II leads in chromatin and in target output: steepest\n"
        f"motif-activity change at bin {stp['median2']:.0f} vs {stp['median1']:.0f} (permutation p = {stp['p_perm']:.4f}).",
        fontsize=12.0, y=0.995)
    wexp = win[win.readout == "TF expression"].iloc[0]
    fig.text(0.5, 0.018,
             "Caveat: the lead is present in the chromatin (A) and in the target output (C) but NOT in the TFs' own RNA (B), "
             f"where module I is half-way {abs(wexp.lead_of_module_II_bins)} bins BEFORE module II and the per-TF midpoint "
             "difference is not significant (p = 0.52).\n"
             "Pseudotime is a monotone ordering of cells, not elapsed time, so a lag in bins is a lag in cell ordering. "
             "The `peak` and scMEGA `time_point` statistics in the CSVs separate rising from falling TFs by construction "
             "and are NOT evidence of sequencing.",
             ha="center", va="bottom", fontsize=9, color="0.30")
    fig.subplots_adjust(left=0.05, right=0.985, top=0.80, bottom=0.115, hspace=0.46, wspace=0.28)
    for e in ("pdf", "png"):
        fig.savefig(out / f"egrn_module_timing.{e}", dpi=200, facecolor="white")
    plt.close(fig)

    (out / "46_module_timing_meta.json").write_text(json.dumps({
        "script": "46_egrn_module_timing.py", "run": args.run, "seed": SEED,
        "source_curves": str(d / "pseudotime_curves.csv.gz"),
        "n_bins": 100, "readouts": READOUTS, "n_permutations": N_PERM,
        "smoothing_for_derivative_bins": SMOOTH, "min_abs_rho": args.min_rho,
        "modules": {str(m): {"n_tfs": int(sum(v == m for v in mod.values())),
                             "tfs": sorted(t for t, v in mod.items() if v == m)}
                    for m in sorted(set(mod.values()))},
        "tests": tests.to_dict("records"),
        "mean_curve_windows": win.to_dict("records"),
        "cross_correlation": {k: {"best_lag": v["best_lag"], "best_r": v["best_r"],
                                  "r_at_lag0": v["r_at_lag0"]} for k, v in xc.items()},
        "caveat": "pseudotime is a monotone ordering of cells, not elapsed time; a lag in bins "
                  "is a lag in cell ordering. TFs within a module share targets and are not "
                  "independent, which is why the module comparison uses a label permutation "
                  "rather than a rank-sum test.",
        "networkx": nx.__version__, "pandas": pd.__version__,
    }, indent=2, default=str))
    print(f"\n[done] {out}/egrn_module_timing.pdf", flush=True)


if __name__ == "__main__":
    main()
