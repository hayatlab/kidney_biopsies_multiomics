#!/usr/bin/env python
"""Parameterised rebuild of the Fig-3j CCI panel: any LR pair, any two cores, boxplots.

Leon (2026-09-25):
  task 5  "use other cores and only show aPT and PT and maybe use boxplots instead of the
           ones there. Left i need Utrecht 2 F5 and right i need Utrecht 2 F2"
  task 6  "see if CD47-SIRPA comes up as a CCI pair and whether it has decent specificity
           for aPT"

Identical to 29_fig3j_app_cd74.py in data and statistics -- same local CellPhoneDB v5 run
(2-hop Voronoi neighbourhood per root cell, 200 permutations, seed 42), same per-core mean,
same PAIRED Wilcoxon over cores with a 2000-sample bootstrap CI on the median paired
difference, and the unpaired Mann-Whitney reported alongside. What changed:

  * senders are aPT and PT only (dPT dropped),
  * the two map cores are named on the command line instead of picked by the eGFR rule,
  * the distribution is a boxplot, not a violin,
  * the pair is a parameter, so CD47-SIRPA uses exactly the same code path as APP-CD74.

Core names are pathology identifiers and are NOT drawn; panels are labelled by disease and
eGFR, and the K-number goes to the meta file only.

Inputs : outputs/local_cci_results/interaction_scores_*.parquet, outputs/geodataframes/
         cell_boundaries.parquet, celltype_concordance/outputs/xenium_labeled.h5ad.
Outputs: figures/xenium/<stem>.{pdf,png,svg}, outputs/<stem>_mean_cci_per_core.csv,
         outputs/<stem>_meta.json.
Several pairs can be drawn in one run (`--pairs`), which matters because the per-root score
parquets are 384 MB and reading them once for five pairs instead of five times is the whole
cost of the script. `--pairs top5` ranks the 35 scored pairs by mean aPT->Immune score and
takes the top N; ligand and receptor are then read from `gene_a`/`gene_b` in the scored data,
not guessed by splitting the pair name.

Env spatialdata.  Usage:
  python scripts/34_cci_panel.py --pair APP_CD74  --ligand APP  --receptor CD74
  python scripts/34_cci_panel.py --pair CD47_SIRPA --ligand CD47 --receptor SIRPA
  python scripts/34_cci_panel.py --pairs top5          # the five strongest aPT->Immune pairs
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import argparse, glob, json, re, sys
from pathlib import Path

import anndata as ad
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import scipy.stats as ss  # noqa: E402
import seaborn as sns  # noqa: E402

HERE = Path(f"{KB_ANALYSIS}/vivien_takeover")
XENIUM = f"{KB_ANALYSIS}/celltype_concordance/outputs/xenium_labeled.h5ad"
TARGET = "Immune"
PANEL_CSV = Path(f"{KB_ANALYSIS}/apt_immune_lr/outputs"
                 "/apt_immune_pairs_panel.csv")
# the four ways of ranking the 35 scored pairs, all computed by apt_immune_lr/scripts/01
RELA_CSV = Path(f"{KB_ANALYSIS}/apt_immune_lr/outputs"
                "/rela_downstream_ligands.csv")
CRITERIA = {"mean_score": "mean aPT\u2192Immune score",
            "pct_cells_p_lt_0.05": "% of aPT roots with CellPhoneDB p < 0.05",
            "rho_factor1_core": "per-core Spearman rho vs Factor 1",
            "min_spec_xenium": "weaker of the two Xenium specificities"}
N_BOOT, SEED = 2000, 42
DCOL = {"CTRL": "#8c564b", "control": "#8c564b", "control_AKI": "#e377c2", "DN": "#9467bd",
        "FSGS": "#1f77b4", "IgAN": "#ff7f0e", "IgA": "#ff7f0e", "MN": "#2ca02c", "TIN": "#d62728"}
# Leon named the two cores for the maps (2026-09-25); matched on the leading tokens so the
# spacing of the TMA_ID ("Utrecht2 F5 ...") does not have to be typed exactly.
DEFAULT_CORES = ["Utrecht2 F5", "Utrecht2 F2"]


def resolve_core(pattern, cores):
    pat = re.sub(r"\s+", r"\\s*", pattern.strip())
    hit = [c for c in cores if re.match(pat, c, flags=re.I)]
    if len(hit) != 1:
        raise SystemExit(f"core pattern {pattern!r} matched {len(hit)} cores: {hit[:6]}")
    return hit[0]


def load_scores(pairs, senders, smoke):
    """Per-root scores for every requested pair at once -- the parquets are read only once."""
    files = sorted(glob.glob(str(HERE / "outputs/local_cci_results/interaction_scores_*.parquet")))
    if smoke:
        files = files[:2]
    cols = ["interacting_pair", "gene_a", "gene_b", "root", "source", "target", "score"]
    out = []
    for i, f in enumerate(files, 1):
        d = pd.read_parquet(f, columns=cols)
        out.append(d[d.interacting_pair.isin(pairs) & (d.target == TARGET)
                     & d.source.isin(senders)])
        print(f"  [scores] {i}/{len(files)} {Path(f).name}", flush=True)
    d = pd.concat(out, ignore_index=True)
    missing = sorted(set(pairs) - set(d.interacting_pair))
    if missing:
        raise SystemExit(f"{missing} -> {TARGET} have no scored rows; are they among the "
                         f"35 scored pairs?")
    d["core"] = d["root"].str.rsplit("_", n=1).str[0]
    print(f"[scores] {len(d):,} rows over {len(pairs)} pairs, {d.core.nunique()} cores, "
          f"senders {sorted(d.source.unique())}", flush=True)
    return d


def rank_pairs(n, senders, smoke):
    """The n strongest pairs by mean aPT->Immune score, from the precomputed summary."""
    f = HERE / "outputs/local_cci_results/cci_summary_all_pairs.parquet"
    summ = pd.read_parquet(f).reset_index()
    r = (summ[(summ.source == senders[0]) & (summ.target == TARGET)]
         .sort_values("mean_score", ascending=False).head(n))
    print(f"[rank] top {n} by mean {senders[0]}\u2192{TARGET} score: "
          + ", ".join(f"{a} ({b:.2f})" for a, b in zip(r.interacting_pair, r.mean_score)),
          flush=True)
    return list(r.interacting_pair)


def paired_stats(wide, a, b, rng):
    s = wide[[a, b]].dropna()
    diff = (s[a] - s[b]).to_numpy()
    boot = [float(np.median(rng.choice(diff, len(diff), replace=True))) for _ in range(N_BOOT)]
    return {"sender_a": a, "sender_b": b, "n_cores": int(len(s)),
            "median_paired_difference": float(np.median(diff)),
            "ci_lo": float(np.percentile(boot, 2.5)), "ci_hi": float(np.percentile(boot, 97.5)),
            "wilcoxon_p_paired": float(ss.wilcoxon(s[a], s[b]).pvalue),
            "mannwhitney_p_unpaired": float(ss.mannwhitneyu(wide[a].dropna(),
                                                            wide[b].dropna()).pvalue)}


def stars(p):
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"


def draw_pair(pair, lig, rec, senders, per, tests, wide, joined, picks, gdf, egfr, cond,
              cells, n_type, rng, fig_dir, out_dir, stem, tag):
    """One Fig-3j style panel: two core maps + the per-core distribution by sender."""
    n_maps = len(picks)
    fig = plt.figure(figsize=(5.2 * n_maps + 6.4, 5.4))
    gs = fig.add_gridspec(1, n_maps + 2, width_ratios=[1] * n_maps + [0.3, 1.2], wspace=0.35)
    vmax = {g: float(np.nanpercentile(
        joined.loc[joined.level2_refined.isin(senders if g == lig else [TARGET]), g], 99))
        for g in (lig, rec)}
    sender_of = {}
    for k, core in enumerate(picks):
        ax = fig.add_subplot(gs[0, k])
        gc = joined[joined.TMA_ID == core]
        gc.plot(ax=ax, facecolor="none", edgecolor="0.85", linewidth=0.12)
        send = max(senders, key=lambda x: n_type.loc[core, x])
        sender_of[core] = send
        gc[gc.level2_refined == send].plot(ax=ax, column=lig, cmap="Oranges", vmin=0,
                                           vmax=vmax[lig], linewidth=0.1, edgecolor="none")
        gc[gc.level2_refined == TARGET].plot(ax=ax, column=rec, cmap="Blues", vmin=0,
                                             vmax=vmax[rec], linewidth=0.1, edgecolor="none")
        ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(f"{cond.get(core)}, eGFR {float(egfr[core]):.0f}\n"
                     f"$\\bf{{{lig}}}$ in {send} to $\\bf{{{rec}}}$ in immune", fontsize=11)
    cax = fig.add_subplot(gs[0, n_maps]); cax.axis("off")
    for g, cmap, x in [(lig, "Oranges", 0.06), (rec, "Blues", 0.56)]:
        sub = cax.inset_axes([x, 0.30, 0.16, 0.42])
        fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0, vmax[g]), cmap=cmap), cax=sub)
        sub.set_title(g, fontsize=9, pad=4); sub.tick_params(labelsize=7)

    ax = fig.add_subplot(gs[0, n_maps + 1])
    order = [x for x in senders if x in set(per.source)]
    sns.boxplot(data=per, x="source", y="score", order=order, color="white", width=0.55,
                showfliers=False, linewidth=0.9, ax=ax,
                boxprops=dict(facecolor="white", edgecolor="0.25"),
                medianprops=dict(color="0.10", linewidth=1.5),
                whiskerprops=dict(color="0.25"), capprops=dict(color="0.25"))
    for dis, g in per.groupby("disease"):
        ax.scatter(pd.Categorical(g["source"], categories=order).codes
                   + rng.uniform(-0.16, 0.16, len(g)),
                   g["score"], s=9, color=DCOL.get(dis, "0.4"), label=dis, zorder=3, alpha=0.9)
    y = per.score.max()
    for i, t in enumerate(tests):
        j = order.index(t["sender_b"])
        yy = y * (1.06 + 0.1 * i)
        ax.plot([0, 0, j, j], [yy, yy * 1.02, yy * 1.02, yy], lw=0.8, c="black")
        ax.text(j / 2, yy * 1.03,
                f"{stars(t['wilcoxon_p_paired'])}  "
                f"({'aPT higher' if t['median_paired_difference'] > 0 else t['sender_b'] + ' higher'})",
                ha="center", va="bottom", fontsize=8.5)
    ax.set_ylim(0, y * 1.32)
    ax.set_xlabel("sender", fontsize=11); ax.set_ylabel("CCI score", fontsize=11)
    ax.set_title(f"{lig}\u2013{rec} with immune as target\n"
                 f"(average per core, {wide.dropna().shape[0]} cores)", fontsize=11)
    ax.legend(frameon=False, fontsize=8, ncol=2, loc="upper right", handletextpad=0.2,
              columnspacing=0.8)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    for e in ("pdf", "png", "svg"):
        fig.savefig(fig_dir / f"{stem}{tag}.{e}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    def kno(c):
        m = re.search(r"\bK\d+\b", c)
        return m.group(0) if m else None

    meta = {"script": "34_cci_panel.py", "pair": pair, "ligand": lig, "receptor": rec,
            "senders": senders, "target": TARGET, "n_cores": int(wide.dropna().shape[0]),
            "tests": tests,
            "plot": "boxplot (median, IQR, 1.5 IQR whiskers) + one point per core",
            "cores_shown": [{"side": ["left", "right"][i], "K_number": kno(c),
                             "condition": cond.get(c), "egfr": float(egfr[c]),
                             "n_cells": int(cells[c]), "ligand_drawn_in": sender_of[c],
                             "n_aPT": int(n_type.loc[c, "aPT"]), "n_PT": int(n_type.loc[c, "PT"]),
                             "n_Immune": int(n_type.loc[c, TARGET])}
                            for i, c in enumerate(picks)],
            "vmax_99th_percentile": vmax, "seed": SEED}
    (out_dir / f"{stem}_meta{tag}.json").write_text(json.dumps(meta, indent=2, default=str))
    print(f"[done] {fig_dir}/{stem}{tag}.pdf", flush=True)
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pair", default="APP_CD74")
    ap.add_argument("--pairs", default="",
                    help="comma-separated list, or 'topN' for the N strongest aPT->Immune pairs")
    ap.add_argument("--ligand", default="APP")
    ap.add_argument("--receptor", default="CD74")
    ap.add_argument("--senders", default="aPT,PT")
    ap.add_argument("--cores", default=",".join(DEFAULT_CORES),
                    help="two TMA_ID prefixes, left then right")
    ap.add_argument("--stem", default="")
    ap.add_argument("--outsub", default="",
                    help="write the panels into figures/xenium/<outsub>/ instead of "
                         "figures/xenium/")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    SENDERS = [x.strip() for x in args.senders.split(",")]
    tag = "_smoke" if args.smoke else ""
    rng = np.random.default_rng(SEED)
    fig_dir = HERE / "figures" / "xenium" / args.outsub if args.outsub \
        else HERE / "figures" / "xenium"
    out_dir = HERE / "outputs"
    fig_dir.mkdir(parents=True, exist_ok=True)

    import importlib.util
    spec = importlib.util.spec_from_file_location("xf", Path(__file__).resolve().parent / "05_xenium_figures.py")
    xf = importlib.util.module_from_spec(spec); spec.loader.exec_module(xf)

    gdf = gpd.read_parquet(HERE / f"outputs/geodataframes/cell_boundaries{tag}.parquet")
    egfr, _, _, _ = xf.corrected_core_egfr(gdf[["TMA_ID", "EGFR"]])
    egfr.index = egfr.index.astype(str)
    cond = gdf.groupby("TMA_ID", observed=True)["condition"].first().astype(str)
    all_cores = sorted(gdf["TMA_ID"].astype(str).unique())
    cells = gdf.groupby("TMA_ID", observed=True).size()
    n_type = gdf.groupby(["TMA_ID", "level2_refined"], observed=True).size().unstack(fill_value=0)

    if args.pairs.lower().startswith("top"):
        pairs = rank_pairs(int(args.pairs[3:]), SENDERS, args.smoke)
    elif args.pairs:
        pairs = [x.strip() for x in args.pairs.split(",")]
    else:
        pairs = [args.pair]

    d = load_scores(pairs, SENDERS, args.smoke)
    lr = (d.groupby("interacting_pair")[["gene_a", "gene_b"]].first()
          if args.pairs else None)

    picks = [resolve_core(x, all_cores) for x in args.cores.split(",")]
    for c in picks:
        print(f"[core] {cond.get(c)}  eGFR {float(egfr.get(c, np.nan)):.0f}  {cells[c]} cells  "
              f"aPT {n_type.loc[c, 'aPT']}  PT {n_type.loc[c, 'PT']}  "
              f"{TARGET} {n_type.loc[c, TARGET]}", flush=True)

    genes = sorted({g for pr in pairs for g in
                    ((lr.loc[pr, "gene_a"], lr.loc[pr, "gene_b"]) if lr is not None
                     else (args.ligand, args.receptor))})
    a = ad.read_h5ad(XENIUM)
    X = a.layers["counts"] if "counts" in a.layers else a.X
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
    lib = np.asarray(X.sum(axis=1)).ravel()
    key = a.obs["TMA_ID"].astype(str).str.strip() + "_" + a.obs["cell_id"].astype(str)
    expr = pd.DataFrame({g: np.log1p(np.asarray(X[:, a.var_names.get_loc(g)].todense()).ravel()
                                     / np.where(lib > 0, lib, 1) * 1e4) for g in genes},
                        index=key.to_numpy())
    del a, X
    joined = gdf.join(expr, how="left")
    frac = joined[genes[0]].notna().mean() * 100
    print(f"[expr] {len(genes)} genes joined onto {frac:.1f} % of cells", flush=True)
    assert frac > 95, "cell ids of the boundaries and the h5ad do not match"

    summary = []
    for pair in pairs:
        lig = lr.loc[pair, "gene_a"] if lr is not None else args.ligand
        rec = lr.loc[pair, "gene_b"] if lr is not None else args.receptor
        stem = (args.stem if args.stem and len(pairs) == 1
                else f"fig3j_{pair.lower()}_{'_'.join(x.lower() for x in SENDERS)}")
        dp = d[d.interacting_pair == pair]
        per = dp.groupby(["core", "source"]).score.mean().reset_index()
        per["disease"] = per["core"].map(cond)
        per["egfr"] = per["core"].map(egfr)
        per.to_csv(out_dir / f"{stem}_mean_cci_per_core{tag}.csv", index=False)
        wide = per.pivot(index="core", columns="source", values="score")
        tests = [paired_stats(wide, SENDERS[0], b, rng) for b in SENDERS[1:] if b in wide]
        print(f"\n=== {pair} ({lig} -> {rec}) ===\n{json.dumps(tests, indent=1)}", flush=True)
        m = draw_pair(pair, lig, rec, SENDERS, per, tests, wide, joined, picks, gdf, egfr,
                      cond, cells, n_type, rng, fig_dir, out_dir, stem, tag)
        summary.append({"pair": pair, "ligand": lig, "receptor": rec, "stem": stem,
                        "n_cores": m["n_cores"],
                        "mean_score_aPT": float(per.loc[per.source == SENDERS[0], "score"].mean()),
                        **{f"{k}_{t['sender_b']}": t[k] for t in tests
                           for k in ("median_paired_difference", "ci_lo", "ci_hi",
                                     "wilcoxon_p_paired")}})
    if len(pairs) > 1:
        sm = pd.DataFrame(summary)
        # why is each pair in the set? rank it under each of the four criteria
        if PANEL_CSV.exists():
            pan = pd.read_csv(PANEL_CSV).set_index("interacting_pair")
            for c in CRITERIA:
                r = pan[c].rank(ascending=False, method="min").astype(int)
                sm[f"rank_{c}"] = sm.pair.map(r)
                sm[c] = sm.pair.map(pan[c])
            sm["selected_by"] = [", ".join(
                k for k in CRITERIA if row[f"rank_{k}"] <= 5) or "\u2014"
                for _, row in sm.iterrows()]
        if RELA_CSV.exists():
            # how much evidence is there that the aPT-side molecule is downstream of RELA?
            # 0-5 lines from apt_immune_lr/scripts/04 (CRISPRa, kB peak gained, aPT RNA,
            # CollecTRI NF-kB target, Hallmark TNFA_SIGNALING_VIA_NFKB)
            ev = pd.read_csv(RELA_CSV).set_index("ligand")["n_evidence"]
            sm["rela_evidence_ligand"] = sm.ligand.map(ev).fillna(0).astype(int)
        sm.to_csv(out_dir / f"fig3j_pairs_summary{tag}.csv", index=False)
        print(f"\n{sm.to_string(index=False)}", flush=True)

        # ---- cross-pair summary: the paired aPT-PT difference, and why each pair is here
        sm2 = sm.sort_values("mean_score_aPT", ascending=False).reset_index(drop=True)
        dcol = f"median_paired_difference_{SENDERS[1]}"
        fig = plt.figure(figsize=(15.0, 0.52 * len(sm2) + 3.6))
        gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1], wspace=0.50)
        ax = fig.add_subplot(gs[0, 0])
        y = np.arange(len(sm2))
        sig = sm2[f"wilcoxon_p_paired_{SENDERS[1]}"] < 0.05
        col = np.where(sm2[dcol] > 0, "#c44e52", "#4c72b0")
        col = np.where(sig, col, "#bbbbbb")
        ax.errorbar(sm2[dcol], y,
                    xerr=[sm2[dcol] - sm2[f"ci_lo_{SENDERS[1]}"],
                          sm2[f"ci_hi_{SENDERS[1]}"] - sm2[dcol]],
                    fmt="none", ecolor="0.5", elinewidth=1.1, capsize=2.5, zorder=2)
        ax.scatter(sm2[dcol], y, s=60, c=col, zorder=3, edgecolors="white", linewidths=0.8)
        ax.axvline(0, color="0.35", lw=1.0, ls="--")
        ax.set_yticks(y)
        ax.set_yticklabels([f"{a}\u2013{b}" for a, b in zip(sm2.ligand, sm2.receptor)],
                           fontsize=9)
        ax.invert_yaxis()
        lo, hi = ax.get_xlim()                    # reserve a margin so the p labels never
        ax.set_xlim(lo, hi + 0.62 * (hi - lo))     # sit on top of a confidence interval
        for i, r in sm2.iterrows():
            ax.annotate(f"p = {r[f'wilcoxon_p_paired_{SENDERS[1]}']:.1e}  ({int(r.n_cores)} cores)",
                        (hi + 0.03 * (hi - lo), i), fontsize=7.5, color="0.35", va="center",
                        ha="left")
        ax.axvspan(hi, ax.get_xlim()[1], color="white", zorder=1)
        ax.set_xlabel(f"median per-core difference, {SENDERS[0]} \u2212 {SENDERS[1]} "
                      f"(95 % bootstrap CI)")
        ax.set_title(f"A  is {SENDERS[0]} really the stronger sender?\n"
                     f"grey = not significant at p < 0.05 (paired Wilcoxon over cores)",
                     fontsize=11.5, fontweight="bold", loc="left")

        ax = fig.add_subplot(gs[0, 1])
        if PANEL_CSV.exists():
            M = np.array([[sm2.loc[i, f"rank_{c}"] for c in CRITERIA] for i in range(len(sm2))],
                         dtype=float)
            im = ax.imshow(M, aspect="auto", cmap="YlGnBu_r", vmin=1, vmax=35)
            for i in range(M.shape[0]):
                for j in range(M.shape[1]):
                    ax.text(j, i, int(M[i, j]), ha="center", va="center", fontsize=8,
                            color="white" if M[i, j] <= 12 else "0.25",
                            fontweight="bold" if M[i, j] <= 5 else "normal")
            ax.set_xticks(range(len(CRITERIA)))
            ax.set_xticklabels([v.replace(" ", "\n") for v in CRITERIA.values()], fontsize=8)
            if "rela_evidence_ligand" in sm2:
                for i, v in enumerate(sm2.rela_evidence_ligand):
                    ax.annotate("\u25cf" * int(v) if v else "\u2013",
                                (len(CRITERIA) - 0.35, i), xycoords=("data", "data"),
                                fontsize=8, ha="left", va="center",
                                color="#c44e52" if v >= 3 else "0.45", annotation_clip=False)
                ax.annotate("RELA\nevidence", (len(CRITERIA) - 0.35, -0.75),
                            fontsize=8, ha="left", va="center", annotation_clip=False)
            ax.set_yticks(y)
            ax.set_yticklabels([f"{a}\u2013{b}" for a, b in zip(sm2.ligand, sm2.receptor)],
                               fontsize=9)
            fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02, label="rank of 35 scored pairs")
        ax.set_title("B  why each pair is in this set\n"
                     "bold = in the top 5 under that criterion", fontsize=11.5,
                     fontweight="bold", loc="left")
        for a in fig.axes[:1]:
            for spn in ("top", "right"):
                a.spines[spn].set_visible(False)
        fig.suptitle(
            f"aPT\u2192Immune ligand\u2013receptor pairs: {len(sm2)} of the 35 scored by the local "
            f"CellPhoneDB v5 run\nranked four different ways; one Fig-3j panel per pair in this "
            f"folder", fontsize=12.5, y=0.99)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        for e in ("pdf", "png"):
            fig.savefig(fig_dir / f"00_overview_pairs{tag}.{e}", dpi=200, bbox_inches="tight",
                        facecolor="white")
        plt.close(fig)
        print(f"[done] {len(pairs)} panels in {fig_dir}; summary in "
              f"{out_dir}/fig3j_pairs_summary{tag}.csv", flush=True)


if __name__ == "__main__":
    main()
