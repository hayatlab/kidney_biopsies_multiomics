#!/usr/bin/env python
"""Filter non-PEC lineage markers out of the cl9-vs-{0,2,4,11} MAST DE, then re-draw the two
volcanos with the 10 niche-20 marker genes force-labelled (+ ~10 additional by |log2FC|).

"non-PEC markers" = the union of the per-lineage CANONICAL marker panels (PANELS) used for the
level2 subclustering, for every lineage EXCEPT PEC, MINUS a protected set (the PEC panel itself +
the 10 niche-20 label genes). This is the user-selected "canonical non-PEC lineages" definition
(2026-07-30): a small, curated, reproducible list rather than a data-driven marker sweep.

Genes removed wherever they occur (both up- and down-regulated) since a lineage marker is non-PEC
regardless of DE direction. Protected genes are never removed even if they also appear in a non-PEC
panel (e.g. CD44/ITGB1/MMP7/VCAM1 are activated-PEC markers that also tag POD-injury / PT-injury).

Outputs (alongside the original DE table, suffixed *_pecfiltered):
  markers_MAST_9_vs_0_2_4_11_pecfiltered.csv   filtered DE table (feeds the GSEA re-run)
  removed_nonpec_genes.csv                     every excluded gene: source lineage(s) + DE stats
  volcano_9_vs_0_2_4_11_pecfiltered.{pdf,png}          standard volcano, niche-20 highlighted
  volcano_topFClabeled_9_vs_0_2_4_11_pecfiltered.{pdf,png}   |log2FC|-ranked labels + niche-20
  volcano_pecfiltered_labeledgenes.csv         the label set actually drawn

Labels use a deterministic gutter layout: every label is pushed into the left (down-regulated) or
right (up-regulated) margin with a leader line to its point, ordered by significance so leader lines
never cross -- this keeps labels clear of the dense point cloud (adjustText only repels labels from
each other, not from background points, so buried labels stay unreadable). Pure pandas/matplotlib.
Env: nichecompass.
"""
from __future__ import annotations
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import ast
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as _fm

for _p in ("/usr/share/fonts/liberation-sans/LiberationSans-Regular.ttf",
           "/usr/share/fonts/liberation-sans/LiberationSans-Bold.ttf"):
    try: _fm.fontManager.addfont(_p)
    except Exception: pass
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]
matplotlib.rcParams["pdf.fonttype"] = 42; matplotlib.rcParams["svg.fonttype"] = "none"

# ---- paths ----
DE_DIR = Path(f"{KB_ANALYSIS}/pec_resubcluster_20260717/"
              "de_9vs0241/de_results/9_vs_0_2_4_11/DE")
CSV = DE_DIR / "markers_MAST_9_vs_0_2_4_11.csv"
# canonical marker PANELS: parsed from the copy of subcluster_sweep.py shipped in this repository
PANELS_SRC = Path(__file__).resolve().parents[2] / "subcluster_sweep" / "subcluster_sweep.py"

# ---- config ----
P_CUT, FC_CUT = 0.05, 0.25
N_ADD = 10                       # curated volcano: additional labels beyond the niche genes (up/down by |lfc|)
N_TOP_LIST = [5, 10, 20]         # additional volcanos: top/bottom N significant genes by log2FC (no Xenium set)
NICHE20 = ["CD44", "ITGB1", "TMSB10", "PDLIM7", "MMP7",      # niche-20 top5 (up in niche 20)
           "CFH", "CDH1", "SLC4A11", "BMP7", "NDRG1"]        # niche-20 bottom5 (down in niche 20)
# shared injury / dediff / EMT / ECM markers that appear in non-PEC panels but are genuine
# crescent-PEC activation signal, not lineage contamination -> protect from removal (user, 2026-07-30)
EXTRA_PROTECT = ["VIM", "SPP1", "FN1", "HAVCR1", "PROM1", "SOX9", "ITGB6", "CD24",
                 "CST3", "POSTN", "CDKN1A", "CDKN2A"]
# lineage/stromal markers the colleague asked to KEEP eligible (un-exclude) for volcano + GSEA
# (2026-07-30). Not highlighted -- just no longer filtered out.
COLLEAGUE_KEEP = ["PDPN", "ACTA2", "TAGLN", "COL1A1", "COL1A2", "SERPINE2", "LRP2", "CUBN", "SLC5A2"]
# niche-20 genes NOT to label on the volcano (colleague: not significant in this contrast)
NICHE_LABEL_DROP = ["CD44", "MMP7", "TMSB10"]


def load_panels(src: Path) -> dict:
    """Extract the PANELS literal from subcluster_sweep.py without importing scanpy."""
    tree = ast.parse(src.read_text())
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == "PANELS" for t in node.targets):
            return ast.literal_eval(node.value)
    raise RuntimeError("PANELS not found in " + str(src))


def build_exclusion(panels: dict):
    """Return (exclusion: {gene: 'lineages'}, protected: set). PEC panel + niche-20 are protected."""
    pec_genes = {g for grp, gs in panels["PEC"].items() for g in gs}
    protected = pec_genes | set(NICHE20) | set(EXTRA_PROTECT) | set(COLLEAGUE_KEEP)
    src_of = {}
    for ct, groups in panels.items():
        if ct == "PEC":
            continue
        for grp, gs in groups.items():
            for g in gs:
                src_of.setdefault(g, set()).add(f"{ct}:{grp}")
    exclusion = {g: ";".join(sorted(l)) for g, l in src_of.items() if g not in protected}
    return exclusion, protected, pec_genes


def main():
    panels = load_panels(PANELS_SRC)
    exclusion, protected, pec_genes = build_exclusion(panels)

    m = pd.read_csv(CSV).dropna(subset=["p_val_adj", "avg_log2FC"]).copy()
    present = set(m["gene"])

    # transparency table: every excluded gene + whether it was in the DE table + its stats
    idx = m.set_index("gene")
    rows = []
    for g, lin in sorted(exclusion.items()):
        r = idx.loc[g] if g in present else None
        rows.append(dict(gene=g, lineages=lin, in_DE_table=g in present,
                         avg_log2FC=float(r["avg_log2FC"]) if r is not None else np.nan,
                         p_val_adj=float(r["p_val_adj"]) if r is not None else np.nan))
    removed = pd.DataFrame(rows)
    removed.to_csv(DE_DIR / "removed_nonpec_genes.csv", index=False)

    n_removed_present = int(removed["in_DE_table"].sum())
    m_filt = m[~m["gene"].isin(exclusion)].copy()
    m_filt.to_csv(DE_DIR / "markers_MAST_9_vs_0_2_4_11_pecfiltered.csv", index=False)

    print(f"PANELS lineages: {sorted(panels)}")
    print(f"exclusion set: {len(exclusion)} canonical non-PEC genes "
          f"({n_removed_present} present in DE table)")
    print(f"protected: {len(protected)} genes (PEC panel + niche-20 + injury/EMT extras + colleague-keep {COLLEAGUE_KEEP})")
    print(f"DE rows: {len(m)} -> {len(m_filt)} after filtering")
    # borderline injury/stress markers that got removed (flag for the user)
    borderline = removed[removed["gene"].isin(
        ["VIM", "HAVCR1", "PROM1", "ITGB6", "SPP1", "SOX9", "CD24", "CST3",
         "CDKN1A", "CDKN2A", "FN1", "POSTN"])]
    if len(borderline):
        print("\n[flag] shared injury/stress/ECM markers removed as non-PEC lineage markers:")
        print(borderline[["gene", "lineages", "in_DE_table", "avg_log2FC", "p_val_adj"]]
              .to_string(index=False))

    # ---- significance masks ----
    m_filt["neglog10padj"] = -np.log10(m_filt["p_val_adj"].clip(lower=np.nextafter(0, 1)))
    sig = m_filt["p_val_adj"] < P_CUT
    up = sig & (m_filt["avg_log2FC"] > FC_CUT)
    dn = sig & (m_filt["avg_log2FC"] < -FC_CUT)
    other_sig = sig & ~up & ~dn

    # ---- label set 1 (curated): significant Xenium genes + top N_ADD/2 up + N_ADD/2 down by |lfc| ----
    niche_present = [g for g in NICHE20 if g in set(m_filt["gene"]) and g not in NICHE_LABEL_DROP]
    half = N_ADD // 2
    add_up = m_filt[up & ~m_filt["gene"].isin(niche_present)].nlargest(half, "avg_log2FC")
    add_dn = m_filt[dn & ~m_filt["gene"].isin(niche_present)].nsmallest(half, "avg_log2FC")
    add_genes = list(add_up["gene"]) + list(add_dn["gene"])
    lab_curated = m_filt[m_filt["gene"].isin(niche_present + add_genes)].copy()
    lab_curated["is_xenium"] = lab_curated["gene"].isin(niche_present)
    lab_curated.sort_values("avg_log2FC", ascending=False)[
        ["gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "is_xenium"]
    ].to_csv(DE_DIR / "volcano_pecfiltered_labeledgenes.csv", index=False)
    print(f"\ncurated labels: {len(niche_present)} Xenium (sig) + {len(add_genes)} by |lfc|")
    print("  Xenium labelled:", niche_present, "| dropped (non-sig):", NICHE_LABEL_DROP)
    print("  additional:", add_genes)

    def annotate(ax, lab_df, band=1.7, shelf_dn=-2.7, shelf_up=2.9,
                 force_text=(0.4, 0.6), expand_text=(1.1, 1.25)):
        # Hybrid: labels whose point sits OUTSIDE the dense central band (|lfc| >= band) stay at their
        # point (adjustText nudges them off local neighbours). Labels INSIDE the band (buried in the
        # cloud) move to a moderate shelf just outside the band, with a leader line -- close, not flung.
        col_of = lambda lfc: "#d62728" if lfc > 0 else "#1f77b4"
        buried = lab_df[lab_df["avg_log2FC"].abs() < band]
        edge = lab_df[lab_df["avg_log2FC"].abs() >= band]

        texts = [ax.text(r["avg_log2FC"], r["neglog10padj"], r["gene"], fontsize=7.5,
                         color=col_of(r["avg_log2FC"]), zorder=6) for _, r in edge.iterrows()]
        pts = m_filt.loc[up | dn, ["avg_log2FC", "neglog10padj"]]
        try:
            from adjustText import adjust_text
            adjust_text(texts, x=pts["avg_log2FC"].to_numpy(), y=pts["neglog10padj"].to_numpy(), ax=ax,
                        arrowprops=dict(arrowstyle="-", color="0.55", lw=0.5),
                        expand_points=(1.6, 1.9), expand_text=expand_text,
                        force_text=force_text, force_points=(0.5, 0.8))
            status = "hybrid"
        except Exception as e:
            status = f"hybrid-noadjust ({type(e).__name__})"

        y0, y1 = ax.get_ylim(); gap = 0.045 * (y1 - y0)
        for sign, shelf, ha in [(1, shelf_up, "left"), (-1, shelf_dn, "right")]:
            d = buried[np.sign(buried["avg_log2FC"]) == sign].sort_values("neglog10padj", ascending=False)
            if len(d) == 0:
                continue
            ys = d["neglog10padj"].to_numpy().astype(float)     # declutter downward to keep min gap
            for i in range(1, len(ys)):
                ys[i] = min(ys[i], ys[i - 1] - gap)
            for (_, r), yy in zip(d.iterrows(), ys):
                ax.annotate(r["gene"], xy=(r["avg_log2FC"], r["neglog10padj"]), xytext=(shelf, yy),
                            ha=ha, va="center", fontsize=7.5, color=col_of(r["avg_log2FC"]), zorder=6,
                            annotation_clip=False,
                            arrowprops=dict(arrowstyle="-", color="0.55", lw=0.5, shrinkA=1, shrinkB=2))
        return status

    def draw(lab_df, title, outstem, figsize=(9.5, 8.0),
             force_text=(0.4, 0.6), expand_text=(1.1, 1.25)):
        fig, ax = plt.subplots(figsize=figsize); fig.patch.set_facecolor("white")
        ax.scatter(m_filt.loc[~sig, "avg_log2FC"], m_filt.loc[~sig, "neglog10padj"], s=6, lw=0,
                   color="#d9d9d9", label="NS")
        ax.scatter(m_filt.loc[other_sig, "avg_log2FC"], m_filt.loc[other_sig, "neglog10padj"], s=7, lw=0,
                   color="#bdbdbd", label=f"padj<{P_CUT}, |lfc|<{FC_CUT}")
        ax.scatter(m_filt.loc[up, "avg_log2FC"], m_filt.loc[up, "neglog10padj"], s=10, lw=0,
                   color="#d62728", label=f"UP in cl9 ({int(up.sum())})")
        ax.scatter(m_filt.loc[dn, "avg_log2FC"], m_filt.loc[dn, "neglog10padj"], s=10, lw=0,
                   color="#1f77b4", label=f"DOWN in cl9 ({int(dn.sum())})")
        ax.axhline(-np.log10(P_CUT), ls="--", c="0.6", lw=0.8)
        ax.axvline(FC_CUT, ls="--", c="0.6", lw=0.8); ax.axvline(-FC_CUT, ls="--", c="0.6", lw=0.8)
        status = annotate(ax, lab_df, force_text=force_text, expand_text=expand_text)
        ax.set_xlabel("avg log2 fold change  (cl9 vs cl0/2/4/11)")
        ax.set_ylabel("-log10 adjusted p (Bonferroni)")
        ax.set_title(title, fontsize=11)
        ax.legend(fontsize=8, frameon=False, loc="upper center", ncol=2)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        fig.tight_layout()
        fig.savefig(f"{outstem}.pdf", bbox_inches="tight", facecolor="white")
        fig.savefig(f"{outstem}.png", dpi=200, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print(f"[{status}] wrote {outstem}.pdf/.png")

    # Plot A: standard volcano (curated labels)
    draw(lab_curated,
         "cl9 vs cl0/2/4/11 (MAST) — non-PEC lineage markers removed\n"
         f"significant Xenium genes + top {N_ADD} by |log2FC| labelled",
         str(DE_DIR / "volcano_9_vs_0_2_4_11_pecfiltered"))
    # Plot B: |log2FC|-ranked title (curated labels)
    draw(lab_curated,
         "cl9 vs cl0/2/4/11 (MAST) — non-PEC markers removed; labels by |log2FC|\n"
         f"significant Xenium genes + top {half} up + {half} down among padj<{P_CUT}",
         str(DE_DIR / "volcano_topFClabeled_9_vs_0_2_4_11_pecfiltered"))
    # Plots C (additional): top-N/bottom-N by log2FC, ignoring the Xenium marker set -- one per N.
    # Canvas + label repulsion scale with label count (2N) so denser sets don't collide.
    top_specs = {5:  ((9.5, 8.0),   (0.5, 0.7), (1.15, 1.30)),
                 10: ((10.5, 9.0),  (0.7, 0.9), (1.30, 1.50)),
                 20: ((12.0, 11.0), (0.9, 1.2), (1.40, 1.70))}
    for n in N_TOP_LIST:
        lab_top = pd.concat([m_filt[up].nlargest(n, "avg_log2FC"),
                             m_filt[dn].nsmallest(n, "avg_log2FC")]).copy()
        lab_top.sort_values("avg_log2FC", ascending=False)[
            ["gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"]
        ].to_csv(DE_DIR / f"volcano_top{n}bottom{n}_labeledgenes.csv", index=False)
        fs, ft, et = top_specs.get(n, ((10.5, 9.0), (0.7, 0.9), (1.3, 1.5)))
        print(f"top{n}/bottom{n} by log2FC: labelling {len(lab_top)} genes")
        draw(lab_top,
             "cl9 vs cl0/2/4/11 (MAST) — non-PEC markers removed\n"
             f"top {n} up + {n} down by log2FC (padj<{P_CUT})",
             str(DE_DIR / f"volcano_top{n}bottom{n}_9_vs_0_2_4_11_pecfiltered"),
             figsize=fs, force_text=ft, expand_text=et)


if __name__ == "__main__":
    main()
