#!/usr/bin/env python
"""
Excel export of the 3-way Venn regions (restricted universe, top200 cut) for sharing.

Reads the outputs of venn_niche20_mofa_factor1.py and writes one workbook with every gene of the
3-way Venn labelled by which region it falls in, annotated with its niche-20 DE stats and its actual
PT / PEC Factor 1 loadings.

env: spatial_env (openpyxl 3.1.5)
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)

import os
import pandas as pd

BASE = f"{KB_ANALYSIS}"
VENN = os.path.join(BASE, "niche20_pec/venn_mofa_factor1")
SRC = os.path.join(VENN, "main_restricted_universe/top200")
DE_CSV = os.path.join(
    BASE, "niche20_pec/DE/niche20_vs_others_pecOnly/markers_MAST_niche20_vs_others_pecOnly.csv"
)
LOADINGS_CSV = os.path.join(BASE, "loadings~hvg~5.csv")
OUT = os.path.join(VENN, "niche20_vs_MOFA_Factor1_3way_top200.xlsx")

REGION_ORDER = ["all three", "PT + niche20", "PEC + niche20", "niche20 only",
                "PT + PEC", "PT only", "PEC only"]


def region_label(sets_str, direction):
    parts = set(sets_str.split("|"))
    has_pt = any(p.startswith("PT_F1") for p in parts)
    has_pec = any(p.startswith("PEC_F1") for p in parts)
    has_n20 = any(p.startswith("niche20") for p in parts)
    if has_pt and has_pec and has_n20:
        return "all three"
    if has_pt and has_n20:
        return "PT + niche20"
    if has_pec and has_n20:
        return "PEC + niche20"
    if has_pt and has_pec:
        return "PT + PEC"
    if has_n20:
        return "niche20 only"
    if has_pt:
        return "PT only"
    return "PEC only"


def build():
    memb = pd.read_csv(os.path.join(SRC, "gene_set_membership.csv"))
    de = pd.read_csv(DE_CSV).set_index("gene")
    L = pd.read_csv(LOADINGS_CSV)
    f1 = L[L["Factor"] == "Factor1"]
    pt = f1[f1["view"] == "PT"].set_index("variable")["Loading"]
    pec = f1[f1["view"] == "PEC"].set_index("variable")["Loading"]
    pt_rank = pt.abs().rank(ascending=False).astype(int)
    pec_rank = pec.abs().rank(ascending=False).astype(int)

    tabs = {}
    for direction in ("UP", "DOWN"):
        s = memb[memb["venn"] == f"3way_{direction}"].copy()
        s["region"] = [region_label(x, direction) for x in s["sets"]]
        s["avg_log2FC"] = s["gene"].map(de["avg_log2FC"])
        s["p_val_adj"] = s["gene"].map(de["p_val_adj"])
        s["pct_niche20"] = s["gene"].map(de["pct.1"])
        s["pct_other_PEC"] = s["gene"].map(de["pct.2"])
        s["PT_F1_loading"] = s["gene"].map(pt)
        s["PT_F1_absrank"] = s["gene"].map(pt_rank)
        s["PEC_F1_loading"] = s["gene"].map(pec)
        s["PEC_F1_absrank"] = s["gene"].map(pec_rank)
        s["region"] = pd.Categorical(s["region"], REGION_ORDER, ordered=True)
        s = s.sort_values(["region", "gene"])
        tabs[direction] = s[["gene", "region", "avg_log2FC", "p_val_adj",
                             "pct_niche20", "pct_other_PEC",
                             "PT_F1_loading", "PT_F1_absrank",
                             "PEC_F1_loading", "PEC_F1_absrank"]]

    counts = pd.concat([
        tabs[d]["region"].value_counts().rename(d).reindex(REGION_ORDER).fillna(0).astype(int)
        for d in ("UP", "DOWN")], axis=1).rename_axis("Venn region").reset_index()

    notes = pd.DataFrame({"": [
        "niche 20 (crescent PEC) DE  vs  MOFA Factor 1 loadings — 3-way Venn regions",
        "Leon Frommer · Hayat Group · generated 2026-08-04",
        "",
        "SCOPE OF THIS WORKBOOK",
        "  Restricted universe (Xenium panel n MOFA view; M = 207 genes present in the 468-gene",
        "  Xenium panel AND in both the PT and PEC MOFA views) and the top200 loading cut.",
        "  The other three variants (all_signed cut, unrestricted universe) exist as figures and",
        "  CSVs in the same folder and can be added on request.",
        "",
        "THE THREE SETS",
        "  niche20 UP/DOWN   MAST crescent-PEC vs other-PEC, |avg_log2FC| > 0.25 & p_val_adj < 0.05.",
        "                    Positive avg_log2FC = higher in the crescent niche.",
        "  PT F1 / PEC F1    The 200 largest-|loading| genes of that MOFA view for Factor 1, split by",
        "                    sign. UP is compared to positive loadings, DOWN to negative loadings.",
        "",
        "READING THE 'region' COLUMN",
        "  all three         in the niche20 signature AND in both views' top-200 Factor 1 loadings",
        "  PT + niche20      niche20 signature + PT top-200 loadings (not PEC)",
        "  PEC + niche20     niche20 signature + PEC top-200 loadings (not PT)",
        "  niche20 only      in the niche20 signature, in NEITHER view's top-200 loadings",
        "  PT + PEC          both views' top-200 loadings, not in the niche20 signature",
        "  PT only / PEC only   that view's top-200 loadings alone",
        "",
        "IMPORTANT — WHAT 'niche20 only' DOES NOT MEAN",
        "  It does NOT mean the gene is absent from Factor 1. Every gene here has a nonzero loading",
        "  in both views; these simply did not clear the top-200 |loading| cut. Some miss it narrowly",
        "  (CASC15 is PT rank 34, NRP1 rank 59, PVT1 rank 91 of 5043). The PT_F1_loading /",
        "  PT_F1_absrank / PEC_F1_* columns are included so this is visible per gene: rank 1 = largest",
        "  |loading| in that view, and anything past the top 200 of a view falls outside its circle.",
        "  26 of the 39 'niche20 only' UP genes carry a sign-matched POSITIVE loading in both views —",
        "  they move into the shared region under the all_signed cut, which is why that triple overlap",
        "  grows from 6 to 44.",
        "",
        "COLUMNS",
        "  avg_log2FC, p_val_adj, pct_niche20 (pct.1), pct_other_PEC (pct.2) — from the MAST DE table.",
        "  PT_F1_loading / PEC_F1_loading   signed Factor 1 loading in that view.",
        "  PT_F1_absrank / PEC_F1_absrank   rank by |loading| within that view (PT n=5043, PEC n=2089).",
        "",
        "SOURCES",
        "  DE       niche20_pec/DE/niche20_vs_others_pecOnly/markers_MAST_niche20_vs_others_pecOnly.csv",
        "  loadings review_analyses/loadings~hvg~5.csv (Factor1, views PT and PEC)",
        "  script   review_analyses/venn_niche20_mofa_factor1.py -> export_venn_tables_xlsx.py",
        "",
        "Species human (GRCh38), HGNC symbols. Enrichment statistics are NOT in this workbook —",
        "see overlap_stats_all_variants.tsv in the parent folder.",
    ]})

    with pd.ExcelWriter(OUT, engine="openpyxl") as xw:
        notes.to_excel(xw, sheet_name="README", index=False, header=False)
        counts.to_excel(xw, sheet_name="region_counts", index=False)
        tabs["UP"].to_excel(xw, sheet_name="3way_UP", index=False)
        tabs["DOWN"].to_excel(xw, sheet_name="3way_DOWN", index=False)

        wb = xw.book
        from openpyxl.styles import Font, Alignment, PatternFill
        from openpyxl.utils import get_column_letter

        hdr_font = Font(bold=True, color="FFFFFF")
        hdr_fill = PatternFill("solid", fgColor="404040")

        ws = wb["README"]
        ws.column_dimensions["A"].width = 105
        for row in ws.iter_rows():
            for c in row:
                c.alignment = Alignment(vertical="top")
        # bold the section headings — every body line in `notes` is indented, headings are not
        for r in range(1, ws.max_row + 1):
            v = ws.cell(row=r, column=1).value
            if v and not v.startswith(" "):
                ws.cell(row=r, column=1).font = Font(bold=True)

        for name in ("region_counts", "3way_UP", "3way_DOWN"):
            ws = wb[name]
            for c in ws[1]:
                c.font = hdr_font
                c.fill = hdr_fill
                c.alignment = Alignment(horizontal="center")
            ws.freeze_panes = "A2"
            ws.auto_filter.ref = ws.dimensions
            for i, col in enumerate(ws.iter_cols(), start=1):
                width = max(len(str(c.value)) if c.value is not None else 0 for c in col)
                ws.column_dimensions[get_column_letter(i)].width = min(max(width + 2, 10), 16)
            if name.startswith("3way"):
                for row in ws.iter_rows(min_row=2):
                    row[2].number_format = "0.000"          # avg_log2FC
                    row[3].number_format = "0.00E+00"       # p_val_adj
                    row[4].number_format = "0.000"
                    row[5].number_format = "0.000"
                    row[6].number_format = "0.000"          # PT loading
                    row[8].number_format = "0.000"          # PEC loading

    print(f"wrote {OUT}")
    print(counts.to_string(index=False))


if __name__ == "__main__":
    build()
