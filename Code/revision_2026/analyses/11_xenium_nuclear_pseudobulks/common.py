"""Shared helpers for the Xenium Zenodo preparation (review_analyses/xenium_zenodo_20261027).

Core identity: the internal table INTERNAL_Xenium_old_vs_new_labels_DO_NOT_SHARE.csv maps the object's TMA_ID (carries pathology
numbers) to `xenium_core_id` (XEN###) and the Zenodo `new_sample_id`. Nothing derived from TMA_ID may leave private/.
"""
from kbpaths import KB_ANALYSIS, KB_DATA, KB_XENIUM_RAW  # repository path config (config/kbpaths.py)
import re, json, os, csv
import numpy as np, pandas as pd
A = f"{KB_ANALYSIS}/xenium_zenodo_20261027"
OUT = f"{KB_DATA}/zenodo_xenium_20261027"           # outputs/ symlink target (Kramann_lab quota)
RAW = f"{KB_XENIUM_RAW}"                                  # original bundles, read-only here
OBJ = f"{KB_ANALYSIS}/celltype_concordance/outputs/xenium_labeled.h5ad"
INTERNAL = f"{KB_ANALYSIS}/INTERNAL_Xenium_old_vs_new_labels_DO_NOT_SHARE.csv"
PREV_KEY = f"{KB_ANALYSIS}/mcfa_xenium_pseudobulk_20261026/outputs/private/sample_key_PRIVATE.csv"
META = f"{KB_ANALYSIS}/celltype_concordance/outputs/per_sample_metadata_updated.csv"
RUNS = ["XETG00229__0041734__Region_1__20250416__104148", "XETG00229__0041734__Region_3__20250416__104148",
        "XETG00229__0042348__Region_1__20250416__104148", "XETG00229__0042348__Region_2__20250416__104148",
        "XETG00229__0042348__Region_3__20250416__104148"]                       # the 5 runs in the analysed object (0041734 Region_2 unused)
PIXEL_SIZE = 0.2125                                                             # um per pixel (experiment.xenium, all runs checked in 01)
BUFFER_UM = 20.0                                                                # polygon buffer around the convex hull of analysed cells
QV_MIN = 20.0                                                                   # 10x default; reproduces the object's counts layer exactly (checked 27 Oct)
norm = lambda s: re.sub(r"\s+", " ", str(s)).strip()
run_short = lambda r: re.sub(r"XETG00229__(\d+)__Region_(\d)__.*", r"\1_R\2", r)

def core_key():
    """Private core key: TMA_norm -> xenium_core_id, sample_id (Zenodo), axis label, prev_sample_id (26 Oct delivery)."""
    m = pd.read_csv(INTERNAL); m["TMA_norm"] = m["old_label"].map(norm)
    k = pd.read_csv(PREV_KEY); k["TMA_norm"] = k["TMA_ID"].map(norm)
    m = m.merge(k[["TMA_norm", "sample_id", "core_pseudonym"]].rename(columns={"sample_id": "prev_sample_id"}), on="TMA_norm", how="left")
    m["grid_position"] = m["new_axis_label"].str.split(" ").str[0]
    assert m.xenium_core_id.is_unique and m.new_sample_id.is_unique and m.prev_sample_id.notna().all(), "key not 1:1"
    return m.rename(columns={"new_sample_id": "sample_id", "figure_row_SuppFig4": "figure_row"})

def object_obs():
    import anndata as ad
    a = ad.read_h5ad(OBJ, backed="r"); o = a.obs.copy(); o["TMA_norm"] = o["TMA_ID"].astype(str).map(norm)
    o["run"] = o["sample"].astype(str); o["x"] = a.obsm["spatial"][:, 0]; o["y"] = a.obsm["spatial"][:, 1]
    k = core_key().set_index("TMA_norm")
    assert set(o.TMA_norm) == set(k.index), "object cores != key cores"
    for c in ["xenium_core_id", "sample_id", "prev_sample_id", "grid_position", "figure_row"]: o[c] = o.TMA_norm.map(k[c])
    return o

def core_polygons(o_run):
    """Convex hull of the analysed cells of each core, buffered by BUFFER_UM; checks the polygons do not overlap."""
    import shapely
    from shapely.geometry import MultiPoint
    polys = {}
    for cid, g in o_run.groupby("xenium_core_id", observed=True):
        polys[cid] = MultiPoint(np.c_[g.x.values, g.y.values]).convex_hull.buffer(BUFFER_UM)
    ids = list(polys); overlaps = []
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            inter = polys[ids[i]].intersection(polys[ids[j]]).area
            if inter > 0: overlaps.append((ids[i], ids[j], round(inter, 1)))
    return polys, overlaps

def assign_points(polys, x, y):
    """Vectorised point-in-polygon with a bbox prefilter; returns an object array of core ids ('' = outside every core)."""
    import shapely
    lab = np.full(len(x), "", dtype=object); hit = np.zeros(len(x), bool)
    for cid, p in polys.items():
        minx, miny, maxx, maxy = p.bounds
        bb = np.where((x >= minx) & (x <= maxx) & (y >= miny) & (y <= maxy))[0]
        if len(bb) == 0: continue
        inside = shapely.contains_xy(p, x[bb], y[bb])
        sel = bb[inside]
        assert not hit[sel].any(), f"point in two polygons ({cid})"
        lab[sel] = cid; hit[sel] = True
    return lab

def old_label_tokens(key):
    """Old TMA labels / Excel names that must not appear in public files; control cores' Excel names equal their new ids and are excluded."""
    new = set(key.sample_id.astype(str))
    return tuple(t for t in list(key.TMA_norm) + list(key.old_sample_name_excel.astype(str)) if t and t != "nan" and t not in new)

PRIVATE_RE = re.compile(r"H\d{2}-\d{3,}|\b\d{2}-\d{4,}\b|T\d{2}-\d{4,}|Rotterdam|patient_id|H\.number|\d{4}-\d{2}-\d{2}(?!T)|Combi |Appollo")
def privacy_scan(paths, extra_tokens=()):
    """Return {file: hits} for text files; h5ad obs columns are scanned by the caller. extra_tokens: the old TMA labels."""
    hits = {}
    for p in paths:
        if p.endswith((".h5ad", ".h5", ".parquet", ".tif", ".zip", ".tar", ".gz")): continue
        txt = open(p, errors="replace").read()
        m = set(PRIVATE_RE.findall(txt)) | {t for t in extra_tokens if t and t in txt}
        if m: hits[os.path.relpath(p, OUT)] = sorted(m)[:5]
    return hits
