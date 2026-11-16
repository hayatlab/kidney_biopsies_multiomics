# review_analyses/celltype_concordance/metrics.py
"""Pure concordance metrics between two label arrays / clusterings."""
import numpy as np
import pandas as pd
from sklearn.metrics import (adjusted_rand_score, normalized_mutual_info_score,
                             homogeneity_completeness_v_measure)


def clustering_scores(labels_true, labels_pred):
    a = np.asarray(labels_true).astype(str)
    b = np.asarray(labels_pred).astype(str)
    h, c, v = homogeneity_completeness_v_measure(a, b)
    return {
        "ari": float(adjusted_rand_score(a, b)),
        "nmi": float(normalized_mutual_info_score(a, b)),
        "homogeneity": float(h), "completeness": float(c), "v_measure": float(v),
    }


def cluster_purity(clusters, labels):
    # Preserve original cluster dtype so callers can filter with df["cluster"]==original_value.
    # Only cast labels to str for stable groupby comparison.
    clusters_arr = np.asarray(clusters)
    labels_arr = np.asarray(labels).astype(str)
    df = pd.DataFrame({"cluster": clusters_arr, "label": labels_arr})
    rows = []
    for cl, g in df.groupby("cluster", observed=True):
        vc = g["label"].value_counts()
        rows.append({"cluster": cl, "dominant_label": vc.index[0],
                     "purity": float(vc.iloc[0] / len(g)), "size": int(len(g))})
    out = pd.DataFrame(rows).sort_values("cluster").reset_index(drop=True)
    overall = float((out["purity"] * out["size"]).sum() / out["size"].sum())
    return out, overall


def crosstab_normalized(a, b, normalize="index"):
    ct = pd.crosstab(pd.Series(np.asarray(a).astype(str), name="row"),
                     pd.Series(np.asarray(b).astype(str), name="col"))
    if normalize == "index":
        ct = ct.div(ct.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    elif normalize == "columns":
        ct = ct.div(ct.sum(axis=0).replace(0, np.nan), axis=1).fillna(0.0)
    return ct


def bootstrap_ari(a, b, n=200, seed=42):
    a = np.asarray(a); b = np.asarray(b)
    rng = np.random.default_rng(seed)
    N = len(a); vals = np.empty(n)
    for i in range(n):
        idx = rng.integers(0, N, N)
        vals[i] = adjusted_rand_score(a[idx].astype(str), b[idx].astype(str))
    return float(vals.mean()), float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))
