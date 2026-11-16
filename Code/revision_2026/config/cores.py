"""Xenium core identities without patient identifiers.

The internal core labels (`TMA_ID` in the analysed objects) contain pathology numbers, so they are never written into
this repository. Scripts refer to individual cores by their public id (`xenium_core_id`, XEN###, the id used in the
de-identified Zenodo object) and read lists of excluded cores from a private file.

- core_label("XEN085") returns the internal label when the private key is available (KB_CORE_KEY, a CSV with columns
  xenium_core_id and old_label); otherwise it returns the public id itself, which is how the core is named in the
  Zenodo object.
- excluded_core_labels(tag) returns the internal labels of excluded cores (KB_EXCLUDED_CORES, a TSV with columns tag and
  label). Without the private file it returns only the entries that carry no identifier ("Unassigned" and the Apollo d0
  A2 control core); the excluded cores are not part of the Zenodo object, so nothing else needs removing there.
"""
import csv
import os

from kbpaths import KB_CORE_KEY

KB_EXCLUDED_CORES = os.environ.get("KB_EXCLUDED_CORES", "")
_PUBLIC_EXCLUDED = {"nichecompass_blacklist": ["Unassigned", "Appollo d0 A2 CTRL"], "dropped_low_cells": []}


def _norm(s) -> str:
    return " ".join(str(s).split())


def core_label(xenium_core_id: str, available=None) -> str:
    """Internal label of a public core id. Pass `available` (e.g. adata.obs["TMA_ID"].unique()) to get the label exactly
    as spelled in that object: a few labels differ from the key only by whitespace."""
    if not KB_CORE_KEY:
        return xenium_core_id
    with open(KB_CORE_KEY, newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if r["xenium_core_id"] == xenium_core_id]
    if len(rows) != 1:
        raise KeyError(f"{xenium_core_id}: {len(rows)} rows in the core key")
    label = rows[0]["old_label"]
    if available is not None:
        hits = [a for a in available if _norm(a) == _norm(label)]
        if len(hits) != 1:
            raise KeyError(f"{xenium_core_id}: {len(hits)} matching labels in the object")
        label = hits[0]
    return label


def excluded_core_labels(tag: str) -> list:
    if not KB_EXCLUDED_CORES:
        return list(_PUBLIC_EXCLUDED[tag])
    with open(KB_EXCLUDED_CORES, newline="") as fh:
        labels = [r["label"] for r in csv.DictReader(fh, delimiter="\t") if r["tag"] == tag]
    if not labels:
        raise KeyError(f"no excluded cores with tag {tag!r} in {KB_EXCLUDED_CORES}")
    return labels
