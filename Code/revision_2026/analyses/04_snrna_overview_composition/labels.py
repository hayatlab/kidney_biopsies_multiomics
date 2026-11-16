"""Shared display labels and colours (contract C3): `<out_root>/meta/labels.json`."""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

from .config import load_config
from .jsonio import dumps_strict, write_json

HEX = re.compile(r"^#[0-9A-Fa-f]{6}$")

LEVEL2_ORDER = ["PT", "aPT", "dPT", "TAL", "DCT", "CNT", "PC", "IC-A", "IC-B",
                "POD", "PEC", "EC", "FIB", "MES", "VSM-P", "Immune"]

DISEASE_SOURCE_MAP = {
    "CTRL": "Control", "control": "Control",
    "IgA": "IgAN", "IgAN": "IgAN",
    "DN": "DN", "FSGS": "FSGS", "MN": "MN", "TIN": "TIN",
    "control_AKI": "control_AKI",
}
DISEASE_ORDER = ["Control", "DN", "FSGS", "IgAN", "MN", "TIN", "control_AKI"]
# Light, distinct hues for the dark UI (after Paul Tol "light"/"bright"). The first choice, Tol "muted",
# made Control (#332288) and FSGS (#117733) near-invisible on the #0a0a0a background.
DISEASE_COLORS = {
    "Control": "#C0C7D0",
    "DN": "#EECC66",
    "FSGS": "#44BB99",
    "IgAN": "#EE6677",
    "MN": "#AA88FF",
    "TIN": "#BBCC33",
    "control_AKI": "#77AADD",
}

PROGRESSION_ORDER = ["Fast Progressor", "stable", "improved", "unclassified"]
# Okabe-Ito (colour-blind safe); grey for "no class". ("not in corrected table" was dropped on 2026-09-14 when
# Xenium cores outside the corrected table got the jm_slope rule.)
PROGRESSION_COLORS = {
    "Fast Progressor": "#D55E00",
    "stable": "#56B4E9",
    "improved": "#009E73",
    "unclassified": "#999999",
}

FACTORS = ["Factor1", "Factor2", "Factor3", "Factor4", "Factor5"]

BASE_KEYS = ("level2", "disease", "progression", "factors")


def display_disease(s: str) -> str:
    """Unified display label (CTRL/control -> Control, IgA/IgAN -> IgAN). Unknown labels raise KeyError."""
    try:
        return DISEASE_SOURCE_MAP[s]
    except (KeyError, TypeError):
        raise KeyError(f"unknown disease label {s!r}; known: {sorted(DISEASE_SOURCE_MAP)}") from None


def level2_colors(palette_path: str | Path | None = None) -> dict[str, str]:
    """level2 -> hex colour from the concordance legend CSV (columns group,color), in LEVEL2_ORDER."""
    path = Path(palette_path) if palette_path is not None else Path(load_config()["sources"]["xenium_celltype_legend"])
    colors: dict[str, str] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames or not {"group", "color"} <= set(reader.fieldnames):
            raise ValueError(f"{path}: expected columns group,color; got {reader.fieldnames}")
        for row in reader:
            group, color = row["group"].strip(), row["color"].strip()
            if not HEX.match(color):
                raise ValueError(f"{path}: invalid colour {color!r} for {group!r}")
            if group in colors:
                raise ValueError(f"{path}: duplicate group {group!r}")
            colors[group] = color
    missing = [g for g in LEVEL2_ORDER if g not in colors]
    unexpected = [g for g in colors if g not in LEVEL2_ORDER]
    if missing or unexpected:
        raise ValueError(f"{path}: palette missing {missing}, unexpected {unexpected}")
    return {g: colors[g] for g in LEVEL2_ORDER}


def _validate_entry(key: str, value) -> None:
    if isinstance(value, dict) and "order" in value:
        order, colors = value["order"], value.get("colors")
        if not isinstance(order, list) or not all(isinstance(o, str) for o in order):
            raise ValueError(f"labels[{key!r}].order must be a list of strings")
        if not isinstance(colors, dict):
            raise ValueError(f"labels[{key!r}].colors must be a mapping")
        missing = [o for o in order if o not in colors]
        if missing:
            raise ValueError(f"labels[{key!r}]: no colour for {missing}")
        bad = [c for c in colors.values() if not isinstance(c, str) or not HEX.match(c)]
        if bad:
            raise ValueError(f"labels[{key!r}]: invalid colours {bad}")
    dumps_strict(value)  # JSON-serialisable, no NaN


def build_labels(extra: dict | None = None, palette_path: str | Path | None = None) -> dict:
    """C3 structure; `extra` adds module-specific entries (e.g. niche_overall, niche_glom, level1)."""
    labels = {
        "level2": {"order": list(LEVEL2_ORDER), "colors": level2_colors(palette_path)},
        "disease": {"order": list(DISEASE_ORDER), "colors": dict(DISEASE_COLORS),
                    "source_map": dict(DISEASE_SOURCE_MAP)},
        "progression": {"order": list(PROGRESSION_ORDER), "colors": dict(PROGRESSION_COLORS)},
        "factors": list(FACTORS),
    }
    for key, value in (extra or {}).items():
        if key in BASE_KEYS:
            raise ValueError(f"extra may not override the shared entry {key!r}")
        _validate_entry(key, value)
        labels[key] = value
    return labels


def write_labels(out_root: str | Path, extra: dict | None = None, palette_path: str | Path | None = None) -> Path:
    """Write `<out_root>/meta/labels.json`, keeping entries previously added by other modules."""
    path = Path(out_root) / "meta" / "labels.json"
    existing = json.loads(path.read_text()) if path.exists() else {}
    merged = build_labels(extra, palette_path)
    for key, value in existing.items():
        if key not in merged:
            merged[key] = value
    return write_json(path, merged, indent=2)
