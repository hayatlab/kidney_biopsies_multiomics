"""Helpers shared by the entry scripts of the Xenium factor / neighbourhood / CCI / cell-typing panels and the snRNA
NELL1 panels (Fig. 3c, 3f-i; ED 4c, 4f; ED 6a; Supp. Fig. 9a-d, 10d, 12a-b, 13a-b).

- REPO, ANALYSES: repository root and analyses/ (the shipped analysis code, never the analysis workspaces).
- load_script(relpath): import analyses/<relpath> as a module by its path (the file names start with a digit); its
  module-level code runs once, exactly as when the script is imported. on_path=True first puts the script's folder on
  sys.path, for scripts that import their siblings (02_xenium_celltyping: `import config as C`, `import metrics as M`).
- CORES / core_display(): the Xenium cores drawn in these panels, public core id -> public sample id. The internal
  core label (config/cores.py::core_label, private key) selects the cells; it contains a pathology number and is never
  written to a title, file name or log line - core_display() gives the public text instead.
- rename(folder, mapping): give files that a shipped function writes under a fixed name their panel name.
"""
import importlib.util
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
ANALYSES = REPO / "analyses"
if str(REPO / "config") not in sys.path:
    sys.path.insert(0, str(REPO / "config"))

CORES = {
    "XEN085": "K16",       # Fig. 3f: TIN, eGFR 5 (Utrecht1 E4)
    "XEN116": "K18",       # Fig. 3i right map, Supp. Fig. 12a: TIN, eGFR 23 (Utrecht2 F2)
    "XEN119": "K82",       # Fig. 3i left map: FSGS, eGFR 63 (Utrecht2 F5)
    "XEN008": "control",   # Supp. Fig. 12b: control core, eGFR 60 (Apollo B5)
}


def core_display(xen):
    """Public text for a core in titles and file names, e.g. 'XEN085 K16'."""
    return f"{xen} {CORES[xen]}"


def load_script(relpath, name=None, on_path=False):
    """Import analyses/<relpath> as a module and return it."""
    path = ANALYSES / relpath
    if on_path and str(path.parent) not in sys.path:
        sys.path.insert(0, str(path.parent))
    name = name or "kb_" + path.stem.replace("-", "_")
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def rename(folder, mapping):
    """os.replace(folder/old, folder/new) for every old -> new in mapping; returns the new paths."""
    folder = Path(folder)
    out = []
    for old, new in mapping.items():
        os.replace(folder / old, folder / new)
        out.append(folder / new)
    return out
