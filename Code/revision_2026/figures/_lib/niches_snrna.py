"""Helpers shared by the entry scripts of the snRNA-seq PEC / immune subcluster panels (Fig. 6g-l, Supp. Fig. 19c-f).

- REPO: repository root; ANALYSES_09: analyses/09_snrna_pec_immune_subclusters (the shipped analysis code)
- load_script(path): import a shipped analysis script as a module by its path (several names start with a digit
  or are not on sys.path); module-level code runs once, exactly as when the script is imported
- CORES: the two Xenium cores shown in Fig. 6g, as public core id -> public sample id (K-number). The internal core
  label is obtained with config/cores.py::core_label (private key) and is never written to titles or file names.
"""
import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
ANALYSES_09 = REPO / "analyses" / "09_snrna_pec_immune_subclusters"
if str(REPO / "config") not in sys.path:
    sys.path.insert(0, str(REPO / "config"))

# Fig. 6g: glomerular niche 20 / crescent core (IgAN) and a TIN core with a "healthy" glomerulus
CORES = {"XEN093": "K50", "XEN116": "K18"}


def load_script(path, name=None):
    """Import the analysis script at `path` (str or Path) as a module and return it."""
    path = Path(path)
    name = name or "kb_" + path.stem.replace("-", "_")
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def core_display(xen):
    """Public text for titles / file names of a core: 'XEN093 K50' and 'K50_XEN093'."""
    k = CORES[xen]
    return f"{xen} {k}", f"{k}_{xen}"
