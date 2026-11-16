"""Helpers shared by the entry scripts of the PT -> aPT snATAC / scMEGA eGRN panels (ED Fig. 3j-m, Supp. Fig. 3a-d).

- SCRIPTS: analyses/10_atac_pt_apt_egrn/ATAC/scripts (the shipped analysis code; file names start with a digit, so they
  are loaded by path)
- load_script(name): import a shipped script as a module. Module-level code runs exactly as when the script is imported
  (for 45-48 it only defines paths, constants and functions; nothing is computed or written at import)
- run_main(name, argv, out): import the script, point its module-level output folder OUT (a release patch whose default
  is the original folder, see build/patch_rules/fig_atac_snrna.py) at `out`, set sys.argv to the original command line
  and call main(). Inputs are still read from the analysis workspace.
"""
import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "analyses" / "10_atac_pt_apt_egrn" / "ATAC" / "scripts"
if str(REPO / "config") not in sys.path:
    sys.path.insert(0, str(REPO / "config"))


def load_script(name):
    """Import analyses/10_atac_pt_apt_egrn/ATAC/scripts/<name> as a module and return it."""
    path = SCRIPTS / name
    modname = "kb_atac_" + path.stem.replace("-", "_")
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    spec.loader.exec_module(mod)
    return mod


def run_main(name, argv, out):
    """Run main() of a shipped script with its outputs redirected to the folder `out` (created if needed)."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    mod = load_script(name)
    mod.OUT = out
    saved = sys.argv
    sys.argv = [str(SCRIPTS / name)] + list(argv)
    try:
        mod.main()
    finally:
        sys.argv = saved
    return mod
