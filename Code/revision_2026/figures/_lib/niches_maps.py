"""Shared helpers for the niche-map and niche-heatmap entry scripts (Fig. 6a-c, Extended Data Fig. 4d-e and g-j,
Supplementary Fig. 13a).

The shipped NicheCompass modules import each other by module name (replot_partA and regenerate_metadata_figures
import `nichecompass_5regions_downstream as ds` and drive it through module globals such as ds.NICHE_KEY), so their
folders are put on sys.path and the modules are imported normally, never loaded twice under different names.
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
A01 = REPO / "analyses" / "01_xenium_preprocessing_niches"
NR = A01 / "nichecompass_replot_level2"


def use_niche_modules() -> None:
    """Put config/, analyses/01_xenium_preprocessing_niches/ and its nichecompass_replot_level2/ on sys.path."""
    for p in (NR, A01, REPO / "config"):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))


def core_targets(figdir: str, panel: str, cores: dict) -> dict:
    """Arguments for replot_boundaries_level2.run(cores=...).

    cores: {public core id (XEN###): sample id (K-number)}. Returns {internal core label: {"label": "XEN093 (K50)",
    "out_base": $KB_FIG_OUT/<figdir>/<panel>_K50_XEN093}}. The internal label (it carries a pathology number) comes
    from the private core key through config/cores.py and is never printed or written; without the key it is the
    public id itself, which is how the core is named in the de-identified Zenodo object.
    """
    use_niche_modules()
    from cores import core_label
    from kbpaths import fig_out
    return {core_label(xen): {"label": f"{xen} ({k})", "out_base": fig_out(figdir, f"{panel}_{k}_{xen}")}
            for xen, k in cores.items()}
