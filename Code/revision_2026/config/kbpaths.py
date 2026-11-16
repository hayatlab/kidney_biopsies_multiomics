"""Data locations used by every Python script in this repository.

Values come from environment variables (set by config/env.sh); the defaults are the Aachen HPC locations where the
analyses were run. All values are plain strings so that they behave exactly like the original hard-coded paths.
"""
import os

_DEFAULTS = {
    "KB_PROJECT": "/data/ep924610/kidney_biopsy",
    "KB_ANALYSIS": "/data/ep924610/kidney_biopsy/review_analyses",
    "KB_DATA": "/data/Kramann_lab/kidney_biopsy/data",
    "KB_RESULTS": "/data/Kramann_lab/kidney_biopsy/results",
    "KB_XENIUM_RAW": "/data/Kramann_lab/From_Sidrah/new_data",
    "CONDA_ROOT": "/data/ep924610/miniforge3",
}

KB_REPO = os.environ.get("KB_REPO", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
KB_PROJECT = os.environ.get("KB_PROJECT", _DEFAULTS["KB_PROJECT"])
KB_ANALYSIS = os.environ.get("KB_ANALYSIS", _DEFAULTS["KB_ANALYSIS"])
KB_DATA = os.environ.get("KB_DATA", _DEFAULTS["KB_DATA"])
KB_RESULTS = os.environ.get("KB_RESULTS", _DEFAULTS["KB_RESULTS"])
KB_XENIUM_RAW = os.environ.get("KB_XENIUM_RAW", _DEFAULTS["KB_XENIUM_RAW"])
CONDA_ROOT = os.environ.get("CONDA_ROOT", _DEFAULTS["CONDA_ROOT"])
KB_FIG_OUT = os.environ.get("KB_FIG_OUT", os.path.join(KB_REPO, "figure_output"))
KB_CORE_KEY = os.environ.get("KB_CORE_KEY", "")


def fig_out(figure: str, filename: str) -> str:
    """Path for a figure panel file, e.g. fig_out("Fig3", "Fig3c_factor1_vs_egfr.pdf"); creates the folder."""
    d = os.path.join(KB_FIG_OUT, figure)
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, filename)
