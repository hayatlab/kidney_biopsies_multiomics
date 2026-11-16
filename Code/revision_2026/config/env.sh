#!/usr/bin/env bash
# Source this file before running any analysis or figure script:   source config/env.sh
# Every data location the code uses is one of the variables below. The defaults are the locations on the
# Aachen HPC where the analyses were run; set a variable before sourcing to point the code elsewhere.

_kb_here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export KB_REPO="${KB_REPO:-$(dirname "${_kb_here}")}"          # this repository

# Inputs and analysis workspaces
export KB_PROJECT="${KB_PROJECT:-/data/ep924610/kidney_biopsy}"                          # project folder (snRNA objects, earlier NicheCompass runs)
export KB_ANALYSIS="${KB_ANALYSIS:-/data/ep924610/kidney_biopsy/review_analyses}"        # analysis workspaces (one sub-folder per analysis, see each analysis README)
export KB_DATA="${KB_DATA:-/data/Kramann_lab/kidney_biopsy/data}"                         # shared processed data (snRNA/snATAC objects, Xenium objects, large outputs)
export KB_RESULTS="${KB_RESULTS:-/data/Kramann_lab/kidney_biopsy/results}"                # shared results (Xenium NicheCompass objects)
export KB_XENIUM_RAW="${KB_XENIUM_RAW:-/data/Kramann_lab/From_Sidrah/new_data}"           # Xenium Onboard Analysis output bundles (read only)

# Outputs of the figure scripts in figures/
export KB_FIG_OUT="${KB_FIG_OUT:-${KB_REPO}/figure_output}"

# Optional private files (never part of the repository; see config/cores.py). Without them the code runs on the
# de-identified Zenodo object, where cores are named by their public id (XEN###).
export KB_CORE_KEY="${KB_CORE_KEY:-}"                 # CSV: xenium_core_id -> internal core label
export KB_EXCLUDED_CORES="${KB_EXCLUDED_CORES:-}"     # TSV: tag, internal label of cores excluded before NicheCompass

# Conda (the analyses use several environments; see envs/README.md)
export CONDA_ROOT="${CONDA_ROOT:-/data/ep924610/miniforge3}"

export PYTHONPATH="${KB_REPO}/config${PYTHONPATH:+:${PYTHONPATH}}"
unset _kb_here
