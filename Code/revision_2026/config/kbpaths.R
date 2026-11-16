# Data locations used by every R script in this repository (see config/env.sh; defaults = Aachen HPC locations).
# Usage at the top of a script:  source(file.path(Sys.getenv("KB_REPO"), "config", "kbpaths.R"))
.kb_get <- function(name, default) { v <- Sys.getenv(name, unset = ""); if (nzchar(v)) v else default }
KB_REPO       <- .kb_get("KB_REPO", NA_character_)
if (is.na(KB_REPO)) stop("KB_REPO is not set: run `source config/env.sh` in the repository first")
KB_PROJECT    <- .kb_get("KB_PROJECT", "/data/ep924610/kidney_biopsy")
KB_ANALYSIS   <- .kb_get("KB_ANALYSIS", "/data/ep924610/kidney_biopsy/review_analyses")
KB_DATA       <- .kb_get("KB_DATA", "/data/Kramann_lab/kidney_biopsy/data")
KB_RESULTS    <- .kb_get("KB_RESULTS", "/data/Kramann_lab/kidney_biopsy/results")
KB_XENIUM_RAW <- .kb_get("KB_XENIUM_RAW", "/data/Kramann_lab/From_Sidrah/new_data")
CONDA_ROOT    <- .kb_get("CONDA_ROOT", "/data/ep924610/miniforge3")
KB_FIG_OUT    <- .kb_get("KB_FIG_OUT", file.path(KB_REPO, "figure_output"))
KB_CORE_KEY   <- .kb_get("KB_CORE_KEY", "")

fig_out <- function(figure, filename) {
  d <- file.path(KB_FIG_OUT, figure)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.path(d, filename)
}
