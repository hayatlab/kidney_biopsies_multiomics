# The two niche runs of analyses/12_xenium_niche_cci_liana_crosstalker (sourced by the scripts of this folder and by the
# figure entry scripts). Needs config/kbpaths.R sourced first.
# Set KB_CCI_OUT to write the LIANA / CrossTalkeR output of a re-run to another folder (default: the folder of the
# published run, next to the input object).
.cci_base <- file.path(KB_RESULTS, "xenium", "sidrah")
CCI_RUNS <- list(
  niche20 = list(                                   # Fig. 6d, 6e: glomerular niche 20 (manual-glomerulus model, Leiden 1.0)
    dir    = file.path(.cci_base, "whole_data_onlyglom_20"),
    object = "xenium_with_niches_glom10.rds",       # Seurat object: the 241 cells of niche 20
    label  = "Niche20",
    report = "overall_glom04_niche20.html"),
  niche1 = list(                                    # Extended Data Fig. 4k: whole-core niche 1 (Leiden 0.4)
    dir    = file.path(.cci_base, "overall_glom_04_level2"),
    object = "xenium_niches_overall_with_glom_niche1.rds",   # Seurat object: the 35,773 cells of niche 1
    label  = "Niche1",
    report = "overall_glom04_niche1.html"))

cci_run <- function(name) {
  if (!name %in% names(CCI_RUNS)) stop("run must be one of: ", paste(names(CCI_RUNS), collapse = ", "))
  r <- CCI_RUNS[[name]]
  r$name <- name
  r$out  <- Sys.getenv("KB_CCI_OUT", unset = file.path(r$dir, "CCI_CrossTalkeR"))   # LIANA result, CrossTalkeR object and report
  r$liana <- file.path(r$out, sub("\\.rds$", "_liana.rds", r$object))
  r
}

# Tables exported for the figure scripts (03_export_tables.R)
CCI_EXPORT <- Sys.getenv("KB_CCI_EXPORT", unset = file.path(KB_ANALYSIS, "sidrah_scripts_20261120", "outputs", "cci"))
