"""Shared paths and constants for the Xenium<->snRNA concordance pipeline."""
from kbpaths import KB_ANALYSIS, KB_DATA, KB_RESULTS  # repository path config (config/kbpaths.py)
from pathlib import Path

SEED = 42

BASE = Path(f"{KB_ANALYSIS}/celltype_concordance")
OUT = BASE / "outputs"
FIG = BASE / "figures"
MET = OUT / "metrics"
MODELS = OUT / "models"

# Read-only inputs
SC_REF = Path(f"{KB_DATA}/charlotte_new_data/adata_20260717.h5ad")
XENIUM = Path(f"{KB_RESULTS}/xenium/leon/new_data/"
              "overall_with_glom/xenium_with_niches_overall_with_glom.h5ad")

# Derived artifacts
HUMANIZED_REF = OUT / "humanized_ref.h5ad"      # snRNA subset to 468 genes, integer counts
XENIUM_LABELED = OUT / "xenium_labeled.h5ad"     # Xenium + predictions + leiden sweep
SHARED_GENES_CSV = OUT / "shared_genes.csv"      # the 468 panel genes, Xenium order

LEIDEN_RES = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0]
SC_BATCH = "sample"
XEN_BATCH = "sample"
LABELS = ["level1", "level2"]

# Existing Xenium cell-typing columns to compare against
XEN_CURRENT = ["cell_type_manual", "predic"]


def ensure_dirs():
    for d in (OUT, FIG, MET, MODELS):
        d.mkdir(parents=True, exist_ok=True)
