"""STEP: regenerate the hallmark "dedup" heatmap + dotplot for every niche version.

Why (2026-07-22): `glom_postprocess.hallmark_pathway_dotplot()` used to run decoupler
ULM straight on `adata.X`, which on the saved NicheCompass objects is RAW COUNTS, while
the "normal" matrixplot (`nichecompass_5regions_downstream.pathway_activity()`) scores
counts -> normalize_total(1e4) -> log1p. The two therefore ranked hallmark pathways
differently and the dedup figure silently dropped pathways the normal figure shows
(on overall 0.4: TNFA_SIGNALING_VIA_NFKB, INFLAMMATORY_RESPONSE, EMT, ANGIOGENESIS were
replaced by KRAS_SIGNALING_UP, XENOBIOTIC_METABOLISM, MTORC1_SIGNALING, APICAL_JUNCTION).
`glom_postprocess.py` now log-normalises first and sorts the per-niche ranking explicitly,
so the dedup view is the same underlying ranking as the normal one.

Per run: ULM is computed ONCE, figures written to the first figures/ dir, then copied to
the version's other figures/ dirs. Regenerated files per dir:
  hallmark_marker_pathways_dedup_heatmap.{pdf,png}
  hallmark_marker_pathways_dotplot.{pdf,png}
  hallmark_score_means.csv
  hallmark_top3_per_niche.csv           (new; which pathways made each niche's top-N)

Usage
-----
python -u nichecompass_replot_level2/regenerate_hallmark_dedup.py --config all
python -u nichecompass_replot_level2/regenerate_hallmark_dedup.py --config manual_glom10 \
        --dest-root <scratch>          # smoke test, writes nothing into the results tree
"""
from __future__ import annotations

import argparse
import shutil
import sys
import traceback
from pathlib import Path

import scanpy as sc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
from regenerate_metadata_figures import CONFIGS          # noqa: E402
import glom_postprocess as gp                            # noqa: E402

N_MARKERS = 3
COPY_FILES = [
    "hallmark_marker_pathways_dedup_heatmap.pdf",
    "hallmark_marker_pathways_dedup_heatmap.png",
    "hallmark_marker_pathways_dedup_heatmap_nostars.pdf",
    "hallmark_marker_pathways_dedup_heatmap_nostars.png",
    "hallmark_marker_pathways_dotplot.pdf",
    "hallmark_marker_pathways_dotplot.png",
    "hallmark_score_means.csv",
    f"hallmark_top{N_MARKERS}_per_niche.csv",
]


def run_config(name: str, cfg: dict, dest_root: str | None = None) -> None:
    print(f"\n{'=' * 70}\n[{name}] {cfg['adata'].name}  niche={cfg['niche_key']}\n{'=' * 70}",
          flush=True)
    # staging mode: compute once into <dest_root>/<config>/figures (see
    # regenerate_metadata_figures.run_config for why), else the real figures/ dirs
    dirs = [Path(dest_root) / name / "figures"] if dest_root else list(cfg["dirs"])
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(cfg["adata"])
    print(f"[load] {adata.shape}", flush=True)
    if cfg["niche_key"] not in adata.obs:
        print(f"  [SKIP] {cfg['niche_key']} not in obs", flush=True)
        return

    top = gp.hallmark_pathway_dotplot(adata, niche_key=cfg["niche_key"],
                                      out_dir=dirs[0], n_markers=N_MARKERS)
    if top:
        for niche in sorted(top, key=lambda x: (0, int(x)) if str(x).isdigit() else (1, x)):
            print(f"    niche {niche}: {top[niche]}", flush=True)

    for d in dirs[1:]:
        for fn in COPY_FILES:
            src = dirs[0] / fn
            if src.exists():
                shutil.copy2(src, d / fn)
        print(f" -> copied to {d}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="all",
                    choices=list(CONFIGS) + ["all"])
    ap.add_argument("--dest-root", default="", help="scratch base for smoke tests")
    args = ap.parse_args()
    todo = list(CONFIGS) if args.config == "all" else [args.config]
    for n in todo:
        try:
            run_config(n, CONFIGS[n], dest_root=args.dest_root or None)
        except Exception as exc:                                   # noqa: BLE001
            print(f"[FAIL] {n}: {exc!r}", flush=True)
            traceback.print_exc()
    print("\n[ALL DONE] " + ", ".join(todo), flush=True)


if __name__ == "__main__":
    main()
