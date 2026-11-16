"""Push a staged figure set into the real results tree, verifying every write.

Why staging exists (2026-07-22): the Kramann_lab group quota is at its 60 TB limit, so
writing figures straight into `/data/Kramann_lab/.../results/xenium/leon/` fails
mid-`savefig` with `OSError(122) Disk quota exceeded` and leaves TRUNCATED png/pdf files
behind. So the regeneration scripts compute into `<stage>/<config>/figures` under the
user's own quota, and this script copies them in afterwards — one short pass, with an
integrity check on every destination file, so a quota failure is loud instead of silent.

Usage
-----
python -u nichecompass_replot_level2/sync_staged_figures.py --stage <stage_root>
python -u nichecompass_replot_level2/sync_staged_figures.py --stage <stage_root> --dry-run
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
from regenerate_metadata_figures import CONFIGS   # noqa: E402


def intact(p: Path) -> tuple[bool, str]:
    """Cheap corruption check: truncated png/pdf are the failure mode here."""
    try:
        if p.stat().st_size == 0:
            return False, "zero bytes"
        if p.suffix == ".png":
            from PIL import Image
            Image.open(p).verify()
        elif p.suffix == ".pdf":
            if b"%%EOF" not in p.read_bytes()[-2048:]:
                return False, "pdf missing %%EOF"
    except Exception as exc:                                       # noqa: BLE001
        return False, f"{type(exc).__name__}: {exc}"
    return True, ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, help="staging root used by --dest-root")
    ap.add_argument("--config", default="all", choices=list(CONFIGS) + ["all"])
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    stage = Path(args.stage)
    todo = list(CONFIGS) if args.config == "all" else [args.config]
    n_copied = n_failed = 0
    failures: list[str] = []

    for name in todo:
        src_dir = stage / name / "figures"
        if not src_dir.is_dir():
            print(f"[{name}] no staged figures at {src_dir} — skipping", flush=True)
            continue
        srcs = sorted(p for p in src_dir.iterdir() if p.is_file())
        for dest_dir in CONFIGS[name]["dirs"]:
            print(f"[{name}] {len(srcs)} files -> {dest_dir}", flush=True)
            if args.dry_run:
                continue
            dest_dir.mkdir(parents=True, exist_ok=True)
            for s in srcs:
                d = dest_dir / s.name
                try:
                    shutil.copy2(s, d)
                except Exception as exc:                           # noqa: BLE001
                    n_failed += 1
                    failures.append(f"{d}: copy failed: {exc!r}")
                    continue
                good, why = intact(d)
                if good:
                    n_copied += 1
                else:
                    n_failed += 1
                    failures.append(f"{d}: {why}")

    print(f"\n[sync] copied+verified {n_copied}, FAILED {n_failed}", flush=True)
    for f in failures:
        print("  FAIL", f, flush=True)
    if failures:
        print("\n[sync] destination files listed above are missing or truncated — "
              "free space on the Kramann_lab quota and re-run this script.", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
