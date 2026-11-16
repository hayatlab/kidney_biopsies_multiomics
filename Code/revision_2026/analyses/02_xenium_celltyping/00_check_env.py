"""Smoke-check environment + inputs. Writes a version block to RUN_LOG.md."""
import sys, platform, datetime, hashlib
import h5py
import config as C


def _decode(x):
    return x.decode() if isinstance(x, bytes) else x


def _var_names(path):
    with h5py.File(path, "r") as f:
        var = f["var"]
        ik = _decode(var.attrs.get("_index", "_index"))
        return [_decode(v) for v in var[ik][:]]


def main():
    C.ensure_dirs()
    import scvi, scanpy, anndata, numpy, sklearn, torch
    assert C.SC_REF.exists(), f"missing snRNA ref: {C.SC_REF}"
    assert C.XENIUM.exists(), f"missing Xenium: {C.XENIUM}"

    xe = _var_names(C.XENIUM)
    sc = set(_var_names(C.SC_REF))
    shared = [g for g in xe if g in sc]   # Xenium order is master
    assert len(shared) == len(xe) == 468, (len(shared), len(xe))

    import csv
    with open(C.SHARED_GENES_CSV, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["gene"])
        for g in shared:
            w.writerow([g])

    block = [
        f"## Env check {datetime.datetime.now().isoformat(timespec='seconds')}",
        f"- python {platform.python_version()} ({sys.executable})",
        f"- scvi-tools {scvi.__version__}, scanpy {scanpy.__version__}, "
        f"anndata {anndata.__version__}, sklearn {sklearn.__version__}, "
        f"numpy {numpy.__version__}, torch {torch.__version__} (cuda={torch.cuda.is_available()})",
        f"- shared genes: {len(shared)} / Xenium {len(xe)} (all panel genes present in snRNA ref)",
        "",
    ]
    with open(C.BASE / "RUN_LOG.md", "a") as fh:
        fh.write("\n".join(block) + "\n")
    print("ENV OK:", block[2])
    print("shared genes written:", C.SHARED_GENES_CSV)


if __name__ == "__main__":
    main()
