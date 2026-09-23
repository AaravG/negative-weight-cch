"""Export a graph to flat binary files for the C++ implementation.

    python export_graph.py NY  cpp/data/NY
    python export_graph.py USA cpp/data/USA      # needs USA_CACHE (see run_nb.py)

Writes meta.txt (n, m), off.i64, tgt.i32, x.f64, y.f64 and w_<metric>.f64.
"""
import sys
from pathlib import Path

import numpy as np

from run_nb import load


def main():
    name = sys.argv[1]
    out = Path(sys.argv[2] if len(sys.argv) > 2 else f"cpp/data/{name}")
    out.mkdir(parents=True, exist_ok=True)
    n, off, tgt, x, y, ws = load(name)
    off.astype(np.int64).tofile(out / "off.i64")
    tgt.astype(np.int32).tofile(out / "tgt.i32")
    x.astype(np.float64).tofile(out / "x.f64")
    y.astype(np.float64).tofile(out / "y.f64")
    z = None
    if name != "USA":
        import dimacs
        for metric in ("ev", "ev_real"):
            try:
                gz = dimacs.load(name, metric, seed=1).z
            except Exception:
                continue
            if gz is not None:
                np.asarray(gz, np.float64).tofile(out / f"z_{metric}.f64")
    if name == "USA":
        import os
        from pathlib import Path as _P
        cache = _P(os.environ.get("USA_CACHE", "data/cache/USA"))
        if (cache / "ev_z.bin").exists():
            z = np.fromfile(cache / "ev_z.bin", np.float64)
    else:
        import dimacs
        g = dimacs.load(name, "ev", seed=1)
        z = np.array(g.z) if g.z is not None else None
    if z is not None:
        z.astype(np.float64).tofile(out / "z.f64")
    names = []
    for w, (wt, p) in ws.items():
        wt.astype(np.float64).tofile(out / f"w_{w}.f64")
        names.append(w)
    (out / "meta.txt").write_text(f"{n}\n{int(off[n])}\n{' '.join(names)}\n")
    print(f"{name}: n={n:,} m={int(off[n]):,} metrics={names} z={z is not None} -> {out}")


if __name__ == "__main__":
    main()
