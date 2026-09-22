"""Battery-constrained, potential-free CCH (Numba) on NY, a region, or the full USA.

    python run_nb_battery.py NY
    python run_nb_battery.py USA       # uses data/rank_USA.npy from run_nb.py

EV-terrain weighting. For battery capacities M = 0.25 D, 1 D, 4 D (D = median
unconstrained energy of random routes): customization time, profile-size
histogram, memory, and queries with b0 = M and M/2 checked against a
label-correcting reference. Writes results/results_nb_battery_<name>.md.
"""
import sys
import time
from pathlib import Path

import numpy as np

import cch_nb as nb
import cch_nb_battery as nbb
from cch_memory import mem_gb
from run_nb import load, log

HERE = Path(__file__).parent


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "NY"
    pairs_n = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    verify = pairs_n > 0          # negative count: time queries, skip the reference
    pairs_n = abs(pairs_n)
    factors = [float(f) for f in sys.argv[3].split(",")] if len(sys.argv) > 3 else [0.25, 1.0, 4.0]
    suffix = "" if len(sys.argv) <= 3 else "_M" + "_".join(sys.argv[3].split(","))
    out = HERE / "results" / f"results_nb_battery_{name}{suffix}.md"
    n, off, tgt, x, y, ws = load(name)
    w = ws["ev"][0]
    del ws
    noff, nbrs, rev = nb.undirected_csr(n, off, tgt)
    rank_file = HERE / "data" / f"rank_{name}.npy"
    if rank_file.exists():
        rank = np.load(rank_file)
    else:
        rank = nb.nested_dissection(n, noff, nbrs, rev, x, y, 64)
        np.save(rank_file, rank)
    del rev
    first, head, parent = nb.symbolic(n, noff, nbrs, rank)
    del noff, nbrs
    e_arc, e_up = nb.edge_arcs(n, off, tgt, rank, first, head)
    log(f"{name}: n={n:,}, shortcut arcs {head.shape[0]:,}, memory {mem_gb()[0]:.1f} GB")

    # D from unconstrained (scalar) distances
    upw, dnw = nb.customize(n, first, head, e_arc, e_up, w)
    df = np.full(n, np.inf)
    dr = np.full(n, np.inf)
    mark = np.zeros(n, np.int64)
    rng = np.random.default_rng(21)
    dists = []
    q = 0
    while len(dists) < 50:
        q += 1
        s, t = int(rng.integers(n)), int(rng.integers(n))
        d = nb.query(rank[s], rank[t], first, head, parent, upw, dnw, df, dr, mark, q)
        if d < np.inf:
            dists.append(d)
    del upw, dnw, df, dr
    D = float(np.median(dists))
    pairs = [(int(rng.integers(n)), int(rng.integers(n))) for _ in range(pairs_n)]

    lines = [f"# Battery-constrained CCH (Numba) on {name}", "",
             f"n = {n:,}; shortcut arcs = {head.shape[0]:,} ({2 * head.shape[0]:,} directions); "
             f"D = {D:,.0f}. Every answer is checked against a label-correcting max-charge search.", "",
             "| Battery M | Customization | Pool pieces | Memory (GB) | Size 0 | Size 1 | 2–5 | 6–10 | 11–50 | 51–100 | >100 | Largest | Query (median) | Infeasible | Correct |",
             "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    soc = np.full(n, -np.inf)
    dsoc = np.full(n, -np.inf)
    mark = np.zeros(n, np.int64)
    for factor in factors:
        M = factor * D
        t0 = time.perf_counter()
        cnt, I1, C1, O1, pst, PI, PC, PO = nbb.battery_customize(n, first, head, e_arc, e_up, w, M)
        t_cust = time.perf_counter() - t0
        h = nbb.size_histogram(cnt)
        mem = mem_gb()[0]
        qt, ok, infeasible, total = [], 0, 0, 0
        for (s, t) in pairs:
            for frac in (1.0, 0.5):
                b0 = frac * M
                q += 1
                t1 = time.perf_counter()
                got = nbb.battery_query(rank[s], rank[t], b0, first, head, parent,
                                        cnt, I1, C1, O1, pst, PI, PC, PO, soc, dsoc, mark, q)
                qt.append(time.perf_counter() - t1)
                total += 1
                if not verify:
                    infeasible += got == -np.inf
                    continue
                want = nbb.battery_reference(n, off, tgt, w, M, s, b0)[t]
                if want == -np.inf:
                    infeasible += 1
                    ok += got == -np.inf
                else:
                    ok += abs(got - want) <= 1e-6 * max(1.0, abs(want))
        nd = cnt.shape[0]
        share = [f"{100 * h[i] / nd:.3f}%" for i in range(7)]
        lines.append(f"| {factor} D | {t_cust:.1f} s | {PI.shape[0]:,} | {mem:.1f} | " + " | ".join(share)
                     + f" | {h[7]} | {1000 * float(np.median(qt)):.2f} ms | {infeasible}/{total} | "+ (f"{ok}/{total} |" if verify else "not checked (reference too slow at this size) |"))
        log(f"{name} M={factor}D: customize {t_cust:.1f}s, sizes {list(h[:7])}, largest {h[7]}, "
            f"memory {mem:.1f} GB, query {1000 * float(np.median(qt)):.2f} ms, "
            f"correct {ok}/{total}" + ("" if verify else " (verification skipped)") + f", infeasible {infeasible}")
        out.write_text("\n".join(lines), encoding="utf-8")
        del cnt, I1, C1, O1, pst, PI, PC, PO
    lines += ["", "Size columns give the share of shortcut directions with that many pieces."]
    out.write_text("\n".join(lines), encoding="utf-8")
    log(f"peak memory {mem_gb()[1]:.1f} GB; wrote {out.name}")


if __name__ == "__main__":
    main()
