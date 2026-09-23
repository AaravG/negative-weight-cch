"""Potential-free CCH (Numba) on a DIMACS graph or region, or the full USA graph.

    python run_nb.py NY                 # DIMACS graph via dimacs.py
    python run_nb.py FLA_R              # region cut out of the USA graph
    python run_nb.py USA                # full USA, from the memory-mapped cache
                                        # built by bigcsr.py (data/cache/USA)

For each weighting (EV terrain, random shift): customization time, a
negative-cycle scan, and N random queries checked against Dijkstra with a
Johnson potential. Writes results/results_nb_<name>.md.
"""
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

import cch_nb as nb
from cch_memory import mem_gb

HERE = Path(__file__).parent
QUERIES = 100


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load(name):
    """Returns n, off, tgt, x, y and {weighting: (w, p or None)}."""
    if name == "USA":
        cache = Path(os.environ.get("USA_CACHE", HERE / "data" / "cache" / "USA"))
        topo = json.loads((cache / "topology.json").read_text())
        n = topo["n"]
        off = np.fromfile(cache / "off.bin", np.int32).astype(np.int64)
        tgt = np.fromfile(cache / "tgt.bin", np.int32)
        x = np.fromfile(cache / "xm.bin", np.float64)
        y = np.fromfile(cache / "ym.bin", np.float64)
        ws = {}
        for w in ("ev", "shifted"):
            fw = np.fromfile(cache / f"{w}_fw.bin", np.float64)
            pf = cache / f"{w}_p.bin"
            p = np.fromfile(pf, np.float64) if pf.exists() else None
            ws[w] = (fw, p)
        return n, off, tgt, x, y, ws
    import dimacs
    ws = {}
    base = None
    for w in ("ev", "shifted", "ev_real"):
        g = dimacs.load(name, w, seed=1)
        if base is None:
            n = g.n
            off = np.zeros(n + 1, np.int64)
            for u in range(n):
                off[u + 1] = off[u] + len(g.adj[u])
            tgt = np.fromiter((v for u in range(n) for v, _ in g.adj[u]), np.int32, off[n])
            base = (n, off, tgt)
            gx = dimacs.load(name, "ev", seed=1) if g.coords is None else g
            x = np.array([c[0] for c in gx.coords])
            y = np.array([c[1] for c in gx.coords])
        wt = np.fromiter((c for u in range(n) for _, c in g.adj[u]), np.float64, off[n])
        ws[w] = (wt, None)
    return n, off, tgt, x, y, ws


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "NY"
    out = HERE / "results" / f"results_nb_{name}.md"
    lines = [f"# Potential-free CCH (Numba) on {name}", ""]

    t0 = time.perf_counter()
    n, off, tgt, x, y, ws = load(name)
    log(f"{name}: n={n:,} arcs={off[n]:,} loaded in {time.perf_counter() - t0:.0f}s, "
        f"memory {mem_gb()[0]:.1f} GB")

    t0 = time.perf_counter()
    noff, nbrs, rev = nb.undirected_csr(n, off, tgt)
    t_und = time.perf_counter() - t0
    log(f"undirected graph: {nbrs.shape[0]:,} slots in {t_und:.0f}s")

    rank_file = HERE / "data" / f"rank_{name}.npy"
    if rank_file.exists():
        rank = np.load(rank_file)
        t_order = float("nan")
        log("ordering loaded from cache")
    else:
        t0 = time.perf_counter()
        rank = nb.nested_dissection(n, noff, nbrs, rev, x, y, 64)
        t_order = time.perf_counter() - t0
        np.save(rank_file, rank)
        log(f"ordering: {t_order:.0f}s")
    assert np.array_equal(np.sort(rank), np.arange(n)), "rank is not a permutation"
    del rev

    t0 = time.perf_counter()
    first, head, parent = nb.symbolic(n, noff, nbrs, rank)
    t_sym = time.perf_counter() - t0
    del noff, nbrs
    depth = np.zeros(n, np.int64)
    for v in range(n - 1, -1, -1):
        if parent[v] >= 0:
            depth[v] = depth[parent[v]] + 1
    d = int(depth.max())
    del depth
    tri = nb.count_triangles(n, first, head)
    t0 = time.perf_counter()
    e_arc, e_up = nb.edge_arcs(n, off, tgt, rank, first, head)
    t_map = time.perf_counter() - t0
    log(f"shortcut arcs {head.shape[0]:,}, triangles {tri:,}, depth {d}; "
        f"symbolic {t_sym:.0f}s, mapping {t_map:.0f}s, memory {mem_gb()[0]:.1f} GB")
    lines += [f"n = {n:,}, original arcs = {off[n]:,}. Ordering: {t_order:.0f} s. "
              f"Symbolic contraction: {t_sym:.0f} s. Shortcut arcs: {head.shape[0]:,}. "
              f"Triangles: {tri:,}. Elimination-tree depth: {d}.", ""]

    rng = np.random.default_rng(7)
    pairs = [(int(rng.integers(n)), int(rng.integers(n))) for _ in range(QUERIES)]
    df = np.full(n, np.inf)
    dr = np.full(n, np.inf)
    mark = np.zeros(n, np.int64)
    for wname, (w, p) in ws.items():
        t0 = time.perf_counter()
        upw, dnw = nb.customize(n, first, head, e_arc, e_up, w)
        t_cust = time.perf_counter() - t0
        t0 = time.perf_counter()
        neg = nb.negative_cycle_arcs(upw, dnw)
        t_neg = time.perf_counter() - t0
        if p is None:
            t0 = time.perf_counter()
            p = nb.spfa_potential(n, off, tgt, w)
            log(f"{wname}: Johnson potential (reference only) {time.perf_counter() - t0:.0f}s")
        q_times, r_times, ok = [], [], 0
        for i, (s, t) in enumerate(pairs):
            t1 = time.perf_counter()
            got = nb.query(rank[s], rank[t], first, head, parent, upw, dnw, df, dr, mark, i + 1)
            q_times.append(time.perf_counter() - t1)
            t1 = time.perf_counter()
            want = nb.dijkstra_potential(n, off, tgt, w, p, s, t)
            r_times.append(time.perf_counter() - t1)
            ok += (got == want) or abs(got - want) <= 1e-6 * max(1.0, abs(want))
        # first call includes JIT compilation; report the median
        qm = float(np.median(q_times)) * 1000
        rm = float(np.median(r_times)) * 1000
        log(f"{wname}: customize {t_cust:.1f}s, negative-cycle arcs {neg} ({t_neg:.2f}s), "
            f"query {qm:.2f} ms vs Dijkstra {rm:.0f} ms, correct {ok}/{QUERIES}")
        lines += [f"## {wname}", "",
                  "| Measure | Value |", "|---|---:|",
                  f"| Customization (compiled) | {t_cust:.1f} s |",
                  f"| Negative-cycle scan | {t_neg:.2f} s ({neg} cycle arcs) |",
                  f"| Median CCH query | {qm:.2f} ms |",
                  f"| Median Johnson + Dijkstra query (compiled, reference) | {rm:.0f} ms |",
                  f"| Correct | {ok}/{QUERIES} |", ""]
        out.write_text("\n".join(lines), encoding="utf-8")
    log(f"peak memory {mem_gb()[1]:.1f} GB; wrote {out.name}")


if __name__ == "__main__":
    main()
