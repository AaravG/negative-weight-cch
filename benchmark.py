"""Benchmark all methods on graphs with negative edges.

Run:  python benchmark.py            (full run)
      python benchmark.py --quick    (small graphs, a few seconds)

Writes results.md and results.csv next to this file.
"""
import csv
import random
import statistics
import sys
import time
from pathlib import Path

import graphs
from heuristics import (ALTPotential, DomainPotential, JohnsonPotential,
                        MaxPotential, naive_euclid)
from parallel import ParallelBidirectional
from preprocess import ALTTables, johnson_potential
from search import astar, astar_reopening, bellman_ford_query, bidirectional_astar

EPS = 1e-6
HERE = Path(__file__).parent


def timed(fn):
    t0 = time.perf_counter()
    out = fn()
    return out, time.perf_counter() - t0


def bench_graph(label, g, queries, landmarks):
    rng = random.Random(42)
    pairs = [(rng.randrange(g.n), rng.randrange(g.n)) for _ in range(queries)]
    print(f"\n== {label}: n={g.n:,} m={g.m:,} negative edges={g.negative_fraction():.0%}")

    truth, bf_times, bf_exp = [], [], []
    for s, t in pairs:
        (d, e), dt = timed(lambda: bellman_ford_query(g, s, t))
        truth.append(d)
        bf_times.append(dt)
        bf_exp.append(e)
    bf_ms = statistics.mean(bf_times) * 1000

    p, t_bf = timed(lambda: johnson_potential(g))
    tables, t_alt = timed(lambda: ALTTables(g, p, k=landmarks, seed=0))
    johnson, alt = JohnsonPotential(p), ALTPotential(tables)
    geometric = g.coords is not None
    dom = DomainPotential(g) if geometric else None
    best = MaxPotential(alt, dom) if geometric else alt

    # (name, preprocessing seconds, query function)
    methods = [("Bellman-Ford (SPFA), no preprocessing", 0.0, None)]
    if geometric:
        h = naive_euclid(g)
        methods.append(("Normal A* (straight-line h, reopening)", 0.0,
                        lambda s, t: astar_reopening(g, s, t, h(t))))
    methods += [
        ("Johnson + Dijkstra", t_bf, lambda s, t: astar(g, s, t, johnson.to_target(t))),
        (f"Johnson + ALT A* ({landmarks} landmarks)", t_bf + t_alt,
         lambda s, t: astar(g, s, t, alt.to_target(t))),
    ]
    if geometric:
        methods += [
            ("Domain-bound A* (no preprocessing)", 0.0,
             lambda s, t: astar(g, s, t, dom.to_target(t))),
            ("Combined max(ALT, domain) A*", t_bf + t_alt,
             lambda s, t: astar(g, s, t, best.to_target(t))),
        ]
    methods.append((f"Bidirectional {best.name}, 1 core", t_bf + t_alt,
                    lambda s, t: bidirectional_astar(g, s, t, best.to_target(t),
                                                     best.from_source(s))))

    rows = []
    for name, pre, fn in methods:
        if fn is None:
            times, exps, dists = bf_times, bf_exp, truth
        else:
            times, exps, dists = [], [], []
            for s, t in pairs:
                (d, e), dt = timed(lambda: fn(s, t))
                times.append(dt)
                exps.append(e)
                dists.append(d)
        rows.append(summarize(label, g, name, pre, times, exps, dists, truth, bf_ms))

    par, t_spawn = timed(lambda: ParallelBidirectional(g, best))
    try:
        par.query(*pairs[0])  # warm-up: workers finish unpickling the graph
        times, exps, dists = [], [], []
        for s, t in pairs:
            (d, e), dt = timed(lambda: par.query(s, t))
            times.append(dt)
            exps.append(e)
            dists.append(d)
    finally:
        par.close()
    rows.append(summarize(label, g, f"Bidirectional {best.name}, 2 cores", t_bf + t_alt,
                          times, exps, dists, truth, bf_ms))
    print(f"   (starting the 2 worker processes took {t_spawn * 1000:.0f} ms, not counted)")
    return rows


def summarize(label, g, name, pre, times, exps, dists, truth, bf_ms):
    ok = sum(abs(a - b) < EPS for a, b in zip(dists, truth))
    q_ms = statistics.mean(times) * 1000
    saved = bf_ms - q_ms
    if pre == 0:
        breakeven = "-"
    elif saved > 0:
        breakeven = f"{pre * 1000 / saved:.1f}"
    else:
        breakeven = "never"
    row = {
        "graph": label, "n": g.n, "method": name,
        "preprocess_s": round(pre, 3),
        "query_ms": round(q_ms, 2),
        "expanded": round(statistics.mean(exps)),
        "expanded_pct": round(100 * statistics.mean(exps) / g.n, 1),
        "speedup_vs_BF": round(bf_ms / q_ms, 1) if q_ms > 0 else 0,
        "correct": f"{ok}/{len(truth)}",
        "breakeven_queries": breakeven,
    }
    print(f"   {name:45s} pre={row['preprocess_s']:7.2f}s  query={row['query_ms']:9.2f}ms  "
          f"expanded={row['expanded']:>7,} ({row['expanded_pct']:5.1f}%)  "
          f"x{row['speedup_vs_BF']:<6} correct={row['correct']}  break-even={breakeven}")
    return row


def write_reports(rows, out="results"):
    with open(HERE / f"{out}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    lines = ["# Benchmark results", ""]
    for label in dict.fromkeys(r["graph"] for r in rows):
        lines += [f"## {label}", "",
                  "| Method | Preprocess (s) | Query (ms) | Nodes expanded | Speedup vs BF | Correct | Break-even queries |",
                  "|---|---:|---:|---:|---:|---:|---:|"]
        for r in rows:
            if r["graph"] == label:
                lines.append(f"| {r['method']} | {r['preprocess_s']} | {r['query_ms']} | "
                             f"{r['expanded']:,} ({r['expanded_pct']}%) | {r['speedup_vs_BF']}x | "
                             f"{r['correct']} | {r['breakeven_queries']} |")
        lines.append("")
    (HERE / f"{out}.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    quick = "--quick" in sys.argv
    out = "results"
    if "--dimacs" in sys.argv:
        import dimacs
        suite = []
        for name in ("NY", "BAY"):
            suite += [(f"DIMACS {name}, EV terrain",
                       lambda name=name: dimacs.load(name, "ev", seed=1)),
                      (f"DIMACS {name}, random potential shift",
                       lambda name=name: dimacs.load(name, "shifted", seed=1))]
        queries, landmarks, out = 20, 16, "results_dimacs"
    elif quick:
        suite = [("EV grid 40x40", lambda: graphs.ev_grid(40, seed=1)),
                 ("Shifted random 2k", lambda: graphs.shifted_random(2000, seed=1))]
        queries, landmarks = 10, 8
    else:
        suite = [("EV grid 100x100", lambda: graphs.ev_grid(100, seed=1)),
                 ("EV grid 200x200", lambda: graphs.ev_grid(200, seed=2)),
                 ("EV road-like 30k", lambda: graphs.ev_geometric(30000, seed=3)),
                 ("Shifted random 20k (no geometry)", lambda: graphs.shifted_random(20000, seed=4))]
        queries, landmarks = 30, 16
    rows = []
    for label, make in suite:
        rows += bench_graph(label, make(), queries, landmarks)
    write_reports(rows, out)
    print(f"\nwrote {out}.md and {out}.csv")


if __name__ == "__main__":
    main()
