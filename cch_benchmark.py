"""CCH-with-negative-weights benchmark on DIMACS NY and BAY.

    python cch_benchmark.py [NY BAY ...]

Writes results_cch.md (rewritten after each graph).
"""
import random
import statistics
import sys
import time
from pathlib import Path

import dimacs
from cch import CCH
from heuristics import ALTPotential, JohnsonPotential
from preprocess import ALTTables, NegativeCycleError, johnson_potential
from search import astar, bidirectional_astar

HERE = Path(__file__).parent
EPS = 1e-6
OUT = None
N_QUERIES = 100
N_CHECK = 30           # queries also run with the (slow) reference methods
LANDMARKS = 16


def timed(fn):
    t0 = time.perf_counter()
    r = fn()
    return r, time.perf_counter() - t0


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def edge_refs(g):
    return [(u, i) for u in range(g.n) for i in range(len(g.adj[u]))]


def bench(name, lines):
    rng = random.Random(11)
    g_ev = dimacs.load(name, "ev", seed=1)
    g_sh = dimacs.load(name, "shifted", seed=1)
    log(f"{name}: n={g_ev.n:,} m={g_ev.m:,}")
    c, t_build = timed(lambda: CCH(g_ev, leaf=64, log=log))
    lines += [f"## {name} (n={g_ev.n:,}, m={g_ev.m:,})", "",
              f"CCH build (metric-independent, done once for both weightings): "
              f"**{t_build:.0f} s** (order {c.t_order:.0f} s, contraction {c.t_contract:.0f} s, "
              f"triangles {c.t_triangles:.0f} s). {c.m:,} arcs, {len(c.tri_a):,} triangles, "
              f"elimination-tree depth {c.depth}.", ""]
    pairs = [(rng.randrange(g_ev.n), rng.randrange(g_ev.n)) for _ in range(N_QUERIES)]

    for wname, g in (("EV terrain", g_ev), ("random shift", g_sh)):
        neg, t_cust = timed(lambda: c.customize(g))
        assert not neg
        p, t_bf = timed(lambda: johnson_potential(g))
        tables, t_alt = timed(lambda: ALTTables(g, p, k=LANDMARKS, seed=0))
        J, A = JohnsonPotential(p), ALTPotential(tables)
        log(f"{name} {wname}: customize {t_cust:.1f}s, Johnson {t_bf:.1f}s, ALT {t_alt:.1f}s")

        rows = {}
        truth = {}
        for label, fn, npairs in (
                ("Johnson + Dijkstra", lambda s, t: astar(g, s, t, J.to_target(t)), N_CHECK),
                (f"Johnson + ALT ({LANDMARKS})", lambda s, t: astar(g, s, t, A.to_target(t)), N_CHECK),
                ("Bidirectional ALT, 1 core",
                 lambda s, t: bidirectional_astar(g, s, t, A.to_target(t), A.from_source(s)), N_CHECK),
                ("**CCH, no potential**", lambda s, t: c.query(s, t), N_QUERIES)):
            times, ok = [], 0
            for s, t in pairs[:npairs]:
                (d, _), dt = timed(lambda: fn(s, t))
                times.append(dt)
                if label.startswith("Johnson + D"):
                    truth[(s, t)] = d
                elif (s, t) in truth:
                    ok += abs(d - truth[(s, t)]) < EPS
            if label.startswith("Johnson + D"):
                ok = len(times)
            rows[label] = (statistics.mean(times) * 1000, f"{ok}/{min(npairs, N_CHECK)}")
        # all CCH answers vs Dijkstra for the remaining pairs too
        extra = 0
        for s, t in pairs[N_CHECK:]:
            d, _ = astar(g, s, t, J.to_target(t))
            extra += abs(c.query(s, t)[0] - d) < EPS
        cch_row = rows["**CCH, no potential**"]
        rows["**CCH, no potential**"] = (cch_row[0], f"{int(cch_row[1].split('/')[0]) + extra}/{N_QUERIES}")

        lines += [f"### {wname} ({g.negative_fraction():.0%} negative edges)", "",
                  "| Preprocessing for this metric | Time |", "|---|---:|",
                  f"| Johnson's Bellman-Ford | {t_bf:.1f} s |",
                  f"| + {LANDMARKS} ALT landmarks | {t_alt:.1f} s |",
                  f"| **CCH customization (includes negative-cycle check)** | **{t_cust:.1f} s** |", "",
                  "| Query method | Mean query (ms) | Correct |", "|---|---:|---:|"]
        base = rows["Johnson + Dijkstra"][0]
        for label, (ms, ok) in rows.items():
            lines.append(f"| {label} | {ms:.2f} ({base / ms:.1f}x) | {ok} |")
        lines.append("")

        # ---- traffic updates: raise k random edges, partial re-customization
        refs = edge_refs(g)
        if not hasattr(c, "by_c"):
            _, t_idx = timed(c.build_update_index)
            lines += [f"Update index (one-time, metric-independent): {t_idx:.1f} s", ""]
        lines += ["| Edges changed | Partial update | Arcs re-evaluated | Full customization | "
                  "Redo Johnson + ALT | Correct after update |",
                  "|---:|---:|---:|---:|---:|---:|"]
        for k in (1, 10, 100, 1000):
            picked = rng.sample(range(len(refs)), k)
            old = [g.adj[refs[e][0]][refs[e][1]][1] for e in picked]
            new = [w + abs(w) * rng.uniform(0.2, 1.0) + 1 for w in old]
            for e, w in zip(picked, new):
                u, i = refs[e]
                g.adj[u][i] = (g.adj[u][i][0], w)
            evaluated, t_part = timed(lambda: c.update_edges(list(zip(picked, new))))
            p2, t_bf2 = timed(lambda: johnson_potential(g))
            _, t_alt2 = timed(lambda: ALTTables(g, p2, k=LANDMARKS, seed=0))
            J2 = JohnsonPotential(p2)
            ok = 0
            for s, t in pairs[:10]:
                ok += abs(c.query(s, t)[0] - astar(g, s, t, J2.to_target(t))[0]) < EPS
            full = CCH.__new__(CCH)
            full.__dict__.update(c.__dict__)
            _, t_full = timed(lambda: full.customize(g))
            same = all(a == b for a, b in zip(full.up, c.up)) and all(
                a == b for a, b in zip(full.down, c.down))
            lines.append(f"| {k} | {t_part * 1000:.1f} ms | {evaluated:,} | {t_full:.1f} s | "
                         f"{t_bf2 + t_alt2:.1f} s | {ok}/10, arcs identical: {same} |")
            log(f"  update k={k}: partial {t_part * 1000:.1f}ms ({evaluated} arcs), full {t_full:.1f}s")
            # restore
            for e, w in zip(picked, old):
                u, i = refs[e]
                g.adj[u][i] = (g.adj[u][i][0], w)
            c.update_edges(list(zip(picked, old)))
        lines.append("")

        # ---- negative cycle detection
        u = next(u for u in range(g.n) if any(v2 != u and any(x == u for x, _ in g.adj[v2])
                                               for v2, _ in g.adj[u]))
        i = next(i for i, (v2, _) in enumerate(g.adj[u]) if any(x == u for x, _ in g.adj[v2]))
        v2, w_old = g.adj[u][i]
        e_idx = refs.index((u, i))
        g.adj[u][i] = (v2, -1e9)
        found_cch, t_neg = timed(lambda: (c.update_edges([(e_idx, -1e9)]), c.has_negative_cycle())[1])
        t0 = time.perf_counter()
        try:
            johnson_potential(g)
            found_bf = False
        except NegativeCycleError:
            found_bf = True
        t_negbf = time.perf_counter() - t0
        g.adj[u][i] = (v2, w_old)
        c.update_edges([(e_idx, w_old)])
        lines += [f"Negative cycle injected: CCH detects it: **{found_cch}** "
                  f"(partial update + scan {t_neg:.2f} s); Bellman-Ford detects it: {found_bf} "
                  f"({t_negbf:.1f} s). After restoring: CCH reports a cycle: {c.has_negative_cycle()}.", ""]
        OUT.write_text("\n".join(lines), encoding="utf-8")


def main():
    global OUT
    names = sys.argv[1:] or ["NY", "BAY"]
    OUT = HERE / "results" / ("results_cch.md" if names == ["NY", "BAY"] else f"results_cch_{'_'.join(names)}.md")
    lines = ["# CCH with negative weights and no potential", "",
             "All CCH answers are checked against Johnson + Dijkstra. "
             "Pure Python; times are single-threaded.", ""]
    for name in names:
        bench(name, lines)
        OUT.write_text("\n".join(lines), encoding="utf-8")
    log("done")


if __name__ == "__main__":
    main()
