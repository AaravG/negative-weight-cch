"""Which CCH acceleration techniques survive negative weights?

Checks, on small graphs with negative edges:
  1. perfect customization  -> does every arc end up holding the true distance?
  2. witness pruning        -> are queries still exact after removing the arcs a
                               witness search would remove?
  3. stall-on-demand        -> counterexample search (expected to fail)

Run:  python test_cch_pruning.py
"""
import random

import graphs
from cch import CCH
from search import bellman_ford_query

EPS = 1e-6


def true_distances(g, nodes):
    return {(s, t): bellman_ford_query(g, s, t)[0] for s in nodes for t in nodes if s != t}


def check_graph(kind, g, rng, report):
    xy = g.coords or [(rng.random(), rng.random()) for _ in range(g.n)]
    c = CCH(g, xy=xy, leaf=8, log=lambda *_: None)
    c.customize(g)

    # ---- 1. perfect customization
    c.perfect_customize()
    bad_arc = 0
    checked = 0
    inv = [0] * c.n
    for v in range(c.n):
        inv[c.rank[v]] = v
    for x in range(c.n):
        for a in range(c.first[x], c.first[x + 1]):
            y = c.head[a]
            u, v = inv[x], inv[y]
            duv, _ = bellman_ford_query(g, u, v)
            dvu, _ = bellman_ford_query(g, v, u)
            checked += 2
            if not (abs(c.up[a] - duv) < EPS or (c.up[a] == duv)):
                bad_arc += 1
            if not (abs(c.down[a] - dvu) < EPS or (c.down[a] == dvu)):
                bad_arc += 1
            if checked > 400:
                break
        if checked > 400:
            break
    report["perfect_checked"] += checked
    report["perfect_bad"] += bad_arc

    # ---- 2. witness pruning (arcs a witness search would drop)
    drop_up, drop_dn = c.prunable_arcs()
    dropped = sum(drop_up) + sum(drop_dn)
    bad_q = 0
    for _ in range(25):
        s, t = rng.randrange(g.n), rng.randrange(g.n)
        truth, _ = bellman_ford_query(g, s, t)
        got = c.query_pruned(s, t, drop_up, drop_dn)
        if not (abs(got - truth) < EPS or got == truth):
            bad_q += 1
        report["prune_queries"] += 1
    report["prune_bad"] += bad_q
    report["prune_dropped"] += dropped
    report["prune_arcs"] += 2 * c.m

    # ---- 3. stalling: does the stall condition ever prune a needed vertex?
    # (a stalled vertex x is one with a cheaper path from a higher neighbour y)
    stalls_that_matter = 0
    for _ in range(25):
        s, t = rng.randrange(g.n), rng.randrange(g.n)
        truth, _ = bellman_ford_query(g, s, t)
        if truth == float("inf"):
            continue
        # forward labels of the elimination-tree sweep
        c.query(s, t)
        df = c._df
        rs = c.rank[s]
        x = rs
        while x >= 0:
            if df[x] < float("inf"):
                for a in range(c.first[x], c.first[x + 1]):
                    y = c.head[a]
                    if df[y] < float("inf") and df[y] + c.down[a] < df[x] - 1e-9:
                        stalls_that_matter += 1   # x would be stalled although its
                        break                     # label is the exact up-distance
            x = c.parent[x]
    report["stall_hits"] += stalls_that_matter
    print(f"  {kind:5s} n={g.n:4d}: perfect-customization mismatches {bad_arc}/{checked}, "
          f"pruning removes {dropped}/{2 * c.m} arcs, wrong queries after pruning {bad_q}/25, "
          f"stallable vertices seen {stalls_that_matter}")


def run():
    rng = random.Random(5)
    report = dict(perfect_checked=0, perfect_bad=0, prune_queries=0, prune_bad=0,
                  prune_dropped=0, prune_arcs=0, stall_hits=0)
    cases = ([("grid", graphs.ev_grid(12, seed=i)) for i in range(3)]
             + [("geo", graphs.ev_geometric(200, seed=i)) for i in range(3)]
             + [("shift", graphs.shifted_random(150, seed=i)) for i in range(3)])
    for kind, g in cases:
        check_graph(kind, g, rng, report)
    print()
    print(f"perfect customization: {report['perfect_bad']} mismatches in "
          f"{report['perfect_checked']} arc distances")
    print(f"witness pruning: {report['prune_bad']} wrong out of {report['prune_queries']} queries; "
          f"{report['prune_dropped']} of {report['prune_arcs']} arc directions removable")
    print(f"stalling: {report['stall_hits']} vertices would be stalled although their sweep "
          f"label is already exact")


if __name__ == "__main__":
    run()
