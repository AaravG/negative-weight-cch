"""Correctness tests for CCH with negative weights (no potentials anywhere).

  * queries match Bellman-Ford on graphs with negative edges
  * negative-cycle detection agrees with Bellman-Ford
  * partial re-customization gives exactly the same arc values as a full one

Run:  python test_cch.py
"""
import random

import graphs
from cch import CCH
from preprocess import NegativeCycleError, johnson_potential
from search import bellman_ford_query

EPS = 1e-6


def edge_list(g):
    return [(u, i) for u in range(g.n) for i in range(len(g.adj[u]))]


def set_weight(g, u, i, w):
    v, _ = g.adj[u][i]
    g.adj[u][i] = (v, w)


def run():
    rng = random.Random(3)
    cases = ([("grid", graphs.ev_grid(15, seed=i)) for i in range(4)]
             + [("geo", graphs.ev_geometric(300, seed=i)) for i in range(4)]
             + [("shift", graphs.shifted_random(200, seed=i)) for i in range(4)])
    checks = 0
    for kind, g in cases:
        xy = g.coords or [(rng.random(), rng.random()) for _ in range(g.n)]
        c = CCH(g, xy=xy, leaf=8, log=lambda *_: None)
        assert not c.customize(g), "false negative cycle"
        for _ in range(25):
            s, t = rng.randrange(g.n), rng.randrange(g.n)
            truth, _ = bellman_ford_query(g, s, t)
            d, scanned = c.query(s, t)
            assert abs(d - truth) < EPS or d == truth, (kind, s, t, d, truth)
            # upward arcs of ancestors stay among ancestors: <= d(d+1)/2 per side
            assert scanned <= c.depth * (c.depth + 1), (scanned, c.depth)
            checks += 2

        # partial updates: raise random edges (traffic), then restore them
        edges = edge_list(g)
        for _ in range(5):
            picked = rng.sample(range(len(edges)), 10)
            old = [g.adj[edges[k][0]][edges[k][1]][1] for k in picked]
            for rnd, weights in (("raise", [w + rng.uniform(1, 50) for w in old]),
                                 ("restore", old)):
                for k, w in zip(picked, weights):
                    set_weight(g, *edges[k], w)
                c.update_edges(list(zip(picked, weights)))
                assert not c.update_created_negative_cycle()
                up, down = list(c.up), list(c.down)
                full = CCH.__new__(CCH)
                full.__dict__.update(c.__dict__)
                full.customize(g)
                assert all(abs(a - b) < EPS or a == b for a, b in zip(up, full.up)), rnd
                assert all(abs(a - b) < EPS or a == b for a, b in zip(down, full.down)), rnd
                c.up, c.down = full.up, full.down
                c.base_up, c.base_down, c.edge_w = full.base_up, full.base_down, full.edge_w
                s, t = rng.randrange(g.n), rng.randrange(g.n)
                assert abs(c.query(s, t)[0] - bellman_ford_query(g, s, t)[0]) < EPS
                checks += 3

        # negative cycle created by a partial update, detected from changed
        # arcs only, then removed again by restoring the weight
        refs = edge_list(g)
        for k, (u, i) in enumerate(refs):
            v, w_old = g.adj[u][i]
            if any(x == u for x, _ in g.adj[v]):
                break
        set_weight(g, u, i, -1e6)
        c.update_edges([(k, -1e6)])
        assert c.update_created_negative_cycle() and c.has_negative_cycle()
        set_weight(g, u, i, w_old)
        c.update_edges([(k, w_old)])
        assert not c.has_negative_cycle()
        checks += 2

        # negative cycle: make one edge and its reverse very negative
        u = rng.randrange(g.n)
        if g.adj[u]:
            v, w_uv = g.adj[u][0]
            back = [i for i, (x, _) in enumerate(g.adj[v]) if x == u]
            set_weight(g, u, 0, -1e6)
            if not back:
                g.adj[v].append((u, 0.0))
                g.radj[u].append((v, 0.0))
                c2 = CCH(g, xy=xy, leaf=8, log=lambda *_: None)
            else:
                c2 = c
            found = bool(c2.customize(g))
            try:
                johnson_potential(g)
                bf_found = False
            except NegativeCycleError:
                bf_found = True
            assert found == bf_found == True, (kind, found, bf_found)
            checks += 1
        print(f"  ok  {kind:5s} n={g.n:4d} arcs={c.m:5d} depth={c.depth}")
    print(f"all {checks} CCH checks passed")


if __name__ == "__main__":
    run()
