"""Correctness checks on many small graphs: every exact method must match Bellman-Ford,
and every potential must be feasible (and a lower bound where it claims to be).

Run:  python test_correctness.py
"""
import random

import graphs
from heuristics import ALTPotential, DomainPotential, JohnsonPotential, MaxPotential
from parallel import ParallelBidirectional
from preprocess import ALTTables, johnson_potential
from search import astar, bellman_ford_query, bidirectional_astar

EPS = 1e-6


def check_potential(g, pot, s, t, lower_bound):
    pt, ps = pot.to_target(t), pot.from_source(s)
    for u in range(g.n):
        for v, w in g.adj[u]:
            assert w + pt(v) - pt(u) >= -EPS, f"{pot.name}: to_target inconsistent"
            assert w + ps(u) - ps(v) >= -EPS, f"{pot.name}: from_source inconsistent"
    if lower_bound:
        for v in random.sample(range(g.n), min(20, g.n)):
            dvt, _ = bellman_ford_query(g, v, t)
            assert pt(v) <= dvt + EPS, f"{pot.name}: not a lower bound"


def run():
    rng = random.Random(1)
    cases = ([("grid", graphs.ev_grid(12, seed=i)) for i in range(4)]
             + [("geo", graphs.ev_geometric(200, seed=i)) for i in range(4)]
             + [("shift", graphs.shifted_random(150, seed=i)) for i in range(4)])
    checked = 0
    for kind, g in cases:
        p = johnson_potential(g)
        alt = ALTPotential(ALTTables(g, p, k=4, seed=0))
        pots = [(JohnsonPotential(p), False), (alt, True)]
        if g.coords is not None:
            dom = DomainPotential(g)
            pots += [(dom, True), (MaxPotential(alt, dom), True)]
        s0, t0 = rng.randrange(g.n), rng.randrange(g.n)
        for pot, lb in pots:
            check_potential(g, pot, s0, t0, lb)
        for _ in range(15):
            s, t = rng.randrange(g.n), rng.randrange(g.n)
            truth, _ = bellman_ford_query(g, s, t)
            for pot, _lb in pots:
                d1, _ = astar(g, s, t, pot.to_target(t))
                d2, _ = bidirectional_astar(g, s, t, pot.to_target(t), pot.from_source(s))
                assert abs(d1 - truth) < EPS, (kind, pot.name, "astar", d1, truth)
                assert abs(d2 - truth) < EPS, (kind, pot.name, "bidir", d2, truth)
                checked += 2
        par = ParallelBidirectional(g, pots[-1][0])
        try:
            for _ in range(15):
                s, t = rng.randrange(g.n), rng.randrange(g.n)
                truth, _ = bellman_ford_query(g, s, t)
                d3, _ = par.query(s, t)
                assert abs(d3 - truth) < EPS, (kind, "parallel", d3, truth)
                checked += 1
        finally:
            par.close()
        print(f"  ok  {kind:5s} n={g.n:4d}  negative edges={g.negative_fraction():.0%}")
    print(f"all {checked} query checks passed")


if __name__ == "__main__":
    run()
