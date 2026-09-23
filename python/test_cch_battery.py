"""Battery-constrained CCH vs. a label-correcting reference.

  * composition/dominance agree with direct simulation along random paths
  * CCH queries equal the reference max-SoC for many (s, t, b0, M)

Run:  python test_cch_battery.py
"""
import random

import graphs
from cch import CCH
from cch_battery import (NEG, BatteryCCH, arc_function, battery_label_correcting,
                         compose, evaluate)

EPS = 1e-6


def simulate(costs, b, M):
    for c in costs:
        if c >= 0:
            if b < c - 1e-12:
                return NEG
            b -= c
        else:
            b = min(M, b - c)
    return b


def test_algebra(rng):
    n = 0
    for _ in range(3000):
        M = rng.uniform(1, 20)
        costs = [rng.uniform(-6, 6) for _ in range(rng.randint(1, 8))]
        f = arc_function(costs[0], M)
        for c in costs[1:]:
            if f is None:
                break
            g = arc_function(c, M)
            f = compose(f, g, M) if g is not None else None
        for _ in range(5):
            b = rng.uniform(0, M)
            want = simulate(costs, b, M)
            got = NEG if f is None else evaluate(f, b)
            assert abs(got - want) < EPS or got == want, (costs, M, b, got, want)
            n += 1
    return n


def run():
    rng = random.Random(9)
    checks = test_algebra(rng)
    print(f"  algebra: {checks} simulated evaluations agree")
    cases = ([("grid", graphs.ev_grid(12, seed=i)) for i in range(3)]
             + [("geo", graphs.ev_geometric(250, seed=i)) for i in range(3)])
    for kind, g in cases:
        c = CCH(g, leaf=8, log=lambda *_: None)
        bc = BatteryCCH(c)
        for M in (3.0, 10.0, 40.0):
            biggest, avg = bc.customize(g, M)
            unreachable = 0
            for _ in range(30):
                s, t = rng.randrange(g.n), rng.randrange(g.n)
                b0 = rng.uniform(0, M)
                want = battery_label_correcting(g, s, b0, M)[t]
                got = bc.query(s, t, b0)
                assert (got == want == NEG) or abs(got - want) < EPS, (kind, M, s, t, b0, got, want)
                unreachable += want == NEG
                checks += 1
            print(f"  ok  {kind:4s} n={g.n} M={M:5.1f}: profile size max {biggest}, "
                  f"mean {avg:.2f}; {unreachable}/30 queries infeasible")
    print(f"all {checks} battery checks passed")


if __name__ == "__main__":
    run()
