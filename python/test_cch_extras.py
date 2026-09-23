"""Further CCH techniques under negative weights: path unpacking and the
tie-based partial update of Dibbelt et al. (2016), section 7.7.

(Parallel customization is tested in the C++ implementation: `cch --parallel`
compares its output against the serial one.)

Run:  python test_cch_extras.py
"""
import random

import graphs
from cch import CCH
from search import bellman_ford_query

EPS = 1e-6


def close(a, b):
    return a == b or abs(a - b) <= EPS * max(1.0, abs(b))


def cases(rng):
    return ([("grid", graphs.ev_grid(12, seed=i)) for i in range(3)]
            + [("geo", graphs.ev_geometric(200, seed=i)) for i in range(2)]
            + [("shift", graphs.shifted_random(150, seed=i)) for i in range(2)])


def test_unpacking(rng):
    tested = good = 0
    for kind, g in cases(rng):
        xy = g.coords or [(rng.random(), rng.random()) for _ in range(g.n)]
        c = CCH(g, xy=xy, leaf=8, log=lambda *_: None)
        c.customize(g)
        for _ in range(20):
            s, t = rng.randrange(g.n), rng.randrange(g.n)
            truth, _ = bellman_ford_query(g, s, t)
            if truth == float("inf"):
                continue
            d, path = c.query_path(s, t)
            tested += 1
            total, ok = 0.0, path[0] == s and path[-1] == t
            for u, v in zip(path, path[1:]):
                w = min([wt for x, wt in g.adj[u] if x == v], default=None)
                if w is None:
                    ok = False
                    break
                total += w
            good += ok and close(total, truth) and close(d, truth)
    print(f"  path unpacking: {good}/{tested} paths are real edge sequences of exactly "
          f"the shortest length")
    return good == tested


def test_tie_updates(rng):
    matched = total = 0
    for kind, g in cases(rng):
        xy = g.coords or [(rng.random(), rng.random()) for _ in range(g.n)]
        c = CCH(g, xy=xy, leaf=8, log=lambda *_: None)
        c.customize(g)
        refs = [(u, i) for u in range(g.n) for i in range(len(g.adj[u]))]
        for trial in range(6):
            picked = rng.sample(range(len(refs)), 8)
            old = [g.adj[refs[k][0]][refs[k][1]][1] for k in picked]
            grow = trial % 2 == 0          # alternate increases and decreases
            new = [w * (1.5 if grow else 0.4) + (5 if grow else -5) for w in old]
            for k, w in zip(picked, new):
                u, i = refs[k]
                g.adj[u][i] = (g.adj[u][i][0], w)
            c.update_edges_tiebased(list(zip(picked, new)))
            tie_up, tie_dn = list(c.up), list(c.down)
            full = CCH.__new__(CCH)
            full.__dict__.update(c.__dict__)
            full.customize(g)
            bad = sum(1 for a, b in zip(tie_up, full.up) if not close(a, b))
            bad += sum(1 for a, b in zip(tie_dn, full.down) if not close(a, b))
            total += 1
            matched += bad == 0
            c.up, c.down = full.up, full.down
            c.base_up, c.base_down, c.edge_w = full.base_up, full.base_down, full.edge_w
    print(f"  tie-based partial updates: {matched}/{total} runs matched full customization")
    return matched == total


def run():
    rng = random.Random(11)
    ok = test_unpacking(rng)
    ok &= test_tie_updates(random.Random(17))
    print("all extra technique checks passed" if ok else "FAILURES (see above)")


if __name__ == "__main__":
    run()
