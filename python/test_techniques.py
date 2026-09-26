"""Which CCH techniques stay exact with negative weights, and which fail?

Every check compares against Bellman-Ford and is repeated on a NON-NEGATIVE
control with the same topology and order, so that failures caused by negative
weights can be told apart from standard CCH conditions.

  1. the counterexamples of the paper (Fig. 1, zero-cycle witness pruning,
     lower-triangle pruning, saturated-integer infinity), run through CCH;
  2. stall-on-demand in the elimination-tree query     -> exact
  3. Dijkstra-based query / distance pruning / early stop -> fail
  4. perfect customization on ALL arcs                  -> exact
  5. witness pruning (upper/intermediate triangles)     -> exact without zero cycles,
     and a lexicographic tie-break repairs zero cycles
  6. the potential read off the customized CCH          -> feasible, equals Bellman-Ford
  7. (--exhaustive) all conservative 3-vertex digraphs with weights in {-2..2}

Run:  python test_techniques.py [--exhaustive]
"""
import heapq
import itertools
import math
import random
import sys
from collections import deque

import graphs
from cch import CCH, INF

NOLOG = lambda *_: None


# ------------------------------------------------------------------ helpers

def make_graph(n, arcs):
    g = graphs.Graph(n)
    for u, v, w in arcs:
        g.add_edge(u, v, float(w))
    return g


def abs_control(g):
    h = graphs.Graph(g.n)
    for u in range(g.n):
        for v, w in g.adj[u]:
            h.add_edge(u, v, abs(w))
    h.coords = g.coords
    return h


def sssp(g, s):
    """Bellman-Ford (queue-based) from s; the graph must be conservative."""
    d = [INF] * g.n
    d[s] = 0.0
    inq = [False] * g.n
    inq[s] = True
    q = deque([s])
    while q:
        u = q.popleft()
        inq[u] = False
        for v, w in g.adj[u]:
            nd = d[u] + w
            if nd < d[v] - 1e-12:
                d[v] = nd
                if not inq[v]:
                    inq[v] = True
                    q.append(v)
    return d


def close(a, b):
    return a == b or (a != INF and b != INF and abs(a - b) <= 1e-6 * max(1.0, abs(b)))


def build(g, order=None, rng=None):
    if order is not None:
        c = CCH(g, log=NOLOG, order=order)
    else:
        rng = rng or random.Random(0)
        xy = g.coords or [(rng.random(), rng.random()) for _ in range(g.n)]
        c = CCH(g, xy=xy, leaf=8, log=NOLOG)
    assert not c.customize(g)
    return c


def ancestors(c, x):
    out = []
    while x >= 0:
        out.append(x)
        x = c.parent[x]
    return out


# --------------------------------------------------------- query variants

def et_query(c, s, t, prune=False, early_stop=False, stall=False):
    """Elimination-tree query over the union of both ancestor paths in rank order.
    prune:      skip relaxing a vertex whose label is >= mu (Buchhold et al.)
    early_stop: stop once every remaining tentative label is >= mu
    stall:      stall-on-demand (strict test, current labels of higher neighbours)"""
    first, head, up, dn = c.first, c.head, c.up, c.down
    A, B = ancestors(c, c.rank[s]), ancestors(c, c.rank[t])
    inA, inB = set(A), set(B)
    df = {x: INF for x in A}
    dr = {x: INF for x in B}
    df[A[0]] = 0.0
    dr[B[0]] = 0.0
    mu = INF
    stalls = 0
    order = sorted(inA | inB)
    for i, x in enumerate(order):
        if early_stop:
            rest = [df[y] for y in order[i:] if y in inA] + [dr[y] for y in order[i:] if y in inB]
            if all(v >= mu for v in rest):
                break
        if x in inA and x in inB:
            mu = min(mu, df[x] + dr[x])
        for side, lab, w_out, w_in, members in ((0, df, up, dn, inA), (1, dr, dn, up, inB)):
            if x not in members or lab[x] == INF:
                continue
            if prune and lab[x] >= mu:
                continue
            if stall and any(lab.get(head[a], INF) + w_in[a] < lab[x]
                             for a in range(first[x], first[x + 1])):
                stalls += 1
                continue
            for a in range(first[x], first[x + 1]):
                y = head[a]
                v = lab[x] + w_out[a]
                if y in members and v < lab[y]:
                    lab[y] = v
    return mu, stalls


def dijkstra_query(c, s, t, label_setting=True, stop=True):
    """Bidirectional Dijkstra-based CCH query on the (unshifted) customized arcs."""
    first, head, up, dn = c.first, c.head, c.up, c.down
    rs, rt = c.rank[s], c.rank[t]
    d = [{rs: 0.0}, {rt: 0.0}]
    q = [[(0.0, rs)], [(0.0, rt)]]
    done = [set(), set()]
    mu = INF
    while q[0] or q[1]:
        tf = q[0][0][0] if q[0] else INF
        tr = q[1][0][0] if q[1] else INF
        if stop and min(tf, tr) >= mu:
            break
        side = 0 if tf <= tr else 1
        dx, x = heapq.heappop(q[side])
        if dx > d[side].get(x, INF) or (label_setting and x in done[side]):
            continue
        done[side].add(x)
        if x in d[1 - side]:
            mu = min(mu, dx + d[1 - side][x])
        w = up if side == 0 else dn
        for a in range(first[x], first[x + 1]):
            y = head[a]
            nd = dx + w[a]
            if label_setting and y in done[side]:
                continue
            if nd < d[side].get(y, INF):
                d[side][y] = nd
                heapq.heappush(q[side], (nd, y))
    return mu


# ------------------------------------------------------ 1. the counterexamples

def check_examples():
    print("1. counterexamples (order = vertex ids)")
    ex = {
        "Ex. 1 Dijkstra-based query": (3, [(0, 1, 0), (0, 2, -1), (1, 2, -2)], 0, 2,
                                       lambda c, s, t: dijkstra_query(c, s, t)),
        "Ex. 1 label-setting, no stop": (3, [(0, 1, 0), (0, 2, -1), (1, 2, -2)], 0, 2,
                                         lambda c, s, t: dijkstra_query(c, s, t, stop=False)),
        "Ex. 2 distance pruning": (4, [(0, 1, 0), (0, 2, 1), (2, 3, -1), (3, 1, -2)], 0, 1,
                                   lambda c, s, t: et_query(c, s, t, prune=True)[0]),
        "Ex. 3 early termination": (3, [(0, 1, -3), (0, 2, -2), (2, 1, -2)], 0, 1,
                                    lambda c, s, t: et_query(c, s, t, early_stop=True)[0]),
    }
    ok = True
    for name, (n, arcs, s, t, method) in ex.items():
        for label, g in (("negative", make_graph(n, arcs)),
                         ("|w| control", abs_control(make_graph(n, arcs)))):
            c = build(g, order=list(range(n)))
            truth = sssp(g, s)[t]
            got, plain = method(c, s, t), c.query(s, t)[0]
            fails = not close(got, truth)
            expect_fail = label == "negative"
            ok &= fails == expect_fail and close(plain, truth)
            print(f"   {name:30s} {label:12s} correct {truth:5g}, variant {got:5g}, "
                  f"ET query {plain:5g}  {'FAILS' if fails else 'exact'}")

    # zero-length cycle: witness pruning removes both mutually witnessing arcs
    for label, arcs in (("negative", [(0, 1, -2), (1, 2, 0), (2, 0, 2)]),
                        ("all-zero", [(0, 1, 0), (1, 2, 0), (2, 0, 0)])):
        g = make_graph(3, arcs)
        c = build(g, order=[0, 1, 2])
        c.perfect_customize()
        drop_up, drop_dn = c.prunable_arcs()
        got, truth = c.query_pruned(0, 1, drop_up, drop_dn), sssp(g, 0)[1]
        ok &= not close(got, truth)
        print(f"   zero-cycle witness pruning    {label:12s} correct {truth:5g}, pruned query {got}")

    # removal via a LOWER triangle breaks queries for positive weights too
    g = make_graph(3, [(1, 0, 1), (0, 2, 1)])
    c = build(g, order=[0, 1, 2])
    c.perfect_customize()
    drop_up = bytearray(c.m)
    drop_dn = bytearray(c.m)
    a = c._find_arc(1, 2)
    drop_up[a] = 1                       # the only up-down representative of 1->0->2
    got = c.query_pruned(1, 2, drop_up, drop_dn)
    ok &= got == INF
    print(f"   lower-triangle pruning        positive     correct 2, pruned query {got}")

    # saturated-integer infinity: INF + (-5) looks finite
    big = 2 ** 31 - 1                    # int32 "infinity" with saturated addition
    v = min(big, big + (-5))             # the triangle computes INF + (-5) without an INF test
    ok &= v < big
    print(f"   saturated int32: INF + (-5) = {v} < INF, so an unreachable vertex looks reachable")
    return ok


# ----------------------------------------------- 2-6. randomized comparisons

def graph_family():
    fam = []
    for i in range(3):
        fam.append(("grid", graphs.ev_grid(10, seed=i)))
        fam.append(("geo", graphs.ev_geometric(150, seed=i)))
        fam.append(("shift", graphs.shifted_random(120, seed=i)))
    return fam


def check_random():
    rng = random.Random(3)
    stats = {k: [0, 0] for k in ("stall", "dijkstra", "prune", "early", "perfect", "witness", "potential")}
    stall_events = [0, 0]
    for kind, g0 in graph_family():
        for ci, g in enumerate((g0, abs_control(g0))):
            c = build(g, rng=random.Random(1))
            D = {}
            def dist(s):
                if s not in D:
                    D[s] = sssp(g, s)
                return D[s]
            for _ in range(40):
                s, t = rng.randrange(g.n), rng.randrange(g.n)
                truth = dist(s)[t]
                got, k = et_query(c, s, t, stall=True)
                stall_events[ci] += k
                stats["stall"][ci] += not close(got, truth)
                stats["dijkstra"][ci] += not close(dijkstra_query(c, s, t), truth)
                stats["prune"][ci] += not close(et_query(c, s, t, prune=True)[0], truth)
                stats["early"][ci] += not close(et_query(c, s, t, early_stop=True)[0], truth)
            # potential read off the CCH: feasible, and equal to Bellman-Ford's
            p = c.potential()
            bf = virtual_source_potential(g)
            bad = sum(1 for u in range(g.n) for v, w in g.adj[u] if w + p[u] - p[v] < -1e-9)
            stats["potential"][ci] += bad + sum(1 for v in range(g.n) if not close(p[v], bf[v]))
            # perfect customization on all arcs
            c.perfect_customize()
            inv = [0] * c.n
            for v in range(c.n):
                inv[c.rank[v]] = v
            for x in range(c.n):
                for a in range(c.first[x], c.first[x + 1]):
                    u, v = inv[x], inv[c.head[a]]
                    stats["perfect"][ci] += (not close(c.up[a], dist(u)[v])) + (not close(c.down[a], dist(v)[u]))
            # witness pruning (upper/intermediate triangles)
            drop_up, drop_dn = c.prunable_arcs()
            for _ in range(40):
                s, t = rng.randrange(g.n), rng.randrange(g.n)
                stats["witness"][ci] += not close(c.query_pruned(s, t, drop_up, drop_dn), dist(s)[t])
    print("2-6. randomized graphs (9 negative-weight graphs, 9 |w| controls; wrong answers)")
    names = {"stall": "stall-on-demand in the ET query", "dijkstra": "Dijkstra-based query (stop rule)",
             "prune": "ET query with distance pruning", "early": "ET query with early termination",
             "perfect": "perfect customization (all arcs)", "witness": "witness pruning (upper/intermediate)",
             "potential": "potential read off the CCH"}
    for k, (neg, ctl) in stats.items():
        print(f"   {names[k]:38s} negative {neg:5d}   control {ctl:5d}")
    print(f"   stall events: negative {stall_events[0]}, control {stall_events[1]}")
    must_be_zero = ("stall", "perfect", "witness", "potential")
    ok = all(stats[k] == [0, 0] for k in must_be_zero)
    ok &= all(stats[k][1] == 0 for k in ("dijkstra", "prune", "early"))   # exact for |w|
    return ok


def check_tiebreak():
    """Witness pruning with zero-length cycles, repaired by a lexicographic tie-break
    (integer weights w' = (n+1) w + 1: every cycle becomes strictly positive)."""
    r = random.Random(12)
    plain = lex = total = 0
    for _ in range(600):
        n = r.randint(3, 7)
        arcs = []
        for _ in range(r.randint(n, 3 * n)):
            u, v = r.randrange(n), r.randrange(n)
            if u != v:
                arcs.append((u, v, r.randint(-3, 3)))
        g = make_graph(n, arcs)
        if _has_negative_cycle(g):
            continue
        D = [sssp(g, s) for s in range(n)]
        order = list(range(n))
        r.shuffle(order)
        K = n + 1
        h = make_graph(n, [(u, v, w * K + 1) for u, v, w in arcs])
        for graph, is_lex in ((g, False), (h, True)):
            c = CCH(graph, log=NOLOG, order=order)
            c.customize(graph)
            c.perfect_customize()
            du, dd = c.prunable_arcs()
            for s in range(n):
                for t in range(n):
                    got = c.query_pruned(s, t, du, dd)
                    if is_lex and got != INF:
                        got = math.floor(got / K)
                    wrong = not close(got, D[s][t])
                    if is_lex:
                        lex += wrong
                    else:
                        plain += wrong
                        total += 1
    print(f"5b. tiny random graphs with zero cycles allowed: witness pruning wrong {plain} of {total}; "
          f"with lexicographic tie-break {lex}")
    return lex == 0


def virtual_source_potential(g):
    """Johnson's potential: Bellman-Ford from a virtual source with 0-arcs."""
    d = [0.0] * g.n
    inq = [True] * g.n
    q = deque(range(g.n))
    while q:
        u = q.popleft()
        inq[u] = False
        for v, w in g.adj[u]:
            if d[u] + w < d[v] - 1e-12:
                d[v] = d[u] + w
                if not inq[v]:
                    inq[v] = True
                    q.append(v)
    return d


def _has_negative_cycle(g):
    d = [0.0] * g.n
    for _ in range(g.n):
        changed = False
        for u in range(g.n):
            for v, w in g.adj[u]:
                if d[u] + w < d[v]:
                    d[v] = d[u] + w
                    changed = True
        if not changed:
            return False
    return True


# ----------------------------------------------------- 7. exhaustive 3 vertices

def check_exhaustive():
    pairs = [(u, v) for u in range(3) for v in range(3) if u != v]
    graphs_seen = 0
    fails = {"stall": 0, "dijkstra": 0, "prune": 0, "early": 0}
    for ws in itertools.product([None, -2, -1, 0, 1, 2], repeat=6):
        arcs = [(u, v, w) for (u, v), w in zip(pairs, ws) if w is not None]
        g = make_graph(3, arcs)
        if _has_negative_cycle(g):
            continue
        graphs_seen += 1
        c = build(g, order=[0, 1, 2])
        for s in range(3):
            D = sssp(g, s)
            for t in range(3):
                truth = D[t]
                fails["stall"] += not close(et_query(c, s, t, stall=True)[0], truth)
                fails["dijkstra"] += not close(dijkstra_query(c, s, t), truth)
                fails["prune"] += not close(et_query(c, s, t, prune=True)[0], truth)
                fails["early"] += not close(et_query(c, s, t, early_stop=True)[0], truth)
    print(f"7. exhaustive: {graphs_seen} conservative 3-vertex digraphs, weights in -2..2, all 9 queries each")
    for k, v in fails.items():
        print(f"   {k:9s} wrong answers: {v}")
    return fails["stall"] == 0


if __name__ == "__main__":
    ok = check_examples()
    ok &= check_random()
    ok &= check_tiebreak()
    if "--exhaustive" in sys.argv:
        ok &= check_exhaustive()
    print("\nall technique checks behave as stated" if ok else "\nUNEXPECTED RESULT")
    sys.exit(0 if ok else 1)
