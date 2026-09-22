"""Battery-constrained energy routing on a CCH, without any potential.

Model (Eisner, Funke & Storandt 2011; Baum et al. 2020)
------------------------------------------------------
Battery capacity M. Driving an arc with energy cost c (negative = recuperation)
maps state of charge (SoC) b to

    c >= 0:  b - c          if b >= c,  else infeasible
    c <  0:  min(M, b - c)                               (battery cannot overfill)

Every path induces a function of the same 3-parameter form

    f(b) = min(out, b - cost)   for b in [in, M],   infeasible for b < in

and composition stays in this family ("f then g"):

    feasible iff out_f >= in_g
    in   = max(in_f, in_g + cost_f)
    cost = cost_f + cost_g
    out  = min(out_g, out_f - cost_g, M - cost)

Why CCH customization stays exact
---------------------------------
Replace (min, +) by (max, o) on monotone functions: "best of two routes" is the
pointwise max, "route then route" is composition. Composition distributes over
max on both sides because all functions are non-decreasing. The analogue of
"no negative cycle" is: every closed walk function satisfies f(b) <= b, which
holds physically (a cycle has non-negative total energy and capping only loses
energy). Under these conditions the proofs of Theorems 1 and 2 of the draft go
through verbatim with ">= in SoC" in place of "<= in length".

Shortcut values are *profiles*: Pareto sets of triples, i.e. upper envelopes.

Query for a fixed initial SoC b0 needs only scalars:
  up:   ancestors of s in increasing rank,  soc[y] = max(soc[y], f_up(soc[x]))
  down: ancestors of t in decreasing rank, soc[x] = max(soc[x], f_down(soc[y]))
The answer is soc[t] after the down pass (-inf = unreachable with this charge).
"""
import time
from array import array
from collections import deque

NEG = float("-inf")
EPS = 1e-9


# ------------------------------------------------------------- functions

def arc_function(c, M):
    if c >= 0:
        return (c, c, M - c) if c <= M else None
    return (0.0, c, M)


def compose(f, g, M):
    """f then g, or None if infeasible."""
    in_f, c_f, o_f = f
    in_g, c_g, o_g = g
    if o_f < in_g - EPS:
        return None
    i = max(in_f, in_g + c_f)
    if i > M + EPS:
        return None
    c = c_f + c_g
    o = min(o_g, o_f - c_g, M - c)
    return (i, c, o)


def evaluate(f, b):
    i, c, o = f
    if b < i - EPS:
        return NEG
    v = b - c
    return o if o < v else v


def dominates(f, g, M):
    """f(b) >= g(b) for every b in [0, M]?

    Both functions are piecewise linear, so the difference is extreme at the
    breakpoints of either function or at the ends of g's feasible range."""
    fi, fc, fo = f
    gi, gc, go = g
    if fi > gi + EPS:
        return False
    lo, hi = gi - EPS, M + EPS
    for b in (gi, M, fo + fc, go + gc):
        if lo <= b <= hi:
            if b < fi - EPS:
                return False
            vf = b - fc
            if fo < vf:
                vf = fo
            vg = b - gc
            if go < vg:
                vg = go
            if vf < vg - EPS:
                return False
    return True


def prune(fs, M):
    """Upper envelope as a minimal Pareto set of triples."""
    fs = sorted(fs, key=lambda f: (f[0], f[1], -f[2]))
    out = []
    for f in fs:
        if any(dominates(g, f, M) for g in out):
            continue
        out = [g for g in out if not dominates(f, g, M)]
        out.append(f)
    return out


def merge_into(profile, h, M):
    """Add function h to a minimal Pareto set (list, modified or replaced).
    Returns the new list. Same result as prune(profile + [h]) up to ties."""
    for g in profile:
        if dominates(g, h, M):
            return profile
    if not profile:
        return [h]
    kept = [g for g in profile if not dominates(h, g, M)]
    kept.append(h)
    return kept


def eval_profile(p, b):
    best = NEG
    for f in p:
        v = evaluate(f, b)
        if v > best:
            best = v
    return best


def compose_profiles(p, q, M):
    res = []
    for f in p:
        for g in q:
            h = compose(f, g, M)
            if h is not None:
                res.append(h)
    return res


# ------------------------------------------------------------- CCH layer

class BatteryCCH:
    """Uses the metric-independent part of a cch.CCH (order, arcs, triangles)."""

    def __init__(self, cch):
        self.c = cch

    def customize(self, g, M):
        """Fast customization: incremental Pareto merge, single-function fast path."""
        c = self.c
        self.M = M
        up, down = self._base(g, M)
        ta, tb, tc = c.tri_a, c.tri_b, c.tri_c
        for i in range(len(ta)):
            a, b, cc = ta[i], tb[i], tc[i]
            da, ub = down[a], up[b]
            if da and ub:
                # y -> x -> z (down[a] then up[b])
                if len(da) == 1 and len(ub) == 1:
                    h = compose(da[0], ub[0], M)
                    if h is not None:
                        up[cc] = merge_into(up[cc], h, M)
                else:
                    tgt = up[cc]
                    for f in da:
                        for g2 in ub:
                            h = compose(f, g2, M)
                            if h is not None:
                                tgt = merge_into(tgt, h, M)
                    up[cc] = tgt
            db, ua = down[b], up[a]
            if db and ua:
                # z -> x -> y (down[b] then up[a])
                if len(db) == 1 and len(ua) == 1:
                    h = compose(db[0], ua[0], M)
                    if h is not None:
                        down[cc] = merge_into(down[cc], h, M)
                else:
                    tgt = down[cc]
                    for f in db:
                        for g2 in ua:
                            h = compose(f, g2, M)
                            if h is not None:
                                tgt = merge_into(tgt, h, M)
                    down[cc] = tgt
        self.up, self.down = up, down
        sizes = [len(p) for p in up] + [len(p) for p in down]
        return max(sizes), sum(sizes) / len(sizes)

    def _base(self, g, M):
        c = self.c
        up = [[] for _ in range(c.m)]
        down = [[] for _ in range(c.m)]
        k = 0
        for u in range(g.n):
            for _, w in g.adj[u]:
                a = c.edge_arc[k]
                k += 1
                if a < 0:
                    continue
                f = arc_function(w, M)
                if f is not None:
                    side = up if c.edge_dir[k - 1] else down
                    side[a] = merge_into(side[a], f, M)
        return up, down

    def customize_reference(self, g, M):
        """Original (slower) customization, kept to validate customize()."""
        c = self.c
        self.M = M
        up = [[] for _ in range(c.m)]
        down = [[] for _ in range(c.m)]
        k = 0
        for u in range(g.n):
            for _, w in g.adj[u]:
                a = c.edge_arc[k]
                k += 1
                if a < 0:
                    continue
                f = arc_function(w, M)
                if f is not None:
                    (up if c.edge_dir[k - 1] else down)[a].append(f)
        for a in range(c.m):
            if len(up[a]) > 1:
                up[a] = prune(up[a], M)
            if len(down[a]) > 1:
                down[a] = prune(down[a], M)
        ta, tb, tc = c.tri_a, c.tri_b, c.tri_c
        for i in range(len(ta)):
            a, b, cc = ta[i], tb[i], tc[i]
            # y -> x -> z  (down[a] then up[b])  and  z -> x -> y
            if down[a] and up[b]:
                new = compose_profiles(down[a], up[b], M)
                if new:
                    up[cc] = prune(up[cc] + new, M)
            if down[b] and up[a]:
                new = compose_profiles(down[b], up[a], M)
                if new:
                    down[cc] = prune(down[cc] + new, M)
        self.up, self.down = up, down
        sizes = [len(p) for p in up] + [len(p) for p in down]
        return max(sizes), sum(sizes) / len(sizes)

    def query(self, s, t, b0):
        """Best SoC at t when leaving s with b0 (NEG if unreachable)."""
        c = self.c
        parent, first, head = c.parent, c.first, c.head
        up, down = self.up, self.down
        soc = {}
        x = c.rank[s]
        soc[x] = b0
        while x >= 0:
            bx = soc.get(x, NEG)
            if bx > NEG:
                for a in range(first[x], first[x + 1]):
                    v = eval_profile(up[a], bx)
                    y = head[a]
                    if v > soc.get(y, NEG):
                        soc[y] = v
            x = parent[x]
        chain = []
        x = c.rank[t]
        while x >= 0:
            chain.append(x)
            x = parent[x]
        on_s = set()
        x = c.rank[s]
        while x >= 0:
            on_s.add(x)
            x = parent[x]
        down_soc = {}
        for x in reversed(chain):          # decreasing rank
            best = soc.get(x, NEG) if x in on_s else NEG
            for a in range(first[x], first[x + 1]):
                y = head[a]
                by = down_soc.get(y, NEG)
                if by > NEG:
                    v = eval_profile(down[a], by)
                    if v > best:
                        best = v
            down_soc[x] = best
        return down_soc[c.rank[t]]


# ------------------------------------------------------------- oracle

def battery_label_correcting(g, s, b0, M):
    """Max SoC at every vertex from s (reference; exact without gain cycles)."""
    soc = [NEG] * g.n
    soc[s] = b0
    inq = [False] * g.n
    dq = deque([s])
    inq[s] = True
    while dq:
        u = dq.popleft()
        inq[u] = False
        bu = soc[u]
        for v, w in g.adj[u]:
            if w >= 0:
                if bu < w - EPS:
                    continue
                nb = bu - w
            else:
                nb = min(M, bu - w)
            if nb > soc[v] + EPS:
                soc[v] = nb
                if not inq[v]:
                    inq[v] = True
                    dq.append(v)
    return soc
