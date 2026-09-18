"""Heuristics expressed as feasible potentials on the ORIGINAL (negative) weights.

For a target t, ``to_target(t)`` returns pi_t with
    w(u,v) + pi_t(v) - pi_t(u) >= 0   for every edge   (consistency)
so A* with key g(v) + pi_t(v) is Dijkstra on non-negative reduced costs and is
exact when t is popped. For a source s, ``from_source(s)`` returns pi_s with
    w(u,v) + pi_s(u) - pi_s(v) >= 0
which is what the backward half of a bidirectional search needs.

Lower bounds (ALT, domain, combined) additionally satisfy pi_t(v) <= d(v,t).
The Johnson potential is feasible but is not a lower bound, so it gives no
pull toward the target: it is the "Johnson + Dijkstra" baseline.
"""
import math

from graphs import ALPHA, BETA_DOWN

INF = float("inf")
_SHRINK = 1 - 1e-9  # guard against float rounding breaking consistency


class JohnsonPotential:
    name = "johnson"

    def __init__(self, p):
        self.p = p

    def to_target(self, t):
        p = self.p
        return lambda v: -p[v]

    def from_source(self, s):
        p = self.p
        return lambda v: p[v]


class ALTPotential:
    name = "alt"

    def __init__(self, tables):
        self.fr = tables.from_L
        self.to = tables.to_L

    def to_target(self, t):
        # d(v,t) >= d(v,L) - d(t,L)   and   d(v,t) >= d(L,t) - d(L,v)
        pairs = [(to, to[t], fr, fr[t]) for fr, to in zip(self.fr, self.to)]
        cache = {}

        def pi(v):
            r = cache.get(v)
            if r is None:
                r = -INF
                for to, to_t, fr, fr_t in pairs:
                    a = to[v] - to_t
                    b = fr_t - fr[v]
                    if a > r:
                        r = a
                    if b > r:
                        r = b
                cache[v] = r
            return r
        return pi

    def from_source(self, s):
        # d(s,v) >= d(L,v) - d(L,s)   and   d(s,v) >= d(s,L) - d(v,L)
        pairs = [(fr, fr[s], to, to[s]) for fr, to in zip(self.fr, self.to)]
        cache = {}

        def pi(v):
            r = cache.get(v)
            if r is None:
                r = -INF
                for fr, fr_s, to, to_s in pairs:
                    a = fr[v] - fr_s
                    b = to_s - to[v]
                    if a > r:
                        r = a
                    if b > r:
                        r = b
                cache[v] = r
            return r
        return pi


class DomainPotential:
    """EV energy bound: every edge costs at least ALPHA*euclid + BETA_DOWN*dz,
    which telescopes to d(v,t) >= ALPHA*|v t| + BETA_DOWN*(z_t - z_v).
    Needs no preprocessing at all."""
    name = "domain"

    def __init__(self, g):
        self.coords, self.z = g.coords, g.z

    def to_target(self, t):
        c, z = self.coords, self.z
        ct, zt = c[t], z[t]
        return lambda v: ALPHA * _SHRINK * math.dist(c[v], ct) + BETA_DOWN * (zt - z[v])

    def from_source(self, s):
        c, z = self.coords, self.z
        cs, zs = c[s], z[s]
        return lambda v: ALPHA * _SHRINK * math.dist(cs, c[v]) + BETA_DOWN * (z[v] - zs)


class MaxPotential:
    """Pointwise max of consistent lower bounds is still a consistent lower bound."""

    def __init__(self, *parts):
        self.parts = parts
        self.name = "max(" + ",".join(p.name for p in parts) + ")"

    def _combine(self, fns):
        cache = {}

        def pi(v):
            r = cache.get(v)
            if r is None:
                r = max(f(v) for f in fns)
                cache[v] = r
            return r
        return pi

    def to_target(self, t):
        return self._combine([p.to_target(t) for p in self.parts])

    def from_source(self, s):
        return self._combine([p.from_source(s) for p in self.parts])


def naive_euclid(g):
    """'Normal A*' heuristic: straight-line energy, ignoring elevation.
    NOT admissible on these graphs (downhill edges cost less than distance)."""
    c = g.coords

    def make(t):
        ct = c[t]
        return lambda v: ALPHA * math.dist(c[v], ct)
    return make
