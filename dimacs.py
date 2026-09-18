"""Load 9th DIMACS Challenge road graphs (USA-road-d.*) and add negative edges.

DIMACS weights are road lengths (all positive), so two weightings are offered:

* ``ev``: real topology and real lengths, plus synthetic smooth terrain.
  cost = len + BETA * dz (BETA_UP uphill, BETA_DOWN downhill), with elevation
  scaled so that steep-ish descents are negative. Coordinates are rescaled so
  straight-line distance never exceeds an edge's length, which keeps the
  domain bound a true lower bound.
* ``shifted``: cost = len + phi(u) - phi(v) with random phi (no geometry used).
"""
import gzip
import math
import random
import statistics
from pathlib import Path

import graphs
from graphs import BETA_DOWN, BETA_UP, Graph

DATA = Path(__file__).parent / "data"


def _read_raw(name):
    coords = {}
    with gzip.open(DATA / f"USA-road-d.{name}.co.gz", "rt") as f:
        for line in f:
            if line.startswith("v"):
                _, i, x, y = line.split()
                coords[int(i) - 1] = (int(x) * 1e-6, int(y) * 1e-6)  # lon, lat
    arcs = []
    with gzip.open(DATA / f"USA-road-d.{name}.gr.gz", "rt") as f:
        for line in f:
            if line.startswith("a"):
                _, u, v, w = line.split()
                arcs.append((int(u) - 1, int(v) - 1, int(w)))
    n = len(coords)
    lat0 = math.radians(statistics.fmean(c[1] for c in coords.values()))
    kx, ky = 111320 * math.cos(lat0), 110540  # degrees -> metres
    xy = [(coords[i][0] * kx, coords[i][1] * ky) for i in range(n)]
    return n, xy, arcs


def load(name, weighting="ev", seed=0, dz_scale=None):
    n, xy, arcs = _read_raw(name)
    rng = random.Random(seed)
    g = Graph(n)
    if weighting == "shifted":
        spread = 3 * statistics.median(w for _, _, w in arcs)
        phi = [rng.uniform(0, spread) for _ in range(n)]
        for u, v, w in arcs:
            if u != v:
                g.add_edge(u, v, w + phi[u] - phi[v])
        return g

    # scale metres -> length units using a low percentile of len/euclid,
    # then lengthen the few edges that are shorter than their straight line
    ratios = sorted(w / d for u, v, w in arcs if (d := math.dist(xy[u], xy[v])) > 1)
    scale = ratios[len(ratios) // 100]
    g.coords = [(x * scale, y * scale) for x, y in xy]
    size = max(max(c[0] for c in g.coords) - min(c[0] for c in g.coords),
               max(c[1] for c in g.coords) - min(c[1] for c in g.coords))
    x0 = min(c[0] for c in g.coords)
    y0 = min(c[1] for c in g.coords)
    terrain = graphs._make_terrain(rng, size, bumps=60, amp=1.0)
    raw_z = [graphs._terrain(x - x0, y - y0, terrain) for x, y in g.coords]
    if dz_scale is None:
        # choose the scale so roughly 10% of edges end up negative
        slopes = sorted(BETA_DOWN * abs(raw_z[v] - raw_z[u]) / w
                        for u, v, w in arcs if w > 0 and raw_z[v] < raw_z[u])
        dz_scale = 1 / slopes[int(len(slopes) * 0.8)]
    g.z = [z * dz_scale for z in raw_z]
    g.lengthened = 0
    for u, v, w in arcs:
        if u == v:
            continue
        e = math.dist(g.coords[u], g.coords[v])
        if w < e:
            w = e
            g.lengthened += 1
        dz = g.z[v] - g.z[u]
        g.add_edge(u, v, w + (BETA_UP if dz > 0 else BETA_DOWN) * dz)
    return g
