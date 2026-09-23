"""Load 9th DIMACS Challenge road graphs (USA-road-d.*) and add negative edges.

DIMACS weights are road lengths (all positive), so two weightings are offered:

* ``ev``: real topology and real lengths, plus synthetic smooth terrain.
  cost = len + BETA * dz (BETA_UP uphill, BETA_DOWN downhill), with elevation
  scaled so that steep-ish descents are negative. Coordinates are rescaled so
  straight-line distance never exceeds an edge's length, which keeps the
  domain bound a true lower bound.
* ``shifted``: cost = len + phi(u) - phi(v) with random phi (no geometry used).
* ``ev_real``: like ``ev`` but with real elevation from the AWS terrain tiles
  (see elevation.py), scaled by CLIMB_PER_METRE: climbing one metre costs as
  much as driving that many metres, so descents steeper than about 3% become
  negative. Sea-floor values (roads sampling water pixels) are clamped to 0.
"""
import gzip
import math
import random
import statistics
from pathlib import Path

import graphs
from graphs import BETA_DOWN, BETA_UP, Graph

DATA = Path(__file__).parent / "data"
CLIMB_PER_METRE = 50.0   # energy of 1 m of climb, in metres of driving


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


# Regions cut out of the full USA graph (lon/lat bounding boxes), largest
# connected component only. Requires data/USA-road-d.USA.{gr,co}.gz.
REGIONS = {
    "FLA_R": (-87.7, 24.4, -79.8, 31.1),    # Florida (~1M nodes)
    "CAL_R": (-124.5, 32.5, -114.0, 42.0),  # California + western Nevada (~2M nodes)
}


def _read_region(bbox):
    lon0, lat0_, lon1, lat1 = bbox
    keep = {}
    pts = []
    with gzip.open(DATA / "USA-road-d.USA.co.gz", "rb") as f:
        for line in f:
            if line[:1] == b"v":
                _, i, x, y = line.split()
                lon, lat = int(x) * 1e-6, int(y) * 1e-6
                if lon0 <= lon <= lon1 and lat0_ <= lat <= lat1:
                    keep[int(i) - 1] = len(pts)
                    pts.append((lon, lat))
    arcs = []
    with gzip.open(DATA / "USA-road-d.USA.gr.gz", "rb") as f:
        for line in f:
            if line[:1] == b"a":
                _, u, v, w = line.split()
                u, v = keep.get(int(u) - 1), keep.get(int(v) - 1)
                if u is not None and v is not None:
                    arcs.append((u, v, int(w)))
    del keep
    # largest connected component (undirected)
    n0 = len(pts)
    nbr = [[] for _ in range(n0)]
    for u, v, _ in arcs:
        nbr[u].append(v)
        nbr[v].append(u)
    comp = [-1] * n0
    best, best_id = 0, -1
    for s0 in range(n0):
        if comp[s0] != -1:
            continue
        comp[s0] = s0
        stack, size = [s0], 0
        while stack:
            x = stack.pop()
            size += 1
            for y in nbr[x]:
                if comp[y] == -1:
                    comp[y] = s0
                    stack.append(y)
        if size > best:
            best, best_id = size, s0
    del nbr
    remap = [-1] * n0
    lonlat = []
    for i in range(n0):
        if comp[i] == best_id:
            remap[i] = len(lonlat)
            lonlat.append(pts[i])
    arcs = [(remap[u], remap[v], w) for u, v, w in arcs if remap[u] >= 0 and remap[v] >= 0]
    return lonlat, arcs


def _read(name):
    if name in REGIONS:
        lonlat, arcs = _read_region(REGIONS[name])
        n = len(lonlat)
        coords = dict(enumerate(lonlat))
    else:
        return _read_raw(name)
    lat0 = math.radians(statistics.fmean(c[1] for c in coords.values()))
    kx, ky = 111320 * math.cos(lat0), 110540
    xy = [(coords[i][0] * kx, coords[i][1] * ky) for i in range(n)]
    return n, xy, arcs


def load(name, weighting="ev", seed=0, dz_scale=None, zoom=11):
    n, xy, arcs = _read(name)
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
    if weighting == "ev_real":
        import elevation
        metres = elevation.sample(raw_lonlat(name), zoom, log=lambda *_: None)
        raw_z = [max(0.0, e) for e in metres]          # clamp sea floor to 0
        dz_scale = CLIMB_PER_METRE * scale             # metres -> length units
    else:
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

def raw_lonlat(name):
    """(lon, lat) per node, in graph order - used by elevation.py."""
    if name in REGIONS:
        lonlat, _ = _read_region(REGIONS[name])
        return lonlat
    coords = {}
    with gzip.open(DATA / f"USA-road-d.{name}.co.gz", "rt") as f:
        for line in f:
            if line.startswith("v"):
                _, i, x, y = line.split()
                coords[int(i) - 1] = (int(x) * 1e-6, int(y) * 1e-6)
    return [coords[i] for i in range(len(coords))]
