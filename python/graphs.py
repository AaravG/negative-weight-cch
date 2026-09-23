"""Graph container and benchmark graph generators.

All generators produce graphs with negative edges but no negative cycles.

* ``ev_grid`` / ``ev_geometric``: electric-vehicle energy model. Driving an edge
  costs ``ALPHA * length + beta * dz`` where ``beta = BETA_UP`` uphill and
  ``BETA_DOWN`` downhill (regenerative braking recovers less than climbing
  costs). Steep downhill edges are negative, every cycle is positive.
* ``shifted_random``: random sparse digraph with no geometry. Weights are
  ``c + phi(u) - phi(v)`` with ``c > 0``, so cycles are positive but many
  edges are negative. Only the ALT heuristic applies here.
"""
import math
import random

ALPHA = 1.0      # energy per unit distance
BETA_UP = 1.0    # energy per unit of climb
BETA_DOWN = 0.6  # energy recovered per unit of descent


class Graph:
    def __init__(self, n):
        self.n = n
        self.adj = [[] for _ in range(n)]   # adj[u]  = [(v, w), ...]
        self.radj = [[] for _ in range(n)]  # radj[v] = [(u, w), ...]
        self.coords = None                  # [(x, y)] for geometric graphs
        self.z = None                       # elevation for EV graphs

    def add_edge(self, u, v, w):
        self.adj[u].append((v, w))
        self.radj[v].append((u, w))

    @property
    def m(self):
        return sum(len(a) for a in self.adj)

    def negative_fraction(self):
        neg = sum(1 for a in self.adj for _, w in a if w < 0)
        return neg / max(1, self.m)


def _terrain(x, y, rng_params):
    """Smooth hilly terrain: sum of a few random bumps."""
    h = 0.0
    for cx, cy, amp, rad in rng_params:
        h += amp * math.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2 * rad * rad))
    return h


def _ev_weight(p, q, zp, zq):
    dz = zq - zp
    beta = BETA_UP if dz > 0 else BETA_DOWN
    return ALPHA * math.dist(p, q) + beta * dz


def _make_terrain(rng, size, bumps=12, amp=None):
    amp = amp if amp is not None else size * 0.6
    return [(rng.uniform(0, size), rng.uniform(0, size),
             rng.uniform(-amp, amp), rng.uniform(size * 0.05, size * 0.2))
            for _ in range(bumps)]


def ev_grid(side, seed=0):
    """side x side 8-connected grid on hilly terrain."""
    rng = random.Random(seed)
    terrain = _make_terrain(rng, side)
    g = Graph(side * side)
    g.coords = [(i % side, i // side) for i in range(g.n)]
    g.z = [_terrain(x, y, terrain) for x, y in g.coords]
    for y in range(side):
        for x in range(side):
            u = y * side + x
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1),
                           (1, 1), (1, -1), (-1, 1), (-1, -1)):
                nx, ny = x + dx, y + dy
                if 0 <= nx < side and 0 <= ny < side:
                    v = ny * side + nx
                    g.add_edge(u, v, _ev_weight(g.coords[u], g.coords[v], g.z[u], g.z[v]))
    return g


def ev_geometric(n, k=6, seed=0):
    """Random points, symmetric k-nearest-neighbour roads, largest component."""
    rng = random.Random(seed)
    size = math.sqrt(n) * 1.0
    pts = [(rng.uniform(0, size), rng.uniform(0, size)) for _ in range(n)]
    # bucket grid for kNN
    cell = 2.0
    buckets = {}
    for i, (x, y) in enumerate(pts):
        buckets.setdefault((int(x // cell), int(y // cell)), []).append(i)
    nbrs = [set() for _ in range(n)]
    for i, (x, y) in enumerate(pts):
        cx, cy, r = int(x // cell), int(y // cell), 1
        while True:
            cand = [j for dx in range(-r, r + 1) for dy in range(-r, r + 1)
                    for j in buckets.get((cx + dx, cy + dy), ()) if j != i]
            if len(cand) >= k or r > 50:
                break
            r += 1
        cand.sort(key=lambda j: math.dist(pts[i], pts[j]))
        for j in cand[:k]:
            nbrs[i].add(j)
            nbrs[j].add(i)
    # largest connected component
    comp, best = [-1] * n, []
    for s in range(n):
        if comp[s] != -1:
            continue
        stack, members = [s], []
        comp[s] = s
        while stack:
            u = stack.pop()
            members.append(u)
            for v in nbrs[u]:
                if comp[v] == -1:
                    comp[v] = s
                    stack.append(v)
        if len(members) > len(best):
            best = members
    remap = {old: new for new, old in enumerate(sorted(best))}
    terrain = _make_terrain(rng, size)
    g = Graph(len(remap))
    g.coords = [None] * g.n
    for old, new in remap.items():
        g.coords[new] = pts[old]
    g.z = [_terrain(x, y, terrain) for x, y in g.coords]
    for old, u in remap.items():
        for oj in nbrs[old]:
            v = remap[oj]
            g.add_edge(u, v, _ev_weight(g.coords[u], g.coords[v], g.z[u], g.z[v]))
    return g


def shifted_random(n, avg_deg=4, seed=0):
    """Random digraph (Hamiltonian cycle + random edges) with potential-shifted weights."""
    rng = random.Random(seed)
    phi = [rng.uniform(0, 300) for _ in range(n)]
    g = Graph(n)
    order = list(range(n))
    rng.shuffle(order)
    edges = set()
    for i in range(n):  # cycle guarantees strong connectivity
        edges.add((order[i], order[(i + 1) % n]))
    while len(edges) < n * avg_deg:
        u, v = rng.randrange(n), rng.randrange(n)
        if u != v:
            edges.add((u, v))
    for u, v in edges:
        g.add_edge(u, v, rng.uniform(1, 100) + phi[u] - phi[v])
    return g
