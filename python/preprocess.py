"""Preprocessing: Johnson potential (Bellman-Ford) and ALT landmark tables."""
import heapq
import random
from collections import deque

INF = float("inf")


class NegativeCycleError(Exception):
    pass


def johnson_potential(g):
    """Bellman-Ford (queue-based) from a virtual source joined to every node by 0-edges.

    Returns p with w(u,v) + p[u] - p[v] >= 0 for every edge.
    """
    n = g.n
    p = [0.0] * n
    pred = [-1] * n
    in_q = [True] * n
    q = deque(range(n))
    pops = 0
    while q:
        u = q.popleft()
        in_q[u] = False
        pops += 1
        if pops % n == 0 and _pred_has_cycle(pred):
            raise NegativeCycleError("negative cycle detected")
        pu = p[u]
        for v, w in g.adj[u]:
            nd = pu + w
            if nd < p[v] - 1e-12:
                p[v] = nd
                pred[v] = u
                if not in_q[v]:
                    in_q[v] = True
                    q.append(v)
    return p


def _pred_has_cycle(pred):
    """A cycle in the predecessor graph is always a negative cycle.
    Checked every n relaxations, so the amortized cost is O(1) per step."""
    state = bytearray(len(pred))  # 0 new, 1 on current walk, 2 done
    for s in range(len(pred)):
        if state[s]:
            continue
        walk = []
        v = s
        while v >= 0 and not state[v]:
            state[v] = 1
            walk.append(v)
            v = pred[v]
        if v >= 0 and state[v] == 1:
            return True
        for x in walk:
            state[x] = 2
    return False


def dijkstra_reweighted(g, p, src, reverse=False):
    """Dijkstra with reduced weights w + p[u] - p[v]; returns TRUE distances."""
    n = g.n
    adj = g.radj if reverse else g.adj
    d = [INF] * n
    d[src] = 0.0
    pq = [(0.0, src)]
    while pq:
        du, u = heapq.heappop(pq)
        if du > d[u]:
            continue
        pu = p[u]
        for v, w in adj[u]:
            # forward edge u->v, or (reverse) original edge v->u
            rw = (w + pu - p[v]) if not reverse else (w + p[v] - pu)
            nd = du + max(0.0, rw)
            if nd < d[v]:
                d[v] = nd
                heapq.heappush(pq, (nd, v))
    # convert reduced distances back to true distances
    ps = p[src]
    if not reverse:
        return [dv - ps + p[v] if dv < INF else INF for v, dv in enumerate(d)]
    return [dv - p[v] + ps if dv < INF else INF for v, dv in enumerate(d)]


class ALTTables:
    """Landmark distances: from_L[i][v] = d(L_i, v), to_L[i][v] = d(v, L_i)."""

    def __init__(self, g, p, k=8, seed=0):
        rng = random.Random(seed)
        self.landmarks = []
        self.from_L, self.to_L = [], []
        # farthest-point selection on reduced (hop-ish) distances
        closest = [INF] * g.n
        cur = rng.randrange(g.n)
        for _ in range(k):
            self.landmarks.append(cur)
            fr = dijkstra_reweighted(g, p, cur)
            to = dijkstra_reweighted(g, p, cur, reverse=True)
            self.from_L.append(fr)
            self.to_L.append(to)
            # next landmark: node farthest (in reduced distance) from all chosen ones
            pl = p[cur]
            best, nxt = -1.0, cur
            for v in range(g.n):
                if fr[v] < INF:
                    closest[v] = min(closest[v], fr[v] + pl - p[v])
                if best < closest[v] < INF:
                    best, nxt = closest[v], v
            cur = nxt
