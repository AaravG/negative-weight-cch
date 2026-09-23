"""Single-pair query algorithms. Each returns (distance, nodes_expanded)."""
import heapq
from collections import deque

INF = float("inf")


def bellman_ford_query(g, s, t):
    """Queue-based Bellman-Ford (SPFA) from s. Ground truth and no-preprocessing baseline."""
    d = [INF] * g.n
    d[s] = 0.0
    in_q = [False] * g.n
    in_q[s] = True
    q = deque([s])
    expanded = 0
    adj = g.adj
    while q:
        u = q.popleft()
        in_q[u] = False
        expanded += 1
        du = d[u]
        for v, w in adj[u]:
            nd = du + w
            if nd < d[v] - 1e-12:
                d[v] = nd
                if not in_q[v]:
                    in_q[v] = True
                    q.append(v)
    return d[t], expanded


def astar(g, s, t, pi):
    """A* with a feasible potential pi (reduced costs >= 0): exact when t is popped."""
    g_cost = {s: 0.0}
    closed = set()
    pq = [(pi(s), s)]
    expanded = 0
    adj = g.adj
    while pq:
        _, u = heapq.heappop(pq)
        if u in closed:
            continue
        if u == t:
            return g_cost[u], expanded
        closed.add(u)
        expanded += 1
        gu = g_cost[u]
        for v, w in adj[u]:
            if v in closed:
                continue
            nd = gu + w
            if nd < g_cost.get(v, INF):
                g_cost[v] = nd
                heapq.heappush(pq, (nd + pi(v), v))
    return INF, expanded


def astar_reopening(g, s, t, h):
    """'Normal' textbook A* with node re-opening and an arbitrary heuristic.
    Stops when t is popped, so with an inadmissible heuristic it can be wrong."""
    g_cost = {s: 0.0}
    pq = [(h(s), 0.0, s)]
    expanded = 0
    adj = g.adj
    while pq:
        _, gu, u = heapq.heappop(pq)
        if gu > g_cost[u]:
            continue
        if u == t:
            return gu, expanded
        expanded += 1
        for v, w in adj[u]:
            nd = gu + w
            if nd < g_cost.get(v, INF):
                g_cost[v] = nd
                heapq.heappush(pq, (nd + h(v), nd, v))
    return INF, expanded


def bidirectional_astar(g, s, t, pi_t, pi_s):
    """Sequential bidirectional A* with the average potential (Ikeda et al.).

    pf = (pi_t - pi_s)/2 for the forward search, pr = -pf for the backward one.
    Both searches then see the same non-negative reduced costs, and we can stop
    when top_f + top_r >= mu.
    """
    if s == t:
        return 0.0, 0

    def pf(v):
        return 0.5 * (pi_t(v) - pi_s(v))

    df, dr = {s: 0.0}, {t: 0.0}
    cf, cr = set(), set()
    qf, qr = [(pf(s), s)], [(-pf(t), t)]
    mu = INF
    expanded = 0
    while qf and qr:
        if qf[0][0] + qr[0][0] >= mu:
            break
        forward = qf[0][0] <= qr[0][0]
        q, d, c, other, adj, sign = ((qf, df, cf, dr, g.adj, 1.0) if forward
                                     else (qr, dr, cr, df, g.radj, -1.0))
        _, u = heapq.heappop(q)
        if u in c:
            continue
        c.add(u)
        expanded += 1
        du = d[u]
        for v, w in adj[u]:
            nd = du + w
            if nd < d.get(v, INF):
                d[v] = nd
                heapq.heappush(q, (nd + sign * pf(v), v))
                ov = other.get(v)
                if ov is not None and nd + ov < mu:
                    mu = nd + ov
    return mu, expanded
