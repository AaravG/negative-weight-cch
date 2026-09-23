"""Inertial-flow graph bisection (Schild & Sommer 2015) for nested dissection.

For a cut direction, the first `balance` fraction of nodes along the
projection are sources and the last fraction are sinks. A minimum *vertex*
cut between them is found with max-flow (Dinic) on the node-split graph with
unit node capacities (sources and sinks have infinite capacity). The cut
vertices form the separator.
"""
from collections import deque

INF_CAP = 1 << 30


def _split(nodes, nbrs, xy, key, balance, st, qid):
    """st: shared stamp array (len n), st[v] == qid marks subgraph membership.
    Returns (side_A, side_B, separator)."""
    order = sorted(nodes, key=key)
    k = max(1, int(len(order) * balance))
    role = {}
    for v in order[:k]:
        role[v] = 1   # source
    for v in order[-k:]:
        role[v] = 2   # sink
    inside = lambda u: st[u] == qid
    # a source adjacent to a sink means no finite vertex cut exists
    for v in order[:k]:
        if any(role.get(u) == 2 and inside(u) for u in nbrs[v]):
            return None
    used = set()     # nodes whose internal (unit) arc carries flow
    net = {}          # net[(u, v)] = net flow on u_out -> v_in

    def neighbours(x):
        v, is_out = x >> 1, x & 1
        out = []
        if is_out:
            for u in nbrs[v]:
                if inside(u):
                    out.append(2 * u)               # out(v) -> in(u), infinite
            if v in used:
                out.append(2 * v)                   # undo internal flow
        else:
            if role.get(v) or v not in used:
                out.append(2 * v + 1)               # in(v) -> out(v)
            for u in nbrs[v]:
                if inside(u) and net.get((u, v), 0) > 0:
                    out.append(2 * u + 1)           # undo u_out -> v_in
        return out

    sources = [2 * v + 1 for v in order[:k]]
    flow = 0
    while True:
        # BFS levels from all sources
        level = {x: 0 for x in sources}
        dq = deque(sources)
        reached = False
        while dq:
            x = dq.popleft()
            if not x & 1 and role.get(x >> 1) == 2:
                reached = True
                continue
            for y in neighbours(x):
                if y not in level:
                    level[y] = level[x] + 1
                    dq.append(y)
        if not reached:
            break
        # blocking flow by iterative DFS
        adj_cache = {}
        ptr = {}
        dead = set()
        for s0 in sources:
            while True:
                path = [s0]
                while path:
                    x = path[-1]
                    if not x & 1 and role.get(x >> 1) == 2:
                        break
                    lst = adj_cache.get(x)
                    if lst is None:
                        lst = adj_cache[x] = [y for y in neighbours(x)
                                              if level.get(y) == level[x] + 1]
                    i = ptr.get(x, 0)
                    while i < len(lst) and lst[i] in dead:
                        i += 1
                    ptr[x] = i
                    if i == len(lst):
                        dead.add(x)
                        path.pop()
                        if path:
                            ptr[path[-1]] = ptr.get(path[-1], 0) + 1
                        continue
                    path.append(lst[i])
                if not path:
                    break
                # augment along path; unit capacity means one unit per path
                for a, b in zip(path, path[1:]):
                    va, vb = a >> 1, b >> 1
                    if va == vb:
                        if a & 1:          # out -> in: cancel internal flow
                            used.discard(va)
                        elif not role.get(va):
                            used.add(va)
                    elif a & 1:            # out(va) -> in(vb)
                        net[(va, vb)] = net.get((va, vb), 0) + 1
                        net[(vb, va)] = net.get((vb, va), 0) - 1
                    else:                  # in(va) -> out(vb): cancel vb->va
                        net[(vb, va)] = net.get((vb, va), 0) - 1
                        net[(va, vb)] = net.get((va, vb), 0) + 1
                flow += 1
                # arcs used by this path may be saturated now: rebuild the
                # cached lists of the path nodes lazily
                for x in path:
                    adj_cache.pop(x, None)
                    ptr.pop(x, None)
    # residual reachability from the sources gives the min cut
    seen = set(sources)
    dq = deque(sources)
    while dq:
        x = dq.popleft()
        for y in neighbours(x):
            if y not in seen:
                seen.add(y)
                dq.append(y)
    sep, A, B = [], [], []
    for v in nodes:
        i_in, i_out = 2 * v in seen, 2 * v + 1 in seen
        if i_in and not i_out:
            sep.append(v)
        elif i_out:
            A.append(v)
        else:
            B.append(v)
    return A, B, sep


def bisect(nodes, nbrs, xy, st, qid, balance=0.25):
    """Best of 4 directions: smallest separator, ties broken by balance."""
    best = None
    for key in (lambda v: xy[v][0], lambda v: xy[v][1],
                lambda v: xy[v][0] + xy[v][1], lambda v: xy[v][0] - xy[v][1]):
        res = _split(nodes, nbrs, xy, key, balance, st, qid)
        if res is None:
            continue
        A, B, sep = res
        score = (len(sep), abs(len(A) - len(B)))
        if best is None or score < best[0]:
            best = (score, A, B, sep)
    if best is None:
        return None
    return best[1], best[2], best[3]
