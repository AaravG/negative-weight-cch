"""Numba-compiled CCH for negative weights, for graphs up to continental size.

Same algorithms as cch.py / inertial_flow.py, rewritten on flat NumPy arrays:
  * undirected CSR with reverse-slot index,
  * nested dissection with inertial-flow bisection (Dinic max-flow on the
    node-split graph, unit vertex capacities), 4 directions per split,
  * symbolic elimination (elimination tree + upward arcs) without Python sets,
  * customization that enumerates lower triangles on the fly (not stored),
  * elimination-tree queries, negative-cycle scan,
  * Dijkstra with a feasible potential, used as the reference answer.
"""
import heapq

import numpy as np
from numba import njit

INF = np.inf


# ---------------------------------------------------------------- graph

@njit(cache=True)
def undirected_csr(n, off, tgt):
    """Symmetric, de-duplicated neighbour lists (no self loops) + reverse slots."""
    deg = np.zeros(n + 1, np.int64)
    for u in range(n):
        for k in range(off[u], off[u + 1]):
            v = tgt[k]
            if u != v:
                deg[u + 1] += 1
                deg[v + 1] += 1
    for i in range(n):
        deg[i + 1] += deg[i]
    pos = deg[:-1].copy()
    nb = np.empty(deg[n], np.int32)
    for u in range(n):
        for k in range(off[u], off[u + 1]):
            v = tgt[k]
            if u != v:
                nb[pos[u]] = v
                pos[u] += 1
                nb[pos[v]] = u
                pos[v] += 1
    # sort and de-duplicate each list
    noff = np.zeros(n + 1, np.int64)
    out = np.empty(deg[n], np.int32)
    w = 0
    for u in range(n):
        seg = np.sort(nb[deg[u]:deg[u + 1]])
        last = -1
        for v in seg:
            if v != last:
                out[w] = v
                w += 1
                last = v
        noff[u + 1] = w
    out = out[:w].copy()
    # reverse slot: position of u in v's list
    rev = np.empty(w, np.int64)
    for u in range(n):
        for k in range(noff[u], noff[u + 1]):
            v = out[k]
            lo, hi = noff[v], noff[v + 1]
            while lo < hi:
                mid = (lo + hi) // 2
                if out[mid] < u:
                    lo = mid + 1
                else:
                    hi = mid
            rev[k] = lo
    return noff, out, rev


# ---------------------------------------------------------------- inertial flow

@njit(cache=True)
def _flow_split(seg, noff, nb, rev, x, y, dirn, member, stamp, role,
                intf, fl, level, lvstamp, ptr, sep_out, a_out, b_out):
    """Min vertex cut between the first and last 25% of `seg` along `dirn`.
    Returns (n_sep, n_a, n_b) written into sep_out / a_out / b_out,
    or (-1, 0, 0) if a source is adjacent to a sink."""
    m = seg.shape[0]
    if dirn == 0:
        key = x[seg]
    elif dirn == 1:
        key = y[seg]
    elif dirn == 2:
        key = x[seg] + y[seg]
    else:
        key = x[seg] - y[seg]
    order = seg[np.argsort(key)]
    k = max(1, m // 4)
    for i in range(m):
        role[order[i]] = 0
    for i in range(k):
        role[order[i]] = 1
        role[order[m - 1 - i]] = 2
    for i in range(k):
        s = order[i]
        for j in range(noff[s], noff[s + 1]):
            v = nb[j]
            if member[v] == stamp and role[v] == 2:
                return -1, 0, 0
    # reset flow state on the subgraph
    for i in range(m):
        v = seg[i]
        intf[v] = 0
        for j in range(noff[v], noff[v + 1]):
            fl[j] = 0
    q = np.empty(2 * m, np.int64)
    stack = np.empty(2 * m + 1, np.int64)
    while True:
        # ---- BFS levels from all sources
        lv = lvstamp[0] + 1
        lvstamp[0] = lv
        qh = 0
        qt = 0
        # level[node] = (lv << 32) | depth, so it never needs clearing
        for i in range(k):
            s = order[i]
            for node in (2 * s, 2 * s + 1):
                level[node] = lv << 32
                q[qt] = node
                qt += 1
        reached = False
        while qh < qt:
            xnode = q[qh]
            qh += 1
            v = xnode >> 1
            d = level[xnode] & 0xFFFFFFFF
            if (xnode & 1) == 0 and role[v] == 2:
                reached = True
                continue
            # enumerate residual arcs
            if xnode & 1:  # out(v)
                if role[v] == 0 and intf[v] == 1:
                    tnode = 2 * v
                    if (level[tnode] >> 32) != lv:
                        level[tnode] = (lv << 32) | (d + 1)
                        q[qt] = tnode
                        qt += 1
                for j in range(noff[v], noff[v + 1]):
                    u = nb[j]
                    if member[u] == stamp:
                        tnode = 2 * u
                        if (level[tnode] >> 32) != lv:
                            level[tnode] = (lv << 32) | (d + 1)
                            q[qt] = tnode
                            qt += 1
            else:  # in(v)
                if role[v] != 0 or intf[v] == 0:
                    tnode = 2 * v + 1
                    if (level[tnode] >> 32) != lv:
                        level[tnode] = (lv << 32) | (d + 1)
                        q[qt] = tnode
                        qt += 1
                for j in range(noff[v], noff[v + 1]):
                    u = nb[j]
                    if member[u] == stamp and fl[j] < 0:
                        tnode = 2 * u + 1
                        if (level[tnode] >> 32) != lv:
                            level[tnode] = (lv << 32) | (d + 1)
                            q[qt] = tnode
                            qt += 1
        if not reached:
            break
        # reset arc pointers of the visited nodes
        for i in range(qt):
            ptr[q[i]] = 0
        # ---- blocking flow: iterative DFS with current-arc pointers
        for i in range(k):
            s = order[i]
            start = 2 * s + 1
            while True:
                sp = 0
                stack[0] = start
                found = False
                while sp >= 0:
                    xnode = stack[sp]
                    v = xnode >> 1
                    if (xnode & 1) == 0 and role[v] == 2:
                        found = True
                        break
                    d = level[xnode] & 0xFFFFFFFF
                    deg = noff[v + 1] - noff[v]
                    advanced = False
                    while ptr[xnode] <= deg:
                        a = ptr[xnode]
                        if a == 0:
                            # internal arc
                            if xnode & 1:
                                ok = role[v] == 0 and intf[v] == 1
                                tnode = 2 * v
                            else:
                                ok = role[v] != 0 or intf[v] == 0
                                tnode = 2 * v + 1
                        else:
                            j = noff[v] + a - 1
                            u = nb[j]
                            if member[u] != stamp:
                                ok = False
                                tnode = 0
                            elif xnode & 1:
                                ok = True
                                tnode = 2 * u
                            else:
                                ok = fl[j] < 0
                                tnode = 2 * u + 1
                        if ok and (level[tnode] >> 32) == lv and (level[tnode] & 0xFFFFFFFF) == d + 1:
                            sp += 1
                            stack[sp] = tnode
                            advanced = True
                            break
                        ptr[xnode] += 1
                    if not advanced:
                        level[xnode] = 0  # dead end
                        sp -= 1
                        if sp >= 0:
                            ptr[stack[sp]] += 1
                if not found:
                    break
                # augment one unit along stack[0..sp]
                for t in range(sp):
                    a_node = stack[t]
                    b_node = stack[t + 1]
                    va = a_node >> 1
                    vb = b_node >> 1
                    if va == vb:
                        if role[va] == 0:
                            intf[va] = 1 if (a_node & 1) == 0 else 0
                    else:
                        j = noff[va] + ptr[a_node] - 1
                        fl[j] += 1
                        fl[rev[j]] -= 1
    # ---- residual reachability from sources: min cut
    lv = lvstamp[0] + 1
    lvstamp[0] = lv
    qh = 0
    qt = 0
    for i in range(k):
        s = order[i]
        for node in (2 * s, 2 * s + 1):
            level[node] = lv << 32
            q[qt] = node
            qt += 1
    while qh < qt:
        xnode = q[qh]
        qh += 1
        v = xnode >> 1
        if xnode & 1:
            if role[v] == 0 and intf[v] == 1:
                tnode = 2 * v
                if (level[tnode] >> 32) != lv:
                    level[tnode] = lv << 32
                    q[qt] = tnode
                    qt += 1
            for j in range(noff[v], noff[v + 1]):
                u = nb[j]
                if member[u] == stamp:
                    tnode = 2 * u
                    if (level[tnode] >> 32) != lv:
                        level[tnode] = lv << 32
                        q[qt] = tnode
                        qt += 1
        else:
            if role[v] != 0 or intf[v] == 0:
                tnode = 2 * v + 1
                if (level[tnode] >> 32) != lv:
                    level[tnode] = lv << 32
                    q[qt] = tnode
                    qt += 1
            for j in range(noff[v], noff[v + 1]):
                u = nb[j]
                if member[u] == stamp and fl[j] < 0:
                    tnode = 2 * u + 1
                    if (level[tnode] >> 32) != lv:
                        level[tnode] = lv << 32
                        q[qt] = tnode
                        qt += 1
    ns = 0
    na = 0
    nb_ = 0
    for i in range(m):
        v = seg[i]
        rin = (level[2 * v] >> 32) == lv
        rout = (level[2 * v + 1] >> 32) == lv
        if rin and not rout:
            sep_out[ns] = v
            ns += 1
        elif rout:
            a_out[na] = v
            na += 1
        else:
            b_out[nb_] = v
            nb_ += 1
    return ns, na, nb_


@njit(cache=True)
def _geometric_split(seg, noff, nb, x, y, member, stamp, side, sep_out, a_out, b_out):
    """Fallback: median split along the wider axis; separator = A-side boundary."""
    m = seg.shape[0]
    xs = x[seg]
    ys = y[seg]
    if xs.max() - xs.min() >= ys.max() - ys.min():
        order = seg[np.argsort(xs)]
    else:
        order = seg[np.argsort(ys)]
    half = m // 2
    for i in range(m):
        side[order[i]] = 0 if i < half else 1
    ns = 0
    na = 0
    nb_ = 0
    for i in range(m):
        v = order[i]
        if side[v] == 1:
            b_out[nb_] = v
            nb_ += 1
            continue
        boundary = False
        for j in range(noff[v], noff[v + 1]):
            u = nb[j]
            if member[u] == stamp and side[u] == 1:
                boundary = True
                break
        if boundary:
            sep_out[ns] = v
            ns += 1
        else:
            a_out[na] = v
            na += 1
    return ns, na, nb_


@njit(cache=True)
def nested_dissection(n, noff, nb, rev, x, y, leaf):
    """Returns rank[v] (0 = eliminated first). Separators get the highest ranks."""
    buf = np.arange(n).astype(np.int64)
    member = np.zeros(n, np.int64)
    role = np.zeros(n, np.int8)
    side = np.zeros(n, np.int8)
    intf = np.zeros(n, np.int8)
    fl = np.zeros(nb.shape[0], np.int32)
    level = np.zeros(2 * n, np.int64)
    ptr = np.zeros(2 * n, np.int64)
    lvstamp = np.zeros(1, np.int64)
    sep_best = np.empty(n, np.int64)
    a_best = np.empty(n, np.int64)
    b_best = np.empty(n, np.int64)
    sep_t = np.empty(n, np.int64)
    a_t = np.empty(n, np.int64)
    b_t = np.empty(n, np.int64)
    # explicit stack of segments [lo, hi) in buf
    cap = 1024
    stk_lo = np.empty(cap, np.int64)
    stk_hi = np.empty(cap, np.int64)
    sp = 0
    stk_lo[0] = 0
    stk_hi[0] = n
    sp = 1
    stamp = 0
    degs = noff[1:] - noff[:-1]
    while sp > 0:
        sp -= 1
        lo = stk_lo[sp]
        hi = stk_hi[sp]
        m = hi - lo
        if m <= 0:
            continue
        seg = buf[lo:hi].copy()
        if m <= leaf:
            order = seg[np.argsort(degs[seg], kind="mergesort")]
            buf[lo:hi] = order
            continue
        stamp += 1
        for i in range(m):
            member[seg[i]] = stamp
        best_ns = -1
        best_na = 0
        best_nb = 0
        best_bal = 0
        for dirn in range(4):
            ns, na, nbb = _flow_split(seg, noff, nb, rev, x, y, dirn, member, stamp, role,
                                      intf, fl, level, lvstamp, ptr, sep_t, a_t, b_t)
            if ns < 0:
                continue
            bal = abs(na - nbb)
            if best_ns < 0 or ns < best_ns or (ns == best_ns and bal < best_bal):
                best_ns, best_na, best_nb, best_bal = ns, na, nbb, bal
                sep_best[:ns] = sep_t[:ns]
                a_best[:na] = a_t[:na]
                b_best[:nbb] = b_t[:nbb]
        if best_ns < 0 or best_na == 0 or best_nb == 0:
            ns, na, nbb = _geometric_split(seg, noff, nb, x, y, member, stamp, side,
                                           sep_best, a_best, b_best)
            best_ns, best_na, best_nb = ns, na, nbb
        # write A | B | S into the segment
        p = lo
        for i in range(best_na):
            buf[p] = a_best[i]
            p += 1
        for i in range(best_nb):
            buf[p] = b_best[i]
            p += 1
        for i in range(best_ns):
            buf[p] = sep_best[i]
            p += 1
        if sp + 2 >= cap:
            cap *= 2
            nl = np.empty(cap, np.int64)
            nh = np.empty(cap, np.int64)
            nl[:sp] = stk_lo[:sp]
            nh[:sp] = stk_hi[:sp]
            stk_lo = nl
            stk_hi = nh
        stk_lo[sp] = lo
        stk_hi[sp] = lo + best_na
        sp += 1
        stk_lo[sp] = lo + best_na
        stk_hi[sp] = lo + best_na + best_nb
        sp += 1
    rank = np.empty(n, np.int64)
    for i in range(n):
        rank[buf[i]] = i
    return rank


# ---------------------------------------------------------------- contraction

@njit(cache=True)
def symbolic(n, noff, nb, rank):
    """Upward arcs of the chordal completion, in rank space.
    Returns first (n+1), head (sorted per row), parent."""
    inv = np.empty(n, np.int64)
    for v in range(n):
        inv[rank[v]] = v
    parent = -np.ones(n, np.int64)
    child_head = -np.ones(n, np.int64)
    child_next = -np.ones(n, np.int64)
    marker = -np.ones(n, np.int64)
    first = np.zeros(n + 1, np.int64)
    cap = max(4 * n, 16)
    head = np.empty(cap, np.int32)
    tmp = np.empty(n, np.int64)
    used = 0
    for x in range(n):
        v = inv[x]
        cnt = 0
        for j in range(noff[v], noff[v + 1]):
            r = rank[nb[j]]
            if r > x and marker[r] != x:
                marker[r] = x
                tmp[cnt] = r
                cnt += 1
        c = child_head[x]
        while c >= 0:
            for a in range(first[c], first[c + 1]):
                r = head[a]
                if r != x and marker[r] != x:
                    marker[r] = x
                    tmp[cnt] = r
                    cnt += 1
            c = child_next[c]
        row = np.sort(tmp[:cnt])
        if used + cnt > cap:
            while used + cnt > cap:
                cap *= 2
            nh = np.empty(cap, np.int32)
            nh[:used] = head[:used]
            head = nh
        for i in range(cnt):
            head[used + i] = row[i]
        used += cnt
        first[x + 1] = used
        if cnt > 0:
            p = row[0]
            parent[x] = p
            child_next[x] = child_head[p]
            child_head[p] = x
    return first, head[:used].copy(), parent


@njit(cache=True)
def _find_arc(first, head, lo, hi):
    a, b = first[lo], first[lo + 1]
    while a < b:
        mid = (a + b) // 2
        if head[mid] < hi:
            a = mid + 1
        else:
            b = mid
    return a


@njit(cache=True)
def edge_arcs(n, off, tgt, rank, first, head):
    """For every original directed arc: shortcut index and direction (1 = up)."""
    m = off[n]
    arc = np.empty(m, np.int64)
    up = np.empty(m, np.int8)
    for u in range(n):
        ru = rank[u]
        for k in range(off[u], off[u + 1]):
            rv = rank[tgt[k]]
            if ru == rv:
                arc[k] = -1
                up[k] = 0
            elif ru < rv:
                arc[k] = _find_arc(first, head, ru, rv)
                up[k] = 1
            else:
                arc[k] = _find_arc(first, head, rv, ru)
                up[k] = 0
    return arc, up


@njit(cache=True)
def count_triangles(n, first, head):
    t = 0
    for x in range(n):
        d = first[x + 1] - first[x]
        t += d * (d - 1) // 2
    return t


# ---------------------------------------------------------------- customization

@njit(cache=True)
def customize(n, first, head, e_arc, e_up, w):
    ma = head.shape[0]
    upw = np.full(ma, INF)
    dnw = np.full(ma, INF)
    for k in range(e_arc.shape[0]):
        a = e_arc[k]
        if a < 0:
            continue
        if e_up[k]:
            if w[k] < upw[a]:
                upw[a] = w[k]
        elif w[k] < dnw[a]:
            dnw[a] = w[k]
    # lower triangles, enumerated on the fly, lowest vertex first
    for x in range(n):
        s, e = first[x], first[x + 1]
        for i in range(s, e):
            y = head[i]
            j = first[y]
            je = first[y + 1]
            d_yx = dnw[i]   # y -> x
            u_xy = upw[i]   # x -> y
            for kk in range(i + 1, e):
                z = head[kk]
                while j < je and head[j] < z:
                    j += 1
                # arc (y, z) is j (chordality)
                v = d_yx + upw[kk]      # y -> x -> z
                if v < upw[j]:
                    upw[j] = v
                v = dnw[kk] + u_xy      # z -> x -> y
                if v < dnw[j]:
                    dnw[j] = v
    return upw, dnw


@njit(cache=True)
def negative_cycle_arcs(upw, dnw):
    c = 0
    for a in range(upw.shape[0]):
        if upw[a] + dnw[a] < -1e-6:
            c += 1
    return c


# ---------------------------------------------------------------- queries

@njit(cache=True)
def query(s_rank, t_rank, first, head, parent, upw, dnw, df, dr, mark, q):
    x = s_rank
    while x >= 0:
        df[x] = INF
        mark[x] = q
        x = parent[x]
    x = t_rank
    while x >= 0:
        dr[x] = INF
        x = parent[x]
    df[s_rank] = 0.0
    dr[t_rank] = 0.0
    x = s_rank
    while x >= 0:
        dx = df[x]
        if dx < INF:
            for a in range(first[x], first[x + 1]):
                v = dx + upw[a]
                y = head[a]
                if v < df[y]:
                    df[y] = v
        x = parent[x]
    best = INF
    x = t_rank
    while x >= 0:
        dx = dr[x]
        if dx < INF:
            if mark[x] == q:
                v = dx + df[x]
                if v < best:
                    best = v
            for a in range(first[x], first[x + 1]):
                v = dx + dnw[a]
                y = head[a]
                if v < dr[y]:
                    dr[y] = v
        x = parent[x]
    return best


@njit(cache=True)
def dijkstra_potential(n, off, tgt, w, p, s, t):
    """Reference: Dijkstra on reduced weights w + p[u] - p[v] (>= 0)."""
    dist = np.full(n, INF)
    done = np.zeros(n, np.bool_)
    dist[s] = 0.0
    h = [(0.0, np.int64(s))]
    while len(h) > 0:
        d, u = heapq.heappop(h)
        if done[u]:
            continue
        done[u] = True
        if u == t:
            break
        pu = p[u]
        for k in range(off[u], off[u + 1]):
            v = np.int64(tgt[k])
            rw = w[k] + pu - p[v]
            if rw < 0.0:
                rw = 0.0
            nd = d + rw
            if nd < dist[v]:
                dist[v] = nd
                heapq.heappush(h, (nd, v))
    if dist[t] == INF:
        return INF
    return dist[t] - p[s] + p[t]


@njit(cache=True)
def spfa_potential(n, off, tgt, w):
    """Johnson potential by queue-based Bellman-Ford (for graphs without a cached one)."""
    p = np.zeros(n)
    inq = np.ones(n, np.bool_)
    q = np.empty(n + 1, np.int64)
    for i in range(n):
        q[i] = i
    head_, tail, size = 0, n % (n + 1), n
    cap = n + 1
    while size > 0:
        u = q[head_]
        head_ = (head_ + 1) % cap
        size -= 1
        inq[u] = False
        pu = p[u]
        for k in range(off[u], off[u + 1]):
            v = tgt[k]
            nd = pu + w[k]
            if nd < p[v] - 1e-9:
                p[v] = nd
                if not inq[v]:
                    inq[v] = True
                    q[tail] = v
                    tail = (tail + 1) % cap
                    size += 1
    return p
