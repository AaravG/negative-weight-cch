"""Numba-compiled battery-constrained CCH (see cch_battery.py for the model).

Profiles are stored per shortcut direction d = 2*arc + (0 up, 1 down):
  * cnt[d] pieces; a single piece lives inline in (I1, C1, O1)[d],
  * profiles with more than one piece live in a shared growable pool
    (PI, PC, PO) at [pst[d], pst[d] + cnt[d]) with capacity pcap[d].
Triangles are enumerated on the fly, as in cch_nb.customize.
"""
import numpy as np
from numba import njit

NEG = -np.inf
EPS = 1e-9


@njit(cache=True, inline="always")
def _dominates(fi, fc, fo, gi, gc, go, M):
    if fi > gi + EPS:
        return False
    lo = gi - EPS
    hi = M + EPS
    for t in range(4):
        if t == 0:
            b = gi
        elif t == 1:
            b = M
        elif t == 2:
            b = fo + fc
        else:
            b = go + gc
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


@njit(cache=True)
def battery_customize(n, first, head, e_arc, e_up, w, M):
    ma = head.shape[0]
    nd = 2 * ma
    cnt = np.zeros(nd, np.int32)
    I1 = np.empty(nd)
    C1 = np.empty(nd)
    O1 = np.empty(nd)
    pst = np.full(nd, -1, np.int64)
    pcap = np.zeros(nd, np.int32)
    pool_cap = max(1 << 20, nd // 8)
    PI = np.empty(pool_cap)
    PC = np.empty(pool_cap)
    PO = np.empty(pool_cap)
    used = 0
    # scratch copies of source profiles (sources may move when the pool grows)
    sa_i = np.empty(1 << 16)
    sa_c = np.empty(1 << 16)
    sa_o = np.empty(1 << 16)
    sb_i = np.empty(1 << 16)
    sb_c = np.empty(1 << 16)
    sb_o = np.empty(1 << 16)

    for k in range(e_arc.shape[0]):
        a = e_arc[k]
        if a < 0:
            continue
        c = w[k]
        if c >= 0:
            if c > M:
                continue
            hi_, hc, ho = c, c, M - c
        else:
            hi_, hc, ho = 0.0, c, M
        d = 2 * a + (0 if e_up[k] else 1)
        # merge (shared code below, duplicated for speed)
        used, PI, PC, PO = _merge(d, hi_, hc, ho, M, cnt, I1, C1, O1, pst, pcap,
                                  PI, PC, PO, used)

    for x in range(n):
        s, e = first[x], first[x + 1]
        for i in range(s, e):
            y = head[i]
            j = first[y]
            je = first[y + 1]
            for kk in range(i + 1, e):
                z = head[kk]
                while j < je and head[j] < z:
                    j += 1
                # y -> x -> z : down(i) then up(kk)  => up(j)
                # z -> x -> y : down(kk) then up(i)  => down(j)
                for side in range(2):
                    if side == 0:
                        da = 2 * i + 1
                        db = 2 * kk
                        dt = 2 * j
                    else:
                        da = 2 * kk + 1
                        db = 2 * i
                        dt = 2 * j + 1
                    na = cnt[da]
                    nb = cnt[db]
                    if na == 0 or nb == 0:
                        continue
                    # copy sources into scratch
                    if na == 1:
                        sa_i[0] = I1[da]
                        sa_c[0] = C1[da]
                        sa_o[0] = O1[da]
                    else:
                        if na > sa_i.shape[0]:
                            sa_i = np.empty(2 * na)
                            sa_c = np.empty(2 * na)
                            sa_o = np.empty(2 * na)
                        p0 = pst[da]
                        for t in range(na):
                            sa_i[t] = PI[p0 + t]
                            sa_c[t] = PC[p0 + t]
                            sa_o[t] = PO[p0 + t]
                    if nb == 1:
                        sb_i[0] = I1[db]
                        sb_c[0] = C1[db]
                        sb_o[0] = O1[db]
                    else:
                        if nb > sb_i.shape[0]:
                            sb_i = np.empty(2 * nb)
                            sb_c = np.empty(2 * nb)
                            sb_o = np.empty(2 * nb)
                        p0 = pst[db]
                        for t in range(nb):
                            sb_i[t] = PI[p0 + t]
                            sb_c[t] = PC[p0 + t]
                            sb_o[t] = PO[p0 + t]
                    for ta in range(na):
                        fi, fc, fo = sa_i[ta], sa_c[ta], sa_o[ta]
                        for tb in range(nb):
                            gi, gc, go = sb_i[tb], sb_c[tb], sb_o[tb]
                            # compose f then g
                            if fo < gi - EPS:
                                continue
                            hi_ = fi if fi > gi + fc else gi + fc
                            if hi_ > M + EPS:
                                continue
                            hc = fc + gc
                            ho = go
                            if fo - gc < ho:
                                ho = fo - gc
                            if M - hc < ho:
                                ho = M - hc
                            used, PI, PC, PO = _merge(dt, hi_, hc, ho, M, cnt, I1, C1, O1,
                                                      pst, pcap, PI, PC, PO, used)
    return cnt, I1, C1, O1, pst, PI[:used].copy(), PC[:used].copy(), PO[:used].copy()


@njit(cache=True)
def _merge(d, hi_, hc, ho, M, cnt, I1, C1, O1, pst, pcap, PI, PC, PO, used):
    c0 = cnt[d]
    if c0 == 0:
        I1[d] = hi_
        C1[d] = hc
        O1[d] = ho
        cnt[d] = 1
        return used, PI, PC, PO
    if c0 == 1:
        gi, gc, go = I1[d], C1[d], O1[d]
        if _dominates(gi, gc, go, hi_, hc, ho, M):
            return used, PI, PC, PO
        if _dominates(hi_, hc, ho, gi, gc, go, M):
            I1[d] = hi_
            C1[d] = hc
            O1[d] = ho
            return used, PI, PC, PO
        # grow to a pool profile of capacity 4
        cap = 4
        if used + cap > PI.shape[0]:
            PI, PC, PO = _grow(PI, PC, PO, used + cap)
        pst[d] = used
        pcap[d] = cap
        PI[used] = gi
        PC[used] = gc
        PO[used] = go
        PI[used + 1] = hi_
        PC[used + 1] = hc
        PO[used + 1] = ho
        cnt[d] = 2
        return used + cap, PI, PC, PO
    p0 = pst[d]
    for t in range(c0):
        if _dominates(PI[p0 + t], PC[p0 + t], PO[p0 + t], hi_, hc, ho, M):
            return used, PI, PC, PO
    # remove pieces dominated by h (compact in place), then append h
    k = 0
    for t in range(c0):
        gi, gc, go = PI[p0 + t], PC[p0 + t], PO[p0 + t]
        if not _dominates(hi_, hc, ho, gi, gc, go, M):
            PI[p0 + k] = gi
            PC[p0 + k] = gc
            PO[p0 + k] = go
            k += 1
    if k == 0:
        I1[d] = hi_
        C1[d] = hc
        O1[d] = ho
        cnt[d] = 1
        return used, PI, PC, PO
    if k == pcap[d]:
        cap = 2 * pcap[d]
        if used + cap > PI.shape[0]:
            PI, PC, PO = _grow(PI, PC, PO, used + cap)
        for t in range(k):
            PI[used + t] = PI[p0 + t]
            PC[used + t] = PC[p0 + t]
            PO[used + t] = PO[p0 + t]
        p0 = used
        pst[d] = p0
        pcap[d] = cap
        used += cap
    PI[p0 + k] = hi_
    PC[p0 + k] = hc
    PO[p0 + k] = ho
    cnt[d] = k + 1
    return used, PI, PC, PO


@njit(cache=True)
def _grow(PI, PC, PO, need):
    cap = PI.shape[0]
    while cap < need:
        cap *= 2
    a = np.empty(cap)
    b = np.empty(cap)
    c = np.empty(cap)
    a[:PI.shape[0]] = PI
    b[:PC.shape[0]] = PC
    c[:PO.shape[0]] = PO
    return a, b, c


@njit(cache=True, inline="always")
def _eval(d, bx, cnt, I1, C1, O1, pst, PI, PC, PO):
    c0 = cnt[d]
    best = NEG
    if c0 == 1:
        i, c, o = I1[d], C1[d], O1[d]
        if bx >= i - EPS:
            v = bx - c
            best = o if o < v else v
    elif c0 > 1:
        p0 = pst[d]
        for t in range(c0):
            i, c, o = PI[p0 + t], PC[p0 + t], PO[p0 + t]
            if bx >= i - EPS:
                v = bx - c
                if o < v:
                    v = o
                if v > best:
                    best = v
    return best


@njit(cache=True)
def battery_query(s_rank, t_rank, b0, first, head, parent, cnt, I1, C1, O1, pst, PI, PC, PO,
                  soc, dsoc, mark, q):
    """Best state of charge at t when leaving s with b0 (NEG if infeasible)."""
    x = s_rank
    while x >= 0:
        soc[x] = NEG
        mark[x] = q
        x = parent[x]
    x = t_rank
    while x >= 0:
        dsoc[x] = NEG
        x = parent[x]
    soc[s_rank] = b0
    x = s_rank
    while x >= 0:
        bx = soc[x]
        if bx > NEG:
            for a in range(first[x], first[x + 1]):
                v = _eval(2 * a, bx, cnt, I1, C1, O1, pst, PI, PC, PO)
                y = head[a]
                if v > soc[y]:
                    soc[y] = v
        x = parent[x]
    # downward pass over t's ancestors in decreasing rank
    chain = np.empty(first.shape[0], np.int64)
    L = 0
    x = t_rank
    while x >= 0:
        chain[L] = x
        L += 1
        x = parent[x]
    for idx in range(L - 1, -1, -1):
        x = chain[idx]
        best = soc[x] if mark[x] == q else NEG
        for a in range(first[x], first[x + 1]):
            by = dsoc[head[a]]
            if by > NEG:
                v = _eval(2 * a + 1, by, cnt, I1, C1, O1, pst, PI, PC, PO)
                if v > best:
                    best = v
        dsoc[x] = best
    return dsoc[t_rank]


@njit(cache=True)
def battery_reference(n, off, tgt, w, M, s, b0):
    """Label-correcting max state of charge from s (exact; no gain cycles)."""
    soc = np.full(n, NEG)
    inq = np.zeros(n, np.bool_)
    cap = n + 1
    q = np.empty(cap, np.int64)
    head_, tail, size = 0, 1, 1
    q[0] = s
    soc[s] = b0
    inq[s] = True
    while size > 0:
        u = q[head_]
        head_ = (head_ + 1) % cap
        size -= 1
        inq[u] = False
        bu = soc[u]
        for k in range(off[u], off[u + 1]):
            c = w[k]
            if c >= 0:
                if bu < c - EPS:
                    continue
                nb = bu - c
            else:
                nb = bu - c
                if nb > M:
                    nb = M
            v = tgt[k]
            if nb > soc[v] + EPS:
                soc[v] = nb
                if not inq[v]:
                    inq[v] = True
                    q[tail] = v
                    tail = (tail + 1) % cap
                    size += 1
    return soc


@njit(cache=True)
def size_histogram(cnt):
    h = np.zeros(8, np.int64)   # 0, 1, 2-5, 6-10, 11-50, 51-100, >100, max
    mx = 0
    for d in range(cnt.shape[0]):
        c = cnt[d]
        if c > mx:
            mx = c
        if c == 0:
            h[0] += 1
        elif c == 1:
            h[1] += 1
        elif c <= 5:
            h[2] += 1
        elif c <= 10:
            h[3] += 1
        elif c <= 50:
            h[4] += 1
        elif c <= 100:
            h[5] += 1
        else:
            h[6] += 1
    h[7] = mx
    return h
