"""Search algorithms on memory-mapped CSR graphs, plus picklable job functions
for a process pool. Same algorithms as search.py / heuristics.py / parallel.py,
rewritten for flat arrays so they fit in memory on the full USA graph."""
import ctypes
import heapq
import math
import sys
import time
from array import array
from collections import deque

from bigcsr import CACHE, CSR
from graphs import BETA_DOWN

INF = float("inf")
SHRINK = 1 - 1e-9
N_LANDMARKS = 8

_graphs = {}
_buffers = {}
_qid = [0]


def low_priority():
    if sys.platform == "win32":
        k = ctypes.windll.kernel32
        k.SetPriorityClass(k.GetCurrentProcess(), 0x4000)  # BELOW_NORMAL


def graph(weighting):
    g = _graphs.get(weighting)
    if g is None:
        g = _graphs[weighting] = CSR(weighting)
    return g


def _bufs(n, key):
    b = _buffers.get(key)
    if b is None or len(b[0]) != n:
        b = _buffers[key] = (array("d", bytes(8 * n)), array("i", bytes(4 * n)),
                             array("i", bytes(4 * n)))
    return b


def _next_qid():
    _qid[0] += 1
    return _qid[0]


def _mapped(fname):
    from bigcsr import _map
    return _map(fname, "d")


# ---------------------------------------------------------------- potentials

def landmark_file(weighting, i, reverse):
    return f"{weighting}_L{i}_{'to' if reverse else 'from'}.bin"


class Potential:
    """kind in {'johnson', 'alt', 'domain', 'max'} (max = max(alt, domain))."""

    def __init__(self, g, kind):
        self.g, self.kind = g, kind
        w = g.weighting
        if kind == "johnson":
            self.p = _mapped(f"{w}_p.bin")
        if kind in ("alt", "max"):
            self.fr = [_mapped(landmark_file(w, i, False)) for i in range(N_LANDMARKS)]
            self.to = [_mapped(landmark_file(w, i, True)) for i in range(N_LANDMARKS)]

    def _alt(self, a, a_key, b, b_key):
        pairs = [(A, A[a_key], B, B[b_key]) for A, B in zip(a, b)]

        def pi(v):
            r = -INF
            for A, Aa, B, Bb in pairs:
                x = A[v] - Aa
                y = Bb - B[v]
                if x > r:
                    r = x
                if y > r:
                    r = y
            return r
        return pi

    def _domain(self, anchor, sign):
        g = self.g
        xm, ym, z, S = g.xm, g.ym, g.z, g.scale * SHRINK
        ax, ay, az = xm[anchor], ym[anchor], z[anchor]
        hyp = math.hypot
        if sign > 0:   # lower bound on d(v, anchor)
            return lambda v: S * hyp(xm[v] - ax, ym[v] - ay) + BETA_DOWN * (az - z[v])
        return lambda v: S * hyp(xm[v] - ax, ym[v] - ay) + BETA_DOWN * (z[v] - az)

    def _cached(self, f):
        cache = {}

        def pi(v):
            r = cache.get(v)
            if r is None:
                r = cache[v] = f(v)
            return r
        return pi

    def to_target(self, t):
        k = self.kind
        if k == "johnson":
            p = self.p
            return lambda v: -p[v]
        if k == "domain":
            return self._domain(t, +1)
        alt = self._alt(self.to, t, self.fr, t)  # d(v,L)-d(t,L), d(L,t)-d(L,v)
        if k == "alt":
            return self._cached(alt)
        dom = self._domain(t, +1)
        return self._cached(lambda v: max(alt(v), dom(v)))

    def from_source(self, s):
        k = self.kind
        if k == "johnson":
            p = self.p
            return lambda v: p[v]
        if k == "domain":
            return self._domain(s, -1)
        alt = self._alt(self.fr, s, self.to, s)  # d(L,v)-d(L,s), d(s,L)-d(v,L)
        if k == "alt":
            return self._cached(alt)
        dom = self._domain(s, -1)
        return self._cached(lambda v: max(alt(v), dom(v)))


# ---------------------------------------------------------------- algorithms

def astar(g, s, t, pi, max_expanded=None, reopen=False):
    d, st, cl = _bufs(g.n, "f")
    q = _next_qid()
    off, tgt, w = g.off, g.tgt, g.fw
    d[s] = 0.0
    st[s] = q
    pq = [(pi(s), 0.0, s)]
    expanded = 0
    pop, push = heapq.heappop, heapq.heappush
    while pq:
        _, gu, u = pop(pq)
        if reopen:
            if gu > d[u]:
                continue
        elif cl[u] == q:
            continue
        if u == t:
            return d[u], expanded
        cl[u] = q
        expanded += 1
        if max_expanded and expanded > max_expanded:
            return None, expanded
        gu = d[u]
        for i in range(off[u], off[u + 1]):
            v = tgt[i]
            if not reopen and cl[v] == q:
                continue
            nd = gu + w[i]
            if st[v] != q or nd < d[v]:
                d[v] = nd
                st[v] = q
                push(pq, (nd + pi(v), nd, v))
    return INF, expanded


def bidirectional(g, s, t, pot):
    if s == t:
        return 0.0, 0
    pt, ps = pot.to_target(t), pot.from_source(s)

    def pf(v):
        return 0.5 * (pt(v) - ps(v))

    q = _next_qid()
    sides = []
    for key, start, off, nbr, w, sign in (("f", s, g.off, g.tgt, g.fw, 1.0),
                                          ("r", t, g.roff, g.rsrc, g.rw, -1.0)):
        d, st, cl = _bufs(g.n, key)
        d[start] = 0.0
        st[start] = q
        sides.append([[(sign * pf(start), start)], d, st, cl, off, nbr, w, sign])
    mu = INF
    expanded = 0
    while sides[0][0] and sides[1][0]:
        if sides[0][0][0][0] + sides[1][0][0][0] >= mu:
            break
        me = 0 if sides[0][0][0][0] <= sides[1][0][0][0] else 1
        pq, d, st, cl, off, nbr, w, sign = sides[me]
        od, ost = sides[1 - me][1], sides[1 - me][2]
        _, u = heapq.heappop(pq)
        if cl[u] == q:
            continue
        cl[u] = q
        expanded += 1
        du = d[u]
        for i in range(off[u], off[u + 1]):
            v = nbr[i]
            nd = du + w[i]
            if st[v] != q or nd < d[v]:
                d[v] = nd
                st[v] = q
                heapq.heappush(pq, (nd + sign * pf(v), v))
                if ost[v] == q and nd + od[v] < mu:
                    mu = nd + od[v]
    return mu, expanded


def spfa_query(g, s, t, time_cap):
    d = array("d", [INF]) * g.n
    inq = bytearray(g.n)
    d[s] = 0.0
    inq[s] = 1
    dq = deque([s])
    off, tgt, w = g.off, g.tgt, g.fw
    expanded = 0
    start = time.perf_counter()
    while dq:
        u = dq.popleft()
        inq[u] = 0
        expanded += 1
        if expanded & 0xFFFFF == 0 and time.perf_counter() - start > time_cap:
            return None, expanded
        du = d[u]
        for i in range(off[u], off[u + 1]):
            v = tgt[i]
            nd = du + w[i]
            if nd < d[v] - 1e-9:
                d[v] = nd
                if not inq[v]:
                    inq[v] = 1
                    dq.append(v)
    return d[t], expanded


# ---------------------------------------------------------------- pool jobs

def job_johnson(weighting):
    low_priority()
    g = graph(weighting)
    t0 = time.perf_counter()
    p = array("d", bytes(8 * g.n))
    inq = bytearray(b"\x01") * g.n
    cnt = array("i", bytes(4 * g.n))
    dq = deque(range(g.n))
    off, tgt, w = g.off, g.tgt, g.fw
    while dq:
        u = dq.popleft()
        inq[u] = 0
        pu = p[u]
        for i in range(off[u], off[u + 1]):
            v = tgt[i]
            nd = pu + w[i]
            if nd < p[v] - 1e-9:
                p[v] = nd
                if not inq[v]:
                    cnt[v] += 1
                    if cnt[v] > g.n:
                        raise RuntimeError("negative cycle")
                    inq[v] = 1
                    dq.append(v)
    with open(CACHE / f"{weighting}_p.bin", "wb") as f:
        p.tofile(f)
    return time.perf_counter() - t0


def job_landmark(weighting, i, L, reverse):
    low_priority()
    g = graph(weighting)
    t0 = time.perf_counter()
    p = _mapped(f"{weighting}_p.bin")
    off, nbr, w = (g.roff, g.rsrc, g.rw) if reverse else (g.off, g.tgt, g.fw)
    d = array("d", [INF]) * g.n
    d[L] = 0.0
    pq = [(0.0, L)]
    pop, push = heapq.heappop, heapq.heappush
    while pq:
        du, u = pop(pq)
        if du > d[u]:
            continue
        pu = p[u]
        for k in range(off[u], off[u + 1]):
            v = nbr[k]
            rw = w[k] + (p[v] - pu if reverse else pu - p[v])
            nd = du + (rw if rw > 0 else 0.0)
            if nd < d[v]:
                d[v] = nd
                push(pq, (nd, v))
    pL = p[L]
    for v in range(g.n):
        if d[v] < INF:
            d[v] = d[v] - p[v] + pL if reverse else d[v] - pL + p[v]
    with open(CACHE / landmark_file(weighting, i, reverse), "wb") as f:
        d.tofile(f)
    return time.perf_counter() - t0


_pots = {}


def job_query(weighting, method, s, t):
    """method: johnson | alt | domain | max | bidir-<kind> | naive"""
    low_priority()
    g = graph(weighting)
    kind = method.split("-")[-1]
    if method != "naive" and (weighting, kind) not in _pots:
        _pots[(weighting, kind)] = Potential(g, kind)
    t0 = time.perf_counter()
    if method == "naive":
        xm, ym, S = g.xm, g.ym, g.scale
        tx, ty = xm[t], ym[t]
        h = lambda v: S * math.hypot(xm[v] - tx, ym[v] - ty)
        dist, exp = astar(g, s, t, h, max_expanded=20_000_000, reopen=True)
    elif method.startswith("bidir"):
        dist, exp = bidirectional(g, s, t, _pots[(weighting, kind)])
    else:
        dist, exp = astar(g, s, t, _pots[(weighting, kind)].to_target(t))
    return dist, exp, time.perf_counter() - t0


def job_bf(weighting, s, t, time_cap):
    low_priority()
    g = graph(weighting)
    t0 = time.perf_counter()
    dist, exp = spfa_query(g, s, t, time_cap)
    return dist, exp, time.perf_counter() - t0


# ---------------------------------------------------------------- 2 cores

def _par_worker(direction, weighting, kind, shared, tasks, results):
    low_priority()
    g = graph(weighting)
    pot = Potential(g, kind)
    fwd = direction == 0
    off, nbr, w = (g.off, g.tgt, g.fw) if fwd else (g.roff, g.rsrc, g.rw)
    my_d, my_st = (shared["df"], shared["sf"]) if fwd else (shared["dr"], shared["sr"])
    ot_d, ot_st = (shared["dr"], shared["sr"]) if fwd else (shared["df"], shared["sf"])
    my_top, ot_top = (shared["tf"], shared["tr"]) if fwd else (shared["tr"], shared["tf"])
    mu, lock, stop, exh = shared["mu"], shared["lock"], shared["stop"], shared["exh"]
    sign = 1.0 if fwd else -1.0
    closed = bytearray(g.n)
    while True:
        task = tasks.get()
        if task is None:
            return
        qid, s, t = task
        pt, ps = pot.to_target(t), pot.from_source(s)
        key = lambda v: sign * 0.5 * (pt(v) - ps(v))
        touched = []
        start = s if fwd else t
        my_d[start] = 0.0
        my_st[start] = qid
        pq = [(key(start), start)]
        expanded = 0
        t0 = time.perf_counter()
        while not stop.value:
            if not pq:
                exh.value = direction + 1
                stop.value = 1
                break
            k, u = pq[0]
            my_top.value = k
            if k + ot_top.value >= mu.value:
                stop.value = 1
                break
            heapq.heappop(pq)
            if closed[u]:
                continue
            closed[u] = 1
            touched.append(u)
            expanded += 1
            du = my_d[u]
            for i in range(off[u], off[u + 1]):
                v = nbr[i]
                nd = du + w[i]
                if my_st[v] != qid or nd < my_d[v]:
                    my_d[v] = nd
                    my_st[v] = qid
                    heapq.heappush(pq, (nd + key(v), v))
                    if ot_st[v] == qid:
                        c = nd + ot_d[v]
                        if c < mu.value:
                            with lock:
                                if c < mu.value:
                                    mu.value = c
        for u in touched:
            closed[u] = 0
        results.put((direction, expanded, time.perf_counter() - t0))


class ParallelCSR:
    def __init__(self, weighting, kind, n):
        import multiprocessing as mp
        ctx = mp.get_context("spawn")
        self.sh = {
            "df": ctx.RawArray("d", n), "dr": ctx.RawArray("d", n),
            "sf": ctx.RawArray("i", n), "sr": ctx.RawArray("i", n),
            "tf": ctx.RawValue("d", 0.0), "tr": ctx.RawValue("d", 0.0),
            "mu": ctx.RawValue("d", INF), "lock": ctx.Lock(),
            "stop": ctx.RawValue("i", 0), "exh": ctx.RawValue("i", 0),
        }
        self.tasks = [ctx.Queue(), ctx.Queue()]
        self.results = ctx.Queue()
        self.procs = [ctx.Process(target=_par_worker, daemon=True,
                                  args=(i, weighting, kind, self.sh, self.tasks[i], self.results))
                      for i in range(2)]
        for p in self.procs:
            p.start()
        self.pot = Potential(graph(weighting), kind)
        self.qid = 0

    def query(self, s, t):
        sh = self.sh
        self.qid += 1
        pt, ps = self.pot.to_target(t), self.pot.from_source(s)
        sh["tf"].value = 0.5 * (pt(s) - ps(s))
        sh["tr"].value = -0.5 * (pt(t) - ps(t))
        sh["mu"].value = INF
        sh["stop"].value = 0
        sh["exh"].value = 0
        t0 = time.perf_counter()
        for q in self.tasks:
            q.put((self.qid, s, t))
        res = [self.results.get() for _ in range(2)]
        elapsed = time.perf_counter() - t0
        exh = sh["exh"].value
        if exh == 1:
            dist = sh["df"][t] if sh["sf"][t] == self.qid else INF
        elif exh == 2:
            dist = sh["dr"][s] if sh["sr"][s] == self.qid else INF
        else:
            dist = sh["mu"].value
        return dist, sum(r[1] for r in res), elapsed

    def close(self):
        for q in self.tasks:
            q.put(None)
        for p in self.procs:
            p.join(timeout=10)
