"""Two-process bidirectional A*: one core searches s->t, the other t->s.

The two processes share, through multiprocessing shared memory:
  * their distance labels (with a per-query stamp, so no O(n) reset),
  * the best meeting cost mu (lock-protected; it only ever decreases),
  * their current smallest queue key (the stopping test is top_f + top_r >= mu).

Each side only ever reads the other's values. A stale read of the other side's
top key is smaller than its true value, and a stale mu is larger than the true
one, so a stale read can only make a side stop later, never too early.
When one side runs out of nodes it has settled everything reachable, so its
label at the far endpoint is exact and is returned directly.
"""
import heapq
import multiprocessing as mp

INF = float("inf")


def _worker(direction, g, potential, shared, tasks, results):
    fwd = direction == 0
    adj = g.adj if fwd else g.radj
    my_d, my_st = (shared["df"], shared["sf"]) if fwd else (shared["dr"], shared["sr"])
    ot_d, ot_st = (shared["dr"], shared["sr"]) if fwd else (shared["df"], shared["sf"])
    my_top, ot_top = (shared["tf"], shared["tr"]) if fwd else (shared["tr"], shared["tf"])
    mu, lock = shared["mu"], shared["lock"]
    stop, exhausted = shared["stop"], shared["exh"]
    sign = 1.0 if fwd else -1.0

    while True:
        task = tasks.get()
        if task is None:
            return
        qid, s, t = task
        pi_t, pi_s = potential.to_target(t), potential.from_source(s)

        def key_pot(v):
            return sign * 0.5 * (pi_t(v) - pi_s(v))

        start = s if fwd else t
        closed = set()
        my_d[start] = 0.0
        my_st[start] = qid
        pq = [(key_pot(start), start)]
        expanded = 0
        while True:
            if stop.value:
                break
            if not pq:
                exhausted.value = direction + 1
                stop.value = 1
                break
            k, u = pq[0]
            my_top.value = k
            if k + ot_top.value >= mu.value:
                stop.value = 1
                break
            heapq.heappop(pq)
            if u in closed:
                continue
            closed.add(u)
            expanded += 1
            du = my_d[u]
            for v, w in adj[u]:
                nd = du + w
                if my_st[v] != qid or nd < my_d[v]:
                    my_d[v] = nd
                    my_st[v] = qid
                    heapq.heappush(pq, (nd + key_pot(v), v))
                    if ot_st[v] == qid:
                        cand = nd + ot_d[v]
                        if cand < mu.value:
                            with lock:
                                if cand < mu.value:
                                    mu.value = cand
        results.put((direction, expanded))


class ParallelBidirectional:
    def __init__(self, g, potential):
        ctx = mp.get_context("spawn")
        n = g.n
        self.shared = {
            "df": ctx.RawArray("d", n), "dr": ctx.RawArray("d", n),
            "sf": ctx.RawArray("i", n), "sr": ctx.RawArray("i", n),
            "tf": ctx.RawValue("d", 0.0), "tr": ctx.RawValue("d", 0.0),
            "mu": ctx.RawValue("d", INF), "lock": ctx.Lock(), "stop": ctx.RawValue("i", 0),
            "exh": ctx.RawValue("i", 0),
        }
        self.tasks = [ctx.Queue(), ctx.Queue()]
        self.results = ctx.Queue()
        self.procs = [ctx.Process(target=_worker,
                                  args=(i, g, potential, self.shared, self.tasks[i], self.results),
                                  daemon=True) for i in range(2)]
        for p in self.procs:
            p.start()
        self.potential = potential
        self.qid = 0

    def query(self, s, t):
        if s == t:
            return 0.0, 0
        sh = self.shared
        self.qid += 1
        pi_t, pi_s = self.potential.to_target(t), self.potential.from_source(s)
        # initial top keys are valid lower bounds on each side's future keys
        sh["tf"].value = 0.5 * (pi_t(s) - pi_s(s))
        sh["tr"].value = -0.5 * (pi_t(t) - pi_s(t))
        sh["mu"].value = INF
        sh["stop"].value = 0
        sh["exh"].value = 0
        for q in self.tasks:
            q.put((self.qid, s, t))
        expanded = sum(self.results.get()[1] for _ in range(2))
        exh = sh["exh"].value
        if exh == 1:  # forward side settled everything reachable from s
            dist = sh["df"][t] if sh["sf"][t] == self.qid else INF
        elif exh == 2:
            dist = sh["dr"][s] if sh["sr"][s] == self.qid else INF
        else:
            dist = sh["mu"].value
        return dist, expanded

    def close(self):
        for q in self.tasks:
            q.put(None)
        for p in self.procs:
            p.join(timeout=5)
