"""Customizable Contraction Hierarchies (CCH) with NEGATIVE edge weights,
without any Bellman-Ford / Johnson potential.

Why this can work
-----------------
CCH splits route planning into three phases:
  1. metric-independent: a nested-dissection node order and the chordal
     "contracted" graph (depends only on the road topology);
  2. customization: for every contracted arc {y,z}, the cost of the best
     y->z and z->y path whose inner nodes are all ranked below y and z.
     It is computed by processing "lower triangles" {x,y,z}, x < y,z, in
     increasing order of x:  up(y,z) = min(up(y,z), cost(y->x) + cost(x->z)).
     This is vertex elimination in the (min,+) semiring. It only needs
     "no negative cycles", not non-negative weights;
  3. query: scan ALL ancestors of s (and of t) in the elimination tree in
     increasing rank order. The upward graph is a DAG processed in
     topological order, so no Dijkstra ordering (and no non-negativity)
     is needed.

Correctness sketch (no negative cycles): take a shortest s-t path and
repeatedly remove its lowest-ranked inner node x, replacing the two arcs by
the contracted arc between x's neighbours. Customization guarantees that
arc is no more expensive, so the result is an up-down path of cost at most
the optimum. Every customized value is the cost of a real walk, and walks
cannot be cheaper than the shortest path without a negative cycle.

Negative cycles: after customization, an arc with up + down < 0 is a
negative cycle (its two endpoints are the two highest-ranked nodes of some
negative closed walk). The converse holds as well, so this is an exact test.
"""
import heapq

import inertial_flow
import time
from array import array

INF = float("inf")


# --------------------------------------------------------- 1. node ordering

def nested_dissection_order(n, nbrs, xy, leaf=64, method="flow"):
    """Recursive geometric bisection. nbrs[v] = iterable of undirected neighbours.
    Returns rank[v] (0 = lowest)."""
    rank = array("i", [-1]) * n
    mark = array("i", [0]) * n   # which side of the current split a node is on
    stamp = [0]
    stack = [(list(range(n)), 0)]  # (nodes, first free rank)
    while stack:
        nodes, lo = stack.pop()
        if len(nodes) <= leaf:
            # small cell: order by degree inside the cell (cheap min-degree)
            nodes.sort(key=lambda v: len(nbrs[v]))
            for i, v in enumerate(nodes):
                rank[v] = lo + i
            continue
        if method == "flow":
            stamp[0] += 1
            for v in nodes:
                mark[v] = stamp[0]
            res = inertial_flow.bisect(nodes, nbrs, xy, mark, stamp[0])
            if res is not None:
                _place(rank, *res, lo, stack)
                continue
        # try 4 cut directions, keep the one with the smallest separator
        best = None
        for key in (lambda v: xy[v][0], lambda v: xy[v][1],
                    lambda v: xy[v][0] + xy[v][1], lambda v: xy[v][0] - xy[v][1]):
            order = sorted(nodes, key=key)
            half = len(order) // 2
            A, B = order[:half], order[half:]
            stamp[0] += 2
            sa, sb = stamp[0], stamp[0] + 1
            for v in A:
                mark[v] = sa
            for v in B:
                mark[v] = sb
            # vertex separator: boundary nodes of the smaller-boundary side
            sepA = [v for v in A if any(mark[u] == sb for u in nbrs[v])]
            sepB = [v for v in B if any(mark[u] == sa for u in nbrs[v])]
            if len(sepB) < len(sepA):
                A, B, sepA, sa = B, A, sepB, sb
            if best is None or len(sepA) < len(best[2]):
                best = (A, B, sepA, sa)
        A, B, sep, _ = best
        stamp[0] += 1
        sa = stamp[0]
        for v in A:
            mark[v] = sa
        for v in sep:
            mark[v] = 0
        A = [v for v in A if mark[v] == sa]
        hi = lo + len(nodes)
        for i, v in enumerate(sep):
            rank[v] = hi - len(sep) + i
        stack.append((A, lo))
        stack.append((B, lo + len(A)))
    return rank


def _place(rank, A, B, sep, lo, stack):
    hi = lo + len(A) + len(B) + len(sep)
    for i, v in enumerate(sep):
        rank[v] = hi - len(sep) + i
    stack.append((A, lo))
    stack.append((B, lo + len(A)))


# --------------------------------------------------------- 2. contraction

class CCH:
    def __init__(self, g, xy=None, leaf=64, method="flow", log=print):
        """g: graphs.Graph (adj lists). xy: coordinates for the ordering."""
        t0 = time.perf_counter()
        n = g.n
        self.n = n
        und = [set() for _ in range(n)]
        for u in range(n):
            for v, _ in g.adj[u]:
                if u != v:
                    und[u].add(v)
                    und[v].add(u)
        xy = xy or g.coords
        rank = nested_dissection_order(n, und, xy, leaf, method)
        self.rank = rank
        self.t_order = time.perf_counter() - t0

        # chordal completion in rank space, via the elimination-tree trick:
        # when x is eliminated, its upper neighbours become neighbours of the
        # lowest of them (which is x's parent).
        t0 = time.perf_counter()
        up = [set() for _ in range(n)]
        for u in range(n):
            ru = rank[u]
            for v in und[u]:
                rv = rank[v]
                if rv > ru:
                    up[ru].add(rv)
        del und
        parent = array("i", [-1]) * n
        for x in range(n):
            U = up[x]
            if U:
                p = min(U)
                parent[x] = p
                if len(U) > 1:
                    up[p].update(U)
                    up[p].discard(p)
        self.parent = parent

        # upward arcs in CSR form, sorted by head rank
        first = array("i", [0]) * (n + 1)
        heads = array("i")
        for x in range(n):
            hs = sorted(up[x])
            heads.extend(hs)
            first[x + 1] = len(heads)
        del up
        self.first, self.head = first, heads
        self.m = len(heads)
        self.t_contract = time.perf_counter() - t0

        # lower triangles {x < y < z}: arcs a=(x,y), b=(x,z), c=(y,z)
        t0 = time.perf_counter()
        tri_a, tri_b, tri_c = array("i"), array("i"), array("i")
        for x in range(n):
            s, e = first[x], first[x + 1]
            for i in range(s, e):
                y = heads[i]
                # arcs of y, to find (y, z) by merging two sorted lists
                j, je = first[y], first[y + 1]
                for k in range(i + 1, e):
                    z = heads[k]
                    while j < je and heads[j] < z:
                        j += 1
                    # chordality guarantees (y, z) exists
                    tri_a.append(i)
                    tri_b.append(k)
                    tri_c.append(j)
        self.tri_a, self.tri_b, self.tri_c = tri_a, tri_b, tri_c
        self.t_triangles = time.perf_counter() - t0

        # original edges mapped onto arcs (for customization)
        t0 = time.perf_counter()
        self.edge_arc = array("i")   # arc id per original edge (in g.adj order)
        self.edge_dir = bytearray()  # 1 = upward (lower -> higher rank)
        for u in range(n):
            ru = rank[u]
            for v, _ in g.adj[u]:
                rv = rank[v]
                if ru == rv:
                    self.edge_arc.append(-1)
                    self.edge_dir.append(0)
                    continue
                lo_, hi_ = (ru, rv) if ru < rv else (rv, ru)
                self.edge_arc.append(self._find_arc(lo_, hi_))
                self.edge_dir.append(1 if ru < rv else 0)
        self.t_map = time.perf_counter() - t0
        self.depth = self._max_depth()
        log(f"  CCH: n={n:,} arcs={self.m:,} triangles={len(tri_a):,} "
            f"etree depth={self.depth} | order {self.t_order:.1f}s, "
            f"contract {self.t_contract:.1f}s, triangles {self.t_triangles:.1f}s")

    def _find_arc(self, lo, hi):
        h, s, e = self.head, self.first[lo], self.first[lo + 1]
        while s < e:  # binary search
            mid = (s + e) // 2
            if h[mid] < hi:
                s = mid + 1
            else:
                e = mid
        assert h[s] == hi
        return s

    def _max_depth(self):
        depth = array("i", [0]) * self.n
        best = 0
        for x in range(self.n - 1, -1, -1):
            p = self.parent[x]
            if p >= 0:
                depth[x] = depth[p] + 1
                if depth[x] > best:
                    best = depth[x]
        return best

    # ----------------------------------------------------- 3. customization

    def base_weights(self):
        """Arc weights from the original edges only (min over parallel edges)."""
        up = array("d", [INF]) * self.m
        down = array("d", [INF]) * self.m
        ea, ed = self.edge_arc, self.edge_dir
        for k, w in enumerate(self.edge_w):
            a = ea[k]
            if a >= 0:
                if ed[k]:
                    if w < up[a]:
                        up[a] = w
                elif w < down[a]:
                    down[a] = w
        return up, down

    def customize(self, g):
        """Full customization from the weights currently in g.
        Returns the list of arcs that witness a negative cycle (empty if none)."""
        self.edge_w = array("d", (w for u in range(self.n) for _, w in g.adj[u]))
        self.base_up, self.base_down = self.base_weights()
        up, down = array("d", self.base_up), array("d", self.base_down)
        ta, tb, tc = self.tri_a, self.tri_b, self.tri_c
        for i in range(len(ta)):
            a, b, c = ta[i], tb[i], tc[i]
            # y -> x -> z   and   z -> x -> y
            v = down[a] + up[b]
            if v < up[c]:
                up[c] = v
            v = down[b] + up[a]
            if v < down[c]:
                down[c] = v
        self.up, self.down = up, down
        return [a for a in range(self.m) if up[a] + down[a] < -1e-9]

    def build_update_index(self):
        """Triangle indices for partial re-customization (built once, lazily)."""
        m = self.m
        ta, tb, tc = self.tri_a, self.tri_b, self.tri_c
        T = len(ta)

        def group(key):
            off = array("i", [0]) * (m + 1)
            for i in range(T):
                off[key[i] + 1] += 1
            for i in range(m):
                off[i + 1] += off[i]
            pos = array("i", off[:-1])
            ids = array("i", [0]) * T
            for i in range(T):
                k = key[i]
                ids[pos[k]] = i
                pos[k] += 1
            return off, ids
        self.by_c = group(tc)                     # lower triangles of an arc
        # triangles in which an arc is a lower side (a or b)
        both = array("i", ta)
        both.extend(tb)
        tid = array("i", range(T))
        tid.extend(range(T))
        off = array("i", [0]) * (m + 1)
        for k in both:
            off[k + 1] += 1
        for i in range(m):
            off[i + 1] += off[i]
        pos = array("i", off[:-1])
        ids = array("i", [0]) * len(both)
        for k, t in zip(both, tid):
            ids[pos[k]] = t
            pos[k] += 1
        self.as_side = (off, ids)
        # rank of the lower endpoint of each arc (processing priority)
        tail = array("i", [0]) * m
        for x in range(self.n):
            for a in range(self.first[x], self.first[x + 1]):
                tail[a] = x
        self.tail = tail


        # original edges grouped by arc, to recompute base values quickly
        ea = self.edge_arc
        off = array("i", [0]) * (m + 1)
        for a in ea:
            if a >= 0:
                off[a + 1] += 1
        for i in range(m):
            off[i + 1] += off[i]
        pos = array("i", off[:-1])
        ids = array("i", [0]) * off[m]
        for k, a in enumerate(ea):
            if a >= 0:
                ids[pos[a]] = k
                pos[a] += 1
        self.arc_edges = (off, ids)

    def update_edges(self, changes):
        """Partial re-customization. changes = [(edge_index, new_weight)], with
        edge indices in g.adj order. Only arcs whose value can change are
        re-evaluated, lowest first. Returns the number of arcs re-evaluated."""
        if not hasattr(self, "by_c"):
            self.build_update_index()
        up, down, bu, bd = self.up, self.down, self.base_up, self.base_down
        ew, ea, ed = self.edge_w, self.edge_arc, self.edge_dir
        aoff, aids = self.arc_edges
        touched = set()
        for k, w in changes:
            ew[k] = w
            if ea[k] >= 0:
                touched.add(ea[k])
        for a in touched:
            nu = nd = INF
            for j in range(aoff[a], aoff[a + 1]):
                k = aids[j]
                if ed[k]:
                    nu = min(nu, ew[k])
                else:
                    nd = min(nd, ew[k])
            bu[a], bd[a] = nu, nd
        ta, tb, tc = self.tri_a, self.tri_b, self.tri_c
        coff, cids = self.by_c
        soff, sids = self.as_side
        tail = self.tail
        heap = [(tail[a], a) for a in touched]
        heapq.heapify(heap)
        queued = set(touched)
        evaluated = 0
        self.last_changed = changed = []
        while heap:
            _, c = heapq.heappop(heap)
            evaluated += 1
            nu, nd = bu[c], bd[c]
            for j in range(coff[c], coff[c + 1]):
                t = cids[j]
                a, b = ta[t], tb[t]
                v = down[a] + up[b]
                if v < nu:
                    nu = v
                v = down[b] + up[a]
                if v < nd:
                    nd = v
            if nu != up[c] or nd != down[c]:
                up[c], down[c] = nu, nd
                changed.append(c)
                for j in range(soff[c], soff[c + 1]):
                    c2 = tc[sids[j]]
                    if c2 not in queued:
                        queued.add(c2)
                        heapq.heappush(heap, (tail[c2], c2))
        return evaluated

    def has_negative_cycle(self):
        up, down = self.up, self.down
        return any(up[a] + down[a] < -1e-9 for a in range(self.m))

    def update_created_negative_cycle(self):
        """Exact test after update_edges(), valid when the metric before the
        update was conservative: an arc whose value did not change still has
        up + down >= 0, so only the changed arcs need checking."""
        up, down = self.up, self.down
        return any(up[a] + down[a] < -1e-9 for a in self.last_changed)

    # ----------------------------------------------------------- 4. query

    def query(self, s, t):
        """Exact s-t distance (original node ids). Returns (dist, arcs_scanned).
        Arcs out of an ancestor of s only lead to ancestors of s, so both
        searches stay on their elimination-tree paths."""
        if not hasattr(self, "_df"):
            n = self.n
            self._df = array("d", [INF]) * n
            self._dr = array("d", [INF]) * n
            self._mark = array("i", [0]) * n
            self._qid = 0
        self._qid += 1
        q = self._qid
        df, dr, mark = self._df, self._dr, self._mark
        parent, first, head = self.parent, self.first, self.head
        up, down = self.up, self.down
        rs, rt = self.rank[s], self.rank[t]
        x = rs
        while x >= 0:
            df[x] = INF
            mark[x] = q
            x = parent[x]
        x = rt
        while x >= 0:
            dr[x] = INF
            x = parent[x]
        df[rs] = 0.0
        dr[rt] = 0.0
        scanned = 0
        x = rs
        while x >= 0:
            dx = df[x]
            if dx < INF:
                e = first[x + 1]
                scanned += e - first[x]
                for a in range(first[x], e):
                    v = dx + up[a]
                    y = head[a]
                    if v < df[y]:
                        df[y] = v
            x = parent[x]
        best = INF
        x = rt
        while x >= 0:
            dx = dr[x]
            if dx < INF:
                if mark[x] == q:
                    v = dx + df[x]
                    if v < best:
                        best = v
                e = first[x + 1]
                scanned += e - first[x]
                for a in range(first[x], e):
                    v = dx + down[a]
                    y = head[a]
                    if v < dr[y]:
                        dr[y] = v
            x = parent[x]
        return best, scanned
