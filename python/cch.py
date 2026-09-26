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
    def __init__(self, g, xy=None, leaf=64, method="flow", log=print, order=None):
        """g: graphs.Graph (adj lists). xy: coordinates for the ordering.
        order: optional explicit rank[v] (0 = lowest), e.g. for hand-built examples."""
        t0 = time.perf_counter()
        n = g.n
        self.n = n
        und = [set() for _ in range(n)]
        for u in range(n):
            for v, _ in g.adj[u]:
                if u != v:
                    und[u].add(v)
                    und[v].add(u)
        if order is not None:
            assert sorted(order) == list(range(n)), "order must be a permutation"
            rank = array("i", order)
        else:
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
        self.loop_edges = [k for k in range(len(self.edge_arc)) if self.edge_arc[k] < 0]
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

    def customize(self, g, stop_on_negative_cycle=False):
        """Full customization from the weights currently in g.
        Returns the list of arcs that witness a negative cycle (empty if none);
        [-1] means a negative self-loop in the input (loops do not enter G+).

        With stop_on_negative_cycle, elimination stops as soon as a cycle becomes
        apparent: when vertex x is eliminated its arcs are final, so a negative
        2-cycle on one of them is a negative cycle of the input."""
        self.edge_w = array("d", (w for u in range(self.n) for _, w in g.adj[u]))
        self.base_up, self.base_down = self.base_weights()
        # a negative self-loop is a negative cycle, but loops never enter G+
        if any(self.edge_w[k] < 0 for k in self.loop_edges):
            self.up, self.down = array("d", self.base_up), array("d", self.base_down)
            self.cycle_arc = -1
            return [-1]
        up, down = array("d", self.base_up), array("d", self.base_down)
        ta, tb, tc = self.tri_a, self.tri_b, self.tri_c
        first = self.first
        x = 0
        for i in range(len(ta)):
            a, b, c = ta[i], tb[i], tc[i]
            if stop_on_negative_cycle:
                # arcs of every vertex below the current one are final by now
                while x < self.n and first[x + 1] <= a:
                    for aa in range(first[x], first[x + 1]):
                        if up[aa] + down[aa] < -1e-9:
                            self.up, self.down = up, down
                            self.cycle_arc = aa
                            return [aa]
                    x += 1
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

    # ------------------------------------------------ path unpacking

    def query_path(self, s, t):
        """Distance plus the shortest path as a list of original vertices.
        Unpacks each shortcut through the lower triangle that realises its value."""
        if not hasattr(self, "_df"):
            self.query(s, t)
        parent, first, head = self.parent, self.first, self.head
        up, down = self.up, self.down
        rs, rt = self.rank[s], self.rank[t]
        n = self.n
        df = [INF] * n
        dr = [INF] * n
        pf = [-1] * n          # arc used to reach this vertex going up from s
        pr = [-1] * n
        df[rs] = 0.0
        dr[rt] = 0.0
        x = rs
        while x >= 0:
            dx = df[x]
            if dx < INF:
                for a in range(first[x], first[x + 1]):
                    y = head[a]
                    if dx + up[a] < df[y]:
                        df[y] = dx + up[a]
                        pf[y] = a
            x = parent[x]
        x = rt
        while x >= 0:
            dx = dr[x]
            if dx < INF:
                for a in range(first[x], first[x + 1]):
                    y = head[a]
                    if dx + down[a] < dr[y]:
                        dr[y] = dx + down[a]
                        pr[y] = a
            x = parent[x]
        best, peak = INF, -1
        x = rt
        anc = set()
        while x >= 0:
            anc.add(x)
            x = parent[x]
        x = rs
        while x >= 0:
            if x in anc and df[x] + dr[x] < best:
                best, peak = df[x] + dr[x], x
            x = parent[x]
        if peak < 0:
            return INF, []
        # shortcut arcs of the up-down path, as (from_rank, to_rank) pairs
        legs = []
        x = peak
        while x != rs:
            a = pf[x]
            v = self._arc_tail(a)
            legs.append((v, x))
            x = v
        legs.reverse()
        x = peak
        while x != rt:
            a = pr[x]
            v = self._arc_tail(a)
            legs.append((x, v))
            x = v
        path = [legs[0][0]] if legs else [rs]
        for v, w in legs:
            self._unpack(v, w, path)
        inv = [0] * n
        for u in range(n):
            inv[self.rank[u]] = u
        return best, [inv[r] for r in path]

    def _arc_tail(self, a):
        if not hasattr(self, "tail"):
            self.build_update_index()
        return self.tail[a]

    def _unpack(self, v, w, out):
        """Append the vertices of the shortest v->w path (excluding v) to out."""
        stack = [(v, w)]
        while stack:
            v, w = stack.pop()
            lo, hi = (v, w) if v < w else (w, v)
            a = self._find_arc(lo, hi)
            cost = self.up[a] if v < w else self.down[a]
            mid = -1
            for u in self._lower_triangle_vertices(a):
                c1 = self.down[self._find_arc(u, v)] if u < v else self.up[self._find_arc(v, u)]
                c2 = self.up[self._find_arc(u, w)] if u < w else self.down[self._find_arc(w, u)]
                if abs(c1 + c2 - cost) <= 1e-9 * max(1.0, abs(cost)):
                    mid = u
                    break
            if mid < 0:
                out.append(w)
            else:
                stack.append((mid, w))
                stack.append((v, mid))

    def _lower_triangle_vertices(self, a):
        if not hasattr(self, "by_c"):
            self.build_update_index()
        coff, cids = self.by_c
        return [self.tail[self.tri_a[cids[j]]] for j in range(coff[a], coff[a + 1])]

    # ------------------------------------------------ acceleration variants

    def perfect_customize(self):
        """Top-down pass over upper/intermediate triangles (CCH survey, step 3).
        Afterwards every arc should hold the true distance between its endpoints.
        Valid with negative weights? -> tested in test_cch_pruning.py."""
        up, down = self.up, self.down
        ta, tb, tc = self.tri_a, self.tri_b, self.tri_c
        tail = getattr(self, "tail", None)
        if tail is None:
            self.build_update_index()
            tail = self.tail
        order = sorted(range(len(ta)), key=lambda i: -tail[ta[i]])
        for i in order:
            a, b, c = ta[i], tb[i], tc[i]      # a=(x,y), b=(x,z), c=(y,z), x lowest
            v = up[b] + down[c]                # x -> z -> y
            if v < up[a]:
                up[a] = v
            v = up[c] + down[b]                # y -> z -> x
            if v < down[a]:
                down[a] = v
            v = up[a] + up[c]                  # x -> y -> z
            if v < up[b]:
                up[b] = v
            v = down[c] + down[a]              # z -> y -> x
            if v < down[b]:
                down[b] = v
        return up, down

    def prunable_arcs(self):
        """Arcs a witness search may remove after perfect customization: the same
        distance is realised through a higher-ranked vertex (upper or
        intermediate triangle), so an up-down path still exists.

        Exact when every cycle has positive length. With a zero-length cycle two
        arcs can witness each other's removal (also for non-negative weights);
        break ties lexicographically, e.g. by integer weights w' = (n+1)w + 1."""
        up, down = self.up, self.down
        ta, tb, tc = self.tri_a, self.tri_b, self.tri_c
        drop_up = bytearray(self.m)
        drop_dn = bytearray(self.m)
        eps = 1e-9
        for i in range(len(ta)):
            a, b, c = ta[i], tb[i], tc[i]
            if abs(up[a] - (up[b] + down[c])) < eps:
                drop_up[a] = 1
            if abs(down[a] - (up[c] + down[b])) < eps:
                drop_dn[a] = 1
            if abs(up[b] - (up[a] + up[c])) < eps:
                drop_up[b] = 1
            if abs(down[b] - (down[c] + down[a])) < eps:
                drop_dn[b] = 1
            # NB: arc c = (y, z) must not be dropped here. Its detour runs through
            # the lower vertex x, which an up-down (elimination-tree) search cannot
            # use, so removing it would break correctness.
        return drop_up, drop_dn

    def query_pruned(self, s, t, drop_up, drop_dn):
        """Elimination-tree query that ignores the arcs marked as removable."""
        up, down = self.up, self.down
        saved_up = [up[a] for a in range(self.m) if drop_up[a]]
        saved_dn = [down[a] for a in range(self.m) if drop_dn[a]]
        for a in range(self.m):
            if drop_up[a]:
                up[a] = INF
            if drop_dn[a]:
                down[a] = INF
        try:
            return self.query(s, t)[0]
        finally:
            i = j = 0
            for a in range(self.m):
                if drop_up[a]:
                    up[a] = saved_up[i]
                    i += 1
                if drop_dn[a]:
                    down[a] = saved_dn[j]
                    j += 1

    def update_edges_tiebased(self, changes):
        """Partial update in the style of Dibbelt et al. (2016), section 7.7:
        propagate a change only to arcs whose value was *realised* by the old
        value (a tie test), instead of re-evaluating every arc above. A decrease
        can also create a new witness where the old value was not tight, so a
        change is propagated as well when the new candidate improves the target.
        Both tests are needed for non-negative weights too. Comparisons are exact:
        every candidate is the same floating-point sum as in customize()."""
        if not hasattr(self, "by_c"):
            self.build_update_index()
        up, down, bu, bd = self.up, self.down, self.base_up, self.base_down
        ew, ea, ed = self.edge_w, self.edge_arc, self.edge_dir
        aoff, aids = self.arc_edges
        ta, tb, tc = self.tri_a, self.tri_b, self.tri_c
        coff, cids = self.by_c
        soff, sids = self.as_side
        tail = self.tail
        touched = set()
        for k, w in changes:
            ew[k] = w
            if ea[k] >= 0:
                touched.add(ea[k])
        for a in touched:
            nu = nd = INF
            for j in range(aoff[a], aoff[a + 1]):
                kk = aids[j]
                if ed[kk]:
                    nu = min(nu, ew[kk])
                else:
                    nd = min(nd, ew[kk])
            bu[a], bd[a] = nu, nd
        heap = [(tail[a], a) for a in touched]
        heapq.heapify(heap)
        queued = set(touched)
        evaluated = 0
        self.last_changed = changed = []
        while heap:
            _, c = heapq.heappop(heap)
            evaluated += 1
            old_up, old_dn = up[c], down[c]
            # recompute this arc from its own lower triangles and input edges
            nu, nd = bu[c], bd[c]
            for j in range(coff[c], coff[c + 1]):
                t = cids[j]
                a, b = ta[t], tb[t]
                nu = min(nu, down[a] + up[b])
                nd = min(nd, down[b] + up[a])
            up[c], down[c] = nu, nd
            if nu == old_up and nd == old_dn:
                continue
            changed.append(c)
            # propagate only where the OLD value was realised (the tie test)
            for j in range(soff[c], soff[c + 1]):
                t = sids[j]
                a, b, e = ta[t], tb[t], tc[t]
                partner = b if a == c else a
                if a == c:
                    old_cu, new_cu = old_dn + up[partner], down[c] + up[partner]
                    old_cd, new_cd = down[partner] + old_up, down[partner] + up[c]
                else:
                    old_cu, new_cu = down[partner] + old_up, down[partner] + up[c]
                    old_cd, new_cd = old_dn + up[partner], down[c] + up[partner]
                # the old value was this arc's witness (it may have to rise), or
                # the new value improves it (it may have to fall)
                realised_up = up[e] == old_cu or new_cu < up[e]
                realised_dn = down[e] == old_cd or new_cd < down[e]
                if (realised_up or realised_dn) and e not in queued:
                    queued.add(e)
                    heapq.heappush(heap, (tail[e], e))
        return evaluated

    def _negative_loop(self):
        return any(self.edge_w[k] < 0 for k in self.loop_edges)

    def has_negative_cycle(self):
        up, down = self.up, self.down
        return self._negative_loop() or any(up[a] + down[a] < -1e-9 for a in range(self.m))

    def update_created_negative_cycle(self):
        """Exact test after update_edges(), valid when the metric before the
        update was conservative: an arc whose value did not change still has
        up + down >= 0, so only the changed arcs need checking."""
        up, down = self.up, self.down
        return self._negative_loop() or any(up[a] + down[a] < -1e-9 for a in self.last_changed)

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

    def potential(self):
        """A feasible potential read off the customized CCH (conservative metric):
        p(v) = min(0, min_u dist(u, v)), the virtual-source distances that
        Johnson's algorithm computes with Bellman-Ford. Every shortest u-v path
        has an up-down representative in G+, so an upward pass over all vertices
        in increasing rank and a downward pass in decreasing rank suffice: O(|E+|)."""
        first, head, up, down = self.first, self.head, self.up, self.down
        p = [0.0] * self.n                 # indexed by rank
        for x in range(self.n):
            px = p[x]
            for a in range(first[x], first[x + 1]):
                v = px + up[a]
                if v < p[head[a]]:
                    p[head[a]] = v
        for x in range(self.n - 1, -1, -1):
            best = p[x]
            for a in range(first[x], first[x + 1]):
                v = p[head[a]] + down[a]
                if v < best:
                    best = v
            p[x] = best
        return [p[self.rank[v]] for v in range(self.n)]
