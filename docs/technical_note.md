# Potential-Free Customizable Contraction Hierarchies for Road Networks with Negative Edge Weights

**Working draft, shared for feedback.** Whether this observation is already known in the CCH community has not yet been confirmed. Corrections and pointers to prior work are very welcome (please open an issue).

*Author: Aarav Gupta*

---

## Abstract

Shortest-path problems with negative edge weights arise in road networks, most prominently in energy-optimal routing for electric vehicles (EVs), where recuperation on downhill segments yields negative energy costs. Existing speed-up techniques for this setting first transform the weights with a feasible *potential* (Johnson's reweighting, or a height-induced potential) so that Dijkstra-based preprocessing and queries become applicable. We point out that *Customizable Contraction Hierarchies* (CCH), whose literature assumes non-negative weights, remain exact for arbitrary *conservative* weights (weights without negative cycles; if the metric is not conservative, shortest paths are undefined and no potential exists either) when the customization is the standard lower-triangle relaxation and queries use the elimination-tree algorithm. No potential and no Bellman–Ford computation are needed. Customization additionally yields an exact negative-cycle test (after an update, only the re-evaluated arcs need checking), and partial re-customization supports traffic-style updates on negative-weight metrics. On the DIMACS New York and San Francisco Bay road networks, with two different negative-weight metrics (10–20% negative arcs), our pure-Python prototype answers all 400 test queries exactly. Queries are 150–620× faster than Johnson + Dijkstra and 30–150× faster than Johnson + ALT A*. A new metric is customized in 0.9–1.6 s, versus 9–14 s for Johnson's Bellman–Ford plus ALT landmarks, and single-arc updates take under a millisecond.

---

## 1 Introduction

Dijkstra's algorithm and virtually all route-planning speed-up techniques built on it (ALT, Contraction Hierarchies, CRP, …) require non-negative arc weights. Several applications violate this. The best known is energy-optimal EV routing: the energy spent on an arc is roughly `α·length + β·Δheight`, and on steep descents regenerative braking makes this negative [Eisner et al. 2011; Baum et al. 2020]. Since every cycle has non-negative total energy (one cannot gain energy by driving in a circle), the weights are *conservative*: no negative cycles.

The standard remedy is **potential shifting**. A potential `p` with `w(u,v) + p(u) − p(v) ≥ 0` is obtained either by Bellman–Ford from a virtual source (Johnson's algorithm) or, for EVs, from vertex heights. The reduced weights are then fed to Dijkstra, A*, or Contraction Hierarchies. This has two costs:

1. Computing a potential for an arbitrary conservative metric needs a Bellman–Ford-type computation, which is expensive on continental graphs and must be redone whenever the metric changes.
2. Domain potentials (height-induced) exist only when the cost model has that structure.

**Customizable Contraction Hierarchies** [Dibbelt, Strasser & Wagner 2016; Bläsius et al. 2025] split preprocessing into a metric-independent phase (a nested-dissection order and a chordal supergraph) and a fast, metric-dependent *customization* phase. Both the original paper and the recent survey define weights as positive or non-negative. The standard CCH query is a bidirectional Dijkstra with a stopping criterion, which is indeed invalid for negative weights.

**Contributions.**

1. We make explicit that basic CCH customization, together with the elimination-tree query, is exact for every conservative metric (Theorem 1, Theorem 2). No potential is required. This is a consequence of classical elimination theory (see below); we give self-contained proofs in CCH terminology.
2. We show that after customization a negative cycle exists **iff** some contracted arc `{v,w}` has `ℓ⁺(v,w) + ℓ⁺(w,v) < 0` (Theorem 3). Negative-cycle detection thus comes for free.
3. We give a partial re-customization procedure that recomputes affected arcs from their lower triangles, and prove it yields exactly the fully customized metric (Theorem 4).
4. We evaluate a prototype on DIMACS road graphs with two negative-weight metrics. We compare against Johnson + Dijkstra, Johnson + ALT A*, bidirectional ALT, and domain-potential A*, including a full-USA study (24M vertices) of the potential-based approaches.

**Relation to classical work (important).** The mathematical core is not new. Algebraic path-problem theory [Carré 1971; Lipton, Rose & Tarjan 1979; Tarjan 1981; Rote 1990] establishes all of the following:

- Shortest paths may have negative arc weights. Rote's running example has arcs of weight −1 and −5 (§2.1).
- The path equations have a solution iff there is no negative cycle (§4).
- Gauss/Gauss–Jordan vertex elimination in any order yields the solution whenever the pivot closures exist, which for shortest paths means no negative cycles (§4.2, Theorems 1–4).
- After eliminating vertices `1…k`, the entry `a_ij` is the best `i–j` path whose intermediate vertices lie in `{1…k}` (§4.3). This is exactly the meaning of a customized CCH shortcut.
- Eliminating a vertex is equivalent to inserting short-cut arcs between its in- and out-neighbours (§4.5, Fig. 3). Nested-dissection orders make this efficient on sparse graphs (§4.5, citing Lipton & Tarjan).
- An LU-type variant satisfies `A* = U*L*` (§4.4), which corresponds to combining an upward and a downward search.
- For shortest paths, the pivot closure reduces to a sign test (§4.4), which is how negative cycles manifest.

CCH customization is precisely this sparse elimination, restricted to the chordal supergraph. What we add is therefore not new mathematics but:

1. the explicit observation that the CCH toolchain inherits this property, which the CCH literature does not state (it assumes weights ≥ 0) and which the negative-weight route-planning literature does not use (it computes potentials first);
2. the CCH-specific consequences: the elimination-tree query as the potential-free query, the 2-cycle form of the negative-cycle test with its update-local variant, and partial re-customization on negative metrics;
3. an implementation and experiments on road networks, including a battery-constrained variant in the (max, ∘) semiring, to which the general theory for ordered semirings (Rote §3–4) applies.

We have not found this application stated in the route-planning literature, but it may be known to experts. Pointers are welcome.

**Scope.** We treat the *unconstrained* shortest-path problem with conservative weights. Energy-optimal EV routing with battery-capacity constraints (a state of charge bounded in `[0, M]`) is not a plain shortest-path problem; Baum et al. handle it with piecewise cost functions. We include a preliminary battery-constrained variant (monotone charge functions under max and composition), validated on small graphs only. Whether its profiles stay small on large networks is open.

---

## 2 Preliminaries

`G = (V, E)` is the undirected graph underlying the road network. A directed weight function is `ℓ : V×V → ℝ ∪ {∞}`, where `ℓ(u,v) = ∞` if `{u,v} ∉ E` or the direction is forbidden. It is **conservative** if every directed cycle has non-negative length. In that case `dist(s,t)` is finite or `∞`, and is attained by a simple path.

**CCH.** A rank order `π` is fixed. Contracting vertices in rank order (completing the upper neighbourhood of each vertex into a clique) gives the chordal supergraph `G⁺ = (V, E⁺)`. For `u ∈ V`, `N↑(u)` denotes its upper neighbourhood, which is a clique in `G⁺`. The **elimination tree** has `parent(u) = min N↑(u)`. Every vertex reachable from `s` in the upward DAG `G↑` is an ancestor of `s` in this tree [Dibbelt et al.]. A **lower triangle** of `{v,w} ∈ E⁺` is a vertex `u` with `u ≺ v`, `u ≺ w`, and `{u,v}, {u,w} ∈ E⁺`.

**Basic customization.** First, `ℓ⁺ ← ℓ` on `E⁺` (shortcuts get `∞`). Then all lower triangles `(u; v, w)` are processed in increasing order of `u`, setting

    ℓ⁺(v,w) ← min(ℓ⁺(v,w), ℓ⁺(v,u) + ℓ⁺(u,w))
    ℓ⁺(w,v) ← min(ℓ⁺(w,v), ℓ⁺(w,u) + ℓ⁺(u,v))

When triangle `(u; v, w)` is processed, the values `ℓ⁺(v,u)`, `ℓ⁺(u,w)`, … are final, because an arc `{u,x}` is changed only by triangles whose lowest vertex is below `u`.

It is the **vertices** that are ordered and eliminated. The triangles are simply grouped by their lowest vertex, and eliminating `u` processes all triangles in which `u` is the lowest (i.e. the middle vertex of `v → u → w`).

**Potentials are irrelevant to the elimination.** Suppose `ℓ` is reweighted with any potential `p`, i.e. `ℓ_p(x,y) = ℓ(x,y) + p(x) − p(y)`. Every customized value is a minimum over lengths of `v–w` walks, and each such length shifts by the same constant `p(v) − p(w)`. So customizing `ℓ_p` yields exactly `ℓ⁺(v,w) + p(v) − p(w)`. Computing a potential first therefore changes nothing but a constant offset per shortcut.

**Negative cycles during elimination.** Elimination itself does not break if negative cycles are present; they simply become apparent. In the algebraic formulation, a pivot closure becomes `−∞`. In CCH terms, some shortcut acquires `ℓ⁺(v,w) + ℓ⁺(w,v) < 0` (Theorem 3). A system can therefore detect a cycle as soon as it appears and reject the offending update.

---

## 3 Correctness for conservative weights

After basic customization, the following hold **for any real weights**, conservative or not:

- **(P0)** Every finite value `ℓ⁺(v,w)` is the length of some `v–w` walk in `G`. This is true by induction: an original arc, or the concatenation of two walks.
- **(P2)** `ℓ⁺ ≤ ℓ` on `E`.
- **(P3)** The lower-triangle inequality `ℓ⁺(v,w) ≤ ℓ⁺(v,u) + ℓ⁺(u,w)` holds for every lower triangle. This follows from the bottom-up processing order alone.

If `ℓ` is conservative, then (P0) gives:

- **(P1)** `ℓ⁺(v,w) ≥ dist(v,w)`.

None of these arguments uses the sign of the weights.

**Theorem 1 (distances are preserved).** Let `ℓ` be conservative. For all `s, t`, the distance in `(G⁺, ℓ⁺)` equals `dist(s,t)`. Moreover, if `dist(s,t) < ∞`, there is an *up-down* `s–t` path in `G⁺` of length `dist(s,t)`.

*Proof.*

1. By (P0) and conservativeness, every `s–t` path in `G⁺` has `ℓ⁺`-length at least `dist(s,t)`. It unpacks to a walk, and walks are no shorter than shortest paths.
2. Conversely, we *construct* an up-down path of length at most `dist(s,t)`. Let `W₀` be a shortest simple `s–t` path in `G`. By (P2), its `ℓ⁺`-length is at most `dist(s,t)`. We build a sequence of `s–t` walks `W₀, W₁, …` in `G⁺`, each with at least one vertex fewer than its predecessor and no greater `ℓ⁺`-length. Each walk is a sequence of consecutive vertices joined by arcs of `E⁺`.
3. Given `Wᵢ`, apply the first applicable rule:
   - **(a) Repeated vertex.** If some vertex occurs twice, delete the closed sub-walk between the two occurrences. By (P0), its `ℓ⁺`-length is the length of a closed walk in `G`, hence `≥ 0`, so the total length does not increase.
   - **(b) Valley.** Otherwise, if some inner vertex `u` has both walk neighbours `v, w` ranked above it, then `v ≠ w` (by (a)) and `v, w ∈ N↑(u)`, which is a clique, so `{v,w} ∈ E⁺`. Replace `v, u, w` by `v, w`. By (P3), the length does not increase.
   - If neither applies, stop.
4. Each step shortens the walk, so the construction terminates. The final walk `W*` is simple and has no inner vertex below both of its neighbours. Its rank sequence therefore has no inner local minimum, i.e. it first increases and then decreases: `W*` is an up-down path.
5. Its length is `≤ dist(s,t)` by construction, and `≥ dist(s,t)` by step 1. ∎

**Theorem 2 (elimination-tree query).** Let `ℓ` be conservative. Define `d↑(x)` by relaxing the upward arcs of the ancestors of `s` in increasing rank order, starting from `d↑(s) = 0`. Define `d↓(x)` analogously from `t`, using the reversed weights `ℓ⁺(y,x)`. Then `min_x d↑(x) + d↓(x) = dist(s,t)`, where `x` ranges over the common ancestors of `s` and `t`.

*Proof.*

1. `G↑` is acyclic, since every arc goes from lower to higher rank. All vertices reachable from `s` in `G↑` are ancestors of `s`, and the arcs leaving an ancestor end at ancestors: they end in `N↑(x)`, whose minimum is `parent(x)`, and `N↑(x) \ {parent(x)} ⊆ N↑(parent(x))`. So the ancestors of `s`, together with their upward arcs, form a sub-DAG, and increasing rank is a topological order of it.
2. Shortest paths in a DAG are computed exactly by relaxing arcs in topological order, whatever the arc signs. So `d↑(x)` is the minimum length of an upward `s–x` path, and similarly for `d↓`.
3. Every up-down path through its peak `x` is counted, since `x` is then a common ancestor. The claim follows from Theorem 1. ∎

**Query cost.** Let `d` be the depth of the elimination tree, measured in arcs. The ancestors of `s` (including `s`) form a tree path of at most `d + 1` vertices, and every upward arc of an ancestor ends at a higher ancestor. So one search scans at most `d(d+1)/2` arcs, a query (two searches) at most `d(d+1)`, and a query costs `O(d²)` time in the worst case, with no logarithmic factor. For nested-dissection orders on graphs with balanced `O(n^α)` separators, `d = O(n^α)` [Dibbelt et al., Lemma 1]. This worst-case bound is the same one that governs the Dijkstra-based CCH query; it does not depend on highway dimension.

> **Why not the usual CCH query?** The standard CCH query is a bidirectional Dijkstra on `G↑` that stops once the smallest queue key exceeds the best distance found so far. Stall-on-demand additionally prunes vertices that are reachable more cheaply from above.
>
> - Both rules assume that extending a path never makes it shorter. With negative arcs, a vertex with a large tentative distance can still lie on the shortest path, so the stopping rule may terminate too early.
> - Such queries would only become valid again after reweighting with a feasible potential, which is exactly the step we avoid.
> - The elimination-tree query has no stopping rule and no priority queue: it relaxes the whole (small) ancestor DAG in topological order. It is therefore the natural potential-free query, and the only one analysed here.

**Theorem 3 (negative-cycle detection).** After basic customization, `ℓ` has a negative cycle **iff** there is an arc `{v,w} ∈ E⁺` with `ℓ⁺(v,w) + ℓ⁺(w,v) < 0`.

*Proof.*

- (⇐) By (P0), `ℓ⁺(v,w) + ℓ⁺(w,v)` is the length of a closed walk. A negative closed walk contains a negative cycle.
- (⇒) Let `C` be a negative simple cycle with at least 2 vertices; loops are excluded from `E`.
  1. By (P2), its `ℓ⁺`-length is `< 0`.
  2. While `|C| ≥ 3`, remove the lowest vertex `u` of `C`. Its two cycle neighbours are distinct upper neighbours of `u`, hence adjacent in `G⁺`. By (P3) — which does not assume conservativeness — the length does not increase.
  3. We end with a 2-cycle `v → w → v` of negative length. ∎

Checking all arcs takes `O(|E⁺|)` time. This is dominated by customization itself, which touches every arc and every triangle.

**Corollary 3a (local check after an update).** Suppose the metric was conservative before an update, and let `D` be the set of arcs whose customized value changed. Then the updated metric has a negative cycle **iff** some `{v,w} ∈ D` has `ℓ⁺(v,w) + ℓ⁺(w,v) < 0`.

*Proof.* By Theorem 3, a negative cycle exists iff some arc has a negative 2-cycle sum. An arc outside `D` has the same values as before, when the sum was `≥ 0` by Theorem 3 applied to the conservative metric. ∎

**Handling a detected cycle.** A negative cycle means shortest paths are undefined, and **no feasible potential exists**, so falling back to potential shifting is impossible. The only sensible reaction is to reject the update: apply the inverse partial update, which restores the previous customization exactly (Theorem 4), or clamp the offending input weights. Our experiments exercise this reject-and-restore path (§5.5).

**Theorem 4 (partial re-customization).** For an arc `c = {v,w}`, define

    F↑_c = min( ℓ(v,w), min over lower triangles (u; v, w) of ℓ⁺(v,u) + ℓ⁺(u,w) )

and `F↓_c` symmetrically for `ℓ⁺(w,v)`. After basic customization, `ℓ⁺(c) = F_c`, evaluated on the final values of arcs whose lower endpoint ranks below that of `c`.

Suppose the input weights of a set `U` of arcs change. Then the following procedure yields exactly the basic customization of the new metric:

1. Keep a priority queue of arcs, keyed by the rank of their lower endpoint, initialised with `U`.
2. Pop an arc `c` and recompute both of its directions from scratch via `F_c`.
3. If a value of `c` changed, then for every triangle `(u; x, y)` in which `c` is one of the two lower sides `{u,x}` or `{u,y}`, enqueue the upper arc `{x,y}` (if it is not already queued).

The theorem makes no assumption about signs or conservativeness: it only states that partial and full basic customization compute the same values.

*Proof sketch.*

1. By induction on the lower-endpoint rank, an arc that is never enqueued has unchanged base weight and unchanged inputs, so its value is unchanged.
2. An enqueued arc is popped only after every arc with a smaller key. All of its inputs have strictly smaller keys, and new insertions always have keys larger than the key being popped. So the inputs are final when the arc is recomputed, and each arc needs to be evaluated only once. ∎

*Cost.* Let `K` be the number of arcs evaluated. The procedure takes

    O( Σ over evaluated arcs c of ( #lower triangles of c + #triangles with c as a lower side ) + K log K )

time. The `log K` term disappears with a bucket queue over ranks. `K` is at most the number of arcs whose lower endpoint ranks at or above the lowest updated arc. In practice it is far smaller (§5.4).

This differs from the partial update of Dibbelt et al. (§7.7), which propagates changes using tie tests (`m(x,z) + w_old = m(y,z)`). Recomputing from lower triangles handles increases and decreases uniformly, including updates that make weights more negative. We compared the result with full customization after every update in all experiments, and it was identical.

### 3.1 Complexity summary

| Phase | Time | Needs a potential? |
|---|---|---|
| Metric-independent (order, contraction, triangles) | order: heuristic; contraction `O(|E⁺| α(n))` [Dibbelt et al.]; triangles `O(#triangles)` | no |
| Customization + negative-cycle test | `O(#triangles + |E⁺|)` | no |
| Update of `K` arcs + local cycle test | see Theorem 4 / Corollary 3a | no |
| Query | `O(d²)`, `d` = elimination-tree depth | no |
| *Johnson potential (baseline)* | `O(nm)` worst case (Bellman–Ford) | computes one |
| *Near-linear negative SSSP (theory)* | `Õ(m log W)` [BNW22, BCF23] | computes one |

The near-linear algorithms compute a single-source solution or a potential per metric. They are not query structures and would have to be rerun after every metric change. CCH customization instead reuses the metric-independent structure, and its cost does not depend on the magnitude of negative weights.

---

## 4 Implementation

- **Language.** Pure Python 3.14, standard library only.
- **Ordering.** Nested dissection with inertial-flow bisection [Schild & Sommer 2015]: min vertex cuts via Dinic max-flow on the node-split graph, sources and sinks being the extreme 25% along 4 directions, keeping the smallest separator. Geometric median bisection is used as a fallback when a source is adjacent to a sink.
- **Contraction.** Chordal completion via the elimination-tree merge (upper neighbours of `u` are merged into `parent(u)`).
- **Triangles.** Enumerated once, in lower-vertex order, by merging sorted neighbour lists.
- **Customization.** A single pass over the triangle arrays.
- **Queries.** Elimination-tree walks with timestamped labels.
- **Partial updates.** Indices from arcs to their lower triangles, and to the triangles they are a lower side of.
- **Infinity arithmetic (porting pitfall).** We use IEEE floats, where `∞ + x = ∞` for every finite `x`. Integer CCH implementations often encode `∞` as a large constant with saturated addition (Dibbelt et al., §7.6). **With negative weights this breaks:** `INF + (−5)` becomes a finite value just below `INF`, which then wins `min` comparisons against genuinely unreachable entries. A port must test operands for `INF` explicitly before adding, and must choose `INF` and the weight range so that sums of up to `n` negative weights cannot overflow.

**Instances.**

- DIMACS 9th Challenge `USA-road-d` NY (264,346 vertices, 733,846 arcs) and BAY (321,270 / 800,172).
- Two metrics:
  1. **EV terrain.** Real road lengths plus a synthetic smooth elevation field; energy = `len + β·Δz` with β = 1 uphill and 0.6 downhill. The elevation is scaled so that 10% of arcs are negative.
  2. **Random shift.** `len + φ(u) − φ(v)` with random `φ`, giving 20% negative arcs and no exploitable geometry.

  Both are conservative by construction.

**Correctness testing.**

- 996 checks on small grid, geometric and random graphs:
  - queries against Bellman–Ford;
  - the `O(d²)` scan bound;
  - negative-cycle detection against Bellman–Ford, both by full customization and by partial update with the local check;
  - removal of the cycle by the inverse update;
  - partial versus full customization (bit-identical).
- On the road graphs: every CCH answer checked against Johnson + Dijkstra, which is itself validated against Bellman–Ford.

---

## 5 Experiments

All times are single-threaded, pure Python, on the same machine. Absolute times are far above those of compiled implementations. Ratios between methods are indicative only: interpreter overhead differs between operation types (heap operations versus flat array scans), so a C++ port would not necessarily preserve them. Arc-scan and node-expansion counts are implementation-independent. Query times are means over random `s–t` pairs (100 pairs for CCH, 30 for the baselines; all CCH answers checked).

### 5.1 Metric-independent preprocessing

| Graph | Order | Contraction + triangles | Arcs in G⁺ | Triangles | Elim. tree depth |
|---|---:|---:|---:|---:|---:|
| NY | 202 s | 2.5 s | 1,529,403 | 10,057,336 | 344 |
| BAY | 183 s | 1.6 s | 1,242,312 | 5,171,378 | 247 |

Our first, purely geometric ordering on NY produced 4.5 M arcs, 299 M triangles and depth 1,236. The flow-based order reduces the number of triangles 30-fold. Order quality dominates everything downstream.

### 5.2 Metric-dependent preprocessing

| Graph / metric | Johnson (Bellman–Ford) | + 16 ALT landmarks | **CCH customization** |
|---|---:|---:|---:|
| NY / EV terrain | 3.0 s | 9.5 s | **1.6 s** |
| NY / random shift | 0.2 s | 8.9 s | **1.6 s** |
| BAY / EV terrain | 3.6 s | 10.6 s | **0.9 s** |
| BAY / random shift | 0.3 s | 10.2 s | **1.0 s** |

The CCH customization time includes the negative-cycle check.

### 5.3 Queries

| Method | NY terrain | NY shift | BAY terrain | BAY shift |
|---|---:|---:|---:|---:|
| Bellman–Ford (SPFA), no preprocessing* | 10.2 s | 10.5 s | 14.0 s | 14.5 s |
| Johnson + Dijkstra | 116.7 ms | 128.0 ms | 183.3 ms | 187.3 ms |
| Johnson + ALT A* (16 landmarks) | 21.4 ms | 23.1 ms | 44.9 ms | 42.7 ms |
| Bidirectional ALT (1 core) | 14.3 ms | 23.7 ms | 36.7 ms | 24.8 ms |
| **CCH, potential-free** | **0.77 ms** | **0.76 ms** | **0.31 ms** | **0.30 ms** |
| Speed-up vs Johnson + Dijkstra | 151× | 169× | 599× | 622× |

All CCH answers were correct (100/100 per metric).

\* Bellman–Ford times come from an earlier run (20 different random pairs, same graphs and metrics, same machine; `results_dimacs.md`). All other rows use identical pairs.

### 5.4 Updates (traffic scenario)

`k` random arcs get their weight increased by 20–100% of `|w|`, plus 1. After each update, 10 queries were checked, and the customized values were compared against a full customization; they were identical in every case.

| Graph / metric | k = 1 | k = 10 | k = 100 | k = 1000 | Full customization | Redo Johnson + ALT |
|---|---:|---:|---:|---:|---:|---:|
| NY / terrain | 0.1 ms | 25 ms | 125 ms | 645 ms | 1.6 s | 11.5–11.8 s |
| NY / shift | 0.3 ms | 16 ms | 138 ms | 627 ms | 1.6 s | 9.4–9.5 s |
| BAY / terrain | 0.6 ms | 9.5 ms | 48 ms | 271 ms | 0.9 s | 13.2–13.4 s |
| BAY / shift | 0.0 ms | 5.6 ms | 30 ms | 249 ms | 0.9–1.0 s | 10.8–11.0 s |

The update index is built once per topology: 4.7 s for NY and 2.6 s for BAY.

### 5.5 Negative cycles

We injected a 2-cycle of weight −10⁹ via a partial update.

- CCH detected it in 0.03–0.07 s (update + scan).
- Bellman–Ford with periodic predecessor-graph cycle checks took 0.1–0.2 s.
- After the arc was restored, CCH correctly reported no cycle.

### 5.6 Context: potential-based approaches at continental scale

On the full `USA-road-d` graph (23.9 M vertices, 58.3 M arcs), we evaluated the potential-based pipeline with memory-mapped arrays (10 queries each, run one at a time):

| Method | EV terrain | Random shift |
|---|---:|---:|
| Johnson (Bellman–Ford) | 34 s | 18 s |
| ALT landmarks (8, 16 Dijkstras in parallel) | 94 s wall | 61 s wall |
| Johnson + Dijkstra query | 11.9 s | 9.3 s |
| Height-potential A* (no preprocessing) | 10.9 s | — |
| Johnson + ALT A* | 3.7 s | 1.9 s |
| Bidirectional ALT, 1 core / 2 cores | 2.8 s / 1.2 s | 1.5 s / 0.8 s |
| "Textbook" A* with straight-line heuristic | 12.1 s, **5/10 wrong** | — |
| Bellman–Ford per query | > 35 min (aborted) | > 35 min (aborted) |

A potential-free CCH for the full USA graph was not built. The prototype's in-memory representation does not scale to this size in Python (see §6).

---

## 6 Limitations and future work

- **Implementation.** Absolute times come from a Python prototype. Published C++ CCH implementations are 1–2 orders of magnitude faster, and FlowCutter/KaHIP orders are better than our inertial-flow order. The time ratios are indicative only (see §5); scan and expansion counts transfer.
- **Negative cycles.** These are detected exactly, but the only possible "recovery" is to reject or undo the offending update (§3). Policies for real systems — clamping recuperation values, or validating updates before applying them — are application decisions.
- **Other baselines.** We compare against potential-based query methods and Bellman–Ford, not against implementations of the near-linear negative-SSSP algorithms. Those solve a different task (one potential or tree per metric) and would still require a query structure afterwards.
- **Scale.** The full-USA CCH still needs to be built, which requires a compact memory layout.
- **Synthetic energy model.** Real elevation data (e.g. SRTM) and a calibrated EV consumption model are needed.
- **Battery constraints.** State-of-charge bounds make path costs non-additive, as noted in the Scope paragraph of §1. Whether CCH customization can be lifted to the bounded cost functions of Baum et al. without a potential is open.
- **Acceleration techniques.** Only basic customization and the elimination-tree query are studied here. Stall-on-demand is invalid with negative weights (§3). Perfect customization, witness pruning and search-space pruning rules have sign-agnostic-looking proofs, but we have not verified them; some of them may not carry over, which would cost part of the usual CCH speed-up.
- **Parallel customization** carries over unchanged (level-synchronous processing does not depend on sign).
- **Turn costs, one-to-many queries.** Not studied.

---

## 7 Related work

- **Negative-weight SSSP.** Classical Bellman–Ford and Goldberg's scaling algorithm, and recent near-linear algorithms [Bernstein, Nanongkai & Wulff-Nilsen 2022; Bringmann, Cassis & Fischer 2023] with engineered implementations [2025]. Tree-depth-parameterized negative-cycle detection [Iwata, Ogasawara & Ohsaka 2017] is closest in spirit, since it also exploits elimination structure.
- **Path algebras and elimination.** Carré (1971); Lipton, Rose & Tarjan (1979); Tarjan (1981); Rote (1990), who gives a self-contained survey covering negative weights, elimination with short-cut arcs, nested dissection and LU-type factorization. This is the theory our correctness results specialize. Gondran & Minoux (*Graphs, Dioids and Semirings*) collect many applications of these algebraic methods.
- **Functions as arc weights.** Time-dependent route planning (e.g. public transport, and time-dependent contraction hierarchies) attaches a function to each arc, such as the earliest arrival for a given departure time. This is the same "monotone functions under max/min and composition" structure that our battery variant uses. In that setting, the size of the stored functions is a known practical concern. For our battery profiles, whether the number of pieces stays small on large networks is **open**; we observed at most 5 pieces on small graphs only.
- **EV routing.** Eisner, Funke & Storandt (AAAI 2011) use Johnson shifting plus CH. Baum, Dibbelt, Pajor & Wagner (Algorithmica 2020) use height-induced potentials, CH, and battery constraints. Recent A*-based work covers resource-constrained search with negative weights (ESA 2025) and profile search (AAAI 2026).
- **Route planning.** ALT [Goldberg & Harrelson 2005]; CH [Geisberger et al. 2008]; CRP [Delling et al. 2011], which has been extended to EV energy; CCH [Dibbelt, Strasser & Wagner 2016]; CCH survey [Bläsius, Buchhold, Wagner, Zeitz & Zündorf 2025]; inertial flow [Schild & Sommer 2015]; FlowCutter [Hamann & Strasser].
- **Bidirectional and parallel A*.** Ikeda et al. (1994); two-thread bidirectional A* implementations.

---

## Reproducibility

All code is in this repository (see the top-level README). Raw results are in `results/`.

## References (to be completed with full bibliographic data)

- Baum, Dibbelt, Pajor, Wagner. Energy-Optimal Routes for Battery Electric Vehicles. Algorithmica 82, 2020.
- Bläsius, Buchhold, Wagner, Zeitz, Zündorf. Customizable Contraction Hierarchies – A Survey. arXiv:2502.10519, 2025.
- Bernstein, Nanongkai, Wulff-Nilsen. Negative-Weight SSSP in Near-Linear Time. FOCS 2022 / JACM.
- Bringmann, Cassis, Fischer. Negative-Weight SSSP in Near-Linear Time: Now Faster! FOCS 2023.
- Carré. An algebra for network routing problems. 1971.
- Dibbelt, Strasser, Wagner. Customizable Contraction Hierarchies. ACM JEA 21, 2016 (arXiv:1402.0402).
- Eisner, Funke, Storandt. Optimal Route Planning for Electric Vehicles in Large Networks. AAAI 2011.
- Goldberg, Harrelson. Computing the shortest path: A* search meets graph theory. SODA 2005.
- Iwata, Ogasawara, Ohsaka. On the Power of Tree-Depth for Fully Polynomial FPT Algorithms. STACS 2018 (arXiv:1710.04376).
- Johnson. Efficient algorithms for shortest paths in sparse networks. JACM 1977.
- Schild, Sommer. On balanced separators in road networks. SEA 2015.
- Tarjan. Fast algorithms for solving path problems. JACM 28(3), 1981.
- Rote. Path problems in graphs. In: Computational Graph Theory, Computing Supplementum 7, Springer, 1990, pp. 155–189. https://page.mi.fu-berlin.de/rote/Papers/pdf/Path+problems+in+graphs.pdf
- Lipton, Rose, Tarjan. Generalized nested dissection. SIAM J. Numer. Anal. 16(2), 1979.
- Gondran, Minoux. Graphs, Dioids and Semirings: New Models and Algorithms. Springer, 2008.
- 9th DIMACS Implementation Challenge – Shortest Paths.
