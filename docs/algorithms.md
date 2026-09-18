# Algorithms and concepts used in this project

A plain-language glossary. Items marked (impl.) are implemented in this repository.

## Basics

- **Negative edge.** An edge with cost < 0. In EV routing, a downhill segment where regenerative braking recovers energy.
- **Negative cycle.** A loop with total cost < 0. Going around it forever makes paths arbitrarily cheap, so shortest paths are undefined.
- **Conservative weights.** Weights without negative cycles. Every method here assumes this. EV energy is conservative: you cannot gain energy by driving in a circle.
- **Dijkstra.** Settles vertices in increasing distance order. Correct only for non-negative weights.
- **Bellman–Ford** (impl.). Relaxes all edges repeatedly (up to n−1 rounds). Handles negative edges, but costs O(nm).
- **SPFA / queue-based Bellman–Ford** (impl.). Only rechecks vertices whose distance just changed. Used as the no-preprocessing baseline.
- **Label-correcting search** (impl.). A vertex's label may be improved many times, as in SPFA. The battery reference search is of this kind.
- **Negative-cycle detection in Bellman–Ford** (impl.). Periodically check whether the predecessor pointers form a cycle; such a cycle is always negative.

## Potentials and A*

- **Potential.** A number `p(v)` per vertex. Reweighting `w'(u,v) = w(u,v) + p(u) − p(v)` changes every `s–t` path by the same constant, so shortest paths are preserved.
- **Feasible potential.** One that makes every `w' ≥ 0`, so Dijkstra becomes applicable. It exists iff the weights are conservative.
- **Johnson's algorithm** (impl.). Bellman–Ford from a virtual source gives a feasible potential; then Dijkstra.
- **A\*** (impl.). Dijkstra ordered by `g(v) + h(v)`, where `h` estimates the remaining distance. Equivalent to Dijkstra with `h` used as a potential.
- **Admissible heuristic.** Never overestimates the remaining distance.
- **Consistent heuristic.** `h(u) ≤ w(u,v) + h(v)` for every edge. A* then never re-expands a vertex.
- **Textbook A\* with a straight-line heuristic** (impl., as a warning). With negative edges, the straight-line energy is not a lower bound, so A* returns wrong answers (about 50% in our tests).

## Goal-directed heuristics

- **ALT: A\*, Landmarks, Triangle inequality** (impl.). Precompute distances to and from a few landmark vertices. Then `dist(v,t) ≥ dist(v,L) − dist(t,L)` and `dist(v,t) ≥ dist(L,t) − dist(L,v)`; take the maximum over all landmarks.
- **Landmark selection** (impl.). "Farthest-point": greedily pick vertices far from those already chosen. "Planar": the farthest vertex in each angular sector around the centre.
- **Domain bound** (impl.). EV-specific lower bound from straight-line distance plus height difference. Needs no preprocessing.
- **Combining bounds** (impl.). The pointwise maximum of consistent lower bounds is still a consistent lower bound.

## Search strategies

- **Bidirectional search** (impl.). A forward search from `s` and a backward search from `t` run together and stop when they provably meet.
- **Average potential** (Ikeda et al. 1994) (impl.). Forward uses `(π_t − π_s)/2`, backward uses its negation. Both searches then see the same non-negative reduced costs.
- **Two-process bidirectional A\*** (impl.). Each direction runs on its own core, sharing distance arrays and the best meeting cost through shared memory.

## Contraction hierarchies

- **Contraction Hierarchies (CH).** Rank vertices by importance and contract them from least to most important, adding shortcuts that preserve distances. Queries only move upward in rank.
- **Shortcut.** An added edge `u–w` standing for a path `u → v → w` through a contracted vertex `v`.
- **Up-down path.** A path whose ranks first increase, then decrease. After (C)CH preprocessing, a shortest path of this shape always exists.
- **CCH** (impl.). CH split into a metric-independent phase (order + shortcuts), customization (shortcut weights), and queries.
- **Chordal supergraph / lower triangle.** The graph after adding all shortcuts. A lower triangle `(u; v, w)` has `u` ranked below `v` and `w`, and gives the rule `ℓ(v,w) ≤ ℓ(v,u) + ℓ(u,w)`.
- **Basic customization** (impl.). Apply that rule to all lower triangles, lowest `u` first.
- **Elimination tree** (impl.). `parent(v)` is the lowest-ranked upper neighbour of `v`. Everything reachable upward from `v` lies on its path to the root.
- **Elimination-tree query** (impl.). Relax the upward arcs of the ancestors of `s` (and, reversed, of `t`) in rank order, then take the best common ancestor. Costs `O(depth²)`.
- **Partial re-customization** (impl.). After some input weights change, recompute only the affected shortcuts from their lower triangles, lowest first.

## Ordering

- **Nested dissection** (impl.). Find a small separator, give it the highest ranks, and recurse on both halves.
- **Vertex separator.** A vertex set whose removal splits the graph into two parts.
- **Inertial flow** (Schild & Sommer 2015) (impl.). Sort vertices along a direction, take the first and last 25% as source and sink, and compute a minimum vertex cut between them. We try 4 directions and keep the smallest cut.
- **Max-flow / min-cut.** The maximum flow between two sets equals the capacity of the smallest cut separating them.
- **Dinic's algorithm** (impl.). A max-flow algorithm: BFS builds level graphs, then augmenting paths are pushed in phases.
- **Node splitting** (impl.). Each vertex `v` becomes `v_in → v_out` with capacity 1, so a minimum cut consists of vertices.

## Algebra

- **(min, +) semiring.** Shortest paths combine paths with `+` and choose between them with `min`.
- **Vertex elimination / Gaussian elimination for path problems** (Carré 1971; Tarjan 1981). Eliminating vertices one by one while updating their neighbours. Exact whenever there is no negative cycle — which is why CCH customization works with negative weights.

## Battery-constrained routing

- **State of charge (SoC).** The battery energy `b ∈ [0, M]`.
- **Charge function** (impl.). An arc maps `b ↦ min(out, b − cost)` if `b ≥ in`, and is infeasible otherwise (Eisner, Funke & Storandt 2011).
- **Composition** (impl.). Driving arc A then arc B gives `B(A(b))`, which has the same 3-parameter form.
- **Profile / Pareto set** (impl.). A shortcut stores the few non-dominated charge functions of its paths; its value is their pointwise maximum.
- **(max, ∘) structure.** The battery analogue of (min, +): `max` chooses the better route, composition chains roads. The "no gain cycles" condition replaces "no negative cycles".

## Mentioned but not implemented

- **Near-linear negative-weight SSSP** (Bernstein–Nanongkai–Wulff-Nilsen 2022; Bringmann–Cassis–Fischer 2023). Theoretical breakthroughs for single-source shortest paths with negative weights.
- **Goldberg–Radzik.** A practical improvement of Bellman–Ford.
- **CRP (Customizable Route Planning).** A multi-level overlay technique with a similar three-phase design.
- **RoutingKit.** An open-source C++ library with a CCH implementation.

## Engineering

- **DIMACS road graphs.** Standard benchmark road networks (NY, BAY, full USA) from the 9th DIMACS Implementation Challenge.
- **CSR (compressed sparse row).** A graph stored as flat offset and neighbour arrays.
- **Memory-mapped files.** Worker processes share one copy of the graph through the OS page cache.
