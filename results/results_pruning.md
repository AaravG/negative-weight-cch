# Which CCH techniques survive negative weights?

Reproduce with `cd python && python test_techniques.py --exhaustive`. Every check is run against
Bellman-Ford and repeated on a non-negative control (same topology and order, weights |w|), so
that effects of negative weights can be told apart from standard CCH conditions.

| Technique | Negative weights (no negative cycles) | Evidence |
|---|---|---|
| Basic customization, elimination-tree query (no distance pruning) | exact | proofs; all tests |
| Stall-on-demand in the elimination-tree query | **exact** | proof; 0 wrong in randomized tests (5,003 stalls) and on all 14,440 conservative 3-vertex digraphs |
| Perfect customization | exact | proof; 0 mismatches on all arcs of 18 graphs |
| Witness pruning (upper/intermediate triangles) | exact if no zero-length cycles | proof; with zero cycles it fails for non-negative weights too (3-vertex example), fixed by a lexicographic tie-break (0 wrong) |
| Removal via lower triangles | invalid, for non-negative weights too | 3-vertex example with positive weights |
| Path unpacking (basic values) | exact | proof |
| Partial updates, negative-cycle test | exact | proofs; `test_cch.py` |
| Level-parallel customization | bit-identical to serial | C++ exact comparison on NY, BAY, USA |
| Potential read off the CCH (two sweeps) | exact, feasible | equals Bellman-Ford on all test graphs |
| Dijkstra-based query (label-setting or stopping rule) | **fails** | 3-vertex counterexample; 1,398 wrong on the exhaustive set; controls 0 |
| Distance pruning of the ET query (Buchhold et al.) | **fails** | 3-vertex counterexample; 560 wrong; controls 0 |
| Early termination of the ET query | **fails** | 3-vertex counterexample; 68 wrong; controls 0 |

Earlier versions of this file called stall-on-demand unsafe. That was wrong: the old test counted
stall events without checking answers, and every stalled label was in fact not a shortest distance.
