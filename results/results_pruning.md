# Which CCH acceleration techniques survive negative weights?

`python/test_cch_pruning.py`, on 9 small graphs with negative edges (grid, geometric,
random-shift), checked against Bellman-Ford.

| Technique | Verdict | Evidence |
|---|---|---|
| **Basic customization** | works | proven (technical note) and tested everywhere |
| **Elimination-tree (sweep) query** | works | ditto |
| **Perfect customization** | **works** | after the top-down pass every arc held the true distance: 0 mismatches in 3,618 arc distances |
| **Witness pruning / arc removal** | **works, with one restriction** | only arcs whose detour runs through a *higher*-ranked vertex (upper or intermediate triangle) may be removed. Then 0 of 225 queries were wrong, and about 60% of arc directions (32,383 of 54,310) were still removable. Dropping arcs whose detour uses a lower-ranked vertex breaks correctness (183 of 225 queries wrong in an earlier run). |
| **Stall-on-demand** | **unsafe** | in 850 cases a vertex would have been stalled although its sweep label was already the exact distance; the stall test assumes a Dijkstra-style order that negative weights destroy |

## Further techniques (`python/test_cch_extras.py`, and `cch --parallel` in C++)

| Technique | Verdict | Evidence |
|---|---|---|
| **Path unpacking** | **works** | 140 of 140 unpacked paths are real edge sequences of exactly the shortest length. Each shortcut is replaced by the lower triangle that realises its value; every shortcut we examined had such a triangle (0 unexplained of 944). |
| **Elimination-tree query on a pruned search graph** | **works** | the pruning test above uses exactly this combination (witness-pruned arcs + sweep query) |
| **Parallel customization by levels** | **correct, but not faster here** | with one thread per vertex of a level and atomic minima, the result is bit-identical to the serial one (0 differing arcs on New York and on the full USA). Speed: 1.6-1.9x on New York (16 threads), 0.9x on the USA (32 threads) - the USA hierarchy has 3,763 levels, so barriers and atomics cost more than they save. Published implementations parallelise *per arc* instead, which needs a stored triangle index (3.15 billion triangles on the USA, which we deliberately do not store). |
| **Tie-based partial updates** (CCH paper, §7.7) | **works, with one addition** | propagating only where the old value was the witness (the equality test alone) missed updates in 26 of 42 runs. Adding the second condition - propagate also where the new value *improves* the target - makes it exact: 42 of 42 runs matched full customization. |

This matches the caveat raised by S. Storandt (personal communication): most of
the CCH speed-up survives, stalling is the casualty.
