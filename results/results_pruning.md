# Which CCH acceleration techniques survive negative weights?

`test_cch_pruning.py`, on 9 small graphs with negative edges (grid, geometric,
random-shift), checked against Bellman-Ford.

| Technique | Verdict | Evidence |
|---|---|---|
| **Basic customization** | works | proven (technical note) and tested everywhere |
| **Elimination-tree (sweep) query** | works | ditto |
| **Perfect customization** | **works** | after the top-down pass every arc held the true distance: 0 mismatches in 3,618 arc distances |
| **Witness pruning / arc removal** | **works, with one restriction** | only arcs whose detour runs through a *higher*-ranked vertex (upper or intermediate triangle) may be removed. Then 0 of 225 queries were wrong, and about 60% of arc directions (32,383 of 54,310) were still removable. Dropping arcs whose detour uses a lower-ranked vertex breaks correctness (183 of 225 queries wrong in an earlier run). |
| **Stall-on-demand** | **unsafe** | in 850 cases a vertex would have been stalled although its sweep label was already the exact distance; the stall test assumes a Dijkstra-style order that negative weights destroy |

This matches the caveat raised by S. Storandt (personal communication): most of
the CCH speed-up survives, stalling is the casualty.
