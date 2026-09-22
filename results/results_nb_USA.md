# Potential-free CCH (Numba) on USA

n = 23,947,347, original arcs = 58,333,344. Ordering: 1421 s. Symbolic contraction: 6 s. Shortcut arcs: 96,738,990. Triangles: 3,152,158,103. Elimination-tree depth: 3771.

## ev

| Measure | Value |
|---|---:|
| Customization (compiled) | 7.2 s |
| Negative-cycle scan | 0.03 s (0 cycle arcs) |
| Median CCH query | 1.01 ms |
| Median Johnson + Dijkstra query (compiled, reference) | 1522 ms |
| Correct | 100/100 |

## shifted

| Measure | Value |
|---|---:|
| Customization (compiled) | 7.3 s |
| Negative-cycle scan | 0.03 s (0 cycle arcs) |
| Median CCH query | 0.98 ms |
| Median Johnson + Dijkstra query (compiled, reference) | 1467 ms |
| Correct | 100/100 |
