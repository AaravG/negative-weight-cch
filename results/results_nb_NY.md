# Potential-free CCH (Numba) on NY

n = 264,346, original arcs = 733,846. Ordering: 12 s. Symbolic contraction: 0 s. Shortcut arcs: 1,531,232. Triangles: 10,164,101. Elimination-tree depth: 343.

## ev

| Measure | Value |
|---|---:|
| Customization (compiled) | 0.1 s |
| Negative-cycle scan | 0.01 s (0 cycle arcs) |
| Median CCH query | 0.04 ms |
| Median Johnson + Dijkstra query (compiled, reference) | 19 ms |
| Correct | 100/100 |

## shifted

| Measure | Value |
|---|---:|
| Customization (compiled) | 0.1 s |
| Negative-cycle scan | 0.00 s (0 cycle arcs) |
| Median CCH query | 0.04 ms |
| Median Johnson + Dijkstra query (compiled, reference) | 17 ms |
| Correct | 100/100 |
