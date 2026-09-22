# CCH with negative weights and no potential

All CCH answers are checked against Johnson + Dijkstra. Pure Python; times are single-threaded.

## FLA_R (n=1,095,184, m=2,772,490)

CCH build (metric-independent, done once for both weightings): **1450 s** (order 1440 s, contraction 3 s, triangles 6 s). 3,958,007 arcs, 14,486,517 triangles, elimination-tree depth 344.

### EV terrain (10% negative edges)

| Preprocessing for this metric | Time |
|---|---:|
| Johnson's Bellman-Ford | 53.6 s |
| + 16 ALT landmarks | 54.6 s |
| **CCH customization (includes negative-cycle check)** | **5.2 s** |

| Query method | Mean query (ms) | Correct |
|---|---:|---:|
| Johnson + Dijkstra | 931.41 (1.0x) | 30/30 |
| Johnson + ALT (16) | 252.42 (3.7x) | 30/30 |
| Bidirectional ALT, 1 core | 280.84 (3.3x) | 30/30 |
| **CCH, no potential** | 1.15 (807.3x) | 100/100 |

Update index (one-time, metric-independent): 13.6 s

| Edges changed | Partial update | Arcs re-evaluated | Full customization | Redo Johnson + ALT | Correct after update |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.1 ms | 14 | 5.0 s | 105.5 s | 10/10, arcs identical: True |
| 10 | 20.4 ms | 1,711 | 5.3 s | 105.8 s | 10/10, arcs identical: True |
| 100 | 74.9 ms | 11,287 | 5.2 s | 106.3 s | 10/10, arcs identical: True |
| 1000 | 648.3 ms | 108,152 | 5.2 s | 105.4 s | 10/10, arcs identical: True |

Negative cycle injected: CCH detects it: **True** (partial update + scan 0.09 s); Bellman-Ford detects it: True (0.9 s). After restoring: CCH reports a cycle: False.

### random shift (20% negative edges)

| Preprocessing for this metric | Time |
|---|---:|
| Johnson's Bellman-Ford | 1.5 s |
| + 16 ALT landmarks | 53.3 s |
| **CCH customization (includes negative-cycle check)** | **5.3 s** |

| Query method | Mean query (ms) | Correct |
|---|---:|---:|
| Johnson + Dijkstra | 994.37 (1.0x) | 30/30 |
| Johnson + ALT (16) | 228.72 (4.3x) | 30/30 |
| Bidirectional ALT, 1 core | 222.13 (4.5x) | 30/30 |
| **CCH, no potential** | 1.15 (867.2x) | 100/100 |

| Edges changed | Partial update | Arcs re-evaluated | Full customization | Redo Johnson + ALT | Correct after update |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.8 ms | 191 | 5.2 s | 54.8 s | 10/10, arcs identical: True |
| 10 | 2.0 ms | 316 | 5.1 s | 55.5 s | 10/10, arcs identical: True |
| 100 | 111.2 ms | 14,394 | 5.2 s | 55.5 s | 10/10, arcs identical: True |
| 1000 | 657.2 ms | 116,102 | 5.1 s | 55.6 s | 10/10, arcs identical: True |

Negative cycle injected: CCH detects it: **True** (partial update + scan 0.10 s); Bellman-Ford detects it: True (1.0 s). After restoring: CCH reports a cycle: False.
