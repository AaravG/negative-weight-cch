# CCH with negative weights and no potential

All CCH answers are checked against Johnson + Dijkstra. Pure Python; times are single-threaded.

## CAL_R (n=1,890,815, m=4,657,742)

CCH build (metric-independent, done once for both weightings): **3275 s** (order 3253 s, contraction 5 s, triangles 13 s). 6,875,035 arcs, 34,189,458 triangles, elimination-tree depth 571.

### EV terrain (10% negative edges)

| Preprocessing for this metric | Time |
|---|---:|
| Johnson's Bellman-Ford | 124.1 s |
| + 16 ALT landmarks | 96.4 s |
| **CCH customization (includes negative-cycle check)** | **11.0 s** |

| Query method | Mean query (ms) | Correct |
|---|---:|---:|
| Johnson + Dijkstra | 1366.42 (1.0x) | 30/30 |
| Johnson + ALT (16) | 472.37 (2.9x) | 30/30 |
| Bidirectional ALT, 1 core | 592.75 (2.3x) | 30/30 |
| **CCH, no potential** | 3.22 (425.0x) | 100/100 |

Update index (one-time, metric-independent): 31.5 s

| Edges changed | Partial update | Arcs re-evaluated | Full customization | Redo Johnson + ALT | Correct after update |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.0 ms | 1 | 11.0 s | 218.0 s | 10/10, arcs identical: True |
| 10 | 50.0 ms | 2,627 | 10.9 s | 218.7 s | 10/10, arcs identical: True |
| 100 | 343.9 ms | 24,704 | 10.6 s | 214.8 s | 10/10, arcs identical: True |
| 1000 | 2255.9 ms | 220,376 | 10.6 s | 217.0 s | 10/10, arcs identical: True |

Negative cycle injected: CCH detects it: **True** (partial update + scan 0.22 s); Bellman-Ford detects it: True (1.5 s). After restoring: CCH reports a cycle: False.

### random shift (20% negative edges)

| Preprocessing for this metric | Time |
|---|---:|
| Johnson's Bellman-Ford | 2.5 s |
| + 16 ALT landmarks | 92.2 s |
| **CCH customization (includes negative-cycle check)** | **10.7 s** |

| Query method | Mean query (ms) | Correct |
|---|---:|---:|
| Johnson + Dijkstra | 1321.85 (1.0x) | 30/30 |
| Johnson + ALT (16) | 324.11 (4.1x) | 30/30 |
| Bidirectional ALT, 1 core | 504.34 (2.6x) | 30/30 |
| **CCH, no potential** | 2.97 (444.7x) | 100/100 |

| Edges changed | Partial update | Arcs re-evaluated | Full customization | Redo Johnson + ALT | Correct after update |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.2 ms | 19 | 10.7 s | 94.6 s | 10/10, arcs identical: True |
| 10 | 59.6 ms | 3,584 | 10.7 s | 95.5 s | 10/10, arcs identical: True |
| 100 | 477.7 ms | 35,056 | 10.8 s | 95.0 s | 10/10, arcs identical: True |
| 1000 | 2272.5 ms | 214,553 | 11.0 s | 95.6 s | 10/10, arcs identical: True |

Negative cycle injected: CCH detects it: **True** (partial update + scan 0.22 s); Bellman-Ford detects it: True (1.7 s). After restoring: CCH reports a cycle: False.
