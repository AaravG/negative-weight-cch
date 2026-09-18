# CCH with negative weights and no potential

All CCH answers are checked against Johnson + Dijkstra. Pure Python; times are single-threaded.

## NY (n=264,346, m=733,846)

CCH build (metric-independent, done once for both weightings): **205 s** (order 202 s, contraction 0 s, triangles 2 s). 1,529,403 arcs, 10,057,336 triangles, elimination-tree depth 344.

### EV terrain (10% negative edges)

| Preprocessing for this metric | Time |
|---|---:|
| Johnson's Bellman-Ford | 3.0 s |
| + 16 ALT landmarks | 9.5 s |
| **CCH customization (includes negative-cycle check)** | **1.6 s** |

| Query method | Mean query (ms) | Correct |
|---|---:|---:|
| Johnson + Dijkstra | 116.65 (1.0x) | 30/30 |
| Johnson + ALT (16) | 21.44 (5.4x) | 30/30 |
| Bidirectional ALT, 1 core | 14.33 (8.1x) | 30/30 |
| **CCH, no potential** | 0.77 (151.4x) | 100/100 |

Update index (one-time, metric-independent): 4.7 s

| Edges changed | Partial update | Arcs re-evaluated | Full customization | Redo Johnson + ALT | Correct after update |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.1 ms | 57 | 1.6 s | 11.5 s | 10/10, arcs identical: True |
| 10 | 25.5 ms | 4,524 | 1.6 s | 11.7 s | 10/10, arcs identical: True |
| 100 | 124.7 ms | 23,340 | 1.6 s | 11.8 s | 10/10, arcs identical: True |
| 1000 | 645.4 ms | 166,538 | 1.6 s | 11.7 s | 10/10, arcs identical: True |

Negative cycle injected: CCH detects it: **True** (partial update + scan 0.07 s); Bellman-Ford detects it: True (0.1 s). After restoring: CCH reports a cycle: False.

### random shift (20% negative edges)

| Preprocessing for this metric | Time |
|---|---:|
| Johnson's Bellman-Ford | 0.2 s |
| + 16 ALT landmarks | 8.9 s |
| **CCH customization (includes negative-cycle check)** | **1.6 s** |

| Query method | Mean query (ms) | Correct |
|---|---:|---:|
| Johnson + Dijkstra | 128.02 (1.0x) | 30/30 |
| Johnson + ALT (16) | 23.12 (5.5x) | 30/30 |
| Bidirectional ALT, 1 core | 23.67 (5.4x) | 30/30 |
| **CCH, no potential** | 0.76 (168.6x) | 100/100 |

| Edges changed | Partial update | Arcs re-evaluated | Full customization | Redo Johnson + ALT | Correct after update |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.3 ms | 116 | 1.6 s | 9.5 s | 10/10, arcs identical: True |
| 10 | 16.4 ms | 2,774 | 1.6 s | 9.5 s | 10/10, arcs identical: True |
| 100 | 138.3 ms | 27,873 | 1.6 s | 9.4 s | 10/10, arcs identical: True |
| 1000 | 627.2 ms | 161,138 | 1.6 s | 9.5 s | 10/10, arcs identical: True |

Negative cycle injected: CCH detects it: **True** (partial update + scan 0.07 s); Bellman-Ford detects it: True (0.2 s). After restoring: CCH reports a cycle: False.

## BAY (n=321,270, m=800,172)

CCH build (metric-independent, done once for both weightings): **185 s** (order 183 s, contraction 1 s, triangles 1 s). 1,242,312 arcs, 5,171,378 triangles, elimination-tree depth 247.

### EV terrain (10% negative edges)

| Preprocessing for this metric | Time |
|---|---:|
| Johnson's Bellman-Ford | 3.6 s |
| + 16 ALT landmarks | 10.6 s |
| **CCH customization (includes negative-cycle check)** | **0.9 s** |

| Query method | Mean query (ms) | Correct |
|---|---:|---:|
| Johnson + Dijkstra | 183.25 (1.0x) | 30/30 |
| Johnson + ALT (16) | 44.93 (4.1x) | 30/30 |
| Bidirectional ALT, 1 core | 36.68 (5.0x) | 30/30 |
| **CCH, no potential** | 0.31 (598.9x) | 100/100 |

Update index (one-time, metric-independent): 2.6 s

| Edges changed | Partial update | Arcs re-evaluated | Full customization | Redo Johnson + ALT | Correct after update |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.6 ms | 193 | 0.9 s | 13.2 s | 10/10, arcs identical: True |
| 10 | 9.5 ms | 2,513 | 0.9 s | 13.3 s | 10/10, arcs identical: True |
| 100 | 47.6 ms | 13,187 | 0.9 s | 13.4 s | 10/10, arcs identical: True |
| 1000 | 270.6 ms | 90,304 | 0.9 s | 13.4 s | 10/10, arcs identical: True |

Negative cycle injected: CCH detects it: **True** (partial update + scan 0.03 s); Bellman-Ford detects it: True (0.1 s). After restoring: CCH reports a cycle: False.

### random shift (20% negative edges)

| Preprocessing for this metric | Time |
|---|---:|
| Johnson's Bellman-Ford | 0.3 s |
| + 16 ALT landmarks | 10.2 s |
| **CCH customization (includes negative-cycle check)** | **1.0 s** |

| Query method | Mean query (ms) | Correct |
|---|---:|---:|
| Johnson + Dijkstra | 187.30 (1.0x) | 30/30 |
| Johnson + ALT (16) | 42.71 (4.4x) | 30/30 |
| Bidirectional ALT, 1 core | 24.77 (7.6x) | 30/30 |
| **CCH, no potential** | 0.30 (621.5x) | 100/100 |

| Edges changed | Partial update | Arcs re-evaluated | Full customization | Redo Johnson + ALT | Correct after update |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.0 ms | 1 | 0.9 s | 10.9 s | 10/10, arcs identical: True |
| 10 | 5.6 ms | 1,561 | 0.9 s | 10.8 s | 10/10, arcs identical: True |
| 100 | 29.9 ms | 9,386 | 0.9 s | 10.8 s | 10/10, arcs identical: True |
| 1000 | 249.4 ms | 85,384 | 1.0 s | 11.0 s | 10/10, arcs identical: True |

Negative cycle injected: CCH detects it: **True** (partial update + scan 0.03 s); Bellman-Ford detects it: True (0.2 s). After restoring: CCH reports a cycle: False.
