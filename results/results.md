# Benchmark results

## EV grid 100x100

| Method | Preprocess (s) | Query (ms) | Nodes expanded | Speedup vs BF | Correct | Break-even queries |
|---|---:|---:|---:|---:|---:|---:|
| Bellman-Ford (SPFA), no preprocessing | 0.0 | 5.48 | 10,004 (100.0%) | 1.0x | 30/30 | - |
| Normal A* (straight-line h, reopening) | 0.0 | 8.16 | 8,590 (85.9%) | 0.7x | 11/30 | - |
| Johnson + Dijkstra | 0.016 | 4.87 | 5,066 (50.7%) | 1.1x | 30/30 | 25.7 |
| Johnson + ALT A* (16 landmarks) | 0.515 | 0.79 | 202 (2.0%) | 7.0x | 30/30 | 109.8 |
| Domain-bound A* (no preprocessing) | 0.0 | 1.66 | 1,265 (12.7%) | 3.3x | 30/30 | - |
| Combined max(ALT, domain) A* | 0.515 | 1.01 | 195 (2.0%) | 5.4x | 30/30 | 115.4 |
| Bidirectional max(alt,domain), 1 core | 0.515 | 1.98 | 245 (2.5%) | 2.8x | 30/30 | 147.3 |
| Bidirectional max(alt,domain), 2 cores | 0.515 | 1.1 | 223 (2.2%) | 5.0x | 30/30 | 117.7 |

## EV grid 200x200

| Method | Preprocess (s) | Query (ms) | Nodes expanded | Speedup vs BF | Correct | Break-even queries |
|---|---:|---:|---:|---:|---:|---:|
| Bellman-Ford (SPFA), no preprocessing | 0.0 | 26.05 | 40,000 (100.0%) | 1.0x | 30/30 | - |
| Normal A* (straight-line h, reopening) | 0.0 | 16.1 | 15,256 (38.1%) | 1.6x | 13/30 | - |
| Johnson + Dijkstra | 0.017 | 21.91 | 18,893 (47.2%) | 1.2x | 30/30 | 4.2 |
| Johnson + ALT A* (16 landmarks) | 2.02 | 3.45 | 1,028 (2.6%) | 7.5x | 30/30 | 89.4 |
| Domain-bound A* (no preprocessing) | 0.0 | 4.22 | 2,877 (7.2%) | 6.2x | 30/30 | - |
| Combined max(ALT, domain) A* | 2.02 | 4.15 | 1,000 (2.5%) | 6.3x | 30/30 | 92.2 |
| Bidirectional max(alt,domain), 1 core | 2.02 | 7.56 | 1,104 (2.8%) | 3.4x | 30/30 | 109.3 |
| Bidirectional max(alt,domain), 2 cores | 2.02 | 3.62 | 1,087 (2.7%) | 7.2x | 30/30 | 90.1 |

## EV road-like 30k

| Method | Preprocess (s) | Query (ms) | Nodes expanded | Speedup vs BF | Correct | Break-even queries |
|---|---:|---:|---:|---:|---:|---:|
| Bellman-Ford (SPFA), no preprocessing | 0.0 | 44.61 | 100,946 (336.5%) | 1.0x | 30/30 | - |
| Normal A* (straight-line h, reopening) | 0.0 | 73.39 | 84,205 (280.7%) | 0.6x | 9/30 | - |
| Johnson + Dijkstra | 0.116 | 15.33 | 13,513 (45.0%) | 2.9x | 30/30 | 3.9 |
| Johnson + ALT A* (16 landmarks) | 1.429 | 3.24 | 709 (2.4%) | 13.8x | 30/30 | 34.5 |
| Domain-bound A* (no preprocessing) | 0.0 | 6.0 | 4,534 (15.1%) | 7.4x | 30/30 | - |
| Combined max(ALT, domain) A* | 1.429 | 3.88 | 699 (2.3%) | 11.5x | 30/30 | 35.1 |
| Bidirectional max(alt,domain), 1 core | 1.429 | 7.06 | 713 (2.4%) | 6.3x | 30/30 | 38.0 |
| Bidirectional max(alt,domain), 2 cores | 1.429 | 3.74 | 744 (2.5%) | 11.9x | 30/30 | 35.0 |

## Shifted random 20k (no geometry)

| Method | Preprocess (s) | Query (ms) | Nodes expanded | Speedup vs BF | Correct | Break-even queries |
|---|---:|---:|---:|---:|---:|---:|
| Bellman-Ford (SPFA), no preprocessing | 0.0 | 17.62 | 40,985 (204.9%) | 1.0x | 30/30 | - |
| Johnson + Dijkstra | 0.011 | 13.49 | 10,162 (50.8%) | 1.3x | 30/30 | 2.7 |
| Johnson + ALT A* (16 landmarks) | 0.93 | 7.57 | 633 (3.2%) | 2.3x | 30/30 | 92.5 |
| Bidirectional alt, 1 core | 0.93 | 5.31 | 333 (1.7%) | 3.3x | 30/30 | 75.5 |
| Bidirectional alt, 2 cores | 0.93 | 2.31 | 289 (1.4%) | 7.6x | 30/30 | 60.7 |
