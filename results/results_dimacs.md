# Benchmark results

## DIMACS NY, EV terrain

| Method | Preprocess (s) | Query (ms) | Nodes expanded | Speedup vs BF | Correct | Break-even queries |
|---|---:|---:|---:|---:|---:|---:|
| Bellman-Ford (SPFA), no preprocessing | 0.0 | 10235.61 | 16,584,555 (6273.8%) | 1.0x | 20/20 | - |
| Normal A* (straight-line h, reopening) | 0.0 | 2283.34 | 2,892,200 (1094.1%) | 4.5x | 10/20 | - |
| Johnson + Dijkstra | 2.444 | 152.07 | 135,169 (51.1%) | 67.3x | 20/20 | 0.2 |
| Johnson + ALT A* (16 landmarks) | 11.89 | 16.97 | 5,773 (2.2%) | 603.3x | 20/20 | 1.2 |
| Domain-bound A* (no preprocessing) | 0.0 | 65.37 | 48,864 (18.5%) | 156.6x | 20/20 | - |
| Combined max(ALT, domain) A* | 11.89 | 20.84 | 5,700 (2.2%) | 491.2x | 20/20 | 1.2 |
| Bidirectional max(alt,domain), 1 core | 11.89 | 22.11 | 3,865 (1.5%) | 462.9x | 20/20 | 1.2 |
| Bidirectional max(alt,domain), 2 cores | 11.89 | 10.02 | 3,455 (1.3%) | 1021.3x | 20/20 | 1.2 |

## DIMACS NY, random potential shift

| Method | Preprocess (s) | Query (ms) | Nodes expanded | Speedup vs BF | Correct | Break-even queries |
|---|---:|---:|---:|---:|---:|---:|
| Bellman-Ford (SPFA), no preprocessing | 0.0 | 10543.73 | 17,074,664 (6459.2%) | 1.0x | 20/20 | - |
| Johnson + Dijkstra | 0.177 | 140.27 | 128,207 (48.5%) | 75.2x | 20/20 | 0.0 |
| Johnson + ALT A* (16 landmarks) | 9.389 | 18.48 | 5,926 (2.2%) | 570.5x | 20/20 | 0.9 |
| Bidirectional alt, 1 core | 9.389 | 16.36 | 3,626 (1.4%) | 644.4x | 20/20 | 0.9 |
| Bidirectional alt, 2 cores | 9.389 | 8.45 | 3,692 (1.4%) | 1248.0x | 20/20 | 0.9 |

## DIMACS BAY, EV terrain

| Method | Preprocess (s) | Query (ms) | Nodes expanded | Speedup vs BF | Correct | Break-even queries |
|---|---:|---:|---:|---:|---:|---:|
| Bellman-Ford (SPFA), no preprocessing | 0.0 | 13983.41 | 22,359,674 (6959.8%) | 1.0x | 20/20 | - |
| Normal A* (straight-line h, reopening) | 0.0 | 286.16 | 417,823 (130.1%) | 48.9x | 12/20 | - |
| Johnson + Dijkstra | 2.779 | 143.54 | 140,441 (43.7%) | 97.4x | 20/20 | 0.2 |
| Johnson + ALT A* (16 landmarks) | 13.427 | 26.85 | 9,632 (3.0%) | 520.7x | 20/20 | 1.0 |
| Domain-bound A* (no preprocessing) | 0.0 | 70.66 | 58,280 (18.1%) | 197.9x | 20/20 | - |
| Combined max(ALT, domain) A* | 13.427 | 33.27 | 9,616 (3.0%) | 420.3x | 20/20 | 1.0 |
| Bidirectional max(alt,domain), 1 core | 13.427 | 22.8 | 4,187 (1.3%) | 613.3x | 20/20 | 1.0 |
| Bidirectional max(alt,domain), 2 cores | 13.427 | 12.27 | 4,346 (1.4%) | 1139.5x | 20/20 | 1.0 |

## DIMACS BAY, random potential shift

| Method | Preprocess (s) | Query (ms) | Nodes expanded | Speedup vs BF | Correct | Break-even queries |
|---|---:|---:|---:|---:|---:|---:|
| Bellman-Ford (SPFA), no preprocessing | 0.0 | 14500.32 | 23,611,873 (7349.5%) | 1.0x | 20/20 | - |
| Johnson + Dijkstra | 0.211 | 154.43 | 145,233 (45.2%) | 93.9x | 20/20 | 0.0 |
| Johnson + ALT A* (16 landmarks) | 11.219 | 36.58 | 11,641 (3.6%) | 396.4x | 20/20 | 0.8 |
| Bidirectional alt, 1 core | 11.219 | 37.71 | 8,542 (2.7%) | 384.6x | 20/20 | 0.8 |
| Bidirectional alt, 2 cores | 11.219 | 13.07 | 6,227 (1.9%) | 1109.0x | 20/20 | 0.8 |
