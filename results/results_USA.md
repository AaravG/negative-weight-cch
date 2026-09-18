# USA road graph results

Started Wed Sep 16 18:17:11 2026. Updated Wed Sep 16 18:22:08 2026.

## ev weighting

n=23,947,347, m=58,333,344, negative edges=10%. Johnson Bellman-Ford: 34.2 s. Landmarks: 93.7 s wall clock (526.0 s summed).

| Method | Queries done | Mean query (s) | Mean nodes expanded | Correct | Speedup vs Johnson+Dijkstra |
|---|---:|---:|---:|---:|---:|
| Domain-bound A* (no preprocessing) | 10 | 16.412 | 9,545,337 | 10/10 | 1.0x |
| Normal A* (straight-line h, reopening) | 10 | 18.021 | 10,421,359 | 5/10 (3 gave up) | 0.9x |
| Johnson + Dijkstra | 10 | 16.983 | 12,832,909 | 10/10 | 1.0x |
| Johnson + ALT A* (8 landmarks) | 10 | 5.685 | 1,330,732 | 10/10 | 3.0x |
| Combined max(ALT, domain) A* | 10 | 5.677 | 1,329,583 | 10/10 | 3.0x |
| Bidirectional max(ALT, domain), 1 core | 10 | 4.004 | 499,533 | 10/10 | 4.2x |
| Bidirectional, 2 cores | 10 | 1.228 | 515,670 | 10/10 | 13.8x |

## shifted weighting

n=23,947,347, m=58,333,344, negative edges=20%. Johnson Bellman-Ford: 17.8 s. Landmarks: 61.0 s wall clock (520.4 s summed).

| Method | Queries done | Mean query (s) | Mean nodes expanded | Correct | Speedup vs Johnson+Dijkstra |
|---|---:|---:|---:|---:|---:|
| Johnson + Dijkstra | 10 | 16.212 | 12,670,325 | 10/10 | 1.0x |
| Johnson + ALT A* (8 landmarks) | 10 | 3.397 | 916,938 | 10/10 | 4.8x |
| Bidirectional ALT, 1 core | 10 | 2.680 | 410,556 | 10/10 | 6.1x |
| Bidirectional, 2 cores | 10 | 0.961 | 423,368 | 10/10 | 16.9x |
