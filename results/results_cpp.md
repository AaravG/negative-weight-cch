# C++ implementation: results and comparison with the classical pipeline

`cpp/cch.cpp`, built with MSVC (/O2). Same algorithms as the Python and Numba
versions. Metrics: `ev` (synthetic terrain), `ev_real` (real elevation from the
AWS terrain tiles, see `elevation.py`), `shifted` (random potential shift).
1,000 random queries per metric, 20 of them verified against Dijkstra with a
Johnson potential. Raw output: `cpp_NY.txt`, `cpp_BAY.txt`.

## New York (264,346 nodes) — preprocessing

| Step | Time |
|---|---:|
| Nested-dissection ordering (inertial flow) | 2.8 s |
| Symbolic contraction + arc mapping | < 0.1 s |
| Customization (per metric, includes the negative-cycle check) | 0.04–0.06 s |

1,531,640 shortcut arcs, 10,168,563 triangles, elimination-tree depth 345.

## Query time (ms), and what each approach needs

| Approach | Potential | NY, ev | NY, ev_real | BAY, ev_real |
|---|---|---:|---:|---:|
| **Ours: no potential, elimination-tree query** | **none** | **0.012** | **0.013** | **0.007** |
| Shifted metric + elimination-tree query | Johnson or height | 0.022 | 0.021 | 0.010 |
| Height-potential CCH + Dijkstra query + stalling (Baum/Eisner style) | height (free) | 0.099 | 0.076 | 0.049 |
| Johnson-shifted CCH + Dijkstra query | Johnson (0.02–0.10 s) | 0.100 | 0.071 | 0.046 |
| Johnson + plain Dijkstra | Johnson | 17 | 16 | 13 |

All approaches agree with each other and with the reference (20/20 per metric).

**Reading the table.** The query speed-up comes from the sweep (elimination-tree)
query, not from avoiding the potential: the same query on a shifted metric is
equally fast. What our approach removes is the potential computation itself,
which matters when the cost model has no height-induced potential (Storandt,
personal communication) or when the metric changes often.

## Full USA road graph (23,947,347 nodes, 58,333,344 arcs), C++

Raw output: `cpp_USA.txt`. One-time ordering 1,204 s (20 min); symbolic contraction
1.4 s; 96,734,624 shortcut arcs; 3,151,504,591 triangles (never stored);
elimination-tree depth 3,762.

| | EV terrain | Random shift |
|---|---:|---:|
| Customization, no potential (includes the negative-cycle check) | **5.8 s** | 5.5 s |
| Negative-cycle scan | 0.03 s | 0.03 s |
| **Our query (sweep, no potential)** | **0.62 ms** | **0.64 ms** |
| Shifted metric + sweep query | 0.61 ms | 0.64 ms |
| Height-potential CCH + Dijkstra query + stalling | 2.69 ms | – |
| Johnson-shifted CCH + Dijkstra query + stalling | 2.58 ms | 2.65 ms |
| Johnson + plain Dijkstra (reference) | 1,577 ms | 1,529 ms |
| Correct | 20/20 | 20/20 |

Johnson's potential, which the potential-free variant does not need, costs
2.56 s (EV terrain) and 1.33 s (random shift) on this graph.

## Real versus synthetic elevation

| Graph | Negative edges (synthetic) | Negative edges (real) | Our query, synthetic | real |
|---|---:|---:|---:|---:|
| New York | 10% | 10.0% | 0.012 ms | 0.013 ms |
| SF Bay Area | 10% | 10.3% | 0.007 ms | 0.007 ms |

Real elevation: AWS Terrain Tiles (SRTM/USGS), zoom 11 (about 60 m per pixel),
nearest-pixel sampling, sea-floor values clamped to 0. Climbing one metre is
charged as 50 metres of driving, so descents steeper than about 3% are negative.
