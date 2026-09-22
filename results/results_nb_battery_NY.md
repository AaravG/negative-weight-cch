# Battery-constrained CCH (Numba) on NY

n = 264,346; shortcut arcs = 1,531,232 (3,062,464 directions); D = 674,674. Every answer is checked against a label-correcting max-charge search.

| Battery M | Customization | Pool pieces | Memory (GB) | Size 0 | Size 1 | 2–5 | 6–10 | 11–50 | 51–100 | >100 | Largest | Query (median) | Infeasible | Correct |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.25 D | 7.0 s | 1,827,148 | 0.4 | 1.516% | 94.601% | 3.401% | 0.320% | 0.148% | 0.010% | 0.004% | 289 | 0.12 ms | 55/60 | 60/60 |
| 1.0 D | 4.4 s | 2,042,804 | 0.4 | 0.050% | 95.943% | 3.483% | 0.349% | 0.161% | 0.010% | 0.004% | 289 | 0.26 ms | 35/60 | 60/60 |
| 4.0 D | 4.8 s | 2,048,304 | 0.4 | 0.000% | 95.990% | 3.485% | 0.350% | 0.161% | 0.010% | 0.004% | 289 | 0.45 ms | 2/60 | 60/60 |

Size columns give the share of shortcut directions with that many pieces.