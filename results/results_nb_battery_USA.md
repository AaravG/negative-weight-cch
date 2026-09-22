# Battery-constrained CCH (Numba) on USA

n = 23,947,347; shortcut arcs = 96,738,990 (193,477,980 directions); D = 18,919,398. Every answer is checked against a label-correcting max-charge search.

| Battery M | Customization | Pool pieces | Memory (GB) | Size 0 | Size 1 | 2–5 | 6–10 | 11–50 | 51–100 | >100 | Largest | Query (median) | Infeasible | Correct |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.25 D | 486.6 s | 29,949,420 | 10.2 | 2.477% | 95.826% | 1.648% | 0.041% | 0.008% | 0.000% | 0.000% | 180 | 2.45 ms | 19/20 | 20/20 |
| 1.0 D | 832.4 s | 32,233,068 | 10.2 | 0.116% | 98.144% | 1.691% | 0.041% | 0.008% | 0.000% | 0.000% | 180 | 15.09 ms | 16/20 | 20/20 |