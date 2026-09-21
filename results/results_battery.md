# Battery-constrained CCH without potentials

Energy metric: DIMACS lengths + synthetic terrain (EV terrain weighting). D = median unconstrained energy of random routes. Every CCH answer is compared with a label-correcting max-SoC search.

## NY

Memory (working set, GB): graph 0.35; after CCH build 0.48 (peak 0.51); after scalar customization + update index 0.67 (peak 0.83). CCH build 258 s. D = 511,678.

| Battery M | Customization | Profile size max / mean | Memory after (GB) | Query (b0 = M) | Query (b0 = M/2) | Infeasible | Correct |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.25 D | 84.0 s | 299 / 1.09 | 1.37 | 2.06 ms | 1.33 ms | 18/20 | 20/20 |
| 1.0 D | 97.9 s | 299 / 1.11 | 1.38 | 4.65 ms | 3.64 ms | 11/20 | 20/20 |
| 4.0 D | 98.9 s | 299 / 1.12 | 1.38 | 6.30 ms | 5.75 ms | 1/20 | 20/20 |

## BAY

Memory (working set, GB): graph 0.79; after CCH build 0.54 (peak 1.43); after scalar customization + update index 0.66 (peak 1.43). CCH build 190 s. D = 1,351,619.

| Battery M | Customization | Profile size max / mean | Memory after (GB) | Query (b0 = M) | Query (b0 = M/2) | Infeasible | Correct |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.25 D | 34.1 s | 222 / 1.07 | 1.12 | 1.24 ms | 1.06 ms | 15/20 | 20/20 |
| 1.0 D | 35.9 s | 222 / 1.07 | 1.12 | 1.85 ms | 1.57 ms | 12/20 | 20/20 |
| 4.0 D | 36.4 s | 222 / 1.07 | 1.12 | 2.24 ms | 2.23 ms | 1/20 | 20/20 |
