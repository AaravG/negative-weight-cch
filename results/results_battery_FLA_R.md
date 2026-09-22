# Battery-constrained CCH without potentials

Energy metric: DIMACS lengths + synthetic terrain (EV terrain weighting). D = median unconstrained energy of random routes. Every CCH answer is compared with a label-correcting max-SoC search.

## FLA_R

Memory (working set, GB): graph 0.98; after CCH build 1.19 (peak 1.80); after scalar customization + update index 1.56 (peak 1.81). CCH build 1482 s. D = 3,575,138.

| Battery M | Customization | Profile size max / mean | Memory after (GB) | Query (b0 = M) | Query (b0 = M/2) | Infeasible | Correct |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.25 D | 42.3 s | 144 / 1.02 | 3.27 | 2.28 ms | 1.82 ms | 20/20 | 20/20 |
| 1.0 D | 46.1 s | 144 / 1.03 | 3.27 | 5.93 ms | 4.00 ms | 11/20 | 20/20 |
| 4.0 D | 46.0 s | 144 / 1.03 | 3.27 | 7.25 ms | 6.85 ms | 1/20 | 20/20 |
