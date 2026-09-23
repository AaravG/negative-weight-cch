# Python implementation

Two layers:

* **reference** (standard library only): readable, used for proofs, tests and
  the smaller experiments;
* **compiled with Numba** (`numpy`, `numba`): the same algorithms on flat
  arrays, used for continental-scale graphs.

## Files

| File | Purpose |
|---|---|
| `cch.py` | CCH with negative weights: ordering, contraction, customization, partial updates (recompute and tie-based), negative-cycle test (also during elimination), elimination-tree query, path unpacking, perfect customization, witness pruning |
| `inertial_flow.py` | nested-dissection bisection via min vertex cuts (Dinic max-flow) |
| `cch_battery.py` | battery-constrained customization and queries (charge-function profiles) |
| `cch_nb.py`, `cch_nb_battery.py` | Numba versions of both, for large graphs |
| `run_nb.py`, `run_nb_battery.py` | drivers for the Numba versions (NY, regions, full USA) |
| `graphs.py`, `dimacs.py`, `elevation.py` | synthetic graphs; DIMACS loader with the `ev`, `ev_real` and `shifted` metrics; real elevation from AWS terrain tiles |
| `preprocess.py`, `heuristics.py`, `search.py`, `parallel.py` | baselines: Bellman–Ford/SPFA, Johnson, ALT, domain bound, A*, bidirectional A*, 2-process bidirectional A* |
| `bigcsr.py`, `bigsearch.py`, `usa_benchmark.py`, `usa_retime.py` | memory-mapped pipeline for the potential-based methods on the full USA |
| `export_graph.py` | writes the binary input for the C++ implementation |
| `benchmark.py`, `cch_benchmark.py`, `battery_benchmark.py`, `battery_profiles.py`, `battery_speed_check.py`, `cch_memory.py` | experiments |
| `test_*.py` | correctness tests |

## Tests

```bash
python test_correctness.py     # A*/ALT/bidirectional vs Bellman–Ford      (1,380 checks)
python test_cch.py             # CCH queries, bounds, updates, cycles         (996)
python test_cch_battery.py     # battery-constrained CCH                   (15,540)
python test_cch_extras.py      # path unpacking, tie-based updates
python test_cch_pruning.py     # perfect customization, witness pruning, stalling
```

## Experiments

```bash
python ../scripts/download_dimacs.py NY BAY   # road data (~14 MB)
python elevation.py NY BAY                    # real elevation tiles (~23 MB)
python benchmark.py --quick                   # synthetic graphs, seconds
python cch_benchmark.py NY BAY                # CCH vs Johnson/ALT (~20 min)
python battery_benchmark.py NY BAY            # battery-constrained CCH
python run_nb.py USA                          # full USA (Numba); see run_nb.py for USA_CACHE
```

Results are written to [`../results/`](../results).
