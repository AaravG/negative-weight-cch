# Potential-free CCH for negative edge weights

Shortest paths on road networks with **negative edge weights**, motivated by energy-optimal routing for electric vehicles, where regenerative braking makes downhill segments cost negative energy.

The core question: can **Customizable Contraction Hierarchies (CCH)** handle negative weights *without* first computing a potential (Johnson / Bellman–Ford reweighting)? This repository contains a working answer, with proofs and experiments, shared as an **open technical log for review**. Whether the observation is already known in the CCH community is not yet confirmed. Pointers to prior work are very welcome — please open an issue.

## Key observation

- **Customization.** CCH processes lower triangles in rank order:

      ℓ⁺(v,w) ← min(ℓ⁺(v,w), ℓ⁺(v,u) + ℓ⁺(u,w))

  This is vertex elimination in the (min, +) algebra, which is exact for any weights *without negative cycles* (Carré 1971, Tarjan 1981).
- **Query.** The elimination-tree query relaxes a small ancestor DAG in rank order. It needs no priority queue and no Dijkstra stopping rule, which are exactly the parts that fail with negative weights.
- **Negative cycles.** A negative cycle exists **iff** some shortcut has `ℓ⁺(v,w) + ℓ⁺(w,v) < 0`. After an update, only the changed shortcuts need checking.
- **Updates.** Partial re-customization recomputes affected shortcuts from their lower triangles, in rank order.

Proofs, complexity and limitations are in [`docs/technical_note.md`](docs/technical_note.md). Every algorithm used here is explained in [`docs/algorithms.md`](docs/algorithms.md).

## Results (pure Python, single-threaded)

All answers are checked against Johnson + Dijkstra, which is itself validated against Bellman–Ford.

**Queries, DIMACS road graphs with 10–20% negative edges**

| Graph | Johnson + Dijkstra | Johnson + ALT A* | **CCH, no potential** |
|---|---:|---:|---:|
| New York (264k nodes) | 117–128 ms | 21–23 ms | **0.77 ms** |
| SF Bay Area (321k nodes) | 183–187 ms | 43–45 ms | **0.30 ms** |

**Preparing a new metric**

| Graph | Johnson's step + 16 ALT landmarks | **CCH customization** |
|---|---:|---:|
| New York | 9.1–12.5 s | **1.6 s** |
| SF Bay Area | 10.5–14.2 s | **0.9–1.0 s** |

The CCH time includes the negative-cycle check.

**Other results**

- **Updates:** a single changed road takes 0.0–0.6 ms; 1,000 changed roads take 0.25–0.65 s. The result matched full re-customization exactly in every test.
- **Full USA graph** (24M nodes, potential-based methods only): Johnson + Dijkstra 9–12 s per query, bidirectional ALT on 2 cores 0.8–1.2 s. Bellman–Ford took more than 35 minutes per query.
- **Textbook A\*** with a straight-line heuristic returns **wrong answers** on about 50% of queries when edges can be negative.

Full tables are in [`results/`](results/). Absolute times are Python times; compiled implementations would be much faster, and the ratios between methods are only indicative.

## Battery constraints (early work)

With a battery of capacity M, each road maps state of charge `b ↦ min(out, b − cost)` for `b ≥ in`. These 3-parameter functions are closed under composition. Together with pointwise max, they satisfy what CCH customization needs, since no loop can gain charge. See `cch_battery.py`.

Validated so far on small graphs only: 15,540 checks against step-by-step simulation and a label-correcting reference search. Large-graph experiments are pending.

## Repository layout

| Files | What |
|---|---|
| `cch.py`, `inertial_flow.py` | CCH with negative weights: ordering (inertial flow + Dinic max-flow), contraction, customization, partial updates, negative-cycle test, elimination-tree query |
| `cch_battery.py` | Battery-constrained customization (profiles of charge functions) and scalar queries |
| `graphs.py`, `dimacs.py` | Synthetic graphs with negative edges; DIMACS loader with EV-terrain and random-shift weightings |
| `preprocess.py`, `heuristics.py`, `search.py`, `parallel.py` | Baselines: Bellman–Ford/SPFA, Johnson, ALT, domain bound, A*, bidirectional A*, 2-process bidirectional A* |
| `bigcsr.py`, `bigsearch.py`, `usa_benchmark.py`, `usa_retime.py` | Memory-mapped pipeline for the full USA graph |
| `test_*.py` | Correctness tests |
| `*benchmark*.py`, `cch_memory.py` | Experiments |
| `docs/`, `results/` | Technical note, algorithm glossary, result tables |

## Running

Requires Python 3.10+ and the standard library only.

```bash
python test_correctness.py        # A*/ALT/bidirectional vs Bellman–Ford
python test_cch.py                # CCH: queries, bounds, updates, negative cycles
python test_cch_battery.py        # battery-constrained CCH
python benchmark.py --quick       # small synthetic benchmark (seconds)

python scripts/download_dimacs.py NY BAY    # ~14 MB into data/
python cch_benchmark.py NY BAY              # CCH vs Johnson/ALT (~20 min)
python battery_benchmark.py NY BAY          # battery-constrained CCH
```

Road data is from the [9th DIMACS Implementation Challenge](http://www.diag.uniroma1.it/challenge9/). It is not redistributed here; the download script fetches it. The elevation used for the EV weighting is **synthetic**. Real elevation data is future work.

## Limitations

- Pure-Python prototype; no C++ implementation yet.
- The inertial-flow ordering is simpler than FlowCutter/KaHIP.
- Synthetic elevation.
- No CCH experiments on the full USA graph yet.
- The battery model has no charging stops.

## Status and feedback

This is research in progress. If you know of prior work, spot an error in the proofs, or work on CH/CCH or EV routing, please open an issue.

The code was developed with the help of an AI coding assistant. All results are checked against independent reference implementations included in this repository.

## License

MIT — see [LICENSE](LICENSE).
