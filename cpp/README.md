# C++ implementation

Single file, C++20, no dependencies (OpenMP optional, for `--parallel`).

## Build

```bat
build.bat                      :: MSVC
```

```bash
g++ -O2 -std=c++20 -fopenmp -o cch cch.cpp      # GCC / Clang
```

## Input

Graphs come from the Python exporter, which writes flat binary files:

```bash
python ../python/export_graph.py NY ../cpp/data/NY
```

`meta.txt` (nodes, arcs, metric names), `off.i64`, `tgt.i32`, `x.f64`, `y.f64`,
`w_<metric>.f64`, and optionally `z_<metric>.f64` (elevation, for the
height-potential baseline).

## Run

```bash
./cch data/NY --queries 1000 --check 20 --shifted-cch
./cch data/USA --order data/USA/order.i32 --parallel 32
```

| Option | Meaning |
|---|---|
| `--order <file>` | read the node order if the file exists, otherwise compute and write it |
| `--queries <n>` | random queries to time (default 1000) |
| `--check <n>` | how many of them to verify against Dijkstra with a Johnson potential |
| `--metric <name>` | run a single metric instead of all |
| `--shifted-cch` | also run the classical baselines: height-potential and Johnson-shifted CCH with Dijkstra queries, with and without stall-on-demand |
| `--parallel <threads>` | also run customization level-parallel and compare it with the serial result |

## What it contains

`cch.cpp` has the whole pipeline: undirected graph construction, nested
dissection with inertial-flow bisection (Dinic max-flow on the node-split
graph), symbolic contraction (elimination tree and shortcut arcs), customization
with triangles enumerated on the fly, the elimination-tree (sweep) query, the
negative-cycle scan, level-parallel customization with atomic minima, and the
baselines listed above.

Results: [`../results/results_cpp.md`](../results/results_cpp.md).
