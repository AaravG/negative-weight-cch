# C++ results

Produced by `cpp/run_all.sh` (medians of 5 repetitions, 1,000 random queries; NY/BAY verify all 1,000 queries, USA 100).
Raw logs: `cpp_NY.txt`, `cpp_BAY.txt`, `cpp_USA.txt`; the USA ordering time (1,204 s) is in `cpp_USA_first_run_with_ordering.txt`.
Tables: `python scripts/summarize_cpp.py`. The paper's Tables 1 and 2 are taken directly from these logs.

Earlier single-run numbers in this file (e.g. 0.62 ms USA query, 5.8 s customization) were measured with a
previous driver and are superseded.
