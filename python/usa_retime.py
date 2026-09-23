"""Re-time the per-query methods one at a time, on a warm page cache, so their
times are comparable with the 2-core stage (which runs on a quiet machine).
Distances are re-checked against the Johnson+Dijkstra answers in results_USA.json.

    python usa_retime.py      (ROAD_GRAPH defaults to USA)
"""
import json
import statistics
import time
from pathlib import Path

import bigcsr
import bigsearch as bs

HERE = Path(__file__).resolve().parent.parent / "results"
METHODS = {"ev": ["johnson", "domain", "alt", "max", "bidir-max", "naive"],
           "shifted": ["johnson", "alt", "bidir-alt"]}


def main():
    path = HERE / f"results_{bigcsr.NAME}.json"
    state = json.loads(path.read_text())
    bs.low_priority()
    lines = ["# Clean re-timing (one query at a time, warm cache)", ""]
    for w, methods in METHODS.items():
        ws = state["weightings"][w]
        pairs = [tuple(map(int, k.split("-"))) for k in ws["results"]["johnson"]]
        truth = {k: r[0] for k, r in ws["results"]["johnson"].items()}
        par = ws["results"].get("par2", {})
        # warm the page cache for the landmark tables and graph
        for s, t in pairs[:1]:
            bs.job_query(w, "alt", s, t)
        lines += [f"## {w}", "", "| Method | Mean query (s) | Mean nodes expanded | Correct |",
                  "|---|---:|---:|---:|"]
        for m in methods:
            times, exps, ok = [], [], 0
            for s, t in pairs:
                d, e, dt = bs.job_query(w, m, s, t)
                times.append(dt)
                exps.append(e)
                ok += d is not None and abs(d - truth[f"{s}-{t}"]) < 1e-6
            row = f"| {m} | {statistics.mean(times):.3f} | {statistics.mean(exps):,.0f} | {ok}/{len(pairs)} |"
            print(f"[{time.strftime('%H:%M:%S')}] {w} {row}", flush=True)
            lines.append(row)
        p2 = bs.ParallelCSR(w, "max" if w == "ev" else "alt", ws["n"])
        try:
            p2.query(*pairs[0])
            res = [(p2.query(s, t), f"{s}-{t}") for s, t in pairs]
        finally:
            p2.close()
        ok = sum(abs(r[0] - truth[k]) < 1e-6 for r, k in res)
        row = (f"| 2 cores | {statistics.mean(r[2] for r, _ in res):.3f} | "
               f"{statistics.mean(r[1] for r, _ in res):,.0f} | {ok}/{len(res)} |")
        print(f"[{time.strftime('%H:%M:%S')}] {w} {row}", flush=True)
        lines.append(row)
        lines.append("")
    with open(HERE / f"results_{bigcsr.NAME}_retime.md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print("wrote clean timings", flush=True)


if __name__ == "__main__":
    main()
