"""Full-USA benchmark, run as many jobs at once as memory allows.

    python usa_benchmark.py [--workers 20] [--queries 10] [--bf-cap-hours 4]
    (set ROAD_GRAPH=NY to test the same pipeline on the small NY graph)

Stages, scheduled as soon as their inputs exist:
  * Bellman-Ford reference queries (2 per weighting, capped) - start immediately
  * domain-bound A* and normal A* queries                     - start immediately
  * Johnson's Bellman-Ford potential                          - start immediately
  * landmark Dijkstras (8 landmarks x 2 directions)           - after Johnson
  * Johnson+Dijkstra queries (the exact ground truth)         - after Johnson
  * ALT, combined and 1-core bidirectional queries            - after landmarks
  * 2-core bidirectional queries (2 dedicated processes)      - after all of the above

All workers run at below-normal priority so the machine stays usable.
Results go to results_<graph>.md / .json, rewritten after every finished job.
Timings of pool jobs are taken while other jobs run, so they include some
memory-bandwidth contention; node-expansion counts are unaffected.
"""
import argparse
import json
import math
import random
import statistics
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path

import bigcsr
import bigsearch as bs

HERE = Path(__file__).resolve().parent.parent / "results"
EPS = 1e-6
LABELS = {
    "bf": "Bellman-Ford (SPFA), no preprocessing",
    "naive": "Normal A* (straight-line h, reopening)",
    "johnson": "Johnson + Dijkstra",
    "alt": f"Johnson + ALT A* ({bs.N_LANDMARKS} landmarks)",
    "domain": "Domain-bound A* (no preprocessing)",
    "max": "Combined max(ALT, domain) A*",
    "bidir-max": "Bidirectional max(ALT, domain), 1 core",
    "bidir-alt": "Bidirectional ALT, 1 core",
    "par2": "Bidirectional, 2 cores",
}


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)


def pick_landmarks(g):
    """'Planar' selection: in each of 8 angular sectors around the centre,
    the node farthest from the centre."""
    xm, ym = g.xm, g.ym
    cx = sum(xm) / g.n
    cy = sum(ym) / g.n
    best = [(-1.0, 0)] * bs.N_LANDMARKS
    sectors = bs.N_LANDMARKS
    for v in range(g.n):
        dx, dy = xm[v] - cx, ym[v] - cy
        k = int((math.atan2(dy, dx) + math.pi) / (2 * math.pi) * sectors) % sectors
        r = dx * dx + dy * dy
        if r > best[k][0]:
            best[k] = (r, v)
    return [v for _, v in best]


class Run:
    def __init__(self, args):
        self.args = args
        self.name = bigcsr.NAME
        self.out = HERE / f"results_{self.name}"
        self.state = {"graph": self.name, "weightings": {}, "started": time.ctime()}

    def save(self):
        (self.out.with_suffix(".json")).write_text(json.dumps(self.state, indent=1))
        lines = [f"# {self.name} road graph results", "",
                 f"Started {self.state['started']}. Updated {time.ctime()}.", ""]
        for w, ws in self.state["weightings"].items():
            lines += [f"## {w} weighting", "",
                      f"n={ws['n']:,}, m={ws['m']:,}, negative edges={ws['neg']:.0%}. "
                      f"Johnson Bellman-Ford: {ws.get('johnson_s', '...')} s. "
                      f"Landmarks: {ws.get('landmarks_wall_s', '...')} s wall clock "
                      f"({ws.get('landmarks_cpu_s', '...')} s summed).", "",
                      "| Method | Queries done | Mean query (s) | Mean nodes expanded | Correct | Speedup vs Johnson+Dijkstra |",
                      "|---|---:|---:|---:|---:|---:|"]
            truth = ws["results"].get("johnson", {})
            base = statistics.mean(r[2] for r in truth.values()) if truth else None
            for m, res in ws["results"].items():
                if not res:
                    continue
                done = [r for r in res.values() if r[0] is not None]
                gave_up = len(res) - len(done)
                mean_t = statistics.mean(r[2] for r in res.values())
                mean_e = statistics.mean(r[1] for r in res.values())
                checked = [(k, r) for k, r in res.items() if k in truth]
                ok = sum(1 for k, r in checked
                         if r[0] is not None and abs(r[0] - truth[k][0]) < EPS)
                corr = f"{ok}/{len(checked)}" + (f" ({gave_up} gave up)" if gave_up else "")
                sp = f"{base / mean_t:.1f}x" if base else "-"
                lines.append(f"| {LABELS[m]} | {len(res)} | {mean_t:.3f} | {mean_e:,.0f} | {corr} | {sp} |")
            lines.append("")
        (self.out.with_suffix(".md")).write_text("\n".join(lines), encoding="utf-8")

    def main(self):
        a = self.args
        log(f"graph {self.name}: converting (skipped if cached)")
        bigcsr.convert(log=log)
        for w in a.weightings:
            bigcsr.make_weights(w, seed=1, log=log)
        g0 = bs.graph(a.weightings[0])
        rng = random.Random(7)
        pairs = [(rng.randrange(g0.n), rng.randrange(g0.n)) for _ in range(a.queries)]
        landmarks = pick_landmarks(g0)
        log(f"landmarks: {landmarks}")
        for w in a.weightings:
            g = bs.graph(w)
            self.state["weightings"][w] = {"n": g.n, "m": g.m, "neg": g.meta["negative_fraction"],
                                           "meta": g.meta, "results": {}}
        bs.low_priority()

        pool = ProcessPoolExecutor(max_workers=a.workers)
        pending = {}

        def submit(tag, fn, *fa):
            pending[pool.submit(fn, *fa)] = tag

        def submit_queries(w, method):
            self.state["weightings"][w]["results"].setdefault(method, {})
            for s, t in pairs:
                submit(("q", w, method, f"{s}-{t}"), bs.job_query, w, method, s, t)

        for w in a.weightings:
            submit(("johnson", w), bs.job_johnson, w)
            self.state["weightings"][w]["results"]["bf"] = {}
            for s, t in pairs[:2]:
                submit(("bf", w, f"{s}-{t}"), bs.job_bf, w, s, t, a.bf_cap_hours * 3600)
            if w == "ev":
                submit_queries(w, "domain")
                submit_queries(w, "naive")

        lm_left = {}
        lm_t0 = {}
        par_done = set()
        while pending:
            done, _ = wait(pending, return_when=FIRST_COMPLETED)
            for f in done:
                tag = pending.pop(f)
                try:
                    res = f.result()
                except Exception as e:  # keep the rest of the run alive
                    log(f"FAILED {tag}: {e!r}")
                    continue
                kind, w = tag[0], tag[1]
                ws = self.state["weightings"][w]
                if kind == "johnson":
                    ws["johnson_s"] = round(res, 1)
                    log(f"{w}: Johnson potential done in {res:.0f}s")
                    submit_queries(w, "johnson")
                    lm_left[w] = 2 * len(landmarks)
                    lm_t0[w] = time.time()
                    ws["landmarks_cpu_s"] = 0.0
                    for i, L in enumerate(landmarks):
                        for rev in (False, True):
                            submit(("lm", w), bs.job_landmark, w, i, L, rev)
                elif kind == "lm":
                    lm_left[w] -= 1
                    ws["landmarks_cpu_s"] = round(ws["landmarks_cpu_s"] + res, 1)
                    if lm_left[w] == 0:
                        ws["landmarks_wall_s"] = round(time.time() - lm_t0[w], 1)
                        log(f"{w}: landmarks done in {ws['landmarks_wall_s']}s wall")
                        for m in (["alt", "max", "bidir-max"] if w == "ev" else ["alt", "bidir-alt"]):
                            submit_queries(w, m)
                elif kind == "bf":
                    ws["results"]["bf"][tag[2]] = list(res)
                    log(f"{w}: Bellman-Ford query {tag[2]} -> {res[0]} in {res[2]:.0f}s")
                else:
                    ws["results"][tag[2]][tag[3]] = list(res)
                    log(f"{w}: {tag[2]} {tag[3]} -> expanded {res[1]:,} in {res[2]:.1f}s")
                self.save()

            # 2-core stage once only Bellman-Ford references remain for a weighting
            for w in a.weightings:
                busy = any(t[1] == w and t[0] != "bf" for t in pending.values())
                if w in par_done or busy or w not in lm_left or lm_left[w] != 0:
                    continue
                par_done.add(w)
                kind = "max" if w == "ev" else "alt"
                log(f"{w}: 2-core bidirectional stage")
                par = bs.ParallelCSR(w, kind, bs.graph(w).n)
                res = self.state["weightings"][w]["results"].setdefault("par2", {})
                try:
                    par.query(*pairs[0])  # warm-up
                    for s, t in pairs:
                        res[f"{s}-{t}"] = list(par.query(s, t))
                        self.save()
                finally:
                    par.close()
        pool.shutdown()
        self.save()
        log(f"finished; wrote {self.out.name}.md")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--queries", type=int, default=10)
    ap.add_argument("--bf-cap-hours", type=float, default=4)
    ap.add_argument("--weightings", nargs="+", default=["ev", "shifted"])
    Run(ap.parse_args()).main()
