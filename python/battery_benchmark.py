"""Battery-constrained, potential-free CCH on DIMACS graphs (+ memory of all phases).

    python battery_benchmark.py NY BAY

Writes results_battery.md.
"""
import random
import statistics
import sys
import time
from pathlib import Path

import dimacs
from cch import CCH
from cch_battery import NEG, BatteryCCH, battery_label_correcting
from cch_memory import mem_gb

HERE = Path(__file__).resolve().parent.parent
OUT = None
PAIRS = 10


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    global OUT
    names = sys.argv[1:] or ["NY", "BAY"]
    OUT = HERE / "results" / ("results_battery.md" if names == ["NY", "BAY"] else f"results_battery_{'_'.join(names)}.md")
    lines = ["# Battery-constrained CCH without potentials", "",
             "Energy metric: DIMACS lengths + synthetic terrain (EV terrain weighting). "
             "D = median unconstrained energy of random routes. Every CCH answer is "
             "compared with a label-correcting max-SoC search.", ""]
    for name in names:
        g = dimacs.load(name, "ev", seed=1)
        m0 = mem_gb()
        t0 = time.perf_counter()
        c = CCH(g, leaf=64, log=log)
        t_build = time.perf_counter() - t0
        m1 = mem_gb()
        c.customize(g)
        c.build_update_index()
        m2 = mem_gb()
        rng = random.Random(21)
        dists = []
        while len(dists) < 50:
            d, _ = c.query(rng.randrange(g.n), rng.randrange(g.n))
            if d < float("inf"):
                dists.append(d)
        D = statistics.median(dists)
        lines += [f"## {name}", "",
                  f"Memory (working set, GB): graph {m0[0]:.2f}; after CCH build {m1[0]:.2f} "
                  f"(peak {m1[1]:.2f}); after scalar customization + update index {m2[0]:.2f} "
                  f"(peak {m2[1]:.2f}). CCH build {t_build:.0f} s. D = {D:,.0f}.", "",
                  "| Battery M | Customization | Profile size max / mean | Memory after (GB) | "
                  "Query (b0 = M) | Query (b0 = M/2) | Infeasible | Correct |",
                  "|---:|---:|---:|---:|---:|---:|---:|---:|"]
        pairs = [(rng.randrange(g.n), rng.randrange(g.n)) for _ in range(PAIRS)]
        for factor in (0.25, 1.0, 4.0):
            M = factor * D
            bc = BatteryCCH(c)
            t0 = time.perf_counter()
            biggest, avg = bc.customize(g, M)
            t_cust = time.perf_counter() - t0
            m3 = mem_gb()
            times = {1.0: [], 0.5: []}
            ok = infeasible = total = 0
            for s, t in pairs:
                for frac in (1.0, 0.5):
                    b0 = frac * M
                    t1 = time.perf_counter()
                    got = bc.query(s, t, b0)
                    times[frac].append(time.perf_counter() - t1)
                    want = battery_label_correcting(g, s, b0, M)[t]
                    ok += (got == want == NEG) or abs(got - want) < 1e-6
                    infeasible += want == NEG
                    total += 1
            lines.append(f"| {factor} D | {t_cust:.1f} s | {biggest} / {avg:.2f} | {m3[0]:.2f} | "
                         f"{statistics.mean(times[1.0]) * 1000:.2f} ms | "
                         f"{statistics.mean(times[0.5]) * 1000:.2f} ms | {infeasible}/{total} | "
                         f"{ok}/{total} |")
            log(f"{name} M={factor}D: customize {t_cust:.1f}s, profiles max {biggest} mean {avg:.2f}, "
                f"correct {ok}/{total}, infeasible {infeasible}")
            del bc
            OUT.write_text("\n".join(lines), encoding="utf-8")
        lines.append("")
        OUT.write_text("\n".join(lines), encoding="utf-8")
    log("done")


if __name__ == "__main__":
    main()
