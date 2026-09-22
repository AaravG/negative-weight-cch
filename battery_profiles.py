"""Distribution of battery-profile sizes on DIMACS graphs.

For each graph and battery capacity: a histogram of profile sizes, where the
large profiles sit in the hierarchy, and whether their pieces are genuinely
different or near-duplicates (which would point at floating-point effects
in the dominance test rather than real complexity).

    python battery_profiles.py NY BAY

Writes results/results_battery_profiles.md.
"""
import random
import statistics
import sys
import time
from pathlib import Path

import dimacs
from cch import CCH
from cch_battery import BatteryCCH

HERE = Path(__file__).parent
BUCKETS = [(0, 0), (1, 1), (2, 5), (6, 10), (11, 50), (51, 100), (101, 10**9)]


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def bucket_name(lo, hi):
    if lo == hi:
        return str(lo)
    return f"{lo}–{hi}" if hi < 10**9 else f"> {lo - 1}"


def spread(profile):
    """Relative spread of the pieces: (max-min)/scale for in, cost and out."""
    ins = [f[0] for f in profile]
    costs = [f[1] for f in profile]
    outs = [f[2] for f in profile]
    scale = max(1.0, max(abs(x) for x in ins + costs + outs))
    return tuple((max(v) - min(v)) / scale for v in (ins, costs, outs))


def main():
    lines = ["# Battery profile sizes", "",
             "Profile = the set of non-dominated charge functions stored on one shortcut "
             "(one direction). Empty profiles are shortcut directions with no feasible path "
             "through lower vertices.", ""]
    for name in sys.argv[1:] or ["NY", "BAY"]:
        g = dimacs.load(name, "ev", seed=1)
        c = CCH(g, leaf=64, log=log)
        c.customize(g)
        rng = random.Random(21)
        dists = []
        while len(dists) < 50:
            d, _ = c.query(rng.randrange(g.n), rng.randrange(g.n))
            if d < float("inf"):
                dists.append(d)
        D = statistics.median(dists)
        tail = [0] * c.m
        for x in range(c.n):
            for a in range(c.first[x], c.first[x + 1]):
                tail[a] = x
        lines += [f"## {name} (n={g.n:,}, {c.m:,} shortcut arcs = {2 * c.m:,} directions)", ""]
        for factor in (0.25, 1.0, 4.0):
            M = factor * D
            bc = BatteryCCH(c)
            t0 = time.perf_counter()
            bc.customize(g, M)
            t = time.perf_counter() - t0
            sizes = [(len(p), a, "up") for a, p in enumerate(bc.up)] + \
                    [(len(p), a, "down") for a, p in enumerate(bc.down)]
            n_dir = len(sizes)
            counts = []
            for lo, hi in BUCKETS:
                k = sum(1 for s, _, _ in sizes if lo <= s <= hi)
                counts.append(k)
            sizes.sort(reverse=True)
            top = sizes[:10]
            big = [s for s in sizes if s[0] > 10]
            # rank of the lower endpoint, as a percentile (100% = top of the hierarchy)
            def pct(a):
                return 100.0 * tail[a] / (c.n - 1)
            big_pct = statistics.median(pct(a) for _, a, _ in big) if big else float("nan")
            total_pieces = sum(s for s, _, _ in sizes)
            big_pieces = sum(s for s, _, _ in big)
            lines += [f"### M = {factor} D ({t:.0f} s customization)", "",
                      "| Profile size | Directions | Share |", "|---|---:|---:|"]
            for (lo, hi), k in zip(BUCKETS, counts):
                lines.append(f"| {bucket_name(lo, hi)} | {k:,} | {100 * k / n_dir:.3f}% |")
            lines += ["",
                      f"- Profiles with more than 10 pieces: **{len(big):,}** of {n_dir:,} directions "
                      f"({100 * len(big) / n_dir:.4f}%), holding {100 * big_pieces / total_pieces:.1f}% of all pieces.",
                      f"- Median position of those large profiles in the hierarchy: "
                      f"{big_pct:.1f}% (0% = least important vertices, 100% = top separator).", "",
                      "Ten largest profiles:", "",
                      "| Size | Lower-endpoint rank percentile | Relative spread of in / cost / out |",
                      "|---:|---:|---|"]
            for s, a, d in top:
                prof = (bc.up if d == "up" else bc.down)[a]
                si, sc, so = spread(prof)
                lines.append(f"| {s} | {pct(a):.1f}% | {si:.2e} / {sc:.2e} / {so:.2e} |")
            lines.append("")
            log(f"{name} M={factor}D: >10 pieces: {len(big)} of {n_dir}; largest {top[0][0]}; "
                f"median position {big_pct:.1f}%")
            del bc
            (HERE / "results" / "results_battery_profiles.md").write_text("\n".join(lines), encoding="utf-8")
    log("done")


if __name__ == "__main__":
    main()
