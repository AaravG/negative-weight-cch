"""Summarize results/cpp_<graph>.txt (output of cpp/cch.exe) as Markdown tables."""
import re
import sys
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"
PAT = {
    "custom": r"customization \(no potential\): ([\d.]+) ms",
    "scan": r"negative-cycle scan: ([\d.]+) ms",
    "et": r"elimination-tree query \(no potential\): ([\d.]+) ms",
    "pcch": r"potential read off the CCH \(two sweeps\): ([\d.]+) ms",
    "spfa": r"Johnson potential \(FIFO Bellman-Ford\): ([\d.]+) ms",
    "pops": r"([\d.]+) queue pops per vertex",
    "verified": r"verified against Dijkstra with Johnson potential: (\d+/\d+)",
    "ref": r"correct \(reference ([\d.]+) ms per query\)",
    "h_custom": r"height-shifted: customization ([\d.]+) ms",
    "h_et": r"height-shifted \+ elimination-tree query: ([\d.]+) ms",
    "h_dij": r"height-shifted \+ Dijkstra-based CCH query: ([\d.]+) ms",
    "h_stall": r"height-shifted \+ Dijkstra-based CCH query \+ stall-on-demand: ([\d.]+) ms",
    "j_custom": r"Johnson-shifted: customization ([\d.]+) ms",
    "j_et": r"Johnson-shifted \+ elimination-tree query: ([\d.]+) ms",
    "j_dij": r"Johnson-shifted \+ Dijkstra-based CCH query: ([\d.]+) ms",
    "j_stall": r"Johnson-shifted \+ Dijkstra-based CCH query \+ stall-on-demand: ([\d.]+) ms",
    "par": r"parallel customization \(\d+ threads, \d+ levels\): ([\d.]+) ms",
    "speedup": r"speed-up ([\d.]+)x",
    "pardiff": r"arcs differing from serial \(exact comparison\): (\d+)",
    "negfrac": r"negative arcs: \d+ \(([\d.]+)%\)",
}


def parse(name):
    out = {}
    for line in (RESULTS / f"cpp_{name}.txt").read_text().splitlines():
        m = re.match(r"\[(\w+)\] (.*)", line)
        if not m:
            continue
        metric, rest = m.groups()
        for key, pat in PAT.items():
            mm = re.search(pat, rest)
            if mm and not (key == "speedup" and "parallel" not in rest):
                out.setdefault(metric, {})[key] = mm.group(1)
    return out


if __name__ == "__main__":
    for name in sys.argv[1:] or ["NY", "BAY", "USA"]:
        try:
            data = parse(name)
        except FileNotFoundError:
            continue
        for metric, d in data.items():
            print(name, metric, d)
