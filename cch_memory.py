"""Peak memory, query arc scans vs. the depth bound, and the changed-arcs-only
negative-cycle check, for the DIMACS graphs. Windows only (peak working set).

    python cch_memory.py NY BAY
"""
import ctypes
import random
import statistics
import sys
import time
from ctypes import wintypes

import dimacs
from cch import CCH


class PMC(ctypes.Structure):
    _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]


def mem_gb():
    pmc = PMC()
    pmc.cb = ctypes.sizeof(PMC)
    k32 = ctypes.WinDLL("kernel32", use_last_error=True)
    k32.GetCurrentProcess.restype = wintypes.HANDLE
    k32.K32GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(PMC), wintypes.DWORD]
    k32.K32GetProcessMemoryInfo.restype = wintypes.BOOL
    if not k32.K32GetProcessMemoryInfo(k32.GetCurrentProcess(), ctypes.byref(pmc), pmc.cb):
        raise OSError(ctypes.get_last_error())
    return pmc.WorkingSetSize / 2**30, pmc.PeakWorkingSetSize / 2**30


def main():
    for name in sys.argv[1:] or ["NY", "BAY"]:
        g = dimacs.load(name, "ev", seed=1)
        base = mem_gb()[0]
        c = CCH(g, leaf=64, log=lambda *_: None)
        after_build = mem_gb()
        c.customize(g)
        c.build_update_index()
        after_all = mem_gb()
        rng = random.Random(5)
        scans = [c.query(rng.randrange(g.n), rng.randrange(g.n))[1] for _ in range(200)]
        refs = [(u, i) for u in range(g.n) for i in range(len(g.adj[u]))]
        k, (u, i) = next((k, r) for k, r in enumerate(refs)
                         if any(x == r[0] for x, _ in g.adj[g.adj[r[0]][r[1]][0]]))
        v, w_old = g.adj[u][i]
        g.adj[u][i] = (v, -1e9)
        c.update_edges([(k, -1e9)])
        t0 = time.perf_counter()
        local = c.update_created_negative_cycle()
        t_local = time.perf_counter() - t0
        t0 = time.perf_counter()
        full = c.has_negative_cycle()
        t_full = time.perf_counter() - t0
        n_changed = len(c.last_changed)
        g.adj[u][i] = (v, w_old)
        c.update_edges([(k, w_old)])
        print(f"{name}: graph loaded {base:.2f} GB | after CCH build {after_build[0]:.2f} GB "
              f"(peak {after_build[1]:.2f}) | after customization + update index "
              f"{after_all[0]:.2f} GB (peak {after_all[1]:.2f})")
        print(f"{name}: depth d={c.depth}, bound d(d+1)={c.depth * (c.depth + 1):,}; arcs scanned per "
              f"query mean {statistics.mean(scans):,.0f}, max {max(scans):,}")
        print(f"{name}: negative cycle via changed arcs only: {local} in {t_local * 1000:.2f} ms "
              f"({n_changed:,} arcs) vs full scan: {full} in {t_full * 1000:.0f} ms", flush=True)


if __name__ == "__main__":
    main()
