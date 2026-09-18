"""Compact, memory-mapped graphs for the full USA road network (24M nodes, 58M arcs).

``convert()`` parses the DIMACS files once and writes flat binary arrays to
data/cache/. ``CSR(weighting)`` memory-maps them, so every worker process
shares the same OS page cache instead of holding its own copy.

Terrain for the ``ev`` weighting is 2-octave value noise (O(n)), because the
bump-sum terrain used for small graphs would be far too slow at this size.
"""
import gzip
import json
import math
import os
import mmap
import random
import time
from array import array
from pathlib import Path

from graphs import BETA_DOWN, BETA_UP

DATA = Path(__file__).parent / "data"
NAME = os.environ.get("ROAD_GRAPH", "USA")
CACHE = DATA / "cache" / NAME


def _save(arr, fname):
    with open(CACHE / fname, "wb") as f:
        arr.tofile(f)


def _load_array(fname, fmt):
    a = array(fmt)
    with open(CACHE / fname, "rb") as f:
        a.frombytes(f.read())
    return a


def _counting_order(keys, n):
    off = array("i", bytes(4 * (n + 1)))
    for k in keys:
        off[k + 1] += 1
    for i in range(n):
        off[i + 1] += off[i]
    pos = array("i", off[:-1])
    order = array("i", bytes(4 * len(keys)))
    for i, k in enumerate(keys):
        order[pos[k]] = i
        pos[k] += 1
    return off, order


def convert(log=print):
    CACHE.mkdir(parents=True, exist_ok=True)
    if (CACHE / "topology.json").exists():
        return
    t0 = time.time()
    with gzip.open(DATA / f"USA-road-d.{NAME}.co.gz", "rb") as f:
        for line in f:
            if line[:1] == b"v":
                _, i, x, y = line.split()
                i = int(i) - 1
                lat = int(y) * 1e-6
                xm[i] = int(x) * 1e-6 * 111320 * math.cos(math.radians(lat))
                ym[i] = lat * 110540
            elif line[:1] == b"p":
                n = int(line.split()[-1])
                xm, ym = array("d", bytes(8 * n)), array("d", bytes(8 * n))
    log(f"  coordinates parsed ({time.time() - t0:.0f}s)")
    with gzip.open(DATA / f"USA-road-d.{NAME}.gr.gz", "rb") as f:
        k = 0
        for line in f:
            if line[:1] == b"a":
                _, u, v, w = line.split()
                us[k] = int(u) - 1
                vs[k] = int(v) - 1
                ls[k] = int(w)
                k += 1
            elif line[:1] == b"p":
                m = int(line.split()[-1])
                us, vs = array("i", bytes(4 * m)), array("i", bytes(4 * m))
                ls = array("i", bytes(4 * m))
    log(f"  arcs parsed: n={n:,} m={m:,} ({time.time() - t0:.0f}s)")
    off, ford = _counting_order(us, n)
    _save(off, "off.bin")
    _save(array("i", (vs[i] for i in ford)), "tgt.bin")
    roff, rord = _counting_order(vs, n)
    _save(roff, "roff.bin")
    _save(array("i", (us[i] for i in rord)), "rsrc.bin")
    for a, fn in ((xm, "xm.bin"), (ym, "ym.bin"), (us, "arc_u.bin"), (vs, "arc_v.bin"),
                  (ls, "arc_len.bin"), (ford, "ford.bin"), (rord, "rord.bin")):
        _save(a, fn)
    (CACHE / "topology.json").write_text(json.dumps({"n": n, "m": m}))
    log(f"  CSR built ({time.time() - t0:.0f}s)")


def _value_noise(xm, ym, cell, seed):
    rng = random.Random(seed)
    x0, y0 = min(xm), min(ym)
    cols = int((max(xm) - x0) / cell) + 2
    rows = int((max(ym) - y0) / cell) + 2
    lat = array("d", (rng.uniform(-1, 1) for _ in range(cols * rows)))
    out = array("d", bytes(8 * len(xm)))
    for i in range(len(xm)):
        gx, gy = (xm[i] - x0) / cell, (ym[i] - y0) / cell
        cx, cy = int(gx), int(gy)
        fx, fy = gx - cx, gy - cy
        fx, fy = fx * fx * (3 - 2 * fx), fy * fy * (3 - 2 * fy)
        b = cy * cols + cx
        top = lat[b] + (lat[b + 1] - lat[b]) * fx
        bot = lat[b + cols] + (lat[b + cols + 1] - lat[b + cols]) * fx
        out[i] = top + (bot - top) * fy
    return out


def make_weights(weighting, seed=1, log=print):
    meta_path = CACHE / f"{weighting}.json"
    if meta_path.exists():
        return
    t0 = time.time()
    us, vs = _load_array("arc_u.bin", "i"), _load_array("arc_v.bin", "i")
    ls = _load_array("arc_len.bin", "i")
    xm, ym = _load_array("xm.bin", "d"), _load_array("ym.bin", "d")
    m = len(us)
    rng = random.Random(seed)
    meta = {}
    if weighting == "shifted":
        sample = sorted(ls[rng.randrange(m)] for _ in range(200_000))
        spread = 3 * sample[len(sample) // 2]
        phi = array("d", (rng.uniform(0, spread) for _ in range(len(xm))))
        W = array("d", (ls[i] + phi[us[i]] - phi[vs[i]] for i in range(m)))
    else:
        ratios = []
        for _ in range(1_000_000):
            i = rng.randrange(m)
            d = math.hypot(xm[us[i]] - xm[vs[i]], ym[us[i]] - ym[vs[i]])
            if d > 1:
                ratios.append(ls[i] / d)
        ratios.sort()
        scale = ratios[len(ratios) // 100]
        z = _value_noise(xm, ym, 5000, seed)
        z2 = _value_noise(xm, ym, 1500, seed + 1)
        for i in range(len(z)):
            z[i] += 0.3 * z2[i]
        slopes = []
        for _ in range(1_000_000):
            i = rng.randrange(m)
            dz = z[vs[i]] - z[us[i]]
            if dz < 0 and ls[i] > 0:
                slopes.append(BETA_DOWN * -dz / ls[i])
        slopes.sort()
        dz_scale = 1 / slopes[int(len(slopes) * 0.8)]
        for i in range(len(z)):
            z[i] *= dz_scale
        _save(z, "ev_z.bin")
        W = array("d", bytes(8 * m))
        lengthened = 0
        for i in range(m):
            u, v = us[i], vs[i]
            e = scale * math.hypot(xm[u] - xm[v], ym[u] - ym[v])
            L = ls[i]
            if L < e:
                L = e
                lengthened += 1
            dz = z[v] - z[u]
            W[i] = L + (BETA_UP if dz > 0 else BETA_DOWN) * dz
        meta.update(scale=scale, dz_scale=dz_scale, lengthened=lengthened)
    meta["negative_fraction"] = sum(1 for x in W if x < 0) / m
    ford, rord = _load_array("ford.bin", "i"), _load_array("rord.bin", "i")
    _save(array("d", (W[i] for i in ford)), f"{weighting}_fw.bin")
    _save(array("d", (W[i] for i in rord)), f"{weighting}_rw.bin")
    meta_path.write_text(json.dumps(meta))
    log(f"  {weighting} weights built ({time.time() - t0:.0f}s): {meta}")


def _map(fname, fmt):
    f = open(CACHE / fname, "rb")
    mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
    return memoryview(mm).cast(fmt)


class CSR:
    def __init__(self, weighting):
        topo = json.loads((CACHE / "topology.json").read_text())
        self.n, self.m = topo["n"], topo["m"]
        self.weighting = weighting
        self.meta = json.loads((CACHE / f"{weighting}.json").read_text())
        self.off, self.tgt = _map("off.bin", "i"), _map("tgt.bin", "i")
        self.roff, self.rsrc = _map("roff.bin", "i"), _map("rsrc.bin", "i")
        self.fw, self.rw = _map(f"{weighting}_fw.bin", "d"), _map(f"{weighting}_rw.bin", "d")
        self.xm, self.ym = _map("xm.bin", "d"), _map("ym.bin", "d")
        self.z = _map("ev_z.bin", "d") if weighting == "ev" else None
        self.scale = self.meta.get("scale")
