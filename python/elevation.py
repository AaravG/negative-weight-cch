"""Real elevation from the public AWS "Terrain Tiles" dataset (Terrarium PNGs).

Source: https://registry.opendata.aws/terrain-tiles/ (SRTM / USGS derived),
no account needed. Tiles are cached under data/terrain/<z>/<x>/<y>.png.

    python elevation.py NY BAY --zoom 11     # download what those graphs need

Elevation (metres) of a Terrarium pixel: (R * 256 + G + B / 256) - 32768.
PNG decoding uses zlib + the standard PNG filters, so no third-party packages.
"""
import argparse
import math
import struct
import sys
import urllib.error
import urllib.request
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

BASE = "https://s3.amazonaws.com/elevation-tiles-prod/terrarium"
ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "terrain"
TILE = 256


def deg2tile(lon, lat, z):
    """Web-Mercator tile coordinates (float, so the fraction gives the pixel)."""
    n = 2 ** z
    x = (lon + 180.0) / 360.0 * n
    lat = max(min(lat, 85.05112878), -85.05112878)
    r = math.radians(lat)
    y = (1.0 - math.log(math.tan(r) + 1 / math.cos(r)) / math.pi) / 2.0 * n
    return x, y


def _decode_png(raw):
    """Minimal PNG reader for 8-bit RGB/RGBA images; returns (w, h, channels, bytes)."""
    assert raw[:8] == b"\x89PNG\r\n\x1a\n", "not a PNG"
    pos = 8
    idat = bytearray()
    w = h = channels = None
    while pos < len(raw):
        (length,) = struct.unpack(">I", raw[pos:pos + 4])
        ctype = raw[pos + 4:pos + 8]
        data = raw[pos + 8:pos + 8 + length]
        pos += 12 + length
        if ctype == b"IHDR":
            w, h, depth, color = struct.unpack(">IIBB", data[:10])
            assert depth == 8 and color in (2, 6), f"unsupported PNG (depth {depth}, color {color})"
            channels = 3 if color == 2 else 4
        elif ctype == b"IDAT":
            idat += data
        elif ctype == b"IEND":
            break
    buf = zlib.decompress(bytes(idat))
    stride = w * channels
    out = bytearray(h * stride)
    prev = bytearray(stride)
    p = 0
    for row in range(h):
        f = buf[p]
        p += 1
        line = bytearray(buf[p:p + stride])
        p += stride
        if f == 1:
            for i in range(channels, stride):
                line[i] = (line[i] + line[i - channels]) & 0xFF
        elif f == 2:
            for i in range(stride):
                line[i] = (line[i] + prev[i]) & 0xFF
        elif f == 3:
            for i in range(stride):
                left = line[i - channels] if i >= channels else 0
                line[i] = (line[i] + ((left + prev[i]) >> 1)) & 0xFF
        elif f == 4:
            for i in range(stride):
                a = line[i - channels] if i >= channels else 0
                b = prev[i]
                c = prev[i - channels] if i >= channels else 0
                pa, pb, pc = abs(b - c), abs(a - c), abs(a + b - 2 * c)
                pr = a if (pa <= pb and pa <= pc) else (b if pb <= pc else c)
                line[i] = (line[i] + pr) & 0xFF
        out[row * stride:(row + 1) * stride] = line
        prev = line
    return w, h, channels, bytes(out)


_tile_cache = {}


def tile_elevations(z, tx, ty, download=True):
    """256x256 list of metres for one tile (cached in memory and on disk)."""
    key = (z, tx, ty)
    if key in _tile_cache:
        return _tile_cache[key]
    path = CACHE / str(z) / str(tx) / f"{ty}.png"
    if not path.exists():
        if not download:
            return None
        path.parent.mkdir(parents=True, exist_ok=True)
        url = f"{BASE}/{z}/{tx}/{ty}.png"
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                path.write_bytes(r.read())
        except urllib.error.HTTPError as e:
            if e.code == 404:          # ocean / outside coverage
                path.write_bytes(b"")
            else:
                raise
    raw = path.read_bytes()
    if not raw:
        elev = [0.0] * (TILE * TILE)
    else:
        w, h, ch, px = _decode_png(raw)
        elev = [0.0] * (w * h)
        for i in range(w * h):
            r, g, b = px[i * ch], px[i * ch + 1], px[i * ch + 2]
            elev[i] = (r * 256 + g + b / 256.0) - 32768.0
    if len(_tile_cache) > 4000:
        _tile_cache.clear()
    _tile_cache[key] = elev
    return elev


def sample(lonlat, z=11, download=True, log=print):
    """Elevation in metres for a list of (lon, lat), nearest-pixel sampling."""
    needed = set()
    coords = []
    for lon, lat in lonlat:
        fx, fy = deg2tile(lon, lat, z)
        tx, ty = int(fx), int(fy)
        px = min(TILE - 1, int((fx - tx) * TILE))
        py = min(TILE - 1, int((fy - ty) * TILE))
        needed.add((tx, ty))
        coords.append((tx, ty, py * TILE + px))
    log(f"elevation: {len(needed):,} tiles at zoom {z} for {len(coords):,} points")
    if download:
        missing = [(tx, ty) for tx, ty in needed
                   if not (CACHE / str(z) / str(tx) / f"{ty}.png").exists()]
        if missing:
            log(f"downloading {len(missing):,} tiles ...")
            with ThreadPoolExecutor(max_workers=8) as ex:
                list(ex.map(lambda t: tile_elevations(z, t[0], t[1]), missing))
    out = [0.0] * len(coords)
    for i, (tx, ty, idx) in enumerate(coords):
        out[i] = tile_elevations(z, tx, ty, download)[idx]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("graphs", nargs="+")
    ap.add_argument("--zoom", type=int, default=11)
    args = ap.parse_args()
    import dimacs
    for name in args.graphs:
        lonlat = dimacs.raw_lonlat(name)
        z = sample(lonlat, args.zoom)
        lo, hi = min(z), max(z)
        mean = sum(z) / len(z)
        print(f"{name}: {len(z):,} nodes, elevation {lo:.0f} to {hi:.0f} m (mean {mean:.0f} m)")
    total = sum(f.stat().st_size for f in CACHE.rglob("*.png"))
    print(f"tile cache: {total / 2**20:.0f} MB in {CACHE}")


if __name__ == "__main__":
    main()
