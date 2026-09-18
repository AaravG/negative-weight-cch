"""Download DIMACS 9th Challenge road graphs (distance graphs + coordinates) into data/.

    python scripts/download_dimacs.py NY BAY        # ~14 MB
    python scripts/download_dimacs.py USA           # ~550 MB

Source: http://www.diag.uniroma1.it/challenge9/download.shtml
"""
import sys
import urllib.request
from pathlib import Path

BASE = "http://www.diag.uniroma1.it/challenge9/data/USA-road-d/"
DATA = Path(__file__).resolve().parent.parent / "data"


def main():
    names = sys.argv[1:] or ["NY", "BAY"]
    DATA.mkdir(exist_ok=True)
    for name in names:
        for kind in ("gr", "co"):
            fname = f"USA-road-d.{name}.{kind}.gz"
            target = DATA / fname
            if target.exists():
                print(f"{fname}: already present")
                continue
            print(f"downloading {fname} ...", flush=True)
            urllib.request.urlretrieve(BASE + fname, target)
            print(f"  {target.stat().st_size / 2**20:.1f} MB")


if __name__ == "__main__":
    main()
