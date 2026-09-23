"""Check that the fast battery customization gives the same profiles as the
reference one, and measure the speed-up.

    python battery_speed_check.py            # small synthetic graphs
    python battery_speed_check.py NY         # also a DIMACS graph
"""
import random
import sys
import time

import graphs
from cch import CCH
from cch_battery import BatteryCCH, eval_profile


def same_profiles(p, q, M, rng):
    """Equal as functions (checked at many charge values) and equal in size."""
    if len(p) != len(q):
        return False
    pts = [0.0, M] + [rng.uniform(0, M) for _ in range(20)]
    pts += [f[0] for f in p] + [f[0] + 1e-7 for f in p]
    return all(abs(eval_profile(p, b) - eval_profile(q, b)) < 1e-6
               or eval_profile(p, b) == eval_profile(q, b) for b in pts if 0 <= b <= M)


def check(name, g, Ms):
    rng = random.Random(4)
    c = CCH(g, leaf=64 if g.n > 5000 else 8, log=lambda *_: None)
    for M in Ms:
        fast, ref = BatteryCCH(c), BatteryCCH(c)
        t0 = time.perf_counter()
        fast.customize(g, M)
        t_fast = time.perf_counter() - t0
        t0 = time.perf_counter()
        ref.customize_reference(g, M)
        t_ref = time.perf_counter() - t0
        bad = sum(not same_profiles(p, q, M, rng) for p, q in zip(fast.up, ref.up))
        bad += sum(not same_profiles(p, q, M, rng) for p, q in zip(fast.down, ref.down))
        print(f"{name} M={M:,.0f}: reference {t_ref:.2f}s, fast {t_fast:.2f}s "
              f"({t_ref / t_fast:.1f}x), differing profiles: {bad}", flush=True)
        assert bad == 0


def main():
    for i in range(3):
        check(f"grid{i}", graphs.ev_grid(15, seed=i), [3.0, 10.0, 40.0])
        check(f"geo{i}", graphs.ev_geometric(300, seed=i), [3.0, 10.0, 40.0])
    for name in sys.argv[1:]:
        import dimacs
        g = dimacs.load(name, "ev", seed=1)
        check(name, g, [511_678.0])
    print("all profiles identical")


if __name__ == "__main__":
    main()
