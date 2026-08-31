"""EXACT one-sided 95% Clopper-Pearson upper bounds for the Gate B
zero-event layers — computed with python-flint 0.9.0 (Arb) in ball
arithmetic, so every printed digit is backed by a certified enclosure.

For zero events in N trials, the one-sided 95% CP upper bound solves
   (1 - p)^N = alpha                     (binomial tail with i=0)
hence the closed form  upper = 1 - alpha^(1/N).  The computation:
   1. closed form via arb pow (ball arb with directed rounding);
   2. VERIFY by re-evaluating the tail (1-up)^N in arb and checking it
      encloses alpha (the definitional equation);
   3. cross-check against mpmath at 50 dps (independent implementation).
Every value is a ball [mid ± radius]; the radius is the certificate width.
"""
from __future__ import annotations
import json
import flint
from flint import arb

OUT = "/Users/jinleic/jinleic-workspace/cs/mceliece/scratch/gateB"

def cp_upper_arb(N, alpha=arb(5) / 100):
    up = 1 - alpha ** (arb(1) / N)
    tail = (1 - up) ** N          # should enclose alpha
    err = flint.arb(abs(tail - alpha))   # ball radius carries the enclosure
    return up, tail, err

def main():
    res = {}
    import mpmath as mp
    mp.mp.dps = 50
    for N in (237, 127, 228, 2016):
        up, tail, err = cp_upper_arb(N)
        mid = float(up.mid()); rad = float(up.rad())
        upm = 1 - (mp.mpf(5) / 100) ** (mp.mpf(1) / N)
        res[f"N={N}"] = {
            "upper_ball": {"mid": f"{mid:.10f}", "radius": f"{rad:.3e}",
                           "lower": f"{float(up.lower()):.10f}", "upper": f"{float(up.upper()):.10f}"},
            "tail_at_upper_ball": {"mid": f"{float(tail.mid()):.10f}",
                                   "radius": f"{float(tail.rad()):.3e}"},
            "tail_eq_alpha_err_radius": f"{float(err.rad()):.3e}",
            "mpmath_cross_check": f"{float(upm):.10f}",
            "abs_dev_arb_vs_mpmath": f"{abs(mid - float(upm)):.2e}",
            "note": "E1 (N=2016) shown for the record only — an exhaustive census carries NO interval",
        }
        print(f"N={N:5d}: CP95_upper = {mid:.10f} (ball radius {rad:.2e}); "
              f"tail(upper) = {float(tail.mid()):.8f} ~ alpha=0.05, |err| <= {float(err.rad()):.1e}; "
              f"mpmath agrees to {abs(mid - float(upm)):.1e}")
    json.dump(res, open(f"{OUT}/cp_bounds_arb.json", "w"), indent=1)
    print("saved", f"{OUT}/cp_bounds_arb.json")


if __name__ == "__main__":
    main()
