"""Independent check of the exact low-temperature series -- WITHOUT the finite-lattice method.

The FLM result in `results/series/sc_lt_free_energy.json` underpins the falsification of Zhang's
claimed low-temperature expansion, so it must be verified by a genuinely different route.

Method used here: a 4x4 x c slab, PERIODIC in the two cross-section directions and with `+`
(frozen up) boundary conditions on the two z faces.  Write

    Xi(a,b,c) = sum_{S subset of the slab} x^{|boundary S|}          (exact integers)

With `+` boundaries the droplet weights are position independent, so the cluster expansion gives

    log Xi(4,4,c) = sum_{clusters} (c - h + 1)^+ W_h ,   h = z-extent of the cluster,

hence the FIRST difference in c isolates the bulk value:

    D(c) := log Xi(4,4,c) - log Xi(4,4,c-1) = sum_{h <= c} W_h = 16 * (phi_LT - 3K) + error.

Two error sources, both controlled:
  * clusters with z-extent h > c.  A connected droplet spanning h layers has surface area at
    least 4h + 2 (a straight column of h cells), so the omitted clusters start at order
    x^{4(c+1)+2}.  For c = 3 that is x^{18}.
  * clusters that wrap the periodic 4x4 cross-section.  The cheapest wrapping object is a
    1x1 tube of length 4 around a periodic direction, of surface area 4 * 4 = 16, so wrapping
    contamination starts at order x^{16}.

Therefore  D(3)/16  equals the bulk low-temperature series  phi_LT - 3K  exactly through x^15.
This route shares NO code with the finite-lattice method: no inclusion-exclusion, no box set, a
different boundary condition and a different lattice geometry.
"""

from __future__ import annotations

import json
from fractions import Fraction

from ising.transfer_matrix import box_broken_bond_poly

ORDER = 15


def series_log(coeffs, order):
    """log(1 + t) composition, exact rationals.  `coeffs[0]` must be 1."""
    assert coeffs[0] == 1, f"Xi(0) = {coeffs[0]} != 1 (the all-up state must be unique)"
    t = [Fraction(c) for c in coeffs[: order + 1]] + [Fraction(0)] * max(0, order + 1 - len(coeffs))
    t[0] = Fraction(0)
    out = [Fraction(0)] * (order + 1)
    power = [Fraction(0)] * (order + 1)
    power[0] = Fraction(1)
    for k in range(1, order + 1):
        new = [Fraction(0)] * (order + 1)
        for i, pi in enumerate(power):
            if pi:
                for j, tj in enumerate(t):
                    if tj and i + j <= order:
                        new[i + j] += pi * tj
        power = new
        sign = Fraction((-1) ** (k + 1), k)
        for i in range(order + 1):
            if power[i]:
                out[i] += sign * power[i]
    return out


def fmt(fr):
    return str(fr) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"


if __name__ == "__main__":
    print("Independent low-temperature series via 4x4xc slabs (periodic cross-section, + faces)")
    logs = {}
    for c in (2, 3):
        coeffs = box_broken_bond_poly((4, 4, c), periodic=(True, True), plus_boundary=True)
        print(f"  Xi(4,4,{c}): {c*16} sites, degree {len(coeffs)-1}, "
              f"first coefficients {coeffs[:8]} ... sum = 2^{c*16}: {sum(coeffs) == 2**(c*16)}")
        logs[c] = series_log(coeffs, ORDER)
    D = [logs[3][i] - logs[2][i] for i in range(ORDER + 1)]
    per_site = [d / 16 for d in D]

    print("\nbulk low-temperature series  phi_LT - 3K  (this route, valid through x^15):")
    terms = [(k, per_site[k]) for k in range(ORDER + 1) if per_site[k]]
    print("   " + "  ".join(f"{fmt(cf)} x^{k}" for k, cf in terms))

    ref = None
    try:
        ref_json = json.load(open("results/series/sc_lt_free_energy.json"))
        d = ref_json.get("data", ref_json)
        for key in ("coefficients", "series", "sc_lt_minus_3K", "phi_lt_minus_3K"):
            if key in d:
                ref = d[key]
                break
        if ref is None:
            for v in d.values():
                if isinstance(v, dict) and all(str(k).isdigit() for k in v):
                    ref = v
                    break
    except Exception as e:
        print(f"\n(could not load FLM reference: {e})")

    ok = True
    if ref is not None:
        print("\ncomparison with the finite-lattice-method series:")
        refd = {int(k): Fraction(str(v)) for k, v in (ref.items() if isinstance(ref, dict)
                                                      else enumerate(ref))}
        for k in range(ORDER + 1):
            a = per_site[k]
            b = refd.get(k, Fraction(0))
            if a or b:
                same = (a == b)
                ok &= same
                print(f"   x^{k:2d}:  slab {fmt(a):>12s}   FLM {fmt(b):>12s}   "
                      f"{'MATCH' if same else 'DIFFER'}")
        print(f"\n{'PASS' if ok else 'FAIL'}: two independent methods agree through x^{ORDER}"
              if ok else f"\nFAIL: methods disagree -- investigate")
    else:
        print("\nFLM reference series not found in the expected shape; printing this route only.")

    out = dict(
        provenance=dict(script="experiments/e26_lt_series_independent.py",
                        method="4x4xc slab, periodic cross-section, plus boundary, first "
                               "difference in c; no finite-lattice method",
                        valid_through=f"x^{ORDER}",
                        error_analysis="omitted z-extent clusters start at x^18; "
                                       "cross-section wrapping starts at x^16"),
        data=dict(series={str(k): fmt(per_site[k]) for k in range(ORDER + 1) if per_site[k]}),
        checks=[dict(name="Xi(0)=1 for both slabs", passed=True,
                     detail="all-up state is the unique zero-energy configuration"),
                dict(name="agrees with FLM series", passed=bool(ok and ref is not None),
                     detail=f"through x^{ORDER}")],
    )
    json.dump(out, open("results/series/lt_independent_check.json", "w"), indent=1)
    print("written results/series/lt_independent_check.json")
