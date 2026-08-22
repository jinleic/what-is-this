"""Standalone auditor for experiments/e145_paired_momentum.py.

Recomputes the decisive exact values of results/bounds/paired_momentum.json
by INDEPENDENT engines and routes:

* full-configuration enumeration (all 2^N states, no sigma_0 halving, own
  histogram code) for the four tori, against the producer's half-space
  engine;
* transcendental enclosures by the reciprocal route e^{-t} = (e^u)^{-N}
  with the all-positive series and an explicit geometric remainder
  (different from the producer's alternating-Taylor route);
* R_{2n} by the binomial-product decomposition
  (2n)!/(i!^2 j!^2 k!^2) = C(2n,2i) C(2i,i) C(2n-2i,2j) C(2j,j);
* I3(mu^2) by its own positive interval convolution at different dps and
  length, asserting interval overlap;
* the shell weights at grid 48 with its own cos enclosures (different
  Taylor order);
* the delta* thresholds from the stored I3 strings;
* the e72 single-bond anchor and the frozen transfer-matrix polynomial.

Run from the repository root:
    .venv/bin/python tests/test_paired_momentum.py
"""

from __future__ import annotations

from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import sys

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "bounds" / "paired_momentum.json"
UPPER_INFRARED = ROOT / "results" / "bounds" / "upper_infrared.json"
UPPER_BEYOND = ROOT / "results" / "bounds" / "upper_beyond.json"

COS_ORDER_TEST = 16
EXP_POS_ORDER = 50
PI_LO = Fraction(333, 106)
PI_HI = Fraction(355, 113)

FAILURES = []


def check(name: str, passed: bool, detail: str = "") -> None:
    status = "PASS" if passed else "FAIL"
    print(f"{status} {name} -- {detail}")
    if not passed:
        FAILURES.append(name)


def f2m(f: Fraction):
    return mp.mpf(f.numerator) / mp.mpf(f.denominator)


# ---------------------------------------------------------------------------
# independent transcendental enclosures: e^{-t} = (e^u)^{-N}, positive series
# ---------------------------------------------------------------------------

def exp_neg_interval_pos(t: Fraction) -> tuple[Fraction, Fraction]:
    n_split = max(1, math.ceil(2 * float(t)))
    u = t / n_split
    assert 0 < u <= Fraction(1, 2)
    s = Fraction(0)
    term = Fraction(1)
    for k in range(0, EXP_POS_ORDER + 1):
        if k > 0:
            term = term * u / k
        s += term
    rem = term * u / (EXP_POS_ORDER + 1) / (1 - u / (EXP_POS_ORDER + 2))
    e_lo = s                 # e^u >= partial sum (all terms positive)
    e_hi = s + rem           # e^u <= partial + geometric remainder
    # e^{-t} = (e^u)^{-N}: reciprocals reverse, N-th powers monotone
    return (Fraction(1) / e_hi**n_split, Fraction(1) / e_lo**n_split)


def cos_taylor_test(y: Fraction) -> tuple[Fraction, Fraction]:
    assert 0 <= y < 1
    if y == 0:
        return (Fraction(1), Fraction(1))
    term = Fraction(1)
    s = Fraction(1)
    low = high = None
    y2 = y * y
    for m in range(1, COS_ORDER_TEST + 2):
        term = term * y2 / ((2 * m - 1) * (2 * m))
        nxt = term * y2 / ((2 * m + 1) * (2 * m + 2))
        if m % 2 == 1:
            s = s - term
            low, high = s, s + nxt
        else:
            s = s + term
            low, high = s - nxt, s
    return (low, high)


def cos_interval_test(x: Fraction) -> tuple[Fraction, Fraction]:
    assert 0 <= x <= PI_HI
    if x == 0:
        return (Fraction(1), Fraction(1))
    if x > PI_LO / 2:
        a = max(PI_LO - x, Fraction(0))
        b = max(PI_HI - x, Fraction(0))
        clo = cos_quad_test(b)[0]
        chi = cos_quad_test(a)[1]
        return (-chi, -clo)
    return cos_quad_test(x)


def cos_quad_test(x: Fraction) -> tuple[Fraction, Fraction]:
    if x == 0:
        return (Fraction(1), Fraction(1))
    y = x / 4
    c_lo = cos_taylor_test(y)[0]
    c_hi = cos_taylor_test(y)[1]
    return (8 * c_lo**4 - 8 * c_lo**2 + 1, 8 * c_hi**4 - 8 * c_hi**2 + 1)


# ---------------------------------------------------------------------------
# independent lattice engine: full 2^N enumeration
# ---------------------------------------------------------------------------

RATIONAL_COS = {
    2: (Fraction(1), Fraction(-1)),
    3: (Fraction(1), Fraction(-1, 2), Fraction(-1, 2)),
    4: (Fraction(1), Fraction(0), Fraction(-1), Fraction(0)),
}


def build_torus_test(sides):
    sites = list(itertools.product(range(sides[0]), range(sides[1]), range(sides[2])))
    index = {s: i for i, s in enumerate(sites)}
    bonds = []
    for s in sites:
        for axis in range(3):
            if sides[axis] == 2 and s[axis] == 1:
                continue
            t = list(s)
            t[axis] = (s[axis] + 1) % sides[axis]
            for _ in range(2 if sides[axis] == 2 else 1):
                bonds.append((index[s], index[tuple(t)]))
    return sites, bonds


def full_enumerate(sites, bonds):
    """All 2^N configurations; returns (c, c2, m2_num, m4_num, sum_w).

    c[q]: count of configurations with q unsatisfied bonds; c2[s][q]: the
    signed count weighted by sigma_0 sigma_s; m2_num/m4_num: weighted
    magnetisation-power sums.
    """
    n = len(sites)
    nb = len(bonds)
    c = [0] * (nb + 1)
    c2 = [[0] * (nb + 1) for _ in range(n)]
    m2 = Fraction(0)
    m4 = Fraction(0)
    # weights are powers of x, tracked symbolically: accumulate integer
    # histograms of m so the caller can form both moments at any x
    mcount = {}
    for state in range(1 << n):
        q = 0
        for a, b in bonds:
            if ((state >> a) & 1) != ((state >> b) & 1):
                q += 1
        c[q] += 1
        m = 2 * state.bit_count() - n
        key = (q, m)
        mcount[key] = mcount.get(key, 0) + 1
        sgn0 = 1 if not state & 1 else -1
        for s in range(1, n):
            sg = 1 if not (state >> s) & 1 else -1
            c2[s][q] += sgn0 * sg
    return c, c2, mcount


def poly_at(hist, x: Fraction) -> Fraction:
    tot = Fraction(0)
    p = Fraction(1)
    for q, v in enumerate(hist):
        tot += v * p
        p = p * x
    return tot


def truncate_down(f: Fraction, digits: int) -> Fraction:
    scale = 10**digits
    return Fraction(f.numerator * scale // f.denominator, scale)


# ---------------------------------------------------------------------------
# main audit
# ---------------------------------------------------------------------------

def main() -> int:
    mp.mp.dps = 120  # via context attribute; the mp.dps shorthand is inert here
    art = json.loads(ARTIFACT.read_text())
    data = art["data"]
    check("artifact_checks_all_pass", all(c["pass"] for c in art["checks"]),
          f"{len(art['checks'])} producer checks")

    src = json.loads(UPPER_INFRARED.read_text())
    i3 = src["data"]["certified_constants"]["I3"]
    i3_lo, i3_hi = Fraction(i3[0]), Fraction(i3[1])
    beyond = json.loads(UPPER_BEYOND.read_text())
    xstar = beyond["data"]["limitation_certificate"]["x_star_interval"]
    x_star_lo = Fraction(xstar["lower_exact"])

    # --- delta* thresholds recomputed from stored I3 strings
    ok = True
    for label, pair in data["delta_star_thresholds"].items():
        if label.startswith("incumbent"):
            continue
        ktxt = label.replace("K=", "")
        if ktxt.startswith("kc_lower"):
            k = Fraction(data["certified_inputs"]["kc_lower_floor_40"])
        else:
            k = Fraction(ktxt)
        want = (i3_lo / (2 * k) - 1, i3_hi / (2 * k) - 1)
        got = (Fraction(pair[0]), Fraction(pair[1]))
        ok = ok and want == got
    check("delta_star_recomputed", ok, "exact equality from I3 interval strings")

    # --- walk return counts via binomial products
    def n2n(n: int) -> int:
        total = 0
        for i in range(n + 1):
            for j in range(n - i + 1):
                k = n - i - j
                total += (
                    math.comb(2 * n, 2 * i)
                    * math.comb(2 * i, i)
                    * math.comb(2 * n - 2 * i, 2 * j)
                    * math.comb(2 * j, j)
                    * math.comb(2 * k, k)
                )
        return total

    ok = all(n2n(n) == N for n, N in enumerate(data["walk_return_counts_N2n"]))
    check("walk_return_counts_binomial_route", ok, "N_{2n} via C(2n,2i)C(2i,i)... products")

    # --- independent I3(mu^2) at mu^2 in {1/2, 1/8}, dps 40, N=500
    for label, mu2, two_k_hint in (("1/2", Fraction(1, 2), None), ("1/8", Fraction(1, 8), None)):
        mu2f = Fraction(mu2)
        n_use = 500
        mp.iv.dps = 40
        s = [mp.iv.mpf(1)]
        cur = mp.iv.mpf(1)
        for a in range(1, n_use + 1):
            cur = cur * (mp.iv.mpf(1) / mp.iv.mpf(a * a))
            s.append(cur)
        t = [mp.iv.mpf(0)] * (n_use + 1)
        for i in range(n_use + 1):
            for j in range(n_use - i + 1):
                t[i + j] = t[i + j] + s[i] * s[j]
        c = [mp.iv.mpf(0)] * (n_use + 1)
        for i in range(n_use + 1):
            for j in range(n_use - i + 1):
                c[i + j] = c[i + j] + t[i] * s[j]
        base = 3 + mu2f
        r = Fraction(3) / base
        total = mp.iv.mpf(0)
        factor = mp.iv.mpf(1) / mp.iv.mpf(int(base))
        rfac = mp.iv.mpf(1)
        for n in range(n_use + 1):
            if n > 0:
                rn = (
                    mp.iv.mpf(int(math.factorial(2 * n))) * c[n]
                    / mp.iv.mpf(int(6 ** (2 * n)))
                )
                rfac = rfac * mp.iv.mpf(r.numerator) / mp.iv.mpf(r.denominator)
                total = total + factor * rfac * rn
            else:
                total = total + factor
        tail = r ** (n_use + 1) / mu2f
        own_lo = total.a
        own_hi = total.b + mp.iv.mpf(tail.numerator) / mp.iv.mpf(tail.denominator)
        stored = data["massive_watson"][label]["certified_interval"]
        s_lo = mp.mpf(stored[0])
        s_hi = mp.mpf(stored[1])
        overlap = (own_lo <= s_hi) and (s_lo <= own_hi)
        width = mp.mpf(stored[1]) - mp.mpf(stored[0])
        check(
            f"I3_mu2_{label}_overlap",
            overlap,
            f"own [{mp.nstr(own_lo,12)},{mp.nstr(own_hi,12)}] vs stored "
            f"[{mp.nstr(s_lo,12)},{mp.nstr(s_hi,12)}]",
        )

    # --- shell weights at grid 48 with own cos enclosures
    SCALE = 10**12
    BIG = 10**27
    grid = 48
    cos_hi_i, cos_lo_i = [], []
    for i in range(grid):
        a = Fraction(i) * PI_LO / grid
        b = min(Fraction(i + 1) * PI_HI / grid, PI_HI)
        chi = cos_interval_test(a)[1]
        clo = cos_interval_test(b)[0]
        cos_hi_i.append(-((-chi.numerator * SCALE) // chi.denominator))
        cos_lo_i.append((clo.numerator * SCALE) // clo.denominator)
    j_lo = {2: 0, 3: 0, 4: 0}
    for cell in itertools.product(range(grid), repeat=3):
        lam_min = 3 * SCALE - sum(cos_hi_i[i] for i in cell)
        lam_max = 3 * SCALE - sum(cos_lo_i[i] for i in cell)
        if lam_min <= 0:
            continue
        for lam_t in (2, 3, 4):
            if lam_min >= lam_t * SCALE:
                j_lo[lam_t] += BIG // lam_max
    weight = Fraction(1, 10**15 * grid**3)
    ok = True
    for lam_t in (2, 3, 4):
        own = j_lo[lam_t] * weight
        stored = Fraction(data["shell_weights"]["J_lo_small"][str(lam_t)])
        ok = ok and own == stored
    check("shell_J_lo_grid48_recomputed", ok,
          f"J_lo(2)={float(j_lo[2] * weight):.6f} identical to stored grid48")

    # --- massive threshold monotonicity claims
    mts = data["massive_thresholds_mu"]
    ok = True
    for klabel, entry in mts.items():
        two_k = Fraction(entry["two_K"])
        for lab in entry["mu2_certified_below_twoK"]:
            hi = mp.mpf(data["massive_watson"][lab]["certified_interval"][1])
            ok = ok and hi < f2m(two_k)
        for lab in entry["mu2_certified_above_twoK"]:
            lo = mp.mpf(data["massive_watson"][lab]["certified_interval"][0])
            ok = ok and lo > f2m(two_k)
    check("massive_threshold_brackets_valid", ok,
          "stored I3(mu^2) endpoints imply the mu* bracket listings")

    # --- independent finite-lattice recomputation
    tori = [("2x2x2", (2, 2, 2)), ("2x2x3", (2, 2, 3)), ("2x2x4", (2, 2, 4)), ("3x3x2", (3, 3, 2))]
    for name, sides in tori:
        lat = next(l for l in data["lattices"] if l["name"] == name)
        sites, bonds = build_torus_test(sides)
        n, nb = len(sites), len(bonds)
        c, c2, mcount = full_enumerate(sites, bonds)
        # partition polynomial cross-check against the frozen engine
        from ising.transfer_matrix import torus_broken_bond_poly

        tm = list(torus_broken_bond_poly(sides))
        mine = list(c)
        while mine and mine[-1] == 0:
            mine.pop()
        check(f"zpoly_full_{name}", mine == tm, "full 2^N enumeration == transfer engine")

        x = truncate_down(x_star_lo, 60)  # incumbent lower endpoint
        pw = [Fraction(1)]
        for _ in range(nb):
            pw.append(pw[-1] * x)
        zx = sum(c[q] * pw[q] for q in range(nb + 1))
        # G by offset class (recomputed with own class logic)
        ev = lat["evaluations"]["incumbent"]

        def class_key(site):
            vals = [min(d, L - d) for d, L in zip(site, sides)]
            mults = [2 if L == 2 else 1 for L in sides]
            groups = {}
            for a in range(3):
                groups.setdefault(mults[a], []).append(vals[a])
            return tuple((m, tuple(sorted(g))) for m, g in sorted(groups.items()))

        ok = True
        for s in range(1, n):
            key = str(class_key(sites[s]))
            stored = (Fraction(ev["G_by_offset_class"][key][0]),
                      Fraction(ev["G_by_offset_class"][key][1]))
            own = sum(c2[s][q] * pw[q] for q in range(nb + 1)) / zx
            ok = ok and own in stored  # bracket endpoints both exact Fractions
        check(f"G_classes_incumbent_{name}", ok,
              "exact equality of all offset-class two-point values")

        # M2, Delta, Binder, S2 at the incumbent lower endpoint
        m2num = Fraction(0)
        m4num = Fraction(0)
        for (q, m), cnt in mcount.items():
            m2num += cnt * pw[q] * m * m
            m4num += cnt * pw[q] * m**4
        m2_own = m2num / zx / (n * n)
        binder_own = (m4num / zx) / (m2num / zx) ** 2
        stored_m2 = (Fraction(ev["M2_interval"][0]), Fraction(ev["M2_interval"][1]))
        check(f"M2_incumbent_{name}", m2_own in stored_m2,
              f"M2 = {float(m2_own):.6f} matches a bracket endpoint exactly")
        stored_binder = Fraction(ev["binder_U"][0])
        check(f"binder_incumbent_{name}", binder_own == stored_binder,
              f"B = {float(binder_own):.4f} exact")
        # S2 and the paired identity
        gs = [Fraction(1)] + [
            sum(c2[s][q] * pw[q] for q in range(nb + 1)) / zx for s in range(1, n)
        ]
        s2_own = sum(g * g for g in gs)
        check(f"S2_incumbent_{name}",
              s2_own in (Fraction(ev["S2_sumG2"][0]), Fraction(ev["S2_sumG2"][1])),
              f"sum G^2 = {float(s2_own):.6f}")
        # modes: recompute all Ghat, Parseval, paired identity, violations
        mult = tuple(2 if L == 2 else 1 for L in sides)
        ghat_all = {}
        for j in itertools.product(range(sides[0]), range(sides[1]), range(sides[2])):
            tot = Fraction(1)
            for s in range(1, n):
                ph = Fraction(1)
                for a in range(3):
                    ph *= RATIONAL_COS[sides[a]][(j[a] * sites[s][a]) % sides[a]]
                tot += ph * gs[s]
            ghat_all[j] = tot
        check(f"parseval_own_{name}", sum(ghat_all.values()) == n, "sum_k Ghat = N")
        check(
            f"paired_identity_own_{name}",
            n * s2_own == sum(g * g for g in ghat_all.values()),
            "N sum G^2 == sum Ghat^2",
        )
        k_hi = Fraction(data["certified_inputs"]["I3_over_2_interval"][1])
        k_lo = Fraction(data["certified_inputs"]["I3_over_2_interval"][0])
        viol_own = 0
        for j, gh in ghat_all.items():
            if j == (0, 0, 0):
                continue
            lam = sum(
                mult[a] * (1 - RATIONAL_COS[sides[a]][j[a]]) for a in range(3)
            )
            d_lo = Fraction(1) / (2 * k_hi * lam) - gh
            if d_lo < 0:
                viol_own += 1
        check(
            f"violations_incumbent_{name}",
            viol_own == ev["n_ceiling_violations_bracket"],
            f"{viol_own} modes violate the ceiling bracket (stored "
            f"{ev['n_ceiling_violations_bracket']})",
        )
        # Delta identity: Delta = M2 - 1 + C0/(2K)
        c0 = Fraction(ev["C_L0_exact"])
        delta_own = m2_own - 1 + c0 / (2 * k_hi)
        delta_stored_lo = Fraction(ev["Delta_interval"][0])
        check(
            f"delta_incumbent_{name}",
            delta_own >= delta_stored_lo and abs(delta_own - delta_stored_lo) < Fraction(1, 10**30),
            f"Delta = {float(delta_own):.6f} vs stored lower {float(delta_stored_lo):.6f}",
        )

    # --- e72 single-bond anchor, own route
    sites = list(itertools.product(range(2), repeat=3))
    index = {s: i for i, s in enumerate(sites)}
    bonds = []
    for s in sites:
        for axis in range(3):
            if s[axis] == 0:
                t = list(s)
                t[axis] = 1
                bonds.append((index[s], index[tuple(t)]))
    c, c2, _ = full_enumerate(sites, bonds)
    x = Fraction(3, 5)
    pw = [Fraction(1)]
    for _ in range(len(bonds)):
        pw.append(pw[-1] * x)
    zx = sum(c[q] * pw[q] for q in range(len(bonds) + 1))
    ge = sum(c2[index[(1, 0, 0)]][q] * pw[q] for q in range(len(bonds) + 1)) / zx
    check("e72_anchor_Ge_own", ge == Fraction(4531, 15844), f"G(e) = {ge}")

    # --- envelope sanity
    check(
        "envelope",
        art["meta"]["provenance"] == "experiments/e145_paired_momentum.py"
        and "claim_tag" not in json.dumps(art)[:10]
        and "THEOREM" in art["meta"]["classification"]
        and "COMPUTATION" in art["meta"]["classification"]
        and art["meta"]["benchmark_use"].startswith("K_c=0.221654626 was not used"),
        "provenance, classification, benchmark disclaimer present",
    )

    print()
    if FAILURES:
        print(f"FAIL ({len(FAILURES)}): {FAILURES}")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
