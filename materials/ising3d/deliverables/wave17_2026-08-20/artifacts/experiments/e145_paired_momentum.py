"""Theorem-gap analysis for the rigorous upper endpoint (paired-momentum defect).

Wave-14 proved the plain Peierls head+tail certificate cannot beat the
incumbent upper endpoint K_c <= I3/2 (proofs/upper_beyond.md), and wave-8
proved the audited two-point constraint class is exactly saturated at I3/2
(proofs/upper_infrared.md, Theorem D).  This experiment delivers the gap
theorem for everything beyond those two walls:

* [THEOREM PM]  On any finite torus the exact Parseval identity gives
      Delta_L(K) = M_L^2(K) - 1 + C_L(0)/(2K)
                 = (1/N) sum_{k!=0} [ 1/(2K*lambda_T(k)) - Ghat_L(k) ],
  an exact identity between the magnetisation and the *paired-momentum
  ceiling defect* (pure algebra; the FSS infrared ceiling is only needed
  to assert Delta_L >= 0 on even tori).  Consequently, for every
  K* < I3/2 the condition
      liminf_{L even} Delta_L(K*) > delta*(K*) := I3/(2K*) - 1
  is exactly the weakest additional inequality that certifies K_c <= K*:
  it is equivalent (up to the certified O(1/L) correction of
  proofs/mag_floor.md) to a uniform magnetisation floor, while no
  two-point consequence can supply it (Theorem D).

* Exact thresholds: delta*(K*) across the certified uncertainty window,
  the massive-ceiling corollary K_c <= I3(mu^2)/2 with I3(mu^2) certified
  by a directed interval series over exact walk return probabilities, and
  the shell corollary eps*(K*, Lambda) = (I3-2K*)/J_lo(Lambda) with
  certified Brillouin-zone shell weights.

* Exact finite-lattice instance: every torus the in-house engine reaches
  (2x2x2, 2x2x3, 2x2x4, 3x3x2, CODE_API doubled-bond convention) is
  enumerated exactly (integer signed histograms in x = e^{-2K}); offset
  classes, all momentum modes, M_L^2, Delta_L, the summed two-point second
  moment, the Binder ratio and the paired-momentum four-point intensities
  <|sighat(k)|^4> are evaluated at certified rational enclosures of the
  incumbent bracket and of benchmark-free rational targets, with every
  identity re-verified exactly.

Run from the repository root:
    .venv/bin/python experiments/e145_paired_momentum.py
"""

from __future__ import annotations

from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import itertools
import sys
import json
import math
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)  # certified x* endpoints exceed 4300 digits
SCRIPT = "experiments/e145_paired_momentum.py"
RESULT = ROOT / "results" / "bounds" / "paired_momentum.json"

UPPER_INFRARED = ROOT / "results" / "bounds" / "upper_infrared.json"
UPPER_BEYOND = ROOT / "results" / "bounds" / "upper_beyond.json"
SAW_UNION4 = ROOT / "results" / "bounds" / "saw_union4.json"

DPS = 50
EXP_ORDER = 44
COS_ORDER = 20
PI_LO = Fraction(333, 106)
PI_HI = Fraction(355, 113)
WALK_ANCHOR_STEPS = 30
CONV_NMAX = 3200
SHELL_GRID = 96
SHELL_GRID_SMALL = 48

RATIONAL_COS = {
    2: (Fraction(1), Fraction(-1)),
    3: (Fraction(1), Fraction(-1, 2), Fraction(-1, 2)),
    4: (Fraction(1), Fraction(0), Fraction(-1), Fraction(0)),
}
RATIONAL_SIN = {
    2: (Fraction(0), Fraction(0)),
    3: (Fraction(0), None, None),
    4: (Fraction(0), Fraction(1), Fraction(0), Fraction(-1)),
}


# ---------------------------------------------------------------------------
# exact walk return counts N_{2n}(0) = sum_{i+j+k=n} (2n)!/(i!^2 j!^2 k!^2)
# ---------------------------------------------------------------------------

def return_counts(n_max: int) -> list[int]:
    out = []
    fact = [math.factorial(m) for m in range(2 * n_max + 1)]
    for n in range(n_max + 1):
        total = 0
        row = fact[2 * n]
        for i in range(n + 1):
            fi = fact[i]
            for j in range(n - i + 1):
                total += row // (fi * fi * fact[j] * fact[j] * fact[n - i - j] ** 2)
        out.append(total)
    return out


# ---------------------------------------------------------------------------
# exact transcendental enclosures (producer route: alternating series)
# ---------------------------------------------------------------------------

def exp_neg_interval(t: Fraction) -> tuple[Fraction, Fraction]:
    """Certified enclosure of exp(-t), t > 0.

    t = N*u, u <= 1/2: the alternating Taylor series of e^{-u} has strictly
    decreasing term magnitudes, so every odd partial sum is a strict lower
    bound and every even partial sum a strict upper bound; N-th powers are
    monotone on positives.
    """
    assert t > 0
    n_split = max(1, math.ceil(2 * float(t)))
    u = t / n_split
    assert 0 < u <= Fraction(1, 2)
    s = Fraction(1)
    term = Fraction(1)
    low = high = None
    for k in range(1, EXP_ORDER + 2):
        term = term * u / k
        s = s + (term if k % 2 == 0 else -term)
        if k % 2 == 1:
            low = s              # partial sum ending in a negative term: lower
        else:
            high = s             # partial sum ending in a positive term: upper
    return (low**n_split, high**n_split)


def _cos_taylor(y: Fraction) -> tuple[Fraction, Fraction]:
    """cos(y) for 0 <= y < 1 (alternating, strictly decreasing terms).

    term_m = y^(2m)/(2m)!: partial sums ending in a negative term are
    strict lower bounds, those ending in a positive term strict upper
    bounds; the drop to the next partial sum is bounded by the next term.
    """
    assert 0 <= y < 1
    if y == 0:
        return (Fraction(1), Fraction(1))
    term = Fraction(1)
    s = Fraction(1)
    low = high = None
    y2 = y * y
    for m in range(1, COS_ORDER + 2):
        term = term * y2 / ((2 * m - 1) * (2 * m))
        if m % 2 == 1:
            s = s - term
            low = s
            high = s + term * y2 / ((2 * m + 1) * (2 * m + 2))
        else:
            s = s + term
            high = s
            low = s - term * y2 / ((2 * m + 1) * (2 * m + 2))
    return (low, high)


def cos_interval(x: Fraction) -> tuple[Fraction, Fraction]:
    """Certified enclosure of cos(x) for 0 <= x <= pi (exact Fractions).

    Upper half via cos(x) = -cos(pi-x); then cos(x) = 8c^4-8c^2+1 with
    c = cos(x/4) in [0, pi/8] where f(c) is increasing (c >= cos(pi/8) >
    1/sqrt(2) whenever x > 0; at x = 0 return exact 1).
    """
    assert 0 <= x <= PI_HI
    if x == 0:
        return (Fraction(1), Fraction(1))
    if x > PI_LO / 2:
        # cos(x) = -cos(pi - x), pi - x in (0, pi/2]
        a = max(PI_LO - x, Fraction(0))
        b = max(PI_HI - x, Fraction(0))
        clo, chi = _cos_quad(b)[0], _cos_quad(a)[1]  # cos decreasing
        return (-chi, -clo)
    return _cos_quad(x)


def _cos_quad(x: Fraction) -> tuple[Fraction, Fraction]:
    """cos(x) for 0 <= x <= pi/2 via c = cos(x/4) in [0, pi/8]."""
    if x == 0:
        return (Fraction(1), Fraction(1))
    y = x / 4
    assert 0 < y <= PI_HI / 8 < 1
    c_lo = _cos_taylor(y)[0]
    c_hi = _cos_taylor(y)[1]
    f_lo = 8 * c_lo**4 - 8 * c_lo**2 + 1
    f_hi = 8 * c_hi**4 - 8 * c_hi**2 + 1
    return (f_lo, f_hi)


# ---------------------------------------------------------------------------
# certified inputs
# ---------------------------------------------------------------------------

def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_certified() -> dict:
    src = json.loads(UPPER_INFRARED.read_text())
    i3 = src["data"]["certified_constants"]["I3"]
    i3h = src["data"]["certified_constants"]["I3_over_2"]
    beyond = json.loads(UPPER_BEYOND.read_text())
    xstar = beyond["data"]["limitation_certificate"]["x_star_interval"]
    saw = json.loads(SAW_UNION4.read_text())
    floor40 = saw["data"]["endpoint"]["floor_40"]
    return {
        "I3": (Fraction(i3[0]), Fraction(i3[1])),
        "I3_over_2": (Fraction(i3h[0]), Fraction(i3h[1])),
        "x_star": (
            Fraction(xstar["lower_exact"]),
            Fraction(xstar["upper_exact"]),
        ),
        "kc_lower_floor_40": Fraction(floor40),
        "sha256": {
            "results/bounds/upper_infrared.json": _sha256(UPPER_INFRARED),
            "results/bounds/upper_beyond.json": _sha256(UPPER_BEYOND),
            "results/bounds/saw_union4.json": _sha256(SAW_UNION4),
        },
    }


# ---------------------------------------------------------------------------
# massive Watson integral I3(mu^2) = <1/(lambda+mu^2)> by directed intervals
# I3(mu^2) = sum_n 3^n R_n / (3+mu^2)^{n+1},  R_n = N_{2n}/6^{2n},
# N_{2n} = (2n)! * [z^n] I0(2 sqrt z)^3 (all-positive series: no cancellation)
# ---------------------------------------------------------------------------

def interval_series_s(n_max: int, dps: int):
    mp.iv.dps = dps
    one = mp.iv.mpf(1)
    s = [one]
    cur = one
    for a in range(1, n_max + 1):
        cur = cur * (one / mp.iv.mpf(a * a))
        s.append(cur)
    return s


def convolve_positive(a, b, n_max):
    out = [mp.iv.mpf(0)] * (n_max + 1)
    for i, ai in enumerate(a):
        if i > n_max:
            break
        for j, bj in enumerate(b[: n_max - i + 1]):
            out[i + j] = out[i + j] + ai * bj
    return out


def massive_watson(mu2: Fraction, c_series, dps: int, n_use: int):
    mp.iv.dps = dps
    base = 3 + mu2
    r = Fraction(3) / base
    total = mp.iv.mpf(0)
    factor = mp.iv.mpf(1) / mp.iv.mpf(int(base))
    rfac = mp.iv.mpf(1)
    rn_iv = mp.iv.mpf(1)
    for n in range(n_use + 1):
        if n > 0:
            rn_iv = (
                mp.iv.mpf(int(math.factorial(2 * n)))
                * c_series[n]
                / mp.iv.mpf(int(6 ** (2 * n)))
            )
            rfac = rfac * mp.iv.mpf(int(r.numerator)) / mp.iv.mpf(int(r.denominator))
        # n = 0: R_0 = 1, rfac = 1, factor = 1/(3+mu^2)
        total = total + factor * rfac * (rn_iv if n > 0 else mp.iv.mpf(1))
    tail = (r ** (n_use + 1)) / mu2  # <= sum_{n>N} r^n/(3+mu^2), since R_n<=1
    hi = total.b + mp.iv.mpf(tail.numerator) / mp.iv.mpf(tail.denominator)

    def _end_str(v, side: str) -> str:
        parts = mp.nstr(v, 30).strip("[]").split(", ")
        return parts[0] if side == "a" else parts[-1]

    return (_end_str(total, "a"), _end_str(hi, "b"), frac_str(tail))


def frac_str(fr: Fraction) -> str:
    return f"{fr.numerator}/{fr.denominator}"


# ---------------------------------------------------------------------------
# Brillouin-zone shell weights, exact rational cells
# ---------------------------------------------------------------------------

def shell_weights(grid: int):
    """Certified BZ shell weights by exact integer cell sums.

    Per interval [i*pi/G, (i+1)*pi/G]: cos_hi = ceiling of the enclosure
    sup to a multiple of 1e-12 (stays >= true sup since cos is decreasing);
    cos_lo = floor (stays <= true inf).  Scaled lambda bounds are integers
    (SCALE = 10^12).  A cell with scaled lambda_min >= Lambda*SCALE lies
    certainly in {lambda >= Lambda} and contributes floor(10^27 /
    lambda_max_scaled)/10^15, itself a lower bound of 1/lambda at every
    point of the cell.
    """
    SCALE = 10**12
    BIG = 10**27
    cos_lo_i, cos_hi_i = [], []
    for i in range(grid):
        a = Fraction(i) * PI_LO / grid          # <= true left end
        b = min(Fraction(i + 1) * PI_HI / grid, PI_HI)  # >= true right end
        chi = cos_interval(a)[1]
        clo = cos_interval(b)[0]
        cos_hi_i.append(-((-chi.numerator * SCALE) // chi.denominator))  # ceil
        cos_lo_i.append((clo.numerator * SCALE) // clo.denominator)  # floor
    j_lo = {2: 0, 3: 0, 4: 0}
    inner_lo = {2: 0, 3: 0, 4: 0}
    sum_all = 0
    for cell in itertools.product(range(grid), repeat=3):
        lam_min = 3 * SCALE - sum(cos_hi_i[i] for i in cell)
        lam_max = 3 * SCALE - sum(cos_lo_i[i] for i in cell)
        if lam_min <= 0:
            continue
        sum_all += BIG // lam_max
        for lam_t in (2, 3, 4):
            if lam_min >= lam_t * SCALE:
                j_lo[lam_t] += BIG // lam_max
            if lam_max <= lam_t * SCALE:
                inner_lo[lam_t] += BIG // lam_min
    weight = Fraction(1, 10**15 * grid**3)
    return (
        {k: v * weight for k, v in j_lo.items()},
        {k: v * weight for k, v in inner_lo.items()},
        sum_all * weight,
    )


# ---------------------------------------------------------------------------
# exact torus engine (integer signed histograms in x = e^{-2K})
# ---------------------------------------------------------------------------

def build_torus(sides):
    sites = list(itertools.product(range(sides[0]), range(sides[1]), range(sides[2])))
    index = {s: i for i, s in enumerate(sites)}
    bonds = []
    for s in sites:
        for axis in range(3):
            if sides[axis] == 2 and s[axis] == 1:
                continue  # length-2 pair emitted from coordinate 0 only
            t = list(s)
            t[axis] = (s[axis] + 1) % sides[axis]
            mult = 2 if sides[axis] == 2 else 1
            for _ in range(mult):
                bonds.append((index[s], index[tuple(t)]))
    return sites, bonds


def offset_class_key(offset, sides):
    """Canonical orbit key: reflect each axis (d -> min(d, L-d)) and permute
    only axes with equal bond multiplicity (length-2 axes are doubled)."""
    vals = [min(d, L - d) for d, L in zip(offset, sides)]
    mults = [2 if L == 2 else 1 for L in sides]
    groups = {}
    for a in range(3):
        groups.setdefault(mults[a], []).append(vals[a])
    return tuple((m, tuple(sorted(g))) for m, g in sorted(groups.items()))


def enumerate_lattice(sites, bonds, integer_modes):
    """sigma_0 = +1 fixed (all observables even). Bit i set <=> spin i = -1."""
    n = len(sites)
    nb = len(bonds)
    bond_pairs = [(min(a, b), max(a, b)) for a, b in bonds]
    c = [0] * (nb + 1)
    c2 = [[0] * (nb + 1) for _ in range(n)]
    m_span = 2 * n + 1
    jm = [0] * ((nb + 1) * m_span)
    i_span = 4 * n * n + 1
    mode_hists = []
    for cvec, svec in integer_modes:
        pos = [i for i in range(n) if cvec[i] > 0]
        neg = [i for i in range(n) if cvec[i] < 0]
        spos = [i for i in range(n) if svec[i] > 0]
        sneg = [i for i in range(n) if svec[i] < 0]
        masks = (
            sum(1 << i for i in pos),
            sum(1 << i for i in neg),
            sum(1 << i for i in spos),
            sum(1 << i for i in sneg),
        )
        mode_hists.append(
            {
                "pos_neg": (len(pos), len(neg), len(spos), len(sneg)),
                "masks": masks,
                "jk": [0] * ((nb + 1) * i_span),
            }
        )
    i_span_used = i_span
    for mask in range(0, 1 << n, 2):  # even masks: bit 0 (= site 0) is +1
        q = 0
        for a, b in bond_pairs:
            q += (mask >> a ^ mask >> b) & 1
        c[q] += 1
        pc = mask.bit_count()
        m = n - 2 * pc
        jm[q * m_span + (m + n)] += 1
        for s in range(1, n):
            c2[s][q] += -1 if (mask >> s) & 1 else 1
        for mh in mode_hists:
            (np_, nn_, nsp, nsn), (mp_, mn_, msp, msn) = mh["pos_neg"], mh["masks"]
            a = (np_ - 2 * (mask & mp_).bit_count()) - (
                nn_ - 2 * (mask & mn_).bit_count()
            )
            b = (nsp - 2 * (mask & msp).bit_count()) - (
                nsn - 2 * (mask & msn).bit_count()
            )
            mh["jk"][q * i_span_used + a * a + b * b] += 1
    return {
        "n": n,
        "nb": nb,
        "c": c,
        "c2": c2,
        "jm": jm,
        "m_span": m_span,
        "mode_hists": mode_hists,
        "i_span": i_span_used,
    }


def _truncate_down(f: Fraction, digits: int) -> Fraction:
    """floor(f * 10^digits)/10^digits <= f (positive f)."""
    scale = 10**digits
    return Fraction(f.numerator * scale // f.denominator, scale)


def _truncate_up(f: Fraction, digits: int) -> Fraction:
    """ceil(f * 10^digits)/10^digits >= f (positive f)."""
    scale = 10**digits
    return Fraction(-((-f.numerator * scale) // f.denominator), scale)


def powers_table(x: Fraction, nb: int) -> list[Fraction]:
    out = [Fraction(1)]
    for _ in range(nb):
        out.append(out[-1] * x)
    return out


def eval_hist(hist, pw) -> Fraction:
    return sum(coef * pw[q] for q, coef in enumerate(hist) if coef)


def cos_phase(j, site, sides) -> Fraction:
    ph = Fraction(1)
    for a in range(3):
        ph *= RATIONAL_COS[sides[a]][(j[a] * site[a]) % sides[a]]
    return ph


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> int:
    mp.mp.dps = DPS  # NB: the `mp.dps = n` shorthand is inert in this venv
    # (mpmath 1.3.0 under CPython 3.14); the context must be set directly.
    checks = []

    def record(name, passed, detail):
        checks.append({"name": name, "pass": bool(passed), "detail": detail})
        print(("PASS " if passed else "FAIL ") + name + " -- " + detail)
        if not passed:
            raise AssertionError(name)

    cert = load_certified()
    i3_lo, i3_hi = cert["I3"]

    # own enclosure of e^{-I3} overlaps the stored x* bracket
    xa = exp_neg_interval(i3_lo)
    xb = exp_neg_interval(i3_hi)
    xs_lo, xs_hi = cert["x_star"]
    record(
        "exp_enclosure_overlaps_stored_xstar",
        xb[1] >= xs_lo and xa[0] <= xs_hi,
        f"own e^-I3 in [{float(xa[0]):.12f},{float(xb[1]):.12f}] vs stored "
        f"[{float(xs_lo):.12f},{float(xs_hi):.12f}]",
    )

    counts = return_counts(WALK_ANCHOR_STEPS)
    partial = Fraction(0)
    for n_, Nn in enumerate(counts):
        partial += Fraction(Nn, 6 ** (2 * n_))
    record(
        "I3_walk_partial_below_certified_lower_endpoint",
        partial / 3 <= i3_hi,
        f"(1/3)sum R_n (n<=2*{WALK_ANCHOR_STEPS}) = {float(partial / 3):.15f} "
        f"<= I3_hi {float(i3_hi):.15f} (positive terms: strict lower bound)",
    )

    # delta*(K*) thresholds
    targets = {
        "incumbent_I3_2_lo": cert["I3_over_2"][0],
        "incumbent_I3_2_hi": cert["I3_over_2"][1],
        "K=63/250": Fraction(63, 250),
        "K=1/4": Fraction(1, 4),
        "K=49/200": Fraction(49, 200),
        "K=6/25": Fraction(6, 25),
        "K=11/50": Fraction(11, 50),
        "K=kc_lower_floor_40": cert["kc_lower_floor_40"],
    }
    delta_star = {}
    for label, k in targets.items():
        delta_star[label] = (
            frac_str(i3_lo / (2 * k) - 1),
            frac_str(i3_hi / (2 * k) - 1),
        )

    # massive Watson integrals
    s_series = interval_series_s(CONV_NMAX, DPS)
    t_series = convolve_positive(s_series, s_series, CONV_NMAX)
    c_series = convolve_positive(t_series, s_series, CONV_NMAX)
    mu_grid = [
        ("1/2", Fraction(1, 2), 400),
        ("1/4", Fraction(1, 4), 800),
        ("1/8", Fraction(1, 8), 1600),
        ("1/16", Fraction(1, 16), 2400),
        ("1/32", Fraction(1, 32), CONV_NMAX),
        ("1/128", Fraction(1, 128), CONV_NMAX),
        ("1/1024", Fraction(1, 1024), CONV_NMAX),
    ]
    massive = {}
    for label, mu2, n_use in mu_grid:
        lo, hi, tail = massive_watson(mu2, c_series, DPS, n_use)
        massive[label] = {
            "mu2": frac_str(mu2),
            "n_used": n_use,
            "certified_interval": [lo, hi],
            "tail_bound": tail,
            "tight": n_use < CONV_NMAX or label in ("1/32",),
        }
    # anchor: interval R_{2n} vs exact multinomial for n<=6
    mp.iv.dps = DPS
    ok_all = True
    for n_ in range(1, 7):
        rn_int = Fraction(counts[n_], 6 ** (2 * n_))
        rn_iv = (
            mp.iv.mpf(int(math.factorial(2 * n_)))
            * c_series[n_]
            / mp.iv.mpf(int(6 ** (2 * n_)))
        )
        val = mp.iv.mpf(rn_int.numerator) / mp.iv.mpf(rn_int.denominator)
        ok_all = ok_all and (rn_iv.a <= val <= rn_iv.b)
    record("interval_Rn_matches_exact", ok_all, "n=1..6 exact multinomial vs interval convolution")

    # massive thresholds: smallest grid mu with certified I3(mu^2) <= 2K*
    mu_thresholds = {}
    for klabel in ("K=63/250", "K=1/4", "K=49/200", "K=6/25"):
        two_k = mp.mpf((2 * targets[klabel]).numerator) / mp.mpf((2 * targets[klabel]).denominator)
        below = [
            lab
            for lab, d in massive.items()
            if mp.mpf(d["certified_interval"][1]) < two_k
        ]
        above = [
            lab
            for lab, d in massive.items()
            if mp.mpf(d["certified_interval"][0]) > two_k
        ]
        mu_thresholds[klabel] = {
            "two_K": frac_str(2 * targets[klabel]),
            "mu2_certified_below_twoK": below,
            "mu2_certified_above_twoK": above,
        }

    # shell weights
    j_lo_big, inner_lo_big, i3_lower_grid = shell_weights(SHELL_GRID)
    j_lo_small, _, i3_lower_small = shell_weights(SHELL_GRID_SMALL)
    record(
        "shell_grid_lower_bound_below_I3",
        i3_lower_grid <= i3_hi and i3_lower_grid <= i3_lower_small + Fraction(1),
        f"grid96 lower sum {float(i3_lower_grid):.6f} <= I3_hi {float(i3_hi):.6f}; "
        f"grid48 {float(i3_lower_small):.6f}",
    )
    shell_eps = {}
    for lam_t in (2, 3, 4):
        jl = j_lo_big[lam_t]
        for klabel in ("K=63/250", "K=1/4", "K=49/200", "K=6/25"):
            k = targets[klabel]
            eps_suf = (i3_lo - 2 * k) / jl  # sufficient if eps > this
            eps_req = (i3_hi - 2 * k) / jl  # provably sufficient only above this
            shell_eps[f"{klabel},Lambda={lam_t}"] = {
                "eps_star_interval": [float(eps_suf), float(eps_req)],
                "eps_star_interval_exact": [frac_str(eps_suf), frac_str(eps_req)],
                "J_lo_exact": frac_str(jl),
            }
    for lam_t in (2, 3):
        a, b = shell_eps[f"K=1/4,Lambda={lam_t}"]["eps_star_interval"]
        record(
            f"shell_eps_usable_Lambda{lam_t}",
            0 < a and b < 1,
            f"eps*(K*=1/4,Lambda={lam_t}) in [{a:.5f},{b:.5f}] within (0,1)",
        )

    # ------------------------------------------------------------------
    # finite lattices
    # ------------------------------------------------------------------
    from ising.transfer_matrix import torus_broken_bond_poly

    tori = [("2x2x2", (2, 2, 2)), ("2x2x3", (2, 2, 3)), ("2x2x4", (2, 2, 4)), ("3x3x2", (3, 3, 2))]
    k_points = {
        "incumbent": (cert["I3_over_2"][0], cert["I3_over_2"][1]),
        "K=1/4": (Fraction(1, 4), Fraction(1, 4)),
        "K=49/200": (Fraction(49, 200), Fraction(49, 200)),
        "K=6/25": (Fraction(6, 25), Fraction(6, 25)),
        "K=kc_lower_floor_40": (cert["kc_lower_floor_40"], cert["kc_lower_floor_40"]),
    }
    lattices_out = []
    for name, sides in tori:
        sites, bonds = build_torus(sides)
        index = {s: i for i, s in enumerate(sites)}
        n, nb = len(sites), len(bonds)
        mult = tuple(2 if L == 2 else 1 for L in sides)
        modes = []
        for j in itertools.product(range(sides[0]), range(sides[1]), range(sides[2])):
            cosv = tuple(RATIONAL_COS[sides[a]][j[a]] for a in range(3))
            sinv = tuple(RATIONAL_SIN[sides[a]][j[a]] for a in range(3))
            lam = sum(mult[a] * (1 - cosv[a]) for a in range(3))
            modes.append({"j": j, "lambda": lam, "integer": all(v is not None for v in sinv)})
        integer_modes = []
        for md in modes:
            if md["integer"] and md["j"] != (0, 0, 0):
                # per-site coefficients of Re/Im sigma_hat(k):
                # e^{-i k.x} = prod_a (cos th_a - i sin th_a) with
                # th_a = k_a x_a; all factors in {0,+-1} for these tori
                cvec = []
                svec = []
                for site in sites:
                    cc = [RATIONAL_COS[sides[a]][(md["j"][a] * site[a]) % sides[a]] for a in range(3)]
                    ss = [RATIONAL_SIN[sides[a]][(md["j"][a] * site[a]) % sides[a]] for a in range(3)]
                    re = cc[0] * cc[1] * cc[2] - cc[0] * ss[1] * ss[2] - ss[0] * cc[1] * ss[2] - ss[0] * ss[1] * cc[2]
                    im = -(ss[0] * cc[1] * cc[2] + cc[0] * ss[1] * cc[2] + cc[0] * cc[1] * ss[2] - ss[0] * ss[1] * ss[2])
                    cvec.append(int(re))
                    svec.append(int(im))
                integer_modes.append((cvec, svec))
        eng = enumerate_lattice(sites, bonds, integer_modes)
        tm_poly = list(torus_broken_bond_poly(tuple(sides)))
        mine_poly = [2 * a for a in eng["c"]]
        while mine_poly and mine_poly[-1] == 0:
            mine_poly.pop()
        record(
            f"zpoly_matches_transfer_{name}",
            mine_poly == tm_poly,
            f"2 x (sigma_0=+1 half histograms, trailing zeros stripped) == "
            f"torus_broken_bond_poly({sides}) ({nb} bonds)",
        )
        classes = {}
        for s in range(1, n):
            classes.setdefault(offset_class_key(sites[s], sides), []).append(s)
        for key, members in classes.items():
            base = eng["c2"][members[0]]
            same = all(eng["c2"][s] == base for s in members[1:])
            record(
                f"orbit_constant_{name}_{key}",
                same,
                f"{len(members)} members share one signed histogram",
            )
        class_rep = {key: members[0] for key, members in classes.items()}

        evaluations = {}
        for klabel, (k_lo, k_hi) in k_points.items():
            if klabel == "incumbent":
                x_lo, x_hi = cert["x_star"]
                x_lo = _truncate_down(x_lo, 60)
                x_hi = _truncate_up(x_hi, 60)
            else:
                ea = exp_neg_interval(2 * k_lo)
                eb = exp_neg_interval(2 * k_hi)
                x_lo, x_hi = min(ea[0], eb[0]), max(ea[1], eb[1])
            pw_lo = powers_table(x_lo, nb)
            pw_hi = powers_table(x_hi, nb)
            zx_lo = eval_hist(eng["c"], pw_lo)
            zx_hi = eval_hist(eng["c"], pw_hi)
            g_lo_vals = {
                s: eval_hist(eng["c2"][s], pw_lo) / zx_lo for s in range(1, n)
            }
            g_hi_vals = {
                s: eval_hist(eng["c2"][s], pw_hi) / zx_hi for s in range(1, n)
            }
            gvals = {}
            for key, s0 in class_rep.items():
                a, b = g_lo_vals[s0], g_hi_vals[s0]
                gvals[str(key)] = [frac_str(min(a, b)), frac_str(max(a, b))]

            def ghat(j, x_end: int) -> Fraction:
                total = Fraction(1)
                gv = g_lo_vals if x_end == 0 else g_hi_vals
                for s in range(1, n):
                    total += cos_phase(j, sites[s], sides) * gv[s]
                return total

            mode_records = []
            n_viol = 0
            n_certain_pos = 0
            for md in modes:
                jj = md["j"]
                gh0, gh1 = ghat(jj, 0), ghat(jj, 1)
                entry = {
                    "j": list(jj),
                    "lambda_T": frac_str(md["lambda"]),
                    "ghat": [frac_str(min(gh0, gh1)), frac_str(max(gh0, gh1))],
                }
                if jj != (0, 0, 0):
                    lam = md["lambda"]
                    d_lo = Fraction(1) / (2 * k_hi * lam) - max(gh0, gh1)
                    d_hi = Fraction(1) / (2 * k_lo * lam) - min(gh0, gh1)
                    entry["defect"] = [frac_str(d_lo), frac_str(d_hi)]
                    if d_lo < 0:
                        n_viol += 1
                    if d_lo > 0:
                        n_certain_pos += 1
                    if min(gh0, gh1) < 0:
                        record(f"ghat_nonneg_{name}_{klabel}_{jj}", False, "GKS")
                mode_records.append(entry)

            # M2, S2, binder, chi
            m2_vals, s2_vals, binder_vals = [], [], []
            for pw, zx in ((pw_lo, zx_lo), (pw_hi, zx_hi)):
                num2 = Fraction(0)
                num4 = Fraction(0)
                for idx, v in enumerate(eng["jm"]):
                    if v:
                        q, moff = divmod(idx, eng["m_span"])
                        m = moff - n
                        w = v * pw[q]
                        num2 += w * m * m
                        num4 += w * m**4
                s2 = Fraction(1)
                for s in range(1, n):
                    g = eval_hist(eng["c2"][s], pw) / zx
                    s2 += g * g
                m2_vals.append(num2 / (zx * n * n))
                s2_vals.append(s2)
                # Binder ratio B = <M^4>/<M^2>^2 (Gaussian ref 3; Lebowitz B<=3)
                binder_vals.append((num4 / zx) / (num2 / zx) ** 2)
            c0 = Fraction(0)
            for md in modes:
                if md["j"] != (0, 0, 0):
                    c0 += 1 / md["lambda"]
            c0 /= n
            delta_lo = min(m2_vals) - 1 + c0 / (2 * k_hi)
            delta_hi = max(m2_vals) - 1 + c0 / (2 * k_lo)

            # paired-mode four-point intensities
            pair_out = {}
            for mh, md in zip(eng["mode_hists"], [m for m in modes if m["integer"] and m["j"] != (0, 0, 0)]):
                jj = md["j"]
                fvals = []
                for pw, zx in ((pw_lo, zx_lo), (pw_hi, zx_hi)):
                    num = Fraction(0)
                    for idx, v in enumerate(mh["jk"]):
                        if v:
                            q, inten = divmod(idx, eng["i_span"])
                            num += v * pw[q] * inten
                    fvals.append(num / zx)
                gh_min = Fraction(min(ghat(jj, 0), ghat(jj, 1)))
                real_mode = all((2 * jj[a]) % sides[a] == 0 for a in range(3))
                gauss = (3 if real_mode else 2) * (n * gh_min) ** 2
                pair_out[str(jj)] = {
                    "F_interval": [frac_str(min(fvals)), frac_str(max(fvals))],
                    "gaussian_ref": frac_str(gauss),
                }

            evaluations[klabel] = {
                "K_interval": [frac_str(k_lo), frac_str(k_hi)],
                "x_interval": [frac_str(x_lo), frac_str(x_hi)],
                "M2_interval": [frac_str(min(m2_vals)), frac_str(max(m2_vals))],
                "chi_per_site": [frac_str(min(m2_vals) * n), frac_str(max(m2_vals) * n)],
                "C_L0_exact": frac_str(c0),
                "Delta_interval": [frac_str(delta_lo), frac_str(delta_hi)],
                "S2_sumG2": [frac_str(min(s2_vals)), frac_str(max(s2_vals))],
                "binder_U": [frac_str(min(binder_vals)), frac_str(max(binder_vals))],
                "G_by_offset_class": gvals,
                "modes": mode_records,
                "n_modes_nonzero": len(modes) - 1,
                "n_ceiling_violations_bracket": n_viol,
                "n_defect_certainly_positive": n_certain_pos,
                "pair_modes": pair_out,
            }

        # exact identities at the incumbent lower x endpoint
        x_chk = _truncate_down(cert["x_star"][0], 60)
        pw = powers_table(x_chk, nb)
        zx = eval_hist(eng["c"], pw)
        gv = {s: eval_hist(eng["c2"][s], pw) / zx for s in range(1, n)}
        tot_parseval = Fraction(1)
        for md in modes:
            gh = Fraction(1)
            for s in range(1, n):
                gh += cos_phase(md["j"], sites[s], sides) * gv[s]
            tot_parseval += gh - (Fraction(1) if md["j"] != (0, 0, 0) else 0) if False else 0
        # recompute cleanly
        tot_parseval = Fraction(0)
        for md in modes:
            gh = Fraction(1)
            for s in range(1, n):
                gh += cos_phase(md["j"], sites[s], sides) * gv[s]
            tot_parseval += gh
        record(
            f"parseval_{name}",
            tot_parseval == n,
            f"sum over all {len(modes)} modes of Ghat = {tot_parseval} = N",
        )
        lhs = Fraction(1) + sum(g * g for g in gv.values())
        rhs = Fraction(0)
        for md in modes:
            gh = Fraction(1)
            for s in range(1, n):
                gh += cos_phase(md["j"], sites[s], sides) * gv[s]
            rhs += gh * gh
        record(
            f"paired_identity_{name}",
            lhs * n == rhs,
            f"N*sum_z G(z)^2 == sum_k Ghat(k)^2 ({float(lhs):.9f})",
        )

        lattices_out.append(
            {
                "name": name,
                "sides": list(sides),
                "n_sites": n,
                "n_bonds": nb,
                "bond_multiplicities": list(mult),
                "z_poly": list(eng["c"]),
                "offset_class_sizes": {str(k): len(v) for k, v in classes.items()},
                "evaluations": evaluations,
            }
        )

    # e72 single-bond 8-site anchor at x = 3/5
    sites = list(itertools.product(range(2), repeat=3))
    index = {s: i for i, s in enumerate(sites)}
    bonds = []
    for s in sites:
        for axis in range(3):
            if s[axis] == 0:
                t = list(s)
                t[axis] = 1
                bonds.append((index[s], index[tuple(t)]))
    eng = enumerate_lattice(sites, bonds, [])
    x_anchor = Fraction(3, 5)
    pw = powers_table(x_anchor, len(bonds))
    zx = eval_hist(eng["c"], pw)
    g_e = eval_hist(eng["c2"][index[(1, 0, 0)]], pw) / zx
    record("e72_anchor_Ge", g_e == Fraction(4531, 15844), f"G(e)={g_e}")
    gh = Fraction(1)
    gv = {s: eval_hist(eng["c2"][s], pw) / zx for s in range(1, 8)}
    for s in range(1, 8):
        gh += cos_phase((1, 0, 0), sites[s], (2, 2, 2)) * gv[s]
    record("e72_anchor_mode", gh == Fraction(280125, 269348), f"Ghat(pi,0,0)={gh}")

    result = {
        "meta": {
            "script": SCRIPT,
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "classification": (
                "[THEOREM] Theorem PM: the paired-momentum ceiling-defect floor "
                "liminf_{L even} Delta_L(K*) > delta*(K*) = I3/(2K*)-1 is the exact "
                "weakest additional inequality class improving the certified upper "
                "endpoint (equivalent up to the certified mag-floor O(1/L) correction "
                "to a uniform magnetisation floor; unreachable by any two-point "
                "consequence by upper_infrared Theorem D). "
                "[COMPUTATION] exact finite-torus instance data (2x2x2/2x2x3/2x2x4/3x3x2, "
                "integer signed histograms in x=e^{-2K} evaluated at certified rational "
                "enclosures), massive-Watson thresholds I3(mu^2) by directed positive "
                "interval series over walk return probabilities with exact rational tail "
                "r^(N+1)/mu^2, and Brillouin-zone shell thresholds eps*(K*,Lambda)="
                "(I3-2K*)/J_lo(Lambda) with exact rational cell bounds."
            ),
            "arithmetic": (
                "finite lattices and thresholds: exact integers/Fractions; exponentials: "
                "exact-Fraction alternating-Taylor enclosures (subdivision u<=1/2, order 44); "
                "I3(mu^2): mpmath.iv directed intervals at dps=50, all-positive series "
                "(no cancellation), tail an exact Fraction upper bound using R_n<=1; "
                "no floating-point value decides any stored comparison"
            ),
            "benchmark_use": (
                "K_c=0.221654626 was not used to select any threshold, target, or "
                "rounding; evaluation points are the certified endpoints (I3/2 bracket, "
                "kc lower floor) and the plain rationals 63/250, 1/4, 49/200, 6/25, 11/50"
            ),
            "provenance": SCRIPT,
            "inputs_sha256": cert["sha256"],
        },
        "data": {
            "certified_inputs": {
                "I3_interval": [str(i3_lo), str(i3_hi)],
                "I3_over_2_interval": [str(v) for v in cert["I3_over_2"]],
                "x_star_interval": [str(v) for v in cert["x_star"]],
                "kc_lower_floor_40": str(cert["kc_lower_floor_40"]),
            },
            "delta_star_thresholds": delta_star,
            "massive_watson": massive,
            "massive_thresholds_mu": mu_thresholds,
            "shell_weights": {
                "grid": SHELL_GRID,
                "J_lo": {str(k): frac_str(v) for k, v in j_lo_big.items()},
                "inner_upper": {str(k): frac_str(v) for k, v in inner_lo_big.items()},
                "I3_grid_lower_sum": frac_str(i3_lower_grid),
                "grid_small": SHELL_GRID_SMALL,
                "J_lo_small": {str(k): frac_str(v) for k, v in j_lo_small.items()},
            },
            "shell_thresholds_eps_star": shell_eps,
            "lattices": lattices_out,
            "walk_return_counts_N2n": counts,
        },
        "checks": checks,
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(result, indent=1))
    print(f"\nwrote {RESULT}")
    print("PASS" if all(ch["pass"] for ch in checks) else "FAIL")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
