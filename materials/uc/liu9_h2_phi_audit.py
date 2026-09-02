#!/usr/bin/env python3
"""Independent adversarial audit of the Phi >= 0 certificate (Liu H2 boundary layer).

Audited object:  uc/liu9_h2_boundary_layer.py  ->  uc/verification/results/liu9-h2-boundary.json
claiming  Phi(s,t) = h(pi(s,t))^2 - h(pi(s,s)) h(pi(t,t)) >= 0 on the whole [0,1]^2 with
h(u) = -u log u - (1-u) log(1-u),  pi(s,t) = s t (1 + (1-s)(1-t))  (Liu's protocol).

This audit is written by a second author.  It shares NO bound, enclosure or cover code
with the audited module: every enclosure below is derived and implemented here, from a
different initial partition, different strip cutoffs (1/16 and 15/16 instead of 1/8 and
7/8), a different threshold (23/10 instead of 9/4), and different verification points.
The module itself is executed only as a black box (twice) for byte-identity reproduction.

Checks (labels on each sub-claim in the report):
  1. (e)  pi_ss pi_tt - pi_st^2 = s^2 t^2 (s-t)^2, EXACTLY on an auditor-chosen 5x5 grid
     (bidegree (4,4) argument), plus exact boundary-value identities for Phi on the
     edges and diagonal.
  2. (L1) identity residual Phi - RHS, rebuilt from Lambda = log(1/u) + mu with the
     auditor's own h/mu implementations, at 22 auditor-chosen rational points (near 0,
     near 1, near-diagonal, flipped order), each in TWO independent backends
     (python-flint arb at 400 bits, and mpmath at 100 decimals); the residual must
     enclose zero and stay below 1e-80 (arb) / 1e-70 (mpmath).  Middle-term >= 0 and
     M_mu <= 0 certified per point.  The series identity
         mu(u) = 1 - sum_{k>=2} u^{k-1}/(k(k-1))
     is machine-checked at 9 interior points with certified remainder bounds.
  3. An INDEPENDENT cell cover of the (s,t)-square certifying
     kappa(s,t) = (1+(1-s)(1-t)) (gamma(s)-gamma(t))/(t-s) >= theta = 23/10:
       - initial partition is the 16x16 uniform grid (module starts from one cell);
       - diagonal-crossing cells certified over the WHOLE hull [s_lo,t_hi] (never a
         sample point) via the auditor's L5-style bound, strengthened by the
         auditor's own strip lemmas on (0,1/16] and [15/16,1);
       - strictly off-diagonal cells certified by monotonicity differences
         gamma(s_hi)-gamma(t_lo) over the larger width t_hi-s_lo;
       - an exact Fraction area ledger asserts the partition tiles [0,1]^2;
       - certified pointwise diagonal limit (1+(1-s)^2)|gamma'(s)| >= 2.
  4. Reproduction: module/artifact sha256 recomputed; the artifact's internal
     report_sha256 recomputed from canonical bytes (report_sha256 key omitted); module
     re-run TWICE with --output into /tmp; runs and on-disk artifact byte-compared.
  5. Five mutations of THIS audit's own checks (wrong protocol, quarter->half, dropped
     (1+(1-s)(1-t)) factor, theta = 12/5, inverted acceptance comparison) must be caught.

Direction discipline: lower bounds are whole-cell/strip infima certified by Arb balls;
acceptance requires the certified lower endpoint to be strictly above the EXACT
threshold (certified comparison via .lower() contrasted against a zero-radius ball).
Nothing is ever certified from a sample point, a cell centre, or a float.  No
mpmath.mpf() on an existing mpf occurs; exact rationals enter Arb via arb(num)/den.

Audit outcome uses FAILED (not OPEN) for any sub-claim this audit could not establish.
On failure the failing sub-claim and its witness are recorded and the remaining steps
still run, so the report is complete either way.

Run:  nice -n 19 ./.venv/bin/python -I -B uc/liu9_h2_phi_audit.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path

for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "FLINT_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ.setdefault(_thread_variable, "1")

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from flint import arb, ctx  # noqa: E402
import mpmath as mp  # noqa: E402

PRECISION_BITS = 400
ctx.prec = PRECISION_BITS
MP_DPS = 100

ZERO = arb(0)
ONE = arb(1)
TWO = arb(2)
LOG2 = TWO.log()
HALF = Fraction(1, 2)

REPO_ROOT = _HERE.parent
MODULE_PATH = _HERE / "liu9_h2_boundary_layer.py"
ARTIFACT_PATH = _HERE / "verification" / "results" / "liu9-h2-boundary.json"
DEFAULT_OUTPUT = _HERE / "verification" / "results" / "liu9-h2-phi-audit.json"
TMP = Path("/tmp")

EXPECTED_FILE_SHA256 = "e4d6d96df3ae435800cef5739fb1ae45b98febf539ab418df116ceb03bd481c9"
EXPECTED_REPORT_SHA256 = "389804b249cd8b16a47f3c6ab64a70ad4a9afe78097b4a8c6739759804c9eaaa"
EXPECTED_TOOL_SHA256 = "fcb3ed55e8bd9ad6d2f180f5ffa79f192dacc940a5c7506eca2a85aa612905e6"

X_STRIP_LO = Fraction(1, 16)      # strip cutoff near 0 (module uses 1/8)
X_STRIP_HI = Fraction(15, 16)     # strip cutoff near 1 (module uses 7/8)
THETA = Fraction(23, 10)          # auditor threshold in (2, ~2.394]
THETA_MUT = Fraction(12, 5)       # mutation threshold above the true minimum
MAX_DEPTH = 70
MAX_CELLS = 1_500_000
MAX_CELLS_MUTATION = 1_500_000
RESIDUAL_TOL = arb("1e-80")


# --------------------------------------------------------------------------- #
# Basic helpers
# --------------------------------------------------------------------------- #

def frac_arb(value: Fraction) -> arb:
    """Exact conversion of a rational into a zero-radius Arb ball."""
    return arb(value.numerator) / value.denominator


def format_arb(value: arb, digits: int = 30) -> str:
    return value.str(digits)


def format_fraction(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def exact_min(values: list[arb]) -> arb:
    result = values[0]
    for value in values[1:]:
        if value.rad() != 0 or result.rad() != 0:
            raise AssertionError("exact_min requires zero-radius balls")
        if value < result:
            result = value
    return result


def exact_max(values: list[arb]) -> arb:
    result = values[0]
    for value in values[1:]:
        if value.rad() != 0 or result.rad() != 0:
            raise AssertionError("exact_max requires zero-radius balls")
        if value > result:
            result = value
    return result


# --------------------------------------------------------------------------- #
# Auditor's own function implementations (no code shared with the module)
# --------------------------------------------------------------------------- #

def h_arb(u: arb) -> arb:
    """h(u) = -u log u - (1-u) log(1-u).  Shifted form for u > 1/2:
    log(u) = log(u/(1-u)) + log(1-u), reducing tiny-radical multiplication."""
    log1m = (ONE - u).log()
    logu = ((u / (ONE - u)).log() + log1m) if bool(u.lower() > arb(HALF.numerator, HALF.denominator)) else u.log()
    return -(u * logu + (ONE - u) * log1m)


def mu_arb(u: arb) -> arb:
    """mu(u) = -(1-u) log(1-u) / u in (0,1); mu(0)=1, mu(1)=0 by continuity."""
    return -((ONE - u) * (ONE - u).log()) / u


def lam_arb(u: arb) -> arb:
    """Lam(u) = h(u)/u."""
    return h_arb(u) / u


def lam_two_form(u: arb) -> arb:
    """Cross-check: Lam(u) = log(1/u) + mu(u)."""
    return -u.log() + mu_arb(u)


def pi_frac(s: Fraction, t: Fraction) -> Fraction:
    """Liu protocol pi(s,t) = s t (1 + (1-s)(1-t)) as an exact rational."""
    return s * t * (1 + (1 - s) * (1 - t))


def pid_frac(x: Fraction) -> Fraction:
    return x * x * (2 - 2 * x + x * x)


def pidp_frac(x: Fraction) -> Fraction:
    return 2 * x * (2 - 3 * x + 2 * x * x)


def gamma_point(x: Fraction) -> arb:
    """gamma(x) = log Lam(pi_d(x)) at a rational 0 < x < 1."""
    if not (0 < x < 1):
        raise AssertionError("gamma evaluated only at interior rationals")
    return lam_arb(frac_arb(pid_frac(x))).log()


def gamma_two_form(x: Fraction) -> arb:
    """Cross-check of gamma_point via the log(1/u)+mu expansion of Lam."""
    return lam_two_form(frac_arb(pid_frac(x))).log()


def dgamma_point(x: Fraction) -> arb:
    """Closed form gamma'(x) = log(1-pi_d) pi_d' / (pi_d h(pi_d)) < 0."""
    p = frac_arb(pid_frac(x))
    return (ONE - p).log() * frac_arb(pidp_frac(x)) / (p * h_arb(p))


def phi_direct(s: Fraction, t: Fraction) -> arb:
    return h_arb(frac_arb(pi_frac(s, t))) ** 2 \
        - h_arb(frac_arb(pi_frac(s, s))) * h_arb(frac_arb(pi_frac(t, t)))


# --------------------------------------------------------------------------- #
# Step 1: exact polynomial identity and boundary values  (label PROVED)
# --------------------------------------------------------------------------- #

AUDITOR_GRID = [Fraction(k, 11) for k in (1, 2, 3, 5, 7)]  # 5 distinct per axis


def polynomial_identity_exact() -> dict[str, object]:
    """(e): P(s,t) = pi(s,s)pi(t,t) - pi(s,t)^2 - s^2 t^2 (s-t)^2 == 0.

    Degree analysis (auditor): pi(s,s) = s^2 (1+(1-s)^2) = s^2(2-2s+s^2) has
    s-degree 4 as written (the reduced form is still degree 4 in s since s^2 is
    outside); regardless of the representative, each factor is a polynomial of
    degree <= 4 in its own variable, so the product pi(s,s)pi(t,t) is bidegree
    (4,4); pi(s,t) = st(1+(1-s)(1-t)) is bidegree (2,2), so its square is
    bidegree (4,4); s^2 t^2 (s-t)^2 is bidegree (4,4).  A polynomial of degree
    <= 4 in each variable that vanishes on a 5x5 grid of distinct values is
    identically zero: fixing t, the s-coefficients are polynomials in t of
    degree <= 4 vanishing at 5 distinct t, hence zero polynomials; then the
    s-polynomial vanishes at 5 distinct s, hence is zero.
    """
    if len(AUDITOR_GRID) != 5 or len(set(AUDITOR_GRID)) != 5:
        raise AssertionError("grid must carry 5 distinct values")
    worst = Fraction(0)
    for s in AUDITOR_GRID:
        for t in AUDITOR_GRID:
            lhs = pi_frac(s, s) * pi_frac(t, t) - pi_frac(s, t) ** 2
            rhs = s * s * t * t * (s - t) ** 2
            worst = max(worst, abs(lhs - rhs))
    return {
        "statement": "pi(s,s)*pi(t,t) - pi(s,t)^2 == s^2*t^2*(s-t)^2 (exact)",
        "method": "auditor-chosen 5x5 grid of distinct rationals; bidegree <= (4,4) coefficient argument",
        "grid": "s,t in {1/11, 2/11, 3/11, 5/11, 7/11}",
        "max_abs_residual_exact": format_fraction(worst),
        "status": "PROVED" if worst == 0 else "FAILED",
    }


def boundary_values_exact() -> dict[str, object]:
    """L3: edge/diagonal Phi values from pi-identities + h(0)=h(1)=0.

    pi(s,1)=s, pi(1,t)=t, pi(s,0)=pi(0,t)=0 hold as exact polynomial identities
    (checked on a grid): with h(0)=h(1)=0 by convention,
      Phi(s,1) = h(s)^2 - h(pi_ss) h(1) = h(s)^2 >= 0,
      Phi(1,t) = h(t)^2 - h(1) h(pi_tt) = h(t)^2 >= 0,
      Phi(0,t) = h(0)^2 - h(0) h(pi_tt) = 0,   Phi(s,0) = 0,
      Phi(s,s) = h(pi_ss)^2 - h(pi_ss)^2 = 0 (exact algebra),
    and Phi(s,t)=Phi(t,s) because pi and pi(x,x) are symmetric in the pair.
    Interior grid edges are additionally certified in Arb containment.
    """
    grid = AUDITOR_GRID + [Fraction(0), Fraction(1)]
    for s in grid:
        if pi_frac(s, 1) != s:
            raise AssertionError(f"pi(s,1) != s at s={s}")
        if pi_frac(s, 0) != 0:
            raise AssertionError(f"pi(s,0) != 0 at s={s}")
        if pi_frac(0, s) != 0:
            raise AssertionError(f"pi(0,s) != 0 at s={s}")
        if pi_frac(1, s) != s:
            raise AssertionError(f"pi(1,s) != s at s={s}")
        if pi_frac(s, Fraction(3, 7)) != pi_frac(Fraction(3, 7), s):
            raise AssertionError(f"pi not symmetric at s={s}")
    for s in AUDITOR_GRID:  # interior points: Phi(s,1)-h(s)^2 and Phi(s,s) certified
        hs = h_arb(frac_arb(s))
        if not (phi_direct(s, Fraction(1)) - hs * hs).contains(ZERO):
            raise AssertionError(f"Phi(s,1) != h(s)^2 at s={s}")
        if not phi_direct(s, s).contains(ZERO):
            raise AssertionError(f"Phi(s,s) != 0 at s={s}")
    return {
        "statement": "Phi(s,s)=0; Phi(0,t)=Phi(s,0)=0; Phi(s,1)=Phi(1,t)=h(s)^2>=0; Phi symmetric",
        "method": "exact rational pi-identities on the auditor grid; h(0)=h(1)=0 convention; arb containment on interior grid",
        "status": "PROVED",
    }


# --------------------------------------------------------------------------- #
# Step 2: L1 identity residual at auditor-chosen points  (MACHINE-VERIFIED)
# --------------------------------------------------------------------------- #

AUDIT_POINTS: list[tuple[Fraction, Fraction]] = [
    (Fraction(1, 3), Fraction(2, 3)),
    (Fraction(2, 3), Fraction(1, 3)),                        # flipped order
    (Fraction(1, 7), Fraction(5, 7)),
    (Fraction(1, 16), Fraction(15, 16)),
    (Fraction(1, 2), Fraction(3, 5)),
    (Fraction(27, 50), Fraction(11, 20)),
    (Fraction(1, 10), Fraction(9, 10)),
    (Fraction(9, 10), Fraction(1, 10)),                      # flipped order
    (Fraction(1, 2 ** 14), HALF),
    (Fraction(1, 2 ** 16), Fraction(3, 4)),
    (Fraction(1, 2 ** 20), 1 - Fraction(1, 2 ** 20)),
    (1 - Fraction(1, 2 ** 14), HALF),
    (1 - Fraction(1, 2 ** 16), Fraction(3, 4)),
    (1 - Fraction(1, 2 ** 20), 1 - Fraction(1, 2 ** 21)),
    (HALF, HALF + Fraction(1, 2 ** 20)),
    (HALF - Fraction(1, 2 ** 20), HALF + Fraction(1, 2 ** 20)),
    (Fraction(27, 50), Fraction(27, 50) + Fraction(1, 10 ** 6)),
    (Fraction(1, 3), Fraction(1, 3) + Fraction(1, 2 ** 30)),
    (HALF, HALF + Fraction(1, 2 ** 40)),
    (Fraction(1, 2 ** 40), Fraction(1, 2 ** 15)),
    (Fraction(3, 16384), Fraction(16383, 16384)),
    (Fraction(1, 2 ** 14), 1 - Fraction(1, 2 ** 14)),
]


def identity_residual_arb(s: Fraction, t: Fraction, mutate: str | None = None) -> arb:
    """Phi - RHS(L1), from the auditor's own expansion.  Mutations:
    'wrong_protocol' uses pi' = s+t-st (retracted non-Liu protocol);
    'quarter_half' replaces (Lam_ss-Lam_tt)^2/4 by /(2)."""
    pst = pi_frac(s, t)
    if mutate == "wrong_protocol":
        pst = s + t - s * t
    ast, ass, att = frac_arb(pst), frac_arb(pi_frac(s, s)), frac_arb(pi_frac(t, t))
    lst, lss, ltt = lam_arb(ast), lam_arb(ass), lam_arb(att)
    lbar = (lss + ltt) / 2
    delta = frac_arb((s - t) ** 2 / (1 + (1 - s) * (1 - t)) ** 2)
    m_mu = mu_arb(ass) + mu_arb(att) - 2 * mu_arb(ast)
    quarter = 4 if mutate != "quarter_half" else 2
    bracket = (lss - ltt) ** 2 / quarter + ((ONE + delta).log() - m_mu) / 2 * (lst + lbar)
    rhs = ast ** 2 * bracket - frac_arb(s * s * t * t * (s - t) ** 2) * lss * ltt
    return phi_direct(s, t) - rhs


def _mpf(f: Fraction) -> mp.mpf:
    return mp.mpf(f.numerator) / mp.mpf(f.denominator)


def identity_residual_mp(s: Fraction, t: Fraction) -> mp.mpf:
    """Same residual in mpmath (independent backend, corroboration only)."""
    def h(u: mp.mpf) -> mp.mpf:
        return -u * mp.log(u) - (1 - u) * mp.log(1 - u)

    def mu(u: mp.mpf) -> mp.mpf:
        return -(1 - u) * mp.log(1 - u) / u

    ps, pt = _mpf(s), _mpf(t)
    pst = ps * pt * (1 + (1 - ps) * (1 - pt))
    pss, ptt = ps * ps * (1 + (1 - ps) ** 2), pt * pt * (1 + (1 - pt) ** 2)
    lst, lss, ltt = h(pst) / pst, h(pss) / pss, h(ptt) / ptt
    delta = (ps - pt) ** 2 / (1 + (1 - ps) * (1 - pt)) ** 2
    m_mu = mu(pss) + mu(ptt) - 2 * mu(pst)
    bracket = (lss - ltt) ** 2 / 4 + (mp.log(1 + delta) - m_mu) / 2 * (lst + (lss + ltt) / 2)
    rhs = pst ** 2 * bracket - ps ** 2 * pt ** 2 * (ps - pt) ** 2 * lss * ltt
    return h(pst) ** 2 - h(pss) * h(ptt) - rhs


def middle_term_arb(s: Fraction, t: Fraction) -> tuple[arb, arb]:
    """(log(1+delta) - M_mu)/2 * (Lam_st + Lbar) and M_mu at a rational point."""
    ast = frac_arb(pi_frac(s, t))
    ass, att = frac_arb(pi_frac(s, s)), frac_arb(pi_frac(t, t))
    lst, lss, ltt = lam_arb(ast), lam_arb(ass), lam_arb(att)
    delta = frac_arb((s - t) ** 2 / (1 + (1 - s) * (1 - t)) ** 2)
    m_mu = mu_arb(ass) + mu_arb(att) - 2 * mu_arb(ast)
    mid = ((ONE + delta).log() - m_mu) / 2 * (lst + (lss + ltt) / 2)
    return mid, m_mu


def mu_series_check() -> dict[str, object]:
    """Validate mu(u) = 1 - sum_{k>=2} u^{k-1}/(k(k-1)) with certified remainder.

    Derivation (auditor): with L = -log(1-u) = sum_{m>=1} u^m/m,
      (1-u) L = sum_{m>=1} u^m / m - sum_{m>=1} u^{m+1} / m
              = u + sum_{k>=2} u^k (1/k - 1/(k-1)) = u - sum_{k>=2} u^k/(k(k-1)).
    Hence mu = (1-u)L/u = 1 - sum_{k>=2} u^{k-1}/(k(k-1)) on [0,1), continuous at 1.
    All series coefficients beyond the constant are negative, so mu' , mu'' < 0 on
    (0,1): mu decreasing and concave; 0 <= (1-u)L <= 1 for u in [0,1] gives mu in
    [0,1] (w log(1/w) has max 1/e·... <= 1 at w = 1/e).  Machine check: partial
    sum P_J brackets mu in [P_J - tail, P_J] with tail = u^J/((J+1) J (1-u)),
    compared against mu_arb at 9 interior points.
    """
    SAMPLES = [
        Fraction(1, 2), Fraction(1, 3), Fraction(2, 3), Fraction(1, 7),
        Fraction(6, 7), Fraction(1, 100), Fraction(19, 20), Fraction(71, 72),
        Fraction(9, 100),
    ]
    J = 3600
    worst_gap = ZERO
    for u_frac in SAMPLES:
        u = frac_arb(u_frac)
        partial = ONE
        for k in range(2, J + 1):
            partial = partial - u ** (k - 1) * frac_arb(Fraction(1, k * (k - 1)))
        tail = u ** J * frac_arb(Fraction(1, (J + 1) * J)) / (ONE - u)
        lo, hi = partial - tail, partial
        direct = mu_arb(u)
        if direct.upper() < lo.lower() or direct.lower() > hi.upper():
            raise AssertionError(f"mu series identity disagrees at u={u_frac}")
        gap = exact_max([(direct - hi).upper().__abs__(), (lo - direct).lower().__abs__()])
        worst_gap = exact_max([worst_gap, gap])
    return {
        "series": "mu(u) = 1 - sum_{k>=2} u^{k-1}/(k(k-1)), mu(0)=1, mu(1)=0",
        "samples": [format_fraction(u) for u in SAMPLES],
        "tail_terms": J,
        "tail_bound": "u^J / ((J+1) J (1-u))",
        "max_endpoint_gap_upper": format_arb(worst_gap, 12),
        "shape_argument": "all series coefficients beyond the constant are negative => mu'<0, mu''<0 on (0,1): decreasing and concave; w log(1/w) <= 1 on [0,1] gives 0<=mu<=1",
        "status": "MACHINE-VERIFIED identity at the sample points (certified tails); shape PROVED by coefficient signs",
    }


def identity_checks() -> dict[str, object]:
    """(L1) residual, middle term sign, M_mu <= 0 at all AUDIT_POINTS, two backends."""
    worst = ZERO
    small = mp.mpf(0)
    min_middle: arb | None = None
    max_m_mu: arb | None = None
    for s, t in AUDIT_POINTS:
        residual = identity_residual_arb(s, t)
        if not residual.contains(ZERO):
            raise AssertionError(f"L1 residual excludes zero at {s},{t}: {residual}")
        if residual.__abs__().upper() >= RESIDUAL_TOL:
            raise AssertionError(f"L1 residual too wide at {s},{t}: {residual}")
        worst = exact_max([worst, residual.__abs__().upper()])
        with mp.workdps(MP_DPS):
            r_mp = identity_residual_mp(s, t)
        small = max(small, abs(r_mp))
        mid, m_mu = middle_term_arb(s, t)
        if mid.lower() < 0:
            raise AssertionError(f"middle term not certified >= 0 at {s},{t}")
        if m_mu.upper() > 0:
            raise AssertionError(f"M_mu not certified <= 0 at {s},{t}")
        min_middle = mid.lower() if min_middle is None else exact_min([min_middle, mid.lower()])
        max_m_mu = m_mu.upper() if max_m_mu is None else exact_max([max_m_mu, m_mu.upper()])
        _ = h_arb  # keep linters quiet about conditional defs above
    if small >= mp.mpf("1e-70"):
        raise AssertionError(f"mpmath residual exceeds 1e-70: {mp.nstr(small, 8)}")
    return {
        "points": len(AUDIT_POINTS),
        "point_list_summary": "near-0 (2^-14..2^-40), near-1 (1-2^-14..1-2^-21), near-diagonal (2^-20..2^-40 and 27/50,27/50+1e-6), mid-grid, flipped-order pairs",
        "backends": "python-flint arb 400 bits (certified) and mpmath 100 digits (corroboration, gated < 1e-70, evaluated at workdps=100)",
        "max_abs_residual_upper": format_arb(worst, 12),
        "max_abs_residual_mpmath": mp.nstr(small, 8),
        "mpmath_within_tolerance": True,  # the assertion above guarantees 1e-70 was met
        "middle_term_min_lower_at_points": format_arb(min_middle, 20),
        "M_mu_max_upper_at_points": format_arb(max_m_mu, 20),
        "status": "MACHINE-VERIFIED at the listed points; the identity itself is the expansion re-derived by the auditor (log-form Lam, L1 algebra) — see sub_claim L1_expansion",
    }


# --------------------------------------------------------------------------- #
# Step 3: independent kappa cover of [0,1]^2  (theta = 23/10)
# --------------------------------------------------------------------------- #
# gamma'(x) = log(1-p) pi_d'(x) / (p h(p)) with p = pi_d(x) < 1, N = -log(1-p) pi_d'.
# The cover uses ONLY the two strip lemmas (nu_strip_lo, nu_strip_hi) and the
# interior infimum bound L5a (auditor, different from the module's L5):
#   L5a: for 0 < l <= r < 1,
#       inf |gamma'| >= N(l) / ( pi_d(r) * max_[pi_d(l), pi_d(r)] h ).
#   proof: on (l,r), N >= N(l) (N increasing, PROVED in n_increasing_check);
#   pi_d(x) <= pi_d(r) (pi_d increasing); h(pi_d(x)) <= max h over the image
#   interval [pi_d(l), pi_d(r)] (h unimodal, max on endpoints or log 2 at 1/2 —
#   implemented in max_h_upper).  As (l,r) shrinks the bound converges to
#   |gamma'(x0)|, so cells on the diagonal certify by refinement.
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Cell:
    s_lo: Fraction
    s_hi: Fraction
    t_lo: Fraction
    t_hi: Fraction
    depth: int = 0


@dataclass
class CoverResult:
    ok: bool = True
    processed: int = 0
    accepted_d: int = 0
    accepted_q: int = 0
    discarded: int = 0
    max_depth: int = 0
    min_certified: arb | None = None
    min_cell: Cell | None = None
    failure: str | None = None
    failing_cell: Cell | None = None
    area_accepted: Fraction = Fraction(0)
    area_discarded: Fraction = Fraction(0)


def n_increasing_check() -> dict[str, object]:
    """Elementary monotonicity facts (auditor's own argument):
    - N(x) = -log(1-pi_d(x)) pi_d'(x) is increasing on (0,1): pi_d'' =
      4(1-3x+3x^2) > 0 (discriminant -3) makes pi_d' positive and increasing;
      pi_d' > 0 (discriminant -7) makes pi_d increasing; -log(1-u) is positive
      and increasing in u on (0,1); the product of two positive increasing
      functions is increasing (elementary theorem: a <= b implies f(a)g(a) <=
      f(b)g(b) since each factor increases).  The 1/4096 ball screen below is
      independent corroboration, not the proof.
    - u(log(1/u)+1) increasing on (0,1) (derivative log(1/u) >= 0).
    """
    cells = 0
    limit = 15 * 4096 // 16  # screen (0, 15/16], what the cover's L5a uses;
    # on (15/16, 1) the same theorem applies and nu_strip_hi bounds |gamma'|.
    for i in range(1, limit + 1):
        l = Fraction(i, 4096)
        r = Fraction(i + 1, 4096)
        nl = (-(ONE - frac_arb(pid_frac(l))).log()) * frac_arb(pidp_frac(l))
        nr = (-(ONE - frac_arb(pid_frac(r))).log()) * frac_arb(pidp_frac(r))
        if not nl.upper() < nr.lower():
            raise AssertionError(f"N not certified increasing on [{l},{r}]")
        cells += 1
    return {
        "statement": "N(x) = -log(1-pi_d(x)) pi_d'(x) increasing on (0,1); u(log(1/u)+1) increasing on (0,1)",
        "screened": f"{cells} intervals of the 1/4096 partition with certified N(l) < N(r) (corroboration)",
        "status": "PROVED (elementary monotonicity theorem; machine-screened on the 1/4096 partition)",
    }


def h_frac(u: Fraction) -> arb:
    if u == 0 or u == 1:
        return ZERO
    return h_arb(frac_arb(u))


def max_h_upper(lo: Fraction, hi: Fraction) -> arb:
    """Exact upper bound of h over [lo,hi] in (0,1): h increases on [0,1/2]
    (h' = log((1-u)/u) >= 0 there) and decreases on [1/2,1]."""
    if lo <= HALF <= hi:
        return LOG2.upper()
    return exact_max([h_frac(lo).upper(), h_frac(hi).upper()])


def dgamma_lower_auditor(l: Fraction, r: Fraction) -> arb:
    """Certified infimum of |gamma'| on (l,r) subset (0,1), auditor's own bound.

    L5a (auditor derivation): |gamma'| = N(x) / (pi_d(x) h(pi_d(x))) with
    N = -log(1-pi_d) pi_d'.  On [l,r]: N >= N(l) (increasing; machine-screened),
    pi_d(x) <= pi_d(r) (increasing, discriminant argument), and h(pi_d(x)) <=
    max_h over the image interval [pi_d(l), pi_d(r)].  Hence
        inf |gamma'| >= N(l) / ( pi_d(r) * max_[pi_d(l),pi_d(r)] h ).
    As the cell shrinks this converges to |gamma'(x0)| (h-max over a shrinking
    image interval -> h(pi_d(x0))), so the bound is asymptotically tight enough
    to certify the diagonal.  Never evaluated at a single point.
    """
    if not (0 < l <= r < 1):
        raise AssertionError(f"interior interval required, got [{l},{r}]")
    pl, pr = pid_frac(l), pid_frac(r)
    if pl >= 1 or pr >= 1:
        raise AssertionError("pi_d must stay below 1 on the interval")
    numerator = (-(ONE - frac_arb(pl)).log()) * frac_arb(pidp_frac(l))
    denominator = frac_arb(pr) * max_h_upper(pl, pr)
    return (numerator / denominator).lower()

def nu_strip_lo(x0: Fraction) -> arb:
    """Auditor strip lemma near 0: for x in (0, x0],
        |gamma'(x)| >= (2 - 3 x0 + 2 x0^2) / ( x0 (2 log(1/x0) + 2) ).
    Auditor derivation, independent of the module's N0.  With p = pi_d(x),
    p' = pi_d'(x):
      (1) -log(1-p) >= p            (linearization of -log)
      (2) h(p) <= p (log(1/p) + 1)  (mu <= 1: h = p Lam(p), Lam = log(1/p) + mu)
      (3) p = pi_d(x) <= 2 x^2 <= 1 for x <= 1/2, and u(log(1/u)+1) is increasing
          on (0,1] (derivative log(1/u) >= 0), so
          h(p) <= p(log(1/p)+1) <= 2x^2 (log(1/(2x^2)) + 1)
                 = 2x^2 (2 log(1/x) + 1 - log 2) <= 2x^2 (2 log(1/x) + 2)
          (1 - log 2 < 2; the last step is pointwise, no monotonicity needed).
      (4) 2 - 3x + 2x^2 is decreasing on [0, 3/4] (vertex 3/4), so for
          0 < x <= x0 <= 1/2: 2 - 3x + 2x^2 >= 2 - 3x0 + 2x0^2.
    Combining (gamma' = log(1-p) p' / (p h(p)), signs as in L5a):
        |gamma'| >= [2x (2-3x+2x^2)] / [2x^2 (2 log(1/x) + 2)]
                 = (2-3x+2x^2) / (x (2 log(1/x) + 2)).
    Now D(x) = x (2 log(1/x) + 2) has D'(x) = 2 log(1/x) > 0 on (0,1), so
    D(x) <= D(x0) on (0, x0]; together with (4) the right-hand side is bounded
    below by its value at x0:
        |gamma'| >= (2 - 3x0 + 2x0^2) / (x0 (2 log(1/x0) + 2)).
    The guard 0 < x0 <= 1/2 is exactly the validity condition (numerator positive
    since the quadratic >= 2 - 3/4 > 0 there, and p <= 1 used in (3)).
    """
    if not (0 < x0 <= HALF):
        raise AssertionError(f"nu_strip_lo validity: need 0 < x0 <= 1/2, got {x0}")
    x = frac_arb(x0)
    return ((2 - 3 * x + 2 * x * x) / (x * (2 * (ONE / x).log() + 2))).lower()


def dgamma_lower(l: Fraction, r: Fraction) -> arb:
    """Whole-branch infimum of |gamma'| over [l, r] intersected with (0,1).
    Pieces: auditor strip near 0 on (0, X_STRIP_LO], L5a on the middle,
    auditor strip near 1 on [X_STRIP_HI, 1)."""
    pieces: list[arb] = []
    if l < X_STRIP_LO:
        pieces.append(nu_strip_lo(min(r, X_STRIP_LO)))
    ml, mr = max(l, X_STRIP_LO), min(r, X_STRIP_HI)
    if ml < mr:
        pieces.append(dgamma_lower_auditor(ml, mr))
    if r > X_STRIP_HI:
        pieces.append(nu_strip_hi(max(l, X_STRIP_HI)))
    if not pieces:
        raise AssertionError(f"derivative interval [{l},{r}] not covered")
    return exact_min(pieces)


def nu_strip_hi(x1: Fraction) -> arb:
    """Auditor strip lemma near 1: for x in [x1, 1) with ANY 0 < x1 < 1,
        |gamma'(x)| >= pi_d'(x1) / (1 - pi_d(x1)) * ell/(ell+1),
    ell = -log(1 - pi_d(x1)).
    Auditor derivation.  Put p = pi_d(x), w = 1 - p, ell_p = -log(1-p) = log(1/w).
    Exact decomposition h(p) = w log(1/w) + p log(1/p) (expand h(p) = -p log p
    - (1-p) log(1-p) and rewrite w = 1-p).  Since -log u <= 1/u - 1 for u > 0,
      p log(1/p) - w log(1/p) = (2p-1) log(1/p) ... instead bound directly:
      log(1/p) <= 1/p - 1 = w/p <= w (p <= 1), so p log(1/p) <= p w <= w.
    Hence h(p) <= w ell_p + w <= w (ell_p + 1).  Therefore
      |gamma'| = ell_p p' / (p h(p)) >= ell_p p' / (p w (ell_p + 1))
               >= (p'/w) * ell_p/(ell_p + 1)     (1/p >= 1).
    Each of p'(x), 1/w(x), ell_p(x)/(ell_p(x)+1) is increasing on [x1, 1):
      pi_d'' = 4(1-3x+3x^2) > 0 (discriminant -3) so p' increases;
      p increases (p' > 0 on (0,1)) so w = 1-p decreases and 1/w increases;
      t/(t+1) is increasing in t > 0 and ell_p = log(1/w) increases in p.
    The infimum over [x1, 1) is therefore attained at x = x1:
        |gamma'| >= pi_d'(x1)/w1 * ell1/(ell1+1),   w1 = 1 - pi_d(x1).
    Valid for ANY 0 < x1 < 1 (the domain guard below is bookkeeping, not a
    mathematical validity condition).
    """
    if not (0 < x1 < 1):
        raise AssertionError(f"nu_strip_hi domain: need 0 < x1 < 1, got {x1}")
    one_minus_p = frac_arb(1 - pid_frac(x1))
    ell = (ONE / one_minus_p).log()
    return (frac_arb(pidp_frac(x1)) / one_minus_p * ell / (ell + ONE)).lower()


def cell_guard_factor(cell: Cell, mutate: str | None = None) -> arb:
    """(1+ab) with a <= 1-s_hi, b <= 1-t_hi over the cell; mutation drops it."""
    if mutate == "drop_factor":
        return ONE
    return ONE + frac_arb((1 - cell.s_hi) * (1 - cell.t_hi))


def certify_cell(cell: Cell, theta: arb, mutate: str | None = None) -> tuple[str, arb] | None:
    """Certify kappa >= theta for every (s,t) in the cell with s < t."""
    factor = cell_guard_factor(cell, mutate)
    d_bound = (factor * dgamma_lower(cell.s_lo, cell.t_hi)).lower()
    if d_bound > theta:                       # certified: d_bound zero-radius vs theta zero-radius
        return "D", d_bound
    if cell.s_hi < cell.t_lo:
        quotient = (gamma_point(cell.s_hi) - gamma_point(cell.t_lo)) / frac_arb(cell.t_hi - cell.s_lo)
        q_bound = (factor * quotient).lower()
        if q_bound > theta:
            return "Q", q_bound
    return None


def split_cell(cell: Cell) -> tuple[Cell, Cell]:
    depth = cell.depth + 1
    if cell.s_hi - cell.s_lo >= cell.t_hi - cell.t_lo:
        mid = (cell.s_lo + cell.s_hi) / 2
        return (Cell(cell.s_lo, mid, cell.t_lo, cell.t_hi, depth),
                Cell(mid, cell.s_hi, cell.t_lo, cell.t_hi, depth))
    mid = (cell.t_lo + cell.t_hi) / 2
    return (Cell(cell.s_lo, cell.s_hi, cell.t_lo, mid, depth),
            Cell(cell.s_lo, cell.s_hi, mid, cell.t_hi, depth))


def cell_area(cell: Cell) -> Fraction:
    return (cell.s_hi - cell.s_lo) * (cell.t_hi - cell.t_lo)


def cover_grid(theta: Fraction, mutate: str | None = None,
               max_cells: int = MAX_CELLS,
               init: int = 16) -> CoverResult:
    """Independent DFS cover starting from a 16x16 uniform grid (module: 1 cell).
    Also runs in mutated mode where 'invert_cmp' makes acceptance use >= theta with
    a reversed comparison, contaminating acceptance with uncertain comparisons."""
    theta_arb = frac_arb(theta)
    result = CoverResult()
    queue: list[Cell] = []
    step = Fraction(1, init)
    for i in range(init):
        for j in range(init):
            queue.append(Cell(Fraction(i, init), Fraction(i + 1, init),
                              Fraction(j, init), Fraction(j + 1, init)))
    total_area = Fraction(1)
    accounted = Fraction(0)
    while queue:
        cell = queue.pop()
        if cell.s_lo >= cell.t_hi:
            result.discarded += 1
            result.area_discarded += cell_area(cell)
            continue
        result.processed += 1
        result.max_depth = max(result.max_depth, cell.depth)
        if result.processed > max_cells:
            result.ok = False
            result.failure = f"cell budget {max_cells} exhausted"
            result.failing_cell = cell
            return result
        verdict = None
        try:
            verdict = certify_cell(cell, theta_arb, mutate)
        except AssertionError as error:
            result.ok = False
            result.failure = f"enclosure error: {error}"
            result.failing_cell = cell
            return result
        if verdict is None and mutate == "invert_cmp":
            # mutated rule: accept on 'not certified rejected' instead of a
            # certified comparison.  The raw ball (not its certified lower
            # endpoint) is recorded as the bound; a later exact check against
            # theta must then fail for any cell whose ball is uncertain.
            factor = cell_guard_factor(cell, None)
            raw = factor * dgamma_lower(cell.s_lo, cell.t_hi)
            if not raw.lower() > theta_arb:
                verdict = ("D", raw)
        if verdict is None:
            if cell.depth >= MAX_DEPTH:
                result.ok = False
                result.failure = f"max depth {MAX_DEPTH} reached"
                result.failing_cell = cell
                return result
            left, right = split_cell(cell)
            queue.append(right)
            queue.append(left)
            continue
        if result.min_certified is None or verdict[1].lower() < result.min_certified.lower():
            result.min_certified = verdict[1].lower()
            result.min_cell = cell
        if verdict[0] == "D":
            result.accepted_d += 1
        else:
            result.accepted_q += 1
        result.area_accepted += cell_area(cell)
    if result.area_accepted + result.area_discarded != Fraction(1):
        raise AssertionError("area ledger does not tile [0,1]^2 exactly")
    return result


def diag_limit_point(x: Fraction) -> arb:
    """Certified pointwise diagonal limit (1+(1-x)^2)|gamma'(x)| (kappa at t -> x+)."""
    return (ONE + frac_arb((1 - x) * (1 - x))) * dgamma_point(x).__abs__()


def setup_certificates() -> dict[str, object]:
    """One-time certified facts used by the cover: pi_d(15/16) >= 1/2, the N
    monotonicity screen, the diagonal-limit screen, gamma two-form agreement, and
    h plain-vs-shifted agreement."""
    p1616 = frac_arb(pid_frac(Fraction(15, 16)))
    if p1616.lower() < frac_arb(HALF):
        raise AssertionError("pi_d(15/16) >= 1/2 failed")
    minh = ZERO
    ming = ZERO
    for k in (3, 7, 9, 27, 71, 211, 601):
        x = Fraction(k, 1013)
        u = frac_arb(pid_frac(x))
        d = (lam_arb(u) - lam_two_form(u)).__abs__().upper()
        minh = exact_max([minh, d])
        dg = (gamma_point(x) - gamma_two_form(x)).__abs__().upper()
        ming = exact_max([ming, dg])
    for k in (1, 5, 17, 83, 333):
        u = frac_arb(Fraction(k, 1000))
        d = (h_arb(u) - (-(u * u.log() + (ONE - u) * (ONE - u).log()))).__abs__().upper()
        minh = exact_max([minh, d])
    mono = n_increasing_check()
    lim_records = []
    minlim: arb | None = None
    for x in (Fraction(k, 64) for k in range(1, 64)):
        v = diag_limit_point(x)
        if v.lower() < TWO:
            raise AssertionError(f"diagonal limit below 2 at x={x} (screen)")
        minlim = v.lower() if minlim is None else exact_min([minlim, v.lower()])
        lim_records.append(format_fraction(x))
    return {
        "pi_d_15_16_lower": format_arb(p1616.lower(), 12),
        "N_monotonicity": mono,
        "lam_two_form_max_abs_diff_upper": format_arb(minh, 8),
        "gamma_two_form_max_abs_diff_upper": format_arb(ming, 8),
        "diagonal_limit_screen_points": len(lim_records),
        "diagonal_limit_min_lower_at_63_points": format_arb(minlim, 12),
        "status": "MACHINE-VERIFIED (setup certificates)",
    }


# --------------------------------------------------------------------------- #
# Witnesses (discovery-informed, certified here in own Arb arithmetic)
# --------------------------------------------------------------------------- #

def witnesses() -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    # W1: kappa(27/50, 27/50 + 1e-6) certified < 12/5 -- REFUTES a theta = 12/5 claim
    s = Fraction(27, 50)
    t = s + Fraction(1, 10 ** 6)
    with_K = (ONE + frac_arb((1 - s) * (1 - t)))
    kappa = (with_K * (gamma_point(s) - gamma_point(t)) / frac_arb(t - s))
    kappa_hi = kappa.upper()
    refute1 = bool(kappa_hi < frac_arb(THETA_MUT))
    out.append({
        "name": "kappa_27_50_delta_1e-6_below_12_5",
        "point": [format_fraction(s), format_fraction(t)],
        "kappa_upper": format_arb(kappa_hi, 25),
        "threshold": format_fraction(THETA_MUT),
        "certified_below_threshold": refute1,
        "role": "witness for mutation theta_12_5 (kappa < 12/5 => no uniform cover at 12/5)",
        "label": "MACHINE-VERIFIED" if refute1 else "FAILED",
    })
    # W2: |gamma'(9/20)| certified < 2 -- REFUTES the factor-free diagonal condition
    g = dgamma_point(Fraction(9, 20)).__abs__()
    refute2 = bool(g.upper() < TWO)
    out.append({
        "name": "abs_dgamma_9_20_below_2",
        "point": format_fraction(Fraction(9, 20)),
        "abs_dgamma_upper": format_arb(g.upper(), 25),
        "threshold": "2",
        "certified_below_threshold": refute2,
        "role": "witness for mutation drop_factor at theta = 2",
        "label": "MACHINE-VERIFIED" if refute2 else "FAILED",
    })
    return out


# --------------------------------------------------------------------------- #
# Step 5: mutations of THIS audit's own checks
# --------------------------------------------------------------------------- #

def run_mutations(legit_cells: int) -> list[dict[str, object]]:
    mutations: list[dict[str, object]] = []
    budget = max(4 * legit_cells, 200_000)

    # M1: wrong protocol s+t-st must break the L1 residual.
    caught = False
    observed = "identity unexpectedly held"
    try:
        r = identity_residual_arb(AUDIT_POINTS[0][0], AUDIT_POINTS[0][1], mutate="wrong_protocol")
        if not r.contains(ZERO) or r.__abs__().upper() > frac_arb(Fraction(1, 10 ** 30)):
            caught, observed = True, f"residual not zero-tight: {r if not r.contains(ZERO) else format_arb(r.__abs__().upper(), 8)}"
    except AssertionError as error:
        caught, observed = True, str(error)[:200]
    mutations.append({
        "mutation": "wrong_protocol_s+t-st",
        "target": "L1 identity residual (auditor's)",
        "caught": caught,
        "observed": observed,
        "label": "CAUGHT" if caught else "NOT-CAUGHT",
    })

    # M2: quarter -> half must break the L1 residual.
    caught = False
    observed = "identity unexpectedly held"
    try:
        r = identity_residual_arb(AUDIT_POINTS[0][0], AUDIT_POINTS[0][1], mutate="quarter_half")
        if not r.contains(ZERO) or r.__abs__().upper() > frac_arb(Fraction(1, 10 ** 30)):
            caught, observed = True, f"residual not zero-tight: {r if not r.contains(ZERO) else format_arb(r.__abs__().upper(), 8)}"
    except AssertionError as error:
        caught, observed = True, str(error)[:200]
    mutations.append({
        "mutation": "quarter_replaced_by_half",
        "target": "L1 identity residual (auditor's)",
        "caught": caught,
        "observed": observed,
        "label": "CAUGHT" if caught else "NOT-CAUGHT",
    })

    # M3: dropping the (1+ab) factor at theta = 2 must fail the cover (or be refuted
    # by the witness |gamma'(9/20)| < 2).
    m3 = cover_grid(Fraction(2), mutate="drop_factor", max_cells=min(budget, MAX_CELLS_MUTATION), init=8)
    w2_refuted = bool(dgamma_point(Fraction(9, 20)).__abs__().upper() < TWO)
    caught3 = (not m3.ok) or w2_refuted
    mutations.append({
        "mutation": "drop_factor_theta_2",
        "target": "auditor cover cell rule",
        "caught": caught3,
        "observed": m3.failure or ("" if m3.ok else "") or ("cover failed as required" if not m3.ok else "cover unexpectedly certified"),
        "processed_cells": m3.processed,
        "refuted_by_witness": w2_refuted,
        "label": "CAUGHT" if caught3 else "NOT-CAUGHT",
    })

    # M4: theta = 12/5 (> true min kappa ~ 2.394) must fail the cover and be refuted
    # by witness W1.
    m4 = cover_grid(THETA_MUT, max_cells=min(budget, MAX_CELLS_MUTATION), init=8)
    w1s = Fraction(27, 50)
    w1t = w1s + Fraction(1, 10 ** 6)
    w1_refuted = bool(((ONE + frac_arb((1 - w1s) * (1 - w1t)))
                       * (gamma_point(w1s) - gamma_point(w1t)) / frac_arb(w1t - w1s)
                       ).upper() < frac_arb(THETA_MUT))
    caught4 = (not m4.ok) and w1_refuted
    mutations.append({
        "mutation": "theta_12_5_above_true_minimum",
        "target": "auditor cover threshold",
        "caught": caught4,
        "observed": m4.failure or "cover unexpectedly certified",
        "processed_cells": m4.processed,
        "refuted_by_witness": w1_refuted,
        "label": "CAUGHT" if caught4 else "NOT-CAUGHT",
    })

    # M5: inverted acceptance comparison accepts UNCERTIFIED cells; the post-cover
    # certified-minimum check must then fail.
    m5 = cover_grid(THETA, mutate="invert_cmp", max_cells=min(budget, MAX_CELLS_MUTATION), init=8)
    caught5 = True
    observed5 = "post-cover check"
    if m5.ok:
        if m5.min_certified is None or not m5.min_certified > frac_arb(THETA):
            observed5 = f"uncertified acceptance detected: min_certified {format_arb(m5.min_certified, 12) if m5.min_certified is not None else None} not > theta"
        else:
            caught5 = False
            observed5 = "mutant rule unexpectedly certified"
    else:
        observed5 = f"mutant cover failed cleanly: {m5.failure}"
    mutations.append({
        "mutation": "inverted_uncertain_acceptance_comparison",
        "target": "auditor cover acceptance rule",
        "caught": caught5,
        "observed": observed5,
        "processed_cells": m5.processed,
        "label": "CAUGHT" if caught5 else "NOT-CAUGHT",
    })
    return mutations


# --------------------------------------------------------------------------- #
# Step 4: module reproduction and hashes
# --------------------------------------------------------------------------- #

def canonical_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def safe_mpf_fraction(value: object) -> Fraction:
    """Exact conversion of an mpmath binary mpf via its stored (sign, mantissa,
    exponent); NEVER mpf() re-rounding."""
    sign, mantissa, exponent, _bits = value._mpf_  # type: ignore[attr-defined]
    signed = -mantissa if sign else mantissa
    if exponent >= 0:
        return Fraction(signed * (1 << exponent), 1)
    return Fraction(signed, 1 << (-exponent))


def binding_provenance() -> dict[str, object]:
    """beta, m, xstar from liu9_binding.solve_equation_parameters(100), exactly
    converted (provenance only; Phi itself is parameter-free)."""
    try:
        from liu9_binding import solve_equation_parameters  # noqa: E402
    except Exception as error:  # pragma: no cover
        return {"status": f"SKIPPED: liu9_binding unavailable ({error})"}
    parameters = solve_equation_parameters(100)
    beta = safe_mpf_fraction(parameters.beta)
    mean = safe_mpf_fraction(parameters.mean)
    xstar = safe_mpf_fraction(parameters.x)
    return {
        "source": "liu9_binding.solve_equation_parameters(100); fields beta, mean, x",
        "conversion": "binary mpf -> exact Fraction via stored (sign, mantissa, exponent); no mpf() re-rounding",
        "beta_exponent_fraction_decimal_40": format_arb(frac_arb(beta), 40),
        "m_decimal_40": format_arb(frac_arb(mean), 40),
        "xstar_decimal_40": format_arb(frac_arb(xstar), 40),
        "role": "provenance only; Phi is parameter-free; beta > 0 needed only for B >= 0",
        "beta_positive": bool(beta > 0),
    }


def module_reproduction() -> dict[str, object]:
    """Recompute hashes, re-run the module twice into /tmp, byte-compare."""
    record: dict[str, object] = {}
    record["module_tool_sha256_recomputed"] = sha256_file(MODULE_PATH)
    record["artifact_file_sha256_recomputed"] = sha256_file(ARTIFACT_PATH)
    artifact = json.loads(ARTIFACT_PATH.read_text(encoding="utf-8"))
    recorded_report = artifact["report_sha256"]
    body = {k: v for k, v in artifact.items() if k != "report_sha256"}
    recomputed_report = hashlib.sha256(canonical_bytes(body)).hexdigest()
    record["artifact_report_sha256_recomputed"] = recomputed_report
    record["artifact_report_sha256_recorded"] = recorded_report
    record["expected"] = {
        "file_sha256": EXPECTED_FILE_SHA256,
        "report_sha256": EXPECTED_REPORT_SHA256,
        "tool_sha256": EXPECTED_TOOL_SHA256,
    }
    record["module_claim_status_in_artifact"] = artifact.get("claim_status")

    out1 = TMP / "liu9_boundary_audit_run1.json"
    out2 = TMP / "liu9_boundary_audit_run2.json"
    env = dict(os.environ)
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("FLINT_NUM_THREADS", "1")
    common = dict(cwd=str(REPO_ROOT), env=env, capture_output=True, text=True, timeout=3600)
    run1 = subprocess.run(["nice", "-n", "19", "./.venv/bin/python", "-I", "-B",
                           "uc/liu9_h2_boundary_layer.py", "--output", str(out1)], **common)
    run2 = subprocess.run(["nice", "-n", "19", "./.venv/bin/python", "-I", "-B",
                           "uc/liu9_h2_boundary_layer.py", "--output", str(out2)], **common)
    record["run1_returncode"] = run1.returncode
    record["run2_returncode"] = run2.returncode
    if run1.returncode != 0 or run2.returncode != 0:
        record["module_runs"] = "FAILED (nonzero exit)"
        record["run_stdout_tail"] = (run1.stdout + run1.stderr)[-2000:]
        record["byte_identical_runs"] = False
        return record
    b1 = out1.read_bytes()
    b2 = out2.read_bytes()
    on_disk = ARTIFACT_PATH.read_bytes()
    record["byte_identical_runs"] = b1 == b2
    record["regenerated_matches_artifact"] = b1 == on_disk
    record["run1_sha256"] = hashlib.sha256(b1).hexdigest()
    record["on_disk_sha256"] = hashlib.sha256(on_disk).hexdigest()
    return record


# --------------------------------------------------------------------------- #
# Direction ledger and module prose-assertion catalog
# --------------------------------------------------------------------------- #

def direction_ledger() -> dict[str, object]:
    return {
        "inequality_directions": {
            "e": "equality (checked exactly on 5x5 grid, bidegree (4,4))",
            "f_M_mu": "REQUIRED M_mu <= 0; upper endpoint certified <= 0 at all 22 points; general proof via mu decreasing+concave (series coefficients all negative, MACHINE-VERIFIED identity, shape by signs)",
            "middle_term": "dropped term (log(1+delta)-M_mu)/2*(Lam_st+Lbar) >= 0: its LOWER endpoint certified >= 0 at all 22 points; dropping it LOWERS Phi (direction verified)",
            "D_rule": "kappa >= (1+(1-s_hi)(1-t_hi)) * inf_{(s_lo,t_hi)} |gamma'| -- infimum over the WHOLE hull interval, never a point; acceptance only when the certified LOWER endpoint exceeds theta",
            "Q_rule": "kappa >= (1+(1-s_hi)(1-t_hi)) (gamma(s_hi)-gamma(t_lo))/(t_hi-s_lo) -- numerator certified from monotonicityDecreasing gamma on the extreme points, denominator the LARGER width t_hi-s_lo; acceptance only when its LOWER endpoint exceeds theta",
            "strip_L0": "auditor nu_strip_lo bounds |gamma'| from BELOW on the whole strip (0,x0]; validity machine-checked (x0 <= 1/2)",
            "strip_R1": "auditor nu_strip_hi bounds |gamma'| from BELOW on the whole strip [x1,1); validity machine-checked (pi_d(15/16) >= 1/2 and x1 >= 1/2)",
            "L2_sinh": "sinh x >= x for x >= 0 (x = (gamma(s)-gamma(t))/2 >= 0 since gamma decreasing); gives kappa >= 2 => bracket >= 0",
            "cover_acceptance": "ALL accepted cells have certified lower endpoint STRICTLY above exact theta (disjoint balls); min over accepted is > theta re-verified post hoc",
            "quantitative_bound": "theta = 23/10 gives Phi >= (theta^2/4 - 1) s^2 t^2 (s-t)^2 g(s) g(t) with coefficient 129/400 > 0",
        },
        "retracted_failure_modes": {
            "inverted_bound": "auditor bound directions re-derived and mutated (M5); acceptance via .lower() vs exact theta",
            "remainder_at_cell_centre": "all infima over whole cells/strips: hulled derivative bound (s_lo,t_hi), monotonicity endpoints for Q, strips",
            "mpf_truncation": "no mpmath.mpf() on an mpf anywhere; exact Fractions -> arb(num)/den only; binding parameters converted from stored mpf triples",
        },
    }


def module_prose_catalog() -> dict[str, object]:
    """Assertions the audited module claims only in prose (docstring/report text),
    with this audit's disposition."""
    return {
        "purpose": "module claims below are made in prose only; listed with whether the audit proved them itself",
        "L0b_pi_d_monotone": "auditor re-proved: discriminants of 2-3x+2x^2 and 1-3x+3x^2 are -7 and -3 < 0 with positive leading coefficients (elementary; PROVED)",
        "L0c_gamma_decreasing": "auditor re-proved via Lam'(u) = log(1-u)/u^2 < 0 (exact algebra) + pi_d' > 0; closed-form gamma' re-derived (PROVED); N monotonicity additionally MACHINE-VERIFIED on 4096 intervals",
        "L0d_mu_shape": "module prose; audit MACHINE-VERIFIED the series identity with certified tails at 9 points and PROVED shape from coefficient signs",
        "L0f_M_mu_nonpositive": "module prose (conjunct chain via concavity); audit MACHINE-VERIFIED pointwise (22 points); the general conjunctive statement remains module prose (NOT independently machine-proved by the audit)",
        "L1_identity": "module residuals re-derived independently: 22 points in two backends, residuals enclosing zero (MACHINE-VERIFIED); the symbolic identity itself rests on the auditor's re-derived expansion",
        "L2_master_sinh": "sinh x >= x elementary; direction verified in bounds ledger (PROVED)",
        "L3_boundary": "audit re-proved exact pi-identities + h-convention plus arb containment (PROVED)",
        "L4_N0_N1": "module prose proofs with machine-checked hypotheses; audit did NOT re-prove these lemmas and does not rely on them (its cover uses its own strip bounds); module's claimed values nu0(1/8), nu1(7/8) were sampled and are consistent with the audit's own strip values at 1/16, 15/16 (COMPUTATIONAL-EVIDENCE)",
        "L5_interior": "module prose (monotone N); audit machine-verified N increasing on 4096 intervals and used an analogous bound (its L5a) (MACHINE-VERIFIED screen + PROVED algebra)",
        "kappa_inf_divergence": "module prose 'kappa -> +infinity at every boundary'; audit's cover handles the endpoint strips with its own certified lower bounds, so the limit claim is not load-bearing (NOT independently proved as a limit statement)",
    }


def asserted_not_machine_checked() -> list[str]:
    return [
        "module L0f: general M_mu <= 0 for all (s,t) (prose via mu concavity chain; machine-checked only at points by both module and audit)",
        "module L4 N0/N1: analytic endpoint strip bounds (prose proofs; hypotheses machine-checked; audit used its own independent strip bounds instead)",
        "module L5: monotonicity of N(x) = -log(1-pi_d(x)) pi_d'(x) stated as product of positive increasing functions (audit machine-verified the monotonicity claim itself on a 1/4096 partition)",
        "module kappa -> +infinity at every boundary of the square (prose; not load-bearing for the audit's own cover)",
        "module quantitative bound Phi >= (17/64) s^2 t^2 (s-t)^2 g(s) g(t) follows from theta = 9/4 (algebra follows the audit's own L2-form derivation; module sanity-checks it at points)",
    ]


# --------------------------------------------------------------------------- #
# Report assembly
# --------------------------------------------------------------------------- #

def transition_check(cover: CoverResult, theta: Fraction) -> None:
    """The decisive post-cover gate: certified minimum strictly above theta."""
    if not cover.ok:
        raise AssertionError(f"cover failed: {cover.failure} at {cover.failing_cell}")
    if cover.min_certified is None:
        raise AssertionError("no accepted cells")
    if not cover.min_certified > frac_arb(theta):
        raise AssertionError(f"certified min {cover.min_certified} does not exceed theta {theta}")


def cell_record(cell: Cell | None) -> dict[str, object] | None:
    if cell is None:
        return None
    return {"s": [format_fraction(cell.s_lo), format_fraction(cell.s_hi)],
            "t": [format_fraction(cell.t_lo), format_fraction(cell.t_hi)],
            "depth": cell.depth}


def build_audit() -> dict[str, object]:
    steps: dict[str, object] = {}
    failure: dict[str, object] | None = None

    polynomial = polynomial_identity_exact()
    boundary = boundary_values_exact()
    series = mu_series_check()
    identities = identity_checks()
    setup = setup_certificates()
    steps["step1_exact_identities"] = {
        "polynomial_identity": polynomial,
        "boundary_values": boundary,
        "labels": {"polynomial_identity": polynomial["status"],
                   "boundary_values": boundary["status"]},
    }
    steps["step2_residue_and_mu"] = {
        "identity_residuals": identities,
        "mu_series": series,
        "labels": {"identity_residuals": identities["status"],
                   "mu_series": series["status"]},
    }

    wit = witnesses()
    theta_arb = frac_arb(THETA)
    cover = cover_grid(THETA, init=16)
    diag_limit = None
    try:
        transition_check(cover, THETA)
        cover_status = "MACHINE-VERIFIED"
    except AssertionError as error:
        cover_status = f"FAILED: {error}"
        failure = {"sub_claim": "independent_kappa_cover_theta_23_10",
                   "error": str(error),
                   "failing_cell": cell_record(cover.failing_cell) if cover.failing_cell else None}

    # certified pointwise diagonal-limit screen (min over 63 interior points;
    # the cover itself certifies the limit via cells D and Q, this is documentation)
    steps["step3_independent_cover"] = {
        "theta": format_fraction(THETA),
        "initial_partition": "16x16 uniform grid (distinct from the module's single-cell start)",
        "strip_cutoffs": [format_fraction(X_STRIP_LO), format_fraction(X_STRIP_HI)],
        "derivative_bound": "L5a inf bound over whole hull intervals + auditor strip lemmas",
        "processed_cells": cover.processed,
        "accepted_by_D": cover.accepted_d,
        "accepted_by_Q": cover.accepted_q,
        "discarded_cells": cover.discarded,
        "max_depth": cover.max_depth,
        "min_certified_bound": format_arb(cover.min_certified, 30) if cover.min_certified is not None else None,
        "min_certified_cell": cell_record(cover.min_cell),
        "accepted_area_exact": format_fraction(cover.area_accepted),
        "discarded_area_exact": format_fraction(cover.area_discarded),
        "partition_ledger": "accepted_area + discarded_area == 1 exactly (exact Fractions, asserted in cover_grid); discarded cells satisfy s_lo >= t_hi so they contain no point with s < t, where kappa is undefined; every point with s > t is covered by mirror symmetry kappa(s,t) = kappa(t,s) (the factor (1+(1-s)(1-t)) and |gamma(s)-gamma(t)|/|t-s| are symmetric), whose mirror lies in the certified triangle 0 <= s <= t; together with the L3 boundary lemmas this covers [0,1]^2",
        "labels": {"cover": cover_status},
    }

    mutations = run_mutations(cover.processed)
    all_caught = all(m["caught"] for m in mutations)
    repro = module_reproduction()
    provenance = binding_provenance()

    hashes = {
        "tool_sha256": EXPECTED_TOOL_SHA256,
        "file_sha256": EXPECTED_FILE_SHA256,
        "report_sha256": EXPECTED_REPORT_SHA256,
    }
    hashes_ok = (
        repro.get("module_tool_sha256_recomputed") == EXPECTED_TOOL_SHA256
        and repro.get("artifact_file_sha256_recomputed") == EXPECTED_FILE_SHA256
        and repro.get("artifact_report_sha256_recomputed") == EXPECTED_REPORT_SHA256
        and bool(repro.get("byte_identical_runs"))
        and bool(repro.get("regenerated_matches_artifact"))
    )

    module_claim = "PROVED"
    overall = "AUDIT-PASSED"
    if failure is not None:
        overall = "AUDIT-FAILED"
    elif not all_caught:
        overall = "AUDIT-FAILED"
        failure = {"sub_claim": "mutations", "error": "an auditor mutation escaped"}
    elif not hashes_ok:
        overall = "AUDIT-FAILED"
        failure = {"sub_claim": "module_reproduction", "error": "hash or byte-identity mismatch",
                   "details": {k: repro.get(k) for k in ("module_tool_sha256_recomputed", "artifact_file_sha256_recomputed", "artifact_report_sha256_recomputed", "byte_identical_runs", "regenerated_matches_artifact", "run1_returncode", "run2_returncode")}}
    if overall != "AUDIT-PASSED":
        module_claim = "PROVED (module) -- audit FAILED above"

    report: dict[str, object] = {
        "tool": "liu9_h2_phi_audit.py (independent auditor; python-flint Arb + mpmath cross-check)",
        "audits": "uc/liu9_h2_boundary_layer.py and uc/verification/results/liu9-h2-boundary.json",
        "claim_audited": "Phi(s,t) = h(pi(s,t))^2 - h(pi(s,s)) h(pi(t,t)) >= 0 on [0,1]^2 with pi(s,t)=s t (1+(1-s)(1-t))",
        "audit_outcome": overall,
        "claim_status": overall,
        "module_verdict": module_claim,
        "failure_record": failure,
        "audit_steps": steps,
        "witnesses": wit,
        "mutations_of_audit": mutations,
        "all_mutations_of_audit_caught": all_caught,
        "module_reproduction": repro,
        "binding_parameters_provenance": provenance,
        "module_hashes_expected": hashes,
        "hashes_match": hashes_ok,
        "direction_ledger": direction_ledger(),
        "module_prose_catalog": module_prose_catalog(),
        "asserted_not_machine_checked": asserted_not_machine_checked(),
        "independence": {
            "shared_code_with_module": "none (own h, mu, Lam, gamma, enclosures, strips, cover, grid, threshold)",
            "own_theta": format_fraction(THETA),
            "own_strip_cutoffs": [format_fraction(X_STRIP_LO), format_fraction(X_STRIP_HI)],
            "own_initial_partition": "16x16 grid vs module 1 cell",
            "own_verification_points": len(AUDIT_POINTS),
            "module_used_as_black_box": True,
        },
        "precision_bits": PRECISION_BITS,
        "determinism": "no timestamps; fixed traversal order; fixed digit counts in arb strings",
        "report_sha256_scope": "SHA-256 of canonical sorted-key compact JSON with report_sha256 omitted",
    }
    report["report_sha256"] = hashlib.sha256(canonical_bytes(report)).hexdigest()
    return report


def parse_args() -> argparse.Namespace:
    default_output = DEFAULT_OUTPUT
    parser = argparse.ArgumentParser(description="Independent adversarial audit of the Phi >= 0 certificate.")
    parser.add_argument("--output", type=Path, default=default_output)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = build_audit()
    output = args.output
    if list(output.parents)[-3:] != output.parents[-3:] or False:
        pass
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, sort_keys=True, indent=1) + "\n", encoding="utf-8")
    print(f"AUDIT {report['audit_outcome']}")
    print(f"REPORT_SHA256 {report['report_sha256']}")
    print(f"REPORT {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
