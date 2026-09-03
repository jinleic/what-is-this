#!/usr/bin/env python3
"""Rigorous Arb certificates: A >= 0 and B >= 0 on [0,1]^2 for the candidate
protocol f(x) = (4/5) x (1-x), at the EXACT candidate balls of
uc/verification/results/liu9-cprime-frontier.json.

Natural logarithms throughout.  Put kappa = (4/5)^2 = 16/25 and

    h(u)      = -u log u - (1-u) log(1-u),        h(0) = h(1) = 0
    pi(s,t)   = s t + f(s) f(t) = s t (1 + kappa (1-s)(1-t))
    q(s)      = pi(s,s) = s^2 (1 + kappa (1-s)^2)
    P2(s,t)   = (1-beta) h(st) - (t h(s) + s h(t)) / (2m)
    phi(s)    = sqrt(beta h(q(s)))
    A(s,t)    = P2(s,t) + phi(s) phi(t)
    B(s,t)    = beta h(pi(s,t)) - phi(s) phi(t)
    Phi(s,t)  = h(pi(s,t))^2 - h(q(s)) h(q(t))          (B >= 0  <=>  Phi >= 0)

beta and m are the frontier's exact candidate balls: x is the unique root of
2x^2 + kappa x^2 (1-x)^2 - 1 in (2/3, 1/sqrt 2), isolated by a 70-digit rational
bracket, p = h(x)/h(x^2), m = p x, and beta is the diagonal stationarity value
of Liu (89)-(90) transcribed to this protocol.  The balls are recomputed here
by the pinned frontier module and compared string-for-string with the pinned
artifact.

A >= 0 is proved with the stratification of uc/liu9_h2_twovar_lemmas.py
(Lemma A2 there, kappa = 1): exact faces, an analytic near-zero strip, an
analytic (1,1) corner, a second-order Taylor certificate at the interior zero
(x,x), and deterministic Arb branch-and-bound on the compact remainder.  The
face coefficient c0 = 1 - beta - 1/(2m) is 0.0196 here against 0.0900 for
Liu's protocol, so the two f-independent analytic strata need epsilon = 1e-40
and delta = 1e-36 instead of 1e-8, and the branch-and-bound must certify
boxes down to those scales.  Two exact rescalings make that cheap:

    Acheck(s,t) = A/(st) = c0 log(1/(st)) + (1-beta) mu(st)
                  - (mu(s)+mu(t))/(2m) + psi(s) psi(t),
        mu(u) = -(1-u) log(1-u)/u,   psi(s) = sqrt(beta Lam(q(s)) (1+kappa(1-s)^2)),
        Lam(u) = h(u)/u = log(1/u) + mu(u);
    Acorner(a,b) = A(1-a,1-b)/u,  u = a+b,  r = ab/u,  alpha = a/u,  W = 1-st = u(1-r):
        Acorner = [c0 + r(1/m - 1 + beta)] log(1/u)
                  + (1-beta)(1-r)(log(1/(1-r)) + mu(W))
                  - [(1-b) alpha (log(1/alpha) + mu(a)) + (1-a)(1-alpha)(log(1/(1-alpha)) + mu(b))]/(2m)
                  + phi(1-a) phi(1-b)/u.

Every term of both rescalings is a product of monotone factors, so exact
endpoint ranges give lower bounds whose loss is first order in the box width
instead of the O(log) cancellation loss of the raw kernel.  Both identities are
verified in Arb at rational points to 1e-60 (the candidate balls have radius
about 1e-70, so residuals of order 1e-68 are the floor).  A box is accepted
when any of the Acheck bound, the Acorner bound, or the centered gradient
bound of liu9_h2_twovar_lemmas.py is certified positive.

Phi >= 0 follows the proof of uc/liu9_h2_boundary_layer.py with kappa carried
through: pi(s,s) pi(t,t) - pi(s,t)^2 = kappa s^2 t^2 (s-t)^2 (exact polynomial
identity), the L1 identity with delta = kappa (s-t)^2/(1+kappa ab)^2, and the
master condition

    Phi >= s^2 t^2 (s-t)^2 g(s) g(t) (kappa_master^2/4 - kappa),
    kappa_master(s,t) = (1 + kappa (1-s)(1-t)) (gamma(s)-gamma(t))/(t-s),

so kappa_master >= theta > 2 sqrt(kappa) = 8/5 suffices; theta = 2 is
certified (numerical minimum 2.275 near s = t = 0.50).  The endpoint bounds
N0, N1 and the interior bound L5 are re-derived with kappa; their proofs are in
the docstrings of nu0, nu1 and dgamma_lower_middle.

Finally the mixture-theorem decomposition
    numerator - (M/m) ehx = sum w_k w_l <nu_k x nu_l, A> + sum w_k <nu_k x nu_k, B> + Var_w <nu_k, phi>
is verified in Arb at random finite mixtures with THIS kernel, so with A >= 0,
B >= 0 the two-protocol inequality numerator >= (M/m) ehx holds for every
finite conditionally i.i.d. mixture under the f = (4/5)x(1-x) protocol; the
passage to general P_U is the Fubini step of liu9_h2_mixture_theorem.py.

Labels: PROVED / MACHINE-VERIFIED.  Run from the repository root:

    nice -n 19 ./.venv/bin/python -I -B uc/liu9_cprime_four_fifths_ab.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import sys
from fractions import Fraction
from pathlib import Path
from typing import Callable, Iterable, Sequence

for _name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "FLINT_NUM_THREADS",
):
    os.environ[_name] = "1"

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from flint import arb, ctx  # noqa: E402

from liu9_cprime_frontier import (  # noqa: E402
    DEFAULT_OUTPUT as FRONTIER_ARTIFACT,
    build_candidate,
    format_arb,
    sha256_file,
)
from liu9_h2_boundary_layer import (  # noqa: E402
    IDENTITY_POINTS,
    Cell,
    CoverResult,
    cell_record,
    exact_max,
    exact_min,
    h_max_upper,
    lam_arb,
    mu_arb,
    split_cell,
)
from liu9_h2_mixture_theorem import Mixture, numerator_ehx  # noqa: E402
from liu9_h2_twovar_lemmas import (  # noqa: E402
    Box,
    CoverStats,
    TaylorCertificate,
    arb_fraction,
    box_record,
    cover_record,
    format_fraction,
    h_interior,
    h_point_fraction,
    h_range_fraction,
    hp,
    hpp,
    hppp,
    hull,
    maximum_abs_lower,
    maximum_abs_upper,
    min_enclosure,
    split_box,
    split_box_dimension,
    sqrt_nonnegative_range,
    stronger_lower,
)

PRECISION_BITS = 400
ctx.prec = PRECISION_BITS

ZERO = arb(0)
ONE = arb(1)
TWO = arb(2)
LOG2 = TWO.log()
SQRT2 = TWO.sqrt()
INV_E = ONE / ONE.exp()

KAPPA = Fraction(16, 25)                  # (4/5)^2
KAPPA_ARB = arb_fraction(KAPPA)
SCALE = Fraction(4, 5)

EPSILON = Fraction(1, 10**40)             # near-zero strip: min(s,t) <= EPSILON
CORNER_DELTA = Fraction(1, 10**36)        # corner square [1-delta,1]^2
LOCAL_HALF = Fraction(16, 1000)           # local square half-width around (x,x)
TAYLOR_BOX_RADIUS = arb("0.03")
MAX_BB_DEPTH = 600
MAX_BB_CELLS = 4_000_000

X0 = Fraction(1, 8)                       # N0 applies on (0, min(t2, X0)]
X1 = Fraction(7, 8)                       # N1 applies on [max(s1, X1), 1)
THETA = Fraction(2)                       # certified kappa_master lower bound, > 8/5
PHI_MAX_DEPTH = 60
PHI_MAX_CELLS = 400_000
RESIDUAL_TOLERANCE = arb("1e-100")        # parameter-free identities (Phi)
BALL_RESIDUAL_TOLERANCE = arb("1e-60")    # identities through the 1e-70-wide candidate balls
MIXTURE_SEED = 20260903

DEFAULT_OUTPUT = HERE / "verification" / "results" / "liu9-cprime-four-fifths-ab.json"

DEPENDENCIES = {
    "frontier_module": {
        "path": "liu9_cprime_frontier.py",
        "sha256": "fcacffc1c33662f12068c4a06054304ba344f56508090e1625178264ac49c8b3",
        "role": "build_candidate: the exact candidate balls (root bracket, p, m, beta, c0)",
    },
    "frontier_artifact": {
        "path": "verification/results/liu9-cprime-frontier.json",
        "sha256": "81ee7a5a76c45c0229b750babd32f2a77aebcdf95f24f2df1fa8af7899257bb0",
        "claim_status": "PROVED",
        "role": "the recorded candidate balls this certificate is compared against string for string",
    },
    "twovar_module": {
        "path": "liu9_h2_twovar_lemmas.py",
        "sha256": "2afe8e242aa5f506c6d54967a133fe246a7e5cb8bad4961b9860e93b86f67163",
        "role": "Lemma A2 machinery: entropy ranges, Taylor and branch-and-bound scaffolding (kappa = 1 there)",
    },
    "boundary_module": {
        "path": "liu9_h2_boundary_layer.py",
        "sha256": "fcb3ed55e8bd9ad6d2f180f5ffa79f192dacc940a5c7506eca2a85aa612905e6",
        "role": "Phi >= 0 machinery: mu, Lam, cell cover scaffolding and identity points (kappa = 1 there)",
    },
    "mixture_module": {
        "path": "liu9_h2_mixture_theorem.py",
        "sha256": "3098a1ca30a0f16582df02e1fb7dfd8fa27970c031167125e696e58ac43dfbab",
        "role": "numerator_ehx with a protocol argument; the (I)-(III) decomposition and its Fubini step",
    },
    "mixture_artifact": {
        "path": "verification/results/liu9-h2-mixture-theorem.json",
        "sha256": "329f7e2d71af8cd78d1a921c72b9ae4d05b71113134eb3341d69932f19ea24b3",
        "claim_status": "PROVED",
        "role": "the kernel-free (I)-(III) identities in the free pairing algebra",
    },
}


# --------------------------------------------------------------------------- #
# Candidate constants: exact balls, cross-checked against the pinned artifact
# --------------------------------------------------------------------------- #

def check_dependencies() -> dict[str, object]:
    out: dict[str, object] = {}
    for name, spec in DEPENDENCIES.items():
        path = (HERE / str(spec["path"])).resolve()
        observed = sha256_file(path)
        if observed != spec["sha256"]:
            raise AssertionError(f"dependency {name} changed: {observed} != {spec['sha256']}")
        record: dict[str, object] = {
            "path": str(path.relative_to(HERE.parent)),
            "sha256": observed,
            "role": spec["role"],
            "match": True,
        }
        if "claim_status" in spec:
            payload = json.loads(path.read_text(encoding="utf-8"))
            if payload.get("claim_status") != spec["claim_status"]:
                raise AssertionError(f"dependency {name} status changed")
            record["claim_status"] = payload["claim_status"]
        out[name] = record
    return out


CANDIDATE = build_candidate("scaled_four_fifths", SCALE)
X = CANDIDATE.x
P = CANDIDATE.p
M = CANDIDATE.m
BETA = CANDIDATE.beta
C0 = CANDIDATE.face_coefficient           # 1 - beta - 1/(2m)
C1 = ONE - ONE / (2 * M)
CORNER_SLOPE = ONE / M - ONE + BETA       # coefficient of r in Acorner's log(1/u) factor
CENTER = 2 * CANDIDATE.root_lo - CANDIDATE.root_hi   # rational, just left of the bracket

if not (BETA > 0 and BETA < 1 and M > 0 and C0 > 0 and C1 > 0 and CORNER_SLOPE > 0):
    raise AssertionError("candidate coefficient signs are not certified")


def artifact_crosscheck() -> dict[str, object]:
    """The balls used here must be the balls the frontier artifact recorded."""
    payload = json.loads(FRONTIER_ARTIFACT.read_text(encoding="utf-8"))
    record = payload["general_f_candidates"]["scaled_four_fifths"]
    equations = record["active_equations"]
    expected = {
        "root_bracket": [str(CANDIDATE.root_lo), str(CANDIDATE.root_hi)],
        "stationarity_beta": format_arb(BETA, 40),
        "m": format_arb(M, 40),
        "one_minus_m": format_arb(CANDIDATE.c, 40),
        "face_coefficient": format_arb(C0, 30),
        "f": CANDIDATE.formula,
    }
    observed = {
        "root_bracket": equations["root_bracket"],
        "stationarity_beta": equations["stationarity_beta"],
        "m": record["m"],
        "one_minus_m": record["one_minus_m"],
        "face_coefficient": record["A_face_check"]["coefficient"],
        "f": record["f"],
    }
    if expected != observed:
        raise AssertionError(f"candidate balls differ from the frontier artifact: {expected} != {observed}")
    if record["remaining_continuum_claim"]["label"] != "OPEN":
        raise AssertionError("frontier artifact no longer records the continuum claim as OPEN")
    return {
        "label": "MACHINE-VERIFIED",
        "artifact": DEPENDENCIES["frontier_artifact"]["path"],
        "compared_fields": sorted(expected),
        "all_equal": True,
        "frontier_open_statement": record["remaining_continuum_claim"]["statement"],
        "frontier_one_minus_m_minus_cprime": record["one_minus_m_minus_cprime"],
    }


# --------------------------------------------------------------------------- #
# Kernel primitives with kappa
# --------------------------------------------------------------------------- #

def pi_frac(s: Fraction, t: Fraction, kappa: Fraction = KAPPA) -> Fraction:
    return s * t * (1 + kappa * (1 - s) * (1 - t))


def pi_arb(s: arb, t: arb, kappa: arb = KAPPA_ARB) -> arb:
    return s * t * (ONE + kappa * (ONE - s) * (ONE - t))


def q_frac(s: Fraction, kappa: Fraction = KAPPA) -> Fraction:
    return s * s * (1 + kappa * (1 - s) ** 2)


def qp_frac(s: Fraction, kappa: Fraction = KAPPA) -> Fraction:
    """q'(s) = 2s (1 + kappa (1-s)(1-2s))."""
    return 2 * s * (1 + kappa * (1 - s) * (1 - 2 * s))


def q_arb(s: arb, kappa: arb = KAPPA_ARB) -> arb:
    return s * s * (ONE + kappa * (ONE - s) ** 2)


def qdiag_derivatives(s: arb, kappa: arb = KAPPA_ARB) -> tuple[arb, arb, arb, arb]:
    """h(q(s)) and its first three s-derivatives."""
    q = q_arb(s, kappa)
    q1 = 2 * s + 2 * kappa * s * (ONE - s) * (ONE - 2 * s)
    q2 = 2 + 2 * kappa * (ONE - 6 * s + 6 * s**2)
    q3 = 2 * kappa * (12 * s - 6)
    return (
        h_interior(q),
        hp(q) * q1,
        hpp(q) * q1**2 + hp(q) * q2,
        hppp(q) * q1**3 + 3 * hpp(q) * q1 * q2 + hp(q) * q3,
    )


def mu_frac(u: Fraction) -> arb:
    """mu(u) = -(1-u) log(1-u)/u on (0,1], mu(1) = 0; mu(0) = 1 by continuity."""
    if u == 0:
        return ONE
    if u == 1:
        return ZERO
    if not 0 < u < 1:
        raise ValueError(f"mu argument outside [0,1]: {u}")
    return mu_arb(arb_fraction(u))


def lam_frac(u: Fraction) -> arb:
    """Lam(u) = h(u)/u on (0,1]; Lam(1) = 0."""
    if not 0 < u <= 1:
        raise ValueError(f"Lam argument outside (0,1]: {u}")
    return h_point_fraction(u) / arb_fraction(u)


def log_inverse_frac(u: Fraction) -> arb:
    if not u > 0:
        raise ValueError(f"log argument not positive: {u}")
    return (ONE / arb_fraction(u)).log()


def entropy_term_frac(alpha: Fraction) -> arb:
    """E(alpha) = alpha log(1/alpha) on [0,1] with E(0) = E(1) = 0."""
    if alpha == 0 or alpha == 1:
        return ZERO
    return arb_fraction(alpha) * log_inverse_frac(alpha)


def entropy_term_max(lo: Fraction, hi: Fraction) -> arb:
    """Upper bound of E on [lo,hi]: E increases on [0,1/e] and decreases after."""
    if not 0 <= lo <= hi <= 1:
        raise ValueError("bad E interval")
    ends = hull(entropy_term_frac(lo), entropy_term_frac(hi)).upper()
    if arb_fraction(hi) < INV_E or arb_fraction(lo) > INV_E:
        return ends
    return INV_E.upper()


def phi_frac(s: Fraction) -> arb:
    if s == 0 or s == 1:
        return ZERO
    return (BETA * h_point_fraction(q_frac(s))).sqrt()


def phi_range(lo: Fraction, hi: Fraction) -> arb:
    """Range of phi on [lo,hi]: q is increasing, h ranges exactly, sqrt monotone."""
    diagonal_entropy = h_range_fraction(q_frac(lo), q_frac(hi))
    return sqrt_nonnegative_range(BETA * diagonal_entropy)


def psi_frac(s: Fraction) -> arb:
    """psi(s) = phi(s)/s = sqrt(beta Lam(q(s)) (1+kappa(1-s)^2)), decreasing on (0,1]."""
    radicand = BETA * lam_frac(q_frac(s)) * (ONE + KAPPA_ARB * arb_fraction((1 - s) ** 2))
    return sqrt_nonnegative_range(radicand)


def p2_point(s: Fraction, t: Fraction) -> arb:
    return (
        (ONE - BETA) * h_point_fraction(s * t)
        - (arb_fraction(t) * h_point_fraction(s) + arb_fraction(s) * h_point_fraction(t)) / (2 * M)
    )


def a_point(s: Fraction, t: Fraction) -> arb:
    return p2_point(s, t) + phi_frac(s) * phi_frac(t)


def r_point(s: Fraction, t: Fraction) -> arb:
    return p2_point(s, t) + BETA * h_point_fraction(pi_frac(s, t))


def b_point(s: Fraction, t: Fraction) -> arb:
    return BETA * h_point_fraction(pi_frac(s, t)) - phi_frac(s) * phi_frac(t)


def acheck_point(s: Fraction, t: Fraction) -> arb:
    """A/(st) through the rescaled identity, at a rational point of (0,1]^2."""
    return (
        C0 * log_inverse_frac(s * t)
        + (ONE - BETA) * mu_frac(s * t)
        - (mu_frac(s) + mu_frac(t)) / (2 * M)
        + psi_frac(s) * psi_frac(t)
    )


def acorner_point(s: Fraction, t: Fraction) -> arb:
    """A/(a+b) through the rescaled identity, a = 1-s, b = 1-t, (a,b) != (0,0)."""
    a, b = 1 - s, 1 - t
    u = a + b
    r = a * b / u
    alpha = a / u
    w = u * (1 - r)
    coefficient = C0 + arb_fraction(r) * CORNER_SLOPE
    positive = (ONE - BETA) * arb_fraction(1 - r) * (log_inverse_frac(1 - r) + mu_frac(w))
    subtracted = (
        arb_fraction(1 - b) * (entropy_term_frac(alpha) + arb_fraction(alpha) * mu_frac(a))
        + arb_fraction(1 - a) * (entropy_term_frac(1 - alpha) + arb_fraction(1 - alpha) * mu_frac(b))
    ) / (2 * M)
    return (
        coefficient * log_inverse_frac(u)
        + positive
        - subtracted
        + phi_frac(s) * phi_frac(t) / arb_fraction(u)
    )


# --------------------------------------------------------------------------- #
# A: identities and the exact interior zero
# --------------------------------------------------------------------------- #

def a_identity_checks() -> dict[str, object]:
    grid = [Fraction(k, 10) for k in range(1, 10)]
    edge = [
        Fraction(1, 10**30), Fraction(1, 10**12), Fraction(1, 1000), Fraction(1, 7),
        Fraction(1, 2), Fraction(9, 10), 1 - Fraction(1, 1000), 1 - Fraction(1, 10**12),
        1 - Fraction(1, 10**30),
    ]
    points = [(s, t) for s in grid for t in grid] + [(s, t) for s in edge for t in edge]
    worst_check = ZERO
    worst_corner = ZERO
    worst_diag = ZERO
    corner_points = 0
    for s, t in points:
        value = a_point(s, t)
        residual_check = value - arb_fraction(s * t) * acheck_point(s, t)
        if not residual_check.contains(0) or not abs(residual_check).upper() < BALL_RESIDUAL_TOLERANCE:
            raise AssertionError(f"Acheck identity fails at {s},{t}: {residual_check}")
        worst_check = exact_max([worst_check, abs(residual_check).upper()])
        if (s, t) != (Fraction(1), Fraction(1)):
            residual_corner = value - arb_fraction(2 - s - t) * acorner_point(s, t)
            if not residual_corner.contains(0) or not abs(residual_corner).upper() < BALL_RESIDUAL_TOLERANCE:
                raise AssertionError(f"Acorner identity fails at {s},{t}: {residual_corner}")
            worst_corner = exact_max([worst_corner, abs(residual_corner).upper()])
            corner_points += 1
        diagonal = a_point(s, s) - r_point(s, s)
        if not diagonal.contains(0):
            raise AssertionError("A(s,s) = R(s,s) fails")
        worst_diag = exact_max([worst_diag, abs(diagonal).upper()])

    worst_face = ZERO
    for u in [Fraction(k, 16) for k in range(17)] + [Fraction(1, 10**20), 1 - Fraction(1, 10**20)]:
        residuals = (
            a_point(Fraction(0), u),
            a_point(u, Fraction(0)),
            a_point(Fraction(1), u) - C0 * h_point_fraction(u),
            a_point(u, Fraction(1)) - C0 * h_point_fraction(u),
            b_point(u, u),
        )
        for residual in residuals:
            if not residual.contains(0):
                raise AssertionError("face identity enclosure misses zero")
            worst_face = exact_max([worst_face, abs(residual).upper()])

    # The candidate equation q(x) + x^2 = 1 (i.e. 2x^2 + f(x)^2 = 1) makes
    # h(q(x)) = h(x^2); with p = h(x)/h(x^2), m = p x this gives A(x,x) = R(x,x) = 0,
    # and beta is defined by d/ds R(s,s) = 0 at x, so A_s = A_t = 0 there.
    quartic_residual = q_arb(X) + X * X - ONE
    at_zero = a_derivatives(X, X)
    if not (quartic_residual.contains(0) and at_zero[0].contains(0)
            and at_zero[1].contains(0) and at_zero[2].contains(0)):
        raise AssertionError("exact interior zero/stationarity enclosure misses zero")
    return {
        "label": "MACHINE-VERIFIED",
        "points": len(points),
        "acheck_identity": "A(s,t) = s t Acheck(s,t)",
        "acheck_max_abs_residual_upper": format_arb(worst_check, 12),
        "acorner_identity": "A(s,t) = (2-s-t) Acorner(s,t)",
        "acorner_points": corner_points,
        "acorner_max_abs_residual_upper": format_arb(worst_corner, 12),
        "diagonal_identity": "A(s,s) = R(s,s) = P2(s,s) + beta h(q(s))",
        "diagonal_max_abs_residual_upper": format_arb(worst_diag, 12),
        "faces": "A(0,t) = A(s,0) = 0; A(1,t) = A(t,1) = c0 h(t); B(s,s) = 0",
        "face_max_abs_residual_upper": format_arb(worst_face, 12),
        "residual_tolerance": "1e-60 (the candidate balls have radius about 1e-70, so Arb cannot cancel them exactly)",
        "interior_zero": {
            "candidate_equation_residual_q(x)+x^2-1": format_arb(quartic_residual, 12),
            "A_at_xx": format_arb(at_zero[0], 12),
            "A_s_at_xx": format_arb(at_zero[1], 12),
            "A_t_at_xx": format_arb(at_zero[2], 12),
            "derivation": (
                "q(x) = 1 - x^2 gives h(q(x)) = h(x^2); p = h(x)/h(x^2) and m = p x give "
                "R(x,x) = 0 = A(x,x); beta is the root of d/ds R(s,s) = 0 at x, and "
                "A_s(s,s) = (1/2) d/ds R(s,s) because phi phi' is half the diagonal derivative of beta h(q)"
            ),
        },
    }


# --------------------------------------------------------------------------- #
# A: derivatives (for the Taylor certificate and the centered gradient bound)
# --------------------------------------------------------------------------- #

def p2_derivatives(s: arb, t: arb) -> tuple[arb, ...]:
    u = s * t
    hu = h_interior(u)
    h1u, h2u, h3u = hp(u), hpp(u), hppp(u)
    hs, ht = h_interior(s), h_interior(t)
    value = (ONE - BETA) * hu - (t * hs + s * ht) / (2 * M)
    ds = (ONE - BETA) * h1u * t - (t * hp(s) + ht) / (2 * M)
    dt = (ONE - BETA) * h1u * s - (s * hp(t) + hs) / (2 * M)
    dss = (ONE - BETA) * h2u * t**2 - t * hpp(s) / (2 * M)
    dtt = (ONE - BETA) * h2u * s**2 - s * hpp(t) / (2 * M)
    dst = (ONE - BETA) * (h2u * s * t + h1u) - (hp(s) + hp(t)) / (2 * M)
    dsss = (ONE - BETA) * h3u * t**3 - t * hppp(s) / (2 * M)
    dttt = (ONE - BETA) * h3u * s**3 - s * hppp(t) / (2 * M)
    dsst = (ONE - BETA) * (h3u * s * t**2 + 2 * h2u * t) - hpp(s) / (2 * M)
    dstt = (ONE - BETA) * (h3u * t * s**2 + 2 * h2u * s) - hpp(t) / (2 * M)
    return value, ds, dt, dss, dst, dtt, dsss, dsst, dstt, dttt


def phi_derivatives(s: arb, kappa: arb = KAPPA_ARB) -> tuple[arb, arb, arb, arb]:
    entropy = qdiag_derivatives(s, kappa)
    w = BETA * entropy[0]
    w1 = BETA * entropy[1]
    w2 = BETA * entropy[2]
    w3 = BETA * entropy[3]
    root = w.sqrt()
    first = w1 / (2 * root)
    second = w2 / (2 * root) - w1**2 / (4 * root**3)
    third = w3 / (2 * root) - 3 * w1 * w2 / (4 * root**3) + 3 * w1**3 / (8 * root**5)
    return root, first, second, third


def a_derivatives(s: arb, t: arb, kappa: arb = KAPPA_ARB) -> tuple[arb, ...]:
    p2 = p2_derivatives(s, t)
    ps = phi_derivatives(s, kappa)
    pt = phi_derivatives(t, kappa)
    return (
        p2[0] + ps[0] * pt[0],
        p2[1] + ps[1] * pt[0],
        p2[2] + ps[0] * pt[1],
        p2[3] + ps[2] * pt[0],
        p2[4] + ps[1] * pt[1],
        p2[5] + ps[0] * pt[2],
        p2[6] + ps[3] * pt[0],
        p2[7] + ps[2] * pt[1],
        p2[8] + ps[1] * pt[2],
        p2[9] + ps[0] * pt[3],
    )


def taylor_certificate(
    derivatives: Callable[[arb, arb], tuple[arb, ...]],
    half_width: Fraction,
    derivative_box_radius: arb = TAYLOR_BOX_RADIUS,
) -> TaylorCertificate:
    """Second-order lower bound around the exact zero (x,x); see liu9_h2_twovar_lemmas."""
    at_zero = derivatives(X, X)
    if not (at_zero[0].contains(0) and at_zero[1].contains(0) and at_zero[2].contains(0)):
        raise AssertionError("interior zero/stationarity enclosure misses zero")
    lambda_s = at_zero[3].lower() - at_zero[4].abs_upper()
    lambda_t = at_zero[5].lower() - at_zero[4].abs_upper()
    lambda_lower = min_enclosure(lambda_s, lambda_t).lower()
    if not lambda_lower > 0:
        raise AssertionError("Hessian lower eigenvalue is not positive")
    derivative_box = hull(X.lower() - derivative_box_radius, X.upper() + derivative_box_radius)
    on_box = derivatives(derivative_box, derivative_box)
    third = tuple(on_box[6:10])
    if not all(value.is_finite() for value in third):
        raise AssertionError("third derivative enclosure is not finite")
    max_third = maximum_abs_upper(third)
    c3 = (2 * SQRT2 * max_third).upper()
    rho = 3 * lambda_lower / c3
    if not (rho > 0 and rho.upper() < derivative_box_radius):
        raise AssertionError("Taylor radius does not lie in derivative box")
    center_error = (X - arb_fraction(CENTER)).abs_upper()
    distance = (SQRT2 * (arb_fraction(half_width) + center_error)).upper()
    coefficient = (lambda_lower / 2 - c3 * distance / 6).lower()
    if not (distance < rho.lower() and coefficient > 0):
        raise AssertionError("chosen local square does not fit positive Taylor ball")
    point_c3_required = maximum_abs_lower((at_zero[6], at_zero[9])).lower()
    return TaylorCertificate(
        value=at_zero[0],
        grad_s=at_zero[1],
        grad_t=at_zero[2],
        h_ss=at_zero[3],
        h_st=at_zero[4],
        h_tt=at_zero[5],
        lambda_lower=lambda_lower,
        third_partial_enclosures=third,
        max_third_upper=max_third,
        c3_upper=c3,
        rho_formula=rho,
        derivative_box_radius=derivative_box_radius,
        chosen_half_width=half_width,
        chosen_distance_upper=distance,
        coefficient_lower=coefficient,
        point_c3_required_lower=point_c3_required,
    )


def taylor_record(certificate: TaylorCertificate) -> dict[str, object]:
    return {
        "value_enclosure": format_arb(certificate.value, 30),
        "gradient_s_enclosure": format_arb(certificate.grad_s, 30),
        "gradient_t_enclosure": format_arb(certificate.grad_t, 30),
        "hessian_ss_enclosure": format_arb(certificate.h_ss, 40),
        "hessian_st_enclosure": format_arb(certificate.h_st, 40),
        "hessian_tt_enclosure": format_arb(certificate.h_tt, 40),
        "lambda_min_lower_enclosure": format_arb(certificate.lambda_lower, 40),
        "third_partial_enclosures": [format_arb(v, 30) for v in certificate.third_partial_enclosures],
        "max_third_partial_upper_enclosure": format_arb(certificate.max_third_upper, 30),
        "c3_upper_enclosure": format_arb(certificate.c3_upper, 30),
        "rho_equals_3lambda_over_c3_enclosure": format_arb(certificate.rho_formula, 30),
        "derivative_box_radius_enclosure": format_arb(certificate.derivative_box_radius, 30),
        "proof_center": "the certified Arb ball bounded by the 70-digit rational root bracket",
        "partition_square_center": format_fraction(CENTER),
        "partition_center_to_x_distance_upper_enclosure": format_arb((X - arb_fraction(CENTER)).abs_upper(), 12),
        "distance_bound_direction": "sqrt(2)*(square half-width + sup|x - partition center|)",
        "chosen_square_half_width": format_fraction(certificate.chosen_half_width),
        "chosen_euclidean_distance_upper_enclosure": format_arb(certificate.chosen_distance_upper, 30),
        "quadratic_minus_remainder_coefficient_lower_enclosure": format_arb(certificate.coefficient_lower, 30),
        "bound": "A >= ||d||^2 (lambda_min/2 - C3 ||d||/6) > 0 on the square minus its center",
    }


# --------------------------------------------------------------------------- #
# A: box lower bounds
# --------------------------------------------------------------------------- #

def scaled_lower(scale_lo: Fraction, scale_hi: Fraction, value_lower: arb) -> arb:
    """Lower bound of scale*value when scale in [scale_lo,scale_hi] >= 0 and value >= value_lower."""
    if value_lower >= 0:
        return (arb_fraction(scale_lo) * value_lower).lower()
    return (arb_fraction(scale_hi) * value_lower).lower()


def acheck_lower(box: Box) -> arb:
    """Certified lower bound of Acheck = A/(st) on a box with s_lo, t_lo > 0.

    log(1/(st)) and mu(st) are decreasing in each coordinate; mu(s), mu(t) are
    decreasing; psi is decreasing and nonnegative, so psi(s)psi(t) >= psi(s_hi)psi(t_hi).
    """
    if not (box.s_lo > 0 and box.t_lo > 0):
        raise ValueError("Acheck needs strictly positive coordinates")
    product_hi = box.s_hi * box.t_hi
    return (
        C0 * log_inverse_frac(product_hi)
        + (ONE - BETA) * mu_frac(product_hi)
        - (mu_frac(box.s_lo) + mu_frac(box.t_lo)) / (2 * M)
        + psi_frac(box.s_hi) * psi_frac(box.t_hi)
    ).lower()


def acorner_lower(box: Box) -> arb | None:
    """Certified lower bound of Acorner = A/(a+b) on a box not touching (1,1).

    With a = 1-s, b = 1-t: u = a+b and W = u - ab increase in a and b; r = ab/u
    increases in a and b; alpha = a/u increases in a and decreases in b; mu is
    decreasing; E(alpha) = alpha log(1/alpha) is bounded through entropy_term_max;
    phi(1-a)phi(1-b) >= 0 uses the exact entropy range on [q(s_lo), q(s_hi)].
    """
    a_lo, a_hi = 1 - box.s_hi, 1 - box.s_lo
    b_lo, b_hi = 1 - box.t_hi, 1 - box.t_lo
    u_lo, u_hi = a_lo + b_lo, a_hi + b_hi
    if u_lo <= 0 or u_hi > 1:
        return None
    r_lo = a_lo * b_lo / (a_lo + b_lo)
    r_hi = a_hi * b_hi / (a_hi + b_hi)
    alpha_lo = a_lo / (a_lo + b_hi) if a_lo > 0 else Fraction(0)
    alpha_hi = a_hi / (a_hi + b_lo) if a_hi > 0 else Fraction(0)
    w_hi = a_hi + b_hi - a_hi * b_hi
    coefficient = C0 + arb_fraction(r_lo) * CORNER_SLOPE
    positive = (ONE - BETA) * arb_fraction(1 - r_hi) * (log_inverse_frac(1 - r_lo) + mu_frac(w_hi))
    subtracted = (
        arb_fraction(1 - b_lo) * (entropy_term_max(alpha_lo, alpha_hi) + arb_fraction(alpha_hi) * mu_frac(a_lo))
        + arb_fraction(1 - a_lo) * (
            entropy_term_max(1 - alpha_hi, 1 - alpha_lo) + arb_fraction(1 - alpha_lo) * mu_frac(b_lo)
        )
    ) / (2 * M)
    phi_product = phi_range(box.s_lo, box.s_hi).lower() * phi_range(box.t_lo, box.t_hi).lower()
    return (
        coefficient * log_inverse_frac(u_hi)
        + positive
        - subtracted
        + phi_product / arb_fraction(u_hi)
    ).lower()


def interval_fraction(lo: Fraction, hi: Fraction) -> arb:
    return hull(arb_fraction(lo), arb_fraction(hi))


def phi_prime_range(lo: Fraction, hi: Fraction) -> arb:
    coordinate = interval_fraction(lo, hi)
    diagonal = interval_fraction(q_frac(lo), q_frac(hi))
    qprime = 2 * coordinate + 2 * KAPPA_ARB * coordinate * (ONE - coordinate) * (ONE - 2 * coordinate)
    return BETA * hp(diagonal) * qprime / (2 * phi_range(lo, hi))


def gradient_lower(box: Box) -> arb | None:
    """Centered first-order lower bound of A (liu9_h2_twovar_lemmas, lemma A2)."""
    if box.s_hi == 1 or box.t_hi == 1:
        return None
    s = interval_fraction(box.s_lo, box.s_hi)
    t = interval_fraction(box.t_lo, box.t_hi)
    st = interval_fraction(box.s_lo * box.t_lo, box.s_hi * box.t_hi)
    h_s = h_range_fraction(box.s_lo, box.s_hi)
    h_t = h_range_fraction(box.t_lo, box.t_hi)
    p2_s = (ONE - BETA) * hp(st) * t - (t * hp(s) + h_t) / (2 * M)
    p2_t = (ONE - BETA) * hp(st) * s - (s * hp(t) + h_s) / (2 * M)
    grad_s = p2_s + phi_prime_range(box.s_lo, box.s_hi) * phi_range(box.t_lo, box.t_hi)
    grad_t = p2_t + phi_range(box.s_lo, box.s_hi) * phi_prime_range(box.t_lo, box.t_hi)
    if not (grad_s.is_finite() and grad_t.is_finite()):
        return None
    midpoint = a_point((box.s_lo + box.s_hi) / 2, (box.t_lo + box.t_hi) / 2)
    half_s = arb_fraction((box.s_hi - box.s_lo) / 2)
    half_t = arb_fraction((box.t_hi - box.t_lo) / 2)
    return (midpoint.lower() - grad_s.abs_upper() * half_s - grad_t.abs_upper() * half_t).lower()


def a_global_lower(box: Box) -> arb:
    """Strongest certified lower bound of A on the box."""
    best = scaled_lower(box.s_lo * box.t_lo, box.s_hi * box.t_hi, acheck_lower(box))
    corner = acorner_lower(box)
    if corner is not None:
        best = stronger_lower(best, scaled_lower(2 - box.s_hi - box.t_hi, 2 - box.s_lo - box.t_lo, corner))
    centered = gradient_lower(box)
    if centered is not None:
        best = stronger_lower(best, centered)
    return best


def split_cover_box(box: Box) -> tuple[Box, Box]:
    """Peel logarithmic near-face strips (liu9_h2_twovar_lemmas.split_cover_box)."""
    lower_coefficient = C0 * log_inverse_frac(box.s_lo * box.t_lo) - ONE / M
    upper_coefficient = C0 * log_inverse_frac(box.s_hi * box.t_hi) - ONE / M
    if lower_coefficient > 0 and not upper_coefficient > 0:
        if box.s_lo * box.t_hi <= box.t_lo * box.s_hi:
            return split_box_dimension(box, "s")
        return split_box_dimension(box, "t")
    return split_box(box)


def cover_boxes(
    roots: Iterable[Box],
    evaluator: Callable[[Box], arb],
    max_depth: int = MAX_BB_DEPTH,
    max_cells: int = MAX_BB_CELLS,
) -> CoverStats:
    stack = list(reversed(list(roots)))
    stats = CoverStats()
    while stack:
        box = stack.pop()
        stats.processed += 1
        stats.max_depth = max(stats.max_depth, box.depth)
        if stats.processed > max_cells:
            raise RuntimeError(f"branch-and-bound exceeded {max_cells} cells at {box}")
        lower = evaluator(box).lower()
        if lower > 0:
            stats.accepted += 1
            if stats.worst_lower is None or lower < stats.worst_lower:
                stats.worst_lower = lower
                stats.worst_box = box
            elif not stats.worst_lower <= lower:
                stats.worst_lower = stats.worst_lower.union(lower).lower()
            continue
        if box.depth >= max_depth:
            raise RuntimeError(
                f"branch-and-bound reached depth cap with nonpositive lower bound {format_arb(lower)} on {box}"
            )
        left, right = split_cover_box(box)
        stats.split += 1
        stack.append(right)
        stack.append(left)
    if stats.accepted == 0 or stats.worst_lower is None:
        raise AssertionError("empty global cover")
    if stats.processed != stats.accepted + stats.split:
        raise AssertionError("branch-and-bound accounting mismatch")
    return stats


def global_root_boxes(local_half: Fraction, corner_delta: Fraction) -> list[Box]:
    local_lo = CENTER - local_half
    local_hi = CENTER + local_half
    corner_start = Fraction(1) - corner_delta
    cuts = [EPSILON, local_lo, local_hi, corner_start, Fraction(1)]
    if cuts != sorted(cuts) or len(set(cuts)) != len(cuts):
        raise AssertionError("invalid stratum partition")
    boxes: list[Box] = []
    for s_lo, s_hi in zip(cuts[:-1], cuts[1:]):
        for t_lo, t_hi in zip(cuts[:-1], cuts[1:]):
            inside_local = s_lo >= local_lo and s_hi <= local_hi and t_lo >= local_lo and t_hi <= local_hi
            inside_corner = s_lo >= corner_start and t_lo >= corner_start
            if inside_local or inside_corner:
                continue
            boxes.append(Box(s_lo, s_hi, t_lo, t_hi))
    if len(boxes) != 14:
        raise AssertionError(f"unexpected global partition size: {len(boxes)}")
    return boxes


def nested_zero_cell(
    evaluator: Callable[[Box], arb],
    initial: Box,
    zero_kind: str,
    depth: int,
) -> tuple[arb, Box]:
    """Shrink a box onto a known exact zero; no sound evaluator can discharge it."""
    box = initial
    enclosure = evaluator(box)
    if enclosure.lower() > 0:
        raise AssertionError("mutant start box was discharged although it contains a zero")
    for _ in range(depth):
        left, right = split_box(box)
        if zero_kind == "x":
            split = left.s_hi if left.s_lo != right.s_lo else left.t_hi
            if X.upper() < arb_fraction(split):
                box = left
            elif X.lower() > arb_fraction(split):
                box = right
            else:
                raise AssertionError("root bracket straddles mutant split")
        elif zero_kind == "corner11":
            box = right
        else:
            raise ValueError(zero_kind)
        enclosure = evaluator(box)
        if enclosure.lower() > 0:
            raise AssertionError("zero-containing mutant cell was incorrectly discharged")
    return enclosure, box


# --------------------------------------------------------------------------- #
# Phi >= 0 with kappa (port of liu9_h2_boundary_layer.py)
# --------------------------------------------------------------------------- #

def gamma_point(x: Fraction) -> arb:
    """gamma(x) = log(h(q(x))/q(x)) at a rational 0 < x < 1."""
    if not (0 < x < 1):
        raise AssertionError("gamma is evaluated only at interior points")
    return lam_arb(arb_fraction(q_frac(x))).log()


def dgamma_point(x: Fraction) -> arb:
    """Closed form gamma'(x) = log(1-q) q'/(q h(q)); sanity samples only."""
    p = arb_fraction(q_frac(x))
    return (ONE - p).log() * arb_fraction(qp_frac(x)) / (p * h_interior(p))


def phi_direct(s: Fraction, t: Fraction) -> arb:
    return h_point_fraction(pi_frac(s, t)) ** 2 - h_point_fraction(q_frac(s)) * h_point_fraction(q_frac(t))


def polynomial_identity_exact(kappa: Fraction = KAPPA) -> dict[str, object]:
    """pi(s,s)pi(t,t) - pi(s,t)^2 == kappa s^2 t^2 (s-t)^2, exactly on a 6x6 grid.

    Both sides have degree <= 4 in s and in t, so vanishing on a 6x6 grid of
    distinct rationals proves the identity.  With a = 1-s, b = 1-t the algebra is
    (1+kappa a^2)(1+kappa b^2) - (1+kappa ab)^2 = kappa (a-b)^2.
    """
    grid = [Fraction(k, 7) for k in range(1, 7)]
    worst = Fraction(0)
    for s in grid:
        for t in grid:
            lhs = pi_frac(s, s) * pi_frac(t, t) - pi_frac(s, t) ** 2
            rhs = kappa * s * s * t * t * (s - t) ** 2
            worst = max(worst, abs(lhs - rhs))
    return {
        "statement": "pi(s,s)*pi(t,t) - pi(s,t)^2 == kappa*s^2*t^2*(s-t)^2 with kappa = 16/25",
        "grid": "s,t in {1/7,...,6/7}",
        "degree_bound": "<= 4 in each variable",
        "max_abs_residual_exact": format_fraction(worst),
        "status": "PROVED" if worst == 0 else "FAILED",
    }


def identity_residual(s: Fraction, t: Fraction, mutate: str | None = None) -> arb:
    """Phi minus the right-hand side of the kappa L1 identity (must enclose zero)."""
    kappa = KAPPA if mutate != "liu_kappa_one_identity" else Fraction(1)
    pst, pss, ptt = pi_frac(s, t), pi_frac(s, s), pi_frac(t, t)
    ast, ass, att = arb_fraction(pst), arb_fraction(pss), arb_fraction(ptt)
    lst, lss, ltt = lam_arb(ast), lam_arb(ass), lam_arb(att)
    delta = arb_fraction(kappa * (s - t) ** 2 / (1 + kappa * (1 - s) * (1 - t)) ** 2)
    m_mu = mu_arb(ass) + mu_arb(att) - 2 * mu_arb(ast)
    quarter = 4 if mutate != "wrong_quarter" else 2
    bracket = (lss - ltt) ** 2 / quarter + (delta.log1p() - m_mu) / 2 * (lst + (lss + ltt) / 2)
    rhs = ast**2 * bracket - arb_fraction(kappa * s * s * t * t * (s - t) ** 2) * lss * ltt
    return phi_direct(s, t) - rhs


def middle_term(s: Fraction, t: Fraction) -> tuple[arb, arb]:
    pst, pss, ptt = pi_frac(s, t), pi_frac(s, s), pi_frac(t, t)
    ast, ass, att = arb_fraction(pst), arb_fraction(pss), arb_fraction(ptt)
    lst, lss, ltt = lam_arb(ast), lam_arb(ass), lam_arb(att)
    delta = arb_fraction(KAPPA * (s - t) ** 2 / (1 + KAPPA * (1 - s) * (1 - t)) ** 2)
    m_mu = mu_arb(ass) + mu_arb(att) - 2 * mu_arb(ast)
    return (delta.log1p() - m_mu) / 2 * (lst + (lss + ltt) / 2), m_mu


def phi_identity_checks(mutate: str | None = None) -> dict[str, object]:
    worst = ZERO
    min_middle: arb | None = None
    max_m_mu: arb | None = None
    for s, t in IDENTITY_POINTS:
        residual = identity_residual(s, t, mutate)
        if not residual.contains(ZERO):
            raise AssertionError(f"L1 identity residual excludes zero at {s},{t}: {residual}")
        if not abs(residual).upper() < RESIDUAL_TOLERANCE:
            raise AssertionError(f"L1 identity residual too wide at {s},{t}: {residual}")
        worst = exact_max([worst, abs(residual).upper()])
        mid, m_mu = middle_term(s, t)
        if not mid.lower() >= 0:
            raise AssertionError(f"middle term not certified nonnegative at {s},{t}")
        if not m_mu.upper() <= 0:
            raise AssertionError(f"M_mu not certified nonpositive at {s},{t}")
        min_middle = mid.lower() if min_middle is None else exact_min([min_middle, mid.lower()])
        max_m_mu = m_mu.upper() if max_m_mu is None else exact_max([max_m_mu, m_mu.upper()])
    return {
        "points": len(IDENTITY_POINTS),
        "max_abs_residual_upper": format_arb(worst, 12),
        "residual_tolerance": "1e-100",
        "middle_term_min_lower_at_points": format_arb(min_middle, 20),
        "M_mu_max_upper_at_points": format_arb(max_m_mu, 20),
        "status": "MACHINE-VERIFIED at the listed rational points; the identity is the algebra of L1 with kappa",
    }


def n0_validity(x0: Fraction) -> bool:
    """x0 < exp(-(1+log(1+kappa))/2), certified in Arb."""
    return bool((arb_fraction(x0).log() + (ONE + (ONE + KAPPA_ARB).log()) / 2).upper() < 0)


def nu0(x0: Fraction, mutate: str | None = None) -> arb:
    """Lower bound of |gamma'| on (0, x0] (N0 with kappa).

    |gamma'(x)| = -log(1-q) q'/(q h(q)).  Numerator: -log(1-q) >= q and
    q'(x) = 2x(1 + kappa(1-x)(1-2x)) >= 2x(1 + kappa(1-3x0)) since (1-x)(1-2x) >= 1-3x.
    Denominator: h(q) <= q(log(1/q)+1) (mu <= 1), q <= (1+kappa)x^2 and u(log(1/u)+1)
    increases on (0,1), so h(q) <= (1+kappa)x^2(2log(1/x) + 1 - log(1+kappa)).
    Hence |gamma'| >= 2(1+kappa(1-3x0)) / ((1+kappa) x (2log(1/x) + 1 - log(1+kappa))),
    and x(2log(1/x)+c) increases for x < exp(-(1+log(1+kappa))/2), so it is at most its
    value at x0.

    `mutate="drop_one_plus_kappa"` omits the (1+kappa) bound on q/x^2 in the
    denominator; the result is then not a lower bound of |gamma'| (mutation).
    """
    if not n0_validity(x0):
        raise AssertionError("N0 hypothesis x0 < exp(-(1+log(1+kappa))/2) is not certified")
    if not x0 <= Fraction(1, 3):
        raise AssertionError("N0 uses 1 + kappa(1-3x0) > 0 through x0 <= 1/3")
    x = arb_fraction(x0)
    one_plus_kappa = ONE + KAPPA_ARB
    scale = ONE if mutate == "drop_one_plus_kappa" else one_plus_kappa
    numerator = 2 * (ONE + KAPPA_ARB * (ONE - 3 * x))
    denominator = scale * x * (2 * (ONE / x).log() + ONE - one_plus_kappa.log())
    return (numerator / denominator).lower()


def nu1(x1: Fraction) -> arb:
    """Lower bound of |gamma'| on [x1, 1) (N1 with kappa).

    q <= 1 and h(q) = h(1-q) <= (1-q)(log(1/(1-q)) + 1); q'' = 2 + 2kappa(1-6x+6x^2)
    >= 2 - kappa > 0, so q' increases and q' >= q'(x1); ell/(ell+1) increases in
    ell = log(1/(1-q)) and 1/(1-q) increases in x.
    """
    if not (0 < x1 < 1):
        raise AssertionError("N1 needs 0 < x1 < 1")
    one_minus_q = arb_fraction(1 - q_frac(x1))
    ell = (ONE / one_minus_q).log()
    return (arb_fraction(qp_frac(x1)) / one_minus_q * ell / (ell + ONE)).lower()


def dgamma_lower_middle(l: Fraction, r: Fraction) -> arb:
    """L5: N(x) = -log(1-q(x)) q'(x) is positive and increasing (q and q' increase),
    so inf_[l,r] |gamma'| >= N(l) / (q(r) max_[q(l),q(r)] h)."""
    if not (0 < l <= r < 1):
        raise AssertionError("middle piece must be strictly interior")
    ql, qr = q_frac(l), q_frac(r)
    numerator = (-(ONE - arb_fraction(ql)).log()) * arb_fraction(qp_frac(l))
    denominator = arb_fraction(qr) * h_max_upper(ql, qr)
    return (numerator / denominator).lower()


def dgamma_lower(l: Fraction, r: Fraction) -> arb:
    if not (0 <= l < r <= 1):
        raise AssertionError("bad derivative interval")
    pieces: list[arb] = []
    if l < X0:
        pieces.append(nu0(min(r, X0)))
    if r > X1:
        pieces.append(nu1(max(l, X1)))
    ml, mr = max(l, X0), min(r, X1)
    if ml < mr:
        pieces.append(dgamma_lower_middle(ml, mr))
    if not pieces:
        raise AssertionError("derivative interval not covered")
    return exact_min(pieces)


def cell_factor(cell: Cell, mutate: str | None) -> arb:
    if mutate == "drop_factor":
        return ONE
    return ONE + KAPPA_ARB * arb_fraction((1 - cell.s_hi) * (1 - cell.t_hi))


def certify_cell(cell: Cell, theta: arb, mutate: str | None) -> tuple[str, arb] | None:
    factor = cell_factor(cell, mutate)
    d_bound = (factor * dgamma_lower(cell.s_lo, cell.t_hi)).lower()
    if d_bound > theta:
        return "D", d_bound
    if cell.s_hi < cell.t_lo:
        quotient = (gamma_point(cell.s_hi) - gamma_point(cell.t_lo)) / arb_fraction(cell.t_hi - cell.s_lo)
        q_bound = (factor * quotient).lower()
        if q_bound > theta:
            return "Q", q_bound
    return None


def cover_triangle(theta: Fraction, mutate: str | None = None, max_cells: int = PHI_MAX_CELLS) -> CoverResult:
    theta_arb = arb_fraction(theta)
    result = CoverResult(ok=True)
    stack = [Cell(Fraction(0), Fraction(1), Fraction(0), Fraction(1))]
    while stack:
        cell = stack.pop()
        if cell.s_lo >= cell.t_hi:
            result.discarded += 1
            continue
        result.processed += 1
        result.max_depth = max(result.max_depth, cell.depth)
        if result.processed > max_cells:
            result.ok = False
            result.failure = f"cell budget {max_cells} exhausted"
            result.failing_cell = cell
            return result
        verdict = certify_cell(cell, theta_arb, mutate)
        if verdict is None:
            if cell.depth >= PHI_MAX_DEPTH:
                result.ok = False
                result.failure = f"max depth {PHI_MAX_DEPTH} reached"
                result.failing_cell = cell
                return result
            left, right = split_cell(cell)
            stack.append(right)
            stack.append(left)
            continue
        kind, bound = verdict
        if kind == "D":
            result.accepted_d += 1
        else:
            result.accepted_q += 1
        if result.min_certified is None or bound < result.min_certified:
            result.min_certified = bound
            result.min_cell = cell
    return result


def kappa_point(s: Fraction, t: Fraction) -> arb:
    return (ONE + KAPPA_ARB * arb_fraction((1 - s) * (1 - t))) * (gamma_point(s) - gamma_point(t)) / arb_fraction(t - s)


def endpoint_samples() -> dict[str, object]:
    samples0 = [Fraction(1, 2**k) for k in (4, 8, 16, 32, 64)]
    samples1 = [1 - Fraction(1, 2**k) for k in (3, 8, 16, 32, 64)]
    bound0, bound1 = nu0(X0), nu1(X1)
    for x in samples0:
        if not abs(dgamma_point(x)).lower() >= bound0:
            raise AssertionError("N0 sample below its bound")
    for x in samples1:
        if not abs(dgamma_point(x)).lower() >= bound1:
            raise AssertionError("N1 sample below its bound")
    return {
        "nu0_at_X0": format_arb(bound0, 20),
        "nu1_at_X1": format_arb(bound1, 20),
        "samples_below_X0": {format_fraction(x): format_arb(abs(dgamma_point(x)), 12) for x in samples0},
        "samples_above_X1": {format_fraction(x): format_arb(abs(dgamma_point(x)), 12) for x in samples1},
        "note": "samples confirm the closed-form gamma' exceeds the analytic bounds; the bounds are proved in nu0/nu1",
    }


def phi_quantitative_samples(theta: Fraction) -> dict[str, object]:
    coefficient = arb_fraction(theta * theta / 4 - KAPPA)
    worst_ratio: arb | None = None
    for s, t in IDENTITY_POINTS:
        gs = lam_arb(arb_fraction(q_frac(s)))
        gt = lam_arb(arb_fraction(q_frac(t)))
        lower = coefficient * arb_fraction(s * s * t * t * (s - t) ** 2) * gs * gt
        phi = phi_direct(s, t)
        if not (phi - lower).lower() >= 0:
            raise AssertionError(f"quantitative bound violated at {s},{t}")
        ratio = (phi / lower).lower()
        worst_ratio = ratio if worst_ratio is None else exact_min([worst_ratio, ratio])
        if not kappa_point(s, t).lower() > arb_fraction(theta):
            raise AssertionError(f"kappa_master sample below theta at {s},{t}")
    return {
        "coefficient_theta^2/4-kappa": format_fraction(theta * theta / 4 - KAPPA),
        "min_ratio_Phi_over_bound_at_points": format_arb(worst_ratio, 15),
        "note": "sanity only; the bound is a theorem once the cover certifies kappa_master >= theta",
    }


# --------------------------------------------------------------------------- #
# Mixture decomposition at THIS kernel
# --------------------------------------------------------------------------- #

def frac_arb(value: Fraction) -> arb:
    return arb_fraction(value)


def law_arb(nu: Sequence[tuple[Fraction, Fraction]]) -> list[tuple[arb, arb]]:
    return [(frac_arb(w), frac_arb(x)) for w, x in nu]


def pair1(nu: list[tuple[arb, arb]], f: Callable[[arb], arb]) -> arb:
    total = ZERO
    for w, x in nu:
        total += w * f(x)
    return total


def pair2(a: list[tuple[arb, arb]], b: list[tuple[arb, arb]], k: Callable[[arb, arb], arb]) -> arb:
    total = ZERO
    for w, x in a:
        for v, y in b:
            total += w * v * k(x, y)
    return total


def h_ball(u: arb) -> arb:
    if u == 0 or u == 1:
        return ZERO
    return h_interior(u)


def phi_ball(x: arb) -> arb:
    value = BETA * h_ball(q_arb(x))
    if value == 0:
        return ZERO
    return value.sqrt()


def kernel_A(x: arb, y: arb) -> arb:
    return (ONE - BETA) * h_ball(x * y) - (y * h_ball(x) + x * h_ball(y)) / (2 * M) + phi_ball(x) * phi_ball(y)


def kernel_B(x: arb, y: arb) -> arb:
    return BETA * h_ball(pi_arb(x, y)) - phi_ball(x) * phi_ball(y)


def nonnegative_form(G: Mixture, mutate: str | None = None) -> tuple[arb, arb, arb, arb]:
    comps = [(frac_arb(w), law_arb(nu)) for w, nu in G]
    F = [pair1(nu, phi_ball) for _, nu in comps]
    sum_a = ZERO
    for wk, nk in comps:
        for wl, nl in comps:
            sum_a += wk * wl * pair2(nk, nl, kernel_A)
    sum_b = ZERO
    if mutate == "mixture_form_B":
        for wk, nk in comps:
            for wl, nl in comps:
                sum_b += wk * wl * pair2(nk, nl, kernel_B)
    else:
        for wk, nk in comps:
            sum_b += wk * pair2(nk, nk, kernel_B)
    mean_f = ZERO
    mean_f2 = ZERO
    for (wk, _), fk in zip(comps, F):
        mean_f += wk * fk
        mean_f2 += wk * fk * fk
    var = ZERO if mutate == "drop_variance" else mean_f2 - mean_f * mean_f
    return sum_a, sum_b, var, sum_a + sum_b + var


def random_mixture(rng: random.Random, components: int, atoms: int) -> Mixture:
    mixture: list[tuple[Fraction, list[tuple[Fraction, Fraction]]]] = []
    weights = [Fraction(rng.randint(1, 20)) for _ in range(components)]
    total = sum(weights)
    for weight in weights:
        masses = [Fraction(rng.randint(1, 9)) for _ in range(atoms)]
        mass_total = sum(masses)
        law = []
        for mass in masses:
            choice = rng.random()
            if choice < 0.08:
                atom = Fraction(0)
            elif choice < 0.16:
                atom = Fraction(1)
            else:
                atom = Fraction(rng.randint(1, 999), 1000)
            law.append((mass / mass_total, atom))
        mixture.append((weight / total, law))
    return mixture


def mixture_decomposition(mutate: str | None = None) -> dict[str, object]:
    rng = random.Random(MIXTURE_SEED)
    worst = ZERO
    trials = 0
    for components in (1, 2, 3, 4, 5):
        for atoms in (1, 2, 3, 4):
            G = random_mixture(rng, components, atoms)
            numerator, ehx = numerator_ehx(G, BETA, pi_arb)
            mean = ZERO
            for w, nu in G:
                for mass, atom in nu:
                    mean += frac_arb(w * mass * atom)
            _, _, _, form = nonnegative_form(G, mutate)
            residual = numerator - mean / M * ehx - form
            if not residual.contains(0):
                raise AssertionError(f"mixture decomposition residual excludes zero at K={components}: {residual}")
            if not abs(residual).upper() < BALL_RESIDUAL_TOLERANCE:
                raise AssertionError(f"mixture decomposition residual too wide at K={components}: {residual}")
            worst = exact_max([worst, abs(residual).upper()])
            trials += 1
    return {
        "label": "MACHINE-VERIFIED",
        "statement": "numerator - (M/m) ehx = sum_kl w_k w_l <nu_k x nu_l, A> + sum_k w_k <nu_k x nu_k, B> + Var_w(<nu_k, phi>)",
        "kernel": "K(x,y) = h(pi(x,y)) with pi = xy + f(x)f(y), f = (4/5)x(1-x)",
        "mixtures": trials,
        "components_up_to": 5,
        "atoms_per_component_up_to": 4,
        "atoms_at_0_and_1_included": True,
        "max_abs_residual_upper": format_arb(worst, 12),
        "residual_tolerance": "1e-60 (the candidate balls have radius about 1e-70)",
        "seed": MIXTURE_SEED,
    }


# --------------------------------------------------------------------------- #
# Mutations
# --------------------------------------------------------------------------- #

def caught(action: Callable[[], object]) -> tuple[bool, str]:
    try:
        action()
    except AssertionError as error:
        return True, str(error)[:200]
    return False, "unexpectedly passed"


def a_mutations(taylor: TaylorCertificate, log_epsilon: arb, strip_margin: arb) -> list[dict[str, object]]:
    mutations: list[dict[str, object]] = []

    local_box = Box(CENTER - LOCAL_HALF, CENTER + LOCAL_HALF, CENTER - LOCAL_HALF, CENTER + LOCAL_HALF)
    local_range, local_cell = nested_zero_cell(a_global_lower, local_box, "x", 14)
    mutations.append({
        "name": "drop_x_local_stratum",
        "status": "FAILED_AS_REQUIRED",
        "sound_local_coefficient_lower_enclosure": format_arb(taylor.coefficient_lower, 30),
        "mutant_depth": local_cell.depth,
        "mutant_zero_cell_lower_bound": format_arb(local_range, 20),
        "mutant_zero_cell": box_record(local_cell),
        "reason": "the nested cell contains A(x,x)=0, so no sound lower bound is positive",
    })

    corner_box = Box(1 - CORNER_DELTA, Fraction(1), 1 - CORNER_DELTA, Fraction(1))
    corner_range, corner_cell = nested_zero_cell(a_global_lower, corner_box, "corner11", 14)
    mutations.append({
        "name": "drop_corner_1_1_stratum",
        "status": "FAILED_AS_REQUIRED",
        "mutant_depth": corner_cell.depth,
        "mutant_zero_cell_lower_bound": format_arb(corner_range, 20),
        "mutant_zero_cell": box_record(corner_cell),
        "reason": "the nested cell contains A(1,1)=0, so no sound lower bound is positive",
    })

    mutant_c3 = (taylor.c3_upper / 16).upper()
    if not mutant_c3 < taylor.point_c3_required_lower:
        raise AssertionError("understated C3 mutant was not refuted at x")
    mutations.append({
        "name": "understate_c3_by_factor_16",
        "status": "FAILED_AS_REQUIRED",
        "sound_c3_upper_enclosure": format_arb(taylor.c3_upper, 30),
        "mutant_c3_enclosure": format_arb(mutant_c3, 30),
        "pointwise_required_c3_lower_enclosure": format_arb(taylor.point_c3_required_lower, 30),
        "reason": "the mutant is below a certified pure-coordinate third derivative at x",
    })

    mutant_strip = (-C0 * log_epsilon - ONE / M).upper()
    if not mutant_strip < 0:
        raise AssertionError("flipped strip-sign mutant unexpectedly passed")
    mutations.append({
        "name": "flip_near_zero_log_term_sign",
        "status": "FAILED_AS_REQUIRED",
        "sound_margin_enclosure": format_arb(strip_margin, 30),
        "mutant_margin_enclosure": format_arb(mutant_strip, 30),
        "reason": "the flipped-sign lower coefficient is certified negative",
    })

    liu_kappa = arb(1)
    was_caught, observed = caught(
        lambda: taylor_certificate(lambda s, t: a_derivatives(s, t, liu_kappa), LOCAL_HALF)
    )
    liu_value = a_derivatives(X, X, liu_kappa)[0]
    if not (was_caught and not liu_value.contains(0)):
        raise AssertionError("Liu-diagonal mutant was not caught")
    mutations.append({
        "name": "liu_kappa_one_diagonal_in_phi",
        "status": "FAILED_AS_REQUIRED",
        "observed": observed,
        "mutant_A_at_xx_enclosure": format_arb(liu_value, 20),
        "reason": "with Liu's diagonal q(s)=s^2(1+(1-s)^2) the candidate point is not a zero of A, so the local stratum cannot start",
    })
    return mutations


def phi_mutations(legit_cells: int) -> list[dict[str, object]]:
    mutations: list[dict[str, object]] = []
    budget = max(4 * legit_cells, 20_000)

    m1 = cover_triangle(Fraction(23, 10), None, max_cells=budget)
    w1s, w1t = Fraction(501, 1000), Fraction(501, 1000) + Fraction(1, 10**6)
    w1 = kappa_point(w1s, w1t)
    refuted1 = bool(w1.upper() < arb_fraction(Fraction(23, 10)))
    if not ((not m1.ok) and refuted1):
        raise AssertionError("theta=23/10 mutant not caught")
    mutations.append({
        "mutation": "theta_2.3_above_true_minimum",
        "theta": "23/10",
        "caught": True,
        "observed": m1.failure,
        "processed_cells": m1.processed,
        "failing_cell": cell_record(m1.failing_cell),
        "refutation_witness": {
            "point": [format_fraction(w1s), format_fraction(w1t)],
            "kappa_master_enclosure": format_arb(w1, 25),
            "certified_below_23/10": refuted1,
        },
    })

    m2 = cover_triangle(THETA, "drop_factor", max_cells=budget)
    w2x = Fraction(432, 1000)
    w2 = abs(dgamma_point(w2x))
    refuted2 = bool(w2.upper() < arb_fraction(THETA))
    if not ((not m2.ok) and refuted2):
        raise AssertionError("drop-factor mutant not caught")
    mutations.append({
        "mutation": "drop_(1+kappa(1-s)(1-t))_factor",
        "theta": format_fraction(THETA),
        "caught": True,
        "observed": m2.failure,
        "processed_cells": m2.processed,
        "failing_cell": cell_record(m2.failing_cell),
        "refutation_witness": {
            "point": format_fraction(w2x),
            "abs_dgamma_enclosure": format_arb(w2, 25),
            "certified_below_theta": refuted2,
        },
        "reason": "the diagonal limit of the factor-free condition is |gamma'(x)| >= theta, which fails at x = 0.432",
    })

    for name, mutate in (
        ("identity_with_liu_kappa_one", "liu_kappa_one_identity"),
        ("identity_quarter_replaced_by_half", "wrong_quarter"),
    ):
        was_caught, observed = caught(lambda mutate=mutate: phi_identity_checks(mutate))
        if not was_caught:
            raise AssertionError(f"{name} not caught")
        mutations.append({"mutation": name, "caught": True, "observed": observed})

    mutant_nu0 = nu0(X0, "drop_one_plus_kappa")
    true_dgamma = abs(dgamma_point(X0))
    refuted3 = bool(mutant_nu0 > true_dgamma.upper())
    if not refuted3:
        raise AssertionError("N0-without-(1+kappa) mutant was not refuted at X0")
    mutations.append({
        "mutation": "N0_without_(1+kappa)_denominator_bound",
        "caught": True,
        "x0": format_fraction(X0),
        "sound_nu0_enclosure": format_arb(nu0(X0), 25),
        "mutant_nu0_enclosure": format_arb(mutant_nu0, 25),
        "true_abs_dgamma_at_x0": format_arb(true_dgamma, 25),
        "certified_above_true_derivative": refuted3,
        "reason": (
            "omitting q <= (1+kappa)x^2 makes the claimed lower bound exceed the exact "
            "|gamma'(x0)|, so it is certified false at x0; the guard n0_validity/x0<=1/3 "
            "is a separate hypothesis check, asserted for the live X0"
        ),
    })
    return mutations


def mixture_mutations() -> list[dict[str, object]]:
    mutations: list[dict[str, object]] = []
    for name in ("drop_variance", "mixture_form_B"):
        was_caught, observed = caught(lambda name=name: mixture_decomposition(name))
        if not was_caught:
            raise AssertionError(f"mixture mutant {name} not caught")
        mutations.append({"mutation": name, "caught": True, "observed": observed})
    return mutations


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #

def canonical_bytes(report: dict[str, object]) -> bytes:
    body = dict(report)
    body.pop("report_sha256", None)
    return (json.dumps(body, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")


def build_report() -> dict[str, object]:
    dependencies = check_dependencies()
    crosscheck = artifact_crosscheck()
    identities = a_identity_checks()

    log_epsilon = log_inverse_frac(EPSILON)
    strip_margin = (C0 * log_epsilon - ONE / M).lower()
    if not strip_margin > 0:
        raise AssertionError("near-zero strip lower bound is not positive")

    delta = arb_fraction(CORNER_DELTA)
    log_corner = (ONE / (2 * delta)).log()
    corner_coefficient = C0 - (ONE - BETA) * delta / 2
    corner_margin = (corner_coefficient * log_corner - (LOG2 + ONE) / (2 * M)).lower()
    if not (corner_coefficient > 0 and corner_margin > 0):
        raise AssertionError("(1,1) corner lower bound failed")

    taylor = taylor_certificate(a_derivatives, LOCAL_HALF)
    cover = cover_boxes(global_root_boxes(LOCAL_HALF, CORNER_DELTA), a_global_lower)
    mutations_a = a_mutations(taylor, log_epsilon, strip_margin)

    polynomial = polynomial_identity_exact()
    if polynomial["status"] != "PROVED":
        raise AssertionError("polynomial identity failed")
    phi_identities = phi_identity_checks()
    endpoints = endpoint_samples()
    if not n0_validity(X0):
        raise AssertionError("X0 outside N0 validity")
    rejected_x0_half, _ = caught(lambda: nu0(Fraction(1, 2)))
    if not rejected_x0_half:
        raise AssertionError("N0 hypothesis guard did not reject x0 = 1/2")
    if not arb_fraction(THETA) > 2 * KAPPA_ARB.sqrt():
        raise AssertionError("theta does not exceed 2 sqrt(kappa)")
    phi_cover = cover_triangle(THETA)
    if not phi_cover.ok:
        raise AssertionError(f"Phi cover failed: {phi_cover.failure} at {cell_record(phi_cover.failing_cell)}")
    if not phi_cover.min_certified > arb_fraction(THETA):
        raise AssertionError("certified minimum does not exceed theta")
    quantitative = phi_quantitative_samples(THETA)
    mutations_phi = phi_mutations(phi_cover.processed)

    decomposition = mixture_decomposition()
    mutations_mixture = mixture_mutations()

    tool_hash = sha256_file(Path(__file__))
    report: dict[str, object] = {
        "artifact": "liu9-cprime-four-fifths-ab",
        "claim": (
            "A(s,t) = P2(s,t) + phi(s)phi(t) >= 0 and B(s,t) = beta h(pi(s,t)) - phi(s)phi(t) >= 0 "
            "on [0,1]^2 for pi(s,t) = st + f(s)f(t), f(x) = (4/5)x(1-x), at the exact candidate balls "
            "(beta, m) of liu9-cprime-frontier.json"
        ),
        "claim_status": "PROVED",
        "label_legend": {
            "PROVED": "a complete certificate: analytic strata plus exhaustive Arb covers with certified lower bounds",
            "MACHINE-VERIFIED": "an exact identity or finite check verified in Arb at the listed points",
        },
        "precision_bits": PRECISION_BITS,
        "natural_log": True,
        "tool": "liu9_cprime_four_fifths_ab.py (python-flint Arb)",
        "tool_sha256": tool_hash,
        "dependencies": dependencies,
        "candidate": {
            "label": "MACHINE-VERIFIED",
            "f": CANDIDATE.formula,
            "kappa": format_fraction(KAPPA),
            "protocol": "pi(s,t) = s t (1 + kappa (1-s)(1-t)); q(s) = pi(s,s) = s^2 (1 + kappa (1-s)^2)",
            "root_polynomial": CANDIDATE.root_polynomial,
            "root_bracket": [str(CANDIDATE.root_lo), str(CANDIDATE.root_hi)],
            "x_enclosure": format_arb(X, 50),
            "p_enclosure": format_arb(P, 50),
            "m_enclosure": format_arb(M, 50),
            "one_minus_m_enclosure": format_arb(CANDIDATE.c, 50),
            "beta_enclosure": format_arb(BETA, 50),
            "c0_equals_1_minus_beta_minus_1_over_2m_enclosure": format_arb(C0, 40),
            "c1_equals_1_minus_1_over_2m_enclosure": format_arb(C1, 40),
            "corner_slope_1_over_m_minus_1_plus_beta_enclosure": format_arb(CORNER_SLOPE, 40),
            "artifact_crosscheck": crosscheck,
        },
        "identities": identities,
        "lemma_A": {
            "claim": "A(s,t) >= 0 on [0,1]^2",
            "status": "PROVED",
            "faces": {
                "s_or_t_zero": "A(0,t) = A(s,0) = 0 exactly (h(0) = 0, phi(0) = 0)",
                "s_or_t_one": "A(1,t) = A(t,1) = c0 h(t) >= 0 with certified c0 > 0 (phi(1) = 0 since q(1) = 1)",
            },
            "near_zero_strip": {
                "epsilon": format_fraction(EPSILON),
                "log_1_over_epsilon_enclosure": format_arb(log_epsilon, 30),
                "lower_coefficient_enclosure": format_arb(strip_margin, 30),
                "bound": (
                    "Dropping phi(s)phi(t) >= 0 gives A >= P2 >= s t (c0 log(1/(st)) - 1/m) "
                    "(h(u) >= u log(1/u) on the positive term, h(u) <= u log(1/u) + u on the subtracted ones). "
                    "If min(s,t) <= epsilon this is >= s t lower_coefficient >= 0; endpoint faces are exact zeros."
                ),
            },
            "corner_1_1": {
                "delta": format_fraction(CORNER_DELTA),
                "log_1_over_2delta_enclosure": format_arb(log_corner, 30),
                "log_coefficient_enclosure": format_arb(corner_coefficient, 30),
                "lower_coefficient_enclosure": format_arb(corner_margin, 30),
                "derivation": (
                    "For a = 1-s, b = 1-t, u = a+b and W = 1-st = u-ab, h(W) >= u(1-delta/2) log(1/u). "
                    "With Aent = a log(1/a) + b log(1/b) <= u(log(1/u)+log 2), the subtracted weighted "
                    "entropies are at most Aent + u.  Dropping phi(s)phi(t) >= 0 gives the bound."
                ),
                "bound": (
                    "For u in (0,2delta], A >= P2 >= u (log_coefficient log(1/u) - (log 2 + 1)/(2m)) "
                    ">= u lower_coefficient >= 0; u = 0 is the exact corner zero."
                ),
            },
            "local_x": taylor_record(taylor),
            "global_branch_and_bound": cover_record(cover),
            "global_lower_bounds": {
                "acheck": (
                    "s_lo t_lo Acheck_lo with Acheck_lo = c0 log(1/(s_hi t_hi)) + (1-beta) mu(s_hi t_hi) "
                    "- (mu(s_lo)+mu(t_lo))/(2m) + psi(s_hi) psi(t_hi)"
                ),
                "acorner": (
                    "u_lo Acorner_lo from the monotone factors of the Acorner identity (docstring of acorner_lower); "
                    "only on boxes with 2 - s_hi - t_hi > 0"
                ),
                "gradient": "A(mid) - sup|A_s| half_s - sup|A_t| half_t on boxes with s_hi < 1 and t_hi < 1",
                "acceptance": "a box is accepted when the strongest of the certified lower bounds is positive",
                "splitting": "the liu9_h2_twovar_lemmas strategy: peel the logarithmic near-face strip, otherwise bisect the longer side",
            },
            "mutations": mutations_a,
        },
        "lemma_Phi": {
            "claim": "Phi(s,t) = h(pi(s,t))^2 - h(q(s)) h(q(t)) >= 0 on [0,1]^2, hence B >= 0 (beta > 0)",
            "status": "PROVED",
            "definitions": {
                "Lam": "h(u)/u = log(1/u) + mu(u)",
                "mu": "-(1-u) log(1-u)/u",
                "g": "Lam(q(x))",
                "gamma": "log g(x)",
                "kappa_master": "(1 + kappa (1-s)(1-t)) (gamma(s)-gamma(t))/(t-s) for 0 < s < t < 1",
            },
            "lemmas": {
                "L0e_polynomial_identity": polynomial,
                "L0b_q_monotone": {
                    "statement": "q'(x) = 2x(1 + kappa(1-x)(1-2x)) >= 2x(1 - kappa/8) > 0 on (0,1]; q''(x) = 2 + 2kappa(1-6x+6x^2) >= 2 - kappa > 0",
                    "proof": "(1-x)(1-2x) >= -1/8 (minimum at x = 3/4) and 1-6x+6x^2 >= -1/2 (minimum at x = 1/2)",
                    "status": "PROVED",
                },
                "L0f_M_mu_nonpositive": {
                    "statement": "mu(q(s)) + mu(q(t)) - 2 mu(pi(s,t)) <= 0",
                    "proof": "L0e gives pi(s,t) <= sqrt(q(s)q(t)) <= (q(s)+q(t))/2; mu decreasing then concave",
                    "status": "PROVED",
                },
                "L1_identity": phi_identities,
                "L1_statement": (
                    "Phi = pi_st^2 [(Lam_ss-Lam_tt)^2/4 + (log(1+delta) - M_mu)/2 (Lam_st + Lbar)] "
                    "- kappa s^2 t^2 (s-t)^2 Lam_ss Lam_tt, delta = kappa (s-t)^2/(1+kappa ab)^2"
                ),
                "L2_master": {
                    "statement": "for 0 < s < t < 1: Phi >= s^2 t^2 (s-t)^2 g(s) g(t) (kappa_master^2/4 - kappa); kappa_master >= 2 sqrt(kappa) implies Phi >= 0",
                    "proof": "drop the nonnegative middle term of L1; g(s)-g(t) = 2 sqrt(g(s)g(t)) sinh((gamma(s)-gamma(t))/2); sinh x >= x",
                    "status": "PROVED",
                },
                "L3_boundary": {
                    "statement": "Phi(s,s) = 0; Phi(0,t) = 0; Phi(s,1) = h(s)^2 >= 0 since pi(s,1) = s and q(1) = 1; Phi symmetric",
                    "status": "PROVED",
                },
                "L4_endpoint_bounds": {
                    "N0": {
                        "statement": "0 < x <= x0 <= 1/3, x0 < exp(-(1+log(1+kappa))/2): |gamma'(x)| >= 2(1+kappa(1-3x0))/((1+kappa) x0 (2log(1/x0)+1-log(1+kappa)))",
                        "x0_used": format_fraction(X0),
                        "validity_certified": n0_validity(X0),
                        "value": endpoints["nu0_at_X0"],
                    },
                    "N1": {
                        "statement": "x1 <= x < 1: |gamma'(x)| >= q'(x1)/(1-q(x1)) ell1/(ell1+1), ell1 = log(1/(1-q(x1)))",
                        "x1_used": format_fraction(X1),
                        "value": endpoints["nu1_at_X1"],
                    },
                    "samples": endpoints,
                    "hypothesis_guards": {
                        "label": "MACHINE-VERIFIED",
                        "asserted_for_X0": "x0 < exp(-(1+log(1+kappa))/2) and x0 <= 1/3",
                        "x0_equals_1/2_rejected": rejected_x0_half,
                        "note": (
                            "these are the hypotheses of the N0 derivation, not a mutation: at "
                            "x0 = 1/2 the displayed bound happens to remain below |gamma'| on "
                            "(0,1/2], so the rejection is a domain guard.  The N0 mutation used "
                            "here instead falsifies the bound itself"
                        ),
                    },
                    "status": "PROVED (proofs in the docstrings of nu0 and nu1)",
                },
                "L5_interior_bound": {
                    "statement": "N(x) = -log(1-q(x)) q'(x) is positive increasing; inf_[l,r] |gamma'| >= N(l)/(q(r) max_[q(l),q(r)] h)",
                    "status": "PROVED",
                },
            },
            "cell_cover": {
                "region": "closed triangle 0 <= s <= t <= 1; only 0 < s < t < 1 needs kappa_master (L3); cells with s_lo >= t_hi are discarded",
                "theta": format_fraction(THETA),
                "two_sqrt_kappa": "8/5",
                "acceptance": "lower endpoint of the Arb bound (D) or (Q) strictly exceeds theta",
                "processed_cells": phi_cover.processed,
                "accepted_by_D": phi_cover.accepted_d,
                "accepted_by_Q": phi_cover.accepted_q,
                "discarded_cells": phi_cover.discarded,
                "max_depth": phi_cover.max_depth,
                "min_certified_bound": format_arb(phi_cover.min_certified, 30),
                "min_certified_cell": cell_record(phi_cover.min_cell),
                "numerical_minimum_of_kappa_master_DISCOVERY": "2.275 near s = t = 0.5015 (float scan, not load-bearing)",
            },
            "quantitative_form": {
                "statement": "for 0 < s < t < 1: Phi >= (theta^2/4 - kappa) s^2 t^2 (s-t)^2 g(s) g(t) = (9/25) s^2 t^2 (s-t)^2 g(s) g(t)",
                "sanity": quantitative,
            },
            "mutations": mutations_phi,
        },
        "lemma_B": {
            "claim": "B(s,t) = beta (h(pi(s,t)) - sqrt(h(q(s)) h(q(t)))) >= 0 on [0,1]^2",
            "status": "PROVED",
            "proof": "beta > 0 (candidate ball) and h >= 0, so Phi >= 0 is equivalent to h(pi_st) >= sqrt(h(q_s) h(q_t))",
            "beta_enclosure": format_arb(BETA, 40),
        },
        "mixture_decomposition_at_candidate_kernel": decomposition,
        "mixture_mutations": mutations_mixture,
        "consequence": {
            "label": "PROVED",
            "statement": (
                "For every finite conditionally i.i.d. mixture Gamma = sum_k w_k nu_k x nu_k of every law mu on [0,1], "
                "under the f = (4/5)x(1-x) protocol, "
                "(1-beta)<mu x mu, h(xy)> + beta sum_k w_k <nu_k x nu_k, h(pi)> >= (M/m) <mu, h> with M = <mu, x>: "
                "the decomposition above is a sum of nonnegative terms (A >= 0, B >= 0 pointwise, variance >= 0). "
                "General P_U follows by the Fubini step of liu9_h2_mixture_theorem.py (bounded continuous kernels)."
            ),
            "union_closed_constant": {
                "value_1_minus_m": format_arb(CANDIDATE.c, 50),
                "exceeds_cprime_by": crosscheck["frontier_one_minus_m_minus_cprime"],
                "route": (
                    "the same route as c' (Liu Proposition 3 with the two-protocol inequality at constant M/m); "
                    "the protocol f = (4/5)x(1-x) satisfies 0 <= f <= min(x,1-x) (frontier general_f_derivation)"
                ),
            },
        },
        "all_mutations_caught": True,
        "determinism": "no timestamps; deterministic depth-first traversals; seeded mixtures; Arb strings at fixed digit counts",
        "report_sha256_scope": "SHA-256 of canonical sorted-key compact JSON with report_sha256 omitted, plus a trailing newline",
    }
    report["report_sha256"] = hashlib.sha256(canonical_bytes(report)).hexdigest()
    return report


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    report = build_report()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, sort_keys=True, indent=1) + "\n", encoding="utf-8")
    cover = report["lemma_A"]["global_branch_and_bound"]
    phi_cover = report["lemma_Phi"]["cell_cover"]
    print(f"LEMMA_A PROVED cells={cover['processed_cells']} depth={cover['max_depth']}")
    print(f"LEMMA_PHI PROVED theta={phi_cover['theta']} cells={phi_cover['processed_cells']} depth={phi_cover['max_depth']}")
    print("LEMMA_B PROVED")
    print(f"ONE_MINUS_M {report['candidate']['one_minus_m_enclosure']}")
    print(f"REPORT {args.output}")
    print(f"REPORT_SHA256 {report['report_sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
