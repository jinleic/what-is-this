#!/usr/bin/env python3
"""The frontier of Liu's two protocols, and discovery beyond it.

Natural logarithms throughout.  Put

    h(u) = -u log u - (1-u) log(1-u),
    pi(s,t) = s t (1 + (1-s)(1-t)),
    N_beta(Gamma) =
      (1-beta)<mu x mu,h(st)>
      + beta int <nu_u x nu_u,h(pi(s,t))> dP_U,
    Gamma = int nu_u x nu_u dP_U,  mu = int nu_u dP_U,
    H(mu) = <mu,h>,  M(mu) = <mu,x>.

For a fixed beta, let m_full(beta) be the infimum of the positive m for which

    N_beta(Gamma) >= (M(mu)/m) H(mu)

holds for every Borel probability law mu on [0,1] and every conditionally
i.i.d. self-coupling Gamma with marginal mu.  (Larger m makes the inequality
weaker.)  The one-component endpoint rays
Gamma=mu x mu, mu=(1-w)delta_0+w delta_x force

    m_full(beta) >= m_end(beta)
      := sup_{0<x<1} x h(x) /
         ((1-beta)h(x^2)+beta h(pi(x,x))).

At Liu's algebraic x*, pi(x*,x*)=1-x*^2, so entropy symmetry makes the value
of this obstruction m*=x*h(x*)/h(x*^2) independent of beta.  Its x-derivative
vanishes only at beta*.  Therefore m_full(beta)>m* for beta != beta*, while the
already-certified mixture theorem gives m_full(beta*)=m*.  Thus beta* is the
unique maximizer of the constant 1-m_full(beta) for this same protocol pair.
No beta-dependent global cover is needed, and the common phrasing
"m(beta)=... at x_beta" is only an endpoint obstruction away from beta*: it is
not an equality for the full measure class without a separate global proof.

The second part derives Liu Example 5's analogous active two-point equations
for K_f(s,t)=h(st+f(s)f(t)), then performs deterministic DISCOVERY scans of the
families requested in the research ledger.  Two candidates above c' are
examined at 400 bits on a rational point lattice and a coarse interval
partition.  These finite diagnostics are COMPUTATIONAL-EVIDENCE only.  In
particular, unresolved interval cells are never promoted to a proof.

The last part transcribes the exact three-protocol obligation from Liu's
Lemma 8 and records it OPEN.  It is not attempted here.

Labels used by this artifact: PROVED / MACHINE-VERIFIED /
COMPUTATIONAL-EVIDENCE / OPEN.

Run from the repository root:

    nice -n 19 ./.venv/bin/python -I -B uc/liu9_cprime_frontier.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Callable, Sequence

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

import mpmath as mp  # noqa: E402
import sympy as sp  # noqa: E402
from flint import arb, ctx  # noqa: E402

from liu9_binding import (  # noqa: E402
    ArbParameters,
    certify_equation_parameters,
    solve_equation_parameters,
)

PRECISION_BITS = 400
MP_DPS = 100
ROOT_DIGITS = 70
POINT_DENOMINATOR = 64
COVER_DENOMINATOR = 64
DEFAULT_OUTPUT = HERE / "verification" / "results" / "liu9-cprime-frontier.json"
ALLOWED_LABELS = {"PROVED", "MACHINE-VERIFIED", "COMPUTATIONAL-EVIDENCE", "OPEN"}
EXPECTED_DIAGNOSTIC_COUNTS = {
    "scaled_four_fifths": {
        "point_lattice": {
            "A": {"certified_nonnegative": 4096, "certified_negative": 0, "unresolved": 0},
            "B": {"certified_nonnegative": 4160, "certified_negative": 0, "unresolved": 0},
        },
        "coarse_partition": {
            "A": {"certified_nonnegative": 456, "certified_negative": 0, "unresolved": 3640},
            "B": {"certified_nonnegative": 2150, "certified_negative": 0, "unresolved": 1946},
        },
    },
    "x_squared_one_minus_x": {
        "point_lattice": {
            "A": {"certified_nonnegative": 3970, "certified_negative": 126, "unresolved": 0},
            "B": {"certified_nonnegative": 4160, "certified_negative": 0, "unresolved": 0},
        },
        "coarse_partition": {
            "A": {"certified_nonnegative": 0, "certified_negative": 0, "unresolved": 4096},
            "B": {"certified_nonnegative": 2232, "certified_negative": 0, "unresolved": 1864},
        },
    },
}

ONE = arb(1)
HALF = Fraction(1, 2)
LOG2 = arb(2).log()

DEPENDENCIES = {
    "mixture_theorem": {
        "path": "verification/results/liu9-h2-mixture-theorem.json",
        "sha256": "329f7e2d71af8cd78d1a921c72b9ae4d05b71113134eb3341d69932f19ea24b3",
        "claim_status": "PROVED",
        "role": "attainment: m_full(beta*) <= m* for every Borel law and conditionally i.i.d. self-coupling",
    },
    "binding_source": {
        "path": "liu9_binding.py",
        "sha256": "feb3a1aaded3cff0805d852c4349ff7142c55e24be4b2da3c0a27f8e7f58cac6",
        "role": "exact root bracket and propagation of Liu (87)-(90)",
    },
    "objective_source": {
        "path": "liu9_objective.py",
        "sha256": "f043b8c99cb26c62804c6f4f29190065a7a876e47d5e7c27c9e6f555c4098c6c",
        "role": "audited transcription of the same two protocol weights",
    },
    "liu_source": {
        "path": "../LIU_H1/literature/pdfs/liu_2023_arxiv_2306.08824v1.pdf",
        "sha256": "e7463fcb0d3fdfefb411289282ab344a6fee718e4f4368c1dd3cc832d0e17d5b",
        "role": "Example 5, (16), Lemma 8 and (29)",
    },
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frac_arb(value: Fraction | int) -> arb:
    value = Fraction(value)
    return arb(value.numerator) / value.denominator


def format_arb(value: arb, digits: int = 40) -> str:
    return value.str(digits)


def h_arb(u: arb) -> arb:
    if u == 0 or u == 1:
        return arb(0)
    return -(u * u.log() + (ONE - u) * (ONE - u).log())


def hp_arb(u: arb) -> arb:
    return ((ONE - u) / u).log()


def h_mp(u: mp.mpf) -> mp.mpf:
    if u == 0 or u == 1:
        return mp.mpf(0)
    return -u * mp.log(u) - (1 - u) * mp.log(1 - u)


def hp_mp(u: mp.mpf) -> mp.mpf:
    return mp.log((1 - u) / u)


def h_float(u: float) -> float:
    if u <= 0.0 or u >= 1.0:
        return 0.0
    return -u * math.log(u) - (1.0 - u) * math.log1p(-u)


def hp_float(u: float) -> float:
    return math.log((1.0 - u) / u)


def check_dependencies() -> dict[str, object]:
    out: dict[str, object] = {}
    for name, spec in DEPENDENCIES.items():
        path = (HERE / str(spec["path"])).resolve()
        observed = sha256_file(path)
        if observed != spec["sha256"]:
            raise AssertionError(
                f"dependency {name} changed: {observed} != {spec['sha256']}"
            )
        record: dict[str, object] = {
            "path": str(path.relative_to(HERE.parent)),
            "sha256": observed,
            "role": spec["role"],
            "match": True,
        }
        if "claim_status" in spec:
            payload = json.loads(path.read_text(encoding="utf-8"))
            observed_status = payload.get("claim_status")
            if observed_status != spec["claim_status"]:
                raise AssertionError(
                    f"dependency {name} status {observed_status!r} != {spec['claim_status']!r}"
                )
            record["claim_status"] = observed_status
        out[name] = record
    return out


def symbolic_identities() -> dict[str, object]:
    x, beta = sp.symbols("x beta", positive=True)
    hh, aa, bb = sp.symbols("h a b", positive=True)
    hp, ap, bp = sp.symbols("hp ap bp", real=True)
    w, p, mm = sp.symbols("w p m", positive=True)

    pi_diag = x**2 * (1 + (1 - x) ** 2)
    quartic = x**4 - 2 * x**3 + 3 * x**2 - 1
    quartic_residual = sp.expand(x**2 + pi_diag - 1 - quartic)

    d = aa + beta * (bb - aa)
    endpoint_residual = sp.factor(
        (w**2 * d - (w * x / mm) * w * hh).subs(mm, x * hh / d)
    )

    dprime = ap + beta * (bp - ap)
    rprime_numerator = (hh + x * hp) * d - x * hh * dprime
    rprime_at_equal_values = rprime_numerator.subs(bb, aa)
    beta_stationary = sp.solve(sp.Eq(rprime_at_equal_values, 0), beta)[0]
    binding_beta = ((hp + hh / x) / (hh / aa) - ap) / (bp - ap)
    binding_beta_residual = sp.factor(beta_stationary - binding_beta)
    frontier_derivative_residual = sp.factor(
        rprime_at_equal_values
        + x * hh * (bp - ap) * (beta - beta_stationary)
    )

    ff = sp.symbols("f", positive=True)
    ffp, zp = sp.symbols("fp zp", real=True)
    general_zprime_residual = sp.expand(
        zp - (2 * x + 2 * ff * ffp)
    ).subs(zp, 2 * x + 2 * ff * ffp)

    D = sp.symbols("D", positive=True)
    Dp = sp.symbols("Dp", real=True)
    G = p**2 * D - p * hh
    fixed_mean_direction = sp.diff(G, p) * (-p / x) + p**2 * Dp - p * hp
    reduced_direction = sp.factor(fixed_mean_direction.subs(p, hh / D))
    expected_direction = hh * (hh * Dp - D * (hp + hh / x)) / D**2
    fixed_mean_residual = sp.simplify(reduced_direction - expected_direction)

    ok = all(
        residual == 0
        for residual in (
            quartic_residual,
            endpoint_residual,
            binding_beta_residual,
            frontier_derivative_residual,
            general_zprime_residual,
            fixed_mean_residual,
        )
    )
    return {
        "label": "MACHINE-VERIFIED",
        "all_residuals_zero": ok,
        "quartic_identity": "x^2 + pi(x,x) - 1 = x^4 - 2x^3 + 3x^2 - 1",
        "quartic_residual": str(quartic_residual),
        "endpoint_ray_identity": "N_beta(mu_w) - (M/m)H(mu_w) = 0 when m = x h(x)/D_beta(x)",
        "endpoint_residual": str(endpoint_residual),
        "binding_beta_identity": "beta*=((h'(x)+h(x)/x)/p-a')/(b'-a') with p=h(x)/h(x^2)",
        "binding_beta_residual": str(binding_beta_residual),
        "frontier_derivative_identity": "given stationarity at beta*, the numerator of d/dx[xh/D_beta] at x* is -x h(x)(b'-a')(beta-beta*)",
        "frontier_derivative_residual": str(frontier_derivative_residual),
        "general_protocol_derivative": "d/dx[x^2+f(x)^2] = 2x+2f(x)f'(x)",
        "general_protocol_derivative_residual": str(general_zprime_residual),
        "fixed_mean_stationarity": "with p=h(x)/D, d[p^2D-ph] along d(px)=0 equals h(x)[h(x)D'-D(h'(x)+h(x)/x)]/D^2",
        "fixed_mean_residual": str(fixed_mean_residual),
    }


def certify_same_pair_frontier(parameters: ArbParameters) -> dict[str, object]:
    x = parameters.x
    hx = h_arb(x)
    hxx = h_arb(x * x)
    pi_diag = x * x * (ONE + (ONE - x) ** 2)
    hpi = h_arb(pi_diag)
    a_prime = 2 * x * hp_arb(x * x)
    pi_prime = 2 * x + 2 * x * (ONE - x) ** 2 - 2 * x * x * (ONE - x)
    b_prime = hp_arb(pi_diag) * pi_prime
    delta = b_prime - a_prime
    slope = -x * hx * delta / (hxx * hxx)
    complement_residual = pi_diag + x * x - ONE
    entropy_residual = hpi - hxx
    stationarity_residual = (
        (ONE - parameters.beta) * a_prime
        + parameters.beta * b_prime
        - (hp_arb(x) + hx / x) / parameters.p
    )

    checks = {
        "quartic_bracket_signs": bool(
            (frac_arb(parameters.root_lo) ** 4
             - 2 * frac_arb(parameters.root_lo) ** 3
             + 3 * frac_arb(parameters.root_lo) ** 2 - 1) < 0
            and
            (frac_arb(parameters.root_hi) ** 4
             - 2 * frac_arb(parameters.root_hi) ** 3
             + 3 * frac_arb(parameters.root_hi) ** 2 - 1) > 0
        ),
        "xstar_in_open_unit_interval": bool(x > 0 and x < 1),
        "beta_star_in_open_unit_interval": bool(parameters.beta > 0 and parameters.beta < 1),
        "bprime_minus_aprime_strictly_negative": bool(delta < 0),
        "endpoint_derivative_coefficient_strictly_positive": bool(slope > 0),
        "complement_residual_contains_zero": bool(complement_residual.contains(0)),
        "entropy_symmetry_residual_contains_zero": bool(entropy_residual.contains(0)),
        "stationarity_residual_contains_zero": bool(stationarity_residual.contains(0)),
    }
    if not all(checks.values()):
        raise AssertionError(f"same-pair frontier check failed: {checks}")

    return {
        "label": "PROVED",
        "theorem": (
            "Let m_full(beta) be the infimum of the positive m for which N_beta(Gamma) >= "
            "(M(mu)/m)H(mu) for every Borel probability law mu and every "
            "conditionally i.i.d. self-coupling Gamma with marginal mu. For "
            "every beta in [0,1], m_full(beta) >= m*. Equality holds at beta*, "
            "and for beta != beta* the inequality is strict. Consequently "
            "beta* is the unique maximizer of 1-m_full(beta), with maximum "
            "c'=1-m*."
        ),
        "endpoint_frontier": (
            "m_end(beta)=sup_{0<x<1} r_beta(x), "
            "r_beta(x)=x h(x)/((1-beta)h(x^2)+beta h(pi(x,x))); "
            "the one-component coupling Gamma=mu x mu gives "
            "m_full(beta)>=m_end(beta)"
        ),
        "proof": [
            "For Gamma=mu_w x mu_w and mu_w=(1-w)delta_0+w delta_x, N_beta=w^2 D_beta(x), H=w h(x), M=wx; hence every admissible m is at least r_beta(x).",
            "At the unique quartic root x*, pi(x*,x*)=1-x*^2, so h(pi)=h(x*^2) and r_beta(x*)=x*h(x*)/h(x*^2)=m* for every beta.",
            "The derivative at x* is r_beta'(x*)=[-x*h(x)*(b'-a')/h(x^2)^2]*(beta-beta*). The bracketed coefficient is strictly positive. Thus beta != beta* gives a nearby x on the increasing side with r_beta(x)>m*.",
            "At beta*, the pinned mixture theorem proves N_beta(Gamma)>=(M/m*)H(mu) for every Borel law and every conditionally i.i.d. self-coupling; the endpoint ray attains equality. Therefore m_full(beta*)=m*.",
        ],
        "premise_correction": (
            "The true endpoint tightness equation is "
            "m=xh(x)/[(1-beta)h(x^2)+beta h(pi(x,x))]. For beta>0, imposing "
            "instead m=xh(x)/h(x^2) forces h(pi(x,x))=h(x^2). Binary-entropy "
            "symmetry and strict monotonicity on either side of 1/2 reduce this "
            "to pi=x^2 (only the endpoint x=1) or pi=1-x^2 (the unique interior "
            "quartic root x*). Thus the prompt's beta-dependent x_beta cannot "
            "occur with the beta-free denominator. For the genuine endpoint "
            "frontier an interior maximizer satisfies "
            "D_beta'(x_beta)/D_beta(x_beta)=1/x_beta+h'(x_beta)/h(x_beta); "
            "away from beta* this parameterizes only m_end(beta), an obstruction, "
            "and does not prove m_full(beta)=m_end(beta). The maximizer theorem "
            "above needs no such unproved equality."
        ),
        "off_beta_exact_frontier": {
            "label": "OPEN",
            "statement": (
                "For beta!=beta*, the exact value of m_full(beta) is not "
                "identified here; the endpoint ray proves only the strict lower "
                "bound m_full(beta)>=m_end(beta)>m*. This does not affect the "
                "unique-maximizer conclusion."
            ),
        },
        "constants_exact_balls": {
            "xstar": format_arb(parameters.x, 50),
            "pstar": format_arb(parameters.p, 50),
            "mstar": format_arb(parameters.mean, 50),
            "cprime": format_arb(parameters.c, 50),
            "betastar": format_arb(parameters.beta, 50),
            "xstar_bracket": {
                "lo": str(parameters.root_lo),
                "hi": str(parameters.root_hi),
                "width": str(parameters.root_hi - parameters.root_lo),
            },
        },
        "arb_checks": checks,
        "arb_values": {
            "pi_diag_plus_xstar_squared_minus_one": format_arb(complement_residual, 12),
            "h_pi_minus_h_xstar_squared": format_arb(entropy_residual, 12),
            "bprime_minus_aprime": format_arb(delta, 30),
            "rprime_coefficient_of_beta_minus_betastar": format_arb(slope, 30),
            "stationarity_residual_at_betastar": format_arb(stationarity_residual, 12),
        },
    }


def bracket_active_root(
    polynomial: Callable[[mp.mpf], mp.mpf], guess: str
) -> tuple[Fraction, Fraction]:
    with mp.workdps(MP_DPS + 20):
        root = mp.findroot(polynomial, mp.mpf(guess))
        scale = 10**ROOT_DIGITS
        lower_numerator = int(mp.floor(root * scale))
    return (
        Fraction(lower_numerator, scale),
        Fraction(lower_numerator + 1, scale),
    )


@dataclass(frozen=True)
class Candidate:
    name: str
    formula: str
    root_lo: Fraction
    root_hi: Fraction
    x: arb
    p: arb
    m: arb
    c: arb
    beta: arb
    face_coefficient: arb
    f_arb: Callable[[arb], arb]
    f_fraction: Callable[[Fraction], Fraction]
    root_polynomial: str
    monotonicity_proof: str


def build_candidate(name: str, scale: Fraction | None = None) -> Candidate:
    if name == "scaled_four_fifths":
        if scale != Fraction(4, 5):
            raise ValueError("scaled candidate requires scale 4/5")
        ell = Fraction(4, 5)

        def f_fraction(u: Fraction) -> Fraction:
            return ell * u * (1 - u)

        def f_arb(u: arb) -> arb:
            return frac_arb(ell) * u * (ONE - u)

        def fp_arb(u: arb) -> arb:
            return frac_arb(ell) * (ONE - 2 * u)

        def polynomial_mp(u: mp.mpf) -> mp.mpf:
            ell_mp = mp.mpf(ell.numerator) / ell.denominator
            ff = ell_mp * u * (1 - u)
            return 2 * u * u + ff * ff - 1

        formula = "f(x)=(4/5)x(1-x)"
        polynomial_text = "2x^2+(16/25)x^2(1-x)^2-1"
        guess = "0.697"
        monotonicity = (
            "partial_s[st+f(s)f(t)] = t[1+(16/25)(1-t)(1-2s)] "
            ">= (9/25)t >= 0; hence the protocol argument and its diagonal "
            "are coordinatewise increasing"
        )
    elif name == "x_squared_one_minus_x":

        def f_fraction(u: Fraction) -> Fraction:
            return u * u * (1 - u)

        def f_arb(u: arb) -> arb:
            return u * u * (ONE - u)

        def fp_arb(u: arb) -> arb:
            return 2 * u - 3 * u * u

        def polynomial_mp(u: mp.mpf) -> mp.mpf:
            ff = u * u * (1 - u)
            return 2 * u * u + ff * ff - 1

        formula = "f(x)=x^2(1-x)"
        polynomial_text = "2x^2+x^4(1-x)^2-1"
        guess = "0.699"
        monotonicity = (
            "if s<=2/3 then f'(s)>=0; if s>2/3 then f'(s)>=-1 and "
            "f(t)<=t, so partial_s[st+f(s)f(t)] >= t-f(t)>0; the diagonal "
            "is therefore increasing"
        )
    else:
        raise ValueError(name)

    lo, hi = bracket_active_root(polynomial_mp, guess)
    x = frac_arb(lo).union(frac_arb(hi))
    f = f_arb(x)
    polynomial = 2 * x * x + f * f - ONE
    derivative = 4 * x + 2 * f * fp_arb(x)
    lo_value = 2 * frac_arb(lo) ** 2 + f_arb(frac_arb(lo)) ** 2 - ONE
    hi_value = 2 * frac_arb(hi) ** 2 + f_arb(frac_arb(hi)) ** 2 - ONE
    broad = frac_arb(Fraction(2, 3)).union(
        frac_arb(Fraction(7071067811865476, 10**16))
    )
    broad_derivative = 4 * broad + 2 * f_arb(broad) * fp_arb(broad)
    if not (lo_value < 0 < hi_value and derivative > 0 and broad_derivative > 0):
        raise AssertionError(f"active root not certified for {name}")

    hx = h_arb(x)
    hxx = h_arb(x * x)
    p = hx / hxx
    m = x * p
    c = ONE - m
    z = x * x + f * f
    a_prime = 2 * x * hp_arb(x * x)
    b_prime = (2 * x + 2 * f * fp_arb(x)) * hp_arb(z)
    beta = ((hp_arb(x) + hx / x) / p - a_prime) / (b_prime - a_prime)
    D = (ONE - beta) * hxx + beta * h_arb(z)
    stationarity = (ONE - beta) * a_prime + beta * b_prime \
        - (hp_arb(x) + hx / x) * D / hx
    complement = z + x * x - ONE
    if not (
        x > 0
        and x < 1
        and p > 0
        and m > 0
        and beta > 0
        and beta < 1
        and complement.contains(0)
        and stationarity.contains(0)
    ):
        raise AssertionError(f"candidate equations failed for {name}")

    return Candidate(
        name=name,
        formula=formula,
        root_lo=lo,
        root_hi=hi,
        x=x,
        p=p,
        m=m,
        c=c,
        beta=beta,
        face_coefficient=ONE - beta - ONE / (2 * m),
        f_arb=f_arb,
        f_fraction=f_fraction,
        root_polynomial=polynomial_text,
        monotonicity_proof=monotonicity,
    )


def h_range(lo: Fraction, hi: Fraction) -> arb:
    """Exact range enclosure of h on a rational interval."""
    if not Fraction(0) <= lo <= hi <= Fraction(1):
        raise AssertionError(f"entropy argument outside [0,1]: {lo}, {hi}")
    value = h_arb(frac_arb(lo)).union(h_arb(frac_arb(hi)))
    if lo <= HALF <= hi:
        value = value.union(LOG2)
    return value


def candidate_point_values(candidate: Candidate, s: Fraction, t: Fraction) -> tuple[arb, arb]:
    sa, ta = frac_arb(s), frac_arb(t)
    hs, ht = h_arb(sa), h_arb(ta)
    fs, ft = candidate.f_arb(sa), candidate.f_arb(ta)
    phi_s = (candidate.beta * h_arb(sa * sa + fs * fs)).sqrt()
    phi_t = (candidate.beta * h_arb(ta * ta + ft * ft)).sqrt()
    p2 = (ONE - candidate.beta) * h_arb(sa * ta) \
        - (ta * hs + sa * ht) / (2 * candidate.m)
    q2 = candidate.beta * h_arb(sa * ta + fs * ft)
    return p2 + phi_s * phi_t, q2 - phi_s * phi_t


def classify_interval(value: arb) -> str:
    if value.lower() >= 0:
        return "certified_nonnegative"
    if value.upper() < 0:
        return "certified_negative"
    return "unresolved"


def point_lattice(candidate: Candidate) -> dict[str, object]:
    counts = {
        "A": {"certified_nonnegative": 0, "certified_negative": 0, "unresolved": 0},
        "B": {"certified_nonnegative": 0, "certified_negative": 0, "unresolved": 0},
    }
    first_negative: dict[str, object] = {}
    minimum_lower: dict[str, arb | None] = {"A": None, "B": None}
    denominator = POINT_DENOMINATOR
    for i in range(denominator + 1):
        s = Fraction(i, denominator)
        for j in range(denominator + 1):
            t = Fraction(j, denominator)
            A, B = candidate_point_values(candidate, s, t)
            for kernel, value, structural_zero in (
                ("A", A, s == 0 or t == 0),
                ("B", B, s == t),
            ):
                if structural_zero:
                    continue
                status = classify_interval(value)
                counts[kernel][status] += 1
                lower = value.lower()
                if minimum_lower[kernel] is None or lower < minimum_lower[kernel]:
                    minimum_lower[kernel] = lower
                if status == "certified_negative" and kernel not in first_negative:
                    first_negative[kernel] = {
                        "s": str(s),
                        "t": str(t),
                        "value": format_arb(value, 30),
                    }
    return {
        "label": "COMPUTATIONAL-EVIDENCE",
        "grid": f"{{0,1/{denominator},...,1}}^2",
        "structural_points_omitted": "A: s=0 or t=0; B: s=t",
        "counts": counts,
        "minimum_lower_endpoint": {
            key: format_arb(value, 20) if value is not None else None
            for key, value in minimum_lower.items()
        },
        "first_certified_negative": first_negative or None,
        "limitation": "finite rational points only; a clean lattice is not a continuum certificate",
    }


def candidate_cell_values(
    candidate: Candidate,
    s_lo: Fraction,
    s_hi: Fraction,
    t_lo: Fraction,
    t_hi: Fraction,
) -> tuple[arb, arb]:
    fs_lo, fs_hi = candidate.f_fraction(s_lo), candidate.f_fraction(s_hi)
    ft_lo, ft_hi = candidate.f_fraction(t_lo), candidate.f_fraction(t_hi)
    protocol_lo = s_lo * t_lo + fs_lo * ft_lo
    protocol_hi = s_hi * t_hi + fs_hi * ft_hi
    diag_s_lo = s_lo * s_lo + fs_lo * fs_lo
    diag_s_hi = s_hi * s_hi + fs_hi * fs_hi
    diag_t_lo = t_lo * t_lo + ft_lo * ft_lo
    diag_t_hi = t_hi * t_hi + ft_hi * ft_hi
    if not (
        protocol_lo <= protocol_hi
        and diag_s_lo <= diag_s_hi
        and diag_t_lo <= diag_t_hi
    ):
        raise AssertionError("monotonicity used by coarse cover was violated")

    s_box = frac_arb(s_lo).union(frac_arb(s_hi))
    t_box = frac_arb(t_lo).union(frac_arb(t_hi))
    hs = h_range(s_lo, s_hi)
    ht = h_range(t_lo, t_hi)
    phi_s = (candidate.beta * h_range(diag_s_lo, diag_s_hi)).sqrt()
    phi_t = (candidate.beta * h_range(diag_t_lo, diag_t_hi)).sqrt()
    p2 = (ONE - candidate.beta) * h_range(s_lo * t_lo, s_hi * t_hi) \
        - (t_box * hs + s_box * ht) / (2 * candidate.m)
    q2 = candidate.beta * h_range(protocol_lo, protocol_hi)
    return p2 + phi_s * phi_t, q2 - phi_s * phi_t


def coarse_partition(candidate: Candidate) -> dict[str, object]:
    counts = {
        "A": {"certified_nonnegative": 0, "certified_negative": 0, "unresolved": 0},
        "B": {"certified_nonnegative": 0, "certified_negative": 0, "unresolved": 0},
    }
    first_negative: dict[str, object] = {}
    denominator = COVER_DENOMINATOR
    for i in range(denominator):
        s_lo, s_hi = Fraction(i, denominator), Fraction(i + 1, denominator)
        for j in range(denominator):
            t_lo, t_hi = Fraction(j, denominator), Fraction(j + 1, denominator)
            A, B = candidate_cell_values(candidate, s_lo, s_hi, t_lo, t_hi)
            for kernel, value in (("A", A), ("B", B)):
                status = classify_interval(value)
                counts[kernel][status] += 1
                if status == "certified_negative" and kernel not in first_negative:
                    first_negative[kernel] = {
                        "cell": [str(s_lo), str(s_hi), str(t_lo), str(t_hi)],
                        "value": format_arb(value, 20),
                    }
    return {
        "label": "COMPUTATIONAL-EVIDENCE",
        "partition": f"{denominator}x{denominator} exact rational rectangles",
        "precision_bits": ctx.prec,
        "monotonicity_used_for_protocol_ranges": candidate.monotonicity_proof,
        "counts": counts,
        "first_certified_negative_cell": first_negative or None,
        "limitation": (
            "Natural interval enclosures leave unresolved cells near zero strata and through "
            "dependency overestimation. Only a certified_negative cell would refute a kernel; "
            "zero such cells does not prove nonnegativity."
        ),
    }


def candidate_record(candidate: Candidate, baseline: ArbParameters) -> dict[str, object]:
    improvement = candidate.c - baseline.c
    face_witness = candidate.face_coefficient * LOG2
    improvement_certified = bool(improvement > 0)
    if not improvement_certified:
        raise AssertionError(f"candidate {candidate.name} does not exceed c'")
    face_negative = bool(candidate.face_coefficient < 0)
    point = point_lattice(candidate)
    cover = coarse_partition(candidate)
    observed_counts = {
        "point_lattice": point["counts"],
        "coarse_partition": cover["counts"],
    }
    expected_counts = EXPECTED_DIAGNOSTIC_COUNTS[candidate.name]
    if observed_counts != expected_counts:
        raise AssertionError(
            f"diagnostic count regression for {candidate.name}: "
            f"{observed_counts} != {expected_counts}"
        )
    return {
        "label": "COMPUTATIONAL-EVIDENCE",
        "name": candidate.name,
        "f": candidate.formula,
        "active_two_point_law": (
            f"(1-p)delta_0+p delta_x with x in {format_arb(candidate.x, 35)}, "
            f"p in {format_arb(candidate.p, 35)}"
        ),
        "active_equations": {
            "label": "MACHINE-VERIFIED",
            "root_polynomial": candidate.root_polynomial,
            "root_bracket": [str(candidate.root_lo), str(candidate.root_hi)],
            "root_bracket_width": str(candidate.root_hi - candidate.root_lo),
            "unique_root_on_2/3_to_1/sqrt(2)": True,
            "tightness": "p=h(x)/h(x^2), m=px",
            "stationarity_beta": format_arb(candidate.beta, 40),
        },
        "m": format_arb(candidate.m, 40),
        "one_minus_m": format_arb(candidate.c, 40),
        "one_minus_m_minus_cprime": format_arb(improvement, 30),
        "candidate_value_exceeds_cprime_certified": True,
        "candidate_value_scope": (
            "The strict scalar comparison concerns only the active endpoint "
            "equations; it is not a certified union-closed bound."
        ),
        "A_face_check": {
            "label": "MACHINE-VERIFIED",
            "identity": "A(s,1)=[1-beta-1/(2m)]h(s), because f(1)=h(1)=0",
            "coefficient": format_arb(candidate.face_coefficient, 30),
            "A_half_one": format_arb(face_witness, 30),
            "conclusion": (
                "A(1/2,1)<0 is a certified negative witness"
                if face_negative
                else "the exact t=1 necessary condition is satisfied"
            ),
        },
        "point_lattice_400_bit": point,
        "coarse_interval_partition_400_bit": cover,
        "diagnostic_counts_regression_asserted": True,
        "verdict": (
            "The A>=0 candidate route has a certified negative witness; "
            "this does not refute the underlying protocol inequality"
            if face_negative
            else "No negative witness appears in the requested coarse A/B diagnostics; "
            "continuum A>=0 and B>=0 remain OPEN"
        ),
        "remaining_continuum_claim": {
            "label": "OPEN",
            "statement": (
                "B>=0 on [0,1]^2 remains unproved; A>=0 is already false, so "
                "this A/B sufficient route cannot certify the candidate"
                if face_negative
                else "A>=0 and B>=0 on [0,1]^2 both remain unproved"
            ),
        },
        "no_certificate_claimed": True,
    }


def active_root_float(f: Callable[[float], float]) -> float:
    lo = 2.0 / 3.0
    hi = 1.0 / math.sqrt(2.0)

    def equation(u: float) -> float:
        return 2.0 * u * u + f(u) ** 2 - 1.0

    flo, fhi = equation(lo), equation(hi)
    if flo > 1e-12 or fhi < -1e-12:
        raise ValueError("no bracket on the Liu complementary branch")
    if abs(flo) <= 1e-15:
        return lo
    for _ in range(100):
        mid = (lo + hi) / 2.0
        if equation(mid) <= 0.0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def active_candidate_float(
    name: str,
    formula: str,
    f: Callable[[float], float],
    fp: Callable[[float], float],
) -> dict[str, object]:
    x = active_root_float(f)
    z = x * x + f(x) ** 2
    hxx = h_float(x * x)
    px = h_float(x) / hxx
    m = x * px
    a_prime = 2.0 * x * hp_float(x * x)
    b_prime = (2.0 * x + 2.0 * f(x) * fp(x)) * hp_float(z)
    beta = ((hp_float(x) + h_float(x) / x) / px - a_prime) / (b_prime - a_prime)
    face = 1.0 - beta - 1.0 / (2.0 * m)
    return {
        "name": name,
        "f": formula,
        "x": format(x, ".17g"),
        "p": format(px, ".17g"),
        "m": format(m, ".17g"),
        "one_minus_m": format(1.0 - m, ".17g"),
        "beta": format(beta, ".17g"),
        "beta_in_open_unit_interval": 0.0 < beta < 1.0,
        "A_face_coefficient": format(face, ".17g"),
        "A_face_necessary_condition_survives": face >= 0.0,
    }


def discovery_scan() -> dict[str, object]:
    baseline = active_candidate_float(
        "baseline_k1",
        "x(1-x)",
        lambda x: x * (1.0 - x),
        lambda x: 1.0 - 2.0 * x,
    )
    baseline_c = float(str(baseline["one_minus_m"]))

    powers = []
    for k in range(1, 9):
        powers.append(
            active_candidate_float(
                f"power_k{k}",
                f"x(1-x)^{k}",
                lambda x, k=k: x * (1.0 - x) ** k,
                lambda x, k=k: (1.0 - x) ** k
                - k * x * (1.0 - x) ** (k - 1),
            )
        )

    squared = active_candidate_float(
        "x_squared_one_minus_x",
        "x^2(1-x)",
        lambda x: x * x * (1.0 - x),
        lambda x: 2.0 * x - 3.0 * x * x,
    )

    mixtures: list[dict[str, object]] = []
    denominator = 40
    for ia in range(denominator + 1):
        for ib in range(denominator + 1 - ia):
            if ia == 0 and ib == 0:
                continue
            a, b = ia / denominator, ib / denominator
            record = active_candidate_float(
                f"mix_{ia}_{ib}",
                f"({ia}/{denominator})x(1-x)+({ib}/{denominator})x(1-x)^2",
                lambda x, a=a, b=b: x * (1.0 - x) * (a + b * (1.0 - x)),
                lambda x, a=a, b=b: (1.0 - 2.0 * x) * (a + b * (1.0 - x))
                - b * x * (1.0 - x),
            )
            record["a"] = f"{ia}/{denominator}"
            record["b"] = f"{ib}/{denominator}"
            mixtures.append(record)

    admissible_improvements = [
        r
        for r in mixtures
        if bool(r["beta_in_open_unit_interval"])
        and float(str(r["one_minus_m"])) > baseline_c
    ]
    face_survivors = [
        r for r in admissible_improvements if bool(r["A_face_necessary_condition_survives"])
    ]
    sort_key = lambda r: (-float(str(r["one_minus_m"])), str(r["a"]), str(r["b"]))
    admissible_improvements.sort(key=sort_key)
    face_survivors.sort(key=sort_key)

    return {
        "label": "COMPUTATIONAL-EVIDENCE",
        "arithmetic": "IEEE binary64 bisection and evaluation; discovery only",
        "general_f_constraint": "0<=f(x)<=min(x,1-x) on [0,1]",
        "joint_active_branch": "2x^2+f(x)^2=1, p=h(x)/h(x^2), m=px, beta from fixed-mean stationarity",
        "baseline": baseline,
        "powers_k_1_through_8": powers,
        "power_result": (
            "k=1 is Liu's baseline; k=2,...,8 have beta>1 on this interior joint-stationary branch and hence are not admissible two-protocol optima"
        ),
        "x_squared_one_minus_x": squared,
        "mixture_domain": (
            "a,b in {0,1/40,...,1}, a>=0, b>=0, a+b<=1; this simplex is a sufficient protocol-validity domain"
        ),
        "mixtures_scanned": len(mixtures),
        "admissible_improvements": len(admissible_improvements),
        "admissible_improvements_surviving_A_face_test": len(face_survivors),
        "top_five_admissible": admissible_improvements[:5],
        "top_five_surviving_A_face_test": face_survivors[:5],
        "promotion_policy": (
            "The scan rows are unpromoted endpoint-law hits. The separate "
            "general_f_candidates records promote f=(4/5)x(1-x), chosen for "
            "one rational parameter and a robust positive A-face margin, and "
            "f=x^2(1-x), chosen as the requested named-family negative control; "
            "both receive the 400-bit diagnostics."
        ),
        "interpretation": (
            "Values 1-m above c' in this binary64 scan are discovery hits, not "
            "bounds or certificate claims. Only promoted candidates below receive "
            "400-bit diagnostics, which still decide no continuum kernel globally "
            "unless they exhibit a certified negative witness."
        ),
    }


def general_f_derivation() -> dict[str, object]:
    return {
        "label": "PROVED",
        "protocol": (
            "In Liu Example 5 use x=1-s. A measurable f with "
            "0<=f(x)<=min(x,1-x) gives K_f(x,y)=h(xy+f(x)f(y)); "
            "the stationarity equations below additionally assume that f is "
            "differentiable at the active point x."
        ),
        "endpoint_law": "mu_w=(1-w)delta_0+w delta_x",
        "definitions": [
            "z_f(x)=x^2+f(x)^2",
            "D_{beta,f}(x)=(1-beta)h(x^2)+beta h(z_f(x))",
            "N=w^2 D_{beta,f}(x), H=w h(x), M=wx",
        ],
        "tightness_equations": [
            "p=h(x)/D_{beta,f}(x)",
            "m=px=x h(x)/D_{beta,f}(x)",
        ],
        "range_on_nontrivial_branch": (
            "From f(x)<=x and f(x)^2=1-2x^2 one gets "
            "x>=1/sqrt(3)>1/2. Then f(x)<=1-x gives x>=2/3, while "
            "f(x)>0 gives x<1/sqrt(2)."
        ),
        "fixed_mean_stationarity": (
            "D'_{beta,f}(x)/D_{beta,f}(x)=1/x+h'(x)/h(x), "
            "where D'=(1-beta)2x h'(x^2)+beta(2x+2f f')h'(x^2+f^2)"
        ),
        "beta_formula": (
            "beta=(C-a')/(b'-a'), a'=2x h'(x^2), "
            "b'=(2x+2f(x)f'(x))h'(x^2+f(x)^2), "
            "C=(h'(x)+h(x)/x)/p"
        ),
        "joint_interior_protocol_weight_stationarity": (
            "If the beta optimum is interior and this single endpoint law is the active branch, "
            "envelope stationarity gives h(x^2+f(x)^2)=h(x^2). For f(x)>0 on the "
            "nontrivial branch this is x^2+f(x)^2=1-x^2, i.e. 2x^2+f(x)^2=1; "
            "then p=h(x)/h(x^2) and m=xh(x)/h(x^2)."
        ),
        "scope_warning": (
            "These equations identify Liu's analogous active two-point branch. They do not prove "
            "that the branch is the global minimizer, nor that A>=0 and B>=0 on the continuum."
        ),
    }


def sawin_open_obligation() -> dict[str, object]:
    return {
        "label": "OPEN",
        "source": (
            "Liu, arXiv:2306.08824v1, Definition 1 and equation (16) on page 4; "
            "Example 5 and (23) on page 6; Lemma 8 and (29) on page 7"
        ),
        "classes": {
            "C1(mu)": "the singleton {mu x mu} (Gilmer's independent protocol)",
            "C2(mu)": (
                "for the conservative certificate obligation used here, all symmetric "
                "self-couplings of mu; this contains every coupling induced by Sawin's "
                "max-entropy protocol, and Liu states that probably no better safe "
                "choice is known"
            ),
            "C3(mu)": (
                "{couplings of mu and mu} intersected with the weak closure of the convex hull "
                "of symmetric rank-one measures (the conditionally-i.i.d. class)"
            ),
        },
        "kernels_in_Liu_variables": {
            "g1(s,t)": "h((1-s)(1-t))",
            "g2(s,r)": "h(max(s,r,min(s+r,1/2)))",
            "g3_example4(s,r)": (
                "h((1-s)(1-r)+a(s)a(r)(min(1-s,1-r)-(1-s)(1-r)))"
            ),
            "g3_example5(s,r)": "h((1-s)(1-r)+f(1-s)f(1-r))",
        },
        "exact_obligation": (
            "To raise the current c', one must exhibit c_new>c', C>1, and fixed "
            "weights w1,w2,w3>=0 with w1+w2+w3=1 such that, for every probability "
            "law mu on [0,1] with E_mu[S]<=c_new,\n"
            "  w1 inf_{P_ST in C1(mu)} E[g1(S,T)]\n"
            "+ w2 inf_{P_SR1 in C2(mu)} E[g2(S,R1)]\n"
            "+ w3 inf_{P_SR2 in C3(mu)} E[g3(S,R2)] >= C E[h(S)].\n"
            "Liu's Lemma 8 uses w1=(1-alpha*)(1-beta), "
            "w2=alpha*(1-beta), w3=beta and g3 from Example 4."
        ),
        "complemented_form_for_current_pair": (
            "With X=1-S, the first and Example-5 terms become h(XY) and "
            "h(XY+f(X)f(Y)), exactly the two kernels certified in the H2 chain. "
            "The Sawin kernel remains minimized over C2(mu), not over product or rank-one mixtures."
        ),
        "literature_scope": (
            "Liu's Lemma 8 proves a non-explicit existence result strictly above "
            "the older c*=0.3823455, relative to the Sawin two-protocol input. "
            "It neither gives explicit (beta,c_new,C) nor proves a value above "
            "the current certified c'."
        ),
        "obstruction": (
            "The existing reduction controls the product term and C3 through conditionally-i.i.d. "
            "rank-one mixtures. Sawin's protocol is not conditionally i.i.d. (Liu notes that its "
            "2x2 kernel matrix is non-PSD for s=t<1/4), and C2(mu) must conservatively include all "
            "symmetric self-couplings. No finite-support/extreme-point reduction or uniform lower "
            "bound for inf_{C2(mu)} E[g2] at means up to c_new is proved in this repository. "
            "That infinite-dimensional symmetric-coupling minimization is the named obstruction."
        ),
        "attempted_here": False,
    }


def run_mutations(
    parameters: ArbParameters,
    robust: Candidate,
    squared: Candidate,
) -> list[dict[str, object]]:
    mutations: list[dict[str, object]] = []

    x = sp.symbols("x")
    wrong_pi = x**2 * (1 + (1 - x))
    quartic = x**4 - 2 * x**3 + 3 * x**2 - 1
    wrong_residual = sp.expand(x**2 + wrong_pi - 1 - quartic)
    mutations.append(
        {
            "mutation": "replace (1-x)^2 by (1-x) in pi(x,x)",
            "caught": wrong_residual != 0,
            "witness": str(wrong_residual),
        }
    )

    bx = parameters.x
    hxx = h_arb(bx * bx)
    pi_diag = bx * bx * (ONE + (ONE - bx) ** 2)
    a_prime = 2 * bx * hp_arb(bx * bx)
    pi_prime = 2 * bx + 2 * bx * (ONE - bx) ** 2 - 2 * bx * bx * (ONE - bx)
    delta = hp_arb(pi_diag) * pi_prime - a_prime
    correct_slope = -bx * h_arb(bx) * delta / (hxx * hxx)
    wrong_slope = -correct_slope
    mutations.append(
        {
            "mutation": "reverse the endpoint-frontier derivative direction",
            "caught": bool(wrong_slope < 0 < correct_slope),
            "witness": format_arb(wrong_slope, 20),
        }
    )

    cx = robust.x
    cf = robust.f_arb(cx)
    cz = cx * cx + cf * cf
    wrong_bprime = 2 * cx * hp_arb(cz)  # deliberately omits 2 f f'
    a_prime = 2 * cx * hp_arb(cx * cx)
    hx = h_arb(cx)
    wrong_stationarity = (
        (ONE - robust.beta) * a_prime
        + robust.beta * wrong_bprime
        - (hp_arb(cx) + hx / cx) / robust.p
    )
    mutations.append(
        {
            "mutation": "drop 2 f(x) f'(x) from the general protocol derivative",
            "caught": not wrong_stationarity.contains(0),
            "witness": format_arb(wrong_stationarity, 20),
        }
    )

    wrong_face = ONE - ONE / (2 * squared.m)  # deliberately drops -beta
    true_face = squared.face_coefficient
    mutations.append(
        {
            "mutation": "drop -beta from the A(s,1) face coefficient",
            "caught": bool(true_face < 0 < wrong_face),
            "witness": {
                "true": format_arb(true_face, 20),
                "mutated": format_arb(wrong_face, 20),
            },
        }
    )

    if not all(bool(item["caught"]) for item in mutations):
        raise AssertionError(f"mutation escaped: {mutations}")
    return mutations


def validate_labels(value: object) -> None:
    if isinstance(value, dict):
        for key, child in value.items():
            if key in {"label", "claim_status"} and child not in ALLOWED_LABELS:
                raise AssertionError(f"unsupported label: {child!r}")
            validate_labels(child)
    elif isinstance(value, list):
        for child in value:
            validate_labels(child)


def canonical_bytes(report: dict[str, object]) -> bytes:
    body = dict(report)
    body.pop("report_sha256", None)
    return (json.dumps(body, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")


def build_report() -> dict[str, object]:
    ctx.prec = max(ctx.prec, PRECISION_BITS)
    mp_parameters = solve_equation_parameters(100)
    parameters = certify_equation_parameters(mp_parameters)
    robust = build_candidate("scaled_four_fifths", Fraction(4, 5))
    squared = build_candidate("x_squared_one_minus_x")

    identities = symbolic_identities()
    if not identities["all_residuals_zero"]:
        raise AssertionError("a symbolic identity failed")

    report: dict[str, object] = {
        "artifact": "liu9-cprime-frontier",
        "claim_status": "PROVED",
        "claim": (
            "For Liu's same two protocols with beta varied over [0,1], beta* is "
            "the unique maximizer of the obtainable full-measure constant: the "
            "maximum is c'=1-m*. General-f values above c' below are discovery "
            "candidates only; the Sawin three-protocol obligation remains OPEN."
        ),
        "precision_bits": ctx.prec,
        "natural_log": True,
        "dependencies": check_dependencies(),
        "symbolic_identities": identities,
        "same_protocol_beta_frontier": certify_same_pair_frontier(parameters),
        "general_f_derivation": general_f_derivation(),
        "general_f_discovery_scan": discovery_scan(),
        "general_f_candidates": {
            "scaled_four_fifths": candidate_record(robust, parameters),
            "x_squared_one_minus_x": candidate_record(squared, parameters),
        },
        "sawin_three_protocol": sawin_open_obligation(),
        "mutations": run_mutations(parameters, robust, squared),
        "all_mutations_caught": True,
        "tool_sha256": sha256_file(Path(__file__)),
        "label_legend": sorted(ALLOWED_LABELS),
    }
    validate_labels(report)
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
    args.output.write_text(
        json.dumps(report, indent=1, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("CPRIME_FRONTIER", report["claim_status"])
    print("same_pair beta* unique maximizer; max c'", report["same_protocol_beta_frontier"]["constants_exact_balls"]["cprime"])
    for name, candidate in report["general_f_candidates"].items():
        print(name, candidate["verdict"], "one_minus_m", candidate["one_minus_m"])
    print("SAWIN_THREE_PROTOCOL", report["sawin_three_protocol"]["label"])
    print("MUTATIONS", len(report["mutations"]), "all caught", report["all_mutations_caught"])
    print("REPORT", args.output)
    print("REPORT_SHA256", report["report_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
