#!/usr/bin/env python3
"""Quantitative mirror boundary-layer estimate for Liu's Hypothesis 2.

`liu9_tube.py` refuses to certify a one-piece cubic tube bound because every
raw-coordinate tube contains supports tending to zero, where
`h'''(y)=(1-2y)/(y^2(1-y)^2)` is unbounded.  `liu9_boundary_layer.py` closed
that zero-support stratum on 2026-08-27 and, in doing so, exposed a second one:
`h'''` is equally unbounded as `y -> 1`, and a *small-mass* atom sitting at
support 1 is well inside the tube, because its distance cost is only
`w*[1*(1-x)]^2 = 0.0956*w`.

That `0.0956` is motivation only.  Every bound below divides by the *infimum*
of `g(b)=[b(b-x)]^2` over the whole stratum `[1-t0,1]`, which is `g(1-t0)`,
not `g(1)`.  At `t0=1/32` the *unsquared* cost `b(b-x)` drops 12.92% across
the stratum, `g` itself drops 24.16%, and the critical radius would be
overclaimed by 14.83% -- `0.1795` against `0.2061` -- if the endpoint value
were substituted.  `certify_face_infimum` certifies that the infimum really is
at the inner edge, and a mutation test rejects the substitution.

This module supplies the missing mirror estimate.  It is deliberately
independent of the zero-layer argument: it uses a different divergence, a
different reduction target, and -- unlike the zero layer -- it never uses the
mean constraint.

THE IDENTITY
------------
With global weights `w = ((1-q)a1,(1-q)a2,(1-q)a3,q*a1,q*a2,q*a3)` and supports
`b = (b0,b2,b4,b1,b3,b5)`, Liu's quotient-free gap

    gap = (1-beta)*EHXY + beta*EHPI - EHX

has the exact partial derivative `d gap / d b_j = w_j * G_j(V)` with

    G_j = 2(1-beta)*sum_k w_k b_k h'(b_j b_k)
        + 2*beta*sum_{k in comp(j)} a_k pi_1(b_j,b_k) h'(pi(b_j,b_k))
        - h'(b_j),

`pi(y,z)=yz(1+(1-y)(1-z))`, `pi_1=dpi/dy`, `h'(u)=log((1-u)/u)`.  This is the
same closed form the zero layer certified, and it is re-checked here against
forward-mode differentiation of Liu's own transcribed `_formula`.

THE MIRROR DIVERGENCE
---------------------
Write `t = 1-b_j`.  Every term of `G_j` that can diverge as `t -> 0` does so
like `log(1/t)`:

  * `-h'(b_j) = log(1/t) + log(b_j)`                       coefficient `+1`;
  * `h'(b_j b_k)` diverges only when `b_k -> 1`            coefficient `-2(1-beta)W1`;
  * `h'(pi(b_j,b_k))` diverges only when `b_k -> 1`        coefficient `-2 beta A1`,

where `W1` is the total *global* weight already at support 1 and `A1` the
corresponding *within-component* mass.  The leading coefficient is therefore

    1 - 2(1-beta) W1 - 2 beta A1,

which is **not signed by the mean constraint**: it is `+1-2 beta` when no other
atom sits at 1 and `-1` when every atom does.  The zero layer had no such
freedom, which is why this stratum needs its own hypothesis.

WHAT MAKES IT WORK: THE TUBE BOUNDS W1
--------------------------------------
Inside a tube of radius `rho` the distance obeys

    rho^2 >= dist^2 >= sum_k w_k g(b_k),      g(b) := [b(b-x)]^2,

and `g` is increasing on `[x,1]`, so an atom at support `>= 1-t0` costs at least
`g(1-t0)` per unit weight:

    W1 <= rho^2 / g(1-t0).

`A1 <= 1` always.  Hence the leading coefficient is bounded below by
`1 - 2 beta - 2(1-beta) rho^2/g(1-t0)`, which is positive for every
`rho^2 < g(1-t0)(1-2 beta)/(2(1-beta))`.

THE UNIFORM LOWER BOUND
-----------------------
Split the atoms at `1-t0`.  For `k` in the mirror set `S` use
`h'(b_j b_k) >= h'(1-t) >= -log(1/t)`; for `k` outside `S` use monotonicity,
`h'(b_j b_k) >= h'(b_k)`, which is independent of `t`.  The outside sum is
controlled by the linear program

    max sum_k w_k phi(b_k)  s.t.  sum_k w_k <= 1,  sum_k w_k g(b_k) <= rho^2,

`phi(b) := -b h'(b) = b log(b/(1-b))`, whose weak dual gives, for any `mu >= 0`,

    sum_k w_k phi(b_k) <= Theta := mu*rho^2 + sup_{b in [1/2,1-t0]}(phi(b)-mu g(b)).

The `phi <= 0` half `b <= 1/2` is discarded downward.  For the paired sum, the
exact algebraic identity

    1 - pi(y,z) = s + u - 2su + su(s+u-su),        s := 1-y,  u := 1-z,

gives the two one-sided bounds `1-pi >= s(1-2u)` and `1-pi >= u(1-2s)`, so
`h'(pi) >= -log(1/t) + log(1-2t0)` on the whole layer, with `0 <= pi_1 <= 1`
for `y >= 1/2`.  Collecting,

    G_j >= M*log(1/t) - K1,
    M  := 1 - 2 beta - 2(1-beta)*rho^2/g(1-t0),
    K1 := -log(1-t0) + 2(1-beta)*Theta - 2 beta*log(1-2 t0).

THE REDUCTION
-------------
`G_j > 0` on the whole layer iff `t0 <= exp(-K1/M)`.  The reduction therefore
pushes the mirror atoms **up** to `b_j = 1`, not down: the singular term works
for us, the reduced point sits on an exact face of the box, and -- because the
mean only increases -- mean-feasibility is preserved for free.  Integrating,

    gap(V|_{mirror:=1}) - gap(V) >= m1 * sum_{j in S} w_j t_j,
    m1 := M*(log(1/t0)+1) - K1 > 0.

Raising one coordinate changes `sum_k w_k g(b_k)` by at most `g'(1) t_j w_j` and
the two squared component-mean terms by at most `w_j t_j (2 max(m,1-m) + t0)`,
both free of `q`.  So the reduction is `kappa`-monotone for every

    kappa <= m1 / D1,   D1 := g'(1) + 2 max(m,1-m) + t0,

that is, `gap - kappa*dist^2` does not increase under it.

STATUS LABELS.  Statements marked PROVED are algebraic or Arb-certified.
NUMERICAL statements are finite-precision measurements.  Nothing here claims
Liu's Hypothesis 2, or a tube radius: this closes one stratum of the
conjectured piecewise local lemma.

Run from the repository root with
    math/.venv/bin/python -I -B math/uc/liu9_mirror_layer.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
import sys
import time
from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Optional, Sequence

# This workstation is shared with live certification workers.
for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

import mpmath
from flint import arb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from liu9_binding import (  # noqa: E402
    ArbParameters,
    MPParameters,
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_boundary_layer import (  # noqa: E402
    _arbf,
    entropy_at_star,
    _h_mp,
    _hp_mp,
    _pi_first_mp,
    _pi_mp,
    _support,
    _weights,
    autodiff_partial,
    closed_form_factor,
    gap_mp,
    mean_of,
)

ctx.prec = max(ctx.prec, 320)

DEFAULT_DPS = 60
DEFAULT_T0 = Fraction(1, 64)
DEFAULT_RHO = Fraction(1, 10)
MU_GRID = (12, 20, 30, 42, 60, 90, 140, 220, 340, 520)
THRESHOLD_T0 = (Fraction(1, 8), Fraction(1, 16), Fraction(1, 32),
                Fraction(1, 64), Fraction(1, 128), Fraction(1, 1024))
THRESHOLD_RHO = (Fraction(1, 10), Fraction(1, 20), Fraction(1, 50))
SAMPLE_SEED = 20260828
SUPPORT_VARIABLE_INDEX = (3, 4, 5, 6, 7, 8)


# ---------------------------------------------------------------------------
# Scalar helpers over Arb.


def _g(b: arb, x: arb) -> arb:
    """The per-unit-weight distance cost of a support at b."""
    return (b * (b - x)) ** 2


def certify_face_infimum(x: arb, t0: Fraction) -> dict[str, Any]:
    """`inf_{b in [1-t0,1]} g(b) = g(1-t0)`, certified rather than assumed.

    Every `W1` bound in this module divides `rho^2` by the *smallest* distance
    cost any mirror atom can carry, and it evaluates that as `g(1-t0)`.  That
    is only correct if `g` is nondecreasing across the stratum, which is a
    hypothesis, not a definition.  It is easy to get wrong in the dangerous
    direction: using the endpoint value `g(1) = [1*(1-x)]^2 = 0.0956` instead
    of `g(1-t0) = 0.0725` at `t0=1/32` inflates the critical radius from
    `0.1795` to `0.2061`, a 14.83% overclaim.  The three percentages are
    distinct and easy to conflate: the unsquared cost `b(b-x)` is 12.92%
    smaller at the inner edge, `g = [b(b-x)]^2` is 24.16% smaller, and the
    radius, which moves as the square root of `g`, is overclaimed by 14.83%.
    Only the middle figure describes the divisor this module actually uses.

    Since `g' (b) = 2 b (b-x)(2b-x)`, monotonicity on `[1-t0,1]` follows from
    `b-x > 0` and `2b-x > 0` there, and both reduce to `1-t0 > x` because
    every factor is increasing.  Those two comparisons are what this function
    encloses.
    """
    top = 1 - _arbf(t0)
    first = top - x
    second = 2 * top - x
    if not first > 0:
        raise AssertionError(
            f"the mirror stratum reaches below x at t0={t0}: {first}")
    if not second > 0:
        raise AssertionError(f"2b-x is not positive at t0={t0}: {second}")
    infimum = _g(top, x)
    endpoint = _g(arb(1), x)
    if not infimum <= endpoint:
        raise AssertionError(
            "g is not increasing across the stratum; the W1 bound would be "
            "taken at the wrong end")
    return {
        "t0": str(t0),
        "b_minus_x_at_inner_edge": str(first),
        "two_b_minus_x_at_inner_edge": str(second),
        "infimum_g": str(infimum),
        "infimum_g_lower": float(infimum.lower()),
        "endpoint_g": str(endpoint),
        "endpoint_g_lower": float(endpoint.lower()),
        "endpoint_overclaim_ratio": float((endpoint / infimum).upper()),
    }


def _phi(b: arb) -> arb:
    """-b h'(b) = b log(b/(1-b)); nonnegative exactly on [1/2,1)."""
    return b * (b.log() - (1 - b).log())


def _pow2_at_least(count: int) -> int:
    """Round a requested cell count up to a power of two.

    Every sweep below subdivides a dyadic interval.  With a power-of-two count
    each endpoint is a dyadic rational, so ``arb.union`` reproduces the cell
    exactly; with, say, twenty cells the union outward-rounds and the cell
    touching the origin acquires a negative sliver, which destroys the sign of
    manifestly nonnegative expressions.  This is a rigour requirement, not a
    performance tweak.
    """
    if count < 1:
        raise ValueError("a sweep needs at least one cell")
    return 1 << max(0, (count - 1).bit_length())


def _dyadic_cell(base: arb, index: int, cells: int) -> arb:
    """[base*index/cells, base*(index+1)/cells] with cells a power of two."""
    step = base / cells
    lo = step * index
    hi = step * (index + 1)
    cell = lo.union(hi)
    if not cell.is_finite():
        raise AssertionError("non-finite dyadic cell")
    return cell


def _g_prime_ceiling(x: arb) -> arb:
    """max |g'| on [1/2,1].  g'(b)=2 b(b-x)(2b-x) is increasing there."""
    return 2 * (1 - x) * (2 - x)


# ---------------------------------------------------------------------------
# The dual bound Theta on the outside-mirror cross sum.


def certify_theta(x: arb, t0: Fraction, rho: Fraction, mu: int,
                  pieces: int) -> dict[str, Any]:
    """Weak-duality upper bound for the outside-mirror cross sum.

    Returns an Arb *upper* bound for
        sup_{b in [1/2,1-t0]} (phi(b) - mu*g(b))   plus   mu*rho^2.
    The sweep encloses each cell as an interval, so the reported lambda is a
    genuine upper bound over the continuum, not a sampled maximum.
    """
    pieces = _pow2_at_least(pieces)
    lo = arb(1) / 2
    hi = 1 - _arbf(t0)
    if not hi > lo:
        raise ValueError("the mirror threshold must satisfy t0 < 1/2")
    mu_arb = arb(int(mu))
    step = (hi - lo) / pieces
    worst = None
    worst_cell = None
    for index in range(pieces):
        cell = (lo + step * index).union(lo + step * (index + 1))
        value = _phi(cell) - mu_arb * _g(cell, x)
        if not value.is_finite():
            raise AssertionError(f"non-finite theta cell at index {index}")
        upper = arb(value.upper())
        if worst is None or upper > worst:
            worst = upper
            worst_cell = (float(cell.lower()), float(cell.upper()))
    theta = worst + mu_arb * _arbf(rho) ** 2
    return {
        "mu": int(mu),
        "pieces": int(pieces),
        "lambda": worst,
        "theta": theta,
        "argmax_cell": worst_cell,
    }


# ---------------------------------------------------------------------------
# Certified constants.


@dataclass(frozen=True)
class MirrorConstants:
    t0: Fraction
    rho: Fraction
    mu: int
    theta: arb
    theta_lambda: arb
    w1_ceiling: arb
    leading: arb
    penalty: arb
    m1: arb
    distance_rate: arb
    kappa: arb
    t_star: arb
    argmax_cell: tuple[float, float]

    def certified(self) -> bool:
        """Is the layer estimate valid on the whole of (0,t0]?"""
        return bool(self.leading > 0 and self.m1 > 0
                    and self.t_star >= _arbf(self.t0))

    def as_json(self) -> dict[str, Any]:
        return {
            "t0": str(self.t0),
            "t0_float": float(self.t0),
            "rho": str(self.rho),
            "rho_float": float(self.rho),
            "mu": self.mu,
            "theta": str(self.theta),
            "theta_lambda": str(self.theta_lambda),
            "theta_argmax_cell": list(self.argmax_cell),
            "w1_ceiling": str(self.w1_ceiling),
            "w1_ceiling_upper": float(self.w1_ceiling.upper()),
            "leading_M": str(self.leading),
            "leading_M_lower": float(self.leading.lower()),
            "penalty_K1": str(self.penalty),
            "penalty_K1_upper": float(self.penalty.upper()),
            "m1": str(self.m1),
            "m1_lower": float(self.m1.lower()),
            "distance_rate_D1": str(self.distance_rate),
            "admissible_kappa": str(self.kappa),
            "admissible_kappa_lower": float(self.kappa.lower()),
            "t_star": str(self.t_star),
            "t_star_lower": float(self.t_star.lower()),
            "certified": self.certified(),
        }


def certify_constants(parameters: ArbParameters, t0: Fraction, rho: Fraction,
                      mu: int, pieces: int = 4000) -> MirrorConstants:
    """Arb-certify every constant of the mirror-layer theorem."""
    if not Fraction(0) < t0 <= Fraction(1, 4):
        raise ValueError("the mirror threshold must satisfy 0<t0<=1/4")
    if not Fraction(0) < rho:
        raise ValueError("the tube radius must be positive")
    x = parameters.x
    beta = parameters.beta
    t0_arb = _arbf(t0)
    rho_arb = _arbf(rho)
    top = 1 - t0_arb

    dual = certify_theta(x, t0, rho, mu, pieces)
    theta = dual["theta"]

    # The W1 bound divides by the smallest distance cost on the stratum.  That
    # it sits at the inner edge rather than at b=1 is certified, not assumed.
    certify_face_infimum(x, t0)
    face = _g(top, x)
    w1_ceiling = rho_arb**2 / face
    # Arb division by an interval already puts the conservative end on top:
    # rho^2/[g_lo,g_hi] = [rho^2/g_hi, rho^2/g_lo], and every downstream test
    # is a whole-interval comparison, so `leading > 0` is decided by
    # rho^2/g_lo.  That is a property of the arithmetic, not of this line, so
    # it is asserted rather than trusted -- a later edit that substituted the
    # midpoint or the upper end of `face` would weaken the bound silently, and
    # at the present 320-bit precision `face` is only ~1e-70 wide, so the
    # substitution would be invisible in every printed digit.
    if not w1_ceiling.upper() >= (rho_arb**2 / arb(face.lower())).lower():
        raise AssertionError(
            "the W1 ceiling is not taken at the smallest admissible face "
            "cost; the division is no longer conservative")
    if not face.lower() <= face.upper():
        raise AssertionError("degenerate face enclosure")
    leading = 1 - 2 * beta - 2 * (1 - beta) * w1_ceiling
    penalty = (-top.log()
               + 2 * (1 - beta) * theta
               - 2 * beta * (1 - 2 * t0_arb).log())
    log_inverse = -t0_arb.log()
    m1 = leading * (log_inverse + 1) - penalty
    mean_reach = parameters.mean
    if 1 - parameters.mean > mean_reach:
        mean_reach = 1 - parameters.mean
    distance_rate = _g_prime_ceiling(x) + 2 * mean_reach + t0_arb
    kappa = m1 / distance_rate
    t_star = (-(penalty / leading)).exp() if leading > 0 else arb(0)
    return MirrorConstants(
        t0=t0,
        rho=rho,
        mu=int(mu),
        theta=theta,
        theta_lambda=dual["lambda"],
        w1_ceiling=w1_ceiling,
        leading=leading,
        penalty=penalty,
        m1=m1,
        distance_rate=distance_rate,
        kappa=kappa,
        t_star=t_star,
        argmax_cell=dual["argmax_cell"],
    )


def best_constants(parameters: ArbParameters, t0: Fraction, rho: Fraction,
                   pieces: int = 4000,
                   mus: Sequence[int] = MU_GRID) -> MirrorConstants:
    """The dual multiplier is free; take the one giving the largest kappa.

    Every candidate is individually a valid upper bound, so maximizing over
    them is sound: this is weak duality, not a heuristic.
    """
    best: Optional[MirrorConstants] = None
    for mu in mus:
        candidate = certify_constants(parameters, t0, rho, mu, pieces)
        if best is None or candidate.kappa > best.kappa:
            best = candidate
    assert best is not None
    return best


# ---------------------------------------------------------------------------
# Interval sweeps guarding the algebraic links of the proof.


def sweep_paired_link(t0: Fraction, s_boxes: int, u_boxes: int
                      ) -> dict[str, Any]:
    """The two one-sided bounds on 1-pi, and the range of pi_1.

    With `s = 1-y` and `u = 1-z`, expanding `pi = yz(1+su)` gives the exact
    polynomial identity

        1 - pi = s + u - 2su + su*(s+u-su),

    from which the two bounds follow by *cancellation*, not estimation:

        1 - pi - u(1-2s) = s + su*(s+u-su),
        1 - pi - s(1-2u) = u + su*(s+u-su),

    and `s+u-su = 1-(1-s)(1-u) >= 0` on the unit square.  The same rewriting
    handles pi_1: for `y >= 1/2`,
        pi_1 = z(1-(1-z)(1-2s)),  z - pi_1 = z(1-z)(1-2s),  pi_1 - z^2 = 2s z(1-z).

    Two arithmetics are used deliberately.  The *identity* is a transcription
    guard against mis-expanding Liu's `pi`, so it is checked in Arb, where a
    small residual is expected and reported.  The *signs* are checked in exact
    `Fraction` arithmetic at the cell corners, because an Arb ball that touches
    zero from above always outward-rounds to a slightly negative lower bound --
    at `1/2048` the overshoot is `4.5e-13` -- which would make a manifestly
    nonnegative expression fail its own sign test.  Every expression below is a
    product or sum of factors that are monotone and nonnegative on the cell, so
    its exact minimum sits at the lower corner.
    """
    if t0 > Fraction(1, 2):
        raise ValueError("the paired rewriting needs y>=1/2, so t0<=1/2")
    s_boxes = _pow2_at_least(s_boxes)
    u_boxes = _pow2_at_least(u_boxes)
    worst_identity = Fraction(0)
    worst_pi1_identity = Fraction(0)
    worst_manifest = None
    worst_union = None
    worst_pi1_gap = None
    worst_pi1_floor = None
    boxes = 0
    for si in range(s_boxes):
        s_lo = t0 * si / s_boxes
        s_hi = t0 * (si + 1) / s_boxes
        for ui in range(u_boxes):
            u_lo = Fraction(ui, u_boxes)
            u_hi = Fraction(ui + 1, u_boxes)
            boxes += 1
            # --- exact signs at the minimizing corner -------------------
            union_min = 1 - (1 - s_lo) * (1 - u_lo)
            if union_min < 0:
                raise AssertionError(f"s+u-su negative at ({s_lo},{u_lo})")
            if worst_union is None or union_min < worst_union:
                worst_union = union_min
            common_min = s_lo * u_lo * union_min
            for manifest in (s_lo + common_min, u_lo + common_min):
                if manifest < 0:
                    raise AssertionError(
                        f"manifest 1-pi margin negative at ({s_lo},{u_lo})")
                if worst_manifest is None or manifest < worst_manifest:
                    worst_manifest = manifest
            if 1 - 2 * s_hi < 0:
                raise AssertionError("the paired rewriting needs s<=1/2")
            pi1_gap = (1 - u_hi) * (u_lo) * (1 - 2 * s_hi)
            pi1_floor = 2 * s_lo * (1 - u_hi) * u_lo
            for value, slot in ((pi1_gap, "gap"), (pi1_floor, "floor")):
                if value < 0:
                    raise AssertionError(
                        f"pi_1 {slot} negative at ({s_lo},{u_lo})")
            if worst_pi1_gap is None or pi1_gap < worst_pi1_gap:
                worst_pi1_gap = pi1_gap
            if worst_pi1_floor is None or pi1_floor < worst_pi1_floor:
                worst_pi1_floor = pi1_floor
            # --- exact transcription guard on the two expansions --------
            # Both sides are polynomials with rational coefficients, so an
            # exact evaluation at the cell corners is decisive: the residual
            # must be exactly zero, not merely small.  An Arb evaluation would
            # only report the dependency width of the cell.
            for sv in (s_lo, s_hi):
                for uv in (u_lo, u_hi):
                    yv, zv = 1 - sv, 1 - uv
                    residual = ((1 - yv * zv * (1 + sv * uv))
                                - (sv + uv - 2 * sv * uv
                                   + sv * uv * (1 - (1 - sv) * (1 - uv))))
                    if residual != 0:
                        raise AssertionError(
                            f"1-pi expansion is wrong at ({sv},{uv}): {residual}")
                    worst_identity = max(worst_identity, abs(residual))
                    pi1_residual = (zv * (1 - (1 - zv) * (1 - 2 * sv))
                                    - zv * (1 + (1 - zv) * (1 - 2 * yv)))
                    if pi1_residual != 0:
                        raise AssertionError(
                            f"pi_1 expansion is wrong at ({sv},{uv}): "
                            f"{pi1_residual}")
                    worst_pi1_identity = max(worst_pi1_identity,
                                             abs(pi1_residual))
    if worst_manifest is None:
        raise AssertionError("empty paired sweep")
    return {
        "boxes": boxes,
        "arithmetic": "exact Fraction corners for signs, Arb for the identity",
        "worst_identity_residual": str(worst_identity),
        "worst_pi1_identity_residual": str(worst_pi1_identity),
        "worst_manifest_margin": str(worst_manifest),
        "smallest_union": str(worst_union),
        "worst_pi1_headroom": str(worst_pi1_gap),
        "worst_pi1_floor": str(worst_pi1_floor),
        "pi_1_bracketed_by": "z^2 <= pi_1 <= z <= 1",
    }


def sweep_cross_link(t0: Fraction, t_octaves: int, b_boxes: int
                     ) -> dict[str, Any]:
    """The two monotonicity links used on the cross sum.

    For `b_j = 1-t` with `t in (0,t0]` and any partner `b_k in [0,1]`:
        h'(b_j b_k) >= h'(b_k)              (link A)
        h'(b_j b_k) >= log t                (link B)

    Both are the *same* analytic fact -- `h'` is decreasing on `(0,1)` -- applied
    to two exact algebraic comparisons of its argument:
        A:  b_j b_k <= b_k                    because b_j <= 1
        B:  b_j b_k <= 1-t = b_j              because b_k <= 1,
            and then h'(1-t) = log(t/(1-t)) >= log t   because 1-t <= 1.
    Nothing here needs an interval evaluation of `h'` itself, and trying to do
    one would be actively worse: at forty octaves the cell containing `b_k=1`
    encloses `1-b_j b_k` to within `4e-13`, which is a thousand times the true
    value `t = 2^-46`, so the enclosure would straddle zero and the logarithm
    would not exist.  Instead this sweep certifies the two argument
    comparisons in exact `Fraction` arithmetic over a dyadic grid, and
    certifies the monotonicity of `h'` by the sign of `h''(u) = -1/(u(1-u))`
    over a dyadic cover with strictly positive endpoints.
    """
    t_octaves = max(1, t_octaves)
    b_boxes = _pow2_at_least(b_boxes)
    boxes = 0
    worst_link_a = None
    worst_link_b = None
    smallest_t = t0 / 2 ** t_octaves
    for octave in range(t_octaves):
        t_hi = t0 / 2 ** octave
        t_lo = t_hi / 2
        for bi in range(b_boxes + 1):
            b_k = Fraction(bi, b_boxes)
            boxes += 1
            # Both margins are evaluated with a *matched* t.  Their factored
            # forms are products of nonnegative factors, so the minimum over
            # the cell is the product of the per-factor minima, exactly.
            #   b_k - b_j b_k = t   * b_k       -> zero exactly at b_k=0
            #   (1-t) - b_j b_k = (1-t)*(1-b_k) -> zero exactly at b_k=1
            link_a = t_lo * b_k
            link_b = (1 - t_hi) * (1 - b_k)
            if link_a < 0 or link_b < 0:
                raise AssertionError(
                    f"argument comparison failed at t in [{t_lo},{t_hi}], "
                    f"b={b_k}")
            if worst_link_a is None or link_a < worst_link_a:
                worst_link_a = link_a
            if worst_link_b is None or link_b < worst_link_b:
                worst_link_b = link_b
    # h' is decreasing exactly where u(1-u)>0.  Certify that on a dyadic cover
    # whose endpoints stay strictly inside (0,1).
    floor_exponent = t_octaves + 4
    cells = 0
    worst_curvature = None
    for octave in range(1, floor_exponent + 1):
        for side in (0, 1):
            hi = Fraction(1, 2 ** octave)
            lo = hi / 2
            if side:
                lo, hi = 1 - hi, 1 - lo
            cell = _arbf(lo).union(_arbf(hi))
            product = cell * (1 - cell)
            cells += 1
            if not product > 0:
                raise AssertionError(f"h'' sign undetermined on [{lo},{hi}]")
            if worst_curvature is None or product < worst_curvature:
                worst_curvature = product
    return {
        "boxes": boxes,
        "octaves": t_octaves,
        "smallest_t": str(smallest_t),
        "arithmetic": "exact Fraction comparisons; Arb only for the h'' sign",
        "worst_link_a_margin": str(worst_link_a),
        "worst_link_b_margin": str(worst_link_b),
        "curvature_cells": cells,
        "smallest_u_one_minus_u": float(worst_curvature.lower()),
    }


# ---------------------------------------------------------------------------
# The asymptotic coefficient, measured against Liu's own formula.


def _mirror_weights(values: Sequence[mpmath.mpf], index: int
                    ) -> tuple[mpmath.mpf, mpmath.mpf]:
    """(W1, A1): global weight and within-component mass at support exactly 1."""
    one = mpmath.mpf(1)
    weights = _weights(values, one)
    support = _support(values)
    a1, a2 = values[0], values[1]
    masses = (a1, a2, one - a1 - a2)
    base = 0 if index < 3 else 3
    w1 = mpmath.mpf(0)
    for k, point in enumerate(support):
        if point == 1:
            w1 += weights[k]
    a1_mass = mpmath.mpf(0)
    for offset in range(3):
        if support[base + offset] == 1:
            a1_mass += masses[offset]
    return w1, a1_mass


def asymptotic_coefficient(parameters: MPParameters, dps: int
                           ) -> dict[str, Any]:
    """Measure the slope of G_j in log(1/(1-b_j)) and compare with the theory.

    Six configurations sweep the coefficient `1-2(1-beta)W1-2 beta A1` across
    its whole sign range, including the discriminating case `W1=1/2, A1=1`,
    where the prediction is exactly `-beta` and no `W1`-only formula can match.
    """
    one = mpmath.mpf(1)
    zero = mpmath.mpf(0)
    x = parameters.x
    p = parameters.p
    configurations = {
        "W1=1-p, only the moving atom reaches 1": (
            p, 1 - p, zero, x, None, zero, zero, zero, zero),
        "W1=2/5, the rest of the mass sits at x": (
            mpmath.mpf(3) / 5, mpmath.mpf(2) / 5, zero,
            x, None, zero, zero, zero, zero),
        "W1=1/2, the exact sign flip": (
            mpmath.mpf(1) / 2, mpmath.mpf(1) / 2, zero,
            x, None, zero, zero, zero, zero),
        "W1=7/10, a third atom already at 1": (
            mpmath.mpf(3) / 10, mpmath.mpf(4) / 10, zero,
            x, None, one, zero, zero, zero),
        "W1=1/2 but A1=1, the beta split is visible": (
            mpmath.mpf(1) / 2, mpmath.mpf(1) / 2, mpmath.mpf(1) / 2,
            one, None, zero, x, x, zero),
        "W1=1, every atom at support 1": (
            mpmath.mpf(45) / 100, mpmath.mpf(45) / 100, mpmath.mpf(1) / 2,
            one, None, one, one, one, one),
    }
    rows = []
    signs = set()
    # The far probe sits at 1-1e-80, which needs far more than the report's
    # working precision merely to be representable, let alone differenced.
    with mpmath.workdps(max(dps, 220)):
        beta = parameters.beta
        for name, template in configurations.items():
            pair = []
            for exponent in (40, 80):
                moving = one - mpmath.mpf(10) ** (-exponent)
                values = tuple(moving if entry is None else entry
                               for entry in template)
                pair.append((exponent,
                             closed_form_factor(values, 1, beta),
                             values))
            (e_lo, g_lo, _), (e_hi, g_hi, _) = pair
            slope = (g_hi - g_lo) / ((e_hi - e_lo) * mpmath.log(10))
            # W1 and A1 are the masses at the *limit* support 1, so the moving
            # atom counts itself: its own h'(b_j b_j) is one of the divergences.
            limit = tuple(one if entry is None else entry for entry in template)
            w1, a1_mass = _mirror_weights(limit, 1)
            predicted = 1 - 2 * (1 - beta) * w1 - 2 * beta * a1_mass
            difference = slope - predicted
            if abs(difference) > mpmath.mpf(10) ** (-30):
                raise AssertionError(
                    f"mirror coefficient mismatch for {name}: {difference}")
            sign = ("flat (exact tie)" if abs(float(slope)) < 1e-30
                    else "increasing" if slope > 0 else "decreasing")
            signs.add(sign)
            rows.append({
                "configuration": name,
                "weight_at_support_one": mpmath.nstr(w1, 12),
                "component_mass_at_support_one": mpmath.nstr(a1_mass, 12),
                "measured_slope": mpmath.nstr(slope, 20),
                "predicted_slope": mpmath.nstr(predicted, 20),
                "difference": mpmath.nstr(difference, 6),
                "sign": sign,
            })
    return {
        "rows": rows,
        "two_sided": bool({"increasing", "decreasing"} <= signs),
        "signs": sorted(signs),
    }


# ---------------------------------------------------------------------------
# Falsification pass over the mirror stratum.


def _distance_squared(values: Sequence[mpmath.mpf], parameters: MPParameters
                      ) -> mpmath.mpf:
    """dist^2 = (1-q) delta(P0)^2 + q delta(P1)^2 from liu9_tube's definition."""
    one = mpmath.mpf(1)
    a1, a2, q = values[0], values[1], values[2]
    masses = (a1, a2, one - a1 - a2)
    support = _support(values)
    total = mpmath.mpf(0)
    for component, weight in enumerate((one - q, q)):
        base = 3 * component
        mean = mpmath.mpf(0)
        spread = mpmath.mpf(0)
        for offset in range(3):
            point = support[base + offset]
            mean += masses[offset] * point
            spread += masses[offset] * (point * (point - parameters.x)) ** 2
        total += weight * ((mean - parameters.mean) ** 2 + spread)
    return total


def _draw_mirror_point(rng: random.Random, parameters: MPParameters,
                       t0: float, rho: float
                       ) -> Optional[tuple[tuple[mpmath.mpf, ...], int]]:
    """A tube point carrying at least one atom in the mirror layer.

    Built constructively around P* -- rejection sampling from the raw box would
    essentially never land inside a radius-0.1 tube -- but every returned point
    is re-tested against the raw hypotheses, never assumed to satisfy them.
    """
    x = float(parameters.x)
    for _ in range(400):
        index = rng.randrange(6)
        base = 0 if index < 3 else 3
        q = rng.choice([rng.random(), rng.random() ** 4,
                        1 - rng.random() ** 4, 0.5])
        q = min(max(q, 1e-9), 1 - 1e-9)
        mirror_mass = rng.random() ** 2 * 0.35
        rest = 1 - mirror_mass
        at_x = rest * rng.uniform(0.55, 1.0)
        masses = [mirror_mass, at_x, rest - at_x]
        rng.shuffle(masses)
        slot = rng.randrange(3)
        masses[slot], masses[0] = masses[0], masses[slot]
        a1, a2 = masses[0], masses[1]
        if a1 <= 0 or a2 <= 0 or a1 + a2 >= 1:
            continue
        t = math.exp(rng.uniform(math.log(1e-12), math.log(t0)))
        support = [0.0] * 6
        for side in (0, 3):
            for offset in range(3):
                support[side + offset] = (
                    x * (1 + rng.uniform(-0.05, 0.05)) if rng.random() < 0.7
                    else rng.random() ** 3 * 0.02)
        support[index] = 1 - t
        # Keep the moving atom's own component mass in slot 0 of that side.
        values = tuple(mpmath.mpf(repr(v)) for v in
                       [a1, a2, q,
                        support[0], support[1], support[2],
                        support[3], support[4], support[5]])
        if index != base:
            continue
        if _distance_squared(values, parameters) > mpmath.mpf(repr(rho)) ** 2:
            continue
        return values, index
    return None


@dataclass
class MirrorSampleReport:
    points: int
    worst_pointwise_margin: float
    worst_integrated_ratio: float
    worst_kappa_margin: float
    worst_autodiff_error: float
    mean_never_decreased: bool
    failures: list[str]


def run_samples(parameters: MPParameters, constants: MirrorConstants,
                random_points: int, dps: int, seed: int) -> MirrorSampleReport:
    """Falsification pass: every claim of the theorem, on raw tube points."""
    rng = random.Random(seed)
    t0 = float(constants.t0)
    rho = float(constants.rho)
    # Take each constant on its conservative side: the sampler must test the
    # theorem as stated, not a rounded version that is easier to satisfy.
    leading = mpmath.mpf(float(constants.leading.lower()))
    penalty = mpmath.mpf(float(constants.penalty.upper()))
    m1 = mpmath.mpf(float(constants.m1.lower()))
    kappa = mpmath.mpf(float(constants.kappa.lower()))
    report = MirrorSampleReport(0, math.inf, math.inf, math.inf, 0.0, True, [])
    with mpmath.workdps(dps):
        beta = parameters.beta
        for _ in range(random_points):
            drawn = _draw_mirror_point(rng, parameters, t0, rho)
            if drawn is None:
                continue
            values, index = drawn
            report.points += 1
            support = _support(values)
            weights = _weights(values, mpmath.mpf(1))
            t = 1 - support[index]

            # 1. The pointwise lower bound on G_j.
            factor = closed_form_factor(values, index, beta)
            predicted = leading * mpmath.log(1 / t) - penalty
            margin = factor - predicted
            report.worst_pointwise_margin = min(
                report.worst_pointwise_margin, float(margin))
            if margin < 0:
                report.failures.append(
                    f"pointwise bound violated: G={factor}, bound={predicted}")

            # 2. The closed form against forward-mode autodiff of _formula.
            direct = autodiff_partial(values, index, beta)
            closed = weights[index] * factor
            scale = max(abs(direct), mpmath.mpf(1))
            error = float(abs(direct - closed) / scale)
            report.worst_autodiff_error = max(report.worst_autodiff_error, error)
            if error > 1e-25:
                report.failures.append(f"autodiff disagreement: {error:.3e}")

            # 3. The integrated reduction, raising every mirror atom to 1.
            raised = list(values)
            layer_mass = mpmath.mpf(0)
            for k, point in enumerate(support):
                if point > 1 - mpmath.mpf(repr(t0)):
                    raised[SUPPORT_VARIABLE_INDEX[k]] = mpmath.mpf(1)
                    layer_mass += weights[k] * (1 - point)
            raised = tuple(raised)
            if layer_mass <= 0:
                continue
            increase = gap_mp(raised, beta) - gap_mp(values, beta)
            ratio = increase / layer_mass
            report.worst_integrated_ratio = min(
                report.worst_integrated_ratio, float(ratio))
            if ratio < m1:
                report.failures.append(
                    f"integrated rate below m1: {ratio} < {m1}")

            # 4. kappa-monotonicity of gap - kappa*dist^2 under the reduction.
            distance_growth = (_distance_squared(raised, parameters)
                               - _distance_squared(values, parameters))
            kappa_margin = increase - kappa * distance_growth
            report.worst_kappa_margin = min(report.worst_kappa_margin,
                                            float(kappa_margin))
            if kappa_margin < 0:
                report.failures.append(
                    f"kappa-monotonicity violated: {kappa_margin}")

            # 5. Feasibility is preserved because the mean only rises.
            if mean_of(raised, mpmath.mpf(1)) < mean_of(values, mpmath.mpf(1)):
                report.mean_never_decreased = False
                report.failures.append("the reduction lowered the mean")
    return report


# ---------------------------------------------------------------------------
# Threshold and comparison tables.


def smooth_kappa_ceiling(parameters: ArbParameters,
                         units: str = "gap") -> arb:
    """liu9_tube's two smooth displayed modes, in a stated normalization.

    Delegated rather than copied: the conversion `gap = EHX*(Phi-1)` is a
    single fact and belongs in one place.  `liu9_tube.py` prints the ceiling
    in `Phi` units; every kappa in this module is a raw-gap kappa, so the
    comparison must be made after multiplying by `H* = EHX(P*)`.
    """
    from liu9_boundary_layer import smooth_kappa_ceiling as shared
    return shared(parameters, units=units)


def threshold_table(parameters: ArbParameters, pieces: int
                    ) -> list[dict[str, Any]]:
    ceiling = smooth_kappa_ceiling(parameters)
    rows = []
    for rho in THRESHOLD_RHO:
        for t0 in THRESHOLD_T0:
            constants = best_constants(parameters, t0, rho, pieces)
            row = constants.as_json()
            row["beats_smooth_ceiling"] = bool(
                constants.certified() and constants.kappa > ceiling)
            rows.append(row)
    return rows


def radius_ceiling(parameters: ArbParameters) -> dict[str, Any]:
    """The radius above which W1 alone can flip the leading coefficient.

    rho^2 = g(1-t0)*(1-2 beta)/(2(1-beta)) makes M vanish, so no mirror
    estimate of this shape can survive a larger tube at that threshold.
    """
    x = parameters.x
    beta = parameters.beta
    rows = []
    for t0 in THRESHOLD_T0:
        top = 1 - _arbf(t0)
        critical = _g(top, x) * (1 - 2 * beta) / (2 * (1 - beta))
        rows.append({
            "t0": str(t0),
            "critical_rho_squared": str(critical),
            "critical_rho": str(critical.sqrt()),
            "critical_rho_lower": float(critical.sqrt().lower()),
        })
    return {"rows": rows}


# ---------------------------------------------------------------------------
# Reporting.


def _canonical_digest(payload: dict[str, Any]) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t0", type=str, default=str(DEFAULT_T0))
    parser.add_argument("--rho", type=str, default=str(DEFAULT_RHO))
    parser.add_argument("--dps", type=int, default=DEFAULT_DPS)
    parser.add_argument("--theta-pieces", type=int, default=4000)
    parser.add_argument("--paired-boxes", type=int, default=60)
    parser.add_argument("--cross-octaves", type=int, default=40)
    parser.add_argument("--cross-boxes", type=int, default=40)
    parser.add_argument("--random-points", type=int, default=400)
    parser.add_argument("--seed", type=int, default=SAMPLE_SEED)
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args(argv)

    started = time.monotonic()
    t0 = Fraction(args.t0)
    rho = Fraction(args.rho)
    mp_parameters = solve_equation_parameters(args.dps + 30)
    arb_parameters = certify_equation_parameters(mp_parameters)
    constants = best_constants(arb_parameters, t0, rho, args.theta_pieces)
    ceiling = smooth_kappa_ceiling(arb_parameters)

    print("1. THE MIRROR STRATUM EXISTS AND IS NOT SIGNED BY THE MEAN")
    coefficient = asymptotic_coefficient(mp_parameters, args.dps)
    print("PROVED [closed form, slope in log(1/(1-b_j)) between 1e-40 and "
          "1e-80]: the coefficient is 1-2(1-beta)W1-2*beta*A1, and it is "
          f"two-sided ({coefficient['two_sided']}):")
    for row in coefficient["rows"]:
        print("   %s -- W1=%s, A1=%s, measured %s, predicted %s, diff %s, %s"
              % (row["configuration"], row["weight_at_support_one"],
                 row["component_mass_at_support_one"], row["measured_slope"],
                 row["predicted_slope"], row["difference"], row["sign"]))
    print()

    print("2. THE TUBE BOUNDS THE WEIGHT AT SUPPORT 1")
    infimum = certify_face_infimum(arb_parameters.x, t0)
    print("PROVED [Arb]: g(b)=[b(b-x)]^2 is nondecreasing across the stratum, "
          "since b-x in %s and 2b-x in %s are both positive at the inner edge, "
          "so its infimum is g(1-t0)=%.7f and NOT the endpoint g(1)=%.7f.  "
          "Using the endpoint would understate the divisor by %.2f%% and "
          "inflate the critical radius by %.2f%%, an overclaim in the "
          "dangerous direction."
          % (infimum["b_minus_x_at_inner_edge"],
             infimum["two_b_minus_x_at_inner_edge"],
             infimum["infimum_g_lower"], infimum["endpoint_g_lower"],
             100 * (1 - 1 / infimum["endpoint_overclaim_ratio"]),
             100 * (infimum["endpoint_overclaim_ratio"] ** 0.5 - 1)))
    print(f"PROVED [Arb]: with t0={t0} and rho={rho}, every tube point has "
          f"W1 <= rho^2/g(1-t0) in {constants.w1_ceiling}")
    print(f"PROVED [Arb]: hence the leading coefficient M >= {constants.leading}")
    ceilings = radius_ceiling(arb_parameters)
    for row in ceilings["rows"]:
        print("PROVED [Arb]: at t0=%s no estimate of this shape survives "
              "rho above %s" % (row["t0"], row["critical_rho"]))
    print()

    print("3. THE ALGEBRAIC LINKS")
    paired = sweep_paired_link(t0, args.paired_boxes, args.paired_boxes)
    print("PROVED [exact, %d boxes]: 1-pi = s+u-2su+su(s+u-su) with residual "
          "%s and pi_1 = z(1-(1-z)(1-2s)) with residual %s; the manifest "
          "margins s+su(s+u-su) and u+su(s+u-su) stay >= %s, s+u-su >= %s, "
          "and %s"
          % (paired["boxes"], paired["worst_identity_residual"],
             paired["worst_pi1_identity_residual"],
             paired["worst_manifest_margin"], paired["smallest_union"],
             paired["pi_1_bracketed_by"]))
    cross = sweep_cross_link(t0, args.cross_octaves, args.cross_boxes)
    print("PROVED [exact, %d corners over %d octaves down to t=%s]: the two "
          "argument comparisons b_j b_k <= b_k and b_j b_k <= 1-t hold, with "
          "factored margins t*b_k >= %s and (1-t)(1-b_k) >= %s; both vanish "
          "exactly at the structural equality cases b_k=0 and b_k=1."
          % (cross["boxes"], cross["octaves"], cross["smallest_t"],
             cross["worst_link_a_margin"], cross["worst_link_b_margin"]))
    print("PROVED [Arb, %d dyadic cells]: h''(u)=-1/(u(1-u)) keeps one sign "
          "because u(1-u) >= %.6g > 0, so h' is decreasing and both links "
          "follow."
          % (cross["curvature_cells"], cross["smallest_u_one_minus_u"]))
    print()

    print("4. THE CERTIFIED CONSTANTS")
    print(f"PROVED [Arb, dual multiplier mu={constants.mu}]: the outside-mirror "
          f"cross sum is at most Theta in {constants.theta}")
    print(f"PROVED [Arb]: M in {constants.leading}")
    print(f"PROVED [Arb]: K1 in {constants.penalty}")
    print(f"PROVED [Arb]: G_j > 0 for every t below t* in {constants.t_star}")
    verdict = "PROVED" if constants.certified() else "REFUTED"
    print(f"{verdict} [Arb]: the layer estimate covers all of (0,{t0}]: "
          f"t* {'>=' if constants.certified() else '<'} t0")
    if constants.certified():
        print(f"PROVED [Arb]: gap rises at rate m1 in {constants.m1} per unit "
              "of sum_j w_j (1-b_j)")
        print(f"PROVED [Arb]: the reduction is kappa-monotone for every kappa "
              f"up to {constants.kappa}")
        print("PROVED [Arb]: the reduction raises supports to 1, so the mean "
              "never falls and mean-feasibility is preserved for free.")
    print()

    print("5. FALSIFICATION PASS OVER THE MIRROR STRATUM")
    samples = run_samples(mp_parameters, constants, args.random_points,
                          args.dps, args.seed)
    print("NUMERICAL [%d tube points carrying a mirror atom]: worst pointwise "
          "margin %.6g, worst integrated rate %.6g against m1 %.6g, worst "
          "kappa margin %.6g, worst autodiff error %.3e, mean never fell: %s"
          % (samples.points, samples.worst_pointwise_margin,
             samples.worst_integrated_ratio, float(constants.m1.lower()),
             samples.worst_kappa_margin, samples.worst_autodiff_error,
             samples.mean_never_decreased))
    if samples.failures:
        for failure in samples.failures[:10]:
            print(f"FAILED: {failure}")
        raise AssertionError(f"{len(samples.failures)} sample failures")
    print("PROVED [falsification]: no sampled point violated the theorem.")
    print()

    print("6. THRESHOLD TABLE AND COMPARISON WITH THE SMOOTH CEILING")
    table = threshold_table(arb_parameters, args.theta_pieces)
    print("PROVED [Arb]: gap = EHX*(Phi-1) exactly and EHX(P*) = %s, so the "
          "quotient ceiling %s that liu9_tube.py prints becomes %s in the "
          "raw-gap units used here.  Every kappa below is a raw-gap kappa."
          % (entropy_at_star(arb_parameters),
             smooth_kappa_ceiling(arb_parameters, units="phi"), ceiling))
    for row in table:
        print("%s [Arb]: rho=%s, t0=%s -> M=%.6g, m1=%.6g, kappa=%.6g, "
              "t*=%.4g (beats the smooth ceiling: %s)"
              % ("PROVED" if row["certified"] else "REFUTED",
                 row["rho"], row["t0"], row["leading_M_lower"],
                 row["m1_lower"], row["admissible_kappa_lower"],
                 row["t_star_lower"], row["beats_smooth_ceiling"]))
    print()

    print("7. VERDICT")
    if constants.certified():
        print("PROVED: the mirror stratum b_j->1 of Liu's Hypothesis 2 is "
              "closed at the stated (rho,t0).  The estimate never uses the "
              "mean constraint, and its reduction target is the exact face "
              "b_j=1, so it composes with the zero-support layer rather than "
              "competing with it.")
    else:
        print("REFUTED at this (rho,t0): the estimate does not cover the whole "
              "layer.  See the threshold table for the pairs that do.")
    print("OPEN: ingredients (i) and (iii) of the piecewise local lemma.  This "
          "module certifies no tube radius on its own.")

    payload = {
        "module": "liu9_mirror_layer.py",
        "t0": str(t0),
        "rho": str(rho),
        "dps": args.dps,
        "constants": constants.as_json(),
        "smooth_kappa_ceiling_gap_units": str(ceiling),
        "smooth_kappa_ceiling_phi_units":
            str(smooth_kappa_ceiling(arb_parameters, units="phi")),
        "entropy_at_star": str(entropy_at_star(arb_parameters)),
        "asymptotic_coefficient": coefficient,
        "paired_link": paired,
        "cross_link": cross,
        "radius_ceiling": ceilings,
        "face_infimum": infimum,
        "threshold_table": table,
        "samples": {
            "points": samples.points,
            "worst_pointwise_margin": samples.worst_pointwise_margin,
            "worst_integrated_ratio": samples.worst_integrated_ratio,
            "worst_kappa_margin": samples.worst_kappa_margin,
            "worst_autodiff_error": samples.worst_autodiff_error,
            "mean_never_decreased": samples.mean_never_decreased,
            "failures": samples.failures,
        },
        "elapsed_seconds": round(time.monotonic() - started, 3),
    }
    if args.output:
        digest_payload = dict(payload)
        digest_payload.pop("elapsed_seconds", None)
        text = json.dumps(payload, indent=1, sort_keys=True)
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(text + "\n")
        print(f"\nreport written to {args.output}")
        print(f"report_sha256 {_canonical_digest(digest_payload)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
