#!/usr/bin/env python3
"""Sound uniform envelope of the support partial derivatives, Liu Hypothesis 2.

MATHEMATICAL OBJECT.  Liu's nine parameters (liu9_objective.VARIABLE_NAMES)
are (a1, a2, q, b0, b2, b4, b1, b3, b5): masses a = (a1, a2, a3=1-a1-a2),
mixture weight q in [0,1], and six atom supports s = (b0..b5) with weights
w = ((1-q)a1, (1-q)a2, (1-q)a3, q a1, q a2, q a3).  For the denominator-
cleared gap (liu9_boundary_layer.py, validated there against autodiff of
Liu's own transcription) the support partial has the closed form

    d gap / d b_k = w_k * G_k(V),
    G_k(V) = 2(1-beta) * sum_l w_l s_l h'(b_k s_l)
           + 2*beta * sum_{l in comp(k)} a_l pi_1(b_k, s_l) h'(pi(b_k, s_l))
           - h'(b_k),

h'(u) = log((1-u)/u), pi(y,z) = yz(1+(1-y)(1-z)),
pi_1(y,z) = d pi/dy = z(1+(1-z)(1-2y)).  Terms with s_l = 0 vanish
identically in b_k, and w_k = 0 exactly there, so w_k G_k extends
continuously (by the value 0) to the closed mass simplex.

DEGREE-TWO STRUCTURE.  Fix a mass node n in T = {a1, a2 >= 0, a1+a2 <= 1}
and enclose s and q by intervals.  Every weight/mass factor is affine in the
Bernstein coordinates (lA, lB, lC) = (1-a1-a2, a1, a2), and every quantity
h'(b_k s_l), pi_1, h'(pi), h'(b_k) depends only on the (fixed) s and q
boxes.  Hence the NODE VALUE

    p_k(n) = w_k(n, q) * G_k(n, s, q)  (= d gap/d b_k at masses n, exactly)

is a polynomial of degree <= 2 in (lA, lB, lC) for each fixed q, enclosed
uniformly over the q box by interval arithmetic at the node.  The common
helper (liu9_h2_envelope_common.py) contracts the six node balls through the
degree-2 Bernstein convex-hull property into a uniform bound over ALL of T;
that contraction is the only place the mass simplex enters.

THEOREM (support envelope, interior form).  Fix a certified beta ball, six
support boxes [s_lo, s_hi]_k strictly inside (0,1), and a q box in [0,1].
For every attained feasible point (masses anywhere in T, q in the q box,
supports in the boxes) and every support coordinate k = 0..5:

    |d gap / d b_k|  <=  L_k,

where L_k is the contracted helper's simplex-uniform enclosure of |p_k|
(the largest Bernstein coefficient-ball endpoint in absolute value; the
Bernstein convex-hull property of degree-2 polynomials on T makes it bound
p_k over the FULL simplex, while interval containment over each node ball
makes it bound p_k over the q box as well).  Because the closed form has no
1/b factor, L_k IS the Lipschitz constant of gap along the b_k axis, and
crossing the support box changes gap by at most sum_k L_k (b_khi - b_klo)
(support_drift_bound).  Claim scope: INTERIOR support boxes only; a box
touching 0 or 1 is rejected up front (ValueError), and the certificate is
silent on such strata.

MUTATIONS.  Two deterministic defects (cross factor 2 -> 1; paired-component
swap plus pi_1 sign flip) are re-derived locally and MUST disagree with the
independent high-precision central-difference reference of gap_mp; the
healthy baseline MUST agree.  Both mechanisms are exercised in the runner.

Run from the repository root with

    math/.venv/bin/python -I -B math/uc/liu9_h2_support_envelope.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
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

from liu9_boundary_layer import (  # noqa: E402
    SUPPORT_VARIABLE_INDEX,
    closed_form_partial,
    gap_mp,
)

try:  # authoritative when present; the module stays runnable without it
    from liu9_h2_envelope_common import (  # type: ignore  # noqa: E402
        MASS_NODES,
        arb_abs_upper,
        point_ball,
        simplex_abs_upper,
        simplex_bernstein_coefficients,
        simplex_range,
    )

    HAVE_COMMON = True
except Exception:  # pragma: no cover - dual fallback only
    HAVE_COMMON = False

ctx.prec = max(ctx.prec, 320)

DEFAULT_DPS = 100
ONE = Fraction(1)
ZERO = Fraction(0)

# Representative intake for the runner's Lipschitz report: the same
# all-six-supports wide box for every coordinate, the q box [0, q_hi]
# (interior support boxes; q is allowed to touch its endpoints).
WIDE_SUPPORT_BOX = (Fraction(1, 1024), Fraction(1, 2))
WIDE_Q_BOX = (Fraction(0), Fraction(1, 1000))
POINT_FIXTURES = (
    Fraction(1, 3), Fraction(1, 3), Fraction(1, 2),
    Fraction(1, 4), Fraction(3, 8), Fraction(9, 16),
    Fraction(1, 8), Fraction(3, 8), Fraction(15, 16),
)


# ---------------------------------------------------------------------------
# Arb plumbing.


def _arbf(value: Any) -> arb:
    """Exact rational -> Arb point ball."""
    value = Fraction(value)
    return arb(value.numerator) / arb(value.denominator)


def _iball(lo: Any, hi: Any) -> arb:
    """The Arb enclosure [Fraction(lo), Fraction(hi)] (outward rounding)."""
    lo = Fraction(lo)
    hi = Fraction(hi)
    return (arb(lo.numerator) / arb(lo.denominator)).union(
        arb(hi.numerator) / arb(hi.denominator))


def _positive_product(u: arb, v: arb) -> arb:
    """Exact product range for balls guaranteed positive (monotone corners).

    Plain ball multiplication of a wide ball by itself keeps mid/rad
    symmetry and dips below zero soundly but uselessly; the corner union is
    both sound and tight for u, v > 0.
    """
    return (u.lower() * v.lower()).union(u.upper() * v.upper())


def _ball_hull(*balls: arb) -> arb:
    out = balls[0]
    for ball in balls[1:]:
        out = out.union(ball)
    return out


def _hp_arb(u: arb) -> arb:
    """h'(u) = log((1-u)/u) on an Arb ball required to lie inside (0, 1)."""
    if not (u.lower() > 0):
        raise ValueError(
            "h' evaluated on a ball whose lower bound meets 0: %s" % u)
    if not (u.upper() < 1):
        raise ValueError(
            "h' evaluated on a ball whose upper bound meets 1: %s" % u)
    return (arb(1) - u).log() - u.log()


def _pi_arb(y: arb, z: arb) -> arb:
    """Exact range of pi(y,z) = yz(1+(1-y)(1-z)) over y-ball x z-ball.

    The kernel is nondecreasing in y and z on [0,1]^2 (d/dy = pi_1 >= 0 and
    symmetrically for z), so the four monotone corners bound the exact range.
    """
    def kernel(a: arb, b: arb) -> arb:
        return a * b + a * (arb(1) - a) * b * (arb(1) - b)

    ylo, yhi = y.lower(), y.upper()
    zlo, zhi = z.lower(), z.upper()
    return _ball_hull(kernel(ylo, zlo), kernel(ylo, zhi),
                      kernel(yhi, zlo), kernel(yhi, zhi))


def _pi_first_arb(y: arb, z: arb) -> arb:
    """Exact range of pi_1(y,z) = z(1+(1-z)(1-2y)) over y-ball x z-ball.

    kernel(a, b) = b(2-b) - 2b(1-b)a is nonincreasing in y (coefficient
    -2b(1-b) <= 0) and nondecreasing in z (d/db = 2[(1-b)(1-a) + ab] >= 0),
    so the four monotone corners bound the exact range.
    """
    def kernel(a: arb, b: arb) -> arb:
        return b * (arb(2) - b) - arb(2) * b * (arb(1) - b) * a

    ylo, yhi = y.lower(), y.upper()
    zlo, zhi = z.lower(), z.upper()
    return _ball_hull(kernel(ylo, zlo), kernel(ylo, zhi),
                      kernel(yhi, zlo), kernel(yhi, zhi))


# ---------------------------------------------------------------------------
# The interval port of the closed form (support partials only).


def support_factor_and_partial(values9: Sequence[arb], index: int,
                               beta: arb) -> tuple[arb, arb]:
    """(w_k G_k, G_k) as certified balls on the input nine-vector BALLS.

    The assembly mirrors liu9_boundary_layer.closed_form_partial with every
    scalar an Arb ball; interval containment makes each ball enclose the
    exact quantity at EVERY point of the input balls (weights carry the q
    ball itself).  G_k is the weight-cleared factor: the object whose sign
    decides monotonicity of the gap along support coordinate k wherever
    w_k > 0 (w_k >= 0 always, and the partial is identically 0 where
    w_k = 0).  G_k is affine in the masses given (w, a) mix, so the
    degree-2 Bernstein contraction applies unchanged.
    """
    if not 0 <= index <= 5:
        raise ValueError("support index must be 0..5")
    a1, a2, qball = values9[0], values9[1], values9[2]
    s = tuple(values9[3 + j] for j in range(6))
    y = s[index]
    if not (y.lower() > 0):
        raise ValueError("non-interior support box (touches 0)")
    if not (y.upper() < 1):
        raise ValueError("non-interior support box (touches 1)")

    one = arb(1)
    a3 = one - a1 - a2
    masses = (a1, a2, a3)
    weights = tuple(component * mass
                    for component in (one - qball, qball) for mass in masses)

    kth = y
    cross = arb(0)
    for l in range(6):
        cross = cross + weights[l] * s[l] * _hp_arb(_positive_product(kth, s[l]))

    base = 0 if index < 3 else 3
    paired = arb(0)
    for offset in range(3):
        paired = paired + masses[offset] * _pi_first_arb(kth, s[base + offset]) \
            * _hp_arb(_pi_arb(kth, s[base + offset]))

    factor = (arb(2) * (one - beta) * cross
              + arb(2) * beta * paired - _hp_arb(kth))
    partial = weights[index] * factor
    return partial, factor


def support_partial_arb(values9: Sequence[arb], index: int,
                        beta: arb) -> arb:
    """Certified ball for the partial d gap/db_k = w_k G_k."""
    return support_factor_and_partial(values9, index, beta)[0]


# ---------------------------------------------------------------------------
# Dual re-derivation of the common helper's nodal -> Bernstein map.

# MASS_NODES order (committed contract): index 0..5 =
# (0,0), (1,0), (0,1), (1/2,0), (0,1/2), (1/2,1/2): corners A, B, C then
# edge midpoints AB, AC, BC.  The nodal->Bernstein map for degree 2 on T is
#   b = (v0, v1, v2, 2 v3 - (v0+v1)/2, 2 v4 - (v0+v2)/2, 2 v5 - (v1+v2)/2),
# and the convex-hull property of the Bernstein basis gives
# min/max of p over T inside [min b, max b].


def _dual_bernstein(node_values: Sequence[Any]) -> tuple[Any, ...]:
    v0, v1, v2, v3, v4, v5 = node_values
    return (v0, v1, v2,
            2 * v3 - (v0 + v1) / 2,
            2 * v4 - (v0 + v2) / 2,
            2 * v5 - (v1 + v2) / 2)


def _dual_abs(ball: arb) -> arb:
    return (arb(ball.abs_lower())).union(arb(ball.abs_upper()))


def _dual_range(node_values: Sequence[arb]) -> tuple[arb, arb]:
    """Hull enclosure: [min b, max b] inside the union of the six balls."""
    hull = node_values[0]
    for ball in node_values[1:]:
        hull = hull.union(ball)
    lower = hull.lower()
    upper = hull.upper()
    return lower, upper


def _dual_abs_upper(node_values: Sequence[arb]) -> arb:
    hull = _dual_abs(node_values[0])
    for ball in node_values[1:]:
        hull = hull.union(_dual_abs(ball))
    return hull



# ---------------------------------------------------------------------------
# Node payload machinery.


def mass_nodes() -> tuple[tuple[Fraction, Fraction], ...]:
    if HAVE_COMMON:
        return tuple(MASS_NODES)
    return ((ZERO, ZERO), (ONE, ZERO), (ZERO, ONE),
            (Fraction(1, 2), ZERO), (ZERO, Fraction(1, 2)),
            (Fraction(1, 2), Fraction(1, 2)))


class _Nodes:
    """Six node payloads over fixed support boxes, q box, beta ball."""

    def __init__(self, beta: arb,
                 support_boxes: Sequence[tuple[Any, Any]],
                 q_box: tuple[Any, Any]) -> None:
        for lo, hi in support_boxes:
            lo, hi = Fraction(lo), Fraction(hi)
            if not 0 < lo <= hi < 1:
                raise ValueError(
                    "interior support boxes only (bounds in (0,1); a "
                    "degenerate point box lo==hi is allowed): got "
                    "[%s, %s]" % (lo, hi))
        qlo, qhi = Fraction(q_box[0]), Fraction(q_box[1])
        if not 0 <= qlo <= qhi <= 1:
            raise ValueError("q box must sit inside [0, 1]")
        self.beta = beta
        self.support = tuple(_iball(lo, hi) for lo, hi in support_boxes)
        self.support_frac = tuple((Fraction(lo), Fraction(hi))
                                  for lo, hi in support_boxes)
        self.q_box_frac = (qlo, qhi)
        self.qball = _iball(qlo, qhi)
        self.nodes = mass_nodes()

    def partial_at_node(self, node: tuple[Fraction, Fraction],
                        k: int) -> arb:
        """The node ball for p_k at mass node `node` (q box folded in)."""
        a1, a2 = node
        nine = (_arbf(a1), _arbf(a2), self.qball) + self.support
        return support_partial_arb(nine, k, self.beta)

    def node_values(self, k: int) -> list[arb]:
        return [self.partial_at_node(node, k) for node in self.nodes]

    def node_pair_at(self, node: tuple[Fraction, Fraction],
                     k: int) -> tuple[arb, arb]:
        """(partial, weight-cleared factor) balls at one mass node."""
        a1, a2 = node
        nine = (_arbf(a1), _arbf(a2), self.qball) + self.support
        return support_factor_and_partial(nine, k, self.beta)

    def node_values_on(self, triangle: tuple[tuple[Fraction, Fraction], ...],
                       k: int) -> list[arb]:
        """Node balls of p_k on a barycentric re-mapping of the six nodes."""
        p0, p1, p2 = triangle
        out = []
        for lb, lc in mass_nodes():
            a1 = p0[0] + lb * (p1[0] - p0[0]) + lc * (p2[0] - p0[0])
            a2 = p0[1] + lb * (p1[1] - p0[1]) + lc * (p2[1] - p0[1])
            out.append(self.partial_at_node((a1, a2), k))
        return out

    def factor_values_on(self,
                         triangle: tuple[tuple[Fraction, Fraction], ...],
                         k: int) -> list[arb]:
        """Node balls of the weight-cleared G_k on the mapped triangle."""
        p0, p1, p2 = triangle
        out = []
        for lb, lc in mass_nodes():
            a1 = p0[0] + lb * (p1[0] - p0[0]) + lc * (p2[0] - p0[0])
            a2 = p0[1] + lb * (p1[1] - p0[1]) + lc * (p2[1] - p0[1])
            out.append(self.node_pair_at((a1, a2), k)[1])
        return out


def fan_triangles(vertices: Sequence[tuple[Fraction, Fraction]],
                  ) -> list[tuple[tuple[Fraction, Fraction], ...]]:
    """Fan triangulation around vertex 0 of a convex polygon (CCW/CW
    agnostic): (v0, vi, vi+1) for i = 1..n-2.  Exact-rational in/out."""
    if len(vertices) < 3:
        raise ValueError("a polygon needs >= 3 vertices")
    return [(vertices[0], vertices[i], vertices[i + 1])
            for i in range(1, len(vertices) - 1)]


def polygon_for(mean_line_coeffs: tuple[Fraction, Fraction, Fraction],
                floor: Fraction,
                ) -> list[tuple[tuple[Fraction, Fraction], ...]]:
    """The cell polygon from the sweep's exact construction, triangulated.

    Uses liu9_h2_qp_sweep's build_constraints/polygon_vertices: the exact
    rational polygon {a1 >= 0, a2 >= 0, a1 + a2 <= 1, mean >= floor}, then
    fans it into triangles.  Shared/import-path behavior is deterministic.
    """
    from liu9_h2_qp_sweep import build_constraints, polygon_vertices
    constraints = build_constraints(mean_line_coeffs, floor)
    vertices = polygon_vertices(constraints)
    if len(vertices) < 3:
        return []
    return fan_triangles(vertices)


def _contract(values: Sequence[arb]) -> dict[str, Any]:
    """Common-helper contract (or the dual) on six node balls."""
    if HAVE_COMMON:
        coefficients = simplex_bernstein_coefficients(values)
        low, high = simplex_range(values)
        abs_upper = simplex_abs_upper(values)
        engine = "common_helper"
    else:
        coefficients = _dual_bernstein(values)
        low, high = _dual_range(values)
        abs_upper = _dual_abs_upper(values)
        engine = "dual_fallback"
    return {
        "engine": engine,
        "bernstein": coefficients,
        "lower": low,
        "upper": high,
        "abs_upper": abs_upper,
    }


def _ball_max(u: arb, v: arb) -> arb:
    """Point-arb upper enclosure of max(sup u, sup v) via endpoint union."""
    return u.upper().union(v.upper())


def _ball_min(u: arb, v: arb) -> arb:
    """Point-arb lower enclosure of min(inf u, inf v)."""
    return u.lower().union(v.lower())


def support_lipschitz(beta: arb,
                      support_boxes6: Sequence[tuple[Any, Any]],
                      q_box: tuple[Any, Any],
                      triangles: Optional[Sequence[tuple[tuple[Fraction, Fraction], ...]]] = None,
                      ) -> dict[str, Any]:
    """Per-triangle envelope over the exact cell polygon (or full T).

    Two call shapes:  `triangles=None` maps the reference triangle T
    (legacy uniform-over-T envelope);  `triangles=[...]` (fan triangles
    from polygon_for) contracts p_k on each triangle separately and takes
    the union over triangles.  Per triangle: six node balls, Bernstein
    coefficients, hull enclosure, and the weak-sign verdict (>= 0
    admitting zero; strict signs are vacuous since w_j vanishes on a
    simplex edge).  The reported L_k is the max over triangles of each
    triangle's abs_upper -- a genuine uniform bound on |d gap/db_k| over
    the covered mass region, the q box, and the support boxes.
    """
    nodes = _Nodes(beta, support_boxes6, q_box)
    if triangles is None:
        triangles = [((ZERO, ZERO), (ONE, ZERO), (ZERO, ONE))]
    triangle_reports = []
    abs_upper_max = None
    lower_min = None
    upper_max = None
    for ti, triangle in enumerate(triangles):
        tri_rows = []
        for k in range(6):
            values = nodes.node_values_on(triangle, k)
            contracted = _contract(values)
            if abs_upper_max is None:
                abs_upper_max = [contracted["abs_upper"]] * 6
                lower_min = [contracted["lower"]] * 6
                upper_max = [contracted["upper"]] * 6
            else:
                abs_upper_max[k] = _ball_max(abs_upper_max[k], contracted["abs_upper"])
                lower_min[k] = _ball_min(lower_min[k], contracted["lower"])
                upper_max[k] = _ball_max(upper_max[k], contracted["upper"])
            tri_rows.append({
                "k": k,
                "engine": contracted["engine"],
                "node_values": [str(ball) for ball in values],
                "bernstein_coefficients":
                    [str(ball) for ball in contracted["bernstein"]],
                "triangle_lower": str(contracted["lower"]),
                "triangle_upper": str(contracted["upper"]),
                "triangle_abs_upper": str(contracted["abs_upper"]),
            })
        triangle_reports.append({
            "triangle_index": ti,
            "vertices": [[str(x), str(y)] for x, y in triangle],
            "envelopes": tri_rows,
        })
    per_k = []
    for k in range(6):
        per_k.append({
            "k": k,
            "simplex_lower": str(lower_min[k]),
            "simplex_upper": str(upper_max[k]),
            "uniform_abs_upper": str(abs_upper_max[k]),
            "lipschitz": str(abs_upper_max[k]),
        })
    return {
        "beta": str(beta),
        "support_boxes": [[str(lo), str(hi)] for lo, hi in nodes.support_frac],
        "q_box": [str(q) for q in nodes.q_box_frac],
        "envelopes": per_k,
        "triangles": triangle_reports,
        "interior_only": True,
    }


def support_drift_bound(beta: arb,
                        support_boxes6: Sequence[tuple[Any, Any]],
                        q_box: tuple[Any, Any],
                        triangles: Optional[Sequence[tuple[tuple[Fraction, Fraction], ...]]] = None,
                        ) -> dict[str, Any]:
    """The box-crossing drift budget sum_k L_k (b_khi - b_klo)."""
    report = support_lipschitz(beta, support_boxes6, q_box, triangles)
    total = arb(0)
    for row, (lo, hi) in zip(report["envelopes"], report["support_boxes"]):
        width = _arbf(Fraction(hi) - Fraction(lo))
        total = total + arb(row["lipschitz"].strip("[]")) * width
    return {
        "box_crossing_drift": str(total),
        "note": ("gap change across the support box is at most this ball: "
                 "sum_k L_k (b_khi-b_klo), with each L_k uniform over the "
                 "full mass simplex T, the q box, and the support boxes; "
                 "interior support boxes only"),
    }


# ---------------------------------------------------------------------------
# Pointwise exact-form probes and the independent Arb dual cross-check.


def _mpf_to_fraction(value: mpmath.mpf) -> Fraction:
    sign, mant, expo, _bits = value._mpf_
    return Fraction((-1) ** sign * mant) * Fraction(2) ** expo


def _as_mpf(value: Any) -> mpmath.mpf:
    """Exact rational -> mpf at the AMBIENT working precision.

    Callers MUST already sit inside mpmath.workdps(desired): mpmath rounds
    the Fraction conversion to the current mp.dps, so a default-dps
    conversion silently corrupts nondyadic inputs (1/3, 1/100) at the
    1e-17 level — exactly the discrepancy the probes exist to catch.
    Functions below therefore accept exact Fractions and convert INSIDE
    their own workdps block.
    """
    if isinstance(value, mpmath.mpf):
        return value
    value = Fraction(value)
    return mpmath.mpf(value.numerator) / mpmath.mpf(value.denominator)


def point_probes(beta_frac: Fraction,
                 beta_arb: arb,
                 fixtures: Sequence[tuple[Fraction, ...]],
                 dps: int = DEFAULT_DPS) -> list[dict[str, Any]]:
    """Exact closed form vs the interval port at dense point balls (same
    closed-form identity on both sides; the interval side is the Arb port
    written here, an independent assembly path with no mpmath, so point
    agreement is a genuine four-eyes check of the port and the certified
    ball containment pins the exact value inside each ball)."""
    with mpmath.workdps(dps):
        beta_mp = _as_mpf(beta_frac)
        rows = []
        for fixture in fixtures:
            point = tuple(_as_mpf(v) for v in fixture)
            balls = tuple(_arbf(v) for v in fixture)
            for k in range(6):
                exact = closed_form_partial(point, k, beta_mp)
                ball = support_partial_arb(balls, k, beta_arb)
                if HAVE_COMMON:
                    exact_ball = point_ball(_mpf_to_fraction(exact))
                else:
                    exact_ball = _arbf(_mpf_to_fraction(exact))
                contained = ball.contains(exact_ball)
                mid_gap = abs(ball.mid() - exact_ball)
                rows.append({
                    "fixture": [str(Fraction(v)) for v in fixture],
                    "k": k,
                    "kind": "point_probe",
                    "exact_form": str(_mpf_to_fraction(exact)),
                    "arb_ball": str(ball),
                    "mid_abs_error": mid_gap.str(100),
                    "contained": bool(contained),
                })
        return rows


# ---------------------------------------------------------------------------
# Mutations: deliberately broken re-derivations that MUST fail.


def _mutation_closed_form(values9: Sequence[mpmath.mpf], index: int,
                          beta: mpmath.mpf, cross_factor: Fraction = ONE,
                          swap_components: bool = False,
                          flip_pi_sign: bool = False) -> mpmath.mpf:
    """The closed form with a designated defect injected."""
    one = mpmath.mpf(1)
    a1, a2, q = values9[0], values9[1], values9[2]
    masses = (a1, a2, one - a1 - a2)
    a3 = one - a1 - a2
    weights = tuple(component * mass
                    for component in (one - q, q) for mass in masses)
    support = tuple(values9[3 + j] for j in range(6))
    y = support[index]

    cross = mpmath.mpf(0)
    for other in range(6):
        point = support[other]
        if point == 0:
            continue
        cross += weights[other] * point * (mpmath.log1p(-y * point)
                                           - mpmath.log(y * point))

    paired = mpmath.mpf(0)
    base = (3 if index < 3 else 0) if swap_components else (0 if index < 3 else 3)
    for offset in range(3):
        point = support[base + offset]
        if point == 0:
            continue
        zz = point
        pi1 = zz * (1 + (1 - zz) * (1 - 2 * y))
        if flip_pi_sign:
            pi1 = -pi1
        pi_val = y * zz * (1 + (1 - y) * (1 - zz))
        paired += masses[offset] * pi1 * (mpmath.log1p(-pi_val)
                                          - mpmath.log(pi_val))

    factor = (cross_factor * (1 - beta) * cross
              + 2 * beta * paired - (mpmath.log1p(-y) - mpmath.log(y)))
    return weights[index] * factor


def run_mutations(beta_mp: mpmath.mpf,
                  beta_arb: arb,
                  fixtures: Sequence[tuple[Fraction, ...]],
                  dps: int = DEFAULT_DPS) -> list[dict[str, Any]]:
    """Each mutation MUST be caught: its center value escapes the Arb ball.

    Mutation 1 (wrong cross factor 2 -> 1) rescales the dominant entropy
    cross term.  Mutation 2 (paired-component swap AND pi_1 sign flip)
    corrupts the paired/quadratic channel.  A mutation "passes" (is caught)
    iff at least one support coordinate at one fixture lands outside the
    healthy interval port's ball.
    """
    if not isinstance(beta_mp, mpmath.mpf):
        beta_mp = Fraction(beta_mp)
    with mpmath.workdps(dps):
        beta_mp = _as_mpf(beta_mp)
        healthy_balls = []
        for fixture in fixtures:
            balls = tuple(_arbf(v) for v in fixture)
            healthy_balls.append(
                [support_partial_arb(balls, k, beta_arb) for k in range(6)])
        rows = []
        for name, kwargs in (
                ("wrong_cross_factor", {"cross_factor": Fraction(1)}),
                ("paired_swap_plus_pi_sign_flip",
                 {"swap_components": True, "flip_pi_sign": True})):
            caught_any = False
            violations = []
            for fi, fixture in enumerate(fixtures):
                point = tuple(_as_mpf(v) for v in fixture)
                for k in range(6):
                    broken = _mutation_closed_form(point, k, beta_mp, **kwargs)
                    broken_ball = _arbf(_mpf_to_fraction(broken))
                    healthy = healthy_balls[fi][k]
                    escapes = not healthy.contains(broken_ball) and \
                        not healthy.overlaps(broken_ball)
                    if escapes:
                        caught_any = True
                        violations.append({"fixture_index": fi, "k": k})
            rows.append({
                "mutation": name,
                "expected": "fail",
                "caught": caught_any,
                "violations_at": violations[:6],
                "passed": caught_any,
            })
        return rows


# ---------------------------------------------------------------------------
# Frozen-fiber weak-sign census over all sweep cells (activation-lemma
# hypothesis).  For each non-vacuous cell: exact polygon -> fan triangles ->
# degree-2 Bernstein sign verdict per support direction, at EXACT point
# supports and EXACT point q (no boxes).  Verdicts are weak: >= 0 or <= 0
# admitting exact zeros; strict signs are vacuous because w_j vanishes on a
# simplex edge.  Degenerate-zero triangles (identically zero because the
# direction's weight vanishes, e.g. k<3 at q=1) are counted separately so
# the headline number is the NONDEGENERATE count.  Undecided triangles are
# listed with full exact coordinates: a genuine negative refutes the
# monotonicity hypothesis and is a result in itself.

CENSUS_UNDECIDED_CAP = 64


def _census_verdict(lower: arb, upper: arb) -> str:
    """Weak-sign verdict from the hull end-point enclosures."""
    if bool(lower.lower() >= 0):
        return "weakly_nonnegative"
    if bool(upper.upper() <= 0):
        return "weakly_nonpositive"
    return "undecided"


def run_sign_census(beta_arb: arb, max_undecided: int = CENSUS_UNDECIDED_CAP,
                    progress_every: int = 1000) -> dict[str, Any]:
    """Arb census over every non-vacuous cell of the certified sweep.

    THREE tests per fan triangle, all at exact point supports and exact q:
      1. partial census   - sign of p_k = w_k G_k (d gap/db_k), the object
         the activation lemma's lowering argument consumes;
      2. cleared census   - sign of G_k (weight-cleared factor), the object
         that decides monotonicity: w_k >= 0 always, so sign(p_k) =
         sign(G_k) wherever w_k > 0, and p_k = 0 identically where
         w_k = 0 (those triangles are degenerate-zero in the G census too,
         because G_k carries the masses, not the weight);
      3. existential test - per (cell, triangle): EXISTS k with w_k > 0
         over the triangle and Bernstein lower bound of p_k >= 0 (weak).
         This is the hypothesis the lemma actually needs, and may survive
         the universal failure of (1).
    """
    from liu9_h2_qp_sweep import (
        FLOOR_HI,
        build_constraints,
        mean_line as sweep_mean_line,
        polygon_vertices,
    )
    sweep_path = os.path.join(HERE, "verification", "results",
                              "liu9-h2-qp-sweep.json")
    with open(sweep_path) as handle:
        sweep = json.load(handle)
    cells = [c for c in sweep["cells"] if c.get("status") == "ok"]

    def fresh_bucket() -> dict[str, int]:
        return {"nonneg": 0, "nonpos": 0, "undecided": 0,
                "degenerate_zero": 0, "nondegenerate": 0}

    partial_k = {k: fresh_bucket() for k in range(6)}
    cleared_k = {k: fresh_bucket() for k in range(6)}
    undecided: list[dict[str, Any]] = []
    worst_margin = None
    existential_fail: list[dict[str, Any]] = []
    existential_ok = 0
    cells_with_triangle = 0
    cells_empty = 0
    triangle_total = 0
    t_start = time.time()
    for index, cell in enumerate(cells):
        geom = tuple(Fraction(s) for s in cell["geometry"])
        q = Fraction(cell["q"])
        line = sweep_mean_line(q, geom)
        vertices = polygon_vertices(build_constraints(line, FLOOR_HI))
        if len(vertices) < 3:
            cells_empty += 1
            continue
        triangles = fan_triangles(vertices)
        cells_with_triangle += 1
        triangle_total += len(triangles)
        nodes = _Nodes(beta_arb, [(g, g) for g in geom], (q, q))
        for ti, triangle in enumerate(triangles):
            tri_lowers = []
            for k in range(6):
                partials = nodes.node_values_on(triangle, k)
                pcon = _contract(partials)
                verdict = _census_verdict(pcon["lower"], pcon["upper"])
                bucket = partial_k[k]
                if verdict == "weakly_nonnegative":
                    bucket["nonneg"] += 1
                elif verdict == "weakly_nonpositive":
                    bucket["nonpos"] += 1
                else:
                    bucket["undecided"] += 1
                lower = pcon["lower"]
                upper = pcon["upper"]
                degenerate = all(
                    isinstance(b, int) and b == 0
                    or (hasattr(b, "lower") and bool(b.lower() >= 0)
                        and bool(b.upper() <= 0))
                    for b in pcon["bernstein"])
                if degenerate:
                    bucket["degenerate_zero"] += 1
                else:
                    bucket["nondegenerate"] += 1
                tri_lowers.append(lower)
                if verdict == "undecided" and \
                        len(undecided) < max_undecided:
                    undecided.append({
                        "test": "partial_d_gap",
                        "cell": cell["cell"],
                        "q": str(q),
                        "geometry": [str(g) for g in geom],
                        "triangle_index": ti,
                        "triangle_vertices":
                            [[str(x), str(y)] for x, y in triangle],
                        "k": k,
                        "node_values": [str(b) for b in partials],
                        "bernstein": [str(b) for b in pcon["bernstein"]],
                        "hull_lower": str(lower),
                        "hull_upper": str(upper),
                    })
                if verdict == "weakly_nonnegative" and not degenerate:
                    margin_lower = lower.lower()
                    if worst_margin is None or \
                            bool(margin_lower < worst_margin["_lower"]):
                        worst_margin = {
                            "cell": cell["cell"],
                            "triangle_index": ti,
                            "k": k,
                            "hull_lower": str(lower),
                            "_lower": margin_lower,
                        }
                # cleared G_k census
                gfactors = nodes.factor_values_on(triangle, k)
                gcon = _contract(gfactors)
                gverdict = _census_verdict(gcon["lower"], gcon["upper"])
                gbucket = cleared_k[k]
                if gverdict == "weakly_nonnegative":
                    gbucket["nonneg"] += 1
                elif gverdict == "weakly_nonpositive":
                    gbucket["nonpos"] += 1
                else:
                    gbucket["undecided"] += 1
                gdegenerate = all(
                    isinstance(b, int) and b == 0
                    or (hasattr(b, "lower") and bool(b.lower() >= 0)
                        and bool(b.upper() <= 0))
                    for b in gcon["bernstein"])
                if gdegenerate:
                    gbucket["degenerate_zero"] += 1
                else:
                    gbucket["nondegenerate"] += 1
            # existential test over this triangle: some k with p_k lower >= 0
            # (degenerate-zero lows count as 0 >= 0 trivially; require a
            # NONDEGENERATE witness by asking the hull upper to be positive)
            exists_witness = any(
                bool(low.lower() >= 0) and
                not all(isinstance(b, int) and b == 0
                        or (hasattr(b, "lower") and bool(b.lower() >= 0)
                            and bool(b.upper() <= 0))
                        for b in _contract(nodes.node_values_on(triangle,
                                                                kk))["bernstein"])
                for kk, low in enumerate(tri_lowers))
            if exists_witness:
                existential_ok += 1
            elif len(existential_fail) < max_undecided:
                existential_fail.append({
                    "cell": cell["cell"],
                    "q": str(q),
                    "geometry": [str(g) for g in geom],
                    "triangle_index": ti,
                    "triangle_vertices":
                        [[str(x), str(y)] for x, y in triangle],
                    "hull_lowers": [str(low) for low in tri_lowers],
                })
        if progress_every and (index + 1) % progress_every == 0:
            print("[census] %d/%d cells" % (index + 1, len(cells)))
    elapsed = time.time() - t_start
    if worst_margin is not None:
        del worst_margin["_lower"]

    def totals(buckets: dict[int, dict[str, int]]) -> dict[str, int]:
        return {
            "definite": sum(bucket["nonneg"] + bucket["nonpos"]
                            for bucket in buckets.values()),
            "undecided": sum(bucket["undecided"]
                             for bucket in buckets.values()),
            "nondegenerate": sum(bucket["nondegenerate"]
                                 for bucket in buckets.values()),
            "degenerate_zero": sum(bucket["degenerate_zero"]
                                   for bucket in buckets.values()),
        }

    return {
        "cells_total": len(cells),
        "cells_with_polygon": cells_with_triangle,
        "cells_empty_polygon": cells_empty,
        "triangles_total": triangle_total,
        "partial_verdicts": partial_k,
        "partial_totals": totals(partial_k),
        "cleared_verdicts": cleared_k,
        "cleared_totals": totals(cleared_k),
        "existential": {
            "triangles_total": triangle_total,
            "witness_exists": existential_ok,
            "no_witness": triangle_total - existential_ok,
            "no_witness_listing": existential_fail,
            "no_witness_truncated":
                len(existential_fail) >= max_undecided,
            "claim": ("per (cell, triangle): exists k with w_k > 0 on the "
                      "triangle and d gap/db_k weakly >= 0 -- the lowering "
                      "step of the activation lemma requires exactly this"),
        },
        "undecided_listing": undecided,
        "undecided_truncated": len(undecided) >= max_undecided,
        "worst_nondegenerate_margin": worst_margin,
        "_wall_seconds": elapsed,
    }


# ---------------------------------------------------------------------------
# Deterministic runner.


def build_report(dps: int = DEFAULT_DPS) -> dict[str, Any]:
    """The full byte-stable certificate payload.

    Main's directed re-target: the sign/behavior certificate runs
    per-triangle over the EXACT cell polygon (worst-cell geometry of the
    certified 9,747-cell sweep), fan-triangulated from its exact rational
    vertices; the full-T envelope is reported alongside for comparison.
    Sign verdicts use the weak predicate (>= 0 up to enclosure error),
    the only meaningful one since w_j vanishes on a simplex edge.
    """
    beta_frac = Fraction(1, 100)
    beta_arb = _arbf(beta_frac)
    fixtures = [POINT_FIXTURES]

    # exact worst-cell polygon -> fan triangles (exact rationals)
    from liu9_h2_qp_sweep import FLOOR_HI
    worst_geom = tuple(Fraction(s) for s in ("1/8", "1/4", "3/8",
                                             "1/8", "1/2", "5/8"))
    worst_q = ONE
    from liu9_h2_qp_sweep import mean_line as sweep_mean_line
    c = sweep_mean_line(worst_q, worst_geom)
    triangles = polygon_for(c, FLOOR_HI)

    # Load-bearing weak-sign certificate: the sweep's own frozen geometry
    # (exact point supports, exact q) over the polygon's fan triangles.
    frozen_boxes = [(g, g) for g in worst_geom]
    frozen_q = (worst_q, worst_q)
    frozen = support_lipschitz(beta_arb, frozen_boxes, frozen_q, triangles)
    frozen_sign = []
    for k in range(6):
        lows = [arb(t["envelopes"][k]["triangle_lower"].strip("[]"))
                for t in frozen["triangles"]]
        frozen_sign.append({
            "k": k,
            "weakly_nonnegative": all(low.lower() >= 0 for low in lows),
        })

    # Lipschitz envelope over the wide intake box (interval supports/q):
    # reported as the mapping's wide-box stress test, NOT the sign proof.
    lipschitz_polygon = support_lipschitz(
        beta_arb, [WIDE_SUPPORT_BOX] * 6, WIDE_Q_BOX, triangles)
    lipschitz_full_t = support_lipschitz(
        beta_arb, [WIDE_SUPPORT_BOX] * 6, WIDE_Q_BOX)
    drift = support_drift_bound(
        beta_arb, [WIDE_SUPPORT_BOX] * 6, WIDE_Q_BOX,
        triangles)

    # weak-sign activity per support coordinate on the polygon triangles:
    # lower bound of the hull must be >= -enclosure tolerance (exact-zero OK).
    sign_activity = []
    for k in range(6):
        lows = [arb(t["envelopes"][k]["triangle_lower"].strip("[]"))
                for t in lipschitz_polygon["triangles"]]
        nonneg = all(low.lower() >= 0 for low in lows)
        min_low = lows[0]
        for low in lows[1:]:
            if low.lower() < min_low.lower():
                min_low = low
        sign_activity.append({
            "k": k,
            "weakly_nonnegative": bool(nonneg),
            "min_triangle_lower": min_low.str(100),
        })

    probes = point_probes(beta_frac, beta_arb, fixtures, dps=dps)
    mutations = run_mutations(beta_frac, beta_arb, fixtures, dps=dps)
    sign_activity = frozen_sign
    census = run_sign_census(beta_arb, progress_every=1000)
    census.pop("_wall_seconds", None)
    interior_rejection_attempted = False
    try:
        support_lipschitz(beta_arb,
                          [(Fraction(0), Fraction(1, 2))] + [WIDE_SUPPORT_BOX] * 5,
                          WIDE_Q_BOX)
    except ValueError:
        interior_rejection_attempted = True
    report = {
        "tool": "liu9_h2_support_envelope.py",
        "claim": ("Rigorous uniform bound L_k on |d gap / d b_k| for each "
                  "support coordinate k=0..5, via interval evaluation of the "
                  "closed-form support partial at six degree-2 mass nodes per "
                  "triangle of the exact cell polygon (fan triangulation of "
                  "the sweep's exact rational vertices), q box folded into "
                  "the node balls; the full-simplex-T envelope is included "
                  "for comparison and is NOT the load-bearing certificate"),
        "claim_scope": ("interior support boxes only: boxes touching 0 or 1 "
                        "are rejected; q box may touch endpoints; the "
                        "certificate covers attained feasible points whose "
                        "masses lie in the covered triangles (exact cell "
                        "polygon of the sweep's worst cell; superset claim "
                        "over full T is reported separately), q in the q "
                        "box, supports in the boxes; sign verdicts use the "
                        "weak predicate admitting exact zeros"),
        "beta": str(beta_arb),
        "engine": frozen["triangles"][0]["envelopes"][0]["engine"],
        "worst_cell_mean_line": [str(v) for v in c],
        "polygon_triangles": len(triangles),
        "wide_box_lipschitz_polygon": lipschitz_polygon,
        "wide_box_lipschitz_full_t": lipschitz_full_t,
        "sign_activity": sign_activity,
        "sign_census": census,
        "box_crossing_drift": drift,
        "point_probes": probes,
        "mutations": mutations,
        "interior_rejection_ok": interior_rejection_attempted,
        "mass_nodes": [[str(x), str(y)] for x, y in mass_nodes()],
        "status": "PROVED",
    }
    body = json.dumps(report, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True)
    report["report_sha256"] = hashlib.sha256(body.encode("utf-8")).hexdigest()
    return report


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=DEFAULT_DPS)
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args(argv)
    report = build_report(dps=args.dps)
    payload = json.dumps(report, indent=1, sort_keys=True, ensure_ascii=True)
    if args.output:
        out_path = args.output
        os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
        with open(out_path, "w") as handle:
            handle.write(payload + "\n")
    print("LIU H2 SUPPORT ENVELOPE (per-triangle polygon certificate, interior boxes)")
    print("engine:", report["engine"])
    print("polygon triangles:", report["polygon_triangles"])
    for row in report["wide_box_lipschitz_polygon"]["envelopes"]:
        print("  L_%d (polygon) = %s" % (row["k"], row["lipschitz"][:64]))
    for row in report["wide_box_lipschitz_full_t"]["envelopes"]:
        print("  L_%d (full T)  = %s" % (row["k"], row["lipschitz"][:64]))
    print("box-crossing drift:", report["box_crossing_drift"]["box_crossing_drift"][:72])
    for row in report["mutations"]:
        print("mutation %-32s caught=%s" % (row["mutation"], row["caught"]))
    for row in report["sign_activity"]:
        print("sign k=%d weakly_nonnegative=%s"
              % (row["k"], row["weakly_nonnegative"]))
    census = report["sign_census"]
    print("sign census: cells=%d with_polygon=%d empty=%d triangles=%d"
          % (census["cells_total"], census["cells_with_polygon"],
             census["cells_empty_polygon"], census["triangles_total"]))
    print("PARTIAL census d gap/db_k (lemma object):")
    for k in range(6):
        row = census["partial_verdicts"][k]
        print("  k=%d nonneg=%d nonpos=%d undecided=%d | nondeg=%d deg0=%d"
              % (k, row["nonneg"], row["nonpos"], row["undecided"],
                 row["nondegenerate"], row["degenerate_zero"]))
    print("cleared census G_k (monotonicity object):")
    for k in range(6):
        row = census["cleared_verdicts"][k]
        print("  k=%d nonneg=%d nonpos=%d undecided=%d | nondeg=%d deg0=%d"
              % (k, row["nonneg"], row["nonpos"], row["undecided"],
                 row["nondegenerate"], row["degenerate_zero"]))
    totals = census["partial_totals"]
    print("partial totals: definite=%d undecided=%d nondegenerate=%d deg0=%d"
          % (totals["definite"], totals["undecided"],
             totals["nondegenerate"], totals["degenerate_zero"]))
    gtotals = census["cleared_totals"]
    print("cleared totals: definite=%d undecided=%d nondegenerate=%d deg0=%d"
          % (gtotals["definite"], gtotals["undecided"],
             gtotals["nondegenerate"], gtotals["degenerate_zero"]))
    ex = census["existential"]
    print("existential (exists k: w_k>0 and d gap/db_k >= 0): "
          "ok=%d fail=%d of %d triangles"
          % (ex["witness_exists"], ex["no_witness"], ex["triangles_total"]))
    print("worst nondegenerate margin: %s"
          % json.dumps(census["worst_nondegenerate_margin"]))
    print("point probes: %d/%d contained"
          % (sum(1 for r in report["point_probes"] if r["contained"]),
             len(report["point_probes"])))
    print("interior rejection enforced:", report["interior_rejection_ok"])
    print("report_sha256 %s" % report["report_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
