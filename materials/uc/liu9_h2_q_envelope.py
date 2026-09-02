#!/usr/bin/env python3
"""Exact q-envelope for Liu's paired r = 3 raw gap (Hypothesis 2, q direction).

WHAT THIS MODULE IS (and is not)
--------------------------------
The live scalar wall needs the inward-q statement: for a frozen support box
(the x-triple and y-triple each in an exact box) and a q box inside [0, 1],
a uniform-over-the-WHOLE-mass-simplex statement built from exact quantity
certificates only.  The pieces, each a separately labelled CLAIM (every
claim is refuted by at least one of the two mutations in run_mutations()):

  CLAIM (Q1)  the raw gap is EXACTLY quadratic in q at every frozen
              (supports, masses) fiber: the (A, B, C) extraction at the
              exact rational q values 0, 1/2, 1 (exact-rational inverse of
              the 3-point Vandermonde, applied inclusion-monotonically to
              certified value balls), validated by residual containment at
              further exact rational q probes and at a second mass fiber;
  CLAIM (Q2)  the q partial d gap / d q is an EXACT arithmetic expression
              obtained by dual-number propagation through the repository's
              single shared transcription liu9_objective._formula (the
              liu9_boundary_layer autodiff pattern, over Arb balls).  No
              finite difference participates in any claim path anywhere;
  CLAIM (Q3)  q_lipschitz: a certified uniform bound
                 max { |d gap/d q(s, q, a)| : s in the support boxes,
                       q in the q box, a in the full mass simplex T }
              from six certified node balls via
              liu9_h2_envelope_common's degree-2 nodal-to-Bernstein map
              (the convex-hull range property on T; inclusion monotone);
  CLAIM (Q4)  q_drift_bound: the uniform drift
                 |gap(q', a, s) - gap(q, a, s)| <= (q' - q) * L
              for all q, q' in the q box and ALL a in T, s in the support
              boxes, by the fundamental theorem on the (C^1, polynomial-
              in-q) gap with the (Q3) constant L;
  CLAIM (Q5)  q_point_minimization: the exact q-quadratic point
              minimization with endpoints plus interior critical point
              over the q box (exact rational algebra on the coefficient
              balls, the certified candidates only), for later
              branch-and-bound composition.

NO FIXED-ARGMIN RADIUS SHORTCUT.  The retracted uc/liu9_h2_box_cover.py
pattern (value-bracket radius around ONE fixed KKT argmin mass point;
SOUNDNESS_ARGUMENT step (vi)) is NOT reused anywhere: nothing here bounds
mass-direction behavior at a single mass.  Mass uniformity comes ONLY from
the Bernstein convex-hull property over the six nodes spanning T (Q3), and
the drift term (Q4) is by construction uniform over T.

PROBES (never claims): the finite-difference agreement check
(exact for polynomials of degree <= 2, so it is simultaneously a strong
consistency probe of the dual engine and a corroborating degree-2
check), and the independent transcription-walker containment check.
Each is labelled PROBE in the report and in the code.

UNIVERSAL ENDPOINT THEOREM (cited, NOT reproved).
  uc/liu9_size_biased.py, certify_qone_kernel, claim_status PROVED:
  F(P) >= 0 for EVERY probability law P with mean >= m, via
  F(P) = M E_{Q x Q} C_M with Q(dx) = x P(dx)/M, C_m >= 0 pointwise and
  dC_M/dM >= 0.  The exact split gap is
      gap = F(P_mix) + q(1-q) c_ch        (proved, r <= 3),
  the channel term VANISHES identically at q in {0, 1}, so gap(q=1) = F(P1)
  and gap(q=0) = F(P0); the mixture at an endpoint fiber IS that component
  law and the mixture mean IS its mean, so the cited theorem discharges
  BOTH endpoint q faces as boundary data.  No sign is claimed for c_ch:
  the certified fact is the domination ratio
      q(1-q)|c| / margin <= 0.838883838221919
  on the 9,747 solved sweep cells (zero danger cells), and the certified
  most-negative channel scalar is c = -0.02755 (attained at a mass vertex,
  on all 1,225 geometries) -- asserting "c_ch >= 0" would be FALSE.
  Endpoint q boxes are one-sided by construction (clipped to [0, 1]); the
  derivative/Lipschitz machinery of this module is valid on the CLOSED
  endpoint boxes as well (the q variable enters the transcription only
  through weight factors, never inside an entropy argument).

INTERIOR DISCIPLINE.  Support boxes touching 0 or 1 are REJECTED
explicitly (ValueError) by every derivative enclosure here: the chain-rule
h' enclosure needs a support argument ball strictly inside (0, 1).  The q
box may touch 0 or 1 (one-sided endpoint handling): the q variable never
enters an entropy argument.

BYTE STABILITY.  No randomness, no timestamps; every report quantity
derives from exact rationals and certified Arb balls at 320 bits.  Two
consecutive runs write byte-identical reports; the report carries the
canonical sha256 of its own body.
"""

from __future__ import annotations

import os
import sys

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_name, "1")

import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Optional, Sequence

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from flint import arb, ctx                     # noqa: E402
ctx.prec = max(ctx.prec, 320)

from liu9_boundary_layer import _GapOps, _formula  # noqa: E402
from liu9_h2_envelope_common import (          # noqa: E402
    MASS_NODES, point_ball, simplex_abs_upper,
    simplex_bernstein_coefficients, simplex_range,
)
from liu9_objective import h_arb               # noqa: E402

OUTPUT_DEFAULT = HERE / "verification/results/liu9-h2-q-envelope.json"

ONE = Fraction(1)
ZERO = Fraction(0)

# Demo fiber (deterministic, byte-stable): k/8-grid x-triple (1/8, 1/4, 3/8)
# and y-triple (3/8, 1/2, 3/4); strictly interior masses (19/20, 1/20);
# q box [1/8, 1/4] (interior, and it contains the fiber's certified interior
# critical point ~0.2135 so the interior branch of Q5 fires); support
# half-width 1/32 keeps every support box strictly inside (0, 1).
DEMO_XTRIPLE = (Fraction(1, 8), Fraction(1, 4), Fraction(3, 8))
DEMO_YTRIPLE = (Fraction(3, 8), Fraction(1, 2), Fraction(3, 4))
DEMO_MASSES = (Fraction(19, 20), Fraction(1, 20))
DEMO_Q_BOX = (Fraction(1, 8), Fraction(1, 4))
DEMO_WH = Fraction(1, 32)

# Probe / mutation constants (exact rationals).
FD_STEP = Fraction(1, 64)          # symmetric FD probe step (exact for deg<=2)
SECOND_FIBER = (Fraction(1, 3), Fraction(1, 3))
Q_PROBES = (Fraction(1, 3), Fraction(2, 3))
MUT_COEF_DELTA = Fraction(1, 10 ** 6)


# ---------------------------------------------------------------------------
# Certified parameter ball (never a hardcoded decimal in a claim path).

def certified_beta_ball() -> arb:
    """The equation-defined beta ball from liu9_binding's certified chain,
    re-derived on every run (deterministic at the pinned working precision);
    plus an exact rational two-endpoint enclosure recorded next to it."""
    from liu9_binding import (certify_equation_parameters,
                              solve_equation_parameters)
    params = certify_equation_parameters(solve_equation_parameters(70))
    return params.beta


def exact_enclosure(beta: arb) -> tuple[Fraction, Fraction]:
    """Exact rational endpoints strictly OUTSIDE the certified ball, so the
    hull point balls contain it.  The endpoint PARSING is exact-digit round
    to nearest with an explicit OUTWARD slack of 2^-100 on each side: the
    lower endpoint moves down, the upper up, and the containment is then
    asserted (raising on any parse/rounding surprise)."""
    def parsed(endpoint: arb) -> Fraction:
        return Fraction(endpoint.str(32).split("+/-")[0].strip("[]").strip())

    slack = Fraction(1, 2 ** 100)
    lo = parsed(beta.lower()) - slack
    hi = parsed(beta.upper()) + slack
    enclosure = point_ball(lo).union(point_ball(hi))
    if not (bool(enclosure.lower() <= beta.lower())
            and bool(beta.upper() <= enclosure.upper())):
        raise ValueError("rational enclosure fails to contain the beta ball")
    return lo, hi


# ---------------------------------------------------------------------------
# The certified raw gap on balls.

class _GapOpsArb(_GapOps):
    """No-quotient ops with ENTRY entropy over interval balls:
    liu9_objective.h_arb is the exact range enclosure of the binary entropy
    over an Arb ball (its own monotone-corner / DIC lemma), so the gap
    evaluates on BALL arguments; every other op is standard interval *."""

    def __init__(self) -> None:
        super().__init__(h_arb, arb(1))


def gap_arb(values9: Sequence[arb], beta: arb) -> arb:
    """Certified two-sided enclosure of the raw gap
       gap = (1-beta) EHXY + beta EHPI - EHX
    over the ball inputs (gap = numerator - ehx, defined even where the
    objective N/D is not).  No support/q ball may touch an entropy-domain
    endpoint unless it stays in [0, 1] (h_arb clips semantically); the
    derivative enclosures below additionally require strict interior."""
    if len(values9) != 9:
        raise ValueError("expected nine variables (a1 a2 q b0 b2 b4 b1 b3 b5)")
    if not all(value.is_finite() for value in values9) or not beta.is_finite():
        raise ValueError("non-finite Arb input")
    a1, a2 = values9[0], values9[1]
    # Mass coordinates are point balls of EXACT simplex rationals (this
    # module's calling convention); their ball sums may protrude beyond 1
    # by enclosure roundoff only, so the semantic simplex test allows a
    # 2^-100 enclosure slack while the exact rationals satisfy the
    # constraint exactly.
    slack = point_ball(Fraction(1, 2 ** 100))
    if not (bool(a1.lower() >= -slack) and bool(a2.lower() >= -slack)
            and bool((a1 + a2).upper() <= 1 + slack)):
        raise ValueError("(a1, a2, 1-a1-a2) must lie in the mass simplex")
    for value in values9[2:]:
        if not (bool(value.lower() >= 0) and bool(value.upper() <= 1)):
            raise ValueError("q and supports must lie in [0, 1]")
    terms = _formula(tuple(values9), beta, _GapOpsArb())
    gap = terms.numerator - terms.ehx
    if not gap.is_finite():
        raise ValueError("non-finite gap enclosure")
    return gap


# ---------------------------------------------------------------------------
# The exact q partial: dual-Arb propagation through the shared formula.

class _DBall:
    """First-order dual number over Arb balls: (value, derivative) with the
    ring laws +, -, * (no division is ever needed: the transcription forms
    no quotient on the claim path -- _GapOps.quotient returns None)."""

    __slots__ = ("value", "derivative")

    def __init__(self, value: Any, derivative: Any = None):
        self.value = arb(value)
        self.derivative = arb(0) if derivative is None else arb(derivative)

    @staticmethod
    def _coerce(other: Any) -> "_DBall":
        return other if isinstance(other, _DBall) else _DBall(arb(other))

    def __add__(self, other: Any) -> "_DBall":
        other = self._coerce(other)
        return _DBall(self.value + other.value,
                      self.derivative + other.derivative)

    __radd__ = __add__

    def __sub__(self, other: Any) -> "_DBall":
        other = self._coerce(other)
        return _DBall(self.value - other.value,
                      self.derivative - other.derivative)

    def __rsub__(self, other: Any) -> "_DBall":
        return self._coerce(other) - self

    def __mul__(self, other: Any) -> "_DBall":
        other = self._coerce(other)
        return _DBall(self.value * other.value,
                      self.derivative * other.value
                      + self.value * other.derivative)

    __rmul__ = __mul__

    def frozen(self) -> "_DBall":
        """The dual with its derivative discarded (value preserved): the
        mutation primitive for dropping a derivative contribution."""
        return _DBall(self.value, arb(0))


def _hp_ball(x: arb) -> arb:
    """Enclosure of h'(y) = log((1-y)/y) on a ball whose closure lies
    strictly inside (0, 1): h' is antimonotone, so the two corner
    evaluations cover the exact range; one union certifies the ball."""
    if not (x.lower() > 0 and x.upper() < 1):
        raise ValueError("h' enclosure needs its argument ball strictly "
                         "inside (0, 1)")
    lo_corner = (arb(1) - x.upper()).log() - x.upper().log()
    hi_corner = (arb(1) - x.lower()).log() - x.lower().log()
    return lo_corner.union(hi_corner)


def _entropy_dual(argument: Any) -> Any:
    """Chain-rule entropy over a dual (or plain ball at derivative zero)."""
    if isinstance(argument, _DBall):
        return _DBall(h_arb(argument.value),
                      _hp_ball(argument.value) * argument.derivative)
    return h_arb(arb(argument))


def q_partial_arb(values9: Sequence[Any], beta: Any,
                  cross_contribution: bool = True) -> arb:
    """CLAIM (Q2).  The EXACT partial derivative d gap / d q of the raw
    gap, enclosed on Arb balls, by dual-number propagation through the
    repository's single shared transcription liu9_objective._formula (the
    liu9_boundary_layer autodiff pattern over Arb balls).  NO finite
    difference is taken anywhere on this path.

    values9 = (a1, a2, q, b0, b2, b4, b1, b3, b5); masses (positions 0, 1,
    with a3 = 1 - a1 - a2) exact simplex rationals or point balls, supports
    (3..8) and q (2) Arb balls, beta the certified ball.  Every
    intermediate step is inclusion-monotone interval arithmetic, so the
    returned ball encloses the point value of transcribed-formula's exact
    partial at EVERY point of the input boxes.

    Domain discipline: support balls strictly inside (0, 1) (the chain-rule
    h' enclosure); the q ball may touch 0 or 1 (one-sided endpoint boxes:
    q never enters an entropy argument, only the weight factors).

    cross_contribution=False runs the MUTANT evaluation used by mutation
    (a): the cross-component derivative contribution (the derivative of the
    weight product w_i * w_j for between-component EHXY pairs, i in the P0
    block and j in the P1 block and symmetrically) is dropped -- the
    product is evaluated with frozen (derivative-free) weights.  The value
    chain is unchanged; only the derivative loses those pieces, so the
    result is NOT the transcription's q partial and the finite-difference
    PROBE must refute agreement.  Implemented by the explicit block walker
    in _gap_q_walker, which is itself cross-validated against this shared
    path when cross_contribution=True (containment, byte-stable)."""
    if cross_contribution:
        duals = tuple(_DBall(value if isinstance(value, arb)
                             else point_ball(value),
                             arb(1) if position == 2 else arb(0))
                      for position, value in enumerate(values9))
        beta_dual = _DBall(beta if isinstance(beta, arb)
                           else point_ball(beta))
        terms = _formula(duals, beta_dual,
                         _GapOps(_entropy_dual, _DBall(arb(1))))
        gap = terms.numerator - terms.ehx
        out = gap.derivative
        if not out.is_finite():
            raise ValueError("non-finite q partial enclosure")
        return out
    return _gap_q_walker(values9, beta, cross_terms=False).derivative


def _gap_q_walker(values9: Sequence[Any], beta: Any,
                  cross_terms: bool = True) -> _DBall:
    """Explicit block-structured dual evaluation of the transcribed raw
    gap, mirroring liu9_objective._formula's summation semantics exactly
    (same 36-term EHXY sum, 9-term per-component PI sums, 6-term EHX sum,
    numerator, no quotient).  cross_terms=False freezes the between-
    component EHXY weight products (both the 01 and the symmetric 10
    block) at their value balls -- the mutation path ONLY; the claim path
    is q_partial_arb through the shared formula."""
    if len(values9) != 9:
        raise ValueError("expected nine variables")
    wraps = tuple(_DBall(value if isinstance(value, arb)
                         else point_ball(value),
                         arb(1) if position == 2 else arb(0))
                  for position, value in enumerate(values9))
    a1, a2, q, b0, b2, b4, b1, b3, b5 = wraps
    one = _DBall(arb(1))
    a3 = one - a1 - a2
    masses = (a1, a2, a3)
    p0, p1 = (b0, b2, b4), (b1, b3, b5)
    qbar = one - q
    support = p0 + p1
    weights = tuple(component * mass
                    for component in (qbar, q) for mass in masses)
    beta_dual = _DBall(beta if isinstance(beta, arb) else point_ball(beta))
    zero = _DBall(arb(0))

    def prod3(weight_i: _DBall, weight_j: _DBall) -> _DBall:
        cross = (weight_i is weights[0] or weight_i is weights[1]
                 or weight_i is weights[2]) and \
                (weight_j is weights[3] or weight_j is weights[4]
                 or weight_j is weights[5])
        if cross_terms or not cross:
            return weight_i * weight_j
        return _DBall(weight_i.value * weight_j.value)  # FROZEN (mutant)

    ehxy = zero
    for i in range(6):
        for j in range(6):
            term = prod3(weights[i], weights[j])
            term = term * _entropy_dual(support[i] * support[j])
            ehxy = ehxy + term

    def component_pi(points: Sequence[_DBall]) -> _DBall:
        out = zero
        for i in range(3):
            for j in range(3):
                x, y = points[i], points[j]
                arg = x * y + x * (one - x) * y * (one - y)
                out = out + masses[i] * masses[j] * _entropy_dual(arg)
        return out

    ehpi = (qbar * component_pi(p0)) + (q * component_pi(p1))
    ehx = zero
    for i in range(6):
        ehx = ehx + weights[i] * _entropy_dual(support[i])
    numerator = ((one - beta_dual) * ehxy) + (beta_dual * ehpi)
    gap = numerator - ehx
    if not gap.value.is_finite() or not gap.derivative.is_finite():
        raise ValueError("non-finite walker enclosure")
    return gap


# ---------------------------------------------------------------------------
# (Q1) the exact v-quadratic at a frozen fiber.

def v_fiber_coefficients(a1: Fraction, a2: Fraction,
                         supports: tuple[Fraction, ...], beta: Any
                         ) -> dict[str, arb]:
    """CLAIM (Q1) machinery.  The three certified coefficient balls of the
    exact q-quadratic at the frozen fiber (supports, (a1, a2)): the exact
    rational q fibers 0, 1/2, 1 feed the shared-formula value enclosure
    gap_arb, and the exact inverse of the 3-point Vandermonde
        V = [[1, 0, 0], [1, 1/2, 1/4], [1, 1, 1]],
        V^-1 = [[1, 0, 0], [-3, 4, -1], [2, -4, 2]]
    is applied inclusion-monotonically to the three certified value balls:
        A = v(0),  C = 2 v(0) - 4 v(1/2) + 2 v(1),  B = v(1) - v(0) - C.
    NO finite difference anywhere.  Monomial order (A, B, C) w.r.t.
    (1, q, q^2):  v(q) = A + B q + C q^2."""
    if len(supports) != 6:
        raise ValueError("expected six supports (b0 b2 b4 b1 b3 b5)")
    masses = (point_ball(a1), point_ball(a2), point_ball(ONE - a1 - a2))
    sballs = tuple(point_ball(s) for s in supports)
    betab = beta if isinstance(beta, arb) else point_ball(beta)
    values = {}
    for q in (ZERO, Fraction(1, 2), ONE):
        values[q] = gap_arb((masses[0], masses[1], point_ball(q), *sballs),
                            betab)
    v0, vh, v1 = values[ZERO], values[Fraction(1, 2)], values[ONE]
    two = point_ball(2)
    three = point_ball(3)
    four = point_ball(4)
    coef_a = v0
    # C = row 3 of V^-1: 2 v(0) - 4 v(1/2) + 2 v(1), i.e. the second
    # difference over the 1/2 step divided by (1/2)^2:
    #     v(1) - 2 v(1/2) + v(0) = C * (1/2)^2.
    coef_c = (two * v1 - four * vh) + two * v0
    # B = row 2 of V^-1: -3 v(0) + 4 v(1/2) - v(1)  (= v(1) - v(0) - C).
    coef_b = (four * vh - three * v0) - v1
    return {"A": coef_a, "B": coef_b, "C": coef_c,
            "_v0": v0, "_vh": vh, "_v1": v1}


def v_quadratic_at(coefficients: dict[str, Any], q: Any) -> arb:
    """Certified v(q) at a q point/ball: inclusion-monotone evaluation of
    the certified coefficient balls (one ball per input ball)."""
    coef_a, coef_b, coef_c = (coefficients["A"], coefficients["B"],
                              coefficients["C"])
    qb = q if isinstance(q, arb) else point_ball(q)
    return coef_a + coef_b * qb + coef_c * (qb * qb)


def q_point_minimization(coefficients: dict[str, Any],
                         q_box: tuple[Fraction, Fraction]) -> dict[str, Any]:
    """CLAIM (Q5).  Exact q-quadratic point minimization with endpoints
    plus interior critical point, over the q box (for later
    branch-and-bound).  The critical point -B/(2C) is a CANDIDATE only
    when the certified coefficient brackets decide it: C's ball strictly
    positive (open upward) and B's ball strictly negative, and the critical
    ball strictly inside the box.  No interval-argmin heuristic is used:
    every candidate is an exact rational or a certified ball argument, and
    the deterministic tie-break is (value lower endpoint, branch name).

    Soundness of the returned lower bound: for the TRUE coefficients
    (inside the certified balls) the box minimum is attained at one of the
    candidate arguments (endpoint or certified-interior critical); each
    candidate's ball covers v(candidate argument) for ALL coefficient
    combos including the true one; hence
        min over box of v >= min over candidates of (value ball lower).
    Non-winning candidates only weaken the bound, never invalidate it."""
    coef_a, coef_b, coef_c = (coefficients["A"], coefficients["B"],
                              coefficients["C"])
    qlo, qhi = q_box
    if qlo > qhi or qlo < 0 or qhi > 1:
        raise ValueError("q box must be a nonempty subinterval of [0, 1]")
    candidates = [
        {"branch": "q_left", "arg": point_ball(qlo),
         "value": v_quadratic_at(coefficients, point_ball(qlo))},
        {"branch": "q_right", "arg": point_ball(qhi),
         "value": v_quadratic_at(coefficients, point_ball(qhi))},
    ]
    critical: Optional[arb] = None
    if bool(coef_c.lower() > 0) and bool(coef_b.lower() < 0):
        neg_half_b = -(coef_b / point_ball(2))
        critical = neg_half_b / coef_c
        if (bool(critical.lower() > point_ball(qlo))
                and bool(critical.upper() < point_ball(qhi))):
            candidates.append({"branch": "q_interior", "arg": critical,
                               "value": v_quadratic_at(coefficients,
                                                       critical)})
    candidates.sort(key=lambda row: (row["value"].lower(), row["branch"]))
    return {"min": candidates[0]["value"], "branch": candidates[0]["branch"],
            "critical": critical, "candidates": candidates,
            "q_box": [str(qlo), str(qhi)]}


def v_residual_check(mass_fibers: Sequence[tuple[Fraction, Fraction]],
        supports: tuple[Fraction, ...], beta: Any,
        probes: Sequence[Fraction] = Q_PROBES) -> dict[str, Any]:
    """CLAIM (Q1) guard.  Deterministic q-degree/residual validation: at
    each exact rational q probe and each exact rational mass fiber, the
    per-fiber extracted v-quadratic and the direct shared-formula gap
    enclosure must intersect (the residual ball contains zero).  This
    certifies exact degree <= 2 in q at every probed fiber; the residual's
    WIDTH is reported as diagnostics only, never as a claim."""
    sballs = tuple(point_ball(s) for s in supports)
    betab = beta if isinstance(beta, arb) else point_ball(beta)
    rows = []
    ok = True
    for a1, a2 in mass_fibers:
        # Each mass fiber has its OWN exact v-quadratic: extract it
        # per fiber (the coefficients at a different fiber say nothing
        # about this one).
        coefficients = v_fiber_coefficients(a1, a2, supports, betab)
        masses = (point_ball(a1), point_ball(a2), point_ball(ONE - a1 - a2))
        for q in probes:
            direct = gap_arb((masses[0], masses[1], point_ball(q), *sballs),
                             betab)
            fit = v_quadratic_at(coefficients, point_ball(q))
            residual = fit - direct
            contains = bool(residual.contains(0))
            ok = ok and contains
            rows.append({
                "mass": [str(a1), str(a2)], "q": str(q),
                "direct": direct.str(26), "quadratic": fit.str(26),
                "residual": residual.str(26),
                "residual_width_float": float(residual.upper()
                                              - residual.lower()),
                "contains_zero": contains,
            })
    return {"probes": rows, "ok": bool(ok),
            "q_degree": "exact degree <= 2 in q at each frozen fiber "
                        "(Q1); residual containment certified at the "
                        "probes, widths are diagnostics only"}


# ---------------------------------------------------------------------------
# (Q3)/(Q4) uniform mass-simplex machinery.

def uniform_q_partial_balls(
        support_boxes: Sequence[tuple[Fraction, Fraction]],
        q_box: tuple[Fraction, Fraction], beta: Any) -> dict[str, Any]:
    """The six certified node balls of d gap / d q at the common
    MASS_NODES (documented A, B, C, AB, AC, BC order): at each node the
    masses are EXACT simplex rationals (a1, a2, a3 = 1 - a1 - a2), while
    the supports and q stay BALLS (the caller's box).  Every returned ball
    encloses the exact partial at every point of (support boxes) x (q box)
    at that node's exact masses.  Boundary mass nodes are admissible
    (masses never enter an entropy argument); support boxes touching 0 or
    1 are REJECTED (the h' enclosure needs strict interior), the q box may
    touch 0 or 1 (one-sided endpoints: q never enters an entropy
    argument)."""
    if len(support_boxes) != 6:
        raise ValueError("expected six support boxes (b0 b2 b4 b1 b3 b5)")
    for lo, hi in support_boxes:
        if not (lo > 0 and hi < 1):
            raise ValueError(
                "support boxes touching 0 or 1 are rejected by this "
                "module's derivative enclosures (strict interior required)")
    if q_box[0] < 0 or q_box[1] > 1 or q_box[0] > q_box[1]:
        raise ValueError("q box must be a nonempty subinterval of [0, 1]")
    sballs = tuple(point_ball(lo).union(point_ball(hi))
                   for lo, hi in support_boxes)
    qball = point_ball(q_box[0]).union(point_ball(q_box[1]))
    betab = beta if isinstance(beta, arb) else point_ball(beta)
    nodes = []
    for a1, a2 in MASS_NODES:
        masses = (point_ball(a1), point_ball(a2), point_ball(ONE - a1 - a2))
        nodes.append(q_partial_arb((masses[0], masses[1], qball, *sballs),
                                   betab))
    return {"nodes": tuple(nodes), "q_ball": qball, "support_balls": sballs,
            "beta": betab}


def q_lipschitz(beta: Any,
                support_boxes: Sequence[tuple[Fraction, Fraction]],
                q_box: tuple[Fraction, Fraction]) -> arb:
    """CLAIM (Q3).  A certified uniform bound L for
        max { |d gap/d q(s, q, a)| : s in the support boxes, q in the q
              box, a in the full mass simplex T },
    computed by routing the six node balls of
    uniform_q_partial_balls through
    liu9_h2_envelope_common.simplex_abs_upper: every degree-<=2 mass
    polynomial in the family { d gap/d q(., ., a) : a in T } (degree 2 in
    the masses: the weight products a_j * w_k and the linear a_j in the
    x-part) has its six node values inside the six node balls, the
    node-to-control map is inclusion-monotone, and the convex-hull range
    property on T bounds the polynomial by the control range; taking the
    absolute-value upper endpoint gives the uniform bound."""
    node_balls = uniform_q_partial_balls(support_boxes, q_box, beta)["nodes"]
    bound = simplex_abs_upper(node_balls)
    if not bound.is_finite():
        raise ValueError("non-finite uniform Lipschitz bound")
    return bound.upper()


def q_drift_bound(beta: Any,
                  support_boxes: Sequence[tuple[Fraction, Fraction]],
                  q_box: tuple[Fraction, Fraction]) -> dict[str, Any]:
    """CLAIM (Q4).  The uniform drift: for EVERY q, q' in the q box, EVERY
    a in the full mass simplex T and EVERY s in the support boxes,
        |gap(q', a, s) - gap(q, a, s)| <= (q' - q) * L,
    with L = q_lipschitz(...) (Q3), by the fundamental theorem: the raw gap
    is a polynomial (hence C^1) in q, its q partial is enclosed by L
    uniformly on (support boxes) x (q box) x T, and the segment between q
    and q' lies inside the q box.  The end q = 0 or 1 case is the one-sided
    endpoint box case (allowed here); gap is continuous there and the
    cited universal endpoint theorem (see module docstring) supplies the
    nonnegativity boundary data, NOT this drift."""
    span = q_box[1] - q_box[0]
    lipschitz = q_lipschitz(beta, support_boxes, q_box)
    drift = point_ball(span) * lipschitz
    return {
        "lipschitz": lipschitz,
        "lipschitz_str": lipschitz.str(20),
        "span": span,
        "drift": drift,
        "drift_str": drift.str(20),
        "statement": (
            "CLAIM (Q4): for every q, q' in the q box, every a in the "
            "full mass simplex T and every s in the support boxes, "
            "|gap(q', a, s) - gap(q, a, s)| <= span * L with L the "
            "certified uniform q-Lipschitz constant (Q3); fundamental "
            "theorem on the polynomial-in-q gap."),
    }


# ---------------------------------------------------------------------------
# Deterministic mutation hooks (the claim runner asserts observed FAIL).

def run_mutations(beta: arb, supports: tuple[Fraction, ...],
                  support_boxes: Sequence[tuple[Fraction, Fraction]],
                  q_box: tuple[Fraction, Fraction]) -> list[dict[str, Any]]:
    """Two mutations, both of which MUST FAIL (their passing would mean the
    claim machinery is insensitive to its own inputs, mirroring the sweep's
    mutation semantics):

    (a) drop_cross_component_q_derivative: the finite-difference PROBE is
        run against the MUTANT q partial (cross_contribution=False: the
        between-component EHXY weight products enter the chain with frozen
        derivatives).  The mutant must NOT agree with the certified FD
        quotient -- i.e. the FD agreement probe genuinely discriminates a
        dropped cross-component contribution.
    (b) coef_perturbed_1e-6: perturb the v-quadratic's C coefficient by
        +1e-6; the value at the same exact rational q point must MOVE and
        the certified q-box minimum (Q5 re-solve) must DIFFER from the
        reference."""
    out: list[dict[str, Any]] = []
    q_mid = (q_box[0] + q_box[1]) / 2
    masses = (point_ball(DEMO_MASSES[0]), point_ball(DEMO_MASSES[1]),
              point_ball(ONE - DEMO_MASSES[0] - DEMO_MASSES[1]))
    mid_balls = tuple(point_ball((lo + hi) / 2) for lo, hi in support_boxes)
    point_inputs = (masses[0], masses[1], point_ball(q_mid), *mid_balls)

    # --- mutation (a): the FD probe must refute the cross-dropped mutant.
    h = FD_STEP
    direct_hi = gap_arb((masses[0], masses[1], point_ball(q_mid + h),
                         *mid_balls), beta)
    direct_lo = gap_arb((masses[0], masses[1], point_ball(q_mid - h),
                         *mid_balls), beta)
    fd = (direct_hi - direct_lo) / point_ball(2 * h)
    true_partial = q_partial_arb(point_inputs, beta)
    mutant = q_partial_arb(point_inputs, beta, cross_contribution=False)
    separated = not bool((fd - mutant).contains(0))
    agreement_true = bool((fd - true_partial).contains(0))
    out.append({
        "mutation": "drop_cross_component_q_derivative",
        "mechanism": (
            "the q partial's cross-component EHXY weight-product "
            "derivatives (01 block and symmetric 10 block) are dropped "
            "(frozen weights); the certified finite-difference probe "
            "must REFUTE agreement -- without the cross pieces the "
            "mutant is not the transcription's q partial"),
        "expected": "FAIL",
        "observed": "FAIL" if separated else "PASSED",
        "fd_ball": fd.str(20), "true_partial": true_partial.str(20),
        "mutant_partial": mutant.str(20),
        "fd_true_agreement_probe": bool(agreement_true),
        "fd_mutant_difference": (fd - mutant).str(20),
        "ok": bool(separated),
    })
    if not agreement_true:
        raise AssertionError("FD probe disagrees with the true q partial; "
                             "the dual engine or the probe is broken")

    # --- mutation (b): perturbed q-quadratic coefficient.
    coefficients = v_fiber_coefficients(DEMO_MASSES[0], DEMO_MASSES[1],
                                        supports, beta)
    perturbed = dict(coefficients)
    perturbed["C"] = coefficients["C"] + point_ball(MUT_COEF_DELTA)
    ref_min = q_point_minimization(coefficients, q_box)
    mut_min = q_point_minimization(perturbed, q_box)
    moved_value = not bool((v_quadratic_at(perturbed, point_ball(q_mid))
                            - v_quadratic_at(coefficients,
                                             point_ball(q_mid))).contains(0))
    moved_min = not bool((mut_min["min"] - ref_min["min"]).contains(0))
    out.append({
        "mutation": "coef_perturbed_1e-6",
        "mechanism": (
            "C -> C + 1e-6 on the certified v-quadratic: the value at "
            "the same exact rational q point must move and the "
            "certified q-box minimum must differ from the reference"),
        "expected": "FAIL",
        "observed": "FAIL" if (moved_value and moved_min) else "PASSED",
        "q_reference": str(q_mid),
        "value_reference": v_quadratic_at(coefficients,
                                          point_ball(q_mid)).str(20),
        "value_mutated": v_quadratic_at(perturbed,
                                        point_ball(q_mid)).str(20),
        "min_reference": ref_min["min"].str(20),
        "min_mutated": mut_min["min"].str(20),
        "branch_reference": ref_min["branch"],
        "branch_mutated": mut_min["branch"],
        "ok": bool(moved_value and moved_min),
    })
    return out


# ---------------------------------------------------------------------------
# The deterministic runner / report (byte-stable; no timestamps).

def build_report() -> dict[str, Any]:
    beta = certified_beta_ball()
    beta_enclosure = exact_enclosure(beta)
    supports = DEMO_XTRIPLE + DEMO_YTRIPLE
    support_boxes = tuple((s - DEMO_WH, s + DEMO_WH) for s in supports)
    q_box = DEMO_Q_BOX

    # ---- CLAIM (Q1): exact v-quadratic at the demo fiber + residual guard.
    coefficients = v_fiber_coefficients(DEMO_MASSES[0], DEMO_MASSES[1],
                                        supports, beta)
    coefficient_record = {key: coefficients[key].str(26)
                          for key in ("A", "B", "C")}
    q1 = v_residual_check([DEMO_MASSES, SECOND_FIBER], supports, beta)

    # ---- PROBES anchoring the dual engine (never claims).
    q_mid = (q_box[0] + q_box[1]) / 2
    masses = (point_ball(DEMO_MASSES[0]), point_ball(DEMO_MASSES[1]),
              point_ball(ONE - DEMO_MASSES[0] - DEMO_MASSES[1]))
    mid_balls = tuple(point_ball((lo + hi) / 2) for lo, hi in support_boxes)
    point_inputs = (masses[0], masses[1], point_ball(q_mid), *mid_balls)
    direct_hi = gap_arb((masses[0], masses[1], point_ball(q_mid + FD_STEP),
                         *mid_balls), beta)
    direct_lo = gap_arb((masses[0], masses[1], point_ball(q_mid - FD_STEP),
                         *mid_balls), beta)
    fd = (direct_hi - direct_lo) / point_ball(2 * FD_STEP)
    dual_partial = q_partial_arb(point_inputs, beta)
    fd_agree = bool((fd - dual_partial).contains(0))
    walker = _gap_q_walker(point_inputs, beta, cross_terms=True)
    walker_agree = bool((walker.derivative - dual_partial).contains(0))
    walker_ok = bool((walker.value - gap_arb(point_inputs, beta)).contains(0))
    probes = {
        "fd_symmetric_quotient": {
            "note": "PROBE: for a degree-<=2 polynomial the symmetric "
                    "difference quotient is EXACT (no truncation term), "
                    "so this is a sharp consistency check of the dual "
                    "engine; it is validation, never a claim",
            "step": str(FD_STEP), "q_mid": str(q_mid),
            "fd_ball": fd.str(26), "dual_partial": dual_partial.str(26),
            "difference": (fd - dual_partial).str(26),
            "agree": fd_agree,
        },
        "walker_containment": {
            "note": "PROBE: the explicit block walker must agree with "
                    "the shared-formula dual path on both value and "
                    "derivative at the same exact point inputs",
            "walker_value": walker.value.str(26),
            "shared_value": gap_arb(point_inputs, beta).str(26),
            "walker_partial": walker.derivative.str(26),
            "shared_partial": dual_partial.str(26),
            "value_agree": walker_ok, "derivative_agree": walker_agree,
        },
    }

    # ---- CLAIM (Q3)/(Q4): uniform Lipschitz and drift over the full T.
    node_balls = uniform_q_partial_balls(support_boxes, q_box, beta)
    node_record = [ball.str(26) for ball in node_balls["nodes"]]
    control_balls = simplex_bernstein_coefficients(node_balls["nodes"])
    range_balls = simplex_range(node_balls["nodes"])
    lipschitz = q_lipschitz(beta, support_boxes, q_box)
    q4 = q_drift_bound(beta, support_boxes, q_box)

    # ---- CLAIM (Q5): exact q-quadratic point minimization.
    q5 = q_point_minimization(coefficients, q_box)
    q5_record = {
        "min": q5["min"].str(26), "branch": q5["branch"],
        "critical": None if q5["critical"] is None
        else q5["critical"].str(20),
        "candidates": [{"branch": row["branch"],
                        "arg": row["arg"].str(20),
                        "value": row["value"].str(26)}
                       for row in q5["candidates"]],
        "q_box": [str(q_box[0]), str(q_box[1])],
    }

    mutations = run_mutations(beta, supports, support_boxes, q_box)
    claims_ok = bool(q1["ok"] and fd_agree and walker_agree and walker_ok)
    mutations_ok = all(row["ok"] for row in mutations)
    claim_status = ("Q-ENVELOPE-CERTIFIED" if claims_ok and mutations_ok
                    else "FAILED")

    report = {
        "tool": "uc/liu9_h2_q_envelope.py",
        "claim": (
            "Q-ENVELOPE for the paired r=3 raw gap: at the frozen exact "
            "rational support boxes and q box, (Q1) the raw gap is exactly "
            "quadratic in q at every frozen mass fiber (residual "
            "containment certified at exact rational probes and a second "
            "mass fiber); (Q2) the q partial is the exact arithmetic "
            "expression through liu9_objective._formula (dual-Arb, no "
            "finite difference in any claim path); (Q3) the uniform bound "
            "max |d gap/d q| over (support boxes) x (q box) x FULL mass "
            "simplex via the common degree-2 Bernstein map; (Q4) the "
            "uniform drift |gap(q')-gap(q)| <= span*L over the same "
            "product; (Q5) the exact q-quadratic point minimization with "
            "endpoints plus interior critical point for "
            "branch-and-bound.  No fixed-argmin radius shortcut is used "
            "anywhere."),
        "claim_status": claim_status,
        "claims": {
            "Q1_exact_v_quadratic": {
                "coefficients": coefficient_record,
                "residual_check": q1,
                "ok": bool(q1["ok"]),
            },
            "Q2_exact_q_partial": {
                "engine": "dual-number propagation through the shared "
                          "transcription liu9_objective._formula over Arb "
                          "balls (liu9_boundary_layer autodiff pattern); "
                          "no finite difference in any claim path",
                "probe_reference": "probes.fd_symmetric_quotient (PROBE)",
                "ok": bool(fd_agree and walker_agree),
            },
            "Q3_uniform_lipschitz": {
                "lipschitz": lipschitz.str(26),
                "lipschitz_float": float(lipschitz.upper()),
                "node_balls": node_record,
                "control_balls": [ball.str(26) for ball in control_balls],
                "control_range": [range_balls[0].str(20),
                                  range_balls[1].str(20)],
                "ok": bool(lipschitz.is_finite()),
            },
            "Q4_uniform_drift": {
                "lipschitz": q4["lipschitz_str"],
                "span": str(q4["span"]),
                "drift": q4["drift_str"],
                "drift_float": float(q4["drift"].upper()),
                "statement": q4["statement"],
                "ok": bool(q4["drift"].is_finite()),
            },
            "Q5_point_minimization": q5_record,
        },
        "probes": probes,
        "mutations": mutations,
        "boundary_theorem": {
            "citation": (
                "UNIVERSAL ENDPOINT THEOREM (cited, not reproved): "
                "uc/liu9_size_biased.py certify_qone_kernel "
                "(claim_status PROVED) proves F(P) >= 0 for EVERY law P "
                "of mean >= m via F(P) = M E_{Q x Q} C_M with "
                "Q(dx) = x P(dx)/M, C_m >= 0 pointwise and dC_M/dM >= 0.  "
                "The exact split gap = F(P_mix) + q(1-q) c_ch (proved, "
                "r <= 3); the channel term vanishes identically at "
                "q in {0, 1}, so gap(q=1) = F(P1) and gap(q=0) = F(P0), "
                "and the endpoint mixture IS that component law with the "
                "mixture mean as its mean: the theorem discharges BOTH "
                "endpoint q faces as cited boundary data."),
            "channel_sign_note": (
                "No sign is claimed for c_ch: asserting c_ch >= 0 would "
                "be FALSE.  The certified facts are the domination ratio "
                "q(1-q)|c|/margin <= 0.838883838221919 on the 9,747 "
                "solved sweep cells (zero danger cells) and the "
                "certified most-negative channel scalar c = -0.02755 "
                "(at a mass vertex, all 1,225 geometries)."),
            "role": "boundary data for the one-sided endpoint q boxes "
                    "(q=0 and q=1 faces); this module's own claims are "
                    "the Q1-Q5 statements above",
        },
        "beta_certified_ball": {
            "ball": beta.str(30),
            "exact_enclosure": [str(beta_enclosure[0]),
                                str(beta_enclosure[1])],
            "note": "re-derived per run from liu9_binding's equation-defined "
                    "chain at 320 bits (deterministic); the exact rational "
                    "enclosure's containment of the certified ball is "
                    "asserted on every run",
        },
        "demo_fiber": {
            "x_triple": [str(v) for v in DEMO_XTRIPLE],
            "y_triple": [str(v) for v in DEMO_YTRIPLE],
            "masses": [str(v) for v in DEMO_MASSES],
            "q_box": [str(v) for v in DEMO_Q_BOX],
            "support_half_width": str(DEMO_WH),
            "support_boxes": [[str(lo), str(hi)]
                              for lo, hi in support_boxes],
            "interior_discipline": (
                "all support boxes strictly inside (0, 1); the q box is "
                "interior in the demo and may touch {0, 1} one-sidedly "
                "in general (q never enters an entropy argument)"),
        },
        "limitations": [
            "BETWEEN-BOX GAPS: the uniform claims are stated ON the "
            "support boxes x q box x mass simplex product; the ambient "
            "parameter space between such boxes is NOT covered here.",
            "The uniform Lipschitz constant (Q3) is a Bernstein "
            "convex-hull over-estimate of the true sup (sound one-sided "
            "slack); the drift (Q4) inherits that slack through the "
            "fundamental theorem, so the drift bound is rigorous but "
            "not sharp.",
            "The interior critical candidate in (Q5) participates only "
            "when the certified coefficient brackets decide BOTH signs "
            "(C strictly positive, B strictly negative) and the critical "
            "ball lies strictly inside the box; otherwise only the "
            "endpoints are candidates -- a documented weakening, never a "
            "heuristic.",
            "All claims are 320-bit Arb ball enclosures of exact "
            "quantities; the printed decimal strings are roundings of "
            "certified balls, and floats appear only as display "
            "diagnostics (float-suffixed keys), never as claim data.",
            "Endpoint q faces (q=0, q=1) are discharged by the CITED "
            "universal size-biased theorem as boundary data; this module "
            "neither reproves nor re-encodes it, and its own drift/"
            "minimization claims make no sign statement about c_ch.",
            "The FD and walker validations are PROBES: they gate the "
            "claim_status deterministically but the claims' soundness "
            "rests only on the exact-arithmetic constructions (Q1-Q5).",
        ],
        "evaluator": (
            "liu9_objective._formula via liu9_boundary_layer._GapOps "
            "(gap = numerator - ehx, no quotient); dual q partials per "
            "liu9_boundary_layer's autodiff pattern; common helper "
            "uc/liu9_h2_envelope_common for the nodal degree-2 Bernstein "
            "map; ctx.prec = 320; beta from liu9_binding's certified "
            "equation-defined chain"),
    }
    body = json.dumps(report, sort_keys=True, indent=1)
    report["report_sha256"] = hashlib.sha256(body.encode("utf-8")).hexdigest()
    return report


def main() -> int:
    report = build_report()
    OUTPUT_DEFAULT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_DEFAULT.write_text(
        json.dumps(report, sort_keys=True, indent=1) + "\n")
    print("claim_status:", report["claim_status"])
    q1 = report["claims"]["Q1_exact_v_quadratic"]
    print("Q1 v-quadratic ok:", q1["ok"],
          "coefficients A, B, C:",
          q1["coefficients"]["A"], q1["coefficients"]["B"],
          q1["coefficients"]["C"])
    q3 = report["claims"]["Q3_uniform_lipschitz"]
    print("Q3 uniform Lipschitz:", q3["lipschitz"])
    q4 = report["claims"]["Q4_uniform_drift"]
    print("Q4 uniform drift over span", q4["span"], ":", q4["drift"])
    q5 = report["claims"]["Q5_point_minimization"]
    print("Q5 box min:", q5["min"], "branch:", q5["branch"],
          "critical:", q5["critical"])
    for row in report["mutations"]:
        print("mutation %-34s %s" % (row["mutation"], row["observed"]))
    print("report_sha256", report["report_sha256"])
    return 0 if report["claim_status"] == "Q-ENVELOPE-CERTIFIED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
