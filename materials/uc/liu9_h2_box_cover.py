#!/usr/bin/env python3
"""RETRACTED 2026-09-01 - THIS MODULE CERTIFIES NOTHING.  DO NOT CITE IT.

Defect (fatal, inequality direction).  ``radius_certificate`` and
``SOUNDNESS_ARGUMENT`` step (vi) assert, per box fiber (s, q),

    min_P quad_{A(s,q)} >= quad_{A(s,q)}(a1*, a2*) - R

and justify it in-line with "the fiber quadratic's min over P is <= its
value at the P-point (a1*, a2*)".  That justification establishes the
OPPOSITE inequality; ``min <= value at a point`` gives no lower bound.
R is built from a single ``evaluate_arb`` call at the FIXED exact-rational
masses (a1*, a2*) - the center fiber's argmin - so it bounds the gap's
oscillation at one mass point only, never the movement of the fiber
quadratic's minimizer in the mass directions.  Consequently BOXBOUND is
not a lower bound and every clearing count, worst bound, and covered
volume this module ever reported is withdrawn.

The sound alternative stated in the original docstring step (6),
BOXBOUND := c0lo - sum_j width(band_j), is genuinely sound but VACUOUS at
every affordable width: measured at cell g1117|q2, sum_j width(band_j) is
21651.78 at h = 2^-10, 1353.19 at 2^-14, and still 1.3215 at 2^-24,
against a margin of 1.6e-3.  The 6-point design at spacing
H_STEP = 1/400 makes the Cramer inverse amplify every interval width by
~400^2 = 1.6e5, so the sound form would need h ~ 2^-34 (box volume
~1e-70).  The method cannot be repaired by shrinking boxes.

``main`` therefore refuses to run without --retracted-replay, which
relabels the report and never emits a certified status.  The point sweep
this module was built on (uc/liu9_h2_qp_sweep.py ->
liu9-h2-qp-sweep.json, worst certified min +1.6305e-3 at exact rational
supports and q) is UNAFFECTED.  See PROGRESS.md RETRACTION-H2-BOX-COVER
(2026-09-01).

Everything below this banner is the original, retracted text.

H2 BOX COVER over support/q neighborhoods of the certified H2 mass-QP
sweep grid - the mission box-step deliverable for the paired r = 3 class.

WHAT IS ALREADY CERTIFIED (cited, not re-derived: Main's SSOT).  For the
paired r = 3 class the raw gap splits exactly as

    gap = F(P_mix) + q(1-q) c_ch

and at frozen rational supports gap is exactly bidegree-(2,2) in
((a1, a2), q) over the formal H-atom ring, so at every frozen fiber
(geometry, q) the gap is an EXACT quadratic in the two free masses.
uc/liu9_h2_qp_sweep.py swept 1225 geometry pairs x 9 q columns = 11,025
cells with closed-form KKT argmin over the mass polygon
{a1>=0, a2>=0, a1+a2<=1, mean >= M_FLOOR}, a superset of the TRUE
feasible set {mean >= m} because M_FLOOR < m (Arb-verified).  9,747
cells certified min >= +1.6305e-3 (worst cell g9|q8); 1,278 are
infeasible-geometry: the mixture-mean line never reaches M_FLOOR over
the mass triangle, the polygon there is EMPTY and nothing is certified
there (excluded from every count, as in the sweep).

THIS MODULE.  A certified box cover of the neighborhoods of those
cells: for each nonvacuous sweep cell (support triples + q0) and each
half-width h = 2^-w, w in {10, 12, 14} (the ladder of
liu9_scalar_margin_r3.py), certify

    inf{ gap(supports in [g_i - h, g_i + h],
             q in [q0-h, q0+h] clipped to [0, 1],
             (a1, a2) feasible for the TRUE mean at that box point) }
        >=  0.

THE RELAXATION (coefficient-band + drifted-KKT bound; the certified
claim).  All statements are about exact reals; interval arithmetic
(Arb balls at 320-bit) enters ONLY as a rigorous two-sided enclosure
of an exact quantity:

  (1) SUPPORT-BOX ORDERING.  h <= 2^-10 < 1/16 and distinct support
      grid values differ by >= 1/8, so on every support box the triple
      stays strictly ordered (x1 < x2 < x3, y1 < y2 < y3) and the
      paired r = 3 closed-form raw gap is a well-defined real number
      for every point of the box product with masses in the simplex.
  (2) EXACT DESIGN VALUES.  The sweep's 6 DESIGN points (CENTER =
      (1/5, 7/10) stepped by the exact rational 1/400: a cross
      pattern) are exact rationals strictly inside the mass simplex;
      at every real fiber (s, q) of the box
      the 6 values V_i(s,q) = gap(s, q, D_i) are exact reals,
      enclosed RIGOROUSLY by interval evaluation (evaluate_arb at the
      box balls: support boxes, q box, exact rational masses D_i, the
      certified beta ball).  evaluate_arb's ball-input contract makes
      each returned ball enclose the exact value at EVERY point of the
      input intervals (monotone corner brackets on the explicit H-atom
      terms plus a monotonically-bracketed h-map enclosure; no
      unverified unions anywhere in this module).
  (3) CRAMER EXTRACTION IDENTITY.  MC = the exact rational 6x6
      Vandermonde-type matrix, rows [1, a1, a2, a1a2, a1^2, a2^2] at
      D, IDENTICAL for every fiber and every box (the design is
      box-independent), exactly invertible over Q (the sweep's
      extraction).  Main's structure theorem says gap(.,.) is an
      exact quadratic in (a1,a2) at EVERY frozen fiber, hence the 6
      design values of the fiber satisfy MC A(s,q) = V(s,q) EXACTLY;
      so A(s,q) = MC V(s,q) componentwise, an EXACT identity per
      fiber.
  (4) COEFFICIENT BAND.  Combining (2)+(3): the rigorous interval
      product MC x (V boxes) is, per coefficient, an Arb ball
      containing {A_j(s,q) : (s,q) in the box}; the six balls form the
      certified BAND.  Ball arithmetic is inclusion-monotone, so each
      band ball also contains the CENTER fiber's exact coefficient
      vector A0 (the cell's Cramer coefficients).
  (5) PER-FIBER FEASIBILITY INCLUSION (the polygon widening).  Fix a
      box point (s, q) and masses a.  The mean over the box satisfies,
      with (dc0, dc1, dc2) the EXACT rational drift bracket (below),
        |mean(s,q,a) - mean(center,a)| <= dc0 + dc1 a1 + dc2 a2
                                       <= drift_total =: DT
      on the mass square [0,1]^2 (mean is multilinear in
      (q, b0..b5); per-coefficient envelopes attained at box corners,
      enumerated exactly).  Hence the TRUE per-point feasible set
        {a : mean(s,q,a) >= m}  subset of
        {a : mean(center,a) >= m - DT}  subset of
        P_box := {a1>=0, a2>=0, a1+a2<=1} AND
                 {c0 + c1 a1 + c2 a2 >= FLOOR_BOX},
      with c = mean_line(q0, geom_center) the frozen center line and
      FLOOR_BOX = FLOOR_HI - DT an exact rational: the last inclusion
      because m - DT > FLOOR_HI - DT = FLOOR_BOX (certified FLOOR_HI
      < m, strict).  P_box is a certified SUPERSET of every true
      per-point feasible set over the box; compared to the cell's
      FLOOR_HI polygon it is WIDENED by the tracked O(h) drift DT.
  (6) THE CHAIN OF INEQUALITIES.  For every box point (s,q):
        inf_a gap >= min over P_box of quad_{A(s,q)}      [per-fiber
            identity: Main's structure theorem, exact]
                  >= min over P_box of quad_{A0}
                     - max over P_box |quad_{A(s,q)} - quad_{A0}|
      and quad_A(x) - quad_B(x) is LINEAR in the coefficients, so over
      P_box subset [0,1]^2 where every monomial feature satisfies
      max |e_j| = 1 EXACTLY (1, a1, a2, a1a2, a1^2, a2^2 on [0,1]^2):
        min over P_box of quad_A >= c0lo - sum_j |A_j - A0_j|
      with c0lo the certified lower bound of the KKT minimum of the
      EXACT center quadratic A0 over P_box (re-derived per (cell,
      width) by the sweep's exact candidate list: polygon vertices +
      interior critical + edge criticals + mean-chord critical; min by
      value lower bound, deterministic name tiebreak; candidates Ball-
      certified at exact/ball points).  Finally |A_j - A0_j| <=
      width(band_j) for every fiber A (step 4), so maximizing over the
      band:
        BOXBOUND := c0lo - sum_j width(band_j)
      is a certified LOWER bound of the infimum of the raw gap over
      the ENTIRE box product (all supports in the boxes, all q in the
      q box, all truly-feasible masses).  BOXBOUND >= 0 (rigorous Arb
      comparison on the bound ball's lower endpoint) certifies the
  (7) WHY A0 LIES IN THE BAND, i.e. why the widths are valid slack.
      The center point belongs to the box, so each center design value
      gap(center, D_i) lies in V_i; MC A0 = V(center) EXACTLY (step 3)
      and the band ball is an inclusion-monotone enclosure of
      MC x (V boxes) (step 4), hence A0 in the band componentwise and
      |A_j(s,q) - A0_j| <= upper_j - lower_j for every fiber.  This is
      the ONLY slack the relaxation adds beyond interval rounding: the
      band widths (O(h) in the box half-width, plus q-box and beta
      width and the 320-bit rounding), SUMMED over the six
      coefficients - a rigorous over-estimate of the true interval-QP
      coupling, sound one-sided.

MEASURE COVERED.  Each cleared (cell, w) claims the box neighborhood;
volume per cell: (2h)^6 * qlen, qlen = 2h for interior q0 and h for
q0 in {0, 1} (clipped).  The report gives the grid fraction
cleared/nonvacuous (and /11,025) and the absolute covered-volume sum
with its fraction of the ambient [0,1]^7 - with the honest note that
adjacent cleared boxes may overlap, so the sum is an UPPER estimate of
the union; the disjoint union is NOT computed here.

LIMITATIONS (restated verbatim in the report):
  * BETWEEN-BOX GAPS: the union of certified boxes does not fill the
    ambient parameter space; between adjacent boxes (up to ~1/8 per
    support direction) the gap remains OPEN in this module.
  * r >= 4 untouched.
  * The certificate speaks of the WIDENED per-cell polygon P_box
    (superset of the true feasible sets); the 1,278 infeasible-geometry
    cells certify nothing (empty FLOOR_HI polygon at the cell center),
    exactly as in the sweep.
  * The band relaxation wastes slack (band widths are SUMMED; the
    joint interval-QP minimum over the band is not solved) - a
    rigorous over-estimate of the relaxation loss, sharpness only.
  * q-endpoint columns (q0 in {0, 1}) use one-sided q boxes of
    width h.
  * beta is a certified ball, not a point; the ambiguity propagates
    rigorously through the value boxes (evaluate_arb ball-input
    contract).

MUTATIONS (at the worst cell g9|q8; both must FAIL, mirroring the
sweep's mutation semantics):
  (a) coef_perturbed_1e-6: perturb the center quadratic's A22 by
      +1e-6; the value at the same exact mean-face point must MOVE
      and the KKT re-solve's certified min must DIFFER.
  (b) drop_mean_chord_constraint: re-solve over the bare triangle;
      the certified minimum must DIP below the constrained min.

BYTE STABILITY.  No randomness, no timestamps; every report quantity
derives from exact rationals and certified Arb balls at fixed
precision.  Writes exactly uc/verification/results/
liu9-h2-box-cover.json; two consecutive runs are byte-identical.
"""

from __future__ import annotations

import os
import sys

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_name, "1")


import hashlib
import json
import mpmath
from fractions import Fraction
from pathlib import Path
from typing import Optional, Sequence

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from flint import arb, ctx                     # noqa: E402
ctx.prec = 320

from liu9_binding import (                     # noqa: E402
    certify_equation_parameters, solve_equation_parameters)
from liu9_boundary_layer import gap_mp         # noqa: E402
from liu9_witness_qp import _pt, exact_inverse  # noqa: E402
from liu9_objective import evaluate_arb        # noqa: E402
from liu9_h2_qp_sweep import (                 # noqa: E402
    FLOOR_HI, MEAN_INDEX, build_constraints, certify_floors,
    extract_quadratic, kkt_candidates, kkt_min, mean_line,
    polygon_vertices, quad_eval, solve_cell, DESIGN)

OUTPUT_DEFAULT = HERE / "verification/results/liu9-h2-box-cover.json"
SWEEP_REPORT = HERE / "verification/results/liu9-h2-qp-sweep.json"
SWEEP_PY_PATH = HERE / "liu9_h2_qp_sweep.py"
SWEEP_PY_SHA256 = ("66305dedf82389597ef082e4d8fa6f6b0b9420a76847c53c9"
                   "cf716f0ec65b5ac")

WIDTH_LADDER = (10, 12, 14)          # half-width exponents, r3 discipline
H_HALF = {w: Fraction(1, 2 ** w) for w in WIDTH_LADDER}
COEF_ORDER = ("A0", "A1", "A2", "A12", "A11", "A22")

WORST_CELL_ID = "g9|q8"
WORST_CELL_GEOM = (Fraction(1, 8), Fraction(1, 4), Fraction(3, 8),
                   Fraction(1, 8), Fraction(1, 2), Fraction(5, 8))
WORST_CELL_Q = Fraction(1)


# ---------------------------------------------------------------------------
# Interval helpers.  Exact-rational endpoint discipline: endpoints are
# exact rationals rendered as point balls; the union of TWO POINT balls
# is an exact hull.  NO wide-ball unions anywhere in this module.

def _iball(lo: Fraction, hi: Fraction) -> arb:
    lo = Fraction(lo)
    hi = Fraction(hi)
    if lo > hi:
        raise ValueError("empty interval")
    if lo == hi:
        return _pt(lo)
    return _pt(lo).union(_pt(hi))


# ---------------------------------------------------------------------------
# Box geometry helpers.

def box_geometry(geom: tuple[Fraction, ...], h: Fraction
                 ) -> tuple[tuple[Fraction, Fraction], ...]:
    """Support boxes (g_i - h, g_i + h); grid supports are interior
    (1/8 .. 7/8 with h <= 2^-10 < 1/8), so no clipping."""
    return tuple((g - h, g + h) for g in geom)


def box_q(q0: Fraction, h: Fraction) -> tuple[Fraction, Fraction]:
    """q box (q0 - h, q0 + h) clipped to [0, 1]: one-sided at the
    endpoint columns q0 in {0, 1}."""
    return (max(Fraction(0), q0 - h), min(Fraction(1), q0 + h))


def box_volume(geom_box, q_box) -> Fraction:
    vol = Fraction(1)
    for lo, hi in geom_box:
        vol *= hi - lo
    return vol * (q_box[1] - q_box[0])


def mean_drift(geom_box, q_box) -> tuple[Fraction, Fraction, Fraction]:
    """Exact rational drift bracket of the mixture-mean line
    coefficients over the box.  mean(s,q,a) = c0 + c1 a1 + c2 a2 with
      c0 = (1-q) b4 + q b5,
      c1 = (1-q)(b0-b4) + q(b1-b5),
      c2 = (1-q)(b2-b4) + q(b3-b5),
    Each c_j is multilinear in (q, its two supports); its JOINT
    envelope over the box is attained at a corner (multilinear =>
    extremes on vertices), so the exact corner enumeration returns
    the exact envelope.  The drift triple satisfies, for EVERY
    (s,q) in the box and every (a1,a2) in [0,1]^2,
      |c_j(s,q) - c_j(center)| <= dc_j,
      |mean(s,q,a) - mean(center,a)| <= dc0 + dc1 a1 + dc2 a2
      <= dc0 + dc1 + dc2  (a1, a2 >= 0 in the triangle,
      the last step by a1 + a2 <= 1)."""
    (b0lo, b0hi), (b2lo, b2hi), (b4lo, b4hi) = geom_box[:3]
    (b1lo, b1hi), (b3lo, b3hi), (b5lo, b5hi) = geom_box[3:]
    qlo, qhi = q_box
    # Each c_j is multilinear in its variables (q, supports), so its
    # joint envelope over the box is attained at a corner; exact
    # corner enumeration (8 for c0, 32 for c1/c2) IS the exact joint
    # envelope.
    def c0_at(q, b4, b5):
        return (1 - q) * b4 + q * b5

    def c1_at(q, b0, b4, b1, b5):
        return (1 - q) * (b0 - b4) + q * (b1 - b5)

    c0_vals = [c0_at(q, b4, b5)
               for q in (qlo, qhi)
               for b4 in (b4lo, b4hi)
               for b5 in (b5lo, b5hi)]
    c1_vals = [c1_at(q, b0, b4, b1, b5)
               for q in (qlo, qhi)
               for b0 in (b0lo, b0hi)
               for b4 in (b4lo, b4hi)
               for b1 in (b1lo, b1hi)
               for b5 in (b5lo, b5hi)]
    c2_vals = [c1_at(q, b2, b4, b3, b5)
               for q in (qlo, qhi)
               for b2 in (b2lo, b2hi)
               for b4 in (b4lo, b4hi)
               for b3 in (b3lo, b3hi)
               for b5 in (b5lo, b5hi)]

    return (max(c0_vals) - min(c0_vals),
            max(c1_vals) - min(c1_vals),
            max(c2_vals) - min(c2_vals))

# ---------------------------------------------------------------------------


def run_mutations(beta: arb, geom: tuple[Fraction, ...], q0: Fraction,
                  params_mean: arb) -> list[dict]:
    """The sweep's two mission mutations, on the worst cell's LIVE
    center data; both must FAIL (their absence of movement would mean
    the certificate machinery is insensitive to its own inputs)."""
    center = solve_cell(beta, geom, q0, FLOOR_HI, params_mean)
    coefd = {k: arb(v) for k, v in center["coefficients"].items()}
    constraints = center["_constraints"]
    # mutation (a): perturb A22, value must move at the mean-face
    # midpoint AND the KKT re-solve must differ
    verts = polygon_vertices(constraints)
    a_m, b_m, c_m = constraints[MEAN_INDEX]
    face = sorted({v[0] for v in verts
                   if a_m * v[0] + b_m * v[1] == c_m})
    x_ref = (face[0] + face[-1]) / 2
    y_ref = (c_m - a_m * x_ref) / b_m
    v_ref = quad_eval(coefd, _pt(x_ref), _pt(y_ref))
    pert = dict(coefd)
    pert["A22"] = coefd["A22"] + _pt(Fraction(1, 10 ** 6))
    v_mut = quad_eval(pert, _pt(x_ref), _pt(y_ref))
    moved_pt = not bool((v_mut - v_ref).contains(0))
    cands = kkt_candidates(pert, constraints)
    best_p, _nc, _nf = kkt_min(cands)
    moved_k = best_p is not None and not bool(
        (best_p["_value"] - center["_min_value_arb"]).contains(0))
    mut_a = {
        "mutation": "coef_perturbed_1e-6",
        "cell": WORST_CELL_ID,
        "mechanism": ("A22 -> A22 + 1e-6 on the exact center "
                      "quadratic: the value at the same exact "
                      "mean-face point must move and the KKT "
                      "re-solve min must differ"),
        "expected": "FAIL",
        "observed": "FAIL" if (moved_pt and moved_k) else "PASSED",
        "reference_point": [str(x_ref), str(y_ref)],
        "value_at_ref": v_ref.str(20),
        "value_mutated_at_ref": v_mut.str(20),
        "value_mutated_ksolve": None if best_p is None
        else best_p["value"],
        "value_reference_ksolve": center["min"],
        "ok": bool(moved_pt and moved_k),
    }
    lean = [c for i, c in enumerate(constraints) if i != MEAN_INDEX]
    cands2 = kkt_candidates(coefd, lean, mean_index=None)
    best2, _nc2, _nf2 = kkt_min(cands2)
    dipped = best2 is not None and bool(
        best2["_value"].lower() < center["_min_value_arb"].lower())
    mut_b = {
        "mutation": "drop_mean_chord_constraint",
        "cell": WORST_CELL_ID,
        "mechanism": ("certify over the bare triangle (no mean floor "
                      "chord); the minimum must dip below the same "
                      "cell's constrained certified min"),
        "expected": "FAIL (dips below the constrained certified min)",
        "observed": "FAIL" if dipped else "PASSED",
        "triangle_only_min": None if best2 is None
        else best2["_value"].str(26),
        "triangle_only_arg": None if best2 is None
        else [best2["a1"], best2["a2"]],
        "constrained_min_lower": center["min_lower"],
        "constrained_min_upper": center["min_upper"],
        "ok": dipped,
    }
    return [mut_a, mut_b]


# ---------------------------------------------------------------------------
# Value boxes at the 6 DESIGN points, coefficient band, per-cell
# certificate.

def design_matrix() -> list[list[Fraction]]:
    """The exact rational 6x6 matrix MC: rows
    [1, a1, a2, a1a2, a1^2, a2^2] at the sweep's 6 DESIGN points -
    the SAME DESIGN list used by design_value_boxes, so the matrix
    matches the evaluation points by construction."""
    return [[Fraction(1), d1, d2, d1 * d2, d1 * d1, d2 * d2]
            for d1, d2 in DESIGN]


def design_value_boxes(beta: arb, geom_box, q_box) -> list[arb]:
    """6 certified value boxes V_i = gap(box balls, D_i)."""
    gballs = tuple(_iball(lo, hi) for lo, hi in geom_box)
    qball = _iball(q_box[0], q_box[1])
    out = []
    for d1, d2 in DESIGN:
        terms = evaluate_arb((_pt(d1), _pt(d2), qball) + gballs, beta)
        out.append(terms.numerator - terms.ehx)
    return out


def center_coefs(beta: arb, geom: tuple[Fraction, ...], q0: Fraction
                 ) -> dict[str, arb]:
    """The cell's exact center-quadratic coefficient balls: the
    sweep's own 6-point extraction at (a1, a2) = CENTER (a fresh
    extract_quadratic run; the DESIGN-extraction identity means this
    is MC x V(center), the band's center fiber)."""
    return extract_quadratic(beta, q0, geom)["_coefd"]


def cramer_band(v_boxes: list[arb],
                inv: list[list[Fraction]]) -> dict[str, arb]:
    """Per-coefficient rigorous band balls: MC x (V boxes), each a sum
    of exact-point x interval-ball products (inclusion-monotone ball
    arithmetic at 320 bits).  Each band ball contains the coefficient
    A_j(s,q) of EVERY fiber of the box (used only for the mutated /
    diagnostic band accounting; the certificate itself is the
    value-bracket radius)."""
    band = {}
    for j, name in enumerate(COEF_ORDER):
        acc = arb(0)
        for i in range(6):
            acc = acc + _pt(inv[j][i]) * v_boxes[i]
        band[name] = acc
    return band


def band_width_sum(band: dict[str, arb]) -> arb:
    acc = arb(0)
    for name in COEF_ORDER:
        acc = acc + (band[name].upper() - band[name].lower())
    return acc


def _cand_frac(s: str) -> Fraction:
    """Exact rational midpoint of a candidate's ball-string."""
    return Fraction(s.split("+/-")[0].strip("[]").strip())


def _project_simplex(a1: Fraction, a2: Fraction
                     ) -> tuple[Fraction, Fraction]:
    """Exact-rational projection onto the OPEN simplex face set that
    evaluate_arb's semantic brackets require: every component mass
    (a1, a2, and the derived a3 = 1 - a1 - a2) must be quoted at an
    exact rational whose 320-bit ball is STRICTLY inside [0, 1] on
    the relevant side.  If a1 + a2 >= 1, pull a2 to 1 - a1 - 2^-90;
    if a component is exactly 0, nudge it to 2^-90 (the KKT candidate
    radii in this run are < 1e-16 >> 2^-90, and every nudge is an
    EXACT rational, so the value box remains a rigorous enclosure of
    gap AT THE EXACT POINT EVALUATED - the point itself is what
    changes, never the rigor)."""
    eps = Fraction(1, 2 ** 90)
    if a1 <= 0:
        a1 = eps
    if a2 <= 0:
        a2 = eps
    if a1 + a2 >= 1:
        a2 = 1 - a1 - eps
        if a2 <= 0:
            a2 = eps
            a1 = 1 - eps - a2
    return a1, a2



def _m_lower_exact(params_arb) -> Fraction:
    """Exact rational BELOW the certified m lower endpoint: midpoint
    string of m.lower() at 50 digits (radius 2.42e-51 there), minus
    1e-42 slack."""
    txt = params_arb.mean.lower().str(50).split("+/-")[0]
    return Fraction(txt.strip("[]").strip()) - Fraction(1, 10 ** 42)


def radius_certificate(beta: arb, geom: tuple[Fraction, ...],
                       q0: Fraction, center_c, coefd: dict[str, arb],
                       h: Fraction, m_lo: Fraction) -> dict:
    """VALUE-BRACKET RADIUS certificate, final sound form (the
    mission's method (a)-(c) combined; full statement in
    SOUNDNESS_ARGUMENT (i)-(ix) plus the TILTED-CHORD note).

    Per (cell, h):
      1. drift over the box: the EXACT rational triple (dc0, dc1,
         dc2) with |mean(s,q,a) - mean(center,a)| <= dc0 + dc1 a1 +
         dc2 a2 for every box fiber and every (a1, a2) >= 0.
      2. TILTED-CHORD superset polygon: every true fiber-feasible a
         (mean(s,q,a) >= m) satisfies
             (c1+dc1) a1 + (c2+dc2) a2 >= m_lo - (c0+dc0),
         a single EXACT half-space; P = triangle AND that half-space
         is a certified SUPERSET of every true per-point feasible
         set - and UNIFORMLY TIGHT (tilt O(h), not the full chord
         shift of the naive floor lowering).
      3. KKT on the EXACT center quadratic over P (the sweep's
         candidate list, re-run): certified minimum ball c0 with
         argmin (a1*, a2*) in P.
      4. value box A = evaluate_arb(support box balls, q box ball,
         exact rational (a1*, a2*)): a rigorous enclosure of
         gap(s, q, a1*, a2*) over ALL fibers of the box.
      5. R = max(A.hi - c0.lower(), c0.lower() - A.lo);
         BOXBOUND = c0.lower() - R; passed iff BOXBOUND >= 0.

    Why sound (the chain): for every box fiber (s, q) and every
    fiber-feasible mass a:  a in P (step 2), so
      gap(s,q,a) = quad_{A(s,q)}(a)  [structure theorem, exact]
                 >= min_P quad_{A(s,q)}  [P superset of fiber set]
                 >= quad_{A(s,q)}(a1*,a2*) - R
                    [the fiber quadratic's min over P is <= its
                     value at the P-point (a1*,a2*), and the fiber
                     center-data oscillation |gap(s,q,x) - c0| at
                     the FIXED masses x = (a1*,a2*) is enclosed by
                     R from the value box A around c0]
                 >= A.lo - R >= c0.lower() - R.
    Infimum over box fibers proves the certificate."""
    geom_box = box_geometry(geom, h)
    q_box = box_q(q0, h)
    drift = mean_drift(geom_box, q_box)
    dc0, dc1, dc2 = drift
    c0c, c1c, c2c = center_c
    tilt1 = c1c + dc1
    tilt2 = c2c + dc2
    rhs = m_lo - c0c - dc0
    constraints = [[Fraction(1), Fraction(0), Fraction(0)],
                   [Fraction(0), Fraction(1), Fraction(0)],
                   [Fraction(-1), Fraction(-1), Fraction(-1)],
                   [tilt1, tilt2, rhs]]
    cands = kkt_candidates(coefd, constraints)
    best, ncand, nfeas = kkt_min(cands)
    if best is None:
        raise RuntimeError("tilted polygon empty at KKT solve")
    c0lo = best["_value"].lower()
    a1x, a2x = _project_simplex(_cand_frac(best["a1"]),
                                _cand_frac(best["a2"]))
    gballs = tuple(_iball(lo, hi) for lo, hi in geom_box)
    qball = _iball(q_box[0], q_box[1])
    terms = evaluate_arb((_pt(a1x), _pt(a2x), qball) + gballs, beta)
    raw = terms.numerator - terms.ehx
    r_hi = raw.upper() - c0lo
    r_lo = c0lo - raw.lower()
    radius = r_hi if r_hi > r_lo else r_lo
    bound = c0lo - radius
    passed = bool(bound.lower() >= 0)
    return {
        "geom_box": [[str(lo), str(hi)] for lo, hi in geom_box],
        "q_box": [str(q_box[0]), str(q_box[1])],
        "mean_drift": [str(v) for v in drift],
        "drift_total": str(dc0 + dc1 + dc2),
        "chord_constraint": [str(tilt1), str(tilt2), str(rhs)],
        "kkt_src": best["name"], "kkt_min": best["value"],
        "kkt_arg": [best["a1"], best["a2"]],
        "n_candidates": ncand, "n_feasible": nfeas,
        "c0_lower": c0lo.str(26),
        "value_box_lo": raw.lower().str(16),
        "value_box_hi": raw.upper().str(16),
        "radius": radius.str(16),
        "radius_f": float(radius.upper()),
        "bound": bound.str(16),
        "bound_f": float(bound.lower()),
        "passed": passed,
        "_bound": bound,
    }


SOUNDNESS_ARGUMENT = """RELAXATION (value-bracket radius + widened-KKT
bound; the certified claim of this module).  (i) At every real fiber
(s, q) of the closed box the support triples stay strictly ordered
(h <= 2^-10 < 1/16 < min gap 1/8), so the paired r = 3 raw gap is
a well-defined real F(s, q; a1, a2).  (ii) Main's structure theorem
(PROVED, cited): at EVERY frozen fiber the gap is an EXACT quadratic
in (a1, a2).  (iii) WIDENED POLYGON: the mixture mean is exact
rational-multilinear; the drift bracket (mean_drift) bounds
|mean(s,q,a) - mean(center,a)| <= DT on the mass square; FLOOR_BOX =
FLOOR_HI - DT (exact rational) with FLOOR_HI < m certified makes the
frozen polygon P_box = {a1>=0, a2>=0, a1+a2<=1, center line >=
FLOOR_BOX} a SUPERSET of every true per-point feasible set {mean >=
m} over the box.  (iv) KKT OVER THE WIDENED POLYGON: the sweep's
exact candidate list (vertices, interior critical, edge criticals,
mean-chord critical) re-solved on the center quadratic A0 over P_box
gives the certified minimum ball c0_box with argmin (a1*, a2*) in
P_box, feasible within the widened discipline.  (v) VALUE BOX AT THE
ARGMIN: evaluate_arb at the box balls (six support balls, q ball) and
the exact rational masses (a1*, a2*) returns a rigorous enclosure
[A.lo, A.hi] of {gap(s, q, a1*, a2*) : (s, q) in box} (ball-input
contract: monotone corner brackets on the H-atom terms plus a
monotonically-bracketed h-map; no unverified unions).  (vi) THE
CHAIN: for every box fiber (s, q) and every mass a FEASIBLE for the
fiber (mean(s,q,a) >= m):
  gap(s,q,a) = quad_{A(s,q)}(a)      [structure theorem, exact]
             >= min_P_box quad_{A(s,q)}
                        [fiber-feasible a in P_box by (iii)]
             >= quad_{A(s,q)}(a1*, a2*) - R
                        [P_box-argmin optimality of the fiber
                         quadratic: min <= value at any P_box point]
             >= A.lo - R             [(v) value box]
             >= c0_box.lower() - R,
with the RADIUS R = max(A.hi - c0_box.lower(), c0_box.lower() - A.lo)
(a rigorous ball difference; the max over the two directions covers
the value-box uncertainty on both sides).  Taking the infimum over
box fibers and fiber-feasible masses:
  inf gap >= c0_box.lower() - R =: BOXBOUND.
(vii) PROMOTION: if the P_box argmin (the midpoint rational of the
certified candidate ball) is FEASIBLE AT THE FIBER (mean(s,q,a*) >=
m), then quad_{A(s,q)}(a*) = gap(s, q, a*) lies in the SAME value
box A, and the chain additionally gives
  min over fiber-feasible a of gap(s,q,a) >= A.lo - R >= BOXBOUND.
For fibers where a* is not fiber-feasible, step (vi) still holds
through the widened polygon: P_box is a certified superset of the
fiber-feasible set, so min_P_box quad_{A(s,q)} <= min over the
fiber-feasible set, hence the inequality direction is PRESERVED
(lower-bounding).  The certified quantity is therefore exactly the
mission's positive-measure statement over each box.  (viii) The
VACUOUS columns are excluded (empty FLOOR_HI polygon at the cell
center), exactly as in the sweep; nothing uses grid structure beyond
the box endpoints being exact rationals.  (ix) Sharpness note: R is
first-order in h (the gap variation over the box at fixed masses),
so BOXBOUND = (cell certified margin) - O(h): the width ladder marks
where the oscillation exhausts the margin."""

MEASURE_NOTE = ("Per cleared (cell, w): box volume (2h)^6 * qlen, "
                "qlen = 2h for interior q0, h for q0 in {0, 1} "
                "(clipped).  The absolute covered-volume sum "
                "OVERCOUNTS the union when adjacent boxes overlap; "
                "the true union fraction is NOT computed here "
                "(disjointization would need cross-cell boundary "
                "bookkeeping over 1225 geometries x 9 columns).  The "
                "grid fractions cleared/9,747 and cleared/11,025 are "
                "exact counts.")

LIMITATIONS_TEXT = [
    "BETWEEN-BOX GAPS: the cover is the union of 9,747 box "
    "neighborhoods per width; the ambient parameter space between "
    "adjacent boxes (up to ~1/8 minus two half-widths per support "
    "direction, ~1/8 minus two half-widths in q) is NOT covered by "
    "this module.",
    "r >= 4 (more than three supports per side) is untouched here.",
    "The certificate speaks of the WIDENED per-cell polygon P_box "
    "(FLOOR_BOX = FLOOR_HI - drift_total, exact rational), a "
    "certified superset of the true per-point feasible sets; the "
    "1,278 infeasible-geometry cells certify nothing (empty "
    "FLOOR_HI polygon at the cell center), exactly as in the sweep.",
    "The band relaxation wastes slack: the certificate now uses the "
    "value-bracket radius (a single point evaluation per box), but "
    "the band accounting module (cramer_band) remains a diagnostic; "
    "the radius is a rigorous over-estimate of the oscillation "
    "affecting sharpness only, never soundness.",
    "q-endpoint columns (q0 in {0, 1}) use one-sided q boxes of "
    "width h after clipping to [0, 1].",
    "beta is the certified ball (interval), not a point; the "
    "ambiguity propagates rigorously through the value boxes "
    "(evaluate_arb ball-input contract).",
    "This is a NEIGHBORHOOD certificate around 9,747 certified grid "
    "cells, NOT a continuum theorem between cells; the absolute "
    "covered-volume sum is an upper estimate (overlapping boxes "
    "double-counted).",
]


def _canonical_digest(payload: dict) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def build_report() -> dict:
    params_mp = solve_equation_parameters(100)
    params_arb = certify_equation_parameters(params_mp)
    beta = params_arb.beta
    floors = certify_floors(params_arb)
    m_lo = _m_lower_exact(params_arb)

    raw = json.loads(SWEEP_REPORT.read_text())
    sweep_meta = raw
    sweep_cells = raw["cells"]
    nonvac = [c for c in sweep_cells if c["nonvacuous"] == "true"]
    assert len(sweep_cells) == 11025
    assert len(nonvac) == 9747

    widths = {w: {"cleared": 0, "of": len(nonvac), "failed_sample": []}
              for w in WIDTH_LADDER}
    worst_bound_f = None
    worst_detail = None
    all_failed_cells = []
    width_worst = {}
    for row in nonvac:
        geom = tuple(Fraction(v) for v in row["geometry"])
        q0 = Fraction(row["q"])
        # live center solve (the sweep's exact machinery), kept as the
        # reproduction sanity: the certificate itself re-solves KKT on
        # the TILTED polygon per width.
        center = solve_cell(beta, geom, q0, FLOOR_HI, params_arb.mean)
        if center["status"] != "ok":
            raise RuntimeError("live center solve not ok at %s"
                               % row["cell"])
        sweep_lower = float(row["min_lower_f"])
        live_lower = float(center["_min_value_arb"].lower())
        tol = 4.1e-6 * max(1.0, abs(sweep_lower))
        if abs(sweep_lower - live_lower) > tol:
            raise RuntimeError("live/sweep min mismatch at %s: %r vs %r"
                               % (row["cell"], sweep_lower, live_lower))
        coefd = center_coefs(beta, geom, q0)
        c = mean_line(q0, geom)
        cell_failed_all = True
        for w in WIDTH_LADDER:
            cert = radius_certificate(beta, geom, q0, c, coefd,
                                      H_HALF[w], m_lo)
            if cert["passed"]:
                widths[w]["cleared"] += 1
                cell_failed_all = False
            elif len(widths[w]["failed_sample"]) < 4:
                widths[w]["failed_sample"].append({
                    "cell": row["cell"], "q": row["q"],
                    "bound": cert["bound"], "bound_f": cert["bound_f"]})
            bf = cert["bound_f"]
            if worst_bound_f is None or bf < worst_bound_f:
                worst_bound_f = bf
                worst_detail = {
                    "cell": row["cell"], "q": row["q"], "width": w,
                    "bound": cert["bound"], "bound_f": cert["bound_f"],
                    "c0_lower": cert["c0_lower"],
                    "radius": cert["radius"],
                    "radius_f": cert["radius_f"],
                    "value_box_lo": cert["value_box_lo"],
                    "value_box_hi": cert["value_box_hi"],
                    "kkt_src": cert["kkt_src"],
                    "kkt_arg": cert["kkt_arg"],
                    "chord_constraint": cert["chord_constraint"],
                    "mean_drift_total": cert["drift_total"],
                    "geometry": row["geometry"]}
            wkey = "2^-%d" % w
            wcur = width_worst.get(wkey)
            if wcur is None or bf < wcur["bound_f"]:
                width_worst[wkey] = {
                    "cell": row["cell"], "q": row["q"],
                    "bound": cert["bound"], "bound_f": cert["bound_f"]}
        if cell_failed_all:
            all_failed_cells.append(row)
    largest = None
    for w in sorted(WIDTH_LADDER, reverse=True):
        if widths[w]["cleared"] == len(nonvac):
            largest = w
            break

    # ESCALATION (mission rule): a cell failing at ALL widths is NOT a
    # counterexample - interval slack can push a certified bound
    # negative.  Re-examine at the tightest width (14) with gap_mp at
    # 200 dps at the box center + the cell's certified KKT argmin:
    # candidate witness (DISCOVERY) with full coordinates, flipping
    # the claim status to NEGATIVE-WITNESS-REPORTED.
    escalations = []
    for row in all_failed_cells:
        geom = tuple(Fraction(v) for v in row["geometry"])
        q0 = Fraction(row["q"])
        center = solve_cell(beta, geom, q0, FLOOR_HI, params_arb.mean)
        beta_mid = mpmath.mpf(
            params_arb.beta.midpoint().str(50).split("+/-")[0]
            .strip("[]").strip())
        with mpmath.workdps(200):
            a1 = mpmath.mpf(
                str(center["min_arg"][0]).split("+/-")[0]
                .strip("[]").strip())
            a2 = mpmath.mpf(
                str(center["min_arg"][1]).split("+/-")[0]
                .strip("[]").strip())
            qm = mpmath.mpf(q0.numerator) / mpmath.mpf(q0.denominator)
            vec = (a1, a2, qm) + tuple(
                mpmath.mpf(Fraction(v).numerator)
                / mpmath.mpf(Fraction(v).denominator) for v in geom)
            g = gap_mp(vec, beta_mid)
        escalations.append({
            "cell": row["cell"], "q": row["q"],
            "geometry": row["geometry"],
            "masses_a1_a2_argmin": center["min_arg"],
            "certified_center_min": center["min"],
            "box_bounds_all_widths_negative": True,
            "gap_mp_200dps": mpmath.nstr(g, 30),
            "gap_mp_negative": bool(g < 0),
            "discovery_only": True,
        })
    mutations = run_mutations(beta, WORST_CELL_GEOM, WORST_CELL_Q,
                              params_arb.mean)

    # measure accounting at the largest fully-clearing width
    covered = Fraction(0)
    if largest is not None:
        h = H_HALF[largest]
        for row in nonvac:
            geom = tuple(Fraction(v) for v in row["geometry"])
            q0 = Fraction(row["q"])
            gbox = box_geometry(geom, h)
            qbox = box_q(q0, h)
            covered += box_volume(gbox, qbox)
    covered_frac = covered / 1

    # claim status
    all_clear = largest is not None
    muts_ok = all(m["ok"] for m in mutations)
    if any(not e["gap_mp_negative"] for e in escalations):
        claim_status = "FAILED"
    elif escalations:
        claim_status = "NEGATIVE-WITNESS-REPORTED"
    elif all_clear and muts_ok:
        claim_status = "ALL-CELLS-BOX-COVER-CERTIFIED"
    else:
        claim_status = "FAILED"

    report = {
        "tool": "uc/liu9_h2_box_cover.py",
        "claim": ("BOX COVER: every nonvacuous cell of the H2 mass-QP "
                  "sweep grid (9,747 of 11,025; support triples from "
                  "the k/8 grid, q from the k/8 grid) certifies gap "
                  ">= 0 on the BOX {supports +- 2^-w} x {q +- 2^-w "
                  "clipped} x {masses feasible for the TRUE mean at "
                  "that box point}, at the largest half-width 2^-%d "
                  "for which ALL cells clear; per-cell certificate = "
                  "value-bracket radius around the cell's exact KKT "
                  "argmin re-solved on the tilted-chord polygon "
                  "(exact rational drift dilation of the mean line; "
                  "see relaxation.soundness_argument)" % (
                      largest if largest is not None else 0,)),
        "claim_status": claim_status,
        "claim_scope": ("per-cell box certificates: ONE certificate = "
                        "one (x-triple, y-triple, q) grid cell's "
                        "neighborhood, NOT a continuum theorem; the "
                        "union of boxes does not fill the ambient "
                        "parameter space"),
        "verdict": ("ALL %d nonvacuous cells box-certified at "
                    "half-width 2^-%d" % (len(nonvac), largest)
                    if all_clear else "PARTIAL: no width clears all "
                    "cells; see box_widths"),
        "feasibility": floors,
        "grid": {
            "lineage": "uc/liu9_h2_qp_sweep.py report (frozen; the "
                       "sweep module sha256 pinned next); 11,025 "
                       "cells, 9,747 nonvacuous, worst cell g9|q8 "
                       "certified min +1.6305e-3",
            "sweep_report_sha256": sweep_meta["report_sha256"],
            "sweep_module_sha256_pinned": SWEEP_PY_SHA256,
            "sweep_module_sha256_live": hashlib.sha256(
                SWEEP_PY_PATH.read_bytes()).hexdigest(),
            "cells_total": len(sweep_cells),
            "cells_nonvacuous": len(nonvac),
            "cells_excluded_vacuous": len(sweep_cells) - len(nonvac),
        },
        "relaxation": {
            "name": "value-bracket radius + tilted-chord KKT bound",
            "soundness_argument": SOUNDNESS_ARGUMENT,
            "mean_drift_note": ("mean is multilinear in (q, "
                                "supports); drift bracket from exact "
                                "corner enumeration; the TILTED-CHORD "
                                "half-space dilation keeps the "
                                "superset polygon uniformly tight (no "
                                "full chord shift)"),
        },
        "box_widths": {
            "2^-%d" % w: {
                "half_width": str(H_HALF[w]),
                "cleared": widths[w]["cleared"],
                "of": widths[w]["of"],
                "failed_sample": widths[w]["failed_sample"],
                "worst_bound": width_worst.get("2^-%d" % w),
            } for w in sorted(WIDTH_LADDER)
        },
        "largest_fully_clearing_half_width": (
            "2^-%d" % largest if largest is not None else None),
        "measure_covered": {
            "note": MEASURE_NOTE,
            **({"largest_fully_clearing_width": "2^-%d" % largest,
                "covered_volume_absolute": str(covered),
                "covered_volume_float": float(covered),
                "ambient_fraction": str(covered_frac),
                "ambient_fraction_float": float(covered_frac),
                "grid_fraction_nonvacuous": "%d/%d" % (
                    widths[largest]["cleared"], len(nonvac)),
                "grid_fraction_total": "%d/%d" % (
                    widths[largest]["cleared"], len(sweep_cells))}
               if largest is not None else {}),
        },
        "worst_box_certificate": worst_detail,
        "escalation": escalations,
        "mutations": mutations,
        "limitations": LIMITATIONS_TEXT,
        "evaluator": ("liu9_objective.evaluate_arb, gap = numerator - "
                      "ehx, ctx.prec = 320; KKT/candidates via "
                      "uc/liu9_h2_qp_sweep.solve_cell - exact same "
                      "code path as the frozen sweep; band via "
                      "exact_inverse Cramer products"),
    }
    report["report_sha256"] = _canonical_digest(report)
    return report


RETRACTION_NOTICE = (
    "RETRACTED 2026-09-01: this module's per-cell certificate is not a "
    "lower bound (inverted inequality in radius_certificate / "
    "SOUNDNESS_ARGUMENT step (vi): the value-bracket radius R is measured "
    "at ONE fixed mass point a*, so it cannot bound the fiber quadratic's "
    "minimum over the mass polygon). Its sound sibling "
    "c0lo - sum_j width(band_j) is vacuous at every affordable width "
    "(21651.78 at 2^-10, 1.3215 at 2^-24, against a 1.6e-3 margin). No "
    "clearing count, worst bound, or covered volume from this module may "
    "be cited. See PROGRESS.md RETRACTION-H2-BOX-COVER (2026-09-01)."
)


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    print(RETRACTION_NOTICE)
    if "--retracted-replay" not in argv:
        print("refusing to run: pass --retracted-replay to reproduce the "
              "retracted numbers for audit (the report is relabelled and "
              "no certified status is emitted)")
        return 3
    report = build_report()
    report["claim_status"] = "RETRACTED"
    report["retraction"] = RETRACTION_NOTICE
    report["claim"] = "RETRACTED - " + report["claim"]
    report["verdict"] = "RETRACTED - " + report["verdict"]
    report["report_sha256"] = _canonical_digest(report)
    OUTPUT_DEFAULT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_DEFAULT.write_text(
        json.dumps(report, indent=1, sort_keys=True) + "\n")
    print("claim_status:", report["claim_status"])
    for key, roww in sorted(report["box_widths"].items()):
        print("  box %-8s cleared %d/%d (RETRACTED, not a bound)"
              % (key, roww["cleared"], roww["of"]))
    print("report_sha256", report["report_sha256"])
    return 3


if __name__ == "__main__":
    raise SystemExit(main())
