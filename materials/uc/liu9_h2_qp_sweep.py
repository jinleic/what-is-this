#!/usr/bin/env python3
"""H2 MASS-QP SWEEP: closed-form KKT argmin per (geometry, q) cell over the
k/8 ordered grid, reusing the certified per-fiber machinery of
uc/liu9_witness_qp.py.

BACKGROUND (proved elsewhere, cited not re-proved here).  For the r = 3
paired class (P0, P1 = shared masses a on the x- and y-supports, mixed
with weight q) the raw gap splits exactly as

    gap = F(P_mix) + q(1-q) c_ch

(liu9-block-copositive-diagonalization.json, liu9-paired-class-c.json),
and at frozen rational supports Main has PROVED gap is exactly
bidegree-(2, 2) in ((a1, a2), q) over the formal H-atom ring: at a frozen
(geometry, q) the gap is an EXACT quadratic in the two free masses.
uc/liu9_witness_qp.py exploited that at ONE fiber (random-111143
supports, q = 65519/65536): 6-point exact-rational quadratic extraction
in (a1, a2), deg-2 probe verification, then a closed-form KKT candidate
list over the mass polygon; every value Arb-certified.  Its certified
per-fiber minimum was +4.5036e-5.

THIS MODULE.  Sweep the SAME per-cell statement over the product grid

  * supports: k/8 grid {1/8, ..., 7/8}; ordered triples of distinct
    values per side (C(7,3) = 35 per side), ordered (x-side, y-side)
    pairs -> 35 x 35 = 1225 geometries (the enumeration of
    uc/liu9_scalar_margin_r3.py, SUPPORT_GRID = 8);
  * q: k/8 grid {0, 1/8, ..., 1}, 9 columns INCLUDING the endpoints;

for 1225 x 9 = 11025 cells.  Per cell: exact-rational quadratic
extraction at CENTER = (1/5, 7/10) with step 1/400, 4 deg-2 probe
verifications, then the KKT candidate list over the polygon

    {a1 >= 0, a2 >= 0, a1 + a2 <= 1} AND {mixture mean >= M_FLOOR}

with candidate classes = polygon vertices (all exact rationals) +
interior 2x2 critical + per-face 1D criticals + the mean-chord
{M = M_FLOOR} slice critical (the witness's E_chord_crit formulas,
re-derived; RB's note: they are already correct there, and the
per-edge generalization used here is checked against them).  min over
feasible candidates = certified per-cell QP min; every value an Arb ball
at an exact/ball point.

SUPERSET DISCIPLINE (as the witness and M_UNDER).  M_FLOOR is the exact
rational FLOOR_HI = (certified midpoint of m's Arb enclosure) - 1e-69,
Arb-verified strictly below the certified m and above the witness's
M_UNDER.  So {mean >= M_FLOOR} is a SUPERSET of the true feasible set
{mean >= m} and the certified per-cell minimum LOWER-bounds the true
mass-worst-case gap.  Certified positive per cell => gap >= 0 on the
TRUE per-cell feasible set.

q = 0 AND q = 1 COLUMNS (exact reduction, per mission).  At q = 0 and
q = 1 the channel term q(1-q)c_ch vanishes exactly and gap = F(P0),
resp. F(P1).  The mean entering the feasible polygon is the P-mix mean,
which at q = 1 IS the P1 mean, so the q = 1 column certifies F(P1) >= 0
on {mean(P1) >= M_FLOOR} (a superset of {mean(P1) >= m}) -- the region
where the universal q = 1 theorem F >= 0 is PUBLISHED (524,800-box
centered Arb cover + dC_M/dM >= 0; GENERAL_R_H2_COMPOSITION.md sec. 1).
The QP extraction and KKT solve run unchanged on those columns.

CLAIM SCOPE (honest labels).  Each cell yields a FINITE, machine-verified
statement about ONE (x-triple, y-triple, q) point: "for all shared
masses in the triangle with mixture mean >= m, gap >= 0".  11025 such
cell claims is a finite cover of grid points, NOT a continuum theorem in
(geometry, q) -- the inter-cell gaps in support and q remain open, and
claim_status says so.  If any certified per-cell min is negative, the
run escalates: gap_mp re-verification at 200 dps at the argmin and the
exact negative witness (geometry, q, masses, gap enclosure) becomes the
top-line report (label DISCOVERY on the mpmath re-check; the Arb cell
ball itself is the certified part).

MUTATIONS (must FAIL, per mission): (a) perturb one quadratic
coefficient by +1e-6 => the KKT min must move; (b) drop the mean-chord
(mean-floor) constraint => the min must dip.  Plus a DISCOVERY-only
agreement scan (mpmath feasible-set grid vs certified min) at the worst
cell.

Runner:  math/.venv/bin/python -I -B math/uc/liu9_h2_qp_sweep.py
(one core, nice -n 19; this workstation shares the live k-collector
process -- nothing here spawns threads, and the report is byte-stable:
no timings, no randomness).  Output:
  uc/verification/results/liu9-h2-qp-sweep.json
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any, Optional

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from flint import arb, ctx  # noqa: E402

CTX_PREC = 320
ctx.prec = max(ctx.prec, CTX_PREC)

import mpmath  # noqa: E402

from liu9_binding import (  # noqa: E402
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_boundary_layer import gap_mp  # noqa: E402
from liu9_objective import evaluate_arb  # noqa: E402
from liu9_witness_qp import M_UNDER, _pt, exact_inverse  # noqa: E402

OUTPUT_DEFAULT = HERE / "verification/results/liu9-h2-qp-sweep.json"

# Certified superset floor: exact rational, M_UNDER < FLOOR_HI < m (Arb).
# The decimal is the certified midpoint of m's Arb enclosure
# (solve+certify at 100 dps) minus 1e-69; the strict comparisons are
# re-verified inside certify_floors() on every run.
M_MID_DEC = ("0.61729091208126497006968709790102738837361856749424241279"
             "90865050809197")
FLOOR_HI = Fraction(M_MID_DEC) - Fraction(1, 10 ** 69)

# Grids.
SUPPORT_GRID = 8
SUPPORT_VALUES = [Fraction(k, SUPPORT_GRID) for k in range(1, SUPPORT_GRID)]
TRIPLES = list(combinations(SUPPORT_VALUES, 3))
Q_GRID = 8
Q_VALUES = [Fraction(k, Q_GRID) for k in range(Q_GRID + 1)]

# Quadratic-extraction design (exact rationals around the cell centre;
# identical construction to the witness, finer step).
CENTER = (Fraction(1, 5), Fraction(7, 10))
H_STEP = Fraction(1, 400)
DESIGN = [
    (CENTER[0], CENTER[1]),
    (CENTER[0] + H_STEP, CENTER[1]),
    (CENTER[0] - H_STEP, CENTER[1]),
    (CENTER[0], CENTER[1] + H_STEP),
    (CENTER[0], CENTER[1] - H_STEP),
    (CENTER[0] + H_STEP, CENTER[1] + H_STEP),
]
VALIDATION_POINTS = [
    (Fraction(3, 20), Fraction(4, 5)),
    (Fraction(7, 20), Fraction(13, 20)),
    (Fraction(3, 25), Fraction(3, 5)),
    (Fraction(2, 5), Fraction(3, 5)),
]
QUAD_EPS = Fraction(1, 10 ** 40)
TOL_ZERO = Fraction(1, 10 ** 12)   # negativity tolerance on certified balls
SCAN_N = 21                        # discovery agreement scan resolution
WITNESS_CAP = 200                  # max stored negative witnesses
TOL_MPMATH = mpmath.mpf(10) ** -12  # mpmath mirror of TOL_ZERO

MEAN_INDEX = 3                     # constraint index of the mean line

Constraint = tuple[Fraction, Fraction, Fraction]


# ---------------------------------------------------------------------------
# Exact-rational cell geometry.

def geometry_pairs() -> list[tuple[Fraction, ...]]:
    """Ordered (x-triple, y-triple) pairs on the k/8 grid: the enumeration
    of liu9_scalar_margin_r3.py (35 x 35 = 1225 sorted distinct triples)."""
    return [tuple(list(xtr) + list(ytr)) for xtr in TRIPLES for ytr in TRIPLES]


def mean_line(q: Fraction, geom: tuple[Fraction, ...]) -> tuple[Fraction,
                                                                Fraction,
                                                                Fraction]:
    """mean(a1, a2) = c0 + c1 a1 + c2 a2 at the frozen (supports, q)."""
    b0, b2, b4, b1, b3, b5 = geom
    c0 = (1 - q) * b4 + q * b5
    c1 = (1 - q) * (b0 - b4) + q * (b1 - b5)
    c2 = (1 - q) * (b2 - b4) + q * (b3 - b5)
    return c0, c1, c2


def build_constraints(c: tuple[Fraction, Fraction, Fraction],
                      floor: Fraction) -> list[Constraint]:
    """The polygon as exact-rational halfplanes A a1 + B a2 >= C:
    triangle {a1 >= 0, a2 >= 0, a1 + a2 <= 1} and the mean floor."""
    c0, c1, c2 = c
    return [(Fraction(1), Fraction(0), Fraction(0)),
            (Fraction(0), Fraction(1), Fraction(0)),
            (Fraction(-1), Fraction(-1), Fraction(-1)),
            (c1, c2, floor - c0)]


def polygon_vertices(constraints: list[Constraint]
                     ) -> list[tuple[Fraction, Fraction]]:
    """Exact rational vertices of the halfplane intersection, sorted.
    A bounded intersection of halfplanes is nonempty iff it has a
    vertex; an empty return means the polygon is empty."""
    verts: list[tuple[Fraction, Fraction]] = []
    n = len(constraints)
    for i in range(n):
        for j in range(i + 1, n):
            a1, b1, c1 = constraints[i]
            a2, b2, c2 = constraints[j]
            det = a1 * b2 - a2 * b1
            if det == 0:
                continue
            x = (c1 * b2 - c2 * b1) / det
            y = (a1 * c2 - a2 * c1) / det
            if all(a * x + b * y >= c for a, b, c in constraints):
                verts.append((x, y))
    return sorted(set(verts))


def polygon_edges(constraints: list[Constraint],
                  verts: list[tuple[Fraction, Fraction]]
                  ) -> list[tuple[int, tuple[Fraction, Fraction],
                                  tuple[Fraction, Fraction]]]:
    """(constraint_index, v_lo, v_hi) consecutive vertex pairs along each
    nondegenerate face."""
    edges = []
    for k, (a, b, c) in enumerate(constraints):
        if a == 0 and b == 0:
            continue
        on = [v for v in verts if a * v[0] + b * v[1] == c]
        on.sort(key=lambda v: -b * v[0] + a * v[1])
        for lo, hi in zip(on, on[1:]):
            if lo != hi:
                edges.append((k, lo, hi))
    return edges


# ---------------------------------------------------------------------------
# Certified gap evaluation and quadratic extraction.

def gap_at(q: Fraction, geom: tuple[Fraction, ...], a1: Fraction,
           a2: Fraction, beta: arb) -> arb:
    b0, b2, b4, b1, b3, b5 = geom
    values = tuple(_pt(v) for v in (a1, a2, q, b0, b2, b4, b1, b3, b5))
    terms = evaluate_arb(values, beta)
    return terms.numerator - terms.ehx


def quad_eval(c: dict[str, arb], x: arb, y: arb) -> arb:
    return (c["A0"] + c["A1"] * x + c["A2"] * y + c["A12"] * x * y
            + c["A11"] * x * x + c["A22"] * y * y)


def extract_quadratic(beta: arb, q: Fraction,
                      geom: tuple[Fraction, ...]) -> dict[str, Any]:
    """Certified 6-point exact-rational quadratic extraction at the frozen
    (geometry, q), with deg-2 probe verification: the witness's staging,
    applied per cell."""
    values = [gap_at(q, geom, a1, a2, beta) for a1, a2 in DESIGN]
    matrix = [[Fraction(1), a1, a2, a1 * a2, a1 * a1, a2 * a2]
              for a1, a2 in DESIGN]
    inv = exact_inverse(matrix)
    coef = []
    for j in range(6):
        acc = arb(0)
        for i in range(6):
            acc = acc + _pt(inv[j][i]) * values[i]
        coef.append(acc)
    names = ("A0", "A1", "A2", "A12", "A11", "A22")
    coefd = dict(zip(names, coef))
    widths = [float(x.upper() - x.lower()) for x in coef]
    max_width = max(widths)
    eps_hi, eps_lo = _pt(QUAD_EPS), -_pt(QUAD_EPS)
    rows = []
    ok = True
    worst = 0.0
    for a1, a2 in VALIDATION_POINTS:
        direct = gap_at(q, geom, a1, a2, beta)
        fit = quad_eval(coefd, _pt(a1), _pt(a2))
        residual = fit - direct
        contains = bool(residual.contains(0))
        radius = float(residual.upper() - residual.lower())
        local_ok = (contains and residual.upper() < eps_hi
                    and residual.lower() > eps_lo)
        ok = ok and local_ok
        worst = max(worst, radius)
        rows.append({"a1": str(a1), "a2": str(a2),
                     "residual": residual.str(10),
                     "residual_radius_float": radius,
                     "contains_zero": contains, "ok": bool(local_ok)})
    return {
        "coefficients": {k: v.str(30) for k, v in coefd.items()},
        "_coefd": coefd,
        "coefficient_widths_float": widths,
        "max_coefficient_width_float": max_width,
        "degree_two_verified": bool(ok),
        "validation": rows,
        "worst_residual_radius_float": worst,
        "validation_threshold": str(QUAD_EPS),
        "ok": bool(ok and max_width < 1e-30),
    }


# ---------------------------------------------------------------------------
# Feasibility of ball points (for critical points that are not exact
# rationals): every constraint slack must have a nonnegative lower bound,
# decided at heightened precision if needed.

def _slack_ok_ball(x: arb, y: arb, constraints: list[Constraint]) -> bool:
    for prec in (CTX_PREC, 600):
        saved = ctx.prec
        ctx.prec = prec
        try:
            ok = True
            for a, b, c in constraints:
                s = _pt(a) * x + _pt(b) * y - _pt(c)
                if s.lower() < 0:
                    ok = False
                    break
            if ok:
                return True
        finally:
            ctx.prec = saved
    return False


# ---------------------------------------------------------------------------
# KKT candidate list over the polygon (the witness's Stage 3, generalized
# to the exact-rational vertex/face structure of the cell polygon).

def _cert(name: str, family: str, x: Any, y: Any, coefd: dict[str, arb],
          constraints: list[Constraint], feasible: bool, mean_face: bool,
          exact_point: bool) -> dict[str, Any]:
    xb = x if isinstance(x, arb) else _pt(x)
    yb = y if isinstance(y, arb) else _pt(y)
    value = quad_eval(coefd, xb, yb)
    return {"name": name, "family": family,
            "a1": xb.str(16), "a2": yb.str(16),
            "value": value.str(26), "value_lower": value.lower().str(26),
            "feasible": bool(feasible), "mean_face": bool(mean_face),
            "exact_point": bool(exact_point), "_value": value}


def _mean_face_vertices(constraints: list[Constraint],
                        mean_index: Optional[int]) -> frozenset[int]:
    """Indices of vertices carrying zero slack on the mean constraint."""
    if mean_index is None:
        return frozenset()
    verts = polygon_vertices(constraints)
    a, b, c = constraints[mean_index]
    return frozenset(i for i, v in enumerate(verts)
                     if a * v[0] + b * v[1] == c)


def kkt_candidates(coefd: dict[str, arb], constraints: list[Constraint],
                   mean_index: Optional[int] = MEAN_INDEX
                   ) -> list[dict[str, Any]]:
    """Every minimum of a smooth function on a compact polygon is attained
    at a vertex, at an interior critical point, or in the relative interior
    of a face at the face-restricted critical point -- the exact classes
    below (the witness's list, instantiated for its polygon, is a special
    case)."""
    cands: list[dict[str, Any]] = []
    verts = polygon_vertices(constraints)
    edges = polygon_edges(constraints, verts)
    mean_verts = _mean_face_vertices(constraints, mean_index)

    for idx, (x, y) in enumerate(verts):
        cands.append(_cert("V%d" % idx, "vertex", x, y, coefd, constraints,
                           True, idx in mean_verts, True))

    det = 4 * coefd["A11"] * coefd["A22"] - coefd["A12"] * coefd["A12"]
    if not det.contains(0):
        xc = (-coefd["A1"] * 2 * coefd["A22"]
              + coefd["A12"] * coefd["A2"]) / det
        yc = (-coefd["A2"] * 2 * coefd["A11"]
              + coefd["A12"] * coefd["A1"]) / det
        inside = _slack_ok_ball(xc, yc, constraints)
        cands.append(_cert("interior_critical", "interior", xc, yc, coefd,
                           constraints, inside, False, False))

    for k, lo, hi in edges:
        d1, d2 = hi[0] - lo[0], hi[1] - lo[1]
        x, y = _pt(lo[0]), _pt(lo[1])
        dd1, dd2 = _pt(d1), _pt(d2)
        q2 = (coefd["A11"] * dd1 * dd1 + coefd["A12"] * dd1 * dd2
              + coefd["A22"] * dd2 * dd2)
        q1 = (2 * coefd["A11"] * x * dd1
              + coefd["A12"] * (x * dd2 + y * dd1)
              + 2 * coefd["A22"] * y * dd2
              + coefd["A1"] * dd1 + coefd["A2"] * dd2)
        if q2.contains(0):
            continue                    # degenerate/linear slice: endpoints
        ustar = -q1 / (2 * q2)
        if (ustar - 1).lower() > 0 or (-ustar).lower() > 0:
            continue                    # critical outside the face
        xr, yr = x + ustar * dd1, y + ustar * dd2
        cands.append(_cert("E%d_crit" % k, "edge", xr, yr, coefd,
                           constraints, _slack_ok_ball(xr, yr, constraints),
                           k == mean_index, False))

    # Mean-chord slice candidates {M = M_FLOOR}: substitute the chord
    # a2 = alpha + beta_line a1 into the quadratic -> a quadratic in a1
    # (witness E_chord_crit coefficients s2c/s1c, re-derived; its
    # endpoints are the mean face's vertices, already V candidates).
    if mean_index is not None:
        a, b, cmean = constraints[mean_index]
        if not (a == 0 and b == 0):
            beta_line = -a / b
            alpha = cmean / b
            ab = _pt(beta_line)
            aa = _pt(alpha)
            face_a1 = sorted({v[0] for v in verts
                              if a * v[0] + b * v[1] == cmean})
            if face_a1:
                s2c = (coefd["A11"] + ab * ab * coefd["A22"]
                       + ab * coefd["A12"])
                s1c = (coefd["A1"] + coefd["A2"] * ab + coefd["A12"] * aa
                       + 2 * coefd["A22"] * aa * ab)
                if not (2 * s2c).contains(0):
                    xstar = -s1c / (2 * s2c)
                    lo1, hi1 = _pt(face_a1[0]), _pt(face_a1[-1])
                    if ((xstar - lo1).lower() >= 0
                            and (hi1 - xstar).lower() >= 0):
                        xr, yr = xstar, aa + ab * xstar
                        cands.append(_cert("MEANCHORD_crit", "mean-chord",
                                           xr, yr, coefd, constraints,
                                           _slack_ok_ball(xr, yr,
                                                          constraints),
                                           True, False))
    return cands


def kkt_min(cands: list[dict[str, Any]]
            ) -> tuple[Optional[dict[str, Any]], int, int]:
    """(argmin candidate, n_candidates, n_feasible): the witness selection
    rule -- min over feasible candidates by value lower bound (name as the
    deterministic tiebreak)."""
    insides = [c for c in cands if c["feasible"]]
    if not insides:
        return None, len(cands), 0
    best = min(insides, key=lambda c: (float(c["_value"].lower()),
                                       c["name"]))
    return best, len(cands), len(insides)


# ---------------------------------------------------------------------------
# Per-cell solve.

def solve_cell(beta: arb, geom: tuple[Fraction, ...], q: Fraction,
               floor: Fraction, mean_ball: arb) -> dict[str, Any]:
    c = mean_line(q, geom)
    constraints = build_constraints(c, floor)
    verts = polygon_vertices(constraints)
    edges = polygon_edges(constraints, verts)
    row: dict[str, Any] = {
        "geometry": [str(v) for v in geom], "q": str(q),
        "mean_line": [str(v) for v in c],
        "polygon_vertices": len(verts), "polygon_faces": len(edges),
    }
    if not verts:
        row.update({"status": "infeasible-geometry", "nonvacuous": "false"})
        return row
    extraction = extract_quadratic(beta, q, geom)
    coefd = extraction["_coefd"]
    cands = kkt_candidates(coefd, constraints)
    best, ncand, nfeas = kkt_min(cands)
    # Nonvacuity of the TRUE feasible set {mean >= m} in this cell: the
    # mean is affine, so its max over the triangle sits at a triangle
    # vertex (0,0), (1,0), (0,1); compare exact balls against the
    # certified m ball.
    tri_max = max(c[0], c[0] + c[1], c[0] + c[2])
    nv_ball = _pt(tri_max) - mean_ball
    if nv_ball.lower() > 0:
        nonvac = "true"
    elif nv_ball.upper() < 0:
        nonvac = "false"
    else:
        nonvac = "undecided"
    row.update({
        "status": "ok",
        "n_candidates": ncand, "n_feasible": nfeas,
        "min_src": None if best is None else best["name"],
        "min_family": None if best is None else best["family"],
        "min_mean_face": None if best is None else best["mean_face"],
        "min_arg": None if best is None else [best["a1"], best["a2"]],
        "min": None if best is None else best["value"],
        "min_lower": None if best is None else best["value_lower"],
        "min_upper": None if best is None
        else best["_value"].upper().str(26),
        "min_lower_f": None if best is None else float(best["_value"].lower()),
        "deg2": extraction["degree_two_verified"],
        "max_coef_width_f": extraction["max_coefficient_width_float"],
        "nonvacuous": nonvac,
        "coefficients": extraction["coefficients"],
        "_constraints": constraints,
        "_min_value_arb": None if best is None else best["_value"],
        "_candidates": cands,
    })
    return row


def certify_floors(params_arb) -> dict[str, Any]:
    return {
        "m_enclosure": params_arb.mean.str(40),
        "m_lower_str": params_arb.mean.lower().str(60),
        "beta": params_arb.beta.str(40),
        "M_UNDER": str(M_UNDER),
        "M_UNDER_lt_m_arb": bool(_pt(M_UNDER) < params_arb.mean),
        "M_FLOOR": str(FLOOR_HI),
        "M_FLOOR_lt_m_arb": bool(_pt(FLOOR_HI) < params_arb.mean),
        "M_UNDER_lt_M_FLOOR_exact": bool(M_UNDER < FLOOR_HI),
        "superset_argument": (
            "M_FLOOR < m (Arb), so each cell polygon {mean >= M_FLOOR} is "
            "a SUPERSET of the true {mean >= m}; the certified per-cell "
            "minimum lower-bounds the true mass-worst-case gap"),
    }


# ---------------------------------------------------------------------------

def _mpf_of_fraction(v: Fraction) -> mpmath.mpf:
    return mpmath.mpf(v.numerator) / mpmath.mpf(v.denominator)


def _ball_mid_text(text: str) -> mpmath.mpf:
    return mpmath.mpf(text.split("+/-")[0].strip().strip("[]"))


def escalate_negative(beta_text: str, row: dict[str, Any]) -> dict[str, Any]:
    """200-dps gap_mp re-verification at a certified-negative cell's
    argmin.  DISCOVERY evidence; the Arb cell ball is the certified part."""
    geom = tuple(Fraction(v) for v in row["geometry"])
    q = Fraction(row["q"])
    a1 = _ball_mid_text(row["min_arg"][0])
    a2 = _ball_mid_text(row["min_arg"][1])
    with mpmath.workdps(200):
        vec = (a1, a2, _mpf_of_fraction(q))
        vec += tuple(_mpf_of_fraction(v) for v in geom)
        g = gap_mp(vec, _ball_mid_text(beta_text))
        return {"cell": row["cell"], "geometry": list(row["geometry"]),
                "q": row["q"],
                "masses_a1_a2_argmin_midpoint": [mpmath.nstr(a1, 20),
                                                 mpmath.nstr(a2, 20)],
                "certified_min": row["min"],
                "certified_min_lower": row["min_lower"],
                "gap_mp_200dps": mpmath.nstr(g, 30),
                "gap_mp_negative": bool(g < 0),
                "discovery_only": True}


def discovery_agreement(beta_text: str, row: dict[str, Any]) -> dict[str, Any]:
    """DISCOVERY-only mpmath feasibility scan at the worst cell: the
    sampled min over feasible points can never undercut the certified
    lower bound (over-claim tripwire)."""
    geom = tuple(Fraction(v) for v in row["geometry"])
    q = Fraction(row["q"])
    c0, c1, c2 = (Fraction(v) for v in row["mean_line"])
    with mpmath.workdps(60):
        beta = _ball_mid_text(beta_text)
        m_floor = _mpf_of_fraction(FLOOR_HI)
        c0m, c1m, c2m = (_mpf_of_fraction(c0), _mpf_of_fraction(c1),
                         _mpf_of_fraction(c2))
        best = None
        nfeas = 0
        for i in range(SCAN_N):
            a1 = mpmath.mpf(i) / mpmath.mpf(SCAN_N - 1)
            for j in range(SCAN_N - i):
                a2 = mpmath.mpf(j) / mpmath.mpf(SCAN_N - 1)
                mean = c0m + c1m * a1 + c2m * a2
                if mean < m_floor:
                    continue
                nfeas += 1
                vec = (a1, a2, _mpf_of_fraction(q))
                vec += tuple(_mpf_of_fraction(v) for v in geom)
                g = gap_mp(vec, beta)
                if best is None or g < best:
                    best, arg = g, (a1, a2)
        return {"cell": row["cell"], "scan_grid_n": SCAN_N,
                "scan_feasible_points": nfeas,
                "scan_min_gap": None if best is None
                else mpmath.nstr(best, 20),
                "scan_arg": None if best is None
                else [mpmath.nstr(arg[0], 16), mpmath.nstr(arg[1], 16)],
                "certified_min_lower": row["min_lower"],
                "consistent": bool(best is None
                                   or best >= _ball_mid_text(
                                       row["min_lower"]) - TOL_MPMATH),
                "certified_tolerance": str(TOL_ZERO),
                "discovery_only": True}


# ---------------------------------------------------------------------------

def _mean_face_critical(coefd: dict[str, arb],
                        constraints: list[Constraint]
                        ) -> Optional[tuple[Fraction, Fraction]]:
    """The rational midpoint of the mean face's a1-range: a deterministic,
    exactly-feasible evaluation point on the mean chord."""
    a, b, c = constraints[MEAN_INDEX]
    if a == 0 and b == 0:
        return None
    verts = polygon_vertices(constraints)
    face = sorted({v[0] for v in verts if a * v[0] + b * v[1] == c})
    if not face:
        return None
    return ((face[0] + face[-1]) / 2, (c - a * (face[0] + face[-1]) / 2) / b)


def run_mutations(worst: dict[str, Any]) -> list[dict[str, Any]]:
    """Mission mutations on a mean-active worst cell.  Both must FAIL.

    (a) coef_perturbed_1e-6: the witness's same-point discriminator --
    evaluating the PERTURBED quadratic at a fixed exact mean-face point
    must move the value -- plus a KKT re-solve whose certified min must
    differ from the reference cell KKT min.
    (b) drop_mean_chord_constraint: re-solve over the bare triangle (no
    mean-floor chord); its certified minimum must dip below the
    constrained certified min of the same cell.
    """
    out = []
    coefd = {k: arb(v) for k, v in worst["coefficients"].items()}
    constraints = worst["_constraints"]
    crit = _mean_face_critical(coefd, constraints)
    assert crit is not None, "mutation cell has an empty mean face"
    x_ref, y_ref = crit
    v_ref = quad_eval(coefd, _pt(x_ref), _pt(y_ref))

    pert = dict(coefd)
    pert["A22"] = coefd["A22"] + _pt(Fraction(1, 10 ** 6))

    v_mut_ref = quad_eval(pert, _pt(x_ref), _pt(y_ref))
    moved_same_point = not bool((v_mut_ref - v_ref).contains(0))

    cands = kkt_candidates(pert, constraints)
    best_p, _nc, _nf = kkt_min(cands)
    assert best_p is not None, "perturbed QP produced no feasible candidate"
    assert worst["_min_value_arb"] is not None
    moved_ksolve = not bool(
        (best_p["_value"] - worst["_min_value_arb"]).contains(0))
    out.append({
        "mutation": "coef_perturbed_1e-6",
        "cell": worst["cell"],
        "mechanism": ("A22 -> A22 + 1e-6: the value at the same exact "
                      "mean-face point must move, and the KKT re-solve's "
                      "certified min must differ from the reference"),
        "expected": "FAIL",
        "observed": ("FAIL" if (moved_same_point and moved_ksolve)
                     else "PASSED"),
        "reference_point": [str(x_ref), str(y_ref)],
        "value_at_ref": v_ref.str(20),
        "value_mutated_at_ref": v_mut_ref.str(20),
        "value_mutated_ksolve": best_p["value"],
        "value_reference_ksolve": worst["min"],
        "ok": bool(moved_same_point and moved_ksolve),
    })
    assert worst["_min_value_arb"] is not None
    lean = [c for i, c in enumerate(constraints) if i != MEAN_INDEX]
    cands = kkt_candidates(coefd, lean, mean_index=None)
    best, _nc, _nf = kkt_min(cands)
    dipped = best is not None and bool(
        best["_value"].lower() < worst["_min_value_arb"].lower())
    out.append({
        "mutation": "drop_mean_chord_constraint",
        "cell": worst["cell"],
        "mechanism": ("certify over the bare triangle (no mean floor "
                      "chord); the minimum must dip below the same "
                      "cell's constrained certified min"),
        "expected": "FAIL (dips below the constrained certified min)",
        "observed": "FAIL" if dipped else "PASSED",
        "triangle_only_min": None if best is None
        else best["_value"].str(26),
        "triangle_only_arg": None if best is None
        else [best["a1"], best["a2"]],
        "constrained_min_lower": worst["min_lower"],
        "constrained_min_upper": worst["min_upper"],
        "ok": dipped,
    })
    return out


# ---------------------------------------------------------------------------

def _canonical_digest(payload: dict[str, Any]) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def _strip_internals(row: dict[str, Any]) -> dict[str, Any]:
    return {k: v for k, v in row.items() if not k.startswith("_")}


def build_report() -> dict[str, Any]:
    params_mp = solve_equation_parameters(100)
    params_arb = certify_equation_parameters(params_mp)
    beta = params_arb.beta
    floors = certify_floors(params_arb)

    cells: list[dict[str, Any]] = []
    deg_failures = 0
    infeasible = 0
    neg_tol = 0                # certified min lower < -TOL_ZERO
    neg_strict = 0             # certified min ball entirely below 0
    neg_witnesses_rows: list[dict[str, Any]] = []
    near_zero = 0              # -TOL_ZERO <= lower < 0: below tolerance
    nonvac_true = 0
    nonvac_false = 0
    nonvac_undecided = 0
    worst = None
    tol = float(TOL_ZERO)
    geoms = geometry_pairs()
    for gi, geom in enumerate(geoms):
        for qi, q in enumerate(Q_VALUES):
            row = solve_cell(beta, geom, q, FLOOR_HI, params_arb.mean)
            row["cell"] = "g%d|q%d" % (gi, qi)
            cells.append(row)
            nv = row.get("nonvacuous")
            if nv == "true":
                nonvac_true += 1
            elif nv == "false":
                nonvac_false += 1
            else:
                nonvac_undecided += 1
            if row["status"] != "ok":
                infeasible += 1
                continue
            if not row["deg2"]:
                deg_failures += 1
            lower = row["min_lower_f"]
            if worst is None or lower < worst["min_lower_f"]:
                worst = row
            upper = float(arb(row["min_upper"]))
            if lower < -tol:
                neg_tol += 1
                if upper < 0:
                    neg_strict += 1
                if len(neg_witnesses_rows) < WITNESS_CAP:
                    neg_witnesses_rows.append(row)
            elif lower < 0:
                near_zero += 1

    escalations = [escalate_negative(params_arb.beta.str(40), r)
                   for r in neg_witnesses_rows]
    agreement = None
    if worst is not None:
        agreement = discovery_agreement(params_arb.beta.str(40), worst)

    chord_active = [r for r in cells if r["status"] == "ok"
                    and r.get("min_mean_face")]
    if chord_active:
        pool = chord_active
    else:
        pool = [r for r in cells if r["status"] == "ok"]
    mut_row = None
    if worst is not None and worst.get("min_mean_face"):
        mut_row = worst
    elif pool:
        mut_row = min(pool, key=lambda r: r["min_lower_f"])
    mutations = run_mutations(mut_row) if mut_row is not None else []

    all_nonneg = (deg_failures == 0 and neg_tol == 0 and near_zero == 0
                  and nonvac_undecided == 0)
    all_trivial = nonvac_false == len(cells)
    if neg_tol > 0:
        verdict = "NEGATIVE-WITNESS-REPORTED"
    elif all_nonneg and not all_trivial and mutations \
            and all(m["ok"] for m in mutations) \
            and floors["M_UNDER_lt_m_arb"] and floors["M_FLOOR_lt_m_arb"]:
        verdict = "ALL-CELLS-MARGIN-CERTIFIED"
    else:
        verdict = "FAILED"
    claim_status = {
        "ALL-CELLS-MARGIN-CERTIFIED": "MACHINE-VERIFIED",
        "NEGATIVE-WITNESS-REPORTED": "CERTIFIED-NEGATIVE-WITNESS",
        "FAILED": "FAILED",
    }[verdict]

    compact = []
    for r in cells:
        if r["status"] != "ok":
            compact.append({"cell": r["cell"], "geometry": r["geometry"],
                            "q": r["q"], "status": r["status"],
                            "nonvacuous": r["nonvacuous"]})
            continue
        compact.append({
            "cell": r["cell"], "geometry": r["geometry"], "q": r["q"],
            "status": "ok", "deg2": r["deg2"],
            "max_coef_width_f": r["max_coef_width_f"],
            "polygon_vertices": r["polygon_vertices"],
            "n_candidates": r["n_candidates"],
            "n_feasible": r["n_feasible"],
            "min_src": r["min_src"], "min_family": r["min_family"],
            "min_mean_face": r["min_mean_face"],
            "min_arg": r["min_arg"],
            "min": r["min"], "min_lower": r["min_lower"],
            "min_upper": r["min_upper"],
            "min_lower_f": r["min_lower_f"],
            "nonvacuous": r["nonvacuous"],
        })

    worst_detail = None
    if worst is not None:
        worst_detail = {
            "cell": worst["cell"], "geometry": worst["geometry"],
            "q": worst["q"], "mean_line": worst["mean_line"],
            "polygon_vertices": worst["polygon_vertices"],
            "polygon_faces": worst["polygon_faces"],
            "min_src": worst["min_src"], "min_family": worst["min_family"],
            "min_mean_face": worst["min_mean_face"],
            "min_arg": worst["min_arg"], "min": worst["min"],
            "min_lower": worst["min_lower"], "min_upper": worst["min_upper"],
            "certified_scope": ("one (geometry, q) point; all shared "
                                "masses with mixture mean >= m"),
        }

    report = {
        "tool": "liu9_h2_qp_sweep.py",
        "claim": ("on the k/8 ordered grid (%d support values, C(7,3)=%d "
                  "sorted triples per side, %d ordered geometry pairs) x "
                  "the k/8 q grid (%d columns including q=0 and q=1), at "
                  "every cell the raw gap is an exact quadratic in the "
                  "free masses (a1, a2) at the frozen (geometry, q); its "
                  "closed-form KKT minimum over the polygon {triangle} "
                  "AND {mixture mean >= M_FLOOR}, M_FLOOR < m (Arb), is "
                  "Arb-certified >= 0, so the gap is >= 0 on the TRUE "
                  "feasible set {mean >= m} of every cell (superset "
                  "argument)" % (len(SUPPORT_VALUES), len(TRIPLES),
                                 len(geoms), len(Q_VALUES))),
        "claim_status": claim_status,
        "verdict": verdict,
        "claim_scope": (
            "per-cell machine-verified finite statements: ONE cell = one "
            "(x-triple, y-triple, q) grid point, NOT a continuum theorem "
            "in (geometry, q); support and q directions between grid "
            "points stay open"),
        "q_endpoint_reduction": (
            "at q=0 and q=1 the channel term q(1-q)c_ch vanishes exactly "
            "(proved diagonalization gap = F + q(1-q)c_ch, r <= 3); the "
            "q=1 column certifies F(P1) >= 0 on {mean(P1) >= M_FLOOR} "
            "superset {mean(P1) >= m}, the region where the universal "
            "q=1 theorem (524,800-box centered Arb cover, dC_M/dM >= 0) "
            "is published"),
        "feasibility": floors,
        "grid": {
            "support_grid": "k/%d, k=1..%d" % (SUPPORT_GRID,
                                               SUPPORT_GRID - 1),
            "support_values": [str(v) for v in SUPPORT_VALUES],
            "triples_per_side": len(TRIPLES),
            "ordered_geometry_pairs": len(geoms),
            "q_grid": "k/%d, k=0..%d (endpoints included)" % (Q_GRID,
                                                              Q_GRID),
            "q_values": [str(v) for v in Q_VALUES],
            "cells_total": len(cells),
            "enumeration_source": ("liu9_scalar_margin_r3.py "
                                   "geometry_pairs() convention "
                                   "(SUPPORT_GRID = 8)"),
        },
        "extraction": {
            "design_points": [[str(a1), str(a2)] for a1, a2 in DESIGN],
            "validation_points": [[str(a1), str(a2)] for a1, a2
                                  in VALIDATION_POINTS],
            "center": [str(CENTER[0]), str(CENTER[1])],
            "h_step": str(H_STEP),
            "degree_check_threshold": str(QUAD_EPS),
        },
        "coverage_accounting": {
            "cells_total": len(cells),
            "cells_solved": len(cells) - infeasible,
            "cells_infeasible": infeasible,
            "cells_certified_nonnegative": (len(cells) - infeasible
                                            - neg_tol - near_zero),
            "cells_negative_below_tolerance": neg_tol,
            "cells_negative_strict_ball": neg_strict,
            "cells_near_zero_straddle": near_zero,
            "cells_degree_check_failures": deg_failures,
            "nonvacuous_true": nonvac_true,
            "nonvacuous_false": nonvac_false,
            "nonvacuous_undecided": nonvac_undecided,
            "worst_certified_min_lower": (None if worst is None
                                          else worst["min_lower"]),
            "worst_cell_id": None if worst is None else worst["cell"],
            "worst_cell_geometry": None if worst is None
            else worst["geometry"],
            "worst_cell_q": None if worst is None else worst["q"],
            "worst_cell_min_arg": None if worst is None
            else worst["min_arg"],
            "worst_cell_min_src": None if worst is None
            else worst["min_src"],
        },
        "kkt_accounting": {
            "candidate_classes": ["vertex", "interior", "edge",
                                  "mean-chord"],
            "total_candidates": sum(r.get("n_candidates", 0)
                                    for r in cells),
            "total_feasible_candidates": sum(r.get("n_feasible", 0)
                                             for r in cells),
            "min_sources": dict(sorted(
                _count_by(r.get("min_src") for r in cells
                          if r["status"] == "ok").items())),
            "min_families": dict(sorted(
                _count_by(r.get("min_family") for r in cells
                          if r["status"] == "ok").items())),
            "mean_face_active_minima": sum(
                1 for r in cells if r.get("min_mean_face")),
        },
        "negative_escalation": {
            "trigger": ("any certified per-cell min with ball lower "
                        "< -1e-12 is re-verified at 200 dps via gap_mp "
                        "at its argmin midpoint"),
            "negative_cells_below_tolerance": neg_tol,
            "witnesses": escalations,
            "witness_cap": WITNESS_CAP,
        },
        "discovery_agreement": agreement,
        "worst_cell_detail": worst_detail,
        "mutations": mutations,
        "cells": compact,
        "evaluator": ("liu9_objective.evaluate_arb (gap = numerator - "
                      "ehx), ctx.prec = %d; mpmath appears only in the "
                      "labelled DISCOVERY sections via "
                      "liu9_boundary_layer.gap_mp" % CTX_PREC),
    }
    report["report_sha256"] = _canonical_digest(report)
    return report


def _count_by(values) -> dict[str, int]:
    out: dict[str, int] = {}
    for v in values:
        key = "none" if v is None else str(v)
        out[key] = out.get(key, 0) + 1
    return out


def main() -> int:
    report = build_report()
    OUTPUT_DEFAULT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_DEFAULT.write_text(
        json.dumps(report, indent=1, sort_keys=True) + "\n")
    ac = report["coverage_accounting"]
    print("H2 MASS-QP SWEEP (k/8 ordered geometry x k/8 q, closed-form "
          "KKT argmin per cell)")
    print("M_FLOOR < m:", report["feasibility"]["M_FLOOR_lt_m_arb"],
          "| M_UNDER < M_FLOOR:",
          report["feasibility"]["M_UNDER_lt_M_FLOOR_exact"])
    print("cells: %d total, %d solved, %d infeasible"
          % (ac["cells_total"], ac["cells_solved"], ac["cells_infeasible"]))
    print("degree-2 failures: %d" % ac["cells_degree_check_failures"])
    print("certified nonneg: %d | negative(<-1e-12): %d (strict ball %d) "
          "| near-zero straddle: %d"
          % (ac["cells_certified_nonnegative"],
             ac["cells_negative_below_tolerance"],
             ac["cells_negative_strict_ball"],
             ac["cells_near_zero_straddle"]))
    print("worst certified min lower: %s at cell %s q=%s src=%s"
          % (ac["worst_certified_min_lower"], ac["worst_cell_id"],
             ac["worst_cell_q"], ac["worst_cell_min_src"]))
    if report["negative_escalation"]["witnesses"]:
        w = report["negative_escalation"]["witnesses"][0]
        print("NEGATIVE WITNESS cell %s q=%s gap_mp_200dps=%s "
              "(discovery)" % (w["cell"], w["q"], w["gap_mp_200dps"]))
    for m in report["mutations"]:
        print("mutation %-26s %s" % (m["mutation"], m["observed"]))
    if report["discovery_agreement"]:
        print("discovery agreement consistent:",
              report["discovery_agreement"]["consistent"])
    print("verdict:", report["verdict"], "| claim_status:",
          report["claim_status"])
    print("report_sha256", report["report_sha256"])
    return 0 if report["claim_status"] in ("MACHINE-VERIFIED",
                                           "CERTIFIED-NEGATIVE-WITNESS"
                                           ) else 2


if __name__ == "__main__":
    raise SystemExit(main())
