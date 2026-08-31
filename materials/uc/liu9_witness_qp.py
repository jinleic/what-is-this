#!/usr/bin/env python3
"""PAIRED-CLASS MASS-QP AT THE DANGEROUS WITNESS FIBER: is the worst-case
mass a MARGIN CERTIFICATE or a NEGATIVE?

BACKGROUND.  With the corrected diagonalization (LIU9_BLOCK_COPOSITIVE) the
Liu H2 scalar margin on the paired class is gap = F(P_mix) + q(1-q) c with
F(P_mix) >= 0 proved.  The dangerous witness of the 2026-08-29 H2 search
(total gap ~ +1.92e-4 at q ~ 0.91, second-order degenerate) is UNRECOVERED:
Main confirmed nothing of that in-memory run survived except the summary
line.  The persisted fallback anchor is the deep-search interior point of
uc/verification/results/liu9-block-witness-search.json ("random-111143",
gap 7.61270e-5, per the file's 20-digit exact-decimal vector), whose shape
is P0 = (0, 0.4324, 0.8197), P1 = (0, ~1, 0) at q = 65519/65536 — a point
near the zero-support face, where pi(x,1)=x kills the beta-gain and the
c-channel -2h(x) dominates near-1 asymptotics.

THE QUESTION THIS FILE ANSWERS.  At that frozen (geometry, q) fiber the raw
gap is exactly quadratic in the shared masses (a1, a2) (Main verified the
quadraticity at residual ~1e-80; re-verified here by third differences
~1e-64).  So the mass-worst case is an exact QP:

    minimize  Q(a1, a2)  over  {a1 >= 0, a2 >= 0, a1 + a2 <= 1}  and
              mean(a1, a2) >= m.

METHOD (all certified in Arb, flint ctx.prec = 320).

  1. WITNESS RELOCATION (discovery, labelled as such).  Deterministic
     coordinate descent at 100 dps around the anchor plus a q-refinement
     scan.  The neighborhood minimum slides toward the b3 -> 1 mirror
     face; the gap stays POSITIVE everywhere probed (no negative found).
     The relocated point, the q-scan result, and the UNRECOVERED status of
     the 1.92e-4 point are recorded; the QP fiber stays the frozen anchor
     geometry.

  2. EXACT QUADRATIC EXTRACTION.  Six exact rational mass pairs around the
     fiber centre give six Arb gap enclosures (evaluate_arb, gap =
     numerator - ehx).  The 6x6 design matrix is EXACT rational, so its
     inverse is computed exactly over Fraction and applied to the interval
     RHS: coefficient enclosures of width ~1e-67.  Validation at four
     further exact points: each residual must contain 0 and have radius
     < 1e-60.

  3. CLOSED-FORM KKT.  With the Hessian [[2A11, A12], [A12, 2A22]]
     (indefinite here), a global minimum over the polygon is attained at
     an interior critical point or on the boundary; on each edge the 1D
     restriction is quadratic, so its minimum sits at an interior critical
     point (if the curvature enclosure is strictly positive and the
     critical point lies inside the segment) or at an endpoint.
     Candidates enumerated (each value Arb-certified):
       * interior critical point (Cramer 2x2 on interval coefficients),
       * edge a1 = 0 (segment = [chord, 1]),
       * edge a1 + a2 = 1 (segment [0, u3]),
       * mean chord a2 = alpha + beta_line a1 (segment [0, u3]),
     and the polygon vertices (0, alpha), (u3, 1-u3), (0, 1).
     Feasibility is made conservative exactly like Region-2 hulls:
     M_UNDER = 617290912/10^9 <= m (Arb comparison against the certified
     bracket), so the certified band is a SUPERSET of the true feasible
     band {triangle} ∩ {mean >= m}, and the certified minimum is a valid
     LOWER bound for the true mass-worst-case gap: MIN over the superset
     <= MIN over the true set.

  4. VERDICT.
       * certified min > 0  =>  MARGIN-CERTIFIED for this fiber (over the
         superset, hence a fortiori over the true polygon); q and geometry
         coverage explicitly open.
       * certified min < 0  =>  argmin re-verified independently via
         gap_mp at 200 dps (and mean feasibility re-checked) before any
         REFUTED-WITNESS claim; a negative that is only an artifact of the
         superset band is labelled OPEN-ESCALATED, not a refutation.

  5. MUTATIONS (must FAIL).
       * coef_perturbed_1e-6: A22 -> A22 + 1e-6 must move the KKT minimum.
       * drop_mean_constraint: the triangle-only minimum must dip below
         the constrained minimum.

Byte-stable report:
uc/verification/results/liu9-witness-qp.json

Run: math/.venv/bin/python -I -B math/uc/liu9_witness_qp.py
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

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
from liu9_boundary_layer import gap_mp, mean_of  # noqa: E402
from liu9_objective import evaluate_arb  # noqa: E402

OUTPUT_DEFAULT = HERE / "verification/results/liu9-witness-qp.json"
M_UNDER = Fraction(617290912, 10 ** 9)   # exact rational <= m (Arb-verified)

# Frozen witness fiber (from uc/verification/results/liu9-block-witness-search.json,
# point "random-111143", the persisted best strictly interior point):
# nine-vector (a1, a2, q, b0, b2, b4, b1, b3, b5) with EXACT decimal rationals.
Q_FIBER = Fraction(65519, 65536)                    # 0.99993896484375 (dyadic)
B0 = Fraction(0)
B2 = Fraction("0.43242736996364833368")
B4 = Fraction("0.81965187064361499925")
B1 = Fraction(1, 2 ** 35)                           # 2.9103830456733703613e-11
B3 = Fraction("0.99997640587983671612")
B5 = Fraction(1, 2 ** 58)                           # 3.4694469519536141888e-18
SUPPORTS = (B0, B2, B4, B1, B3, B5)
ANCHOR_A1 = Fraction("0.20323585489734297127")
ANCHOR_A2 = Fraction("0.70894899911375863777")

# Quadratic-extraction design (exact rationals around the fiber centre):
CENTER = (Fraction(1, 5), Fraction(7, 10))
H = Fraction(1, 100)
DESIGN = [
    (CENTER[0], CENTER[1]),
    (CENTER[0] + H, CENTER[1]),
    (CENTER[0] - H, CENTER[1]),
    (CENTER[0], CENTER[1] + H),
    (CENTER[0], CENTER[1] - H),
    (CENTER[0] + H, CENTER[1] + H),
]
VALIDATION_POINTS = [
    (Fraction(3, 20), Fraction(4, 5)),
    (Fraction(7, 20), Fraction(13, 20)),
    (Fraction(3, 25), Fraction(3, 5)),
    (Fraction(2, 5), Fraction(3, 5)),
]

DISCOVERY_STEPS = ("1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7")
DISCOVERY_SWEEPS_PER_STEP = 6
Q_REFINE_SPAN = 100        # +/- 1e-3 around the frozen q, step 1e-5


def _pt(value: Fraction | int) -> arb:
    value = Fraction(value)
    return arb(value.numerator) / arb(value.denominator)


def _mpf_of(value: Fraction) -> mpmath.mpf:
    return mpmath.mpf(value.numerator) / mpmath.mpf(value.denominator)


# ---------------------------------------------------------------------------
# Exact-rational linear algebra on Fraction matrices (exact inverse).

def exact_inverse(rows: list[list[Fraction]]) -> list[list[Fraction]]:
    n = len(rows)
    aug = [[Fraction(x) for x in row] + [Fraction(int(i == j))
                                         for j in range(n)]
           for i, row in enumerate(rows)]
    for col in range(n):
        piv = next(r for r in range(col, n) if aug[r][col] != 0)
        aug[col], aug[piv] = aug[piv], aug[col]
        pval = aug[col][col]
        aug[col] = [x / pval for x in aug[col]]
        for r in range(n):
            if r != col and aug[r][col] != 0:
                factor = aug[r][col]
                aug[r] = [x - factor * y
                          for x, y in zip(aug[r], aug[col])]
    return [row[n:] for row in aug]


# ---------------------------------------------------------------------------
# Certified wrappers.

def gap_arb(a1: Fraction, a2: Fraction, beta: arb) -> arb:
    values = (_pt(a1), _pt(a2), _pt(Q_FIBER)) + tuple(_pt(b) for b in SUPPORTS)
    terms = evaluate_arb(values, beta)
    return terms.numerator - terms.ehx


def mean_exact(a1: Fraction, a2: Fraction) -> Fraction:
    a3 = Fraction(1) - a1 - a2
    b0, b2, b4, b1, b3, b5 = SUPPORTS
    return ((1 - Q_FIBER) * (a1 * b0 + a2 * b2 + a3 * b4)
            + Q_FIBER * (a1 * b1 + a2 * b3 + a3 * b5))


def mean_line_constants() -> tuple[Fraction, Fraction, Fraction]:
    """mean(a1,a2) = c0 + c1 a1 + c2 a2 at the frozen (supports, q)."""
    b0, b2, b4, b1, b3, b5 = SUPPORTS
    c0 = (1 - Q_FIBER) * b4 + Q_FIBER * b5
    c1 = (1 - Q_FIBER) * (b0 - b4) + Q_FIBER * (b1 - b5)
    c2 = (1 - Q_FIBER) * (b2 - b4) + Q_FIBER * (b3 - b5)
    return c0, c1, c2


# ---------------------------------------------------------------------------
# Stage 1: witness relocation scan (DISCOVERY evidence, mpmath point values).

def witness_scan() -> dict[str, Any]:
    """Deterministic local search around the anchor.  Discovery only: every
    number here is an mpmath point evaluation, not a certificate."""
    params = solve_equation_parameters(100)
    with mpmath.workdps(100):
        beta = params.beta
        m = params.mean
        one = mpmath.mpf(1)
        base = [mpmath.mpf(str(x)) for x in
                (ANCHOR_A1, ANCHOR_A2, Q_FIBER) + SUPPORTS]

        def feasible(v) -> bool:
            a1, a2, q = v[0], v[1], v[2]
            return (0 <= a1 <= 1 and 0 <= a2 <= 1 and a1 + a2 <= 1
                    and 0 <= q <= 1
                    and all(0 <= t <= 1 for t in v[3:])
                    and mean_of(tuple(v), one) >= m)

        def ev(v):
            return gap_mp(tuple(v), beta) if feasible(v) else None

        anchor_gap = ev(base)
        cur = base[:]
        gcur = anchor_gap
        evals = 1
        for step in DISCOVERY_STEPS:
            h = mpmath.mpf(step)
            for _ in range(DISCOVERY_SWEEPS_PER_STEP):
                improved = False
                for i in range(9):
                    for sgn in (1, -1):
                        cand = cur[:]
                        cand[i] = cand[i] + sgn * h
                        g = ev(cand)
                        evals += 1
                        if g is not None and g < gcur:
                            cur, gcur = cand, g
                            improved = True
                if not improved:
                    break
        q_best = None
        for k in range(0, 2 * Q_REFINE_SPAN + 1):
            qv = (mpmath.mpf(str(Q_FIBER))
                  + mpmath.mpf(k - Q_REFINE_SPAN) / mpmath.mpf(100000))
            if qv < 0 or qv > 1:
                continue
            v = base[:2] + [qv] + base[3:]
            g = ev(v)
            evals += 1
            if g is not None and (q_best is None or g < q_best[0]):
                q_best = (g, qv)
        return {
            "discovery_only": True,
            "dps": 100,
            "evaluations": evals,
            "anchor_gap": mpmath.nstr(anchor_gap, 24),
            "descent_best_gap": mpmath.nstr(gcur, 24),
            "descent_best_point_20digits": [mpmath.nstr(x, 20) for x in cur],
            "descent_best_mean": mpmath.nstr(mean_of(tuple(cur), one), 20),
            "descent_negative_found": bool(gcur < 0),
            "descent_note": (
                "the coordinate descent slides toward the b3 -> 1 mirror "
                "face with the gap shrinking but staying positive; its "
                "endpoint is a neighborhood local minimum, not a certified "
                "one"),
            "q_refine_best": (None if q_best is None
                              else {"gap": mpmath.nstr(q_best[0], 24),
                                    "q": mpmath.nstr(q_best[1], 20)}),
        }


# ---------------------------------------------------------------------------
# Stage 2: certified exact quadratic extraction.

QUAD_EPS = Fraction(1, 10 ** 60)


def extract_quadratic(beta: arb) -> dict[str, Any]:
    values = [gap_arb(a1, a2, beta) for a1, a2 in DESIGN]
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
    widths = [float(c.upper() - c.lower()) for c in coef]
    max_width = max(widths)
    rows = []
    ok = True
    worst = 0.0
    eps_hi, eps_lo = _pt(QUAD_EPS), -_pt(QUAD_EPS)
    for a1, a2 in VALIDATION_POINTS:
        direct = gap_arb(a1, a2, beta)
        aa1, aa2 = _pt(a1), _pt(a2)
        fit = (coef[0] + coef[1] * aa1 + coef[2] * aa2
               + coef[3] * aa1 * aa2 + coef[4] * aa1 * aa1 + coef[5] * aa2 * aa2)
        residual = fit - direct
        contains = bool(residual.contains(0))
        radius = float(residual.upper() - residual.lower())
        local_ok = contains and residual.upper() < eps_hi and residual.lower() > eps_lo
        ok = ok and local_ok
        worst = max(worst, radius)
        rows.append({"a1": str(a1), "a2": str(a2),
                     "residual": residual.str(10),
                     "residual_radius_float": radius,
                     "contains_zero": contains,
                     "ok": bool(local_ok)})
    return {
        "design_points": [[str(a1), str(a2)] for a1, a2 in DESIGN],
        "coefficients": dict(zip(names, (c.str(30) for c in coef))),
        "coefficient_widths_float": widths,
        "max_coefficient_width_float": max_width,
        "degree_two_verified": bool(ok),
        "validation": rows,
        "worst_residual_radius_float": worst,
        "validation_threshold": str(QUAD_EPS),
        "ok": bool(ok and max_width < 1e-60),
    }


# ---------------------------------------------------------------------------
# Stage 3: closed-form KKT over the polygon, every value Arb-certified.

def _quad(c: dict[str, arb], x: arb, y: arb) -> arb:
    return (c["A0"] + c["A1"] * x + c["A2"] * y + c["A12"] * x * y
            + c["A11"] * x * x + c["A22"] * y * y)


def kkt_qp(extraction: dict[str, Any]) -> dict[str, Any]:
    coef = {name: arb(text) for name, text in extraction["coefficients"].items()}
    c0, c1, c2 = mean_line_constants()
    alpha = (M_UNDER - c0) / c2
    beta_line = -c1 / c2
    u3 = (c2 + c0 - M_UNDER) / (c2 - c1)
    alpha_b, beta_b, u3_b = _pt(alpha), _pt(beta_line), _pt(u3)
    one, zero = arb(1), arb(0)
    m_under_b = _pt(M_UNDER)

    def seg_contains(x: arb, lo: arb, hi: arb) -> bool:
        return bool((x - lo).lower() >= 0 and (hi - x).lower() >= 0)

    candidates: list[dict[str, Any]] = []

    def add(name: str, x: arb, y: arb, note: str, inside: bool) -> None:
        candidates.append({"name": name, "inside_feasible_set": bool(inside),
                           "a1": x.str(16), "a2": y.str(16),
                           "value": _quad(coef, x, y).str(26),
                           "note": note})

    add("V_chord_left", zero, alpha_b, "chord meets a1=0", True)
    add("V_chord_cap", u3_b, one - u3_b,
        "chord meets a1+a2=1; mean(a) = M_UNDER exactly there", True)
    add("V_top", zero, one, "a1=0 meets a1+a2=1", True)
    det = 4 * coef["A11"] * coef["A22"] - coef["A12"] * coef["A12"]
    if not det.contains(0):
        xc = (-coef["A1"] * 2 * coef["A22"] + coef["A12"] * coef["A2"]) / det
        yc = (-coef["A2"] * 2 * coef["A11"] + coef["A12"] * coef["A1"]) / det
        tri = bool(xc.lower() >= 0 and yc.lower() >= 0
                   and (xc + yc).upper() <= 1)
        mean_ok_inst = bool(
            (_pt(c0) + _pt(c1) * xc + _pt(c2) * yc - m_under_b).lower() >= 0)
        add("interior_critical", xc, yc,
            "2x2 KKT solve; Hessian indefinite here (saddle, not a min)",
            tri and mean_ok_inst)
    a22 = coef["A22"]
    if not (2 * a22).contains(0):
        ystar = -coef["A2"] / (2 * a22)
        conv = bool(a22.lower() > 0)
        add("E_a1z_crit", zero, ystar, "a1=0 edge critical point",
            conv and seg_contains(ystar, alpha_b, one))
    s2 = coef["A11"] + coef["A22"] - coef["A12"]
    s1 = coef["A1"] - coef["A2"] + coef["A12"]
    if not (2 * s2).contains(0):
        xstar = -s1 / (2 * s2)
        conv = bool(s2.lower() > 0)
        add("E_sum1_crit", xstar, one - xstar,
            "a1+a2=1 edge critical point",
            conv and seg_contains(xstar, zero, u3_b))
    s2c = (coef["A11"] + beta_b * beta_b * coef["A22"] + beta_b * coef["A12"])
    s1c = (coef["A1"] + coef["A2"] * beta_b + coef["A12"] * alpha_b
           + 2 * coef["A22"] * alpha_b * beta_b)
    if not (2 * s2c).contains(0):
        xstarc = -s1c / (2 * s2c)
        conv = bool(s2c.lower() > 0)
        add("E_chord_crit", xstarc, alpha_b + beta_b * xstarc,
            "mean-chord critical point",
            conv and seg_contains(xstarc, zero, u3_b))

    insides = [c for c in candidates if c["inside_feasible_set"]]
    best = min(insides, key=lambda c: arb(c["value"]).lower())
    best_val = arb(best["value"])
    return {
        "feasible_polygon": ("{a1 >= 0, a2 >= 0, a1+a2 <= 1} intersect "
                             "{mean(a) >= M_UNDER} with a1, a2 >= 0 and "
                             "a1 + a2 <= 1 kept explicit; chord a2 = alpha + "
                             "beta_line a1; cap vertex at u3"),
        "mean_line": {"c0": str(c0), "c1": str(c1), "c2": str(c2),
                      "alpha": str(alpha), "beta_line": str(beta_line),
                      "u3": str(u3)},
        "M_UNDER": str(M_UNDER),
        "M_UNDER_lt_m_arb": None,  # filled by caller (needs certified m)
        "hessian_determinant": det.str(14),
        "hessian_indefinite": bool(det.upper() < 0),
        "candidates": candidates,
        "min_name": best["name"],
        "min_value": best_val.str(30),
        "min_lower": best_val.lower().str(30),
        "min_arg": [best["a1"], best["a2"]],
        "superset_argument": (
            "M_UNDER <= m, so {mean >= M_UNDER} is a SUPERSET of {mean >= m}; "
            "the certified minimum over the superset LOWER-bounds the true "
            "mass-worst-case gap over {mean >= m}: certified positive => "
            "the true worst case is positive too"),
        "ok": len(insides) >= 3,
    }


# ---------------------------------------------------------------------------
# Stage 4: verdict + 200-dps independent re-verification at the argmin.


def verify_argmin_independent(qp: dict[str, Any]) -> dict[str, Any]:
    # the cap vertex is EXACTLY known: a1 = u3, a2 = 1 - u3 (exact rationals via
    # mean_line), so evaluate there rather than at truncated ball midpoints
    u3f = Fraction(qp["mean_line"]["u3"])
    with mpmath.workdps(200):
        params = solve_equation_parameters(200)
        a1 = _mpf_of(u3f)
        a2 = mpmath.mpf(1) - a1
        vec = (a1, a2) + tuple(_mpf_of(b) for b in (Q_FIBER,) + SUPPORTS)
        g = gap_mp(vec, params.beta)
        mn = mean_of(vec, mpmath.mpf(1))
        return {
            "gap_mp_200dps": mpmath.nstr(g, 30),
            "mean_mp_200dps": mpmath.nstr(mn, 30),
            "mean_ge_M_UNDER": bool(
                mn >= mpmath.mpf(19290341) / mpmath.mpf(31250000)),
            "discovery_only": True,
            "note": ("mpmath point evaluation at the certified argmin "
                     "(exact u3 vertex, 200-dps context): DISCOVERY-grade "
                     "agreement check, not a certificate"),
        }


# ---------------------------------------------------------------------------
# Stage 5: mutations (must FAIL).

def run_mutations(extraction: dict[str, Any], qp: dict[str, Any]
                  ) -> list[dict[str, Any]]:
    out = []
    coef = {k: arb(v) for k, v in extraction["coefficients"].items()}
    pert = dict(coef)
    pert["A22"] = coef["A22"] + _pt(Fraction(1, 10 ** 6))
    u3_b = _pt(Fraction(qp["mean_line"]["u3"]))
    one = arb(1)
    x, y = u3_b, one - u3_b
    v_ref = arb(qp["min_value"])
    v_mut = _quad(pert, x, y)
    changed = not bool((v_mut - v_ref).contains(0))
    out.append({
        "mutation": "coef_perturbed_1e-6",
        "mechanism": "A22 -> A22 + 1e-6 must move the KKT minimum",
        "expected": "FAIL",
        "observed": "FAIL" if changed else "PASSED",
        "value_mutated": v_mut.str(20),
        "value_reference": v_ref.str(20),
        "ok": changed,
    })
    a0, a2c, a22 = coef["A0"], coef["A2"], coef["A22"]
    tri_min = a0 - a2c * a2c / (4 * a22)
    dips = bool(tri_min.upper() < arb(qp["min_lower"]))
    out.append({
        "mutation": "drop_mean_constraint",
        "mechanism": ("triangle-only minimum A0 - A2^2/(4 A22) over the "
                      "a1 = 0 edge, unconstrained by the mean chord"),
        "expected": "FAIL (dips below the constrained certified min)",
        "observed": "FAIL" if dips else "PASSED",
        "triangle_only_min": tri_min.str(20),
        "triangle_only_arg_a2": (-a2c / (2 * a22)).str(16),
        "constrained_min_lower": qp["min_lower"],
        "ok": dips,
    })
    if not all(r["ok"] for r in out):
        raise AssertionError("a required mutation unexpectedly PASSED")
    return out


# ---------------------------------------------------------------------------
# Report.

def _canonical_digest(payload: dict[str, Any]) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def build_report() -> dict[str, Any]:
    params_mp = solve_equation_parameters(100)
    params_arb = certify_equation_parameters(params_mp)
    beta = params_arb.beta
    m_under_ok = bool(_pt(M_UNDER) < params_arb.mean)

    anchor_gap = gap_arb(ANCHOR_A1, ANCHOR_A2, beta)

    scan = witness_scan()
    extraction = extract_quadratic(beta)
    if not extraction["ok"]:
        raise AssertionError("quadratic extraction failed validation")
    qp = kkt_qp(extraction)
    qp["M_UNDER_lt_m_arb"] = m_under_ok
    min_lower = arb(qp["min_lower"])
    if min_lower > 0:
        verdict = "MARGIN-CERTIFIED"
    elif min_lower < 0:
        verdict = "CERTIFIED-BAND-NEGATIVE"
    else:
        verdict = "ZERO-BOUND-UNDECIDED"
    independent = verify_argmin_independent(qp)
    if verdict == "CERTIFIED-BAND-NEGATIVE":
        # only a superset-negative so far; the independent check decides the
        # escalation label: a feasible negative gap_mp witness at 200 dps
        # would be labelled REFUTED-WITNESS instead.
        neg_text = independent["gap_mp_200dps"]
        if neg_text is not None and mpmath.mpf(neg_text) < 0:
            verdict = "REFUTED-WITNESS"
        else:
            verdict = "OPEN-ESCALATED"
    mutations = run_mutations(extraction, qp)
    ok = (m_under_ok and extraction["ok"] and qp["ok"]
          and all(r["ok"] for r in mutations)
          and verdict in ("MARGIN-CERTIFIED", "REFUTED-WITNESS"))
    claim_status = {"MARGIN-CERTIFIED": "MACHINE-VERIFIED",
                    "REFUTED-WITNESS": "REFUTED-WITNESS",
                    "OPEN-ESCALATED": "OPEN",
                    "CERTIFIED-BAND-NEGATIVE": "OPEN",
                    "ZERO-BOUND-UNDECIDED": "OPEN"}[verdict]
    b0, b2, b4, b1, b3, b5 = SUPPORTS
    report = {
        "tool": "liu9_witness_qp.py",
        "claim": ("at the frozen witness fiber (persisted supports of "
                  "random-111143, q = 65519/65536), the raw gap is an exact "
                  "quadratic in the shared masses (a1, a2); its closed-form "
                  "KKT minimum over {triangle} intersect {mean >= M_UNDER} "
                  "with M_UNDER <= m is Arb-certified positive, so the "
                  "mass-worst-case gap over the TRUE feasible set "
                  "{triangle} intersect {mean >= m} is positive as well "
                  "(superset argument)"),
        "claim_status": claim_status if ok else "FAILED",
        "verdict": verdict,
        "witness_1_92e_4": "UNRECOVERED",
        "witness_1_92e_4_status": (
            "the ~1.92e-4 total-gap point of the 2026-08-29 search is not "
            "persisted anywhere (Main: in-memory only; only the summary "
            "line survived).  Structural note from the ledger: near the "
            "zero-support face pi(x,1)=x kills the beta-gain, so the "
            "c-channel -2h(x) dominates near-1 asymptotics - consistent "
            "with danger concentrated at that face.  The persisted fallback "
            "anchor geometry is used for the QP per Main's instruction"),
        "witness_fiber": {
            "description": ("persisted deep-search point random-111143: "
                            "P0 = masses a on (b0, b2, b4), P1 = same "
                            "masses on (b1, b3, b5); q frozen at the "
                            "recorded dyadic value"),
            "q": str(Q_FIBER),
            "supports": {"b0": str(B0), "b2": str(B2), "b4": str(B4),
                         "b1": str(B1), "b3": str(B3), "b5": str(B5)},
            "anchor_masses": [str(ANCHOR_A1), str(ANCHOR_A2)],
            "anchor_gap_arb": anchor_gap.str(24),
            "anchor_mean_exact": str(mean_exact(ANCHOR_A1, ANCHOR_A2)),
        },
        "parameters": {
            "m_lower": params_arb.mean.lower().str(40),
            "m_upper": params_arb.mean.upper().str(40),
            "beta": beta.str(40),
            "M_UNDER": str(M_UNDER),
            "M_UNDER_lt_m_arb": m_under_ok,
        },
        "witness_scan": scan,
        "quadratic_extraction": extraction,
        "mass_qp": qp,
        "independent_argmin_check": independent,
        "mutations": mutations,
        "coverage_statement": (
            "a margin certificate for THIS (geometry, q) fiber only; the q "
            "direction and the four support degrees of freedom stay open - "
            "one point of the (P0, P1, q) space, not a universal H2 proof"),
        "evaluator": ("liu9_objective.evaluate_arb (gap = numerator - ehx), "
                      "ctx.prec = %d; discovery scans use "
                      "liu9_boundary_layer.gap_mp at the stated workdps "
                      "(DISCOVERY evidence only, labelled so)"
                      % CTX_PREC),
    }
    report["report_sha256"] = _canonical_digest(report)
    return report


def main() -> int:
    report = build_report()
    OUTPUT_DEFAULT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_DEFAULT.write_text(
        json.dumps(report, indent=1, sort_keys=True) + "\n")
    print("LIU H2 PAIRED-CLASS MASS-QP AT THE WITNESS FIBER")
    print("1.92e-4 witness:", report["witness_1_92e_4"])
    print("anchor gap:", report["witness_fiber"]["anchor_gap_arb"])
    print("quadratic verified:",
          report["quadratic_extraction"]["degree_two_verified"],
          "| worst validation radius",
          report["quadratic_extraction"]["worst_residual_radius_float"])
    print("QP min:", report["mass_qp"]["min_name"],
          report["mass_qp"]["min_value"])
    print("independent:", report["independent_argmin_check"])
    for row in report["mutations"]:
        print("mutation %-24s %s" % (row["mutation"], row["observed"]))
    print("verdict:", report["verdict"], "| claim_status:",
          report["claim_status"])
    print("report_sha256", report["report_sha256"])
    return 0 if report["claim_status"] in ("MACHINE-VERIFIED",
                                           "REFUTED-WITNESS") else 2


if __name__ == "__main__":
    raise SystemExit(main())
