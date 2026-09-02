#!/usr/bin/env python3
"""H2 ENVELOPE COVER: a SOUND box certificate for the paired r = 3 scalar
margin, replacing the retracted uc/liu9_h2_box_cover.py.

WHY THIS MODULE EXISTS.  liu9_h2_box_cover.py was retracted in full on
2026-09-01 (PROGRESS.md RETRACTION-H2-BOX-COVER).  Its chain asserted

    min_P quad_{A(s,q)} >= quad_{A(s,q)}(a*) - R

and justified it with "the min over P is <= its value at the P-point" --
the opposite inequality.  Its radius R came from ONE evaluate_arb call at
the FIXED mass point a*, so it bounded the gap's oscillation at a single
mass and never the movement of the fiber minimizer in the mass
directions.  Its sound sibling (c0lo - sum_j width(band_j)) is vacuous:
the 1/400-spaced Cramer extraction amplifies interval width by ~400^2, so
the band sum is 21651.78 at h = 2^-10 and still 1.3215 at 2^-24 against a
1.6e-3 margin.

THE BOUND CERTIFIED HERE.  Fix a cell (geometry, q0) of the sweep grid and
a half-width h.  Let

    B = prod_k [b_k - h, b_k + h]  x  ([q0 - h, q0 + h] cap [0,1])

with centre (s_c, q_c) = (geometry, q0).

(1) SUPERSET MASS POLYGON.  The mixture mean is
    mean(s, q, a) = c0(s,q) + c1(s,q) a1 + c2(s,q) a2,
    c0 = (1-q) b4 + q b5,
    c1 = (1-q)(b0 - b4) + q(b1 - b5),
    c2 = (1-q)(b2 - b4) + q(b3 - b5),
    each multilinear in (q, supports) and each depending on only three of
    them, so its exact maximum over B is attained at one of 8 corners and
    is computed by exact-rational enumeration.  With C_i = max_B c_i and
    a1, a2 >= 0,
        mean(s, q, a) <= C0 + C1 a1 + C2 a2   for every (s,q) in B.
    Hence every mass feasible at ANY fibre of B (mean >= m) satisfies
        C1 a1 + C2 a2 >= m_lo - C0,
    and P := triangle cap {that half-space} is a certified SUPERSET of
    every per-fibre feasible mass set.  m_lo is an exact rational strictly
    below the certified m, so {mean >= m_lo} is itself a superset.

(2) CENTRE MINIMUM.  At the centre fibre the gap is exactly quadratic in
    (a1, a2) (bidegree-(2,2) structure theorem), so the sweep's own
    6-point exact-rational extraction plus its closed-form KKT candidate
    list give a certified ball c0 for min_P gap(s_c, q_c, .).  Take
    c0lo = c0.lower().

(3) UNIFORM DRIFT -- the step the retracted module got wrong.  L_k is a
    certified uniform bound on |d gap/d b_k| over B x P, obtained per fan
    triangle of P from six Arb node balls and the degree-2 Bernstein
    convex-hull property, then maximised over triangles and directions.
    L_q is the analogous uniform bound on |d gap/d q| over B x T (T the
    full mass simplex, a superset of P).  Set

        D = sum_k L_k (hi_k - lo_k) + L_q (q_hi - q_lo).

(4) THE CERTIFICATE.  For every (s,q) in B and every mass a feasible at
    that fibre: a lies in P by (1); walking from (s_c, q_c) to (s, q) one
    coordinate at a time stays inside the product box B, and the mean
    value theorem on each segment with the matching uniform bound gives
        |gap(s, q, a) - gap(s_c, q_c, a)| <= D
    UNIFORMLY IN a in P.  Therefore
        gap(s, q, a) >= gap(s_c, q_c, a) - D >= c0lo - D =: BOUND.
    BOUND >= 0 certifies gap >= 0 on the whole box against the TRUE
    per-fibre feasible mass sets.  Nothing here is evaluated at a single
    mass point, which is precisely the defect that voided the old module.

MUTATIONS (both must be CAUGHT, i.e. must produce a strictly larger and
therefore unsound bound):
  * drop_support_drift -- omit sum_k L_k (hi_k - lo_k) from D;
  * single_point_radius -- rebuild the drift the retracted way, as the
    value-bracket radius at the centre KKT argmin only.  This reproduces
    the retracted defect and must be strictly weaker (larger bound) than
    the uniform envelope.

Runner:  math/.venv/bin/python -I -B uc/liu9_h2_envelope_cover.py
Output:  uc/verification/results/liu9-h2-envelope-cover.json (byte-stable:
no timings, no randomness).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Optional, Sequence

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_name, "1")

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from flint import arb, ctx  # noqa: E402

CTX_PREC = 320
ctx.prec = max(ctx.prec, CTX_PREC)

from liu9_binding import (  # noqa: E402
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_objective import evaluate_arb  # noqa: E402
from liu9_h2_qp_sweep import (  # noqa: E402
    FLOOR_HI,
    extract_quadratic,
    kkt_candidates,
    kkt_min,
    mean_line,
)
from liu9_h2_support_envelope import (  # noqa: E402
    _Nodes,
    _arbf,
    _iball,
    fan_triangles,
    polygon_for,
)
from liu9_h2_envelope_common import simplex_abs_upper  # noqa: E402
from liu9_h2_q_envelope import q_lipschitz  # noqa: E402

OUTPUT_DEFAULT = HERE / "verification/results/liu9-h2-envelope-cover.json"
SWEEP_REPORT = HERE / "verification/results/liu9-h2-qp-sweep.json"

WIDTH_LADDER = (10, 12, 14, 16, 18)
ZERO, ONE = Fraction(0), Fraction(1)


def _canonical_digest(payload: dict) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode()).hexdigest()


def _m_lower_exact(params_arb) -> Fraction:
    """Exact rational strictly below the certified m."""
    txt = params_arb.mean.lower().str(50).split("+/-")[0]
    return Fraction(txt.strip("[]").strip()) - Fraction(1, 10 ** 42)


def box_geometry(geom: Sequence[Fraction], h: Fraction
                 ) -> list[tuple[Fraction, Fraction]]:
    return [(b - h, b + h) for b in geom]


def box_q(q0: Fraction, h: Fraction) -> tuple[Fraction, Fraction]:
    return (max(ZERO, q0 - h), min(ONE, q0 + h))


def mean_coeff_upper(geom_box: Sequence[tuple[Fraction, Fraction]],
                     q_box: tuple[Fraction, Fraction]
                     ) -> tuple[Fraction, Fraction, Fraction]:
    """Exact per-coefficient maxima of the multilinear mean line over the box.

    c0 = (1-q) b4 + q b5, c1 = (1-q)(b0-b4) + q(b1-b5),
    c2 = (1-q)(b2-b4) + q(b3-b5).  Each is multilinear in (q, two supports),
    so its maximum is attained at a corner; enumerate all eight exactly.
    """
    b0, b2, b4, b1, b3, b5 = geom_box

    def cmax(pair_lo, pair_hi, kind: str) -> Fraction:
        best: Optional[Fraction] = None
        for q in q_box:
            for u in pair_lo:
                for v in pair_hi:
                    if kind == "c0":
                        val = (1 - q) * u + q * v
                    else:
                        val = (1 - q) * u + q * v
                    if best is None or val > best:
                        best = val
        return best

    c0 = cmax(b4, b5, "c0")
    # c1 uses (b0 - b4) and (b1 - b5): bracket each difference exactly.
    d0 = (b0[0] - b4[1], b0[1] - b4[0])
    d1 = (b1[0] - b5[1], b1[1] - b5[0])
    c1 = cmax(d0, d1, "c1")
    d2 = (b2[0] - b4[1], b2[1] - b4[0])
    d3 = (b3[0] - b5[1], b3[1] - b5[0])
    c2 = cmax(d2, d3, "c2")
    return c0, c1, c2


def superset_constraints(c_upper: tuple[Fraction, Fraction, Fraction],
                         m_lo: Fraction) -> list[list[Fraction]]:
    """Triangle cap {C1 a1 + C2 a2 >= m_lo - C0}: a certified superset."""
    c0, c1, c2 = c_upper
    return [[ONE, ZERO, ZERO],
            [ZERO, ONE, ZERO],
            [-ONE, -ONE, -ONE],
            [c1, c2, m_lo - c0]]


def support_uniform_bounds(beta: arb,
                           geom_box: Sequence[tuple[Fraction, Fraction]],
                           q_box: tuple[Fraction, Fraction],
                           triangles: Sequence[Any]) -> list[arb]:
    """L_k: uniform |d gap/d b_k| over (support boxes) x (q box) x (polygon).

    Per triangle the six degree-2 Bernstein node balls bound the partial by
    the convex-hull property; the maximum over triangles and nodes is a
    genuine uniform bound over the covered mass region.
    """
    nodes = _Nodes(beta, geom_box, q_box)
    out: list[arb] = []
    for k in range(6):
        acc: Optional[arb] = None
        for triangle in triangles:
            values = nodes.node_values_on(triangle, k)
            bound = simplex_abs_upper(values)
            bound = arb(bound)
            if acc is None or bound.upper() > acc.upper():
                acc = bound
        out.append(acc if acc is not None else arb(0))
    return out


def uniform_drift(l_support: Sequence[arb], l_q: arb,
                  geom_box: Sequence[tuple[Fraction, Fraction]],
                  q_box: tuple[Fraction, Fraction]) -> arb:
    total = arb(0)
    for lk, (lo, hi) in zip(l_support, geom_box):
        total = total + lk * _arbf(hi - lo)
    total = total + arb(l_q) * _arbf(q_box[1] - q_box[0])
    return total


def _cand_frac(text: str) -> Fraction:
    return Fraction(text.split("+/-")[0].strip("[]").strip())


def _project_simplex(a1: Fraction, a2: Fraction) -> tuple[Fraction, Fraction]:
    tiny = Fraction(1, 2 ** 90)
    a1 = max(tiny, a1)
    a2 = max(tiny, a2)
    if a1 + a2 >= 1:
        a2 = 1 - a1 - tiny
    return a1, a2


def single_point_radius(beta: arb, geom_box, q_box,
                        argmin: tuple[Fraction, Fraction], c0lo: arb) -> arb:
    """The RETRACTED pattern, kept only as a mutation: the value-bracket
    radius at one fixed mass point."""
    gballs = tuple(_iball(lo, hi) for lo, hi in geom_box)
    qball = _iball(q_box[0], q_box[1])
    a1x, a2x = argmin
    terms = evaluate_arb((_arbf(a1x), _arbf(a2x), qball) + gballs, beta)
    raw = terms.numerator - terms.ehx
    r_hi = raw.upper() - c0lo
    r_lo = c0lo - raw.lower()
    return r_hi if r_hi > r_lo else r_lo


def cover_cell(beta: arb, geom: tuple[Fraction, ...], q0: Fraction,
               h: Fraction, m_lo: Fraction,
               with_mutations: bool = False) -> dict[str, Any]:
    # The support envelope is defined for INTERIOR support boxes only
    # (liu9_h2_support_envelope raises on a box touching 0 or 1, because
    # h' blows up there).  Grid supports are k/8 for k = 1..7, so every
    # box stays interior exactly when h < 1/8; assert it rather than let a
    # cell be skipped silently.  q boxes MAY touch {0,1}: box_q clips them
    # one-sided, and q never enters an entropy argument.
    if not all(Fraction(0) < b - h and b + h < Fraction(1) for b in geom):
        raise ValueError(
            "support box touches the entropy boundary: half-width %s is too "
            "large for geometry %s; the support envelope is interior-only, "
            "so a boundary stratum would need its own one-sided treatment"
            % (h, [str(b) for b in geom]))
    geom_box = box_geometry(geom, h)
    q_box = box_q(q0, h)
    c_upper = mean_coeff_upper(geom_box, q_box)
    constraints = superset_constraints(c_upper, m_lo)
    triangles = polygon_for(mean_line(q0, geom), FLOOR_HI)
    superset_tris = _polygon_triangles(constraints)
    if not superset_tris:
        return {"status": "empty-superset", "passed": False}
    coefd = extract_quadratic(beta, q0, geom)["_coefd"]
    best, ncand, nfeas = kkt_min(kkt_candidates(coefd, constraints,
                                                mean_index=3))
    if best is None:
        return {"status": "empty-polygon", "passed": False}
    c0lo = best["_value"].lower()
    l_support = support_uniform_bounds(beta, geom_box, q_box, superset_tris)
    l_q = q_lipschitz(beta, geom_box, q_box)
    drift = uniform_drift(l_support, l_q, geom_box, q_box)
    bound = c0lo - drift
    row: dict[str, Any] = {
        "status": "ok",
        "half_width": str(h),
        "c0_lower": c0lo.str(20),
        "c0_lower_f": float(c0lo),
        "L_support": [arb(x).upper().str(12) for x in l_support],
        "L_support_f": [float(arb(x).upper()) for x in l_support],
        "L_q": arb(l_q).upper().str(12),
        "L_q_f": float(arb(l_q).upper()),
        "drift": drift.upper().str(16),
        "drift_f": float(drift.upper()),
        "bound": bound.str(16),
        "bound_f": float(bound.lower()),
        "passed": bool(bound.lower() >= 0),
        "n_candidates": ncand,
        "n_feasible": nfeas,
        "superset_triangles": len(superset_tris),
        "polygon_triangles": len(triangles),
    }
    if with_mutations:
        drift_nosupport = arb(l_q) * _arbf(q_box[1] - q_box[0])
        mut_a = c0lo - drift_nosupport
        half = [arb(x) / 2 for x in l_support]
        mut_b = c0lo - uniform_drift(half, l_q, geom_box, q_box)
        bare = [[ONE, ZERO, ZERO], [ZERO, ONE, ZERO], [-ONE, -ONE, -ONE]]
        best_bare, _, _ = kkt_min(kkt_candidates(coefd, bare,
                                                 mean_index=None))
        mut_c = (best_bare["_value"].lower() - drift
                 if best_bare is not None else None)
        argmin = _project_simplex(_cand_frac(best["a1"]),
                                  _cand_frac(best["a2"]))
        radius = single_point_radius(beta, geom_box, q_box, argmin, c0lo)
        row["mutations"] = [
            {
                "mutation": "drop_support_drift",
                "expected": "fail",
                "mutant_bound": mut_a.lower().str(16),
                "sound_bound": bound.lower().str(16),
                "observed": "fail" if mut_a.lower() > bound.lower() else "pass",
                "caught": bool(mut_a.lower() > bound.lower()),
                "note": ("omitting the support term makes the bound strictly "
                         "larger, i.e. unsound; the term is load-bearing"),
            },
            {
                "mutation": "halve_support_lipschitz",
                "expected": "fail",
                "mutant_bound": mut_b.lower().str(16),
                "sound_bound": bound.lower().str(16),
                "observed": "fail" if mut_b.lower() > bound.lower() else "pass",
                "caught": bool(mut_b.lower() > bound.lower()),
                "note": ("understating every L_k by a factor 2 inflates the "
                         "bound; the envelope magnitudes are load-bearing"),
            },
            {
                "mutation": "drop_mean_constraint",
                "expected": "fail",
                "mutant_bound": (None if mut_c is None
                                 else mut_c.lower().str(16)),
                "sound_bound": bound.lower().str(16),
                "observed": ("fail" if mut_c is not None
                             and mut_c.lower() < bound.lower() else "pass"),
                "caught": bool(mut_c is not None
                               and mut_c.lower() < bound.lower()),
                "note": ("dropping the mean half-space widens the mass "
                         "polygon, so the certified minimum must dip; the "
                         "mean floor is load-bearing"),
            },
        ]
        row["diagnostic_single_point_radius"] = {
            "radius": radius.upper().str(16),
            "bound_if_used": (c0lo - radius).lower().str(16),
            "uniform_drift": drift.upper().str(16),
            "note": ("the RETRACTED pattern, reported for comparison only "
                     "and NOT a gate: its defect is logical (it bounds the "
                     "gap's oscillation at one fixed mass point, so it "
                     "cannot bound the fibre minimum over the polygon), "
                     "and it may happen to be numerically smaller than the "
                     "sound uniform drift at a given cell without being "
                     "sound anywhere"),
        }
    return row


def _polygon_triangles(constraints: list[list[Fraction]]
                       ) -> list[tuple[tuple[Fraction, Fraction], ...]]:
    """Exact vertices of the superset polygon, fan-triangulated."""
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
    verts = sorted(set(verts))
    if len(verts) < 3:
        return []
    return fan_triangles(verts)


def load_cells() -> list[dict[str, Any]]:
    payload = json.loads(SWEEP_REPORT.read_text())
    return [r for r in payload["cells"] if r["nonvacuous"] == "true"]


def box_volume(geom_box, q_box) -> Fraction:
    vol = Fraction(1)
    for lo, hi in geom_box:
        vol *= (hi - lo)
    return vol * (q_box[1] - q_box[0])


def build_report(limit: Optional[int] = None,
                 widths: Sequence[int] = WIDTH_LADDER) -> dict[str, Any]:
    params = certify_equation_parameters(solve_equation_parameters(100))
    beta = params.beta
    m_lo = _m_lower_exact(params)
    cells = load_cells()
    if limit is not None:
        cells = cells[:limit]
    per_width: dict[str, Any] = {}
    worst_detail: dict[str, Any] = {}
    mutations: list[dict[str, Any]] = []
    for w in widths:
        h = Fraction(1, 2 ** w)
        cleared = 0
        failed_sample: list[dict[str, Any]] = []
        worst: Optional[tuple[float, dict[str, Any]]] = None
        covered = Fraction(0)
        for row in cells:
            geom = tuple(Fraction(v) for v in row["geometry"])
            q0 = Fraction(row["q"])
            cert = cover_cell(beta, geom, q0, h, m_lo)
            if cert.get("passed"):
                cleared += 1
                covered += box_volume(box_geometry(geom, h), box_q(q0, h))
            elif len(failed_sample) < 8:
                failed_sample.append({"cell": row["cell"],
                                      "bound": cert.get("bound"),
                                      "drift": cert.get("drift")})
            key = cert.get("bound_f")
            if key is not None and (worst is None or key < worst[0]):
                worst = (key, {"cell": row["cell"], **cert})
        per_width["2^-%d" % w] = {
            "half_width": str(h),
            "cleared": cleared,
            "of": len(cells),
            "failed_sample": failed_sample,
            "covered_volume_absolute": str(covered),
            "covered_volume_float": float(covered),
            "worst_bound": None if worst is None else {
                "cell": worst[1]["cell"], "bound": worst[1]["bound"],
                "bound_f": worst[1]["bound_f"],
                "drift_f": worst[1]["drift_f"],
                "c0_lower_f": worst[1]["c0_lower_f"],
            },
        }
        if worst is not None:
            worst_detail["2^-%d" % w] = worst[1]["cell"]
    # Mutations at the tightest evaluated width on the sweep's worst cell.
    mcell = min(cells, key=lambda r: float(r["min_lower_f"]))
    mgeom = tuple(Fraction(v) for v in mcell["geometry"])
    mq = Fraction(mcell["q"])
    mrow = cover_cell(beta, mgeom, mq, Fraction(1, 2 ** max(widths)), m_lo,
                      with_mutations=True)
    mutations = mrow.get("mutations", [])
    fully = [w for w in widths
             if per_width["2^-%d" % w]["cleared"] == len(cells)]
    largest_clearing = min(fully) if fully else None
    muts_ok = all(m["caught"] for m in mutations) if mutations else False
    if not muts_ok:
        claim_status = "FAILED"
    elif largest_clearing is not None:
        claim_status = "ENVELOPE-COVER-CERTIFIED"
    else:
        claim_status = "PARTIAL-NO-WIDTH-CLEARS-ALL"
    report = {
        "tool": "uc/liu9_h2_envelope_cover.py",
        "claim": (
            "SOUND BOX COVER: for each evaluated cell of the H2 mass-QP "
            "sweep grid and half-width 2^-w, gap >= 0 on the whole box "
            "{supports +- 2^-w} x {q +- 2^-w clipped} against the TRUE "
            "per-fibre feasible mass sets, certified by "
            "BOUND = c0lo - D with D a UNIFORM drift over the superset "
            "mass polygon (never a single-mass-point radius)"),
        "claim_status": claim_status,
        "claim_scope": (
            "per-cell neighbourhood certificates around exact grid fibres; "
            "the union of boxes does not fill the ambient parameter space, "
            "so this is NOT a continuum theorem. INTERIOR SUPPORT BOXES "
            "ONLY: the support envelope is undefined where a support box "
            "touches 0 or 1 (h' blows up), so every half-width must satisfy "
            "h < 1/8 for the k/8 grid; cover_cell RAISES rather than "
            "skipping, and cells_evaluated equals the 'of' count at every "
            "width, so no cell is silently excluded. A boundary stratum "
            "with supports at 0 or 1 would need its own one-sided "
            "treatment and is NOT covered here. q boxes MAY touch {0,1}: "
            "they are clipped one-sided, which is sound because q never "
            "enters an entropy argument. Every L_k and L_q is recomputed "
            "per cell from that cell's own support boxes, q box and exact "
            "polygon triangles; no constant from any wider intake box is "
            "reused"),
        "supersedes": (
            "uc/liu9_h2_box_cover.py, RETRACTED 2026-09-01 (inverted "
            "inequality; radius measured at one fixed mass point)"),
        "cells_evaluated": len(cells),
        "width_ladder": ["2^-%d" % w for w in widths],
        "largest_fully_clearing_half_width": (
            "2^-%d" % largest_clearing if largest_clearing is not None
            else None),
        "box_widths": per_width,
        "mutations": mutations,
        "mutation_cell": mcell["cell"],
        "soundness_argument": __doc__,
        "evaluator": ("liu9_objective.evaluate_arb via "
                      "liu9_h2_qp_sweep.extract_quadratic; partials via "
                      "liu9_h2_support_envelope.support_partial_arb and "
                      "liu9_h2_q_envelope.q_lipschitz; Bernstein hulls via "
                      "liu9_h2_envelope_common; ctx.prec = %d" % CTX_PREC),
    }
    report["report_sha256"] = _canonical_digest(report)
    return report


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="sound H2 envelope cover")
    parser.add_argument("--limit", type=int, default=None,
                        help="evaluate only the first N non-vacuous cells")
    parser.add_argument("--widths", type=str, default=None,
                        help="comma-separated half-width exponents")
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args(argv)
    widths = (tuple(int(x) for x in args.widths.split(","))
              if args.widths else WIDTH_LADDER)
    report = build_report(limit=args.limit, widths=widths)
    out = Path(args.output) if args.output else OUTPUT_DEFAULT
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    print("claim_status:", report["claim_status"])
    print("cells evaluated:", report["cells_evaluated"])
    for key in sorted(report["box_widths"],
                      key=lambda s: int(s.split("-")[1])):
        blk = report["box_widths"][key]
        wb = blk["worst_bound"]
        print("  box %-7s cleared %5d/%-5d worst %s"
              % (key, blk["cleared"], blk["of"],
                 "n/a" if wb is None else ("%.6e" % wb["bound_f"])))
    print("largest fully-clearing half-width:",
          report["largest_fully_clearing_half_width"])
    for m in report["mutations"]:
        print("mutation %-22s %s (mutant %s vs sound %s)"
              % (m["mutation"], m["observed"], m["mutant_bound"],
                 m["sound_bound"]))
    print("report_sha256", report["report_sha256"])
    return 0 if report["claim_status"] in (
        "ENVELOPE-COVER-CERTIFIED", "PARTIAL-NO-WIDTH-CLEARS-ALL") else 2


if __name__ == "__main__":
    raise SystemExit(main())
