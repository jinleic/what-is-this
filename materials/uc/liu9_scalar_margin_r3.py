#!/usr/bin/env python3
"""R = 3 paired class: certified Arb gap positivity over a full geometry
sweep, with honest interval-width accounting.

CONFIGURATION.  A three-atom paired law is (a1, a2, a3) shared masses on
supports (x1, x2, x3) for the first component and (y1, y2, y3) for the
second, mixed with weight q.  Feasibility for Hypothesis 2 requires the
mixture mean to be at least m = 0.61729091208126497...

EXACT FEASIBILITY.  Instead of solving mean = m (m is irrational) this
module pins the mean to the exact rational

    M_TARGET = 6173/10000  >  m      (verified by an Arb comparison)

by taking the exact rational pinch value

    q0 = (M_TARGET - <a, x>) / (<a, y> - <a, x>)

whenever it lands in [0, 1].  Every retained configuration therefore has
mean exactly M_TARGET >= m as a matter of rational arithmetic, and each
gap is a certified Arb enclosure at an exact rational point.

WHAT IS CERTIFIED.
  1. Point certificates: over EVERY ordered geometry pair on the k/8 grid
     (35 x 35 = 1225 pairs) and every mass point on the k/8 simplex grid,
     the gap enclosure at the exact pinch configuration is strictly
     positive.  Sign certified by Arb, not by floats.
  2. Neighbourhood certificates: around each such point, boxes of half
     width 2^-10, 2^-12 and 2^-14 in (a1, a2, q) are evaluated as
     intervals; the report records, per width, how many boxes have
     nonnegative gap lower bound.  These are positive-measure rigorous
     statements.

WHAT IS NOT CERTIFIED, MEASURED HONESTLY.  A naive interval branch and
bound over the whole (a1, a2, q) unit box does NOT close.  On the root box
the gap enclosure has radius above 10 because interval arithmetic treats
the 36 mass-weight products as independent and so loses the constraint
a1 + a2 + a3 = 1.  The diagnosis section runs that branch and bound on one
geometry with a fixed budget and reports the unresolved count instead of
hiding it.  Closing the full simplex needs the mean-contraction and gauge
machinery of the tube modules, not this file.

This file replaces a retracted version that scanned ONE fixed geometry
with mpmath point evaluations, carried a cell variable (a1v) that was
computed and never used -- so its advertised 768 cells were 48 distinct
evaluations repeated -- and still labelled itself MACHINE-VERIFIED.

Run: math/.venv/bin/python -I -B math/uc/liu9_scalar_margin_r3.py
"""
from __future__ import annotations

import hashlib
import heapq
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

from flint import arb, ctx

from liu9_binding import certify_equation_parameters, solve_equation_parameters
from liu9_objective import evaluate_arb, mean_value

ctx.prec = 320

OUTPUT_DEFAULT = HERE / "verification/results/liu9-scalar-margin-r3.json"
SUPPORT_GRID = 8
MASS_GRID = 8
M_TARGET = Fraction(6173, 10000)
WIDTH_LADDER = (10, 12, 14)
DIAGNOSIS_BUDGET = 1500

Box = tuple[tuple[Fraction, Fraction], ...]


def _pt(value: Fraction) -> arb:
    value = Fraction(value)
    return arb(value.numerator) / arb(value.denominator)


def _ball(lo: Fraction, hi: Fraction) -> arb:
    return _pt(lo) if lo == hi else _pt(lo).union(_pt(hi))


def gap_enclosure(values: tuple[arb, ...], parameters,
                  negate: bool = False) -> arb:
    terms = evaluate_arb(values, parameters.beta)
    return terms.ehx - terms.numerator if negate else terms.numerator - terms.ehx


def geometry_pairs() -> list[tuple[Fraction, ...]]:
    grid = [Fraction(k, SUPPORT_GRID) for k in range(1, SUPPORT_GRID)]
    triples = list(combinations(grid, 3))
    return [tuple(list(x) + list(y)) for x in triples for y in triples]


def mass_points() -> list[tuple[Fraction, Fraction]]:
    out = []
    for i in range(1, MASS_GRID):
        for j in range(1, MASS_GRID - i):
            out.append((Fraction(i, MASS_GRID), Fraction(j, MASS_GRID)))
    return out


def pinch_q(geom: tuple[Fraction, ...], a1: Fraction,
            a2: Fraction) -> Optional[Fraction]:
    x1, x2, x3, y1, y2, y3 = geom
    a3 = 1 - a1 - a2
    xa = a1 * x1 + a2 * x2 + a3 * x3
    ya = a1 * y1 + a2 * y2 + a3 * y3
    if ya == xa:
        return None
    q0 = (M_TARGET - xa) / (ya - xa)
    if q0 < 0 or q0 > 1:
        return None
    return q0


def point_and_box_certificates(parameters, negate: bool = False,
                               widths: tuple[int, ...] = WIDTH_LADDER,
                               m_target: Fraction = M_TARGET) -> dict[str, Any]:
    global M_TARGET
    saved = M_TARGET
    M_TARGET = m_target
    try:
        geoms = geometry_pairs()
        masses = mass_points()
        points = 0
        nonpositive = []
        worst: Optional[arb] = None
        worst_where = None
        width_cleared = {w: 0 for w in widths}
        width_failed_sample: dict[int, list[Any]] = {w: [] for w in widths}
        geometries_with_points = 0
        for geom in geoms:
            local = 0
            for a1, a2 in masses:
                q0 = pinch_q(geom, a1, a2)
                if q0 is None:
                    continue
                points += 1
                local += 1
                vals = (_pt(a1), _pt(a2), _pt(q0)) + tuple(
                    _pt(g) for g in geom)
                gap = gap_enclosure(vals, parameters, negate=negate)
                if worst is None or gap.lower() < worst.lower():
                    worst = gap
                    worst_where = {"geometry": [str(g) for g in geom],
                                   "a1": str(a1), "a2": str(a2),
                                   "q0": str(q0)}
                if not gap > 0 and len(nonpositive) < 8:
                    nonpositive.append({
                        "geometry": [str(g) for g in geom],
                        "a1": str(a1), "a2": str(a2), "q0": str(q0),
                        "gap": gap.str(16)})
                for w in widths:
                    half = Fraction(1, 2 ** w)
                    box = (_ball(a1 - half, a1 + half),
                           _ball(a2 - half, a2 + half),
                           _ball(max(Fraction(0), q0 - half),
                                 min(Fraction(1), q0 + half))) + tuple(
                        _pt(g) for g in geom)
                    bgap = gap_enclosure(box, parameters, negate=negate)
                    if bgap >= 0:
                        width_cleared[w] += 1
                    elif len(width_failed_sample[w]) < 4:
                        width_failed_sample[w].append({
                            "geometry": [str(g) for g in geom],
                            "a1": str(a1), "a2": str(a2), "q0": str(q0),
                            "gap": bgap.str(12)})
            if local:
                geometries_with_points += 1
        return {
            "geometry_pairs_enumerated": len(geoms),
            "geometries_with_feasible_points": geometries_with_points,
            "mass_grid_points": len(masses),
            "feasible_points": points,
            "points_with_nonpositive_gap": len(nonpositive),
            "nonpositive_sample": nonpositive,
            "worst_point_gap": worst.str(20) if worst is not None else None,
            "worst_point": worst_where,
            "box_widths": {
                "2^-%d" % w: {
                    "cleared": width_cleared[w],
                    "of": points,
                    "failed_sample": width_failed_sample[w],
                } for w in widths},
            "all_points_positive": bool(points > 0 and not nonpositive),
        }
    finally:
        M_TARGET = saved


# --- honest diagnosis: naive branch and bound on the full mass box ----------

def _widest(box: Box) -> int:
    return max(range(len(box)), key=lambda i: box[i][1] - box[i][0])


def _bisect(box: Box, index: int) -> tuple[Box, Box]:
    lo, hi = box[index]
    mid = (lo + hi) / 2
    left, right = list(box), list(box)
    left[index] = (lo, mid)
    right[index] = (mid, hi)
    return tuple(left), tuple(right)


def branch_and_bound_diagnosis(parameters, geom: tuple[Fraction, ...],
                               budget: int = DIAGNOSIS_BUDGET) -> dict[str, Any]:
    root: Box = ((Fraction(0), Fraction(1)),) * 3
    heap: list[tuple[float, int, Box]] = [(0.0, 0, root)]
    serial = 0
    processed = cleared = infeasible = unresolved = splits = 0
    root_radius = None
    while heap:
        _, _, box = heapq.heappop(heap)
        if box[0][0] + box[1][0] > 1:
            continue
        processed += 1
        vals = (_ball(*box[0]), _ball(*box[1]), _ball(*box[2])) + tuple(
            _pt(g) for g in geom)
        gap = gap_enclosure(vals, parameters)
        if root_radius is None:
            root_radius = arb(gap.rad()).str(8)
        if mean_value(vals) < parameters.mean:
            infeasible += 1
            continue
        if gap >= 0:
            cleared += 1
            continue
        if processed + len(heap) >= budget:
            unresolved += 1
            continue
        index = _widest(box)
        left, right = _bisect(box, index)
        splits += 1
        for child in (left, right):
            serial += 1
            width = child[_widest(child)][1] - child[_widest(child)][0]
            heapq.heappush(heap, (-float(width), serial, child))
    return {
        "geometry": [str(g) for g in geom],
        "budget": budget,
        "processed": processed,
        "splits": splits,
        "cleared": cleared,
        "mean_infeasible": infeasible,
        "unresolved": unresolved,
        "root_box_gap_radius": root_radius,
        "closes": bool(unresolved == 0),
        "reason": ("interval arithmetic treats the 36 mass-weight products as "
                   "independent, losing a1 + a2 + a3 = 1, so the wide-box gap "
                   "enclosure is orders of magnitude wider than the gap "
                   "itself"),
    }


def run_mutations(parameters, baseline: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    # 1. Feasibility target below m: the guard must reject it.
    bad_target = Fraction(61, 100)
    guard_ok = bool(_pt(bad_target) > parameters.mean)
    out.append({
        "mutation": "mean_target_below_m",
        "mechanism": "M_TARGET = 61/100, which is below m",
        "expected": "FAIL",
        "observed": "FAIL" if not guard_ok else "PASSED",
        "detail": {"arb_comparison_M_target_gt_m": guard_ok}})
    # 2. Gap sign convention flipped: every point must lose positivity.
    flipped = point_and_box_certificates(parameters, negate=True,
                                         widths=(12,))
    out.append({
        "mutation": "gap_sign_flipped",
        "mechanism": "evaluate ehx - numerator instead of numerator - ehx",
        "expected": "FAIL",
        "observed": "FAIL" if not flipped["all_points_positive"] else "PASSED",
        "detail": {"points": flipped["feasible_points"],
                   "nonpositive": flipped["points_with_nonpositive_gap"]}})
    # 3. Box half width 2^-6: the interval cover must stop clearing.
    coarse = point_and_box_certificates(parameters, widths=(6,))
    cleared = coarse["box_widths"]["2^-6"]["cleared"]
    out.append({
        "mutation": "box_half_width_2^-6",
        "mechanism": "neighbourhood boxes 64 times wider than the working "
                     "width",
        "expected": "FAIL",
        "observed": "FAIL" if cleared < coarse["feasible_points"] else "PASSED",
        "detail": {"cleared": cleared, "of": coarse["feasible_points"]}})
    for row in out:
        row["ok"] = bool(row["observed"] == "FAIL")
    return out


def build_report() -> dict[str, Any]:
    parameters = certify_equation_parameters(solve_equation_parameters(100))
    feasibility_guard = bool(_pt(M_TARGET) > parameters.mean)
    certificates = point_and_box_certificates(parameters)
    diagnosis = branch_and_bound_diagnosis(
        parameters, geometry_pairs()[len(geometry_pairs()) // 2])
    muts = run_mutations(parameters, certificates)
    tight = certificates["box_widths"]["2^-%d" % WIDTH_LADDER[-1]]
    ok = (feasibility_guard
          and certificates["all_points_positive"]
          and tight["cleared"] == certificates["feasible_points"]
          and all(r["ok"] for r in muts))
    report = {
        "tool": "liu9_scalar_margin_r3.py",
        "claim": (
            "for every ordered three-atom geometry on the k/%d support grid "
            "and every mass point on the k/%d simplex grid, the exact "
            "rational pinch configuration with mean M_TARGET = %s >= m has "
            "strictly positive gap (Arb-certified), and the gap stays "
            "nonnegative on boxes of half width 2^-%d around every such "
            "point" % (SUPPORT_GRID, MASS_GRID, M_TARGET, WIDTH_LADDER[-1])),
        "claim_status": "MACHINE-VERIFIED finite" if ok else "FAILED",
        "scope_limits": (
            "finite and pointwise-plus-neighbourhood: exact rational supports "
            "on a k/%d grid, exact rational masses on a k/%d grid, mean pinned "
            "to M_TARGET; NOT a theorem about all three-atom paired laws, and "
            "NOT a cover of the (a1, a2, q) simplex -- see the diagnosis "
            "section" % (SUPPORT_GRID, MASS_GRID)),
        "feasibility": {
            "M_TARGET": str(M_TARGET),
            "arb_comparison_M_TARGET_gt_m": feasibility_guard,
            "m_enclosure": parameters.mean.str(24),
            "mechanism": ("q0 chosen as an exact rational so the mixture mean "
                          "equals M_TARGET identically"),
        },
        "certificates": certificates,
        "full_simplex_branch_and_bound_diagnosis": diagnosis,
        "mutations": muts,
        "evaluator": ("liu9_objective.evaluate_arb, gap = numerator - ehx; "
                      "ctx.prec = %d" % ctx.prec),
        "supersedes": ("the retracted point-scan version of this file, whose "
                       "a1v cell variable was never used and whose "
                       "MACHINE-VERIFIED label was overstated"),
    }
    report["report_sha256"] = hashlib.sha256(json.dumps(
        report, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return report


def main() -> int:
    report = build_report()
    OUTPUT_DEFAULT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_DEFAULT.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    c = report["certificates"]
    print("M_TARGET > m:", report["feasibility"]["arb_comparison_M_TARGET_gt_m"])
    print("geometries %d (with feasible points %d), feasible points %d"
          % (c["geometry_pairs_enumerated"],
             c["geometries_with_feasible_points"], c["feasible_points"]))
    print("all point gaps positive:", c["all_points_positive"],
          "| worst point gap:", (c["worst_point_gap"] or "n/a")[:26])
    for key, row in sorted(c["box_widths"].items()):
        print("  box %-8s cleared %d/%d" % (key, row["cleared"], row["of"]))
    d = report["full_simplex_branch_and_bound_diagnosis"]
    print("full-simplex B&B: closes=%s unresolved=%d root gap radius %s"
          % (d["closes"], d["unresolved"], d["root_box_gap_radius"]))
    for row in report["mutations"]:
        print("mutation %-26s %s" % (row["mutation"], row["observed"]))
    print("claim_status:", report["claim_status"])
    print("report_sha256", report["report_sha256"])
    return 0 if report["claim_status"] != "FAILED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
