#!/usr/bin/env python3
"""REGION 2, CURVE I (the q = 1 edge): certified Arb minimum over the
feasible mass band.

RETRACTION FIRST.  The previous version of this file was wrong, in three
separate ways, and its PROVED label was unsound:

  * The three-point coefficient extraction was algebraically incorrect.
    For a quadratic g(a) = A2 a^2 + A1 a + A0 sampled at a = 0, 1/2, 1 the
    second difference is A2/2, not A2, so the correct inversion is

        A0 = v0,  A2 = 2 (v1 - 2 v_half + v0),  A1 = (v1 - v0) - A2,

    while the old file used A2 = v1 - 2 v_half + v0 and
    A1 = v_half - A0 - A2/4, i.e. half of both true coefficients.  Every
    edge and vertex value it "certified" was therefore computed from the
    wrong polynomial.  The mutation `legacy_three_point_formula` below
    encodes this bug as a permanent regression test.
  * `curve2_pinch_scan` was dead code: its inner expression was guarded by
    `if False`, its accumulator was always None, and it returned a
    hardcoded placeholder dict.  Deleted; curve II lives in
    `uc/liu9_scalar_margin_region2_curve2.py`, which is real.
  * The advertised 4 x 4 x 4 x 4 geometry grid was 10 distinct problems
    repeated 10 times each: at q = 1 the first component carries weight
    zero, so the gap does not depend on the x supports at all.  That
    independence is now VERIFIED (not assumed) and only y pairs are
    enumerated.

WHAT THIS FILE NOW CERTIFIES.  At q = 1 the raw gap restricted to a
two-atom paired law is a polynomial of degree 2 in the shared mass a1,

    gapA(a1) = A2 a1^2 + A1 a1 + A0,

with coefficients recovered from three exact rational samples and then
VERIFIED at further exact points (the interpolation residual must contain
zero as an Arb interval).  Feasibility `mean = a1 y1 + (1-a1) y2 >= m` is
enforced conservatively through an exact rational `M_UNDER <= m`, so the
band certified is a superset of the true feasible band.  The minimum over
that band is bounded below by interval subdivision of the certified
quadratic -- no vertex or convexity case analysis, so no sign-of-A2
fragility.

Byte-stable report:
uc/verification/results/liu9-scalar-margin-region2-hulls.json
Run: math/.venv/bin/python -I -B math/uc/liu9_scalar_margin_region2hulls.py
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

from flint import arb, ctx

from liu9_binding import certify_equation_parameters, solve_equation_parameters
from liu9_objective import evaluate_arb

ctx.prec = 320

OUTPUT_DEFAULT = (
    HERE / "verification/results/liu9-scalar-margin-region2-hulls.json")
Y_GRID = 16
M_UNDER = Fraction(617290912, 10 ** 9)   # exact rational <= m, Arb-verified
BAND_CELLS = 64
X_REFERENCE = (Fraction(1, 8), Fraction(1, 4))
X_ALTERNATE = (Fraction(3, 8), Fraction(5, 8))


def _pt(value: Fraction) -> arb:
    value = Fraction(value)
    return arb(value.numerator) / arb(value.denominator)


def _ball(lo: Fraction, hi: Fraction) -> arb:
    return _pt(lo) if lo == hi else _pt(lo).union(_pt(hi))


def gap_q1(a1: Fraction, ypair: tuple[Fraction, Fraction], parameters,
           xpair: tuple[Fraction, Fraction] = X_REFERENCE) -> arb:
    """Arb enclosure of the raw gap at q = 1 for masses (a1, 1-a1)."""
    y1, y2 = ypair
    x1, x2 = xpair
    values = (_pt(a1), _pt(1 - a1), _pt(Fraction(1)),
              _pt(x1), _pt(x2), _pt(Fraction(1, 2)),
              _pt(y1), _pt(y2), _pt(Fraction(1, 2)))
    terms = evaluate_arb(values, parameters.beta)
    return terms.numerator - terms.ehx


def extract_quadratic(ypair, parameters, legacy: bool = False):
    """Certified (A0, A1, A2) from exact samples at a1 = 0, 1/2, 1."""
    v0 = gap_q1(Fraction(0), ypair, parameters)
    vh = gap_q1(Fraction(1, 2), ypair, parameters)
    v1 = gap_q1(Fraction(1), ypair, parameters)
    a0 = v0
    if legacy:                      # the retracted, incorrect inversion
        a2 = v1 - 2 * vh + v0
        a1c = vh - a0 - a2 / 4
    else:
        a2 = 2 * (v1 - 2 * vh + v0)
        a1c = (v1 - v0) - a2
    return a0, a1c, a2


def degree_check(ypair, parameters, legacy: bool = False) -> dict[str, Any]:
    """The recovered quadratic must reproduce further exact samples."""
    a0, a1c, a2 = extract_quadratic(ypair, parameters, legacy=legacy)
    rows = []
    ok = True
    for probe in (Fraction(1, 4), Fraction(3, 4), Fraction(7, 8),
                  Fraction(1, 16)):
        pa = _pt(probe)
        predicted = a2 * pa * pa + a1c * pa + a0
        direct = gap_q1(probe, ypair, parameters)
        residual = predicted - direct
        contains = bool(residual.contains(0))
        ok = ok and contains
        rows.append({"a1": str(probe), "residual": residual.str(8),
                     "contains_zero": contains})
    return {"degree_two_verified": bool(ok), "probes": rows}


def feasible_band(ypair: tuple[Fraction, Fraction]) -> Optional[
        tuple[Fraction, Fraction]]:
    """Conservative superset of {a1 in [0,1] : a1 y1 + (1-a1) y2 >= m}."""
    y1, y2 = ypair
    if max(y1, y2) < M_UNDER:
        return None
    if y1 == y2:
        return (Fraction(0), Fraction(1)) if y1 >= M_UNDER else None
    root = (M_UNDER - y2) / (y1 - y2)
    if y1 > y2:
        lo = max(Fraction(0), root)
        hi = Fraction(1)
    else:
        lo = Fraction(0)
        hi = min(Fraction(1), root)
    return (lo, hi) if lo <= hi else None


def band_minimum_lower(a0: arb, a1c: arb, a2: arb,
                       band: tuple[Fraction, Fraction],
                       cells: int = BAND_CELLS) -> tuple[arb, str]:
    """Interval lower bound of the quadratic over the band by subdivision."""
    lo, hi = band
    worst: Optional[arb] = None
    worst_cell = ""
    for k in range(cells):
        cl = lo + (hi - lo) * Fraction(k, cells)
        ch = lo + (hi - lo) * Fraction(k + 1, cells)
        a = _ball(cl, ch)
        value = a2 * a * a + a1c * a + a0
        if worst is None or value.lower() < worst.lower():
            worst = arb(value.lower())
            worst_cell = "[%s, %s]" % (cl, ch)
    return worst, worst_cell


def certify_curve1(parameters, legacy: bool = False,
                   ignore_feasibility: bool = False,
                   m_under: Fraction = M_UNDER) -> dict[str, Any]:
    global M_UNDER
    saved = M_UNDER
    M_UNDER = m_under
    try:
        grid = [Fraction(k, Y_GRID) for k in range(1, Y_GRID)]
        pairs = list(combinations(grid, 2))
        rows = []
        worst: Optional[arb] = None
        worst_pair = None
        negative = 0
        degree_failures = 0
        skipped = 0
        for ypair in pairs:
            band = ((Fraction(0), Fraction(1)) if ignore_feasibility
                    else feasible_band(ypair))
            if band is None:
                skipped += 1
                continue
            check = degree_check(ypair, parameters, legacy=legacy)
            if not check["degree_two_verified"]:
                degree_failures += 1
            a0, a1c, a2 = extract_quadratic(ypair, parameters, legacy=legacy)
            lower, cell = band_minimum_lower(a0, a1c, a2, band)
            if not lower >= 0:
                negative += 1
            if worst is None or lower.lower() < worst.lower():
                worst = lower
                worst_pair = {"y1": str(ypair[0]), "y2": str(ypair[1]),
                              "band": [str(band[0]), str(band[1])],
                              "cell": cell}
            rows.append({
                "y1": str(ypair[0]), "y2": str(ypair[1]),
                "band": [str(band[0]), str(band[1])],
                "A0": a0.str(12), "A1": a1c.str(12), "A2": a2.str(12),
                "band_min_lower": lower.str(16),
                "nonnegative": bool(lower >= 0),
                "degree_two_verified": check["degree_two_verified"],
            })
        return {
            "y_pairs_enumerated": len(pairs),
            "y_pairs_certified": len(rows),
            "y_pairs_infeasible": skipped,
            "negative_bands": negative,
            "degree_check_failures": degree_failures,
            "worst_band_min_lower": worst.str(20) if worst is not None else None,
            "worst_pair": worst_pair,
            "all_nonnegative": bool(rows and negative == 0),
            "rows": rows,
        }
    finally:
        M_UNDER = saved


def x_independence_check(parameters) -> dict[str, Any]:
    """At q = 1 the x supports carry weight zero; verify, do not assume."""
    rows = []
    ok = True
    for ypair in ((Fraction(1, 2), Fraction(7, 8)),
                  (Fraction(3, 16), Fraction(11, 16)),
                  (Fraction(5, 8), Fraction(15, 16))):
        for a1 in (Fraction(3, 8), Fraction(7, 8)):
            ref = gap_q1(a1, ypair, parameters, xpair=X_REFERENCE)
            alt = gap_q1(a1, ypair, parameters, xpair=X_ALTERNATE)
            residual = ref - alt
            contains = bool(residual.contains(0))
            ok = ok and contains
            rows.append({"y": [str(ypair[0]), str(ypair[1])], "a1": str(a1),
                         "residual": residual.str(8),
                         "contains_zero": contains})
    return {"x_independence_verified": bool(ok), "probes": rows}


def run_mutations(parameters) -> list[dict[str, Any]]:
    out = []
    legacy = certify_curve1(parameters, legacy=True)
    out.append({
        "mutation": "legacy_three_point_formula",
        "mechanism": ("the retracted inversion A2 = v1 - 2 v_half + v0, "
                      "A1 = v_half - A0 - A2/4"),
        "expected": "FAIL",
        "observed": ("FAIL" if legacy["degree_check_failures"] > 0
                     else "PASSED"),
        "detail": {"degree_check_failures": legacy["degree_check_failures"],
                   "of": legacy["y_pairs_certified"]}})
    wide = certify_curve1(parameters, ignore_feasibility=True)
    out.append({
        "mutation": "band_widened_to_unit_interval",
        "mechanism": "drop the mean >= m band and certify over all of [0,1]",
        "expected": "FAIL",
        "observed": "FAIL" if not wide["all_nonnegative"] else "PASSED",
        "detail": {"negative_bands": wide["negative_bands"],
                   "of": wide["y_pairs_certified"]}})
    high = certify_curve1(parameters, m_under=Fraction(99, 100))
    out.append({
        "mutation": "mean_floor_raised_to_99_percent",
        "mechanism": "M_UNDER = 99/100 leaves almost no feasible y pair",
        "expected": "FAIL",
        "observed": ("FAIL" if high["y_pairs_certified"]
                     < certify_curve1(parameters)["y_pairs_certified"]
                     else "PASSED"),
        "detail": {"certified_pairs": high["y_pairs_certified"]}})
    for row in out:
        row["ok"] = bool(row["observed"] == "FAIL")
    return out


def build_report() -> dict[str, Any]:
    parameters = certify_equation_parameters(solve_equation_parameters(100))
    m_under_ok = bool(_pt(M_UNDER) < parameters.mean)
    independence = x_independence_check(parameters)
    curve1 = certify_curve1(parameters)
    muts = run_mutations(parameters)
    ok = (m_under_ok
          and independence["x_independence_verified"]
          and curve1["all_nonnegative"]
          and curve1["degree_check_failures"] == 0
          and all(r["ok"] for r in muts))
    report = {
        "tool": "liu9_scalar_margin_region2hulls.py",
        "claim": ("at q = 1, for every ordered y pair on the k/%d grid whose "
                  "feasible mass band is nonempty, the certified quadratic "
                  "gapA has nonnegative Arb lower bound over a conservative "
                  "superset of that band" % Y_GRID),
        "claim_status": "PROVED" if ok else "FAILED",
        "retraction": (
            "supersedes the previous version of this file, whose three-point "
            "coefficient inversion was algebraically wrong (A2 and A1 both "
            "half their true values), which contained a dead placeholder "
            "curve-II scan, and whose geometry count was inflated tenfold by "
            "x supports that do not enter at q = 1"),
        "feasibility": {
            "M_UNDER": str(M_UNDER),
            "arb_comparison_M_UNDER_le_m": m_under_ok,
            "m_enclosure": parameters.mean.str(24),
            "direction": ("M_UNDER <= m, so the certified band is a superset "
                          "of the true feasible band"),
        },
        "x_independence": independence,
        "curve1": curve1,
        "band_minimisation": ("interval subdivision into %d cells; no "
                              "convexity or vertex case analysis"
                              % BAND_CELLS),
        "mutations": muts,
        "evaluator": ("liu9_objective.evaluate_arb, gap = numerator - ehx; "
                      "ctx.prec = %d" % ctx.prec),
    }
    report["report_sha256"] = hashlib.sha256(json.dumps(
        report, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return report


def main() -> int:
    report = build_report()
    OUTPUT_DEFAULT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_DEFAULT.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    print("M_UNDER <= m:", report["feasibility"]["arb_comparison_M_UNDER_le_m"])
    print("x-independence verified:",
          report["x_independence"]["x_independence_verified"])
    c = report["curve1"]
    print("y pairs: %d enumerated, %d certified, %d infeasible"
          % (c["y_pairs_enumerated"], c["y_pairs_certified"],
             c["y_pairs_infeasible"]))
    print("degree-2 failures: %d | negative bands: %d | worst lower: %s"
          % (c["degree_check_failures"], c["negative_bands"],
             (c["worst_band_min_lower"] or "n/a")[:26]))
    for row in report["mutations"]:
        print("mutation %-34s %s" % (row["mutation"], row["observed"]))
    print("claim_status:", report["claim_status"])
    print("report_sha256", report["report_sha256"])
    return 0 if report["claim_status"] == "PROVED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
