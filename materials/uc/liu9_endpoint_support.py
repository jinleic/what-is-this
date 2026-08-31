#!/usr/bin/env python3
"""Certified endpoint-aware support theorem for Liu H2's q=1 blocker."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path
from typing import Any, Optional, Sequence

for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

import mpmath
from flint import ctx

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from liu9_binding import certify_equation_parameters, solve_equation_parameters
from liu9_boundary_layer import gap_mp
from liu9_qone_face import (
    _single_zero_margin, active_compact, best_qone_endpoint_lower,
    canonical_qone_box, contract_qone_mean,
)

ctx.prec = max(ctx.prec, 320)

OUTPUT_DEFAULT = HERE / "verification/results/liu9-endpoint-support.json"


def _canonical_digest(payload: dict[str, Any]) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode()).hexdigest()


def certify() -> dict[str, Any]:
    mp_parameters = solve_equation_parameters(100)
    parameters = certify_equation_parameters(mp_parameters)
    blocker = canonical_qone_box((
        (13 / 32, 7 / 16),
        (0.0, 1 / 32),
        (1 / 32, 1 / 16),
        (13 / 32, 7 / 16),
        (15 / 16, 1.0),
    ))
    uncontracted, uncontracted_method = best_qone_endpoint_lower(
        blocker, parameters.mean, parameters.beta)
    contracted, contraction_count = contract_qone_mean(
        blocker, parameters.mean)
    if contracted is None:
        raise AssertionError("the exact blocker unexpectedly became infeasible")
    lower, method = best_qone_endpoint_lower(
        contracted, parameters.mean, parameters.beta)
    if not lower > 0:
        raise AssertionError("endpoint support bound did not close the blocker")
    if not uncontracted < 0:
        raise AssertionError("mean-contractor mutation no longer discriminates")

    with mpmath.workdps(100):
        a1 = mpmath.mpf(13) / 32
        a2 = mpmath.mpf(0)
        b1 = mpmath.mpf(1) / 16
        b3 = mpmath.mpf(13) / 32
        a3 = 1 - a1
        b5 = (mp_parameters.mean - a1 * b1) / a3
        point = (
            a1, a2, mpmath.mpf(1),
            mpmath.mpf("0.5"), mpmath.mpf("0.5"), mpmath.mpf("0.5"),
            b1, b3, b5,
        )
        point_gap = gap_mp(point, mp_parameters.beta)
        point_mean = a1 * b1 + a3 * b5
    if point_mean < mp_parameters.mean or point_gap <= 0:
        raise AssertionError("sanctioned raw-gap endpoint control failed")

    sharp_zero = _single_zero_margin(parameters.mean, parameters.beta)
    report: dict[str, Any] = {
        "tool": "liu9_endpoint_support.py",
        "claim_status": "PROVED",
        "scope": (
            "exact q=1 feasible subset of the blocker in "
            "(a1,a2,b1,b3,b5); raw gap only"
        ),
        "blocker": active_compact(blocker),
        "contracted_blocker": active_compact(contracted),
        "mean_contractions": contraction_count,
        "endpoint_lower": str(lower),
        "endpoint_lower_float": float(lower.lower()),
        "endpoint_method": method,
        "uncontracted_mutation_lower": str(uncontracted),
        "uncontracted_mutation_method": uncontracted_method,
        "sharp_single_zero_margin_y0_1_16": str(sharp_zero),
        "sanctioned_gap_mp_control": {
            "values": [mpmath.nstr(value, 80) for value in point],
            "mean": mpmath.nstr(point_mean, 80),
            "raw_gap": mpmath.nstr(point_gap, 80),
        },
        "theorem": (
            "mean deficit contracts the blocker to b1>0.0579 and b5>0.9968; "
            "the factored s*log(s) endpoint extension plus a nonpositive mean "
            "multiplier gives raw gap >= endpoint_lower > 0"
        ),
        "quotient_used": False,
    }
    report["report_sha256"] = _canonical_digest(report)
    return report


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT_DEFAULT)
    args = parser.parse_args(argv)
    report = certify()
    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    print("LIU H2 Q=1 ENDPOINT SUPPORT THEOREM")
    print("PROVED raw gap >= %.12g on the exact feasible blocker"
          % report["endpoint_lower_float"])
    print("report_sha256 %s" % report["report_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
