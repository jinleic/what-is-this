#!/usr/bin/env python3
"""Localize the residual frontier of the Liu-H2 tube complement.

`liu9_tube.py` measures how much of the nine-parameter complement
`{dist >= rho}` an interval branch-and-bound can clear.  At `rho=0.1` a 780 s
single-core run processed 866,551 boxes and still left 76 boxes on the heap,
against 71 after 26,247 boxes in the earlier short run.  A frontier that does
not shrink under a 33x budget increase is not a volume problem, so the useful
question is *where* those boxes live.

This script reruns the same complement branch-and-bound with the same rules
(`liu9_tube.complement_box_bound`, `distance_squared_box`, `mean_corner_range`,
identical priority order and bisection) but keeps the residual heap and reports
its geometry: how many residual boxes straddle the tube boundary, how many
touch the degenerate mixture endpoints `q in {0,1}` where the transverse
curvature `C q (1-q)` vanishes, which support coordinates are pinned at the
entropy-singular corner `0`, and the distribution of gap lower bounds.

Status: NUMERICAL/diagnostic.  Every individual discard inside the run is
certified by exact dyadic geometry plus Arb, but no verdict here is a proof
about Hypothesis 2.  The stop is by box budget, which is deterministic, so the
reported frontier is reproducible.

Usage:
    python -I -B uc/liu9_residual.py --rho 0.1 --box-budget 866551 \
        [--output PATH]
"""

from __future__ import annotations

import argparse
import datetime as _datetime
import hashlib
import heapq
import json
import os
import sys
import time
import uuid
from fractions import Fraction
from pathlib import Path
from typing import Optional, Sequence

for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

from flint import arb, ctx  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from liu9_binding import (  # noqa: E402
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_objective import arb_fraction  # noqa: E402
from liu9_pilot import Box, ONE, ZERO, mean_corner_range  # noqa: E402
from liu9_tube import (  # noqa: E402
    _bisect,
    complement_box_bound,
    distance_squared_box,
)

ctx.prec = max(ctx.prec, 320)

COORDINATES = ("a1", "a2", "q", "b0", "b2", "b4", "b1", "b3", "b5")
SUPPORT_INDICES = tuple(range(3, 9))


def collect_residual(
    rho: Fraction,
    parameters,
    lambdas: Sequence[Fraction],
    box_budget: int,
    seconds: float,
) -> dict[str, object]:
    rho_squared = arb_fraction(rho * rho)
    root: Box = tuple((0.0, 1.0) for _ in range(9))
    heap: list = []
    serial = 0
    processed = 0
    split = 0
    mean_infeasible = 0
    tube_excluded = 0
    objective_cleared = 0
    maximum_depth = 0
    started = time.monotonic()
    deadline = started + seconds
    stop_reason = "complete"

    def classify(box: Box, depth: int) -> None:
        nonlocal serial, processed, mean_infeasible, tube_excluded
        nonlocal objective_cleared, maximum_depth
        processed += 1
        maximum_depth = max(maximum_depth, depth)
        mean_range = mean_corner_range(box)
        if mean_range is None or mean_range.upper() < parameters.mean.lower():
            mean_infeasible += 1
            return
        distance_range = distance_squared_box(box, parameters)
        if distance_range.upper() < rho_squared:
            tube_excluded += 1
            return
        straddles_tube = not distance_range.lower() >= rho_squared
        bound = complement_box_bound(box, parameters, lambdas)
        if bound.gap_lower >= ZERO or bound.objective_lower >= ONE:
            objective_cleared += 1
            return
        serial += 1
        geometry_priority = (
            float(distance_range.upper() - rho_squared)
            if straddles_tube else bound.key
        )
        heapq.heappush(
            heap,
            (0 if straddles_tube else 1, geometry_priority, bound.key,
             serial, depth, box, bound, distance_range),
        )

    classify(root, 0)
    while heap:
        if processed + 2 > box_budget:
            stop_reason = "box_budget"
            break
        if time.monotonic() >= deadline:
            stop_reason = "wall_clock"
            break
        _, _, _, _, depth, box, _, _ = heapq.heappop(heap)
        split += 1
        left, right = _bisect(box)
        classify(left, depth + 1)
        classify(right, depth + 1)

    elapsed = time.monotonic() - started
    residual = []
    for entry in heap:
        straddling = entry[0] == 0
        depth = entry[4]
        box = entry[5]
        bound = entry[6]
        distance_range = entry[7]
        widths = [upper - lower for lower, upper in box]
        residual.append({
            "box": [[lower, upper] for lower, upper in box],
            "depth": depth,
            "distance_squared_lower": float(distance_range.lower()),
            "distance_squared_upper": float(distance_range.upper()),
            "gap_lower": float(bound.gap_lower),
            "maximum_width": max(widths),
            "method": bound.method,
            "objective_lower": float(bound.objective_lower),
            "q_high": box[2][1],
            "q_low": box[2][0],
            "straddles_tube": straddling,
            "supports_touching_zero": [
                COORDINATES[index] for index in SUPPORT_INDICES
                if box[index][0] <= 0.0
            ],
        })
    residual.sort(key=lambda item: item["gap_lower"])
    return {
        "elapsed_seconds": round(elapsed, 3),
        "maximum_depth": maximum_depth,
        "mean_infeasible": mean_infeasible,
        "objective_cleared": objective_cleared,
        "processed": processed,
        "residual": residual,
        "residual_count": len(heap),
        "split": split,
        "stop_reason": stop_reason,
        "tube_excluded": tube_excluded,
    }


def summarize(result: dict[str, object], rho: Fraction) -> dict[str, object]:
    residual = result["residual"]
    count = len(residual)
    if not count:
        return {"residual_count": 0}
    straddling = sum(1 for item in residual if item["straddles_tube"])
    q_at_zero = sum(1 for item in residual if item["q_low"] <= 0.0)
    q_at_one = sum(1 for item in residual if item["q_high"] >= 1.0)
    q_full = sum(
        1 for item in residual if item["q_low"] <= 0.0 and item["q_high"] >= 1.0
    )
    q_interior = sum(
        1 for item in residual if item["q_low"] > 0.0 and item["q_high"] < 1.0
    )
    zero_support = {}
    per_coordinate = {COORDINATES[index]: 0 for index in SUPPORT_INDICES}
    for item in residual:
        key = ",".join(item["supports_touching_zero"]) or "(none)"
        zero_support[key] = zero_support.get(key, 0) + 1
        for name in item["supports_touching_zero"]:
            per_coordinate[name] += 1
    any_zero_support = sum(
        1 for item in residual if item["supports_touching_zero"]
    )
    return {
        "any_support_touching_zero": any_zero_support,
        "distinct_zero_support_patterns": len(zero_support),
        "gap_lower_max": max(item["gap_lower"] for item in residual),
        "gap_lower_min": min(item["gap_lower"] for item in residual),
        "maximum_width_max": max(item["maximum_width"] for item in residual),
        "maximum_width_min": min(item["maximum_width"] for item in residual),
        "q_at_one_only": q_at_one - q_full,
        "q_at_zero_only": q_at_zero - q_full,
        "q_full_range": q_full,
        "q_interior": q_interior,
        "residual_count": count,
        "straddling_tube_boundary": straddling,
        "support_touching_zero_by_coordinate": per_coordinate,
        "wholly_outside_tube": count - straddling,
        "zero_support_patterns": dict(
            sorted(zero_support.items(), key=lambda pair: -pair[1])
        ),
    }


def canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def write_report(path: Path, report: dict[str, object]) -> None:
    path = path.expanduser().resolve()
    path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    payload = json.dumps(report, indent=2, sort_keys=True, ensure_ascii=True) + "\n"
    descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    try:
        view = memoryview(payload.encode("ascii"))
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise RuntimeError(f"short report write: {path}")
            view = view[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    os.replace(temporary, path)


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rho", default="0.1")
    parser.add_argument("--box-budget", type=int, default=866_551)
    parser.add_argument("--seconds", type=float, default=3000.0)
    parser.add_argument("--prec", type=int, default=320)
    parser.add_argument("--dps", type=int, default=90)
    parser.add_argument(
        "--lambdas", nargs="+", default=("0", "-0.4", "-0.8", "-0.9", "-1.0", "-1.2")
    )
    parser.add_argument("--output", type=Path, default=None)
    arguments = parser.parse_args(argv)

    rho = Fraction(arguments.rho)
    if not 0 < rho < 1:
        parser.error("--rho must lie strictly between zero and one")
    lambdas = tuple(Fraction(value) for value in arguments.lambdas)
    if any(value > 0 for value in lambdas):
        parser.error("all feasibility-shift lambdas must be <= 0")
    if arguments.box_budget < 101:
        parser.error("--box-budget must be at least 101")

    ctx.prec = arguments.prec
    mp_parameters = solve_equation_parameters(arguments.dps + 20)
    arb_parameters = certify_equation_parameters(mp_parameters)

    print("LIU H2 COMPLEMENT RESIDUAL LOCALIZER")
    print("NUMERICAL [scope]: diagnostic geometry of the surviving frontier;")
    print("   individual discards are Arb-certified, the verdict is not a proof.")
    print(
        "NUMERICAL [budget]: rho=%s, box budget=%d, time cap=%.0fs, prec=%d bits."
        % (arguments.rho, arguments.box_budget, arguments.seconds, arguments.prec)
    )
    result = collect_residual(
        rho, arb_parameters, lambdas, arguments.box_budget, arguments.seconds
    )
    summary = summarize(result, rho)
    print(
        "NUMERICAL [run]: processed=%d split=%d mean-infeasible=%d tube-out=%d "
        "obj-clear=%d residual=%d depth=%d elapsed=%.1fs stop=%s"
        % (
            result["processed"], result["split"], result["mean_infeasible"],
            result["tube_excluded"], result["objective_cleared"],
            result["residual_count"], result["maximum_depth"],
            result["elapsed_seconds"], result["stop_reason"],
        )
    )
    if summary["residual_count"]:
        print(
            "NUMERICAL [frontier]: straddling tube boundary=%d, wholly outside=%d;"
            % (summary["straddling_tube_boundary"], summary["wholly_outside_tube"])
        )
        print(
            "   q pinned at 0 only = %d, at 1 only = %d, full [0,1] = %d, interior = %d;"
            % (
                summary["q_at_zero_only"], summary["q_at_one_only"],
                summary["q_full_range"], summary["q_interior"],
            )
        )
        print(
            "   boxes with a support coordinate touching 0 = %d (%d distinct patterns);"
            % (
                summary["any_support_touching_zero"],
                summary["distinct_zero_support_patterns"],
            )
        )
        print(
            "   gap lower bound in [%.6g, %.6g]; widest coordinate in [%.6g, %.6g]."
            % (
                summary["gap_lower_min"], summary["gap_lower_max"],
                summary["maximum_width_min"], summary["maximum_width_max"],
            )
        )
        print(
            "   per-coordinate zero contact: %s"
            % summary["support_touching_zero_by_coordinate"]
        )
    else:
        print("NUMERICAL [frontier]: empty; the complement closed at this budget.")

    report = {
        "claim_status": "NUMERICAL",
        "finished_utc": _datetime.datetime.now(_datetime.timezone.utc).isoformat(),
        "lambdas": [str(value) for value in lambdas],
        "limitations": [
            "Diagnostic only: no verdict here is a proof about Liu Hypothesis 2.",
            "The local tube lemma remains unproved, so the complement run is counterfactual.",
            "Stop is by deterministic box budget; a different budget yields a different frontier.",
        ],
        "precision_bits": arguments.prec,
        "report_type": "liu9_complement_residual",
        "rho": str(rho),
        "run": {key: value for key, value in result.items() if key != "residual"},
        "residual": result["residual"],
        "summary": summary,
    }
    report["report_sha256"] = hashlib.sha256(
        canonical_bytes({k: v for k, v in report.items()})
    ).hexdigest()
    if arguments.output is not None:
        write_report(arguments.output, report)
        print("NUMERICAL [report]: %s" % arguments.output)


if __name__ == "__main__":
    main()
