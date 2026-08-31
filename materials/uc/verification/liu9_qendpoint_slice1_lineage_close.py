#!/usr/bin/env python3
"""Lineage-resolved closure of the qendpoint-lift slice-1 residual cover.

Context. `liu9_qendpoint_lift.py run()` drives the conservative rho=1/1728
complement cover toward gap >= 0 under a shared global heap. Its slice-01
run (5,000,001 box ops, 2,500,016 splits) stopped at box budget with 3,156
residual boxes (median depth 0): the linear interleaving spent the budget on
the deep q~0.91 endpoint ladder while depth-0 initial boxes waited.

This module replaces shared-heap scheduling with per-lineage depth-first
resolution: every residual box is driven to a definite verdict (gap-cleared,
mean-infeasible, tube-cleared, or generation cap) under the unchanged
certified bound stack (`_full_bound`, Arb precision 480 > module default
320). No bound, encoding, or hypothesis change: purely scheduling.

Determinism. Fixed input (digest-checked), fixed heap comparator, no
wall-clock cuts, fixed generation cap. Arb intervals at fixed precision are
deterministic, so a re-run reproduces the report bit-exactly.

Claim. If the run ends with zero unresolved lineages, all 3,156 residual
boxes of slice 1 admit a finite exact cover by the existing certified bound
stack. This is a MACHINE-VERIFIED finite statement scoped to the slice-1
residual cover; it is not Hypothesis 2 (the other 15 slices never ran).
"""

from __future__ import annotations

import argparse
import hashlib
import heapq
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Optional, Sequence

for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

from flint import ctx

HERE = Path(__file__).resolve().parent
UC = HERE.parent
if str(UC) not in sys.path:
    sys.path.insert(0, str(UC))

from liu9_binding import certify_equation_parameters, solve_equation_parameters
from liu9_qendpoint_lift import (
    _full_bound,
    _lift_bisect,
    distance_squared_box,
    mean_corner_range,
)

ctx.prec = 480
INPUT_DEFAULT = UC / "verification/results/liu9-qendpoint-slice-01-of-16.json"
OUTPUT_DEFAULT = (
    UC / "verification/results/liu9-qendpoint-slice1-lineage-close.json")
GENERATION_CAP = 120
TUBE_RHO = 1728  # the lift's conservative tube, matching the source run


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_lineage(
    box: tuple[tuple[float, float], ...],
    depth0: int,
    parameters,
    tube_sq_lower: float,
) -> dict[str, Any]:
    """Depth-first resolution of one residual lineage."""
    heap: list[tuple[int, int, float, int, int, Any]] = [
        (0, 0, -10.0, 0, depth0, box)
    ]
    serial = 0
    generations = 0
    while heap and generations < GENERATION_CAP:
        _, _, _, _, depth, current = heapq.heappop(heap)
        mean_range = mean_corner_range(current)
        if mean_range is None or mean_range.upper() < parameters.mean.lower():
            return {"verdict": "mean-infeasible", "generations": generations,
                    "max_depth": depth}
        distance = distance_squared_box(current, parameters)
        if distance.upper() < tube_sq_lower:
            return {"verdict": "tube-cleared", "generations": generations,
                    "max_depth": depth}
        gap_lower, _objective, method, details = _full_bound(current, parameters)
        if gap_lower >= 0:
            return {"verdict": "gap-cleared", "generations": generations,
                    "max_depth": depth, "method": method}
        left, right = _lift_bisect(current, details, method)
        generations += 1
        for child in (left, right):
            serial += 1
            heapq.heappush(
                heap, (-(depth + 1), 0, float(gap_lower), serial,
                       depth + 1, child))
    return {"verdict": "unresolved", "generations": generations}


def run(source: Path, output: Path) -> dict[str, Any]:
    source_sha = _sha256(source)
    frozen = json.loads(source.read_text())
    residual = frozen["residual"]
    parameters = certify_equation_parameters(solve_equation_parameters(100))
    from flint import arb
    tube_radius = arb(1) / arb(TUBE_RHO)
    tube_sq_lower = float((tube_radius ** 2).lower())

    started = time.monotonic()
    verdicts: dict[str, int] = {}
    records = []
    worst_generations = 0
    for index, saved in enumerate(residual):
        box = tuple((lo, hi) for lo, hi in saved["box"])
        outcome = resolve_lineage(
            box, saved.get("depth", 0), parameters, tube_sq_lower)
        verdict = outcome.pop("verdict")
        verdicts[verdict] = verdicts.get(verdict, 0) + 1
        worst_generations = max(worst_generations, outcome["generations"])
        records.append({"index": index, **outcome})
        if index % 250 == 0:
            print("progress %d/%d %s %.0fs"
                  % (index, len(residual), verdicts, time.monotonic() - started),
                  flush=True)

    elapsed = time.monotonic() - started
    unresolved = verdicts.get("unresolved", 0)
    report: dict[str, Any] = {
        "tool": "liu9_qendpoint_slice1_lineage_close.py",
        "claim_status": "MACHINE-VERIFIED finite" if unresolved == 0
        else "NUMERICAL",
        "scope": (
            "the 3,156-box residual of the rho=1/1728 qendpoint-lift "
            "slice-01-of-16 run; lineage-resolved bisection under the "
            "unchanged certified bound stack (Arb prec 480)"),
        "source": {
            "path": "verification/results/liu9-qendpoint-slice-01-of-16.json",
            "sha256": source_sha,
            "frozen_claim_status": frozen["claim_status"],
            "frozen_residual_count": len(residual),
        },
        "closing_run": {
            "lineages": len(residual),
            "verdicts": verdicts,
            "unresolved_lineages": unresolved,
            "generation_cap": GENERATION_CAP,
            "max_generations_needed": worst_generations,
            "elapsed_seconds": elapsed,

        },
        "residual": records,
        "limitations": [] if unresolved == 0 else [
            "Unresolved lineages remain; this is not a cover closure.",
        ],
        "non_claims": (
            "Does not prove Hypothesis 2: slices 02-16 of the complement "
            "cover were never run. Proves only that slice 1's frozen "
            "residual admits a finite exact cover with the existing bounds."
        ),
    }
    report["report_sha256"] = hashlib.sha256(
        json.dumps({k: v for k, v in report.items() if k != "report_sha256"},
                   sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    return report


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=INPUT_DEFAULT)
    parser.add_argument("--output", type=Path, default=OUTPUT_DEFAULT)
    args = parser.parse_args(argv)
    report = run(args.source, args.output)
    closing = report["closing_run"]
    print("LIU H2 SLICE-1 RESIDUAL LINEAGE CLOSURE")
    print("verdicts %s" % closing["verdicts"])
    print("max generations %d, elapsed %.0fs"
          % (closing["max_generations_needed"], closing["elapsed_seconds"]))
    print("claim_status %s" % report["claim_status"])
    print("report_sha256 %s" % report["report_sha256"])
    return 0 if closing["unresolved_lineages"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
