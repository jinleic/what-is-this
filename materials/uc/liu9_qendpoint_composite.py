#!/usr/bin/env python3
"""Aggregate disjoint q-endpoint certificate slices into the global theorem."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Optional, Sequence

HERE = Path(__file__).resolve().parent
SLICE_COUNT = 16
DEFAULT_INPUTS = tuple(
    HERE / (
        f"verification/results/liu9-qendpoint-slice-{index:02d}"
        f"-of-{SLICE_COUNT:02d}.json")
    for index in range(SLICE_COUNT)
)
OUTPUT_DEFAULT = HERE / "verification/results/liu9-qendpoint-composite.json"
SOURCE = HERE / "verification/results/liu9-complement-residual-rho1-1728-100k.json"

LOCAL_TUBE = HERE / "verification/results/liu9-piecewise-tube.json"
EXPECTED_COVERAGE = frozenset({
    "interior_smooth_chart",
    "mirror_support_boundary",
    "q_degenerate_endpoints",
    "q_seam_complete",
    "quotient_gauges_removed",
    "simultaneous_second_order_insertions",
    "strict_mean_half_space",
    "zero_support_boundary",
})
EXPECTED_SOURCE_SHA256 = (
    "494afc93c37eefb07cd0a8e0c982746db020e1387e0c93ab7b731fadd4524b69")
EXPECTED_LOCAL_SHA256 = (
    "4651ad293d5205a9ba416a3474fbc0b4ee6e553c762d4790f2ce7ccc83612742")

def _canonical_digest(payload: dict[str, Any]) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode()).hexdigest()


def _load_verified(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    claimed = payload.pop("report_sha256", None)
    observed = _canonical_digest(payload)
    payload["report_sha256"] = claimed
    if claimed != observed:
        raise AssertionError(f"report digest mismatch: {path}")
    return payload


def _box(saved: dict[str, Any]) -> tuple[tuple[float, float], ...]:
    return tuple((float(lo), float(hi)) for lo, hi in saved["box"])


def _active_key(box) -> tuple[tuple[float, float], ...]:
    return tuple(box[index] for index in (0, 1, 6, 7, 8))


def aggregate(inputs: Sequence[Path]) -> dict[str, Any]:
    source = _load_verified(SOURCE)
    source_boxes = [_box(saved) for saved in source["residual"]]
    active_groups = sorted({_active_key(box) for box in source_boxes})
    if len(source_boxes) != 49767 or len(active_groups) != 230:
        raise AssertionError("source residual identity changed")
    if (
        source.get("report_type") != "liu9_complement_residual"
        or source.get("report_sha256") != EXPECTED_SOURCE_SHA256
    ):
        raise AssertionError("conservative complement source is not authenticated")

    reports = [_load_verified(path) for path in inputs]
    slice_counts = {report["source"]["slice_count"] for report in reports}
    local_tube = _load_verified(LOCAL_TUBE)
    coverage = local_tube.get("coverage", {})
    if (
        local_tube.get("tool") != "liu9_tube.py"
        or local_tube.get("claim_status") != "PROVED"
        or local_tube.get("tube_radius") != "1/1701"
        or set(coverage) != EXPECTED_COVERAGE
        or any(coverage[key] is not True for key in EXPECTED_COVERAGE)
    ):
        raise AssertionError("proved rho=1/1701 local tube is unavailable")
    if local_tube["report_sha256"] != EXPECTED_LOCAL_SHA256:
        raise AssertionError("local tube digest changed")
    q_degenerate = local_tube.get("component_reports", {}).get("q_degenerate", {})
    if q_degenerate.get("report_sha256") != (
        "353131d9124516841fe32da2a8f09829d401f9304890a604d91ff544ca693107"
    ):
        raise AssertionError("near-maximal q-degenerate seam is not authenticated")
    if len(slice_counts) != 1:
        raise AssertionError("slice counts disagree")
    slice_count = slice_counts.pop()
    if slice_count != len(reports):
        raise AssertionError("incomplete slice set")
    indices = {report["source"]["slice_index"] for report in reports}
    if indices != set(range(slice_count)):
        raise AssertionError("slice indices do not partition the cover")

    expected_total = 0
    total_processed = total_split = total_gap = total_mean = total_tube = 0
    component_digests = []
    for report, path in zip(reports, inputs):
        source_meta = report["source"]
        index = source_meta["slice_index"]
        selected = set(active_groups[index::slice_count])
        expected_boxes = sum(
            _active_key(box) in selected for box in source_boxes)
        if source_meta["slice_source_boxes"] != expected_boxes:
            raise AssertionError(f"slice {index} source-box count mismatch")
        if source_meta["active_groups_selected"] != len(selected):
            raise AssertionError(f"slice {index} active-group count mismatch")
        if source_meta["report_sha256"] != source["report_sha256"]:
            raise AssertionError(f"slice {index} source digest mismatch")
        if (
            report.get("tool") != "liu9_qendpoint_lift.py"
            or report.get("rho") != "1/1728"
        ):
            raise AssertionError(f"slice {index} semantic identity mismatch")
        if report["claim_status"] != "PROVED" or report["residual"]:
            raise AssertionError(f"slice {index} is not complete")
        if report["run"]["stop_reason"] != "complete":
            raise AssertionError(f"slice {index} did not terminate by exhaustion")
        expected_total += expected_boxes
        total_processed += report["run"]["processed"]
        total_split += report["run"]["split"]
        total_gap += report["run"]["gap_cleared"]
        total_mean += report["run"]["mean_infeasible"]
        total_tube += report["run"]["tube_cleared"]
        component_digests.append({
            "path": str(path.relative_to(HERE)),
            "slice_index": index,
            "report_sha256": report["report_sha256"],
        })
    if expected_total != len(source_boxes):
        raise AssertionError("slice partition does not cover every source box")

    component_digests.sort(key=lambda item: item["slice_index"])
    composite: dict[str, Any] = {
        "tool": "liu9_qendpoint_composite.py",
        "claim_status": "PROVED",
        "source": {
            "path": str(SOURCE.relative_to(HERE)),
            "report_sha256": source["report_sha256"],
            "residual_boxes": len(source_boxes),
            "active_groups": len(active_groups),
        },
        "partition": {
            "slice_count": slice_count,
            "covered_source_boxes": expected_total,
            "components": component_digests,
        },
        "totals": {
            "processed": total_processed,
            "split": total_split,
            "gap_cleared": total_gap,
            "mean_infeasible": total_mean,
            "tube_cleared": total_tube,
            "residual": 0,
        },
        "local_tube": {
            "path": str(LOCAL_TUBE.relative_to(HERE)),
            "report_sha256": local_tube["report_sha256"],
            "tube_radius": local_tube["tube_radius"],
            "q_degenerate_report_sha256": q_degenerate["report_sha256"],
        },
        "theorem": (
            "the proved rho=1/1701 local tube and the rigorously cleared "
            "rho=1/1728 complement cover imply Liu Hypothesis 2 unconditionally"
        ),
        "raw_gap_claim": "gap >= 0 on every exact mean-feasible point",
        "quotient_used_at_entropy_zero": False,
    }
    composite["report_sha256"] = _canonical_digest(composite)
    return composite


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="*", type=Path, default=DEFAULT_INPUTS)
    parser.add_argument("--output", type=Path, default=OUTPUT_DEFAULT)
    args = parser.parse_args(argv)
    report = aggregate(args.inputs or DEFAULT_INPUTS)
    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    print("LIU H2 Q-ENDPOINT COMPOSITE CERTIFICATE")
    print("PROVED source boxes=%d residual=0 processed=%d"
          % (report["source"]["residual_boxes"], report["totals"]["processed"]))
    print("report_sha256 %s" % report["report_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
