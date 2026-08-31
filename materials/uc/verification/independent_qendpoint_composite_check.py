#!/usr/bin/env python3
"""Independent authentication and cover check for the Liu H2 composite."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
SOURCE = RESULTS / "liu9-complement-residual-rho1-1728-100k.json"
LOCAL = RESULTS / "liu9-piecewise-tube.json"
ENDPOINT = RESULTS / "liu9-endpoint-support-independent.json"
SLICE_COUNT = 16
SLICES = tuple(
    RESULTS / (
        f"liu9-qendpoint-slice-{index:02d}-of-{SLICE_COUNT:02d}.json")
    for index in range(SLICE_COUNT)
)
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
EXPECTED_ENDPOINT_SHA256 = (
    "7b97e2d4e3936c248ef70f76431fbcab082186cf781a810341c9e59317221837")
OUTPUT = RESULTS / "liu9-qendpoint-composite-independent.json"


def _digest(payload: dict[str, Any]) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode()).hexdigest()


def _load(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    claimed = payload.pop("report_sha256", None)
    observed = _digest(payload)
    payload["report_sha256"] = claimed
    if claimed != observed:
        raise AssertionError(f"digest mismatch: {path.name}")
    return payload


def _box(saved: dict[str, Any]):
    return tuple((float(lo), float(hi)) for lo, hi in saved["box"])


def _active(box):
    return tuple(box[index] for index in (0, 1, 6, 7, 8))


def main() -> int:
    source = _load(SOURCE)
    boxes = [_box(saved) for saved in source["residual"]]
    groups = sorted({_active(box) for box in boxes})
    if len(boxes) != 49767 or len(groups) != 230:
        raise AssertionError("authenticated source cover has unexpected shape")
    if (
        source.get("report_type") != "liu9_complement_residual"
        or source.get("report_sha256") != EXPECTED_SOURCE_SHA256
    ):
        raise AssertionError("conservative complement source changed")

    local = _load(LOCAL)
    coverage = local.get("coverage", {})
    if (
        local.get("tool") != "liu9_tube.py"
        or local.get("claim_status") != "PROVED"
        or local.get("tube_radius") != "1/1701"
        or set(coverage) != EXPECTED_COVERAGE
        or any(coverage[key] is not True for key in EXPECTED_COVERAGE)
    ):
        raise AssertionError("local tube theorem is unavailable")
    if local["report_sha256"] != EXPECTED_LOCAL_SHA256:
        raise AssertionError("local tube digest changed")
    verified_local_components = 0
    for component in local["component_reports"].values():
        report = _load(RESULTS / component["path"])
        if report["report_sha256"] != component["report_sha256"]:
            raise AssertionError("local component binding mismatch")
        verified_local_components += 1

    endpoint = _load(ENDPOINT)
    if (
        endpoint.get("claim_status") != "PASS"
        or endpoint.get("report_sha256") != EXPECTED_ENDPOINT_SHA256
    ):
        raise AssertionError("independent endpoint arithmetic did not pass")

    reports = [_load(path) for path in SLICES]
    indices = {report["source"]["slice_index"] for report in reports}
    if indices != set(range(SLICE_COUNT)):
        raise AssertionError("slice index set is incomplete")
    covered = 0
    processed = 0
    for report in reports:
        meta = report["source"]
        index = meta["slice_index"]
        selected = set(groups[index::SLICE_COUNT])
        expected = sum(_active(box) in selected for box in boxes)
        if (
            meta["slice_count"] != SLICE_COUNT
            or meta["slice_source_boxes"] != expected
            or meta["active_groups_selected"] != len(selected)
            or meta["report_sha256"] != source["report_sha256"]
        ):
            raise AssertionError(f"slice {index} partition metadata failed")
        if (
            report.get("tool") != "liu9_qendpoint_lift.py"
            or report.get("rho") != "1/1728"
        ):
            raise AssertionError(f"slice {index} semantic identity failed")
        if (
            report.get("claim_status") != "PROVED"
            or report.get("residual")
            or report["run"].get("stop_reason") != "complete"
        ):
            raise AssertionError(f"slice {index} is incomplete")
        covered += expected
        processed += report["run"]["processed"]
    if covered != len(boxes):
        raise AssertionError("slice cover is not exact")

    result = {
        "tool": "independent_qendpoint_composite_check.py",
        "claim_status": "PASS",
        "source_report_sha256": source["report_sha256"],
        "local_tube_report_sha256": local["report_sha256"],
        "independent_endpoint_report_sha256": endpoint["report_sha256"],
        "verified_local_components": verified_local_components,
        "slice_count": SLICE_COUNT,
        "active_groups": len(groups),
        "covered_source_boxes": covered,
        "processed": processed,
        "residual": 0,
        "conclusion": (
            "authenticated local tube plus exact disjoint complement slices "
            "cover every mean-feasible point"
        ),
    }
    result["report_sha256"] = _digest(result)
    OUTPUT.write_text(json.dumps(result, indent=1, sort_keys=True) + "\n")
    print("INDEPENDENT LIU H2 COMPOSITE CHECK")
    print("PASS slices=%d source_boxes=%d residual=0"
          % (SLICE_COUNT, covered))
    print("report_sha256 %s" % result["report_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
