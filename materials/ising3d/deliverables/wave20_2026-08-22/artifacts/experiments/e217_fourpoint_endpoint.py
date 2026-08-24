#!/usr/bin/env python3
"""Assemble the exact mode-resolved four-point method-limitation artifact."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
EXPERIMENTS = ROOT / "experiments"
if str(EXPERIMENTS) not in sys.path:
    sys.path.insert(0, str(EXPERIMENTS))

import e215_fourpoint_kernel as kernel_component  # noqa: E402
import e216_fourpoint_certificate as certificate_component  # noqa: E402

SCRIPT = "experiments/e217_fourpoint_endpoint.py"
OUTPUT = ROOT / "results" / "bounds" / "upper_fourpoint.json"
INPUTS = (
    ROOT / "proofs" / "upper_endpoint4.md",
    ROOT / "proofs" / "kc_upper_peierls.md",
    ROOT / "proofs" / "upper_infrared.md",
    ROOT / "proofs" / "mag_floor.md",
    ROOT / "experiments" / "e93_mag_floor.py",
    ROOT / "experiments" / "e215_fourpoint_kernel.py",
    ROOT / "experiments" / "e216_fourpoint_certificate.py",
    ROOT / "experiments" / "e217_fourpoint_endpoint.py",
)
COVERAGE = (
    "RATIONAL_CHALLENGE_BELOW_WATSON",
    "EXACT_RATIONAL_TORUS_KERNELS",
    "MODE_RESOLVED_PHYSICAL_MAPPING",
    "FOURPOINT_MOMENT_ROWS",
    "DETERMINISTIC_LIFT_PROJECTION_EQUALITY",
    "EXACT_L12_ETA_COUNTERCERTIFICATE",
    "ALL_SIZE_VANISHING_FLOOR",
    "THERMODYNAMIC_DIRECTION_AND_SCOPE",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _record(
    checks: list[dict[str, Any]], name: str, condition: bool, detail: str
) -> None:
    passed = bool(condition)
    if not passed:
        raise AssertionError(f"{name}: {detail}")
    checks.append({"name": name, "passed": passed, "detail": detail})


def main() -> int:
    kernel_data, kernel_checks = kernel_component.build_kernel_data()
    certificate_data, certificate_checks = certificate_component.build_certificate_data(
        kernel_data
    )

    checks = [
        {
            "name": f"kernel::{row['name']}",
            "passed": bool(row["passed"]),
            "detail": row["detail"],
        }
        for row in kernel_checks
    ]
    checks.extend(
        {
            "name": f"certificate::{row['name']}",
            "passed": bool(row["passed"]),
            "detail": row["detail"],
        }
        for row in certificate_checks
    )

    kernel_sides = [row["side"] for row in kernel_data["sides"]]
    lift_sides = [row["side"] for row in certificate_data["finite_lifts"]]
    _record(
        checks,
        "cross::kernel_lift_sides",
        kernel_sides == lift_sides == list(kernel_component.SIDES),
        f"kernel and four-point lifts share sides {kernel_sides}",
    )
    _record(
        checks,
        "cross::coverage_unique",
        len(COVERAGE) == len(set(COVERAGE)),
        ",".join(COVERAGE),
    )
    _record(
        checks,
        "cross::all_component_checks_pass",
        all(row["passed"] for row in checks),
        f"{len(checks)} exact component and cross checks pass",
    )

    provenance: dict[str, dict[str, Any]] = {}
    for path in INPUTS:
        stat = path.stat()
        provenance[str(path.relative_to(ROOT))] = {
            "sha256": _sha256(path),
            "size_bytes": stat.st_size,
        }

    meta = {
        "producer": SCRIPT,
        "components": [kernel_component.SCRIPT, certificate_component.SCRIPT],
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "provenance": provenance,
        "environment": {
            "python": sys.version,
            "executable": sys.executable,
            "implementation": platform.python_implementation(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "hash_seed": os.environ.get("PYTHONHASHSEED", "unset"),
        },
        "arithmetic": (
            "decisive values use integers, Fraction, or exact quadratic fields; "
            "finite Green kernels stored in the artifact are rational"
        ),
        "benchmark_policy": (
            "K_c=0.221654626 has no selection, fitting, comparison, or validation "
            "role; K=6/25 is a pre-fixed plain rational challenge"
        ),
    }
    data = {
        "headline": (
            "[THEOREM] MR4-power-simplex-L1 has exactly the same p0 optimum as "
            "the inherited two-point relaxation and its certifiable uniform floor "
            "vanishes; it cannot improve the Watson endpoint."
        ),
        "classification": {
            "endpoint": "[UNRESOLVED] the certified upper endpoint remains I3/2",
            "negative_result": (
                "[THEOREM] exact method limitation for MR4-power-simplex-L1"
            ),
            "finite_data": (
                "[COMPUTATION] exact rational Green kernels and exact algebraic "
                "mode lifts on L=4,6,8,12"
            ),
        },
        "coverage": list(COVERAGE),
        "kernel_certificates": kernel_data,
        "fourpoint_certificate": certificate_data,
        "outcome": {
            "improved_endpoint": False,
            "incumbent": "I3/2",
            "challenge": "K=6/25",
            "challenge_strictly_below_incumbent": True,
            "method_class": certificate_component.CLASS_NAME,
            "exact_optimum_statement": (
                "inf_{(p,Q) in MR4} p0 = inf_{p in two-point relaxation} p0 "
                "at every finite L and K"
            ),
            "thermodynamic_limitation": (
                "for every 0<K<=I3/2 the class infimum has a feasible sequence "
                "with p0=O(1/L), so it cannot certify a uniform positive floor"
            ),
            "counterexample_scope": certificate_data["scope"]["not_covered"],
        },
    }
    payload = {"meta": meta, "data": data, "checks": checks}
    _record(
        checks,
        "cross::artifact_shape",
        set(payload) == {"meta", "data", "checks"},
        "top-level artifact envelope is exactly meta/data/checks",
    )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    for check in checks:
        print(f"PASS {check['name']} -- {check['detail']}")
    print(f"wrote {OUTPUT}")
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
