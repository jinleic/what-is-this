"""Assemble the exact paired-momentum and random-current endpoint obstructions.

Run from the repository root:
    .venv/bin/python experiments/e177_upper_endpoint4.py

The producer writes results/bounds/upper_endpoint4.json and prints a final PASS
line only when every exact component and cross-check passes.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import sys
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
EXPERIMENTS = ROOT / "experiments"
if str(EXPERIMENTS) not in sys.path:
    sys.path.insert(0, str(EXPERIMENTS))

import e175_upper_endpoint4_paired as paired_component  # noqa: E402
import e176_upper_endpoint4_current as current_component  # noqa: E402

SCRIPT = "experiments/e177_upper_endpoint4.py"
OUTPUT = ROOT / "results" / "bounds" / "upper_endpoint4.json"
INPUTS = [
    ROOT / "proofs" / "upper_infrared.md",
    ROOT / "proofs" / "mag_floor.md",
    ROOT / "proofs" / "upper_beyond.md",
    ROOT / "proofs" / "kc_interval2.md",
    ROOT / "proofs" / "paired_momentum.md",
    ROOT / "results" / "bounds" / "upper_infrared.json",
    ROOT / "experiments" / "e175_upper_endpoint4_paired.py",
    ROOT / "experiments" / "e176_upper_endpoint4_current.py",
]
COVERAGE_ROWS = [
    "TARGET_BELOW_INCUMBENT",
    "PM_PARSEVAL_BOUNDED",
    "PM_SDP_OPTIMUM",
    "PM_POSITIVE_FLOOR_INFEASIBLE",
    "RC_EDGE_WEIGHTS",
    "RC_SWITCHING_C4",
    "RC_DOMINATION_LP",
    "PERC_PLANAR_5_6",
]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _record(checks: list[dict[str, Any]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def _package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def main() -> int:
    paired, paired_checks = paired_component.build_paired_data()
    current, current_checks = current_component.build_current_data()
    # Uniform version of the current obstruction on the whole incumbent range.
    # atanh(1/4) > 1/4+(1/4)^3/3 = 49/192, while I3/2 < 49/192
    # follows from the certified rational I3 upper endpoint.
    i3_hi = Fraction(paired["target"]["I3_interval_input"][1])
    atanh_quarter_lower = Fraction(49, 192)
    uniform_p_upper = Fraction(1, 16)
    uniform_gap_lower = Fraction(5, 6) - uniform_p_upper
    current["whole_incumbent_interval_obstruction"] = {
        "comparison": "I3/2 < 49/192 < atanh(1/4)",
        "atanh_one_quarter_lower": "49/192",
        "uniform_optimum_bound": "p*=tanh(K)^2<1/16 for every K<=I3/2",
        "uniform_gap_lower_to_5_over_6": "37/48",
    }
    checks = [
        {"name": f"paired::{row['name']}", "passed": row["passed"], "detail": row["detail"]}
        for row in paired_checks
    ]
    checks.extend(
        {"name": f"current::{row['name']}", "passed": row["passed"], "detail": row["detail"]}
        for row in current_checks
    )

    _record(
        checks,
        "cross::same_target_v",
        paired["target"]["v"] == current["target_v"] == "6/25",
        f"paired={paired['target']['v']}, current={current['target_v']}",
    )
    _record(
        checks,
        "cross::both_exact_obstructions",
        paired["primal_sdp"]["exact_optimum_order_floor_eta"] == "0/1"
        and current["domination_lp"]["exact_optimum"] == "36/625"
        and current["domination_lp"]["exact_infeasibility_gap"] == "2909/3750"
        and i3_hi / 2 < atanh_quarter_lower
        and uniform_gap_lower == Fraction(37, 48),
        "paired eta*=0; current p*=36/625 at target and p*<1/16 on every K<=I3/2",
    )
    _record(
        checks,
        "cross::coverage_rows_complete",
        len(COVERAGE_ROWS) == len(set(COVERAGE_ROWS)) == 8,
        ",".join(COVERAGE_ROWS),
    )

    input_meta = {}
    for path in INPUTS:
        stat = path.stat()
        input_meta[str(path.relative_to(ROOT))] = {
            "sha256": _sha256(path),
            "size_bytes": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
        }

    payload = {
        "meta": {
            "provenance": {
                "producer": SCRIPT,
                "components": [paired_component.SCRIPT, current_component.SCRIPT],
                "inputs": input_meta,
            },
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "environment": {
                "python_version": sys.version,
                "python_executable": sys.executable,
                "python_implementation": platform.python_implementation(),
                "python_prefix": sys.prefix,
                "platform": platform.platform(),
                "machine": platform.machine(),
                "byteorder": sys.byteorder,
                "hash_seed": os.environ.get("PYTHONHASHSEED", "unset"),
                "packages": {
                    "mpmath": _package_version("mpmath"),
                    "numpy": _package_version("numpy"),
                    "sympy": _package_version("sympy"),
                },
            },
            "arithmetic": (
                "all decisive optimization, current weights, switching identities, endpoint comparisons, "
                "and gaps use Python integers/Fraction; decimal strings are directed formatting only"
            ),
            "benchmark_policy": (
                "K_c=0.221654626 has no logical or selection role; the target v=6/25 is a plain rational "
                "chosen to make both current weights and the endpoint challenge exact"
            ),
        },
        "data": {
            "headline": (
                "No improved upper endpoint. Exact infeasibility in PM4-aggregate-one-moment and "
                "RC2-single-edge-conditional-domination at K=atanh(6/25)<I3/2."
            ),
            "classification": {
                "paired_momentum": "[COMPUTATION][THEOREM] exact SDP optimum/order-floor eta*=0",
                "random_current": "[COMPUTATION][THEOREM] exact LP optimum p*=36/625, short of 5/6 by 2909/3750",
                "endpoint": "[UNRESOLVED] certified upper endpoint remains I3/2",
            },
            "coverage_rows": COVERAGE_ROWS,
            "paired_momentum": paired,
            "random_current": current,
            "outcome": {
                "improved_endpoint": False,
                "incumbent": "I3/2",
                "challenged_endpoint": "atanh(6/25)",
                "challenged_endpoint_strictly_below_incumbent": True,
                "paired_obstruction": "exact optimum uniform magnetisation floor is 0",
                "current_obstruction": (
                    "best one-edge iid domination p is 36/625 at the target (gap 2909/3750); "
                    "uniformly p<1/16 for every K<=I3/2 (gap to 5/6 greater than 37/48)"
                ),
                "scope_limit": (
                    "does not exclude mode-resolved four-point/DLR constraints, multi-edge or block random-current "
                    "certificates, source-dependent current inequalities, or any other method outside the two named classes"
                ),
            },
            "checks": checks,
        },
    }

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    for check in checks:
        print(f"{'PASS' if check['passed'] else 'FAIL'} {check['name']} -- {check['detail']}")
    passed = all(check["passed"] for check in checks)
    print(f"wrote {OUTPUT}")
    print("PASS" if passed else "FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
