#!/usr/bin/env python3
"""Integrate the c8 closed formula and c10 connected-object obstruction."""

from __future__ import annotations

import gc
import hashlib
import json
import platform
import resource
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from e221_c8_closed import build_c8_closed
from e222_c10_structure import build_c10_structure

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "interlayer" / "c10_structure.json"
CPU_BUDGET_SECONDS = 900.0
RSS_CAP_BYTES = 1_900_000_000
SOURCE_PATHS = (
    ROOT / "experiments" / "e221_c8_closed.py",
    ROOT / "experiments" / "e222_c10_structure.py",
    Path(__file__),
    ROOT / "results" / "interlayer" / "c8_series.json",
    ROOT / "results" / "interlayer" / "c6_closed.json",
    ROOT / "proofs" / "interlayer_c8.md",
    ROOT / "proofs" / "c6_closed.md",
    ROOT / "proofs" / "interlayer_allorders.md",
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def record(
    checks: list[dict[str, object]], name: str, passed: bool, detail: str
) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def guard(started: float, stage: str) -> None:
    elapsed = time.process_time() - started
    rss = max_rss_bytes()
    if elapsed >= CPU_BUDGET_SECONDS:
        raise RuntimeError(
            f"process-time budget exceeded at {stage}: {elapsed:.3f}s"
        )
    if rss >= RSS_CAP_BYTES:
        raise RuntimeError(f"RSS cap exceeded at {stage}: {rss}")


def run_certificate() -> tuple[dict[str, object], list[dict[str, object]]]:
    started = time.process_time()
    if not all(path.is_file() for path in SOURCE_PATHS):
        missing = [str(path) for path in SOURCE_PATHS if not path.is_file()]
        raise FileNotFoundError(", ".join(missing))

    c8, c8_checks = build_c8_closed()
    guard(started, "c8 closed reconstruction")
    gc.collect()
    c10, c10_checks = build_c10_structure()
    guard(started, "c10 partition inventory")

    checks = list(c8_checks) + list(c10_checks)
    c8_profiles = c8["connected_formula"]["profiles"]
    c10_profiles = c10["connected_formula"]["profiles"]
    record(
        checks,
        "integrated c8 and c10 composition counts",
        len(c8_profiles) == 8 and len(c10_profiles) == 16,
        "2^3=8 and 2^4=16",
    )
    record(
        checks,
        "integrated top-connected square obstruction",
        c8["w8_sector"]["w8_square_coefficient"] == "1/40320"
        and c10["w10_sector"]["w10_square_coefficient"] == "1/3628800",
        "1/8! and 1/10!",
    )
    record(
        checks,
        "integrated connected alphabets",
        {kind for row in c8_profiles for key in row["atom_profile_histogram"] for kind in key.split(",")}
        == {"G", "U", "W", "W8"}
        and max(row["maximum_connected_rank"] for row in c10_profiles) == 10,
        "c8 uses G/U/W/W8; c10 reaches W10",
    )
    record(
        checks,
        "source inputs present",
        all(path.is_file() for path in SOURCE_PATHS),
        f"{len(SOURCE_PATHS)} exact inputs",
    )

    elapsed = time.process_time() - started
    rss = max_rss_bytes()
    record(
        checks,
        "producer process-time budget",
        elapsed < CPU_BUDGET_SECONDS,
        f"process_time_seconds={elapsed:.6f}",
    )
    record(
        checks,
        "producer RSS cap",
        rss < RSS_CAP_BYTES,
        f"peak_rss_bytes={rss}",
    )
    if not all(bool(row["passed"]) for row in checks):
        failed = [str(row["name"]) for row in checks if not row["passed"]]
        raise AssertionError("integrated certificate failed: " + ", ".join(failed))

    data = {
        "headline": {
            "tag": "[THEOREM]",
            "statement": (
                "[THEOREM] The eighth direct interlayer coefficient has an exact all-v "
                "finite connected-correlation formula over its eight composition profiles, "
                "and W8 is forced with square coefficient 1/8!. At tenth order there are "
                "exactly sixteen composition profiles; the same derivation forces W10 with "
                "square coefficient 1/10!. More generally C_(2m)^2/(2m)! is unavoidable in "
                "the universal partition-cumulant algebra."
            ),
            "acceptance_result": (
                "[THEOREM] Exact method-limitation obstruction: previously solved connected "
                "objects cannot polynomially close the next even interlayer order without the "
                "new top-rank connected layer datum."
            ),
            "finite_evidence_scope": (
                "[COMPUTATION] Every already available c8 formal coefficient through v^12 "
                "is checked exactly. [UNRESOLVED] No c10 FLM coefficients beyond the all-order "
                "v^0 and v^2 columns are claimed."
            ),
        },
        "c8_closed": c8,
        "c10_structure": c10,
        "resources": {
            "tag": "[COMPUTATION]",
            "process_time_seconds": round(elapsed, 6),
            "peak_rss_bytes": rss,
            "new_flm_boxes_launched": 0,
        },
    }
    return data, checks


def make_artifact(
    data: dict[str, object], checks: list[dict[str, object]]
) -> dict[str, object]:
    return {
        "meta": {
            "experiment": "e223_c10_certificate",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "command": ".venv/bin/python experiments/e223_c10_certificate.py",
            "working_directory": str(Path.cwd()),
            "python_executable": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "arithmetic": (
                "exact integers and fractions.Fraction; hashes are integrity metadata; "
                "process time and RSS are noncertifying resource observations"
            ),
            "budgets": {
                "process_time_seconds": CPU_BUDGET_SECONDS,
                "rss_cap_bytes": RSS_CAP_BYTES,
            },
            "source_sha256": {
                str(path.relative_to(ROOT)): sha256_file(path) for path in SOURCE_PATHS
            },
        },
        "data": data,
        "checks": checks,
    }


def main() -> int:
    data, checks = run_certificate()
    artifact = make_artifact(data, checks)
    ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    ARTIFACT.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {ARTIFACT.relative_to(ROOT)}")
    print(
        json.dumps(
            {
                "checks": len(checks),
                "process_time_seconds": data["resources"]["process_time_seconds"],
                "peak_rss_bytes": data["resources"]["peak_rss_bytes"],
            },
            sort_keys=True,
        )
    )
    print("PASS e223 interlayer c8 closed formula and c10 obstruction")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
