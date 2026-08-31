#!/usr/bin/env python3
"""Validate quant exploration artifact shape without third-party dependencies."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


ALLOWED_STATUSES = {"validation-pass", "falsified", "inconclusive", "blocked"}
EXPECTED_RETURNS_HEADER = [
    "date",
    "strategy_return",
    "benchmark_return",
    "position",
    "turnover",
    "cost",
]


def load_json(path: Path) -> Any:
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


def trial_records(payload: Any) -> list[Any] | None:
    if isinstance(payload, list):
        return payload
    if not isinstance(payload, dict):
        return None
    for key in ("trials", "variants", "records"):
        value = payload.get(key)
        if isinstance(value, list):
            return value
        if isinstance(value, dict):
            return list(value.values())
    dev_metrics = payload.get("dev_metrics_at_5bps")
    if isinstance(dev_metrics, dict):
        table = dev_metrics.get("table")
        if isinstance(table, list):
            return table
        if isinstance(table, dict):
            return list(table.values())
    return None


def validate_directory(directory: Path, required_fields: list[str]) -> dict[str, Any]:
    errors: list[str] = []
    warnings: list[str] = []
    strategy_path = directory / "strategy.json"
    trials_path = directory / "trials.json"
    returns_path = directory / "returns.csv"
    run_path = directory / "run.py"

    try:
        strategy = load_json(strategy_path)
    except (OSError, json.JSONDecodeError) as exc:
        return {"directory": str(directory), "errors": [f"strategy.json: {exc}"], "warnings": []}

    if not isinstance(strategy, dict):
        return {"directory": str(directory), "errors": ["strategy.json must be an object"], "warnings": []}

    missing = sorted(field for field in required_fields if field not in strategy)
    if missing:
        errors.append(f"missing strategy fields: {', '.join(missing)}")

    if strategy.get("schema_version") != 1:
        errors.append("schema_version must equal 1")
    if strategy.get("status") not in ALLOWED_STATUSES:
        errors.append(f"invalid status: {strategy.get('status')!r}")
    if not isinstance(strategy.get("holdout_accessed"), bool):
        errors.append("holdout_accessed must be boolean")

    count = strategy.get("trial_count")
    if not isinstance(count, int) or isinstance(count, bool) or not 0 <= count <= 12:
        errors.append("trial_count must be an integer from 0 through 12")

    if not trials_path.is_file():
        errors.append("trials.json is missing")
    else:
        try:
            trial_payload = load_json(trials_path)
            trials = trial_records(trial_payload)
            explicit_zero = (
                count == 0
                and isinstance(trial_payload, dict)
                and trial_payload.get("trial_count") == 0
            )
            if trials is None:
                warnings.append("trials.json is descriptive rather than an enumerable variant list")
            elif not explicit_zero and isinstance(count, int) and len(trials) != count:
                errors.append(f"trial_count={count} but trials.json contains {len(trials)} records")
        except (OSError, json.JSONDecodeError) as exc:
            errors.append(f"trials.json: {exc}")

    status = strategy.get("status")
    returns_applicable = strategy.get("returns_applicable", True)
    if not isinstance(returns_applicable, bool):
        errors.append("returns_applicable must be boolean when present")
        returns_applicable = True
    empirical = status in {"validation-pass", "falsified"} or returns_path.exists()
    if empirical and not run_path.is_file():
        errors.append("empirical result is missing run.py")
    if empirical and returns_applicable and not returns_path.is_file():
        errors.append("empirical result is missing returns.csv")
    if not returns_applicable and returns_path.exists():
        errors.append("returns_applicable=false but returns.csv exists")

    if returns_path.is_file():
        try:
            with returns_path.open(newline="", encoding="utf-8") as stream:
                reader = csv.reader(stream)
                header = next(reader, None)
                row_count = sum(1 for _ in reader)
            if header != EXPECTED_RETURNS_HEADER:
                errors.append(f"returns.csv header must be {','.join(EXPECTED_RETURNS_HEADER)}")
            if row_count == 0:
                errors.append("returns.csv has no observations")
        except OSError as exc:
            errors.append(f"returns.csv: {exc}")
            row_count = 0
    else:
        row_count = 0

    if strategy.get("holdout_accessed") and directory.name != "frozen-cross-asset-trend-audit":
        warnings.append("holdout_accessed=true; candidate is ineligible for historical confirmation")

    return {
        "directory": str(directory),
        "id": strategy.get("id"),
        "status": status,
        "trial_count": count,
        "return_rows": row_count,
        "errors": errors,
        "warnings": warnings,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    contract = load_json(args.root / "contract.json")
    required_fields = contract["required_files"]["strategy.json"]["required_fields"]
    directories = sorted(path.parent for path in args.root.glob("*/strategy.json"))
    results = [validate_directory(directory, required_fields) for directory in directories]
    report = {
        "schema_version": 1,
        "directories_found": len(results),
        "valid_directories": sum(not item["errors"] for item in results),
        "error_count": sum(len(item["errors"]) for item in results),
        "warning_count": sum(len(item["warnings"]) for item in results),
        "results": results,
    }
    rendered = json.dumps(report, indent=2, sort_keys=False) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0 if report["error_count"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
