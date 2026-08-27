#!/usr/bin/env python3
"""Aggregate eight freshly generated clean-room slice reports.

This is a self-consistency and external-pin validator, not an execution
attestation. Repository reports use public SHA-256 values; an independent
recipient must run the slice replays rather than trust copied reports.

Run from ``math`` with ``python -I -B`` after all eight fresh slice runs.
"""

from __future__ import annotations

import argparse
import datetime as _datetime
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from clean_room_replay import (
    CAMPAIGN_DIRECTORY,
    DEFAULT_LOCK,
    VerificationError,
    _write_json_atomic,
    canonical_bytes,
    expected_input_fingerprints,
    file_fingerprint,
    load_and_validate_lock,
    load_json_file,
    object_sha256,
    require_isolated_interpreter,
    validate_environment,
    verifier_fingerprints,
)

SCHEMA = "uc-clean-room-composite-report-v1"
EMPTY_SHA256 = hashlib.sha256(b"").hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def validate_slice(
    report: dict[str, object],
    index: int,
    lock: dict[str, object],
    tools: dict[str, object],
) -> dict[str, object]:
    require(report.get("report_sha256") == object_sha256(report, "report_sha256"),
            f"slice {index} report object hash mismatch")
    require(report.get("schema") == "uc-clean-room-replay-report-v1",
            f"slice {index} report schema mismatch")
    require(report.get("report_type") == "clean_room_replay",
            f"slice {index} report type mismatch")
    require(report.get("campaign_id") == lock["campaign_id"],
            f"slice {index} campaign mismatch")
    require(report.get("slice") == index, f"slice {index} index mismatch")
    require(report.get("claim_status") == "MACHINE-VERIFIED",
            f"slice {index} is not MACHINE-VERIFIED")
    require(report.get("outcome") == "PASS", f"slice {index} did not PASS")
    require(report.get("errors") == [], f"slice {index} contains errors")

    expected_inputs = expected_input_fingerprints(lock, index)
    inputs = report.get("inputs")
    require(type(inputs) is dict, f"slice {index} has no input map")
    for name in ("source_before", "source_after", "stage_before", "stage_after"):
        require(inputs.get(name) == expected_inputs,
                f"slice {index} {name} differs from external pins")

    require(report.get("verifier_inputs_before") == tools,
            f"slice {index} verifier inputs differ from current bytes")
    require(report.get("verifier_inputs_after") == tools,
            f"slice {index} verifier inputs drifted")
    before = report.get("environment_before")
    after = report.get("environment_after")
    require(type(before) is dict and before == after,
            f"slice {index} environment drifted")
    validate_environment(before, lock)

    command = report.get("command")
    require(type(command) is dict, f"slice {index} has no command record")
    require(command.get("termination") == "exit" and command.get("exit_code") == 0,
            f"slice {index} process did not exit 0")
    require(command.get("stderr_utf8") == "" and command.get("stderr_sha256") == EMPTY_SHA256,
            f"slice {index} emitted stderr")
    stdout = command.get("stdout_utf8")
    require(type(stdout) is str,
            f"slice {index} stdout is absent")
    require(hashlib.sha256(stdout.encode("utf-8")).hexdigest() == command.get("stdout_sha256"),
            f"slice {index} stdout hash mismatch")
    lines = stdout.splitlines()
    require(len(lines) == 1, f"slice {index} stdout is not one JSON line")
    try:
        parsed = json.loads(lines[0])
    except json.JSONDecodeError as exc:
        raise VerificationError(f"slice {index} stdout is not JSON: {exc}") from exc

    expected_tallies = lock["slices"][index]["expected_replay_tallies"]
    require(parsed == expected_tallies,
            f"slice {index} stdout tallies differ from external pins")
    replay = report.get("replay")
    require(type(replay) is dict, f"slice {index} replay record is absent")
    require(replay.get("implementation") == "snapshot/cert3_replay.py",
            f"slice {index} replay implementation differs")
    require(replay.get("implementation_sha256") ==
            lock["snapshot"]["cert3_replay.py"]["sha256"],
            f"slice {index} replay implementation hash differs")
    require(replay.get("shares_frozen_arithmetic") is True,
            f"slice {index} arithmetic-sharing disclosure is absent")
    require(replay.get("tallies") == expected_tallies,
            f"slice {index} replay tallies differ")

    return {
        "attempt_id": report.get("attempt_id"),
        "elapsed_seconds": command.get("elapsed_seconds"),
        "report_sha256": report["report_sha256"],
        "slice": index,
        "stdout_sha256": command["stdout_sha256"],
        "tallies": expected_tallies,
    }


def run(args: argparse.Namespace) -> int:
    require_isolated_interpreter()
    lock_path = args.lock.expanduser().resolve()
    lock, lock_fingerprint = load_and_validate_lock(lock_path)
    reports_root = args.reports.expanduser().resolve()
    tools = verifier_fingerprints(lock_path)
    slices = []
    for index in range(8):
        path = reports_root / f"slice-{index}" / "latest.json"
        report = load_json_file(path, maximum_size=8 * 1024 * 1024)
        slices.append(validate_slice(report, index, lock, tools))

    total_processed = sum(item["tallies"]["processed"] for item in slices)
    total_elapsed = sum(float(item["elapsed_seconds"]) for item in slices)
    require(total_processed == 488_465_854,
            "composite processed tally is not 488,465,854")
    require(all(item["tallies"]["residual"] == 0 for item in slices),
            "composite contains residual events")

    result: dict[str, object] = {
        "campaign_directory": CAMPAIGN_DIRECTORY,
        "campaign_id": lock["campaign_id"],
        "claim_status": "MACHINE-VERIFIED",
        "finished_utc": _datetime.datetime.now(_datetime.timezone.utc).isoformat(),
        "limitations": [
            "Every arithmetic replay uses the same frozen rule implementation as the producer.",
            "This aggregate validates public report bytes; it is not an authenticated execution attestation.",
        ],
        "lock": lock_fingerprint,
        "outcome": "PASS",
        "report_type": "clean_room_composite",
        "schema": SCHEMA,
        "slices": slices,
        "total_elapsed_seconds": total_elapsed,
        "total_processed": total_processed,
        "verifier": file_fingerprint(Path(__file__).resolve()),
        "verifier_inputs": tools,
    }
    result["report_sha256"] = object_sha256(result, "report_sha256")
    output = args.output.expanduser().resolve()
    output.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    _write_json_atomic(output, result)
    print(canonical_bytes({
        "claim_status": "MACHINE-VERIFIED",
        "outcome": "PASS",
        "report": str(output),
        "report_sha256": result["report_sha256"],
        "total_processed": total_processed,
    }).decode("ascii"))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Aggregate eight fresh clean-room replay reports")
    parser.add_argument("--reports", type=Path, required=True,
                        help="directory containing slice-0 through slice-7")
    parser.add_argument("--output", type=Path, required=True,
                        help="composite JSON output path")
    parser.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    return parser


def main() -> int:
    try:
        return run(build_parser().parse_args())
    except (VerificationError, OSError) as exc:
        print(canonical_bytes({
            "claim_status": "FAILED",
            "error": str(exc),
            "outcome": "FAIL_CLOSED",
        }).decode("ascii"), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
