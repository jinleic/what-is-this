#!/usr/bin/env python3
"""Create the externally publishable raw-hash lock for secure replay reports."""

from __future__ import annotations

import argparse
import datetime as _datetime
import hashlib
import json
from pathlib import Path
import sys

EXPECTED_CAMPAIGN_ID = "20260818T212601Z_425f109c15b64a6198785c6cebbbdaab"
EXPECTED_CAMPAIGN_LOCK_SHA256 = "d254a23a4df7cb3a4c6aae2453883324fba2381c1890c4f011d61c7971b84b48"
EXPECTED_VERIFIER_SHA256 = "e67058df1347545c39b6c2d27c53e758d68f4ff338e42225acff51836aa91a9a"
EXPECTED_RUNTIME_SEAL_RAW_SHA256 = "58c49334a92b2d9e8b02b672d9e67da01b772a626b0ec29e88347d1cc7f92a0c"


class LockError(RuntimeError):
    """A fail-closed report-lock error."""


def canonical_bytes(value):
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("ascii")


def object_hash(value):
    reduced = {key: item for key, item in value.items() if key != "report_sha256"}
    return hashlib.sha256(canonical_bytes(reduced)).hexdigest()


def read_json(path, label):
    path = Path(path)
    if not path.is_file() or path.is_symlink():
        raise LockError("%s must be a regular non-symlink file" % label)
    raw = path.read_bytes()
    try:
        value = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise LockError("%s is malformed JSON" % label) from exc
    return raw, value


def create_lock(
    report_directory, launch_manifest_path, live_runtime_attestation_path
):
    report_directory = Path(report_directory)
    launch_raw, launch = read_json(launch_manifest_path, "secure launch manifest")
    if (not isinstance(launch, dict)
            or launch.get("schema") != "uc-secure-independent-arithmetic-launch-v1"
            or launch.get("campaign_id") != EXPECTED_CAMPAIGN_ID
            or launch.get("campaign_lock_sha256") != EXPECTED_CAMPAIGN_LOCK_SHA256
            or launch.get("verifier_sha256") != EXPECTED_VERIFIER_SHA256
            or launch.get("runtime_seal_raw_sha256") != EXPECTED_RUNTIME_SEAL_RAW_SHA256
            or launch.get("resume_used") is not False
            or launch.get("initial_trace_offset") != 0):
        raise LockError("secure launch manifest does not match pinned byte-zero run")
    live_raw, live = read_json(
        live_runtime_attestation_path, "live runtime attestation"
    )
    if (not isinstance(live, dict)
            or live.get("schema") != "uc-live-replay-runtime-attestation-v1"
            or live.get("outcome") != "PASS"
            or live.get("claim_status") != "MACHINE-VERIFIED"
            or live.get("process_count") != 8
            or live.get("runtime_seal_raw_sha256")
            != EXPECTED_RUNTIME_SEAL_RAW_SHA256
            or live.get("launch_manifest_raw_sha256")
            != hashlib.sha256(launch_raw).hexdigest()
            or live.get("report_sha256") != object_hash(live)):
        raise LockError("live runtime attestation does not bind the secure launch")

    reports = []
    for index in range(8):
        path = report_directory / ("slice-%d.json" % index)
        raw, report = read_json(path, "slice %d report" % index)
        if (not isinstance(report, dict)
                or report.get("schema") != "uc-independent-arithmetic-report-v1"
                or report.get("outcome") != "PASS"
                or report.get("claim_status") != "MACHINE-VERIFIED"
                or report.get("campaign_id") != EXPECTED_CAMPAIGN_ID
                or report.get("slice") != index
                or report.get("initial_trace_offset") != 0
                or report.get("resume_used") is not False
                or report.get("arithmetic_complete") is not True
                or report.get("trace_eof") is not True
                or report.get("pending_stack") != 0
                or report.get("finite_event_limit") is not None
                or report.get("report_sha256") != object_hash(report)):
            raise LockError("slice %d is not a sealed byte-zero PASS" % index)
        inputs = report.get("input_sha256", {})
        runtime_seal = report.get("runtime_seal", {})
        if (inputs.get("verifier") != EXPECTED_VERIFIER_SHA256
                or inputs.get("lock") != EXPECTED_CAMPAIGN_LOCK_SHA256
                or runtime_seal.get("raw_sha256") != EXPECTED_RUNTIME_SEAL_RAW_SHA256):
            raise LockError("slice %d input pins differ" % index)
        reports.append({
            "slice": index,
            "path": path.name,
            "size": len(raw),
            "raw_sha256": hashlib.sha256(raw).hexdigest(),
            "report_sha256": report["report_sha256"],
        })

    result = {
        "schema": "uc-independent-report-lock-v1",
        "claim_status": "OPERATOR-RECORDED",
        "created_utc": _datetime.datetime.now(_datetime.timezone.utc).isoformat(),
        "campaign_id": EXPECTED_CAMPAIGN_ID,
        "campaign_lock_sha256": EXPECTED_CAMPAIGN_LOCK_SHA256,
        "verifier_sha256": EXPECTED_VERIFIER_SHA256,
        "runtime_seal_raw_sha256": EXPECTED_RUNTIME_SEAL_RAW_SHA256,
        "launch_manifest": {
            "path": Path(launch_manifest_path).name,
            "size": len(launch_raw),
            "raw_sha256": hashlib.sha256(launch_raw).hexdigest(),
        },
        "live_runtime_attestation": {
            "path": Path(live_runtime_attestation_path).name,
            "size": len(live_raw),
            "raw_sha256": hashlib.sha256(live_raw).hexdigest(),
            "report_sha256": live["report_sha256"],
        },
        "reports": reports,
        "limitations": [
            "This public digest lock is an external publication anchor, not a secret-key signature or hardware attestation.",
            "Acceptance still trusts the local supervisor, operating system, interpreter, and hash implementation not to forge execution evidence."
        ],
    }
    result["report_sha256"] = object_hash(result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description="Lock eight secure UC replay report bytes")
    parser.add_argument("--reports", required=True)
    parser.add_argument("--launch-manifest", required=True)
    parser.add_argument("--live-runtime-attestation", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    try:
        result = create_lock(
            args.reports,
            args.launch_manifest,
            args.live_runtime_attestation,
        )
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("xb") as stream:
            stream.write(canonical_bytes(result) + b"\n")
    except (LockError, OSError, ValueError, TypeError) as exc:
        print("report lock FAIL: %s" % exc, file=sys.stderr)
        return 1
    print(canonical_bytes({
        "claim_status": result["claim_status"],
        "output": str(output),
        "report_sha256": result["report_sha256"],
    }).decode("ascii"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
