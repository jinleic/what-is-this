#!/usr/bin/env python3
"""Exercise the certificate pipeline with real Kissat, XZ, and drat-trim."""
from __future__ import annotations

import json
import re
import shutil
import sys
import time
from pathlib import Path

if str(Path(__file__).resolve().parent) not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parent))

import certify_unsat

WORKSPACE = Path(__file__).resolve().parents[3]
DECIDE = Path(__file__).resolve().parent
SMOKE = DECIDE / "certifier-real-smoke"
CASES = (
    ("proof_nonup", DECIDE / "xz_pipeline_nonup_smoke.cnf"),
    ("input_up", DECIDE / "xz_pipeline_input_up_smoke.cnf"),
)


def write_once(path: Path, payload: bytes) -> None:
    if path.exists():
        if path.is_symlink() or path.read_bytes() != payload:
            raise certify_unsat.Refused(f"preserved smoke artifact differs: {path}")
        return
    certify_unsat.atomic_write(path, payload)


def build_decision() -> dict:
    source = json.loads((DECIDE / "manifest.json").read_bytes())
    instances = []
    for case, path in CASES:
        record = certify_unsat.file_record(path, WORKSPACE)
        variables, clauses, _maximum = certify_unsat.inspect_dimacs(path)
        instances.append(
            {
                **record,
                "clauses": clauses,
                "cube": case,
                "variables": variables,
            }
        )
    return {
        "case_count": len(CASES),
        "case_order": [case for case, _path in CASES],
        "instances": instances,
        "schema": certify_unsat.DECISION_SCHEMA,
        "solver": source["solver"],
    }


def validate_record(record: dict) -> Path:
    if not isinstance(record, dict) or set(record) != {"bytes", "path", "sha256"}:
        raise certify_unsat.Refused("malformed smoke artifact record")
    path = Path(record["path"])
    if not path.is_absolute():
        path = WORKSPACE / path
    actual = certify_unsat.file_record(path, WORKSPACE)
    if actual != record:
        raise certify_unsat.Refused(f"smoke artifact changed: {path}")
    return path


def corrupt_control(certifier: certify_unsat.Certifier, proof_row: dict) -> dict:
    assert certifier.certificate_sha is not None
    control = SMOKE / "negative-controls" / certifier.certificate_sha
    control.mkdir(parents=True, exist_ok=True)
    receipt_path = control / "result.json"
    if receipt_path.exists():
        receipt = json.loads(receipt_path.read_bytes())
        if (
            certify_unsat.canonical_json(receipt) != receipt_path.read_bytes()
            or receipt.get("schema") != "kobon-certifier-corrupt-proof-control/1"
            or receipt.get("certificate_manifest_sha256") != certifier.certificate_sha
            or receipt.get("outcome") != "XZ_INTEGRITY_FAILED"
        ):
            raise certify_unsat.Refused("invalid preserved corrupt-proof control")
        validate_record(receipt["corrupted_proof"])
        validate_record(receipt["xz_integrity"]["log"])
        return receipt

    proof = certifier._artifact_path(proof_row["proof"])
    corrupted = control / "corrupted.drat.xz"
    data = bytearray(proof.read_bytes())
    if not data:
        raise certify_unsat.Refused("cannot corrupt an empty proof")
    data[-1] ^= 1
    write_once(corrupted, bytes(data))
    log_paths = {
        "xz_integrity": control / "xz-integrity.log",
        "input_up": control / "input-up.log",
        "xz_stream": control / "xz-stream.log",
        "checker": control / "checker.log",
    }
    if any(path.exists() for path in log_paths.values()):
        raise certify_unsat.Refused("partial corrupt-proof control is preserved")

    certifier.driver_lock.acquire()
    certifier._install_signals()
    try:
        integrity, input_up, proof_check, verdict, outcome = certifier._verify_proof(
            "proof_nonup", corrupted, log_paths
        )
    finally:
        certifier._terminate_all()
        certifier._restore_signals()
        certifier.driver_lock.release()
    if (
        outcome != "XZ_INTEGRITY_FAILED"
        or integrity is None
        or integrity["exit_code"] == 0
        or input_up is not None
        or proof_check is not None
        or verdict is not None
    ):
        raise certify_unsat.Refused(f"corrupt proof was not rejected safely: {outcome}")
    receipt = {
        "certificate_manifest_sha256": certifier.certificate_sha,
        "corrupted_proof": certify_unsat.file_record(corrupted, WORKSPACE),
        "outcome": outcome,
        "schema": "kobon-certifier-corrupt-proof-control/1",
        "source_proof": proof_row["proof"],
        "xz_integrity": integrity,
    }
    certify_unsat.atomic_write(receipt_path, certify_unsat.canonical_json(receipt))
    return receipt


def main() -> int:
    SMOKE.mkdir(parents=True, exist_ok=True)
    decision = build_decision()
    decision_payload = (json.dumps(decision, indent=2, sort_keys=True) + "\n").encode()
    decision_path = SMOKE / "manifest.json"
    write_once(decision_path, decision_payload)
    decision_sha = certify_unsat.sha256_bytes(decision_payload)

    required_tools = {name: shutil.which(name) for name in ("xz", "pgrep", "ps")}
    if any(path is None for path in required_tools.values()):
        raise certify_unsat.Refused("XZ, pgrep, and ps are required for the real smoke")
    config = certify_unsat.Config(
        workspace=WORKSPACE,
        here=SMOKE,
        manifest_path=decision_path,
        checker_path=WORKSPACE / "scratch/kobon-audit/tools/drat-trim/drat-trim",
        driver_lock_path=SMOKE / ".certify-unsat.lock",
        host_lock_path=WORKSPACE / "scratch/.host-heavy-job.lock",
        expected_decision_sha256=decision_sha,
        expected_cases=tuple(case for case, _path in CASES),
        xz_path=Path(required_tools["xz"]),
        pgrep_path=Path(required_tools["pgrep"]),
        ps_path=Path(required_tools["ps"]),
        disk_floor_bytes=certify_unsat.default_config().disk_floor_bytes,
        memory_limit_bytes=8 << 30,
        disk_poll_seconds=0.1,
        process_poll_seconds=0.05,
        solver_timeout_seconds=60.0,
        checker_timeout_seconds=60.0,
        resource_wait_seconds=28800.0,
        kill_grace_seconds=3.0,
    )
    certifier = certify_unsat.Certifier(config, manifest_builder=build_decision)
    if certifier.run(
        [case for case, _path in CASES], retry_nonterminal=True
    ) != certify_unsat.EXIT_OK:
        raise certify_unsat.Refused("real smoke did not produce two terminal outcomes")

    rows = {row["case"]: row for row in certifier.ledger["rows"] if row["terminal"]}
    proof_row = rows.get("proof_nonup")
    input_row = rows.get("input_up")
    if proof_row is None or proof_row["outcome"] != "CERTIFIED_UNSAT":
        raise certify_unsat.Refused("non-unit formula lacked a proof-consuming verdict")
    if proof_row["evidence_class"] != "PROOF_CONSUMING_BINARY_DRAT_XZ":
        raise certify_unsat.Refused("non-unit formula has the wrong evidence class")
    if input_row is None or input_row["outcome"] != "VERIFIED_INPUT_UP_UNSAT":
        raise certify_unsat.Refused("input-UP formula was not classified separately")
    if input_row["evidence_class"] != "INPUT_UNIT_PROPAGATION_UNSAT":
        raise certify_unsat.Refused("input-UP formula was incorrectly proof-certified")

    checker_log = certifier._artifact_path(proof_row["proof_check"]["checker_log"])
    checker_data = checker_log.read_bytes()
    normalized_checker = checker_data.replace(b"\r", b"\n")
    match = re.search(
        rb"^c ([1-9][0-9]*) of [0-9]+ lemmas in core",
        normalized_checker,
        re.MULTILINE,
    )
    if b"s VERIFIED" not in checker_data or match is None:
        raise certify_unsat.Refused("real checker did not consume nonzero proof lemmas")
    negative = corrupt_control(certifier, proof_row)

    assert certifier.namespace is not None and certifier.ledger_path is not None
    receipt = {
        "certificate_manifest": certify_unsat.file_record(
            certifier.namespace / "certificate-manifest.json", WORKSPACE
        ),
        "certificate_manifest_sha256": certifier.certificate_sha,
        "checked_utc": certify_unsat.utc(time.time_ns()),
        "checks": {
            "corrupt_proof_outcome": negative["outcome"],
            "input_up_evidence_class": input_row["evidence_class"],
            "input_up_outcome": input_row["outcome"],
            "nonunit_core_lemmas": int(match.group(1)),
            "nonunit_evidence_class": proof_row["evidence_class"],
            "nonunit_outcome": proof_row["outcome"],
        },
        "decision_manifest": certify_unsat.file_record(decision_path, WORKSPACE),
        "negative_control": certify_unsat.file_record(
            SMOKE / "negative-controls" / certifier.certificate_sha / "result.json",
            WORKSPACE,
        ),
        "schema": "kobon-certifier-real-tool-smoke/1",
        "status": "PASS",
        "terminal_row_sha256": {
            case: certify_unsat.sha256_bytes(certify_unsat.canonical_json(row))
            for case, row in sorted(rows.items())
        },
    }
    receipt_path = SMOKE / f"result.{certifier.certificate_sha}.json"
    payload = certify_unsat.canonical_json(receipt)
    if receipt_path.exists():
        preserved = json.loads(receipt_path.read_bytes())
        if certify_unsat.canonical_json(preserved) != receipt_path.read_bytes():
            raise certify_unsat.Refused("preserved real-smoke receipt is not canonical")
        preserved_without_time = {key: value for key, value in preserved.items() if key != "checked_utc"}
        receipt_without_time = {key: value for key, value in receipt.items() if key != "checked_utc"}
        if preserved_without_time != receipt_without_time:
            raise certify_unsat.Refused("preserved real-smoke receipt differs")
        receipt = preserved
        payload = certify_unsat.canonical_json(receipt)
    else:
        certify_unsat.atomic_write(receipt_path, payload)
    print(payload.decode().rstrip())
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except certify_unsat.Refused as error:
        print(f"REFUSED: {error}", file=sys.stderr)
        raise SystemExit(certify_unsat.EXIT_REFUSED)
