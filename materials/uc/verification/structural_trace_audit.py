#!/usr/bin/env python3
"""Audit cert3-trace-v1 structure without importing proof arithmetic.

The auditor is standard-library-only. It independently transcribes the one-byte
opcode format and checks hashes, exact dyadic roots, DFS topology/termination,
and result tallies. It does not import or evaluate the generator, collector,
replayer, interval predicates, flint, or mpmath; consequently this is not a
second interval proof.
"""

from __future__ import annotations

import argparse
import datetime as _datetime
import hashlib
import json
import os
from pathlib import Path
import stat
import sys
import time
import uuid
from fractions import Fraction


REPORT_SCHEMA = "uc-structural-trace-audit-report-v1"
LOCK_RAW_SHA256 = "d254a23a4df7cb3a4c6aae2453883324fba2381c1890c4f011d61c7971b84b48"
CODE_SHA256 = "2f23a58ebdb8b14275284e2ffefca7862373f7baf5d354a06e607cff8d78ce05"
LAUNCH_OBJECT_SHA256 = "6414affcaa54bd5e37e9f7bd351c6ed4cf675167b10b55d86cce1436ebdf65cc"
LAUNCH_RAW_SHA256 = "63bb6fbe1cfbcb6b7d4c50d7b66882255e5e778ce604c92ecd1abdfa1ea8337e"
COLLECTOR_SHA256 = "95d09322a7821f41850ef1ce77345e5e93eb4fa570cc0e7f58f817716cf681ec"
CAMPAIGN_ID = "20260818T212601Z_425f109c15b64a6198785c6cebbbdaab"
CAMPAIGN_DIRECTORY = (
    "cert3_20260818T212601Z_425f109c15b64a6198785c6cebbbdaab_"
    "2f23a58ebdb8"
)
HERE = Path(__file__).resolve().parent
DEFAULT_LOCK = HERE / "campaign-lock.json"
DEFAULT_CAMPAIGN = HERE.parent / "campaigns" / CAMPAIGN_DIRECTORY
SCHEMA_PATH = HERE / "report.schema.json"

# Independently transcribed cert3-trace-v1 contract.
TERMINALS = {
    0: "mean_infeasible",
    1: "corner_infeasible",
    2: "corner",
    3: "ratio",
    4: "center",
    5: "center_mixed",
    6: "center_mixed_swap",
    7: "center_w",
    8: "face",
    9: "residual",
}
SPLITS = {16 + coordinate: f"split_{coordinate}" for coordinate in range(5)}
VALID_OPCODES = frozenset(TERMINALS) | frozenset(SPLITS)
REPLAY_TALLIES = (
    "processed",
    "infeasible",
    "corner",
    "ratio",
    "center",
    "center_mixed",
    "center_mixed_swap",
    "center_w",
    "face",
    "residual",
    "split",
)
SNAPSHOT_NAMES = frozenset(
    {
        "arbcore.py",
        "bound_kkt.py",
        "bridge_uc.py",
        "cert2.py",
        "cert3.py",
        "cert3_collect.py",
        "cert3_par.py",
        "cert3_replay.py",
        "decomposition.py",
        "diag_exhaust.py",
        "entropy.py",
        "lemma_rh_proof.py",
        "margin_lemma.py",
        "reduction.py",
        "thmB3_proof.py",
    }
)


class AuditError(RuntimeError):
    pass


def now_utc() -> str:
    return _datetime.datetime.now(_datetime.timezone.utc).isoformat()


def canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def object_hash(value: dict[str, object], omitted: str) -> str:
    return hashlib.sha256(
        canonical_bytes({key: item for key, item in value.items() if key != omitted})
    ).hexdigest()


def no_duplicate_keys(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise AuditError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def require_isolation() -> None:
    if not (
        sys.flags.isolated == 1
        and sys.flags.no_user_site == 1
        and sys.flags.ignore_environment == 1
        and sys.flags.dont_write_bytecode == 1
    ):
        raise AuditError("invoke with an isolated no-bytecode interpreter: python -I -B")


def open_regular(path: Path) -> tuple[int, os.stat_result]:
    flags = os.O_RDONLY
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        raise AuditError(f"cannot open required regular file {path}: {exc}") from exc
    before = os.fstat(descriptor)
    if not stat.S_ISREG(before.st_mode) or stat.S_ISLNK(path.lstat().st_mode):
        os.close(descriptor)
        raise AuditError(f"required input is not a non-symlink regular file: {path}")
    return descriptor, before


def unchanged(before: os.stat_result, after: os.stat_result) -> bool:
    return (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    ) == (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    )


def fingerprint(path: Path) -> dict[str, object]:
    descriptor, before = open_regular(path)
    digest = hashlib.sha256()
    size = 0
    try:
        while True:
            chunk = os.read(descriptor, 4 * 1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
            size += len(chunk)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if not unchanged(before, after) or size != before.st_size:
        raise AuditError(f"input changed while hashing: {path}")
    return {"sha256": digest.hexdigest(), "size": size}


def load_json(path: Path) -> dict[str, object]:
    descriptor, before = open_regular(path)
    try:
        if before.st_size > 2 * 1024 * 1024:
            raise AuditError(f"JSON input is unexpectedly large: {path}")
        payload = bytearray()
        while len(payload) < before.st_size:
            chunk = os.read(descriptor, before.st_size - len(payload))
            if not chunk:
                raise AuditError(f"JSON input was truncated: {path}")
            payload.extend(chunk)
        if os.read(descriptor, 1):
            raise AuditError(f"JSON input grew while reading: {path}")
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if not unchanged(before, after):
        raise AuditError(f"JSON input changed while reading: {path}")
    try:
        value = json.loads(
            payload.decode("utf-8"), object_pairs_hook=no_duplicate_keys
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise AuditError(f"cannot parse JSON input {path}: {exc}") from exc
    if type(value) is not dict:
        raise AuditError(f"JSON input is not an object: {path}")
    return value


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def is_dyadic(value: Fraction) -> bool:
    denominator = value.denominator
    return denominator > 0 and denominator & (denominator - 1) == 0

def load_pinned_lock(path: Path) -> tuple[dict[str, object], dict[str, object]]:
    """Hash and parse one stable no-follow descriptor."""
    descriptor, before = open_regular(path)
    try:
        if before.st_size != 12590:
            raise AuditError("external campaign pin manifest has the wrong size")
        payload = bytearray()
        while len(payload) < before.st_size:
            chunk = os.read(descriptor, before.st_size - len(payload))
            if not chunk:
                raise AuditError("external campaign pin manifest was truncated")
            payload.extend(chunk)
        if os.read(descriptor, 1):
            raise AuditError("external campaign pin manifest grew while reading")
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if not unchanged(before, after):
        raise AuditError("external campaign pin manifest changed while reading")
    lock_fingerprint = {
        "sha256": hashlib.sha256(payload).hexdigest(),
        "size": len(payload),
    }
    if lock_fingerprint != {"sha256": LOCK_RAW_SHA256, "size": 12590}:
        raise AuditError("external campaign pin manifest differs from verifier-pinned bytes")
    try:
        value = json.loads(
            payload.decode("utf-8"), object_pairs_hook=no_duplicate_keys
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise AuditError(f"cannot parse external campaign pin manifest: {exc}") from exc
    if type(value) is not dict:
        raise AuditError("external campaign pin manifest is not an object")
    return value, lock_fingerprint



def load_lock(path: Path) -> tuple[dict[str, object], dict[str, object]]:
    lock, lock_fingerprint = load_pinned_lock(path)
    if lock.get("schema") != "uc-campaign-lock-v1":
        raise AuditError("unsupported campaign lock schema")
    if lock.get("campaign_id") != CAMPAIGN_ID or lock.get("code_sha256") != CODE_SHA256:
        raise AuditError("campaign identity/code digest differs from external pins")
    if lock.get("campaign_directory_name") != CAMPAIGN_DIRECTORY:
        raise AuditError("campaign lock directory identity differs from external pin")
    launch_pin = lock.get("launch")
    if launch_pin != {
        "object_sha256": LAUNCH_OBJECT_SHA256,
        "path": "launch.json",
        "raw_sha256": LAUNCH_RAW_SHA256,
        "size": 4339,
    }:
        raise AuditError("launch pin differs from externally captured identity")
    snapshot = lock.get("snapshot")
    if type(snapshot) is not dict or set(snapshot) != SNAPSHOT_NAMES:
        raise AuditError("snapshot pin inventory is not the exact 15-file set")
    if snapshot["cert3_collect.py"].get("sha256") != COLLECTOR_SHA256:
        raise AuditError("collector digest differs from external pin")
    slices = lock.get("slices")
    if type(slices) is not list or len(slices) != 8:
        raise AuditError("external pin manifest does not contain eight slices")
    for index, item in enumerate(slices):
        if type(item) is not dict or item.get("slice") != index:
            raise AuditError("pinned slices are not ordered exactly 0..7")
        expected_lo = fraction_text(Fraction(1, 2) + Fraction(index, 16))
        expected_hi = fraction_text(Fraction(1, 2) + Fraction(index + 1, 16))
        if item.get("w_lo") != expected_lo or item.get("w_hi") != expected_hi:
            raise AuditError(f"pinned slice {index} root differs from exact partition")
        result = item.get("result")
        trace = item.get("trace")
        tallies = item.get("expected_replay_tallies")
        if type(result) is not dict or type(trace) is not dict or type(tallies) is not dict:
            raise AuditError(f"pinned slice {index} is malformed")
        run_id = item.get("run_id")
        if result.get("path") != f"result_slice{index}_{run_id}.json":
            raise AuditError(f"pinned slice {index} result filename is malformed")
        if trace.get("path") != f"trace_slice{index}_{run_id}.bin":
            raise AuditError(f"pinned slice {index} trace filename is malformed")
        if set(tallies) != set(REPLAY_TALLIES):
            raise AuditError(f"pinned slice {index} structural tallies are incomplete")
        if any(type(value) is not int or value < 0 for value in tallies.values()):
            raise AuditError(f"pinned slice {index} structural tally is invalid")
        if item.get("expected_completion_tallies") != {
            "budget_boxes": 0,
            "budget_time": 0,
            "stack": 0,
        }:
            raise AuditError(f"pinned slice {index} completion state is not empty")
    return lock, lock_fingerprint


def validate_inventory(campaign: Path, lock: dict[str, object]) -> list[str]:
    info = campaign.lstat()
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
        raise AuditError("campaign path must be a real directory, not a symlink")
    expected_results = {item["result"]["path"] for item in lock["slices"]}
    expected_traces = {item["trace"]["path"] for item in lock["slices"]}
    entries = list(campaign.iterdir())
    actual_results = {
        entry.name for entry in entries
        if entry.name.startswith("result_") and entry.name.endswith(".json")
    }
    actual_traces = {
        entry.name for entry in entries
        if entry.name.startswith("trace_") and entry.name.endswith(".bin")
    }
    if actual_results != expected_results or actual_traces != expected_traces:
        raise AuditError("result/trace inventory differs from exact external pins")
    accepted = {"launch.json", "snapshot"} | expected_results | expected_traces
    return sorted(entry.name for entry in entries if entry.name not in accepted)


def validate_launch(campaign: Path, lock: dict[str, object]) -> dict[str, object]:
    observed = fingerprint(campaign / "launch.json")
    if observed != {"sha256": LAUNCH_RAW_SHA256, "size": 4339}:
        raise AuditError("raw launch file differs from external pin")
    launch = load_json(campaign / "launch.json")
    if launch.get("launch_sha256") != object_hash(launch, "launch_sha256"):
        raise AuditError("launch canonical object digest is invalid")
    if (
        launch.get("schema") != "cert3-campaign-v2"
        or launch.get("launch_sha256") != LAUNCH_OBJECT_SHA256
        or launch.get("campaign_id") != CAMPAIGN_ID
        or launch.get("code_sha256") != CODE_SHA256
        or launch.get("nslices") != 8
    ):
        raise AuditError("launch identity differs from external Campaign I pins")
    runs = launch.get("runs")
    if type(runs) is not list or len(runs) != 8:
        raise AuditError("launch does not contain exactly eight runs")
    for index, (run, pinned) in enumerate(zip(runs, lock["slices"])):
        if type(run) is not dict or run.get("slice") != index:
            raise AuditError("launch runs are not ordered exactly 0..7")
        for key, value in {
            "result_filename": pinned["result"]["path"],
            "run_id": pinned["run_id"],
            "trace_filename": pinned["trace"]["path"],
            "w_hi": pinned["w_hi"],
            "w_lo": pinned["w_lo"],
        }.items():
            if run.get(key) != value:
                raise AuditError(f"launch slice {index} field {key} differs from pin")
        root = ((Fraction(0), Fraction(1)),) * 4 + (
            (Fraction(run["w_lo"]), Fraction(run["w_hi"])),
        )
        if any(not (is_dyadic(lo) and is_dyadic(hi) and lo < hi) for lo, hi in root):
            raise AuditError(f"launch slice {index} root is not exact dyadic")
    return launch


def validate_result(
    campaign: Path,
    launch: dict[str, object],
    pinned: dict[str, object],
) -> tuple[dict[str, object], dict[str, object]]:
    path = campaign / pinned["result"]["path"]
    observed = fingerprint(path)
    expected = {
        "sha256": pinned["result"]["raw_sha256"],
        "size": pinned["result"]["size"],
    }
    if observed != expected:
        raise AuditError(f"slice {pinned['slice']} result raw file differs from external pin")
    result = load_json(path)
    if result.get("record_sha256") != object_hash(result, "record_sha256"):
        raise AuditError(f"slice {pinned['slice']} result object digest is invalid")
    if result.get("record_sha256") != pinned["result"]["record_sha256"]:
        raise AuditError(f"slice {pinned['slice']} result object differs from external pin")
    worker = {
        key: value for key, value in result.items()
        if key not in ("record_sha256", "result_schema", "supervisor_committed", "worker_exit_code")
    }
    if worker.get("worker_record_sha256") != object_hash(worker, "worker_record_sha256"):
        raise AuditError(f"slice {pinned['slice']} embedded worker digest is invalid")
    index = pinned["slice"]
    run = launch["runs"][index]
    for key, value in {
        "campaign_id": CAMPAIGN_ID,
        "code_sha256": CODE_SHA256,
        "launch_sha256": LAUNCH_OBJECT_SHA256,
        "result_filename": run["result_filename"],
        "run_id": run["run_id"],
        "slice": index,
        "trace_filename": run["trace_filename"],
        "w_hi": run["w_hi"],
        "w_lo": run["w_lo"],
    }.items():
        if result.get(key) != value:
            raise AuditError(f"slice {index} result field {key} differs from launch/pin")
    for key, value in {
        "complete": True,
        "result_schema": "cert3-result-v2",
        "returned_normally": True,
        "supervisor_committed": True,
        "trace_schema": "cert3-trace-v1",
        "verdict": "COMPLETE",
        "worker_exit_code": 0,
        "worker_schema": "cert3-worker-v2",
    }.items():
        if result.get(key) != value:
            raise AuditError(f"slice {index} result field {key} is not {value!r}")
    if result.get("trace_sha256") != pinned["trace"]["sha256"]:
        raise AuditError(f"slice {index} result trace digest differs from external pin")
    if result.get("trace_size") != pinned["trace"]["size"]:
        raise AuditError(f"slice {index} result trace size differs from external pin")
    tallies = result.get("tallies")
    if type(tallies) is not dict:
        raise AuditError(f"slice {index} tallies are absent")
    if any(type(value) is not int for key, value in tallies.items() if key != "residual_geometry"):
        raise AuditError(f"slice {index} has a non-integer result tally")
    if tallies.get("residual_geometry") != {}:
        raise AuditError(f"slice {index} residual geometry is nonempty")
    for key, value in pinned["expected_completion_tallies"].items():
        if tallies.get(key) != value:
            raise AuditError(f"slice {index} completion tally {key} is nonzero")
    return result, observed


def scan_trace(path: Path, chunk_size: int) -> dict[str, object]:
    descriptor, before = open_regular(path)
    digest = hashlib.sha256()
    counts = [0] * 21
    pending = 1
    maximum_pending = 1
    processed = 0
    try:
        while True:
            chunk = os.read(descriptor, chunk_size)
            if not chunk:
                break
            digest.update(chunk)
            for opcode in chunk:
                event_number = processed + 1
                if pending == 0:
                    raise AuditError(f"{path.name}: event {event_number} follows DFS termination")
                if opcode <= 9:
                    pending -= 1
                elif 16 <= opcode <= 20:
                    pending += 1
                    maximum_pending = max(maximum_pending, pending)
                else:
                    raise AuditError(f"{path.name}: event {event_number} has unknown opcode {opcode}")
                counts[opcode] += 1
                processed += 1
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if not unchanged(before, after) or processed != before.st_size:
        raise AuditError(f"trace changed while scanning: {path}")
    if pending != 0:
        raise AuditError(f"{path.name}: EOF left {pending} DFS nodes pending")
    tallies = {
        "processed": processed,
        "infeasible": counts[0] + counts[1],
        "corner": counts[2],
        "ratio": counts[3],
        "center": counts[4],
        "center_mixed": counts[5],
        "center_mixed_swap": counts[6],
        "center_w": counts[7],
        "face": counts[8],
        "residual": counts[9],
        "split": sum(counts[16:21]),
    }
    terminal_count = sum(counts[:10])
    if terminal_count != tallies["split"] + 1:
        raise AuditError(f"{path.name}: full-binary-tree leaf/split invariant failed")
    return {
        "fingerprint": {"sha256": digest.hexdigest(), "size": processed},
        "maximum_pending_nodes": maximum_pending,
        "opcode_counts": {str(code): counts[code] for code in sorted(VALID_OPCODES)},
        "split_by_coordinate": {str(i): counts[16 + i] for i in range(5)},
        "tallies": tallies,
        "terminal_count": terminal_count,
    }


def write_json(path: Path, value: dict[str, object]) -> None:
    path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    if path.parent.is_symlink():
        raise AuditError("output parent must not be a symlink")
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    payload = json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True).encode("ascii") + b"\n"
    descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    try:
        view = memoryview(payload)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise AuditError(f"short report write: {path}")
            view = view[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    os.replace(temporary, path)


def apply_controls(nice_increment: int, cpu_limit: int) -> None:
    if not 0 <= nice_increment <= 19 or cpu_limit < 0:
        raise AuditError("nice must be 0..19 and CPU limit must be nonnegative")
    if cpu_limit:
        import resource
        resource.setrlimit(resource.RLIMIT_CPU, (cpu_limit, cpu_limit + 60))
    if nice_increment:
        os.nice(nice_increment)


def run(args: argparse.Namespace) -> int:
    require_isolation()
    apply_controls(args.nice, args.cpu_limit_seconds)
    if not 4096 <= args.chunk_size <= 64 * 1024 * 1024:
        raise AuditError("chunk size must be 4096..67108864")
    lock_path = args.lock.expanduser().resolve()
    lock, lock_fingerprint = load_lock(lock_path)
    verifier_before = {
        "campaign_lock": lock_fingerprint,
        "report.schema.json": fingerprint(SCHEMA_PATH),
        "structural_trace_audit.py": fingerprint(Path(__file__).resolve()),
    }
    campaign = args.campaign.expanduser().resolve()
    ignored = validate_inventory(campaign, lock)
    launch = validate_launch(campaign, lock)
    selected = list(range(8)) if not args.slice else sorted(set(args.slice))
    if any(index not in range(8) for index in selected):
        raise AuditError("each --slice must be in 0..7")
    output = None if args.output is None else args.output.expanduser().resolve()
    if output is not None and (output == campaign or campaign in output.parents):
        raise AuditError("output must not be inside the frozen campaign")

    started_utc = now_utc()
    started = time.monotonic()
    before: dict[str, dict[str, object]] = {
        "launch.json": {"sha256": LAUNCH_RAW_SHA256, "size": 4339}
    }
    slice_reports = []
    for index in selected:
        pinned = lock["slices"][index]
        result, result_hash = validate_result(campaign, launch, pinned)
        before[pinned["result"]["path"]] = result_hash
        scanned = scan_trace(campaign / pinned["trace"]["path"], args.chunk_size)
        expected_trace = {
            "sha256": pinned["trace"]["sha256"],
            "size": pinned["trace"]["size"],
        }
        if scanned["fingerprint"] != expected_trace:
            raise AuditError(f"slice {index} trace differs from external pin")
        before[pinned["trace"]["path"]] = scanned["fingerprint"]
        if scanned["tallies"] != pinned["expected_replay_tallies"]:
            raise AuditError(f"slice {index} trace tallies differ from external pins")
        result_tallies = result["tallies"]
        for name, value in scanned["tallies"].items():
            if result_tallies.get(name) != value:
                raise AuditError(f"slice {index} trace/result tally {name} differs")
        if result_tallies.get("processed") != result.get("trace_size"):
            raise AuditError(f"slice {index} processed tally differs from trace size")
        unchecked = sorted(
            set(result_tallies)
            - set(REPLAY_TALLIES)
            - {"budget_boxes", "budget_time", "residual_geometry", "stack"}
        )
        slice_reports.append(
            {
                "claim_status": "MACHINE-VERIFIED",
                "dyadic_midpoint_closure": True,
                "maximum_pending_nodes": scanned["maximum_pending_nodes"],
                "opcode_counts": scanned["opcode_counts"],
                "result_tallies_match": True,
                "root_exact": [["0", "1"]] * 4 + [[pinned["w_lo"], pinned["w_hi"]]],
                "run_id": pinned["run_id"],
                "slice": index,
                "split_by_coordinate": scanned["split_by_coordinate"],
                "structural_tallies": scanned["tallies"],
                "terminal_count": scanned["terminal_count"],
                "topology": {
                    "final_pending_nodes": 0,
                    "initial_pending_nodes": 1,
                    "leaf_count_equals_split_count_plus_one": True,
                    "no_event_after_termination": True,
                },
                "trace": scanned["fingerprint"],
                "unchecked_result_tallies": unchecked,
            }
        )

    before = dict(sorted(before.items()))
    after = {name: fingerprint(campaign / name) for name in before}
    if before != after:
        raise AuditError("source inputs changed between scan and final hash")
    verifier_after = {
        "campaign_lock": fingerprint(lock_path),
        "report.schema.json": fingerprint(SCHEMA_PATH),
        "structural_trace_audit.py": fingerprint(Path(__file__).resolve()),
    }
    if verifier_before != verifier_after:
        raise AuditError("verifier, schema, or external pin bytes changed during audit")
    report: dict[str, object] = {
        "$schema": str(SCHEMA_PATH),
        "arithmetic_rules_checked": False,
        "campaign_id": CAMPAIGN_ID,
        "claim_status": "MACHINE-VERIFIED",
        "finished_utc": now_utc(),
        "ignored_campaign_non_inputs": ignored,
        "imports_frozen_arithmetic": False,
        "inputs": {"source_after": after, "source_before": before},
        "limitations": [
            "This audit checks only bytes, hashes, valid opcodes, exact-root metadata, DFS topology/termination, and equality with committed structural tallies.",
            "It does not evaluate any interval discharge predicate and is not a second interval proof.",
            "Auxiliary producer tallies listed in unchecked_result_tallies are hash-locked but are not encoded by cert3-trace-v1.",
        ],
        "opcode_contract": {
            "split": {str(key): value for key, value in SPLITS.items()},
            "terminal": {str(key): value for key, value in TERMINALS.items()},
        },
        "outcome": "PASS",
        "report_type": "structural_trace_audit",
        "resource_controls": {
            "chunk_size": args.chunk_size,
            "cpu_limit_seconds": args.cpu_limit_seconds,
            "nice_increment": args.nice,
            "parallel_workers": 1,
        },
        "schema": REPORT_SCHEMA,
        "slices": slice_reports,
        "started_utc": started_utc,
        "total_elapsed_seconds": time.monotonic() - started,
        "verifier_inputs_after": verifier_after,
        "verifier_inputs_before": verifier_before,
    }
    report["report_sha256"] = object_hash(report, "report_sha256")
    if output is None:
        print(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=True))
    else:
        write_json(output, report)
        print(canonical_bytes({
            "claim_status": "MACHINE-VERIFIED",
            "outcome": "PASS",
            "report": str(output),
            "report_sha256": report["report_sha256"],
            "slices": selected,
        }).decode("ascii"))
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Structurally audit cert3-trace-v1 without proof arithmetic"
    )
    result.add_argument("--campaign", type=Path, default=DEFAULT_CAMPAIGN)
    result.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    result.add_argument("--slice", type=int, action="append", help="repeat; omit for all")
    result.add_argument("--output", type=Path)
    result.add_argument("--nice", type=int, default=10)
    result.add_argument("--cpu-limit-seconds", type=int, default=0)
    result.add_argument("--chunk-size", type=int, default=4 * 1024 * 1024)
    return result


def main() -> int:
    try:
        return run(parser().parse_args())
    except (AuditError, OSError) as exc:
        print(canonical_bytes({
            "claim_status": "FAILED", "error": str(exc), "outcome": "FAIL_CLOSED"
        }).decode("ascii"), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
