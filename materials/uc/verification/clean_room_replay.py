#!/usr/bin/env python3
"""Clean-room, per-slice replay driver for the frozen Campaign I traces.

This driver is intentionally standard-library-only.  It validates a pinned
campaign lock, stages only the launch/result/trace inputs and the frozen source
snapshot in a fresh read-only tree, then invokes the frozen cert3_replay.py in
an isolated Python subprocess.  It never reads campaign replay markers,
collector output, worker scratch records, bytecode, or import caches.
"""

from __future__ import annotations

import argparse
import datetime as _datetime
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import signal
import stat
import subprocess
import sys
import tempfile
import time
import uuid
from fractions import Fraction


REPORT_SCHEMA = "uc-clean-room-replay-report-v1"
LOCK_SCHEMA = "uc-campaign-lock-v1"
LOCK_RAW_SHA256 = "d254a23a4df7cb3a4c6aae2453883324fba2381c1890c4f011d61c7971b84b48"
CAMPAIGN_CODE_SHA256 = (
    "2f23a58ebdb8b14275284e2ffefca7862373f7baf5d354a06e607cff8d78ce05"
)
CAMPAIGN_LAUNCH_OBJECT_SHA256 = (
    "6414affcaa54bd5e37e9f7bd351c6ed4cf675167b10b55d86cce1436ebdf65cc"
)
CAMPAIGN_LAUNCH_RAW_SHA256 = (
    "63bb6fbe1cfbcb6b7d4c50d7b66882255e5e778ce604c92ecd1abdfa1ea8337e"
)
CAMPAIGN_COLLECTOR_SHA256 = (
    "95d09322a7821f41850ef1ce77345e5e93eb4fa570cc0e7f58f817716cf681ec"
)
CAMPAIGN_DIRECTORY = (
    "cert3_20260818T212601Z_425f109c15b64a6198785c6cebbbdaab_"
    "2f23a58ebdb8"
)
EXPECTED_MODULE_FILE_SHA256 = {
    "flint": "2e5f8f1768d14eccd7961353c635195bd557f297edff8b8de09e2d211f03ec2d",
    "mpmath": "b241584d2c1fc0304b0a1015ea923749d7b0800411dd406dcab7c82bf25d9fe8",
    "numpy": "09295a80660f17925ae23765ce8cbd7ff7ceae968d5f2f89349f1cb74c0b9e11",
    "scipy": "39ccae300a4739cc53719bd3b39ce3f8f66e736c6eb672006e338a4c5cc3dc76",
    "sympy": "4e9476348ba105feab28d82f5bcf6cdba2e3e84de6e059bbfe7a13728c0a4ab0",
}
HERE = Path(__file__).resolve().parent
DEFAULT_LOCK = HERE / "campaign-lock.json"
DEFAULT_CAMPAIGN = HERE.parent / "campaigns" / CAMPAIGN_DIRECTORY
SCHEMA_PATH = HERE / "report.schema.json"

EXECUTABLE_FILES = frozenset(
    {
        "arbcore.py",
        "bound_kkt.py",
        "cert2.py",
        "cert3.py",
        "cert3_par.py",
        "diag_exhaust.py",
        "entropy.py",
    }
)
REVIEW_FILES = frozenset(
    {
        "bridge_uc.py",
        "cert3_collect.py",
        "cert3_replay.py",
        "decomposition.py",
        "lemma_rh_proof.py",
        "margin_lemma.py",
        "reduction.py",
        "thmB3_proof.py",
    }
)
SNAPSHOT_FILES = EXECUTABLE_FILES | REVIEW_FILES
REPLAY_TALLY_NAMES = (
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
HEX64 = re.compile(r"[0-9a-f]{64}\Z")
RUN_ID = re.compile(r"[0-9a-f]{32}\Z")

# This probe runs under the same isolated interpreter and sanitized environment
# used for arithmetic replay.  It imports only external dependencies, before
# the frozen replayer has a chance to add its snapshot to sys.path.
ENVIRONMENT_PROBE = r'''
import hashlib
import json
import os
import platform
import site
import stat
import sys

names = ("flint", "mpmath", "numpy", "scipy", "sympy")
modules = {}
for name in names:
    module = __import__(name)
    filename = getattr(module, "__file__", None)
    digest = None
    size = None
    if filename:
        filename = os.path.realpath(filename)
        try:
            info = os.stat(filename, follow_symlinks=False)
            if stat.S_ISREG(info.st_mode):
                hasher = hashlib.sha256()
                with open(filename, "rb") as handle:
                    while True:
                        chunk = handle.read(1024 * 1024)
                        if not chunk:
                            break
                        hasher.update(chunk)
                digest = hasher.hexdigest()
                size = info.st_size
        except OSError:
            pass
    modules[name] = {
        "file": filename,
        "file_sha256": digest,
        "file_size": size,
        "version": getattr(module, "__version__", "unknown"),
    }
launch_environment = {
    "flint": modules["flint"]["version"],
    "mpmath": modules["mpmath"]["version"],
    "numpy": modules["numpy"]["version"],
    "platform": platform.platform(),
    "python": platform.python_version(),
    "python_implementation": platform.python_implementation(),
    "scipy": modules["scipy"]["version"],
    "sympy": modules["sympy"]["version"],
}
record = {
    "base_prefix": sys.base_prefix,
    "dependencies": modules,
    "executable": os.path.realpath(sys.executable),
    "flags": {
        "dont_write_bytecode": sys.flags.dont_write_bytecode,
        "ignore_environment": sys.flags.ignore_environment,
        "isolated": sys.flags.isolated,
        "no_user_site": sys.flags.no_user_site,
        "safe_path": getattr(sys.flags, "safe_path", None),
    },
    "launch_environment": launch_environment,
    "prefix": sys.prefix,
    "sys_path": list(sys.path),
    "user_site_enabled": site.ENABLE_USER_SITE,
}
print(json.dumps(record, sort_keys=True, separators=(",", ":")))
'''


class VerificationError(RuntimeError):
    """A fail-closed input, environment, staging, or replay rejection."""


def utc_now() -> str:
    return _datetime.datetime.now(_datetime.timezone.utc).isoformat()


def canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def object_sha256(value: dict[str, object], omitted_key: str) -> str:
    reduced = {key: item for key, item in value.items() if key != omitted_key}
    return hashlib.sha256(canonical_bytes(reduced)).hexdigest()


def _reject_duplicate_keys(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise VerificationError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _open_regular(path: Path) -> tuple[int, os.stat_result]:
    flags = os.O_RDONLY
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        raise VerificationError(f"cannot open required regular file {path}: {exc}") from exc
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise VerificationError(f"required input is not a regular file: {path}")
        path_info = path.lstat()
        if stat.S_ISLNK(path_info.st_mode):
            raise VerificationError(f"required input is a symbolic link: {path}")
        return descriptor, before
    except Exception:
        os.close(descriptor)
        raise


def file_fingerprint(path: Path) -> dict[str, object]:
    descriptor, before = _open_regular(path)
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
    identity_before = (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    )
    identity_after = (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    )
    if identity_before != identity_after or size != before.st_size:
        raise VerificationError(f"input changed while it was hashed: {path}")
    return {"sha256": digest.hexdigest(), "size": size}


def load_json_file(path: Path, maximum_size: int = 2 * 1024 * 1024) -> dict[str, object]:
    descriptor, before = _open_regular(path)
    try:
        if before.st_size > maximum_size:
            raise VerificationError(f"JSON input is unexpectedly large: {path}")
        chunks: list[bytes] = []
        remaining = before.st_size
        while remaining:
            chunk = os.read(descriptor, min(remaining, 1024 * 1024))
            if not chunk:
                raise VerificationError(f"JSON input was truncated while reading: {path}")
            chunks.append(chunk)
            remaining -= len(chunk)
        extra = os.read(descriptor, 1)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if extra or (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    ) != (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    ):
        raise VerificationError(f"JSON input changed while it was read: {path}")
    try:
        value = json.loads(
            b"".join(chunks).decode("utf-8"),
            object_pairs_hook=_reject_duplicate_keys,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise VerificationError(f"cannot parse JSON input {path}: {exc}") from exc
    if type(value) is not dict:
        raise VerificationError(f"JSON input must be an object: {path}")
    return value


def require_isolated_interpreter() -> None:
    failures = []
    if sys.flags.isolated != 1:
        failures.append("sys.flags.isolated is not 1 (invoke Python with -I)")
    if sys.flags.no_user_site != 1:
        failures.append("user site is enabled (invoke Python with -I/-s)")
    if sys.flags.ignore_environment != 1:
        failures.append("Python environment variables are active (invoke with -I/-E)")
    if sys.flags.dont_write_bytecode != 1:
        failures.append("bytecode writes are enabled (invoke Python with -B)")
    if failures:
        raise VerificationError("; ".join(failures))


def _safe_relative(value: object, label: str) -> str:
    if type(value) is not str or not value:
        raise VerificationError(f"{label} must be a nonempty relative path")
    path = PurePosixPath(value)
    if path.is_absolute() or any(part in ("", ".", "..") for part in path.parts):
        raise VerificationError(f"{label} is not a safe relative path: {value!r}")
    return value


def _require_digest(value: object, label: str) -> str:
    if type(value) is not str or HEX64.fullmatch(value) is None:
        raise VerificationError(f"{label} is not lowercase SHA-256 hexadecimal")
    return value


def _require_nonnegative_int(value: object, label: str) -> int:
    if type(value) is not int or value < 0:
        raise VerificationError(f"{label} is not a nonnegative integer")
    return value


def _fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _valid_run_id(value: object) -> bool:
    if type(value) is not str or RUN_ID.fullmatch(value) is None:
        return False
    try:
        parsed = uuid.UUID(hex=value)
    except (ValueError, AttributeError):
        return False
    return parsed.hex == value and parsed.version == 4

def _load_pinned_lock(path: Path) -> tuple[dict[str, object], dict[str, object]]:
    """Hash and parse the same no-follow file descriptor, avoiding TOCTOU."""
    descriptor, before = _open_regular(path)
    try:
        if before.st_size != 12590:
            raise VerificationError(
                f"campaign lock size is {before.st_size}, expected 12590"
            )
        payload = bytearray()
        while len(payload) < before.st_size:
            chunk = os.read(descriptor, before.st_size - len(payload))
            if not chunk:
                raise VerificationError("campaign lock was truncated while reading")
            payload.extend(chunk)
        if os.read(descriptor, 1):
            raise VerificationError("campaign lock grew while reading")
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    ) != (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    ):
        raise VerificationError("campaign lock changed while reading")
    fingerprint = {
        "sha256": hashlib.sha256(payload).hexdigest(),
        "size": len(payload),
    }
    if fingerprint != {"sha256": LOCK_RAW_SHA256, "size": 12590}:
        raise VerificationError(
            "campaign lock differs from verifier-pinned bytes: "
            f"expected {LOCK_RAW_SHA256}/12590, got "
            f"{fingerprint['sha256']}/{fingerprint['size']}"
        )
    try:
        value = json.loads(
            payload.decode("utf-8"), object_pairs_hook=_reject_duplicate_keys
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise VerificationError(f"cannot parse campaign lock {path}: {exc}") from exc
    if type(value) is not dict:
        raise VerificationError("campaign lock must be a JSON object")
    return value, fingerprint



def load_and_validate_lock(path: Path) -> tuple[dict[str, object], dict[str, object]]:
    lock, fingerprint = _load_pinned_lock(path)
    if lock.get("schema") != LOCK_SCHEMA:
        raise VerificationError("unsupported campaign lock schema")
    if lock.get("campaign_directory_name") != CAMPAIGN_DIRECTORY:
        raise VerificationError("campaign lock names a different frozen campaign")
    if lock.get("campaign_id") != (
        "20260818T212601Z_425f109c15b64a6198785c6cebbbdaab"
    ):
        raise VerificationError("campaign lock has an unexpected campaign ID")
    if lock.get("code_sha256") != CAMPAIGN_CODE_SHA256:
        raise VerificationError("campaign lock has an unexpected code digest")

    launch_lock = lock.get("launch")
    if type(launch_lock) is not dict:
        raise VerificationError("campaign lock launch entry is not an object")
    if _safe_relative(launch_lock.get("path"), "lock launch path") != "launch.json":
        raise VerificationError("campaign lock launch path is not launch.json")
    if launch_lock.get("raw_sha256") != CAMPAIGN_LAUNCH_RAW_SHA256:
        raise VerificationError("campaign lock has an unexpected raw launch digest")
    if launch_lock.get("object_sha256") != CAMPAIGN_LAUNCH_OBJECT_SHA256:
        raise VerificationError("campaign lock has an unexpected launch object digest")
    if launch_lock.get("size") != 4339:
        raise VerificationError("campaign lock has an unexpected raw launch size")

    snapshot = lock.get("snapshot")
    if type(snapshot) is not dict or set(snapshot) != SNAPSHOT_FILES:
        raise VerificationError("campaign lock snapshot inventory is not the pinned 15-file set")
    for name, entry in snapshot.items():
        if type(entry) is not dict:
            raise VerificationError(f"lock snapshot entry {name} is not an object")
        if _safe_relative(entry.get("path"), f"lock snapshot path {name}") != f"snapshot/{name}":
            raise VerificationError(f"lock snapshot path mismatch for {name}")
        _require_digest(entry.get("sha256"), f"lock snapshot digest {name}")
        _require_nonnegative_int(entry.get("size"), f"lock snapshot size {name}")
    if snapshot["cert3_collect.py"]["sha256"] != CAMPAIGN_COLLECTOR_SHA256:
        raise VerificationError("campaign lock has an unexpected frozen collector digest")

    slices = lock.get("slices")
    if type(slices) is not list or len(slices) != 8:
        raise VerificationError("campaign lock must contain exactly eight slices")
    for position, item in enumerate(slices):
        if type(item) is not dict or item.get("slice") != position:
            raise VerificationError("campaign lock slices must be ordered exactly 0..7")
        run_id = item.get("run_id")
        if not _valid_run_id(run_id):
            raise VerificationError(f"lock slice {position} has an invalid UUID4 run ID")
        expected_lo = _fraction_text(Fraction(1, 2) + Fraction(position, 16))
        expected_hi = _fraction_text(Fraction(1, 2) + Fraction(position + 1, 16))
        if item.get("w_lo") != expected_lo or item.get("w_hi") != expected_hi:
            raise VerificationError(f"lock slice {position} has the wrong exact dyadic root")
        result = item.get("result")
        trace = item.get("trace")
        if type(result) is not dict or type(trace) is not dict:
            raise VerificationError(f"lock slice {position} result/trace entry is malformed")
        expected_result = f"result_slice{position}_{run_id}.json"
        expected_trace = f"trace_slice{position}_{run_id}.bin"
        if _safe_relative(result.get("path"), "lock result path") != expected_result:
            raise VerificationError(f"lock slice {position} result filename mismatch")
        if _safe_relative(trace.get("path"), "lock trace path") != expected_trace:
            raise VerificationError(f"lock slice {position} trace filename mismatch")
        _require_digest(result.get("raw_sha256"), "lock result raw digest")
        _require_digest(result.get("record_sha256"), "lock result object digest")
        _require_digest(trace.get("sha256"), "lock trace digest")
        _require_nonnegative_int(result.get("size"), "lock result size")
        _require_nonnegative_int(trace.get("size"), "lock trace size")
        if trace.get("schema") != "cert3-trace-v1":
            raise VerificationError(f"lock slice {position} has an unexpected trace schema")
        expected = item.get("expected_replay_tallies")
        if type(expected) is not dict or set(expected) != set(REPLAY_TALLY_NAMES):
            raise VerificationError(f"lock slice {position} replay tallies are incomplete")
        for name in REPLAY_TALLY_NAMES:
            _require_nonnegative_int(expected[name], f"lock slice {position} tally {name}")
        completion = item.get("expected_completion_tallies")
        if completion != {"budget_boxes": 0, "budget_time": 0, "stack": 0}:
            raise VerificationError(f"lock slice {position} completion tallies are not zero")
    return lock, fingerprint


def expected_input_fingerprints(
    lock: dict[str, object], slice_index: int
) -> dict[str, dict[str, object]]:
    launch = lock["launch"]
    snapshot = lock["snapshot"]
    item = lock["slices"][slice_index]
    result: dict[str, dict[str, object]] = {
        launch["path"]: {
            "sha256": launch["raw_sha256"],
            "size": launch["size"],
        }
    }
    for entry in snapshot.values():
        result[entry["path"]] = {
            "sha256": entry["sha256"],
            "size": entry["size"],
        }
    result[item["result"]["path"]] = {
        "sha256": item["result"]["raw_sha256"],
        "size": item["result"]["size"],
    }
    result[item["trace"]["path"]] = {
        "sha256": item["trace"]["sha256"],
        "size": item["trace"]["size"],
    }
    return dict(sorted(result.items()))


def validate_campaign_inventory(campaign: Path, lock: dict[str, object]) -> list[str]:
    try:
        campaign_lstat = campaign.lstat()
    except OSError as exc:
        raise VerificationError(f"cannot inspect campaign directory {campaign}: {exc}") from exc
    if stat.S_ISLNK(campaign_lstat.st_mode) or not stat.S_ISDIR(campaign_lstat.st_mode):
        raise VerificationError("campaign path must be a real directory, not a symlink")
    snapshot_dir = campaign / "snapshot"
    try:
        snapshot_lstat = snapshot_dir.lstat()
    except OSError as exc:
        raise VerificationError(f"cannot inspect snapshot directory: {exc}") from exc
    if stat.S_ISLNK(snapshot_lstat.st_mode) or not stat.S_ISDIR(snapshot_lstat.st_mode):
        raise VerificationError("snapshot must be a real directory, not a symlink")
    actual_snapshot = {entry.name for entry in snapshot_dir.iterdir()}
    if actual_snapshot != SNAPSHOT_FILES:
        raise VerificationError(
            "snapshot inventory differs from the pinned sources; "
            f"extra={sorted(actual_snapshot - SNAPSHOT_FILES)} "
            f"missing={sorted(SNAPSHOT_FILES - actual_snapshot)}"
        )

    expected_results = {item["result"]["path"] for item in lock["slices"]}
    expected_traces = {item["trace"]["path"] for item in lock["slices"]}
    entries = list(campaign.iterdir())
    actual_results = {
        entry.name
        for entry in entries
        if entry.name.startswith("result_") and entry.name.endswith(".json")
    }
    actual_traces = {
        entry.name
        for entry in entries
        if entry.name.startswith("trace_") and entry.name.endswith(".bin")
    }
    if actual_results != expected_results:
        raise VerificationError(
            "result inventory differs from launch/lock; "
            f"extra={sorted(actual_results - expected_results)} "
            f"missing={sorted(expected_results - actual_results)}"
        )
    if actual_traces != expected_traces:
        raise VerificationError(
            "trace inventory differs from launch/lock; "
            f"extra={sorted(actual_traces - expected_traces)} "
            f"missing={sorted(expected_traces - actual_traces)}"
        )
    accepted_top_level = {"launch.json", "snapshot"} | expected_results | expected_traces
    return sorted(entry.name for entry in entries if entry.name not in accepted_top_level)


def hash_input_set(
    root: Path, expected: dict[str, dict[str, object]]
) -> dict[str, dict[str, object]]:
    actual: dict[str, dict[str, object]] = {}
    for relative, fingerprint in expected.items():
        observed = file_fingerprint(root / relative)
        if observed != fingerprint:
            raise VerificationError(
                f"input fingerprint mismatch for {relative}: "
                f"expected {fingerprint}, got {observed}"
            )
        actual[relative] = observed
    return actual


def combined_source_hash(per_file: dict[str, str]) -> str:
    digest = hashlib.sha256()
    for name in sorted(per_file):
        digest.update(name.encode("utf-8"))
        digest.update(b"\0")
        digest.update(per_file[name].encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def validate_launch_and_result(
    campaign: Path, lock: dict[str, object], slice_index: int
) -> tuple[dict[str, object], dict[str, object]]:
    launch = load_json_file(campaign / lock["launch"]["path"])
    if launch.get("schema") != "cert3-campaign-v2":
        raise VerificationError("launch schema is not cert3-campaign-v2")
    if launch.get("launch_sha256") != object_sha256(launch, "launch_sha256"):
        raise VerificationError("launch canonical object SHA-256 mismatch")
    if launch.get("launch_sha256") != lock["launch"]["object_sha256"]:
        raise VerificationError("launch object SHA-256 differs from campaign lock")
    if launch.get("campaign_id") != lock["campaign_id"]:
        raise VerificationError("launch campaign ID differs from campaign lock")
    if launch.get("code_sha256") != lock["code_sha256"]:
        raise VerificationError("launch code SHA-256 differs from campaign lock")
    if launch.get("environment") != lock.get("expected_environment"):
        raise VerificationError("launch environment record differs from campaign lock")
    if launch.get("nslices") != 8 or type(launch.get("nslices")) is not int:
        raise VerificationError("launch nslices is not integer 8")
    parameters = launch.get("parameters")
    if type(parameters) is not dict:
        raise VerificationError("launch parameters is not an object")
    if parameters.get("nslices") != 8 or parameters.get("work_prec_bits") != 80:
        raise VerificationError("launch replay precision/slice parameters are unexpected")
    if parameters.get("t_decimal") != launch.get("t_decimal"):
        raise VerificationError("launch target differs from parameters target")
    if parameters.get("level_offset") != launch.get("level_offset"):
        raise VerificationError("launch level offset differs from parameters offset")
    try:
        delta = Fraction(launch["level_offset"])
        target = Fraction(launch["t_decimal"])
    except (KeyError, TypeError, ValueError, ZeroDivisionError) as exc:
        raise VerificationError(f"launch target/offset is not exact rational text: {exc}") from exc
    reduced_target = target - delta
    if not (
        delta > 0
        and 0 <= reduced_target <= Fraction(3, 2)
        and reduced_target * reduced_target - 3 * reduced_target + 1 <= 0
    ):
        raise VerificationError("launch target does not dominate exact psi plus offset")

    executable = launch.get("executable_files")
    review = launch.get("review_files")
    if type(executable) is not dict or set(executable) != EXECUTABLE_FILES:
        raise VerificationError("launch executable source inventory is unexpected")
    if type(review) is not dict or set(review) != REVIEW_FILES:
        raise VerificationError("launch review source inventory is unexpected")
    locked_snapshot = lock["snapshot"]
    for name in SNAPSHOT_FILES:
        launch_digest = executable.get(name, review.get(name))
        if launch_digest != locked_snapshot[name]["sha256"]:
            raise VerificationError(f"launch/lock source digest mismatch for {name}")
    if combined_source_hash(executable) != launch.get("code_sha256"):
        raise VerificationError("launch combined executable code digest mismatch")

    runs = launch.get("runs")
    if type(runs) is not list or len(runs) != 8:
        raise VerificationError("launch does not contain exactly eight runs")
    for position, (run, locked) in enumerate(zip(runs, lock["slices"])):
        if type(run) is not dict or run.get("slice") != position:
            raise VerificationError("launch runs are not ordered exactly 0..7")
        expected = {
            "result_filename": locked["result"]["path"],
            "run_id": locked["run_id"],
            "slice": position,
            "trace_filename": locked["trace"]["path"],
            "w_hi": locked["w_hi"],
            "w_lo": locked["w_lo"],
        }
        for key, value in expected.items():
            if run.get(key) != value:
                raise VerificationError(f"launch slice {position} field {key} differs from lock")
        if type(run.get("time_budget_s")) is not int or run["time_budget_s"] <= 0:
            raise VerificationError(f"launch slice {position} time budget is invalid")

    locked_slice = lock["slices"][slice_index]
    run = runs[slice_index]
    result_path = campaign / locked_slice["result"]["path"]
    result = load_json_file(result_path)
    if result.get("result_schema") != "cert3-result-v2":
        raise VerificationError("result schema is not cert3-result-v2")
    if result.get("worker_schema") != "cert3-worker-v2":
        raise VerificationError("worker schema is not cert3-worker-v2")
    if result.get("record_sha256") != object_sha256(result, "record_sha256"):
        raise VerificationError("result canonical object SHA-256 mismatch")
    if result.get("record_sha256") != locked_slice["result"]["record_sha256"]:
        raise VerificationError("result object SHA-256 differs from campaign lock")
    worker_view = {
        key: value
        for key, value in result.items()
        if key
        not in (
            "record_sha256",
            "result_schema",
            "supervisor_committed",
            "worker_exit_code",
        )
    }
    if worker_view.get("worker_record_sha256") != object_sha256(
        worker_view, "worker_record_sha256"
    ):
        raise VerificationError("embedded worker object SHA-256 mismatch")
    common = {
        "campaign_id": launch["campaign_id"],
        "code_sha256": launch["code_sha256"],
        "launch_sha256": launch["launch_sha256"],
        "level_offset": launch["level_offset"],
        "nslices": 8,
        "result_filename": run["result_filename"],
        "run_id": run["run_id"],
        "slice": slice_index,
        "t_decimal": launch["t_decimal"],
        "time_budget_s": run["time_budget_s"],
        "trace_filename": run["trace_filename"],
        "w_hi": run["w_hi"],
        "w_lo": run["w_lo"],
    }
    for key, value in common.items():
        if result.get(key) != value:
            raise VerificationError(f"result field {key} differs from launch")
    for key in ("complete", "returned_normally", "supervisor_committed"):
        if result.get(key) is not True:
            raise VerificationError(f"result field {key} is not true")
    if result.get("trace_schema") != "cert3-trace-v1":
        raise VerificationError("result trace_schema is not cert3-trace-v1")
    if result.get("verdict") != "COMPLETE":
        raise VerificationError("result verdict is not COMPLETE")
    if type(result.get("worker_exit_code")) is not int or result["worker_exit_code"] != 0:
        raise VerificationError("result worker_exit_code is not integer 0")
    if result.get("trace_sha256") != locked_slice["trace"]["sha256"]:
        raise VerificationError("result trace digest differs from campaign lock")
    if result.get("trace_size") != locked_slice["trace"]["size"]:
        raise VerificationError("result trace size differs from campaign lock")
    tallies = result.get("tallies")
    if type(tallies) is not dict:
        raise VerificationError("result tallies is not an object")
    for name, value in tallies.items():
        if name != "residual_geometry" and type(value) is not int:
            raise VerificationError(f"result tally {name} is not an integer")
    for name, expected_value in locked_slice["expected_replay_tallies"].items():
        if tallies.get(name) != expected_value:
            raise VerificationError(f"result replay tally {name} differs from campaign lock")
    for name, expected_value in locked_slice["expected_completion_tallies"].items():
        if tallies.get(name) != expected_value:
            raise VerificationError(f"result completion tally {name} is not zero")
    if tallies.get("residual_geometry") != {}:
        raise VerificationError("result residual_geometry is not empty")
    if tallies.get("processed") != result.get("trace_size"):
        raise VerificationError("result processed tally differs from trace size")
    return launch, result


def _copy_regular_verified(source: Path, target: Path, expected: dict[str, object]) -> None:
    source_fd, before = _open_regular(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    target_fd = os.open(target, flags, 0o600)
    digest = hashlib.sha256()
    copied = 0
    try:
        while True:
            chunk = os.read(source_fd, 4 * 1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
            copied += len(chunk)
            view = memoryview(chunk)
            while view:
                written = os.write(target_fd, view)
                if written <= 0:
                    raise VerificationError(f"short write while staging {source}")
                view = view[written:]
        os.fsync(target_fd)
        after = os.fstat(source_fd)
    finally:
        os.close(source_fd)
        os.close(target_fd)
    if (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    ) != (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    ):
        raise VerificationError(f"source changed while staging {source}")
    observed = {"sha256": digest.hexdigest(), "size": copied}
    if observed != expected:
        raise VerificationError(
            f"source fingerprint changed while staging {source}: {observed}"
        )
    os.chmod(target, 0o444)


def stage_input_set(
    source_root: Path,
    stage_root: Path,
    expected: dict[str, dict[str, object]],
) -> Path:
    staged_campaign = stage_root / "campaign"
    staged_campaign.mkdir(mode=0o700, parents=True)
    for relative, fingerprint in expected.items():
        _copy_regular_verified(
            source_root / relative, staged_campaign / relative, fingerprint
        )
    os.chmod(staged_campaign / "snapshot", 0o555)
    os.chmod(staged_campaign, 0o555)
    return staged_campaign


def validate_stage_inventory(staged_campaign: Path, expected: dict[str, object]) -> None:
    expected_top = {PurePosixPath(name).parts[0] for name in expected}
    actual_top = {entry.name for entry in staged_campaign.iterdir()}
    if actual_top != expected_top:
        raise VerificationError(
            f"read-only stage top-level inventory changed: {sorted(actual_top)}"
        )
    expected_snapshot = {
        PurePosixPath(name).name for name in expected if name.startswith("snapshot/")
    }
    actual_snapshot = {entry.name for entry in (staged_campaign / "snapshot").iterdir()}
    if actual_snapshot != expected_snapshot:
        raise VerificationError(
            f"read-only stage snapshot inventory changed: {sorted(actual_snapshot)}"
        )


def make_stage_removable(staged_campaign: Path) -> None:
    try:
        os.chmod(staged_campaign, 0o700)
    except OSError:
        return
    snapshot = staged_campaign / "snapshot"
    if snapshot.exists():
        try:
            os.chmod(snapshot, 0o700)
        except OSError:
            pass


def sanitized_environment(runtime_root: Path) -> dict[str, str]:
    home = runtime_root / "home"
    temporary = runtime_root / "tmp"
    cwd = runtime_root / "cwd"
    for directory in (home, temporary, cwd):
        directory.mkdir(mode=0o700, parents=True, exist_ok=False)
    environment = {
        "HOME": str(home),
        "PATH": "/usr/bin:/bin:/usr/sbin:/sbin",
        "TMPDIR": str(temporary),
        "PYTHONDONTWRITEBYTECODE": "1",
        "PYTHONNOUSERSITE": "1",
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
        "BLIS_NUM_THREADS": "1",
    }
    for name in ("LANG", "LC_ALL", "LC_CTYPE", "TZ", "SYSTEM_VERSION_COMPAT"):
        value = os.environ.get(name)
        if value is not None:
            environment[name] = value
    return environment


def capture_environment(environment: dict[str, str], cwd: Path) -> dict[str, object]:
    command = [sys.executable, "-I", "-B", "-s", "-E", "-c", ENVIRONMENT_PROBE]
    completed = subprocess.run(
        command,
        cwd=cwd,
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        timeout=60,
    )
    if completed.returncode != 0:
        raise VerificationError(
            "isolated environment probe failed with exit code "
            f"{completed.returncode}: {completed.stderr.decode('utf-8', 'replace')}"
        )
    if completed.stderr:
        raise VerificationError(
            "isolated environment probe emitted unexpected stderr: "
            + completed.stderr.decode("utf-8", "replace")
        )
    lines = completed.stdout.splitlines()
    if len(lines) != 1:
        raise VerificationError("isolated environment probe did not emit one JSON line")
    try:
        record = json.loads(
            lines[0].decode("ascii"), object_pairs_hook=_reject_duplicate_keys
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise VerificationError(f"cannot parse isolated environment probe: {exc}") from exc
    if type(record) is not dict:
        raise VerificationError("isolated environment probe output is not an object")
    expected_flags = {
        "dont_write_bytecode": 1,
        "ignore_environment": 1,
        "isolated": 1,
        "no_user_site": 1,
        "safe_path": True,
    }
    if record.get("flags") != expected_flags:
        raise VerificationError(
            f"child interpreter isolation flags are unexpected: {record.get('flags')}"
        )
    if record.get("user_site_enabled") is not False:
        raise VerificationError("child interpreter reports an enabled user site")
    return record


def validate_environment(
    observed: dict[str, object], lock: dict[str, object]
) -> None:
    if observed.get("launch_environment") != lock.get("expected_environment"):
        raise VerificationError(
            "isolated interpreter/platform/dependency versions differ from launch: "
            f"expected {lock.get('expected_environment')}, "
            f"got {observed.get('launch_environment')}"
        )
    dependencies = observed.get("dependencies")
    if type(dependencies) is not dict or set(dependencies) != {
        "flint",
        "mpmath",
        "numpy",
        "scipy",
        "sympy",
    }:
        raise VerificationError("environment dependency inventory is unexpected")
    for name, item in dependencies.items():
        if type(item) is not dict or not item.get("file"):
            raise VerificationError(f"dependency {name} has no resolved module file")
        digest = _require_digest(
            item.get("file_sha256"), f"dependency {name} module-file digest"
        )
        if digest != EXPECTED_MODULE_FILE_SHA256[name]:
            raise VerificationError(
                f"dependency {name} module bytes differ from the audited environment"
            )


def verifier_fingerprints(lock_path: Path) -> dict[str, object]:
    if not SCHEMA_PATH.exists():
        raise VerificationError(f"missing report schema: {SCHEMA_PATH}")
    return {
        "campaign_lock": file_fingerprint(lock_path),
        "clean_room_replay.py": file_fingerprint(Path(__file__).resolve()),
        "report.schema.json": file_fingerprint(SCHEMA_PATH),
    }


def _child_setup(nice_increment: int, cpu_limit: int, address_space_gib: int) -> None:
    import resource

    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    if cpu_limit:
        hard = cpu_limit + 60
        resource.setrlimit(resource.RLIMIT_CPU, (cpu_limit, hard))
    if address_space_gib:
        limit = address_space_gib * (1 << 30)
        resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    if nice_increment:
        os.nice(nice_increment)


def run_replayer(
    command: list[str],
    environment: dict[str, str],
    cwd: Path,
    nice_increment: int,
    cpu_limit: int,
    wall_timeout: int,
    address_space_gib: int,
) -> dict[str, object]:
    started = time.monotonic()
    process = subprocess.Popen(
        command,
        cwd=cwd,
        env=environment,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,
        preexec_fn=lambda: _child_setup(
            nice_increment, cpu_limit, address_space_gib
        ),
    )
    termination = "exit"
    try:
        stdout, stderr = process.communicate(
            timeout=None if wall_timeout == 0 else wall_timeout
        )
    except subprocess.TimeoutExpired:
        termination = "wall-timeout"
        os.killpg(process.pid, signal.SIGTERM)
        try:
            stdout, stderr = process.communicate(timeout=10)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            stdout, stderr = process.communicate()
    except KeyboardInterrupt:
        termination = "interrupted"
        os.killpg(process.pid, signal.SIGTERM)
        try:
            stdout, stderr = process.communicate(timeout=10)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            stdout, stderr = process.communicate()
    elapsed = time.monotonic() - started
    return {
        "elapsed_seconds": elapsed,
        "exit_code": process.returncode,
        "stderr_sha256": hashlib.sha256(stderr).hexdigest(),
        "stderr_utf8": stderr.decode("utf-8", "replace"),
        "stdout_sha256": hashlib.sha256(stdout).hexdigest(),
        "stdout_utf8": stdout.decode("utf-8", "replace"),
        "termination": termination,
    }


def parse_replay_output(
    execution: dict[str, object], expected: dict[str, int]
) -> tuple[dict[str, int] | None, list[str]]:
    problems: list[str] = []
    if execution["termination"] != "exit":
        problems.append(f"replayer termination was {execution['termination']}")
    if execution["exit_code"] != 0:
        problems.append(f"replayer exit code was {execution['exit_code']}, not 0")
    if execution["stderr_utf8"]:
        problems.append("replayer emitted unexpected stderr")
    lines = execution["stdout_utf8"].splitlines()
    replayed: dict[str, int] | None = None
    if len(lines) != 1:
        problems.append("replayer stdout was not exactly one JSON line")
    else:
        try:
            value = json.loads(lines[0], object_pairs_hook=_reject_duplicate_keys)
        except (json.JSONDecodeError, VerificationError) as exc:
            problems.append(f"replayer stdout JSON could not be parsed: {exc}")
        else:
            if type(value) is not dict:
                problems.append("replayer stdout JSON was not an object")
            elif set(value) != set(REPLAY_TALLY_NAMES):
                problems.append("replayer stdout tally inventory was unexpected")
            elif any(type(value[name]) is not int for name in REPLAY_TALLY_NAMES):
                problems.append("replayer stdout contained a non-integer tally")
            else:
                replayed = {name: value[name] for name in REPLAY_TALLY_NAMES}
                if replayed != expected:
                    problems.append("replayed tallies differed from locked result tallies")
    return replayed, problems


def _append_jsonl(path: Path, value: dict[str, object]) -> None:
    flags = os.O_WRONLY | os.O_APPEND | os.O_CREAT
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags, 0o600)
    try:
        os.write(descriptor, canonical_bytes(value) + b"\n")
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _write_json_atomic(path: Path, value: dict[str, object]) -> None:
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    payload = json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True).encode(
        "ascii"
    ) + b"\n"
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    descriptor = os.open(temporary, flags, 0o600)
    try:
        view = memoryview(payload)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise VerificationError(f"short write for report {path}")
            view = view[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    os.replace(temporary, path)
    directory_fd = os.open(path.parent, os.O_RDONLY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)


def _inside(path: Path, parent: Path) -> bool:
    return path == parent or parent in path.parents


def prepare_output_directory(output: Path, campaign: Path, slice_index: int) -> Path:
    output = output.expanduser().resolve()
    if _inside(output, campaign):
        raise VerificationError("output directory must not be inside the frozen campaign")
    output.mkdir(mode=0o700, parents=True, exist_ok=True)
    if output.is_symlink() or not output.is_dir():
        raise VerificationError("output path must be a real directory")
    slice_dir = output / f"slice-{slice_index}"
    slice_dir.mkdir(mode=0o700, exist_ok=True)
    if slice_dir.is_symlink() or not slice_dir.is_dir():
        raise VerificationError("per-slice output path must be a real directory")
    return slice_dir




def run(args: argparse.Namespace) -> int:
    require_isolated_interpreter()
    lock_path = args.lock.expanduser().resolve()
    lock, _ = load_and_validate_lock(lock_path)
    if not 0 <= args.slice < 8:
        raise VerificationError("--slice must be an integer in 0..7")
    if args.resume:
        raise VerificationError(
            "--resume is disabled: repository reports are publicly forgeable "
            "without an external attestation key; run a fresh replay"
        )
    if not 0 <= args.nice <= 19:
        raise VerificationError("--nice must be in 0..19")
    for name in ("cpu_limit_seconds", "wall_timeout_seconds", "address_space_gib"):
        if getattr(args, name) < 0:
            raise VerificationError(f"--{name.replace('_', '-')} must be nonnegative")

    campaign = args.campaign.expanduser().resolve()
    ignored_non_inputs = validate_campaign_inventory(campaign, lock)
    expected = expected_input_fingerprints(lock, args.slice)
    source_before = hash_input_set(campaign, expected)
    launch, result = validate_launch_and_result(campaign, lock, args.slice)
    tools_before = verifier_fingerprints(lock_path)
    slice_dir = prepare_output_directory(args.output_dir, campaign, args.slice)
    journal = slice_dir / "events.jsonl"
    attempt_id = uuid.uuid4().hex
    started_utc = utc_now()
    started_monotonic = time.monotonic()
    _append_jsonl(
        journal,
        {
            "attempt_id": attempt_id,
            "event": "START",
            "preflight_only": args.preflight_only,
            "slice": args.slice,
            "utc": started_utc,
        },
    )

    staging_parent: Path | None = None
    if args.staging_root is not None:
        staging_parent = args.staging_root.expanduser().resolve()
        if _inside(staging_parent, campaign):
            raise VerificationError("staging root must not be inside the frozen campaign")
        staging_parent.mkdir(mode=0o700, parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(
        prefix=f"uc-campaign-i-slice-{args.slice}-",
        dir=None if staging_parent is None else staging_parent,
    ) as temporary_name:
        temporary_root = Path(temporary_name)
        runtime_root = temporary_root / "runtime"
        runtime_root.mkdir(mode=0o700)
        environment = sanitized_environment(runtime_root)
        runtime_cwd = runtime_root / "cwd"
        environment_before = capture_environment(environment, runtime_cwd)
        validate_environment(environment_before, lock)

        latest_report = slice_dir / "latest.json"

        staged_campaign: Path | None = None
        execution: dict[str, object] | None = None
        replayed: dict[str, int] | None = None
        errors: list[str] = []
        try:
            staged_campaign = stage_input_set(
                campaign, temporary_root / "staging", expected
            )
            stage_before = hash_input_set(staged_campaign, expected)
            validate_stage_inventory(staged_campaign, expected)
            if not args.preflight_only:
                locked_slice = lock["slices"][args.slice]
                command = [
                    sys.executable,
                    "-I",
                    "-B",
                    "-s",
                    "-E",
                    str(staged_campaign / "snapshot" / "cert3_replay.py"),
                    str(staged_campaign / "launch.json"),
                    str(staged_campaign / locked_slice["result"]["path"]),
                    str(staged_campaign / locked_slice["trace"]["path"]),
                ]
                execution = run_replayer(
                    command,
                    environment,
                    runtime_cwd,
                    args.nice,
                    args.cpu_limit_seconds,
                    args.wall_timeout_seconds,
                    args.address_space_gib,
                )
                replayed, replay_problems = parse_replay_output(
                    execution, locked_slice["expected_replay_tallies"]
                )
                errors.extend(replay_problems)
            stage_after = hash_input_set(staged_campaign, expected)
            validate_stage_inventory(staged_campaign, expected)
            source_after = hash_input_set(campaign, expected)
            environment_after = capture_environment(environment, runtime_cwd)
            validate_environment(environment_after, lock)
            if source_before != source_after:
                errors.append("source input hashes changed during the attempt")
            if stage_before != stage_after:
                errors.append("read-only staged input hashes changed during the attempt")
            if environment_before != environment_after:
                errors.append("interpreter/platform/dependency environment changed during attempt")
            tools_after = verifier_fingerprints(lock_path)
            if tools_before != tools_after:
                errors.append("verifier, schema, or external pin bytes changed during attempt")
        finally:
            if staged_campaign is not None:
                make_stage_removable(staged_campaign)

        elapsed_seconds = time.monotonic() - started_monotonic
        if errors:
            outcome = "FAIL"
            claim_status = "FAILED"
        elif args.preflight_only:
            outcome = "PREFLIGHT_PASS"
            claim_status = "OPEN"
        else:
            outcome = "PASS"
            claim_status = "MACHINE-VERIFIED"
        locked_slice = lock["slices"][args.slice]
        report: dict[str, object] = {
            "$schema": str(SCHEMA_PATH),
            "attempt_id": attempt_id,
            "campaign_id": lock["campaign_id"],
            "claim_status": claim_status,
            "command": execution,
            "environment_after": environment_after,
            "environment_before": environment_before,
            "environment_matches_launch": True,
            "errors": errors,
            "finished_utc": utc_now(),
            "ignored_campaign_non_inputs": ignored_non_inputs,
            "inputs": {
                "source_after": source_after,
                "source_before": source_before,
                "stage_after": stage_after,
                "stage_before": stage_before,
            },
            "limitations": [
                "Arithmetic replay uses the frozen cert3_replay.py and the same frozen interval-arithmetic modules used by the producer; it is not an independently implemented interval proof.",
                "Repository-local reports are self-consistent records, not authenticated execution attestations; fresh replay is required.",
            ],
            "outcome": outcome,
            "report_type": "clean_room_replay",
            "resource_controls": {
                "address_space_gib": args.address_space_gib,
                "cpu_limit_seconds": args.cpu_limit_seconds,
                "nice_increment": args.nice,
                "thread_limits": {
                    name: environment[name]
                    for name in (
                        "BLIS_NUM_THREADS",
                        "MKL_NUM_THREADS",
                        "NUMEXPR_NUM_THREADS",
                        "OMP_NUM_THREADS",
                        "OPENBLAS_NUM_THREADS",
                        "VECLIB_MAXIMUM_THREADS",
                    )
                },
                "wall_timeout_seconds": args.wall_timeout_seconds,
            },
            "root_exact": [
                ["0", "1"],
                ["0", "1"],
                ["0", "1"],
                ["0", "1"],
                [locked_slice["w_lo"], locked_slice["w_hi"]],
            ],
            "run_id": locked_slice["run_id"],
            "schema": REPORT_SCHEMA,
            "slice": args.slice,
            "started_utc": started_utc,
            "total_elapsed_seconds": elapsed_seconds,
            "verifier_inputs_after": tools_after,
            "verifier_inputs_before": tools_before,
            "replay": None
            if args.preflight_only
            else {
                "implementation": "snapshot/cert3_replay.py",
                "implementation_sha256": lock["snapshot"]["cert3_replay.py"]["sha256"],
                "shares_frozen_arithmetic": True,
                "tallies": replayed,
            },
        }
        report["report_sha256"] = object_sha256(report, "report_sha256")
        report_path = (
            slice_dir / f"preflight-{attempt_id}.json"
            if args.preflight_only
            else latest_report
        )
        _write_json_atomic(report_path, report)
        _append_jsonl(
            journal,
            {
                "attempt_id": attempt_id,
                "claim_status": claim_status,
                "event": "FINISH",
                "outcome": outcome,
                "report": str(report_path),
                "report_sha256": report["report_sha256"],
                "slice": args.slice,
                "utc": report["finished_utc"],
            },
        )
        print(
            canonical_bytes(
                {
                    "claim_status": claim_status,
                    "outcome": outcome,
                    "report": str(report_path),
                    "report_sha256": report["report_sha256"],
                    "slice": args.slice,
                }
            ).decode("ascii")
        )
        return 0 if not errors else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Replay one frozen Campaign I trace in a fresh isolated stage"
    )
    parser.add_argument("--slice", type=int, required=True, help="slice index 0..7")
    parser.add_argument(
        "--campaign", type=Path, default=DEFAULT_CAMPAIGN, help="frozen campaign directory"
    )
    parser.add_argument(
        "--lock", type=Path, default=DEFAULT_LOCK, help="pinned campaign lock"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True, help="durable per-slice report directory"
    )
    parser.add_argument(
        "--staging-root", type=Path, help="optional filesystem for fresh temporary copies"
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="disabled: fresh replay is required for proof acceptance",
    )
    parser.add_argument(
        "--preflight-only",
        action="store_true",
        help="stage, hash, and probe without starting arithmetic replay",
    )
    parser.add_argument(
        "--nice", type=int, default=10, help="positive child nice increment (default: 10)"
    )
    parser.add_argument(
        "--cpu-limit-seconds",
        type=int,
        default=0,
        help="per-child RLIMIT_CPU in seconds; 0 disables",
    )
    parser.add_argument(
        "--wall-timeout-seconds",
        type=int,
        default=0,
        help="terminate the child after this wall time; 0 disables",
    )
    parser.add_argument(
        "--address-space-gib",
        type=int,
        default=0,
        help="optional child RLIMIT_AS in GiB; 0 disables",
    )
    return parser


def main() -> int:
    try:
        return run(build_parser().parse_args())
    except (VerificationError, OSError) as exc:
        print(
            canonical_bytes(
                {
                    "claim_status": "FAILED",
                    "error": str(exc),
                    "outcome": "FAIL_CLOSED",
                }
            ).decode("ascii"),
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
