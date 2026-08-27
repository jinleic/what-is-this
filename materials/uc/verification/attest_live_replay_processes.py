#!/usr/bin/env python3
"""Compare live secure replay process images with the content-addressed seal."""

from __future__ import annotations

import argparse
import datetime as _datetime
import hashlib
import json
from pathlib import Path
import stat
import subprocess
import sys

SCHEMA = "uc-live-replay-runtime-attestation-v1"


class AttestationError(RuntimeError):
    """A fail-closed live-process inspection error."""


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


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        while True:
            block = stream.read(1024 * 1024)
            if not block:
                return digest.hexdigest()
            digest.update(block)


def load_json(path, label):
    path = Path(path)
    if not path.is_file() or path.is_symlink():
        raise AttestationError("%s must be a regular file" % label)
    raw = path.read_bytes()
    try:
        value = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise AttestationError("%s is malformed JSON" % label) from exc
    return raw, value


def inspect_process(pid, slice_index, sealed_files):
    command_result = subprocess.run(
        ["/bin/ps", "-p", str(pid), "-o", "command="],
        capture_output=True,
        text=True,
        timeout=30,
    )
    if command_result.returncode != 0 or not command_result.stdout.strip():
        raise AttestationError("secure slice %d process is not live" % slice_index)
    command = command_result.stdout.strip()
    required_fragments = (
        "independent_arithmetic_replay_secure.py",
        "--runtime-seal",
        "runtime-environment-seal-secure.json",
        "--slice %d" % slice_index,
    )
    if any(fragment not in command for fragment in required_fragments):
        raise AttestationError("secure slice %d command differs" % slice_index)
    if "--resume" in command:
        raise AttestationError("secure slice %d process uses resume" % slice_index)

    lsof_result = subprocess.run(
        ["/usr/sbin/lsof", "-p", str(pid), "-Fn"],
        capture_output=True,
        text=True,
        timeout=30,
    )
    if lsof_result.returncode != 0:
        raise AttestationError("lsof failed for secure slice %d" % slice_index)
    current_type = None
    loaded = []
    for line in lsof_result.stdout.splitlines():
        if line.startswith("f"):
            current_type = line[1:]
        elif line.startswith("n") and current_type == "txt":
            loaded.append(line[1:])

    checked = []
    system_unsealed = []
    for text in sorted(set(loaded)):
        path = Path(text)
        if text.startswith(("/usr/lib/", "/System/")):
            system_unsealed.append(text)
            continue
        try:
            resolved = path.resolve(strict=True)
            metadata = resolved.lstat()
        except OSError as exc:
            raise AttestationError("cannot resolve loaded image %s" % text) from exc
        if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISREG(metadata.st_mode):
            raise AttestationError("loaded image is not regular: %s" % resolved)
        sealed = sealed_files.get(str(resolved))
        digest = sha256_file(resolved)
        if (sealed is not None
                and (metadata.st_size != sealed.get("size")
                     or digest != sealed.get("sha256"))):
            raise AttestationError("loaded image bytes differ from seal: %s" % resolved)
        checked.append({
            "path": str(resolved),
            "size": metadata.st_size,
            "sha256": digest,
            "present_in_prelaunch_seal": sealed is not None,
        })
    return {
        "slice": slice_index,
        "pid": pid,
        "command": command,
        "loaded_non_system": checked,
        "system_unsealed": system_unsealed,
        "lsof_sha256": hashlib.sha256(lsof_result.stdout.encode()).hexdigest(),
    }


def attest(runtime_seal_path, launch_manifest_path):
    seal_raw, seal = load_json(runtime_seal_path, "runtime seal")
    launch_raw, launch = load_json(launch_manifest_path, "secure launch manifest")
    if (seal.get("schema") != "uc-runtime-environment-seal-v1"
            or seal.get("outcome") != "PASS"
            or seal.get("report_sha256") != object_hash(seal)):
        raise AttestationError("runtime seal envelope differs")
    if (launch.get("schema") != "uc-secure-independent-arithmetic-launch-v1"
            or launch.get("resume_used") is not False
            or launch.get("initial_trace_offset") != 0
            or launch.get("runtime_seal_raw_sha256")
            != hashlib.sha256(seal_raw).hexdigest()):
        raise AttestationError("secure launch manifest envelope differs")
    sealed_files = {
        item["path"]: item for item in seal.get("files", [])
        if isinstance(item, dict) and isinstance(item.get("path"), str)
    }
    workers = launch.get("workers")
    if not isinstance(workers, list) or len(workers) != 8:
        raise AttestationError("secure launch manifest needs eight workers")
    process_reports = []
    for worker in workers:
        process_reports.append(
            inspect_process(worker["pid"], worker["slice"], sealed_files)
        )
    result = {
        "schema": SCHEMA,
        "claim_status": "MACHINE-VERIFIED",
        "outcome": "PASS",
        "created_utc": _datetime.datetime.now(_datetime.timezone.utc).isoformat(),
        "runtime_seal_raw_sha256": hashlib.sha256(seal_raw).hexdigest(),
        "runtime_seal_report_sha256": seal["report_sha256"],
        "launch_manifest_raw_sha256": hashlib.sha256(launch_raw).hexdigest(),
        "processes": process_reports,
        "process_count": len(process_reports),
        "limitations": [
            "This is an operating-system process inspection, not hardware-backed attestation.",
            "macOS dyld-shared-cache system images are listed but not hashed.",
            "Loaded non-system images absent from the prelaunch seal are hashed and identified in this mid-run report rather than silently accepted.",
            "A malicious kernel, lsof, ps, interpreter, or hash implementation could forge this report."
        ],
    }
    result["report_sha256"] = object_hash(result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description="Attest live secure UC replay images")
    parser.add_argument("--runtime-seal", required=True)
    parser.add_argument("--launch-manifest", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    try:
        result = attest(args.runtime_seal, args.launch_manifest)
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("xb") as stream:
            stream.write(canonical_bytes(result) + b"\n")
    except (AttestationError, OSError, subprocess.SubprocessError, ValueError, TypeError) as exc:
        print("live runtime attestation FAIL: %s" % exc, file=sys.stderr)
        return 1
    print(canonical_bytes({
        "claim_status": result["claim_status"],
        "outcome": result["outcome"],
        "process_count": result["process_count"],
        "report": str(output),
        "report_sha256": result["report_sha256"],
    }).decode("ascii"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
