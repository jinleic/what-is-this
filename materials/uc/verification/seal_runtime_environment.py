#!/usr/bin/env python3
"""Content-address the observed CPython/Arb verifier runtime.

This is an integrity manifest, not a sandbox or execution attestation.  It
loads the standalone independent verifier, initializes its Arb covers, hashes
every imported module file, and hashes every regular file shipped by the
python-flint distribution.
"""

from __future__ import annotations

import argparse
import datetime as _datetime
import hashlib
import importlib.metadata as importlib_metadata
import importlib.util
import json
import os
from pathlib import Path
import platform
import stat
import subprocess
import sys

SCHEMA = "uc-runtime-environment-seal-v1"
THREAD_VARIABLES = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "PYTHONHASHSEED",
    "PYTHONNOUSERSITE",
)


class SealError(RuntimeError):
    """A fail-closed sealing error."""


def canonical_bytes(value):
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("ascii")


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        while True:
            block = stream.read(1024 * 1024)
            if not block:
                return digest.hexdigest()
            digest.update(block)


def regular_file(path):
    path = Path(path)
    try:
        metadata = path.lstat()
    except OSError as exc:
        raise SealError("cannot stat %s" % path) from exc
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISREG(metadata.st_mode):
        raise SealError("expected a regular non-symlink file: %s" % path)
    return metadata


def file_record(path, origins):
    path = Path(path).resolve(strict=True)
    metadata = regular_file(path)
    return {
        "path": str(path),
        "size": metadata.st_size,
        "sha256": sha256_file(path),
        "origins": sorted(origins),
    }


def load_verifier(path):
    path = Path(path).resolve(strict=True)
    before = set(sys.modules)
    spec = importlib.util.spec_from_file_location("uc_independent_sealed", path)
    if spec is None or spec.loader is None:
        raise SealError("cannot construct independent-verifier import spec")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module._ensure_global_caps()
    return module, before


def module_files(before):
    records = {}
    for name, module in sorted(sys.modules.items()):
        primary = getattr(module, "__file__", None)
        candidates = []
        if isinstance(primary, str):
            candidates.append((primary, "module-file"))
        cached = getattr(module, "__cached__", None)
        if isinstance(cached, str) and Path(cached).exists():
            candidates.append((cached, "module-cache"))
        spec = getattr(module, "__spec__", None)
        origin = getattr(spec, "origin", None)
        if (isinstance(origin, str)
                and origin not in {"built-in", "frozen"}
                and Path(origin).exists()):
            candidates.append((origin, "module-origin"))
        for candidate, origin_kind in candidates:
            try:
                resolved = Path(candidate).resolve(strict=True)
                metadata = resolved.lstat()
            except OSError as exc:
                raise SealError("cannot resolve imported module file %s" % candidate) from exc
            if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISREG(metadata.st_mode):
                raise SealError("imported module path is not a regular file: %s" % resolved)
            origins = records.setdefault(resolved, set())
            origins.add(
                "newly_imported_module" if name not in before else "startup_module"
            )
            origins.add("module:%s" % name)
            origins.add(origin_kind)
    return records


def distribution_files(records):
    try:
        distribution = importlib_metadata.distribution("python-flint")
    except importlib_metadata.PackageNotFoundError as exc:
        raise SealError("python-flint distribution metadata is unavailable") from exc
    for item in distribution.files or ():
        path = Path(distribution.locate_file(item))
        try:
            resolved = path.resolve(strict=True)
            metadata = resolved.lstat()
        except OSError as exc:
            raise SealError(
                "cannot resolve python-flint distribution file %s" % path
            ) from exc
        if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISREG(metadata.st_mode):
            raise SealError(
                "python-flint distribution path is not a regular file: %s"
                % resolved
            )
        records.setdefault(resolved, set()).add("python-flint-distribution")
    return distribution.version


def otool_dependencies(paths):
    tool = Path("/usr/bin/otool")
    if not tool.is_file():
        raise SealError("/usr/bin/otool is unavailable")
    output = {}
    for path in sorted(paths):
        if path.suffix not in {".so", ".dylib"} and path != Path(sys.executable).resolve():
            continue
        completed = subprocess.run(
            [str(tool), "-L", str(path)],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
        if completed.returncode != 0:
            raise SealError("otool failed for %s: %s" % (path, completed.stderr.strip()))
        output[str(path)] = completed.stdout.splitlines()[1:]
    return output

def extend_macho_dependencies(records):
    while True:
        dependency_map = otool_dependencies(records)
        added = False
        for loader_text, lines in dependency_map.items():
            loader = Path(loader_text)
            for line in lines:
                dependency = line.strip().split(" (", 1)[0]
                candidate = None
                if dependency.startswith("@loader_path/"):
                    candidate = loader.parent / dependency[len("@loader_path/"):]
                elif dependency.startswith("/"):
                    candidate = Path(dependency)
                if candidate is None or not candidate.exists():
                    continue
                resolved = candidate.resolve(strict=True)
                if str(resolved).startswith(("/usr/lib/", "/System/")):
                    continue
                if resolved not in records:
                    records[resolved] = {"resolved-macho-dependency"}
                    added = True
                else:
                    records[resolved].add("resolved-macho-dependency")
        if not added:
            return dependency_map


def object_hash(value):
    reduced = {key: item for key, item in value.items() if key != "report_sha256"}
    return hashlib.sha256(canonical_bytes(reduced)).hexdigest()


def seal(verifier_path):
    verifier_path = Path(verifier_path).resolve(strict=True)
    verifier, before = load_verifier(verifier_path)
    records = module_files(before)
    records.setdefault(verifier_path, set()).add("independent-verifier")
    flint_version = distribution_files(records)

    executable_link = Path(sys.executable)
    executable_target = executable_link.resolve(strict=True)
    records.setdefault(executable_target, set()).add("cpython-executable")
    source = Path(__file__).resolve(strict=True)
    records.setdefault(source, set()).add("environment-sealer")
    dependency_map = extend_macho_dependencies(records)

    file_records = [file_record(path, origins) for path, origins in records.items()]
    file_records.sort(key=lambda item: item["path"])
    result = {
        "schema": SCHEMA,
        "claim_status": "MACHINE-VERIFIED",
        "outcome": "PASS",
        "created_utc": _datetime.datetime.now(_datetime.timezone.utc).isoformat(),
        "runtime": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "executable_link": str(executable_link),
            "executable_target": str(executable_target),
            "prefix": sys.prefix,
            "base_prefix": sys.base_prefix,
            "python_flint": flint_version,
            "arb_context_precision_after_cover": verifier.ctx.prec,
        },
        "thread_environment": {
            name: os.environ.get(name) for name in THREAD_VARIABLES
        },
        "independent_verifier": {
            "path": str(verifier_path),
            "sha256": sha256_file(verifier_path),
            "size": verifier_path.stat().st_size,
        },
        "files": file_records,
        "file_count": len(file_records),
        "otool_dependencies": dependency_map,
        "sealer": {
            "path": str(source),
            "sha256": sha256_file(source),
            "size": source.stat().st_size,
        },
        "limitations": [
            "This manifest records observed bytes; it is not a sandbox or authenticated execution attestation.",
            "System libraries supplied through the macOS dyld shared cache, the kernel, firmware, and hardware are identified by platform metadata rather than copied or hashed here.",
            "A malicious interpreter or kernel could forge the manifest.",
            "Reproduction requires restoring every listed file or a byte-identical environment and then rerunning the verifier.",
        ],
    }
    result["report_sha256"] = object_hash(result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description="Seal the UC verifier runtime bytes")
    parser.add_argument("--verifier", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    try:
        result = seal(args.verifier)
        destination = Path(args.output)
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("xb") as stream:
            stream.write(canonical_bytes(result) + b"\n")
    except (SealError, OSError, subprocess.SubprocessError, ValueError) as exc:
        print("runtime seal FAIL: %s" % exc, file=sys.stderr)
        return 1
    print(
        canonical_bytes(
            {
                "claim_status": result["claim_status"],
                "file_count": result["file_count"],
                "outcome": result["outcome"],
                "report": str(destination),
                "report_sha256": result["report_sha256"],
            }
        ).decode("ascii")
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
