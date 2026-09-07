#!/usr/bin/env python3
"""Bind the proof-free n=12 target-39 campaign to exact CNFs, solver,
exec wrapper, and supervisor implementation."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
CNF_DIR = WORKSPACE / "scratch" / "kobon" / "n12"
ENGINE = WORKSPACE / "math" / "kobon" / "engine.py"
RUNNER = HERE / "run_cubes.sh"
IMPL = HERE / "supervise.py"
OUTPUT = HERE / "manifest.json"
CUBES = (
    "simple", "par_noc", "par_conc",
    "c0", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9",
    "c10", "c11",
)


class ManifestFailure(RuntimeError):
    """A required campaign input is absent or malformed."""


def _atomic_write_text(path: Path, text: str) -> None:
    """Whole-file JSON write: sibling temp + fsync + replace + dir fsync."""
    tmp = path.with_name(f"{path.name}.tmp-{os.getpid()}")
    with tmp.open("w", encoding="utf-8") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)
    fd = os.open(path.parent, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def _manifest_path(path: Path) -> str:
    try:
        return str(path.relative_to(WORKSPACE))
    except ValueError:
        return str(path)


def _record(path: Path) -> dict:
    if not path.is_file():
        raise ManifestFailure(f"missing file: {path}")
    return {
        "path": _manifest_path(path),
        "bytes": path.stat().st_size,
        "sha256": _sha256(path),
    }


def _dimacs(path: Path) -> tuple[int, int]:
    with path.open("rb") as stream:
        for line in stream:
            if line.startswith(b"c"):
                continue
            fields = line.split()
            if len(fields) == 4 and fields[:2] == [b"p", b"cnf"]:
                return int(fields[2]), int(fields[3])
            raise ManifestFailure(f"invalid DIMACS header in {path}")
    raise ManifestFailure(f"missing DIMACS header in {path}")


def build() -> dict:
    solver_text = shutil.which("kissat")
    if solver_text is None:
        raise ManifestFailure("kissat is not on PATH")
    solver = Path(solver_text).resolve()
    version = subprocess.run(
        [str(solver), "--version"], check=True,
        capture_output=True, text=True).stdout.strip()
    instances = []
    for cube in CUBES:
        path = CNF_DIR / f"kobon_n12_t39_proof_{cube}.cnf"
        record = _record(path)
        variables, clauses = _dimacs(path)
        record.update({
            "cube": cube,
            "variables": variables,
            "clauses": clauses,
        })
        instances.append(record)
    return {
        "schema": "kobon-n12-decision-manifest/2",
        "claim_scope": "K_gen(12), broad convention, relaxed target-39 15-case split",
        "case_count": len(instances),
        "case_order": list(CUBES),
        "evidence_contract": {
            "SAT": "Boolean candidate only; requires extraction, straightening, and exact-rational verify_selection(minimum=39)",
            "all_UNSAT_without_proofs": "discovery evidence only",
            "manifest_binding": "Hashes the current engine and existing CNFs independently; does not prove the engine emitted those CNFs. Byte-regeneration/provenance audit is required before certification",
            "theorem_grade_upper_bound": "all cases require complete DRAT, drat-trim s VERIFIED, and CNF provenance/cover audit",
        },
        "solver": {**_record(solver), "version": version},
        "engine": _record(ENGINE),
        "runner": _record(RUNNER),
        "implementation": _record(IMPL),
        "instances": instances,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    text = json.dumps(build(), indent=2, sort_keys=True) + "\n"
    if args.check:
        if not OUTPUT.is_file() or OUTPUT.read_text() != text:
            raise ManifestFailure("manifest is absent or stale")
    else:
        _atomic_write_text(OUTPUT, text)
    print(json.dumps({
        "case_count": len(CUBES),
        "manifest_sha256": hashlib.sha256(text.encode()).hexdigest(),
        "status": "OK",
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
