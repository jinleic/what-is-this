#!/usr/bin/env python3
"""Independently replay the committed source-zero VeriPB/CakePB certificate."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parents[1]
DEFAULT_MANIFEST = ROOT / "data" / "involution_f5_w_dfs_certificate_source0.json"
DEFAULT_VERIPB = WORKSPACE / "scratch/r55-cert-tools/veripb/target/release/veripb"
DEFAULT_CAKEPB = WORKSPACE / "scratch/r55-cert-tools/cakepb/cake_pb"
DEFAULT_WORK = WORKSPACE / "scratch/r55-involution-f5/dfs-proof-check"
VERIPB_SHA256 = "2744ee038573246c3b88cddabb03e8986a8ce1f2b1a7bb1167faaef3d0fa3860"
CAKEPB_SHA256 = "1938a5685d8cb8c66d50eef67d186beb97299e08d7be75e3565e8b3e3b1ead76"


class CheckViolation(RuntimeError):
    """The manifest, proof payload, or pinned checker failed."""


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise CheckViolation(f"duplicate JSON member {key!r}")
        result[key] = value
    return result


def load_document(path: Path) -> dict:
    try:
        document = json.loads(
            path.read_text(), object_pairs_hook=_unique_object,
            parse_constant=lambda token: (_ for _ in ()).throw(
                CheckViolation(f"nonstandard JSON constant {token}")))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise CheckViolation(f"cannot parse manifest: {error}") from error
    if type(document) is not dict:
        raise CheckViolation("manifest root must be an object")
    return document


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def _require_file(path: Path, size: int, digest: str, label: str) -> None:
    if not path.is_file() or path.stat().st_size != size:
        raise CheckViolation(f"{label} size/path mismatch")
    if _sha256(path) != digest:
        raise CheckViolation(f"{label} SHA-256 mismatch")


def _proof_policy(path: Path) -> None:
    conclusions = []
    with path.open() as stream:
        for line_number, line in enumerate(stream, 1):
            stripped = line.lstrip()
            if not stripped or stripped.startswith("*"):
                continue
            command = stripped.split(None, 1)[0].rstrip(";")
            if command in {"del", "d", "deld", "wiplvl"}:
                raise CheckViolation(
                    f"unchecked deletion {command!r} at line {line_number}")
            tokens = stripped.replace(";", " ;").split()
            if tokens and tokens[0] == "conclusion":
                conclusions.append(tokens[1] if len(tokens) > 1 else "")
    if conclusions != ["UNSAT"]:
        raise CheckViolation(f"proof conclusions {conclusions!r} != ['UNSAT']")


def _extract(record: dict, destination: Path, label: str) -> None:
    relative = record.get("relative_path")
    if type(relative) is not str or not relative.startswith("r55/data/"):
        raise CheckViolation(f"invalid {label} relative path")
    archive = ROOT / relative.removeprefix("r55/")
    _require_file(
        archive, record.get("gzip_bytes"), record.get("gzip_sha256"),
        f"{label} gzip")
    with gzip.open(archive, "rb") as source, destination.open("wb") as target:
        shutil.copyfileobj(source, target, length=1 << 20)
    _require_file(
        destination, record.get("uncompressed_bytes"),
        record.get("uncompressed_sha256"), label)


def verify(manifest_path: Path, veripb: Path, cakepb: Path,
           work: Path, run_checkers: bool = True) -> dict:
    document = load_document(manifest_path)
    if (document.get("schema_version") != 1
            or document.get("source_index") != 0
            or document.get("disposition") !=
            "F5_W_DFS_CERTIFICATE_SOURCE_0_VERIFIED"):
        raise CheckViolation("unexpected certificate identity")
    if document.get("coverage", {}).get("campaign_complete") is not False:
        raise CheckViolation("pilot manifest must not claim full coverage")
    generator = document.get("generator")
    if (type(generator) is not dict
            or generator.get("relative_source_path") !=
            "r55/src/involution_f5_w_dfs_certificate.cpp"):
        raise CheckViolation("invalid DFS generator provenance")
    generator_source = ROOT / "src" / "involution_f5_w_dfs_certificate.cpp"
    generator_binary = ROOT / "src" / "involution_f5_w_dfs_certificate_c"
    if (_sha256(generator_source) != generator.get("source_sha256")
            or _sha256(generator_binary) != generator.get("binary_sha256")):
        raise CheckViolation("DFS generator source/binary hash mismatch")
    work.mkdir(parents=True, exist_ok=True)
    formula = work / "source0.opb"
    proof = work / "source0.pbp"
    kernel = work / "source0.kernel.pbp"
    artifacts = document.get("artifacts")
    if type(artifacts) is not dict:
        raise CheckViolation("missing proof artifacts")
    _extract(artifacts.get("formula", {}), formula, "formula")
    _extract(artifacts.get("veripb_proof", {}), proof, "VeriPB proof")
    _extract(artifacts.get("cakepb_kernel_proof", {}), kernel, "CakePB proof")
    _proof_policy(proof)
    _proof_policy(kernel)
    for path, digest, label in (
            (veripb, VERIPB_SHA256, "VeriPB"),
            (cakepb, CAKEPB_SHA256, "CakePB")):
        if not path.is_file():
            raise CheckViolation(f"missing pinned {label} binary")
        _require_file(path, path.stat().st_size, digest, label)
    if run_checkers:
        verified = subprocess.run(
            [str(veripb), "--force-checked-deletion", str(formula), str(proof)],
            check=False, capture_output=True, text=True)
        if (verified.returncode != 0
                or "s VERIFIED UNSATISFIABLE" not in verified.stdout.splitlines()):
            raise CheckViolation(
                f"VeriPB rejected certificate: {verified.stdout} {verified.stderr}")
        environment = os.environ.copy()
        environment.update({
            "CML_HEAP_SIZE": "65536",
            "CML_STACK_SIZE": "16384",
        })
        checked = subprocess.run(
            [str(cakepb), str(formula), str(kernel)],
            check=False, capture_output=True, text=True, env=environment)
        if (checked.returncode != 0
                or "s VERIFIED UNSATISFIABLE" not in checked.stdout.splitlines()):
            raise CheckViolation(
                f"CakePB rejected certificate: {checked.stdout} {checked.stderr}")
    return {
        "source_index": 0,
        "veripb_unsat": True,
        "cakepb_unsat": True,
        "unchecked_deletion_present": False,
        "campaign_complete": False,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", nargs="?", type=Path,
                        default=DEFAULT_MANIFEST)
    parser.add_argument("--veripb", type=Path, default=DEFAULT_VERIPB)
    parser.add_argument("--cakepb", type=Path, default=DEFAULT_CAKEPB)
    parser.add_argument("--work", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--no-run-checkers", action="store_true")
    args = parser.parse_args(argv)
    result = verify(
        args.manifest, args.veripb, args.cakepb, args.work,
        run_checkers=not args.no_run_checkers)
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
