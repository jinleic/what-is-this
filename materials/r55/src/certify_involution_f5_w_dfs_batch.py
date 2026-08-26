#!/usr/bin/env python3
"""Sequentially certify all 550 W-negative support representatives."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import involution_f5_signed_completion as completion
import involution_f5_w_square_census as w_census

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parents[1]
W_CENSUS = ROOT / "data" / "involution_f5_w_square_census.json"
SOURCE0_MANIFEST = ROOT / "data" / "involution_f5_w_dfs_certificate_source0.json"
DEFAULT_OUTPUT = ROOT / "data" / "involution_f5_w_dfs_certificates.json"
DEFAULT_PROOF_DIR = ROOT / "data" / "involution_f5_w_dfs_proofs"
DEFAULT_WORK = WORKSPACE / "scratch/r55-involution-f5/w-dfs-batch"
GENERATOR_SOURCE = ROOT / "src" / "involution_f5_w_dfs_certificate.cpp"
GENERATOR = ROOT / "src" / "involution_f5_w_dfs_certificate_c"
VERIPB = WORKSPACE / "scratch/r55-cert-tools/veripb/target/release/veripb"
CAKEPB = WORKSPACE / "scratch/r55-cert-tools/cakepb/cake_pb"
W_CENSUS_SHA256 = "e99963ed0494b90a7cd3cb117113b317132b711e0592867d27353e76e8a11a6a"
GENERATOR_SOURCE_SHA256 = "565cc32f1925d813e094ee862ba24b2da2004ae889e9f7c3b0202b70eb6fb0ca"
GENERATOR_SHA256 = "1eef7e4274611a53b034f78f0ebb588f4f65d9d9b4e82c348d773650f10dad76"
VERIPB_SHA256 = "2744ee038573246c3b88cddabb03e8986a8ce1f2b1a7bb1167faaef3d0fa3860"
CAKEPB_SHA256 = "1938a5685d8cb8c66d50eef67d186beb97299e08d7be75e3565e8b3e3b1ead76"
ROWS = tuple(range(20))


class BatchViolation(RuntimeError):
    """A proof generator, checker, pin, or persisted record failed."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def _record(path: Path) -> dict:
    return {"bytes": path.stat().st_size, "sha256": _sha256(path)}


def _require_pin(path: Path, digest: str, label: str) -> None:
    if not path.is_file() or _sha256(path) != digest:
        raise BatchViolation(f"{label} pin mismatch")


def _proof_policy(path: Path) -> None:
    conclusions = []
    with path.open() as stream:
        for line_number, line in enumerate(stream, 1):
            stripped = line.lstrip()
            if not stripped or stripped.startswith("*"):
                continue
            command = stripped.split(None, 1)[0].rstrip(";")
            if command in {"del", "d", "deld", "wiplvl"}:
                raise BatchViolation(
                    f"unchecked deletion {command!r} at {path}:{line_number}")
            tokens = stripped.replace(";", " ;").split()
            if tokens and tokens[0] == "conclusion":
                conclusions.append(tokens[1] if len(tokens) > 1 else "")
    if conclusions != ["UNSAT"]:
        raise BatchViolation(f"{path} does not conclude UNSAT exactly once")


def _run(command: list[str], timeout: float, env=None) -> tuple[subprocess.CompletedProcess, float]:
    started = time.monotonic()
    process = subprocess.run(
        command, check=False, capture_output=True, text=True,
        timeout=timeout, env=env)
    return process, time.monotonic() - started


def _gzip(source: Path, destination: Path) -> dict:
    with source.open("rb") as input_stream, gzip.GzipFile(
            filename=str(destination), mode="wb", compresslevel=1,
            mtime=0) as output_stream:
        shutil.copyfileobj(input_stream, output_stream, length=1 << 20)
    return _record(destination)


def _seed_source_zero() -> dict:
    source = json.loads(SOURCE0_MANIFEST.read_text())
    if (source.get("source_index") != 0
            or source.get("checkers", {}).get("veripb", {}).get(
                "verified_conclusion") != "UNSAT"
            or source.get("checkers", {}).get("cakepb", {}).get(
                "verified_conclusion") != "UNSAT"):
        raise BatchViolation("source-zero certificate seed is invalid")
    return {
        "source_index": 0,
        "status": "VERIPB_CAKEPB_VERIFIED",
        "seed_manifest": "r55/data/involution_f5_w_dfs_certificate_source0.json",
        "formula": source["artifacts"]["formula"],
        "veripb_proof": source["artifacts"]["veripb_proof"],
        "cakepb_kernel_proof": source["artifacts"]["cakepb_kernel_proof"],
        "generator": source["generator"],
        "timing": {
            "veripb_verify_seconds":
                source["checkers"]["veripb"]["verification_wall_seconds"],
            "veripb_elaborate_seconds":
                source["checkers"]["veripb"]["elaboration_wall_seconds"],
            "cakepb_seconds":
                source["checkers"]["cakepb"]["verification_wall_seconds"],
        },
    }


def _write_progress(output: Path, records: list[dict], target: int,
                    started: float) -> None:
    document = {
        "schema_version": 1,
        "campaign_id": "involution_f5_w_dfs_certificates",
        "disposition": (
            "F5_W_NEGATIVE_550_VERIPB_CAKEPB_CERTIFIED"
            if len(records) == target else
            f"F5_W_CERTIFICATE_BATCH_IN_PROGRESS_{len(records)}_OF_{target}"
        ),
        "input_dependency": {
            "relative_path": "r55/data/involution_f5_w_square_census.json",
            "bytes": W_CENSUS.stat().st_size,
            "sha256": W_CENSUS_SHA256,
        },
        "toolchain": {
            "generator_source_sha256": GENERATOR_SOURCE_SHA256,
            "generator_binary_sha256": GENERATOR_SHA256,
            "veripb_commit": completion.VERIPB_COMMIT,
            "veripb_binary_sha256": VERIPB_SHA256,
            "cakepb_commit": completion.CAKEPB_COMMIT,
            "cakepb_binary_sha256": CAKEPB_SHA256,
            "unchecked_deletion_allowed": False,
            "cakepb_memory_profiles_mib": [
                {"heap": 65536, "stack": 4096},
                {"heap": 90112, "stack": 4096}
            ],
        },
        "coverage": {
            "target_w_negative_support_orbits": target,
            "certified_support_orbits": len(records),
            "remaining_support_orbits": target - len(records),
            "campaign_complete": len(records) == target,
        },
        "batch_wall_seconds_so_far": round(time.monotonic() - started, 6),
        "records": records,
        "strongly_regular_fixed_five_branch_only": True,
        "general_ramsey_bound_claimed": False,
    }
    output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")


def run(output: Path, proof_dir: Path, work: Path,
        timeout_per_step: float) -> dict:
    for path, digest, label in (
            (W_CENSUS, W_CENSUS_SHA256, "W census"),
            (GENERATOR_SOURCE, GENERATOR_SOURCE_SHA256, "DFS source"),
            (GENERATOR, GENERATOR_SHA256, "DFS binary"),
            (VERIPB, VERIPB_SHA256, "VeriPB"),
            (CAKEPB, CAKEPB_SHA256, "CakePB")):
        _require_pin(path, digest, label)
    census = json.loads(W_CENSUS.read_text())
    targets = [record for record in census["census"]["candidates"]
               if not record["w_square_satisfiable"]]
    if len(targets) != 550 or targets[0]["source_index"] != 0:
        raise BatchViolation("unexpected W-negative target set")
    proof_dir.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)
    records = [_seed_source_zero()]
    if output.is_file():
        previous = json.loads(output.read_text())
        prior_records = previous.get("records")
        if type(prior_records) is list and prior_records:
            records = prior_records
    completed = {record["source_index"] for record in records}
    started = time.monotonic()
    _write_progress(output, records, len(targets), started)
    dump = work / "current.dom"
    formula = work / "current.opb"
    proof = work / "current.pbp"
    kernel = work / "current.kernel.pbp"
    rows = ",".join(map(str, ROWS))
    for ordinal, target in enumerate(targets):
        source_index = target["source_index"]
        if source_index in completed:
            continue
        instance = completion.load_instance(source_index)
        dump_record = w_census.write_domain_dump(instance, dump)
        generated, generator_wall = _run(
            [str(GENERATOR), str(dump), rows, str(formula), str(proof)],
            timeout_per_step)
        if generated.returncode != 0:
            raise BatchViolation(
                f"DFS generator failed for {source_index}: "
                f"{generated.stdout} {generated.stderr}")
        generator_stats = json.loads(generated.stdout.strip())
        _proof_policy(proof)
        elaborated, elaborate_wall = _run(
            [str(VERIPB), "--force-checked-deletion",
             "--elaborate", str(kernel), str(formula), str(proof)],
            timeout_per_step)
        if (elaborated.returncode != 0
                or "s VERIFIED UNSATISFIABLE" not in
                elaborated.stdout.splitlines()):
            raise BatchViolation(
                f"VeriPB elaboration failed for {source_index}")
        _proof_policy(kernel)
        cake_wall = 0.0
        cake_profile = None
        checked = None
        memory_profiles = (
            (90112,) if formula.stat().st_size > 250_000_000
            else (65536, 90112)
        )
        for heap_mib in memory_profiles:
            cake_environment = os.environ.copy()
            cake_environment.update({
                "CML_HEAP_SIZE": str(heap_mib),
                "CML_STACK_SIZE": "4096",
            })
            attempt, attempt_wall = _run(
                [str(CAKEPB), str(formula), str(kernel)],
                timeout_per_step, env=cake_environment)
            cake_wall += attempt_wall
            if (attempt.returncode == 0
                    and "s VERIFIED UNSATISFIABLE" in
                    attempt.stdout.splitlines()):
                checked = attempt
                cake_profile = {"heap_mib": heap_mib, "stack_mib": 4096}
                break
            if "CakeML heap space exhausted" not in attempt.stderr:
                raise BatchViolation(
                    f"CakePB rejected source {source_index}: "
                    f"{attempt.stdout} {attempt.stderr}")
        if checked is None:
            raise BatchViolation(
                f"CakePB exhausted every memory profile for {source_index}")
        proof_archive = proof_dir / f"w{source_index:03d}.pbp.gz"
        kernel_archive = proof_dir / f"w{source_index:03d}.kernel.pbp.gz"
        record = {
            "source_index": source_index,
            "ordinal_in_w_negative_set": ordinal,
            "status": "VERIPB_CAKEPB_VERIFIED",
            "row_domains_sha256": dump_record["row_domains_sha256"],
            "generator": generator_stats,
            "formula": _record(formula),
            "veripb_proof": {
                **_record(proof),
                "gzip": _gzip(proof, proof_archive),
                "relative_gzip_path":
                    f"r55/data/involution_f5_w_dfs_proofs/{proof_archive.name}",
            },
            "cakepb_kernel_proof": {
                **_record(kernel),
                "gzip": _gzip(kernel, kernel_archive),
                "relative_gzip_path":
                    f"r55/data/involution_f5_w_dfs_proofs/{kernel_archive.name}",
            },
            "timing": {
                "generator_seconds": round(generator_wall, 6),
                "veripb_elaborate_seconds": round(elaborate_wall, 6),
                "cakepb_seconds": round(cake_wall, 6),
                "cakepb_memory_profile": cake_profile,
            },
            "unchecked_deletion_present": False,
        }
        records.append(record)
        completed.add(source_index)
        _write_progress(output, records, len(targets), started)
        print(json.dumps({
            "certified": len(records),
            "target": len(targets),
            "source_index": source_index,
            "cakepb_seconds": round(cake_wall, 6),
        }, sort_keys=True), file=sys.stderr, flush=True)
    return json.loads(output.read_text())


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--proof-dir", type=Path, default=DEFAULT_PROOF_DIR)
    parser.add_argument("--work", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--timeout-per-step", type=float, default=900.0)
    args = parser.parse_args(argv)
    document = run(
        args.output, args.proof_dir, args.work, args.timeout_per_step)
    print(json.dumps({
        "disposition": document["disposition"],
        "certified_support_orbits":
            document["coverage"]["certified_support_orbits"],
        "remaining_support_orbits":
            document["coverage"]["remaining_support_orbits"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
