#!/usr/bin/env python3
"""Generate and check a W-square UNSAT certificate with no unchecked deletion."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import subprocess
import shutil
import sys
import time
from pathlib import Path

import involution_f5_signed_completion as completion

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parents[1]
DEFAULT_ROUNDINGSAT = (
    WORKSPACE / "scratch/r55-cert-tools/roundingsat/"
    "build-nodelete/roundingsat"
)
DEFAULT_VERIPB = (
    WORKSPACE / "scratch/r55-cert-tools/veripb/target/release/veripb"
)
DEFAULT_CAKEPB = WORKSPACE / "scratch/r55-cert-tools/cakepb/cake_pb"
ROUNDINGSAT_LOGGER_SOURCE = (
    WORKSPACE / "scratch/r55-cert-tools/roundingsat/src/Logger.cpp"
)
DEFAULT_WORK = WORKSPACE / "scratch/r55-involution-f5/certificates"
PINNED_ROUNDINGSAT_SHA256 = (
    "f18950b389989edef5b8e3d80468aa34c22b5f3f7a6d1ac678f70fdd6f6ee7bf"
)
PINNED_ROUNDINGSAT_LOGGER_SHA256 = (
    "e907f35cba70994e2225b780794ee9777b2ce6df9112a89dd2b2023f202d0b34"
)
PINNED_VERIPB_SHA256 = (
    "2744ee038573246c3b88cddabb03e8986a8ce1f2b1a7bb1167faaef3d0fa3860"
)
PINNED_CAKEPB_SHA256 = (
    "1938a5685d8cb8c66d50eef67d186beb97299e08d7be75e3565e8b3e3b1ead76"
)


class CertificateViolation(RuntimeError):
    """The solver, proof, pinned checker, or certificate contract failed."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def _file_record(path: Path) -> dict:
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": _sha256(path),
    }
def _require_pinned_file(path: Path, expected_sha256: str, label: str) -> None:
    if not path.is_file():
        raise CertificateViolation(f"missing {label} file {path}")
    observed = _sha256(path)
    if observed != expected_sha256:
        raise CertificateViolation(
            f"{label} SHA-256 {observed} != pinned {expected_sha256}")




def _run(command: list[str], *, timeout: float | None) -> tuple[subprocess.CompletedProcess, float]:
    started = time.monotonic()
    process = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    return process, time.monotonic() - started


def _assert_checked_deletion_policy(
        path: Path, *, allow_checked_deletion: bool) -> None:
    """Reject unchecked deletion grammar, independent of whitespace."""
    unchecked_commands = {"del", "d", "deld", "wiplvl"}
    with path.open("rt", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            stripped = line.lstrip()
            if not stripped or stripped.startswith("*"):
                continue
            command = stripped.split(None, 1)[0].rstrip(";")
            if command in unchecked_commands:
                raise CertificateViolation(
                    f"proof contains unchecked deletion command {command!r} "
                    f"at line {line_number}")
            if command == "delc" and not allow_checked_deletion:
                raise CertificateViolation(
                    f"proof unexpectedly contains checked deletion at line "
                    f"{line_number}")


def _require_unsat_conclusion(path: Path, label: str) -> None:
    conclusions = []
    with path.open("rt", encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if not stripped or stripped.startswith("*"):
                continue
            tokens = stripped.replace(";", " ;").split()
            if tokens and tokens[0] == "conclusion":
                conclusions.append(tokens[1] if len(tokens) > 1 else "")
    if conclusions != ["UNSAT"]:
        raise CertificateViolation(
            f"{label} conclusions {conclusions!r} != ['UNSAT']")


def certify(source_index: int, *, roundingsat: Path, veripb: Path,
            cakepb: Path, work: Path, timeout: float | None,
            core_rows: tuple[int, ...] | None = None) -> dict:
    for label, path, digest in (
            ("RoundingSat", roundingsat, PINNED_ROUNDINGSAT_SHA256),
            ("VeriPB", veripb, PINNED_VERIPB_SHA256),
            ("CakePB", cakepb, PINNED_CAKEPB_SHA256),
            ("RoundingSat deletion-free Logger.cpp",
             ROUNDINGSAT_LOGGER_SOURCE,
             PINNED_ROUNDINGSAT_LOGGER_SHA256)):
        _require_pinned_file(path, digest, label)
    work.mkdir(parents=True, exist_ok=True)
    instance = completion.load_instance(source_index)
    rows = tuple(range(20)) if core_rows is None else core_rows
    if (type(rows) is not tuple or not rows or len(set(rows)) != len(rows)
            or any(type(row) is not int or not 0 <= row < 20
                   for row in rows)):
        raise CertificateViolation("invalid W certificate core rows")
    row_set = set(rows)
    for permutation in completion.residual_generators(instance):
        if {permutation[row] for row in row_set} != row_set:
            raise CertificateViolation(
                "W certificate core is not residual-symmetry closed")
    formula, registry = completion.build_signed_formula(
        instance,
        include_t_square=False,
        w_square_rows=rows,
        include_w_domain_tables=True,
        w_domain_table_rows=rows,
        include_ramsey=False,
    )
    stem = f"w-square-{source_index:03d}"
    # Formula and raw RoundingSat proof are verified before the next sequential
    # candidate overwrites them.  Per-candidate CakePB kernel proofs persist.
    opb = work / "current.opb"
    proof = work / "current.pbp"
    kernel = work / "current.kernel.pbp"
    kernel_archive = work / f"{stem}.kernel.pbp.gz"
    manifest = work / f"{stem}.json"
    opb.write_text(formula.to_opb(registry))

    solver, solver_wall = _run([
        str(roundingsat),
        "--verbosity=0",
        f"--proof-log={proof}",
        str(opb),
    ], timeout=timeout)
    if solver.returncode != 20 or "s UNSATISFIABLE" not in solver.stdout:
        raise CertificateViolation(
            f"RoundingSat did not prove UNSAT (exit {solver.returncode}): "
            f"{solver.stdout[-1000:]} {solver.stderr[-1000:]}")
    _assert_checked_deletion_policy(
        proof, allow_checked_deletion=False)
    _require_unsat_conclusion(proof, "RoundingSat proof")

    verifier, verifier_wall = _run([
        str(veripb),
        "--force-checked-deletion",
        "--elaborate", str(kernel),
        str(opb), str(proof),
    ], timeout=timeout)
    verifier_lines = verifier.stdout.splitlines()
    if (verifier.returncode != 0
            or "s VERIFIED UNSATISFIABLE" not in verifier_lines):
        raise CertificateViolation(
            f"VeriPB did not verify UNSAT: {verifier.stdout[-1000:]} "
            f"{verifier.stderr[-1000:]}")
    _assert_checked_deletion_policy(
        kernel, allow_checked_deletion=True)
    _require_unsat_conclusion(kernel, "elaborated kernel proof")

    cake, cake_wall = _run(
        [str(cakepb), str(opb), str(kernel)], timeout=timeout)
    if cake.returncode != 0:
        raise CertificateViolation(
            f"CakePB rejected the kernel proof: {cake.stdout[-1000:]} "
            f"{cake.stderr[-1000:]}")
    kernel_record = _file_record(kernel)
    with kernel.open("rb") as source, gzip.GzipFile(
            filename=str(kernel_archive), mode="wb",
            compresslevel=1, mtime=0) as target:
        shutil.copyfileobj(source, target, length=1 << 20)
    record = {
        "schema_version": 1,
        "campaign_id": "involution_f5_w_square_certificate",
        "source_index": source_index,
        "claim": "W_SQUARE_UNSAT_FOR_THIS_RAMSEY_SUPPORT_REPRESENTATIVE",
        "formula": {
            **_file_record(opb),
            "variables": registry.variable_count,
            "constraints": formula.constraint_count,
            "direct_w_square_row_core": list(rows),
            "row_domain_extension_used": True,
            "symmetry_breaking_used": False,
        },
        "roundingsat": {
            "base_commit": completion.ROUNDINGSAT_COMMIT,
            "binary": _file_record(roundingsat),
            "deletion_free_logger_source":
                _file_record(ROUNDINGSAT_LOGGER_SOURCE),
            "deletion_free_logger_patch": True,
            "proof": _file_record(proof),
            "wall_seconds": round(solver_wall, 6),
            "unchecked_deletion_present": False,
        },
        "veripb": {
            "commit": completion.VERIPB_COMMIT,
            "binary": _file_record(veripb),
            "wall_seconds": round(verifier_wall, 6),
            "force_checked_deletion": True,
            "verified_conclusion": "UNSAT",
            "checked_kernel_deletion_allowed": True,
            "accepted": True,
        },
        "cakepb": {
            "commit": completion.CAKEPB_COMMIT,
            "binary": _file_record(cakepb),
            "kernel_proof_uncompressed": kernel_record,
            "kernel_proof_gzip": _file_record(kernel_archive),
            "wall_seconds": round(cake_wall, 6),
            "accepted": True,
            "verified_kernel_conclusion": "UNSAT",
        },
        "projection_or_orbit_count_used": False,
        "general_ramsey_bound_claimed": False,
    }
    manifest.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-index", type=int, required=True)
    parser.add_argument("--roundingsat", type=Path, default=DEFAULT_ROUNDINGSAT)
    parser.add_argument("--veripb", type=Path, default=DEFAULT_VERIPB)
    parser.add_argument("--cakepb", type=Path, default=DEFAULT_CAKEPB)
    parser.add_argument("--work", type=Path, default=DEFAULT_WORK)
    parser.add_argument(
        "--core-rows",
        help="comma-separated residual-symmetry-closed W row core")
    parser.add_argument("--timeout", type=float)
    args = parser.parse_args(argv)
    core_rows = (
        tuple(int(value) for value in args.core_rows.split(","))
        if args.core_rows else None
    )
    record = certify(
        args.source_index,
        roundingsat=args.roundingsat,
        veripb=args.veripb,
        cakepb=args.cakepb,
        work=args.work,
        timeout=args.timeout,
        core_rows=core_rows,
    )
    print(json.dumps({
        "source_index": record["source_index"],
        "claim": record["claim"],
        "proof_bytes": record["roundingsat"]["proof"]["bytes"],
        "kernel_proof_bytes":
            record["cakepb"]["kernel_proof_uncompressed"]["bytes"],
        "kernel_proof_gzip_bytes":
            record["cakepb"]["kernel_proof_gzip"]["bytes"],
        "veripb_accepted": record["veripb"]["accepted"],
        "cakepb_accepted": record["cakepb"]["accepted"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
