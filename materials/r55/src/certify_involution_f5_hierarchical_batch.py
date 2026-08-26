#!/usr/bin/env python3
"""Certify the 155 deeper signed/Ramsey-negative support representatives."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import struct
import subprocess
import sys
import time
from pathlib import Path

import involution_f5_signed_completion as completion
import involution_f5_t_square_census as t_census
import involution_f5_w_square_census as w_census

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parents[1]
W_CENSUS = ROOT / "data" / "involution_f5_w_square_census.json"
SIGNED_CENSUS = ROOT / "data" / "involution_f5_signed_square_census.json"
RAMSEY_CENSUS = ROOT / "data" / "involution_f5_ramsey_square_census.json"
DEFAULT_OUTPUT = ROOT / "data" / "involution_f5_hierarchical_certificates.json"
DEFAULT_PROOF_DIR = ROOT / "data" / "involution_f5_hierarchical_proofs"
DEFAULT_WORK = WORKSPACE / "scratch/r55-involution-f5/hierarchical-batch"
W_TERMINAL_SOURCE = ROOT / "src" / "involution_f5_w_terminal_certificate.cpp"
W_TERMINAL = ROOT / "src" / "involution_f5_w_terminal_certificate_c"
T_DFS_SOURCE = ROOT / "src" / "involution_f5_t_dfs_certificate.cpp"
T_DFS = ROOT / "src" / "involution_f5_t_dfs_certificate_c"
T_RAMSEY_SOURCE = ROOT / "src" / "involution_f5_t_ramsey_certificate.cpp"
T_RAMSEY = ROOT / "src" / "involution_f5_t_ramsey_certificate_c"
VERIPB = WORKSPACE / "scratch/r55-cert-tools/veripb/target/release/veripb"
CAKEPB = WORKSPACE / "scratch/r55-cert-tools/cakepb/cake_pb"
PINS = {
    W_CENSUS: "e99963ed0494b90a7cd3cb117113b317132b711e0592867d27353e76e8a11a6a",
    SIGNED_CENSUS: "01a805f5a783d572966785cedb5c31c89b69ff986ed6bdc77bf3a5dc4c3b8150",
    RAMSEY_CENSUS: "6069740d3c5339dc1d734e10f70d07be53f7046ddcbf07c37974fab322f359c1",
    W_TERMINAL_SOURCE: "bea4993c06dd2c0274481cb21937751e2cff38b1d3511ad28542bd9923af4e1d",
    W_TERMINAL: "651342c43352fcb67d650221304eb83b6ec3863df023e8069fa8fa16e3c2202b",
    T_DFS_SOURCE: "66fff245e8222eabd4406643aef9c9543524423a5645eefb0e952bf5aa82bfe5",
    T_DFS: "940856b07861d999e0d8d2d871a6d50519271e3c533d4fb72597d2f0d3e974ec",
    T_RAMSEY_SOURCE: "6f2cee305c9ebd50a129894754291b603ff66471d7c83e70527002b234e364d6",
    T_RAMSEY: "8a8ff0a8fe223d710adc98366813d42c17be598a39ddc7c39891fda37c097b28",
    VERIPB: "2744ee038573246c3b88cddabb03e8986a8ce1f2b1a7bb1167faaef3d0fa3860",
    CAKEPB: "1938a5685d8cb8c66d50eef67d186beb97299e08d7be75e3565e8b3e3b1ead76",
}
ROWS = tuple(range(20))
ROWS_TEXT = ",".join(map(str, ROWS))


class HierarchyViolation(RuntimeError):
    """A modular proof, terminal binding, checker, or pin failed."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def _record(path: Path) -> dict:
    return {"bytes": path.stat().st_size, "sha256": _sha256(path)}


def _require_pins() -> None:
    for path, digest in PINS.items():
        if not path.is_file() or _sha256(path) != digest:
            raise HierarchyViolation(f"pin mismatch: {path}")


def _proof_policy(path: Path) -> None:
    conclusions = []
    with path.open() as stream:
        for line_number, line in enumerate(stream, 1):
            stripped = line.lstrip()
            if not stripped or stripped.startswith("*"):
                continue
            command = stripped.split(None, 1)[0].rstrip(";")
            if command in {"del", "d", "deld", "wiplvl"}:
                raise HierarchyViolation(
                    f"unchecked deletion at {path}:{line_number}")
            tokens = stripped.replace(";", " ;").split()
            if tokens and tokens[0] == "conclusion":
                conclusions.append(tokens[1] if len(tokens) > 1 else "")
    if conclusions != ["UNSAT"]:
        raise HierarchyViolation(f"non-UNSAT proof conclusion: {path}")


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


def _checked_bundle(formula: Path, proof: Path, kernel: Path,
                    timeout: float, large: bool) -> dict:
    _proof_policy(proof)
    elaborated, elaborate_wall = _run(
        [str(VERIPB), "--force-checked-deletion", "--elaborate", str(kernel),
         str(formula), str(proof)], timeout)
    if (elaborated.returncode != 0
            or "s VERIFIED UNSATISFIABLE" not in
            elaborated.stdout.splitlines()):
        raise HierarchyViolation("VeriPB elaboration failed")
    _proof_policy(kernel)
    profiles = (
        ((90112,) if formula.stat().st_size > 250_000_000
         else (65536, 90112))
        if large else (8192, 32768, 65536)
    )
    cake_wall = 0.0
    cake_profile = None
    for heap_mib in profiles:
        environment = os.environ.copy()
        environment.update({
            "CML_HEAP_SIZE": str(heap_mib),
            "CML_STACK_SIZE": "4096",
        })
        checked, attempt_wall = _run(
            [str(CAKEPB), str(formula), str(kernel)], timeout, env=environment)
        cake_wall += attempt_wall
        if (checked.returncode == 0
                and "s VERIFIED UNSATISFIABLE" in checked.stdout.splitlines()):
            cake_profile = {"heap_mib": heap_mib, "stack_mib": 4096}
            break
        if "CakeML heap space exhausted" not in checked.stderr:
            raise HierarchyViolation(
                f"CakePB rejected proof: {checked.stdout} {checked.stderr}")
    if cake_profile is None:
        raise HierarchyViolation("CakePB exhausted every memory profile")
    return {
        "formula": _record(formula),
        "proof": _record(proof),
        "kernel": _record(kernel),
        "veripb_elaborate_seconds": round(elaborate_wall, 6),
        "cakepb_seconds": round(cake_wall, 6),
        "cakepb_memory_profile": cake_profile,
        "unchecked_deletion_present": False,
    }


def _put_u32(stream, value: int) -> None:
    stream.write(struct.pack(">I", value))


def _read_u32(data: bytes, position: int) -> tuple[int, int]:
    if position + 4 > len(data):
        raise HierarchyViolation("truncated terminal sidecar")
    return struct.unpack(">I", data[position:position + 4])[0], position + 4


def _write_t_input(instance: completion.CompletionInstance, source_index: int,
                   w_trits: str, path: Path, ramsey: bool) -> None:
    with path.open("wb") as stream:
        stream.write(b"TRAM1" if ramsey else b"TDOM1")
        _put_u32(stream, source_index)
        if ramsey:
            stream.write(bytes(instance.masks))
        stream.write(bytes(instance.internal))
        stream.write(w_trits.encode("ascii"))


def _parse_w_terminals(path: Path, source_index: int) -> list[dict]:
    data = path.read_bytes()
    if data[:6] != b"WTERM1":
        raise HierarchyViolation("bad W terminal magic")
    position = 6
    observed, position = _read_u32(data, position)
    count, position = _read_u32(data, position)
    if observed != source_index:
        raise HierarchyViolation("W terminal source mismatch")
    terminals = []
    for terminal_index in range(count):
        if position + 190 > len(data):
            raise HierarchyViolation("truncated W terminal")
        w_trits = data[position:position + 190].decode("ascii")
        position += 190
        formula_id, position = _read_u32(data, position)
        terminals.append({
            "terminal_index": terminal_index,
            "w_half_trits": w_trits,
            "formula_id": formula_id,
        })
    if position != len(data):
        raise HierarchyViolation("trailing W terminal bytes")
    return terminals


def _parse_t_terminals(
        path: Path, instance: completion.CompletionInstance,
        w_trits: str, source_index: int) -> list[dict]:
    data = path.read_bytes()
    if data[:6] != b"TTERM1":
        raise HierarchyViolation("bad T terminal magic")
    position = 6
    observed, position = _read_u32(data, position)
    count, position = _read_u32(data, position)
    if observed != source_index:
        raise HierarchyViolation("T terminal source mismatch")
    terminals = []
    for terminal_index in range(count):
        if position + 190 > len(data):
            raise HierarchyViolation("truncated T terminal")
        t_trits = data[position:position + 190].decode("ascii")
        position += 190
        formula_id, position = _read_u32(data, position)
        if position + 6 > len(data):
            raise HierarchyViolation("truncated Ramsey witness")
        kind = chr(data[position])
        vertices = tuple(data[position + 1:position + 6])
        position += 6
        states = t_census.relation_states(w_trits, t_trits)
        graph = completion.verify_state_assignment(instance, states)
        if kind == "K":
            witness = graph["first_K5"]
        elif kind == "I":
            witness = graph["first_I5"]
        else:
            raise HierarchyViolation("invalid Ramsey witness kind")
        if witness is None or tuple(witness) != vertices:
            raise HierarchyViolation("Ramsey terminal witness mismatch")
        terminals.append({
            "terminal_index": terminal_index,
            "formula_id": formula_id,
            "kind": kind,
            "vertices": list(vertices),
        })
    if position != len(data):
        raise HierarchyViolation("trailing T terminal bytes")
    return terminals


def _archive_bundle(prefix: str, proof_dir: Path, proof: Path,
                    kernel: Path, checked: dict) -> dict:
    proof_archive = proof_dir / f"{prefix}.pbp.gz"
    kernel_archive = proof_dir / f"{prefix}.kernel.pbp.gz"
    return {
        **checked,
        "proof_gzip": _gzip(proof, proof_archive),
        "proof_relative_path":
            f"r55/data/involution_f5_hierarchical_proofs/{proof_archive.name}",
        "kernel_gzip": _gzip(kernel, kernel_archive),
        "kernel_relative_path":
            f"r55/data/involution_f5_hierarchical_proofs/{kernel_archive.name}",
    }


def _write_progress(output: Path, records: list[dict], target: int,
                    started: float) -> None:
    document = {
        "schema_version": 1,
        "campaign_id": "involution_f5_hierarchical_certificates",
        "disposition": (
            "F5_DEEPER_155_VERIPB_CAKEPB_CERTIFIED"
            if len(records) == target else
            f"F5_HIERARCHICAL_CERTIFICATES_IN_PROGRESS_{len(records)}_OF_{target}"
        ),
        "coverage": {
            "target_deeper_support_orbits": target,
            "certified_support_orbits": len(records),
            "remaining_support_orbits": target - len(records),
            "campaign_complete": len(records) == target,
        },
        "toolchain": {
            "pins": {str(path.relative_to(WORKSPACE)): digest
                     for path, digest in PINS.items()},
            "veripb_commit": completion.VERIPB_COMMIT,
            "cakepb_commit": completion.CAKEPB_COMMIT,
            "unchecked_deletion_allowed": False,
        },
        "batch_wall_seconds_so_far": round(time.monotonic() - started, 6),
        "records": records,
        "strongly_regular_fixed_five_branch_only": True,
        "general_ramsey_bound_claimed": False,
    }
    output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")


def run(output: Path, proof_dir: Path, work: Path,
        timeout_per_step: float) -> dict:
    _require_pins()
    signed_document = json.loads(SIGNED_CENSUS.read_text())
    ramsey_document = json.loads(RAMSEY_CENSUS.read_text())
    signed_records = {
        record["source_index"]: record
        for record in signed_document["census"]["candidates"]
    }
    ramsey_records = {
        record["source_index"]: record
        for record in ramsey_document["census"]["candidates"]
    }
    targets = sorted(signed_records)
    if len(targets) != 155 or len(ramsey_records) != 108:
        raise HierarchyViolation("unexpected deeper target set")
    proof_dir.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)
    records = []
    if output.is_file():
        previous = json.loads(output.read_text())
        if type(previous.get("records")) is list:
            records = previous["records"]
    completed = {record["source_index"] for record in records}
    started = time.monotonic()
    _write_progress(output, records, len(targets), started)
    dump = work / "current.wdom"
    w_formula = work / "current.w.opb"
    w_proof = work / "current.w.pbp"
    w_kernel = work / "current.w.kernel.pbp"
    w_terminals_path = work / "current.w.term"
    t_input = work / "current.tdom"
    t_formula = work / "current.t.opb"
    t_proof = work / "current.t.pbp"
    t_kernel = work / "current.t.kernel.pbp"
    t_terminals_path = work / "current.t.term"
    for target_ordinal, source_index in enumerate(targets):
        if source_index in completed:
            continue
        instance = completion.load_instance(source_index)
        ramsey_mode = source_index in ramsey_records
        dump_record = w_census.write_domain_dump(instance, dump)
        generated, w_generate_wall = _run(
            [str(W_TERMINAL), str(dump), ROWS_TEXT, str(w_formula),
             str(w_proof), str(w_terminals_path)], timeout_per_step)
        if generated.returncode != 0:
            raise HierarchyViolation(
                f"W terminal generator failed: {generated.stdout} {generated.stderr}")
        w_generator = json.loads(generated.stdout.strip())
        w_terminals = _parse_w_terminals(w_terminals_path, source_index)
        expected_w = (ramsey_records[source_index]["w_completions_tested"]
                      if ramsey_mode else
                      signed_records[source_index]["w_completions_tested"])
        if len(w_terminals) != expected_w:
            raise HierarchyViolation("W terminal count mismatch")
        t_records = []
        signed_terminal_total = 0
        for w_terminal in w_terminals:
            w_index = w_terminal["terminal_index"]
            w_trits = w_terminal["w_half_trits"]
            _write_t_input(
                instance, source_index, w_trits, t_input, ramsey_mode)
            generator = T_RAMSEY if ramsey_mode else T_DFS
            command = [str(generator), str(t_input), ROWS_TEXT,
                       str(t_formula), str(t_proof)]
            if ramsey_mode:
                command.append(str(t_terminals_path))
            generated_t, t_generate_wall = _run(command, timeout_per_step)
            if generated_t.returncode != 0:
                raise HierarchyViolation(
                    f"T generator failed for {source_index}/{w_index}: "
                    f"{generated_t.stdout} {generated_t.stderr}")
            t_generator = json.loads(generated_t.stdout.strip())
            t_terminals = []
            if ramsey_mode:
                t_terminals = _parse_t_terminals(
                    t_terminals_path, instance, w_trits, source_index)
                if len(t_terminals) != t_generator[
                        "terminal_signed_completions"]:
                    raise HierarchyViolation("T terminal count mismatch")
                signed_terminal_total += len(t_terminals)
            checked_t = _checked_bundle(
                t_formula, t_proof, t_kernel, timeout_per_step, large=False)
            prefix = f"s{source_index:03d}-w{w_index:03d}-t"
            archived_t = _archive_bundle(
                prefix, proof_dir, t_proof, t_kernel, checked_t)
            t_records.append({
                "w_terminal_index": w_index,
                "w_formula_id": w_terminal["formula_id"],
                "w_half_trits_sha256":
                    hashlib.sha256(w_trits.encode()).hexdigest(),
                "generator": t_generator,
                "ramsey_terminals": t_terminals,
                "certificate": archived_t,
                "generator_seconds": round(t_generate_wall, 6),
            })
        if ramsey_mode and signed_terminal_total != ramsey_records[
                source_index]["signed_completions_tested"]:
            raise HierarchyViolation("signed terminal aggregate mismatch")
        checked_w = _checked_bundle(
            w_formula, w_proof, w_kernel, timeout_per_step, large=True)
        archived_w = _archive_bundle(
            f"s{source_index:03d}-w", proof_dir,
            w_proof, w_kernel, checked_w)
        terminal_copy = proof_dir / f"s{source_index:03d}.wterm"
        shutil.copyfile(w_terminals_path, terminal_copy)
        record = {
            "source_index": source_index,
            "target_ordinal": target_ordinal,
            "kind": "RAMSEY_NEGATIVE" if ramsey_mode else "SIGNED_SQUARE_NEGATIVE",
            "status": "VERIPB_CAKEPB_VERIFIED",
            "row_domains_sha256": dump_record["row_domains_sha256"],
            "w_generator": w_generator,
            "w_terminal_sidecar": {
                **_record(terminal_copy),
                "relative_path":
                    f"r55/data/involution_f5_hierarchical_proofs/{terminal_copy.name}",
            },
            "t_certificates": t_records,
            "w_certificate": archived_w,
            "w_generator_seconds": round(w_generate_wall, 6),
            "unchecked_deletion_present": False,
        }
        records.append(record)
        completed.add(source_index)
        _write_progress(output, records, len(targets), started)
        print(json.dumps({
            "certified": len(records),
            "target": len(targets),
            "source_index": source_index,
            "w_terminals": len(w_terminals),
            "signed_terminals": signed_terminal_total,
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
