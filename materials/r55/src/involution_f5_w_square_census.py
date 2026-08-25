#!/usr/bin/env python3
"""Run the exact W-row clique gate over all 705 Ramsey-support survivors."""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import struct
import subprocess
import sys
import time
from pathlib import Path

import involution_f5_signed_completion as completion

SCHEMA_VERSION = 1
CAMPAIGN_ID = "involution_f5_w_square_census"
ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parents[1]
DEFAULT_ENGINE_SOURCE = ROOT / "src" / "involution_f5_w_square.cpp"
DEFAULT_ENGINE = ROOT / "src" / "involution_f5_w_square_c"
DEFAULT_WORK = WORKSPACE / "scratch" / "r55-involution-f5" / "w-square-current.dom"
DEFAULT_OUTPUT = ROOT / "data" / "involution_f5_w_square_census.json"
CONTROL_SOURCE_INDEX = 393


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()




def write_domain_dump(instance: completion.CompletionInstance,
                      path: Path) -> dict:
    """Write one deterministic binary W-domain table for the C++ engine."""
    path.parent.mkdir(parents=True, exist_ok=True)
    row_counts = []
    domain_digest = hashlib.sha256()
    with path.open("wb") as stream:
        stream.write(b"WDOM1")
        stream.write(struct.pack(">I", instance.source_index))
        stream.write(bytes(instance.masks))
        stream.write(bytes(instance.internal))
        for row in range(20):
            domains = completion.enumerate_w_row_domains(instance, row)
            row_counts.append(len(domains))
            stream.write(struct.pack(">I", len(domains)))
            domain_digest.update(bytes((row,)))
            domain_digest.update(struct.pack(">I", len(domains)))
            for domain in domains:
                code = completion.w_domain_code(domain)
                encoded = struct.pack(">I", code)
                stream.write(encoded)
                domain_digest.update(encoded)
    return {
        "row_domain_counts": row_counts,
        "total_row_domains": sum(row_counts),
        "row_domains_sha256": domain_digest.hexdigest(),
        "dump_bytes": path.stat().st_size,
        "dump_sha256": _sha256(path),
    }


def run_engine(engine: Path, dump: Path, timeout: float | None = None) -> dict:
    started = time.monotonic()
    process = subprocess.run(
        [str(engine), str(dump)],
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    wall = time.monotonic() - started
    if process.returncode not in (0, 20):
        raise completion.CompletionViolation(
            f"W engine failed with {process.returncode}: {process.stderr.strip()}")
    lines = [line for line in process.stdout.splitlines() if line.strip()]
    if len(lines) != 1:
        raise completion.CompletionViolation("W engine did not emit one JSON record")
    try:
        record = json.loads(lines[0])
    except json.JSONDecodeError as error:
        raise completion.CompletionViolation("W engine emitted invalid JSON") from error
    expected_sat = process.returncode == 0
    if type(record) is not dict or record.get("w_square_satisfiable") is not expected_sat:
        raise completion.CompletionViolation("W engine status/exit-code mismatch")
    record["driver_wall_seconds"] = round(wall, 6)
    return record


def verify_w_witness(instance: completion.CompletionInstance,
                     trits: str) -> None:
    if type(trits) is not str or len(trits) != 190 or set(trits) - set("012"):
        raise completion.CompletionViolation("invalid W half-trit witness")
    rows = [[0] * 20 for _ in range(20)]
    position = 0
    for left in range(20):
        for right in range(left + 1, 20):
            value = int(trits[position]) - 1
            position += 1
            rows[left][right] = rows[right][left] = value
    for row in range(20):
        completion.verify_w_row_domain(instance, row, tuple(rows[row]))
    columns = instance.r_columns
    for left in range(20):
        for right in range(20):
            value = 0
            for index in range(20):
                left_value = (1 - 2 * instance.internal[left]
                              if index == left else 2 * rows[left][index])
                right_value = (1 - 2 * instance.internal[right]
                               if index == right else 2 * rows[right][index])
                value += left_value * right_value
            gram = sum(columns[left][vertex] * columns[right][vertex]
                       for vertex in range(5))
            if value != 45 * (left == right) - 2 - 2 * gram:
                raise completion.CompletionViolation("W witness fails the square equation")


def _deterministic_record(record: dict) -> dict:
    return {key: value for key, value in record.items()
            if key not in ("wall_seconds", "driver_wall_seconds")}


def build_analysis(engine: Path = DEFAULT_ENGINE, work: Path = DEFAULT_WORK,
                   *, timeout_per_candidate: float | None = None) -> dict:
    if not engine.is_file():
        raise completion.CompletionViolation(
            f"missing W engine {engine}; compile {DEFAULT_ENGINE_SOURCE.name}")
    records = []
    exact_digest = hashlib.sha256()
    total_row_domains = 0
    total_nodes = 0
    total_checks = 0
    total_wall = 0.0

    control = dataclasses.replace(
        completion.load_instance(CONTROL_SOURCE_INDEX), restrictions=())
    control_domains = write_domain_dump(control, work)
    control_result = run_engine(engine, work, timeout_per_candidate)
    if not control_result["w_square_satisfiable"]:
        raise completion.CompletionViolation("unfiltered published control lacks a W witness")
    verify_w_witness(control, control_result["w_half_trits"])
    control_record = {
        "source_index": CONTROL_SOURCE_INDEX,
        "ramsey_relation_restrictions_applied": False,
        **control_domains,
        **control_result,
    }

    for ordinal, instance in enumerate(completion.load_all_instances()):
        domain_record = write_domain_dump(instance, work)
        result = run_engine(engine, work, timeout_per_candidate)
        if result.get("source_index") != instance.source_index:
            raise completion.CompletionViolation("W engine source-index mismatch")
        if result["w_square_satisfiable"]:
            verify_w_witness(instance, result.get("w_half_trits"))
        record = {
            "ordinal": ordinal,
            "source_index": instance.source_index,
            "orbit_size": instance.orbit_size,
            **domain_record,
            **result,
        }
        deterministic = _deterministic_record(record)
        exact_digest.update(json.dumps(
            deterministic, sort_keys=True, separators=(",", ":")).encode())
        exact_digest.update(b"\n")
        total_row_domains += domain_record["total_row_domains"]
        total_nodes += result["nodes"]
        total_checks += result["compatibility_checks"]
        total_wall += result["driver_wall_seconds"]
        records.append(record)
        print(json.dumps({
            "ordinal": ordinal,
            "source_index": instance.source_index,
            "w_square_satisfiable": result["w_square_satisfiable"],
            "nodes": result["nodes"],
            "wall_seconds": result["driver_wall_seconds"],
        }, sort_keys=True), file=sys.stderr, flush=True)

    satisfiable = [record for record in records
                   if record["w_square_satisfiable"]]
    if satisfiable:
        disposition = (
            "F5_W_SQUARE_CENSUS_HAS_"
            f"{len(satisfiable)}_SURVIVING_SUPPORT_ORBITS"
        )
    else:
        disposition = "F5_RAMSEY_FIXED5_W_SQUARE_EMPTY_OVER_705_D5_ORBITS"
    return {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "disposition": disposition,
        "input_dependencies": {
            "ramsey_frontier": {
                "relative_path": "r55/data/involution_f5_ramsey_filter.json",
                "bytes": completion.FRONTIER_BYTES,
                "sha256": completion.FRONTIER_SHA256,
            },
            "support_census": {
                "relative_path": "r55/data/involution_f5_support_census.json",
                "bytes": completion.SUPPORT_BYTES,
                "sha256": completion.SUPPORT_SHA256,
            },
        },
        "engine": {
            "relative_source_path": "r55/src/involution_f5_w_square.cpp",
            "source_sha256": _sha256(DEFAULT_ENGINE_SOURCE),
            "binary_sha256": _sha256(engine),
            "algorithm": (
                "exact 20-partite clique over complete linear W row domains; "
                "pair compatibility is symmetry plus the fixed W-row Gram entry"
            ),
        },
        "unfiltered_positive_control": control_record,
        "census": {
            "input_support_orbits": len(records),
            "input_labelled_supports": sum(record["orbit_size"]
                                            for record in records),
            "w_square_satisfiable_support_orbits": len(satisfiable),
            "w_square_unsatisfiable_support_orbits": len(records) - len(satisfiable),
            "total_row_domains": total_row_domains,
            "total_search_nodes": total_nodes,
            "total_compatibility_checks": total_checks,
            "measured_driver_wall_seconds": round(total_wall, 6),
            "deterministic_records_sha256": exact_digest.hexdigest(),
            "candidates": records,
        },
        "downstream": {
            "t_square_candidate_support_orbits": len(satisfiable),
            "graph_reconstruction_candidate_support_orbits": len(satisfiable),
            "negative_certificate_status": "PENDING_VERIPB_CAKEPB",
            "general_ramsey_bound_claimed": False,
        },
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", type=Path, default=DEFAULT_ENGINE)
    parser.add_argument("--work", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timeout-per-candidate", type=float)
    args = parser.parse_args(argv)
    document = build_analysis(
        args.engine, args.work,
        timeout_per_candidate=args.timeout_per_candidate)
    args.output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "disposition": document["disposition"],
        "w_square_satisfiable_support_orbits":
            document["census"]["w_square_satisfiable_support_orbits"],
        "w_square_unsatisfiable_support_orbits":
            document["census"]["w_square_unsatisfiable_support_orbits"],
        "output": str(args.output),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
