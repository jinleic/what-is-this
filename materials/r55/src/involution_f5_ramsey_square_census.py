#!/usr/bin/env python3
"""Exhaust signed square completions for K5/I5 over all 108 survivors."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import involution_f5_signed_completion as completion
import involution_f5_t_square_census as t_census
import involution_f5_w_square_census as w_census

SCHEMA_VERSION = 1
CAMPAIGN_ID = "involution_f5_ramsey_square_census"
ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parents[1]
DEFAULT_SIGNED_CENSUS = ROOT / "data" / "involution_f5_signed_square_census.json"
DEFAULT_ENGINE_SOURCE = ROOT / "src" / "involution_f5_ramsey_square.cpp"
DEFAULT_ENGINE = ROOT / "src" / "involution_f5_ramsey_square_c"
DEFAULT_WORK = WORKSPACE / "scratch/r55-involution-f5/ramsey-square-current.dom"
DEFAULT_OUTPUT = ROOT / "data" / "involution_f5_ramsey_square_census.json"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def run_engine(engine: Path, dump: Path,
               timeout: float | None = None) -> dict:
    started = time.monotonic()
    process = subprocess.run(
        [str(engine), str(dump)], check=False,
        capture_output=True, text=True, timeout=timeout)
    wall = time.monotonic() - started
    if process.returncode not in (0, 20):
        raise completion.CompletionViolation(
            f"Ramsey engine failed with {process.returncode}: {process.stderr}")
    lines = [line for line in process.stdout.splitlines() if line.strip()]
    if len(lines) != 1:
        raise completion.CompletionViolation(
            "Ramsey engine did not emit one JSON record")
    result = json.loads(lines[0])
    if result.get("ramsey_square_satisfiable") is not (process.returncode == 0):
        raise completion.CompletionViolation(
            "Ramsey engine status/exit-code mismatch")
    result["driver_wall_seconds"] = round(wall, 6)
    return result


def build_analysis(signed_path: Path = DEFAULT_SIGNED_CENSUS,
                   engine: Path = DEFAULT_ENGINE,
                   work: Path = DEFAULT_WORK,
                   *, timeout_per_candidate: float | None = None) -> dict:
    raw = signed_path.read_bytes()
    signed_document = json.loads(raw)
    if signed_document.get("census", {}).get(
            "signed_square_satisfiable_support_orbits") != 108:
        raise completion.CompletionViolation("unexpected signed survivor count")
    if not engine.is_file():
        raise completion.CompletionViolation(f"missing Ramsey engine {engine}")
    records = []
    total_wall = 0.0
    for candidate in signed_document["census"]["candidates"]:
        if not candidate["signed_square_satisfiable"]:
            continue
        source_index = candidate["source_index"]
        instance = completion.load_instance(source_index)
        domain_record = w_census.write_domain_dump(instance, work)
        result = run_engine(engine, work, timeout_per_candidate)
        if result.get("source_index") != source_index:
            raise completion.CompletionViolation(
                "Ramsey engine source-index mismatch")
        record = {
            "source_index": source_index,
            "orbit_size": instance.orbit_size,
            "row_domains_sha256": domain_record["row_domains_sha256"],
            **result,
        }
        if result["ramsey_square_satisfiable"]:
            w_census.verify_w_witness(instance, result["w_half_trits"])
            states = t_census.relation_states(
                result["w_half_trits"], result["t_half_trits"])
            graph_check = completion.verify_state_assignment(instance, states)
            if (not graph_check["srg_45_22_10_11"]
                    or not graph_check["ramsey_good"]):
                raise completion.CompletionViolation(
                    "Ramsey witness failed independent graph verification")
            record["graph_check"] = graph_check
        if (result["signed_completions_with_k5"]
                > result["signed_completions_tested"]
                or result["signed_completions_with_i5"]
                > result["signed_completions_tested"]):
            raise completion.CompletionViolation(
                "invalid obstruction counters from Ramsey engine")
        total_wall += result["driver_wall_seconds"]
        records.append(record)
        print(json.dumps({
            "source_index": source_index,
            "ramsey_square_satisfiable": result["ramsey_square_satisfiable"],
            "signed_completions_tested": result["signed_completions_tested"],
            "wall_seconds": result["driver_wall_seconds"],
        }, sort_keys=True), file=sys.stderr, flush=True)
    ramsey = [record for record in records
              if record["ramsey_square_satisfiable"]]
    disposition = (
        "F5_RAMSEY_SIGNED_SQUARE_EMPTY_OVER_108_SUPPORT_ORBITS"
        if not ramsey else
        f"F5_RAMSEY_SIGNED_SQUARE_HAS_{len(ramsey)}_SUPPORT_ORBITS"
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "disposition": disposition,
        "input_dependency": {
            "relative_path": "r55/data/involution_f5_signed_square_census.json",
            "bytes": len(raw),
            "sha256": hashlib.sha256(raw).hexdigest(),
        },
        "engine": {
            "relative_source_path": "r55/src/involution_f5_ramsey_square.cpp",
            "source_sha256": _sha256(DEFAULT_ENGINE_SOURCE),
            "binary_sha256": _sha256(engine),
            "algorithm": (
                "exhaustive W-row clique search; exhaustive T-row clique "
                "search modulo a spanning-tree switching gauge; exact K5/I5 "
                "check on every signed completion"
            ),
            "switching_gauge_preserves_ramsey_property": True,
        },
        "census": {
            "signed_square_survivor_support_orbits": len(records),
            "ramsey_good_support_orbits": len(ramsey),
            "ramsey_bad_support_orbits": len(records) - len(ramsey),
            "signed_completions_tested": sum(
                record["signed_completions_tested"] for record in records),
            "signed_completions_with_k5": sum(
                record["signed_completions_with_k5"] for record in records),
            "signed_completions_with_i5": sum(
                record["signed_completions_with_i5"] for record in records),
            "measured_driver_wall_seconds": round(total_wall, 6),
            "candidates": records,
        },
        "claim": {
            "complete_over_all_w_solutions": True,
            "complete_over_t_switching_classes": True,
            "all_d5_images_covered_by_explicit_prior_expansion": True,
            "negative_certificate_status": "PENDING_VERIPB_CAKEPB",
            "strongly_regular_fixed_five_branch_only": True,
            "general_ramsey_bound_claimed": False,
        },
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--signed-census", type=Path, default=DEFAULT_SIGNED_CENSUS)
    parser.add_argument("--engine", type=Path, default=DEFAULT_ENGINE)
    parser.add_argument("--work", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timeout-per-candidate", type=float)
    args = parser.parse_args(argv)
    document = build_analysis(
        args.signed_census, args.engine, args.work,
        timeout_per_candidate=args.timeout_per_candidate)
    args.output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "disposition": document["disposition"],
        "ramsey_good_support_orbits":
            document["census"]["ramsey_good_support_orbits"],
        "signed_completions_tested":
            document["census"]["signed_completions_tested"],
        "output": str(args.output),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
