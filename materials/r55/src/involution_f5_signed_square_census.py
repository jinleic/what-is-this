#!/usr/bin/env python3
"""Complete W/T square census over the 155 W-square support survivors."""

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
CAMPAIGN_ID = "involution_f5_signed_square_census"
ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parents[1]
DEFAULT_W_CENSUS = ROOT / "data" / "involution_f5_w_square_census.json"
DEFAULT_ENGINE_SOURCE = ROOT / "src" / "involution_f5_signed_square.cpp"
DEFAULT_ENGINE = ROOT / "src" / "involution_f5_signed_square_c"
DEFAULT_WORK = WORKSPACE / "scratch/r55-involution-f5/signed-square-current.dom"
DEFAULT_OUTPUT = ROOT / "data" / "involution_f5_signed_square_census.json"


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
            f"signed engine failed with {process.returncode}: {process.stderr}")
    lines = [line for line in process.stdout.splitlines() if line.strip()]
    if len(lines) != 1:
        raise completion.CompletionViolation(
            "signed engine did not emit one JSON record")
    result = json.loads(lines[0])
    if result.get("signed_square_satisfiable") is not (process.returncode == 0):
        raise completion.CompletionViolation(
            "signed engine status/exit-code mismatch")
    result["driver_wall_seconds"] = round(wall, 6)
    return result


def build_analysis(w_path: Path = DEFAULT_W_CENSUS,
                   engine: Path = DEFAULT_ENGINE,
                   work: Path = DEFAULT_WORK,
                   *, timeout_per_candidate: float | None = None) -> dict:
    raw = w_path.read_bytes()
    w_document = json.loads(raw)
    if w_document.get("census", {}).get(
            "w_square_satisfiable_support_orbits") != 155:
        raise completion.CompletionViolation("unexpected W survivor count")
    if not engine.is_file():
        raise completion.CompletionViolation(f"missing signed engine {engine}")
    records = []
    total_wall = 0.0
    for candidate in w_document["census"]["candidates"]:
        if not candidate["w_square_satisfiable"]:
            continue
        source_index = candidate["source_index"]
        instance = completion.load_instance(source_index)
        domain_record = w_census.write_domain_dump(instance, work)
        result = run_engine(engine, work, timeout_per_candidate)
        if result.get("source_index") != source_index:
            raise completion.CompletionViolation(
                "signed engine source-index mismatch")
        record = {
            "source_index": source_index,
            "orbit_size": instance.orbit_size,
            "row_domains_sha256": domain_record["row_domains_sha256"],
            **result,
        }
        if result["signed_square_satisfiable"]:
            w_census.verify_w_witness(instance, result["w_half_trits"])
            states = t_census.relation_states(
                result["w_half_trits"], result["t_half_trits"])
            graph_check = completion.verify_state_assignment(instance, states)
            if not graph_check["srg_45_22_10_11"]:
                raise completion.CompletionViolation(
                    "signed witness does not reconstruct the required SRG")
            record["graph_check"] = graph_check
        total_wall += result["driver_wall_seconds"]
        records.append(record)
        print(json.dumps({
            "source_index": source_index,
            "signed_square_satisfiable": result["signed_square_satisfiable"],
            "w_completions_tested": result["w_completions_tested"],
            "wall_seconds": result["driver_wall_seconds"],
        }, sort_keys=True), file=sys.stderr, flush=True)
    signed = [record for record in records
              if record["signed_square_satisfiable"]]
    ramsey_good_first = [record for record in signed
                         if record["graph_check"]["ramsey_good"]]
    return {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "disposition": (
            "F5_SIGNED_SQUARE_CENSUS_HAS_"
            f"{len(signed)}_SURVIVING_SUPPORT_ORBITS"
        ),
        "input_dependency": {
            "relative_path": "r55/data/involution_f5_w_square_census.json",
            "bytes": len(raw),
            "sha256": hashlib.sha256(raw).hexdigest(),
        },
        "engine": {
            "relative_source_path": "r55/src/involution_f5_signed_square.cpp",
            "source_sha256": _sha256(DEFAULT_ENGINE_SOURCE),
            "binary_sha256": _sha256(engine),
            "algorithm": (
                "exhaustive W-row clique search; each complete W invokes an "
                "exact T-row clique search modulo a positive spanning-tree "
                "switching gauge"
            ),
            "switching_gauge_theorem": (
                "conjugating T by a diagonal sign matrix preserves T^2 and "
                "can make every edge of a support spanning tree positive"
            ),
        },
        "census": {
            "w_square_survivor_support_orbits": len(records),
            "signed_square_satisfiable_support_orbits": len(signed),
            "signed_square_unsatisfiable_support_orbits":
                len(records) - len(signed),
            "first_signed_witness_ramsey_good": len(ramsey_good_first),
            "measured_driver_wall_seconds": round(total_wall, 6),
            "candidates": records,
        },
        "claim": {
            "complete_over_all_w_solutions": True,
            "complete_t_existence_for_each_w_support": True,
            "complete_ramsey_filter_over_all_signed_solutions": False,
            "general_ramsey_bound_claimed": False,
        },
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--w-census", type=Path, default=DEFAULT_W_CENSUS)
    parser.add_argument("--engine", type=Path, default=DEFAULT_ENGINE)
    parser.add_argument("--work", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timeout-per-candidate", type=float)
    args = parser.parse_args(argv)
    document = build_analysis(
        args.w_census, args.engine, args.work,
        timeout_per_candidate=args.timeout_per_candidate)
    args.output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "disposition": document["disposition"],
        "signed_square_satisfiable_support_orbits":
            document["census"]["signed_square_satisfiable_support_orbits"],
        "first_signed_witness_ramsey_good":
            document["census"]["first_signed_witness_ramsey_good"],
        "output": str(args.output),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
