#!/usr/bin/env python3
"""Test exact T-square completion of the first W witness for each W survivor."""

from __future__ import annotations

import argparse
import hashlib
import json
import struct
import subprocess
import sys
import time
from pathlib import Path

import involution_f5_signed_completion as completion
import involution_f5_w_square_census as w_census

SCHEMA_VERSION = 1
CAMPAIGN_ID = "involution_f5_t_square_census"
ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parents[1]
DEFAULT_W_CENSUS = ROOT / "data" / "involution_f5_w_square_census.json"
DEFAULT_ENGINE_SOURCE = ROOT / "src" / "involution_f5_t_square.cpp"
DEFAULT_ENGINE = ROOT / "src" / "involution_f5_t_square_c"
DEFAULT_WORK = WORKSPACE / "scratch/r55-involution-f5/t-square-current.dom"
DEFAULT_OUTPUT = ROOT / "data" / "involution_f5_t_square_census.json"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def write_t_input(instance: completion.CompletionInstance, w_trits: str,
                  path: Path) -> dict:
    w_census.verify_w_witness(instance, w_trits)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as stream:
        stream.write(b"TDOM1")
        stream.write(struct.pack(">I", instance.source_index))
        stream.write(bytes(instance.internal))
        stream.write(w_trits.encode("ascii"))
    return {"bytes": path.stat().st_size, "sha256": _sha256(path)}


def run_engine(engine: Path, path: Path,
               timeout: float | None = None) -> dict:
    started = time.monotonic()
    process = subprocess.run(
        [str(engine), str(path)],
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    wall = time.monotonic() - started
    if process.returncode not in (0, 20):
        raise completion.CompletionViolation(
            f"T engine failed with {process.returncode}: {process.stderr}")
    lines = [line for line in process.stdout.splitlines() if line.strip()]
    if len(lines) != 1:
        raise completion.CompletionViolation("T engine did not emit one JSON row")
    result = json.loads(lines[0])
    if result.get("t_square_satisfiable") is not (process.returncode == 0):
        raise completion.CompletionViolation("T engine status/exit-code mismatch")
    result["driver_wall_seconds"] = round(wall, 6)
    return result


def relation_states(w_trits: str, t_trits: str) -> dict[tuple[int, int], str]:
    if (type(w_trits) is not str or type(t_trits) is not str
            or len(w_trits) != 190 or len(t_trits) != 190):
        raise completion.CompletionViolation("signed witnesses must have 190 trits")
    states = {}
    position = 0
    for left in range(20):
        for right in range(left + 1, 20):
            w_value = int(w_trits[position]) - 1
            t_value = int(t_trits[position]) - 1
            position += 1
            if (w_value, t_value) == (1, 0):
                state = "00"
            elif (w_value, t_value) == (-1, 0):
                state = "11"
            elif (w_value, t_value) == (0, 1):
                state = "01"
            elif (w_value, t_value) == (0, -1):
                state = "10"
            else:
                raise completion.CompletionViolation(
                    "W/T witnesses do not have complementary support")
            states[left, right] = state
    return states


def build_analysis(w_path: Path = DEFAULT_W_CENSUS,
                   engine: Path = DEFAULT_ENGINE,
                   work: Path = DEFAULT_WORK,
                   *, timeout_per_candidate: float | None = None) -> dict:
    raw = w_path.read_bytes()
    w_document = json.loads(raw)
    if w_document.get("census", {}).get(
            "w_square_satisfiable_support_orbits") != 155:
        raise completion.CompletionViolation("unexpected W-square survivor count")
    if not engine.is_file():
        raise completion.CompletionViolation(f"missing T engine {engine}")
    records = []
    for candidate in w_document["census"]["candidates"]:
        if not candidate["w_square_satisfiable"]:
            continue
        source_index = candidate["source_index"]
        instance = completion.load_instance(source_index)
        w_trits = candidate["w_half_trits"]
        input_record = write_t_input(instance, w_trits, work)
        result = run_engine(engine, work, timeout_per_candidate)
        if result.get("source_index") != source_index:
            raise completion.CompletionViolation("T engine source-index mismatch")
        record = {
            "source_index": source_index,
            "orbit_size": instance.orbit_size,
            "w_half_trits": w_trits,
            "t_input": input_record,
            **result,
        }
        if result["t_square_satisfiable"]:
            states = relation_states(w_trits, result["t_half_trits"])
            graph_check = completion.verify_state_assignment(instance, states)
            if not graph_check["srg_45_22_10_11"]:
                raise completion.CompletionViolation(
                    "signed W/T witness does not reconstruct the required SRG")
            record["graph_check"] = graph_check
        records.append(record)
        print(json.dumps({
            "source_index": source_index,
            "t_square_satisfiable": result["t_square_satisfiable"],
            "wall_seconds": result["driver_wall_seconds"],
        }, sort_keys=True), file=sys.stderr, flush=True)
    t_sat = [record for record in records if record["t_square_satisfiable"]]
    ramsey_good = [record for record in t_sat
                   if record["graph_check"]["ramsey_good"]]
    return {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "disposition": (
            "F5_FIRST_W_T_SQUARE_DIAGNOSTIC_"
            f"{len(t_sat)}_OF_155_T_COMPLETABLE"
        ),
        "input_dependency": {
            "relative_path": "r55/data/involution_f5_w_square_census.json",
            "bytes": len(raw),
            "sha256": hashlib.sha256(raw).hexdigest(),
        },
        "engine": {
            "relative_source_path": "r55/src/involution_f5_t_square.cpp",
            "source_sha256": _sha256(DEFAULT_ENGINE_SOURCE),
            "binary_sha256": _sha256(engine),
            "algorithm": (
                "exact 20-partite clique over all 2^11 signed T row domains; "
                "pair compatibility is symmetry plus row orthogonality"
            ),
        },
        "census": {
            "w_square_survivor_support_orbits": len(records),
            "first_w_with_t_square_completion": len(t_sat),
            "first_w_without_t_square_completion": len(records) - len(t_sat),
            "first_w_t_completions_ramsey_good": len(ramsey_good),
            "candidates": records,
        },
        "claim": {
            "complete_over_all_w_solutions": False,
            "role": "diagnostic first-W pass before complete signed PB enumeration",
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
        "first_w_with_t_square_completion":
            document["census"]["first_w_with_t_square_completion"],
        "first_w_t_completions_ramsey_good":
            document["census"]["first_w_t_completions_ramsey_good"],
        "output": str(args.output),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
