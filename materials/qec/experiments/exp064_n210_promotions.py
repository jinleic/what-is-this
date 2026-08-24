"""EXP-064: replayed exact certificates for adaptive n=210 frontier points.

Targets are exact results already produced by EXP-060's parity/duality ratchet.
Each promotion rebuilds the constructor, replays the decisive one-side UNSAT
cap with raw kissat404, verifies the minimum-weight witness through NumPy and
bitset paths, and checks the exact BB d_X=d_Z isometry.

Run:
  python experiments/exp064_n210_promotions.py promote --target k24d4 --force
  python experiments/exp064_n210_promotions.py validate --target k24d4
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


E56 = _load("exp056_for_exp064", "exp056_odd_distance.py")
E58 = _load("exp058_for_exp064", "exp058_odd_frontier.py")

SCHEMA = "exp064-n210-frontier-distance-v1"
SOLVER = "kissat404"
RESULT_DIR = ROOT / "results" / "partial_runs" / "exp060_n210"
TARGETS: dict[str, dict[str, Any]] = {
    "k24d4": {
        "certificate": "results/certificates/exp064_210_24_4_distance.json",
        "result": "no_reference_000.json",
        "name": "exp060-15x7-frontier-210-24-4",
        "ell": 15, "m": 7, "expected_k": 24, "expected_d": 4,
        "lower_cap": 3,
        "A": [[0, 0], [0, 1], [0, 3]],
        "B": [[0, 0], [1, 0], [4, 0]],
    },
    "k14d12": {
        "certificate": "results/certificates/exp064_210_14_12_distance.json",
        "result": "undecided_000.json",
        "name": "exp060-15x7-frontier-210-14-12",
        "ell": 15, "m": 7, "expected_k": 14, "expected_d": 12,
        "lower_cap": 11,
        "A": [[0, 0], [1, 1], [4, 3]],
        "B": [[0, 0], [1, 2], [4, 6]],
    },
    "k10d16": {
        "certificate": "results/certificates/exp064_210_10_16_distance.json",
        "result": "undecided_001.json",
        "name": "exp060-15x7-frontier-210-10-16",
        "ell": 15, "m": 7, "expected_k": 10, "expected_d": 16,
        "lower_cap": 15,
        "A": [[0, 0], [1, 1], [2, 3]],
        "B": [[0, 0], [1, 6], [11, 2]],
    },
}


def certificate_path(key: str) -> Path:
    return ROOT / TARGETS[key]["certificate"]


def problem_of(key: str) -> dict[str, Any]:
    target = {field: value for field, value in TARGETS[key].items()
              if field not in {"certificate", "result"}}
    return E56.build_problem(target)


def source_result(key: str) -> dict[str, Any]:
    return json.loads((RESULT_DIR / TARGETS[key]["result"]).read_text(encoding="utf-8"))


def target_payload(key: str, witness_support: list[int]) -> dict[str, Any]:
    target = {field: value for field, value in TARGETS[key].items()
              if field not in {"certificate", "result"}}
    target["witness_support"] = witness_support
    return target


def exact_witness(result: dict[str, Any]) -> tuple[int, list[int]]:
    distance = int(result["exact_distance"])
    for decision in result["decisions"]:
        if decision.get("status") == "SAT" and int(decision.get("weight", -1)) == distance:
            return distance, [int(index) for index in decision["support"]]
    upper = result["upper"]
    if int(upper["weight"]) == distance:
        return distance, [int(index) for index in upper["support"]]
    raise RuntimeError("EXP-064 result has no exact witness")


def identity(key: str, problem: dict[str, Any], target: dict[str, Any]) -> dict[str, Any]:
    return {
        "target_key": key,
        "target": target,
        "n": problem["n"], "k": problem["k"],
        "HX_sha256": E56.matrix_sha256(problem["HX"]),
        "HZ_sha256": E56.matrix_sha256(problem["HZ"]),
    }


def protocol(key: str, problem: dict[str, Any]) -> dict[str, Any]:
    target = TARGETS[key]
    parity = E56.kernel_weight_parity_proof(problem, int(target["lower_cap"]))
    return {
        "encoding_version": E58.ENCODING_VERSION,
        "solver": SOLVER, "replay_solver": SOLVER,
        "decomposition": "monolithic", "side": "z",
        "requested_cap": int(target["lower_cap"]),
        "effective_even_cap": int(parity["effective_even_cap"]),
        "conflict_budget": 0, "replay_required": True,
    }


def validate_exact_certificate_payload(payload: dict[str, Any]) -> None:
    key = payload.get("identity", {}).get("target_key")
    if key not in TARGETS:
        raise RuntimeError("EXP-064 certificate target is unknown")
    result = source_result(key)
    distance, support = exact_witness(result)
    problem = problem_of(key)
    target = target_payload(key, support)
    expected_protocol = protocol(key, problem)
    cap = int(expected_protocol["effective_even_cap"])
    instance = E58.css_side_instance(
        problem["HX"], problem["HZ"], "z", block_length=problem["block"]
    )
    digest = E58.decision_cnf_digest(instance, cap)
    vector = np.zeros(problem["n"], dtype=np.uint8)
    vector[np.asarray(support, dtype=int)] = 1
    parity = E56.kernel_weight_parity_proof(problem, int(TARGETS[key]["lower_cap"]))
    duality = E56.bb_duality_proof(problem)
    upper = E56.verify_upper_witness(problem, vector)
    try:
        lower = payload["lower_bound"]
        valid = bool(
            payload.get("schema") == SCHEMA
            and payload.get("identity") == identity(key, problem, target)
            and payload.get("protocol") == expected_protocol
            and payload.get("protocol_sha256") == E56.canonical_json_sha256(expected_protocol)
            and payload.get("parity") == parity
            and payload.get("duality") == {k: v for k, v in duality.items() if k != "permutation"}
            and payload.get("upper") == upper
            and E58._decision_valid(lower["initial"], SOLVER, digest, cap)
            and E58._decision_valid(lower["replay"], SOLVER, digest, cap)
            and payload.get("verdict") == {
                "classification": "CERTIFIED_EXACT",
                "d": distance, "d_X": distance, "d_Z": distance, "exact": True,
            }
        )
    except (KeyError, TypeError, ValueError):
        valid = False
    if not valid:
        raise RuntimeError("EXP-064 exact certificate does not validate")


def promote(key: str, *, force: bool) -> int:
    path = certificate_path(key)
    if path.exists() and not force:
        payload = json.loads(path.read_text(encoding="utf-8"))
        validate_exact_certificate_payload(payload)
        print(json.dumps(payload["verdict"], indent=1))
        return 0
    result = source_result(key)
    if result.get("exact") is not True:
        raise RuntimeError("EXP-064 source ratchet result is not exact")
    distance, support = exact_witness(result)
    problem = problem_of(key)
    target = target_payload(key, support)
    E56.atomic_write_json(path, {
        "schema": SCHEMA,
        "identity": identity(key, problem, target),
        "verdict": {"classification": "RUNNING", "exact": False},
    })
    expected_protocol = protocol(key, problem)
    cap = int(expected_protocol["effective_even_cap"])
    initial = next(
        decision for decision in result["decisions"]
        if decision["status"] == "UNSAT" and int(decision["weight_cap"]) == cap
    )
    instance = E58.css_side_instance(
        problem["HX"], problem["HZ"], "z", block_length=problem["block"]
    )
    replay = E58.decide_raw(instance, cap, SOLVER)
    parity = E56.kernel_weight_parity_proof(problem, int(TARGETS[key]["lower_cap"]))
    duality = E56.bb_duality_proof(problem)
    vector = np.zeros(problem["n"], dtype=np.uint8)
    vector[np.asarray(support, dtype=int)] = 1
    upper = E56.verify_upper_witness(problem, vector)
    exact = bool(
        replay["status"] == "UNSAT"
        and initial["cnf_sha256"] == replay["cnf_sha256"]
        and upper["weight"] == distance
        and cap + 2 == distance
        and parity["valid"] and duality["valid"]
    )
    payload = {
        "schema": SCHEMA,
        "identity": identity(key, problem, target),
        "protocol": expected_protocol,
        "protocol_sha256": E56.canonical_json_sha256(expected_protocol),
        "parity": parity,
        "duality": {k: v for k, v in duality.items() if k != "permutation"},
        "upper": upper,
        "lower_bound": {"initial": initial, "replay": replay},
        "verdict": ({
            "classification": "CERTIFIED_EXACT",
            "d": distance, "d_X": distance, "d_Z": distance, "exact": True,
        } if exact else {"classification": "REPLAY_INCOMPLETE", "exact": False}),
    }
    E56.atomic_write_json(path, payload)
    if not exact:
        return 1
    validate_exact_certificate_payload(payload)
    print(json.dumps(payload["verdict"], indent=1))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("promote", "validate"))
    parser.add_argument("--target", choices=tuple(TARGETS), required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.command == "promote":
        return promote(args.target, force=args.force)
    payload = json.loads(certificate_path(args.target).read_text(encoding="utf-8"))
    validate_exact_certificate_payload(payload)
    print(json.dumps(payload["verdict"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
