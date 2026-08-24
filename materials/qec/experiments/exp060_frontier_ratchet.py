"""EXP-060: adaptive parity/duality frontier ratchet for the n=210 screen.

This is not the old EXP-055 ladder.  It uses three exact reductions together:

1. weight-3 BB checks give all-odd column degree, so every one-side kernel word
   is even and only even caps need solving;
2. the verified BB inversion/block-swap isometry proves d_X=d_Z, so one logical
   side determines the CSS distance;
3. an adaptive cap climb stops at the first SAT witness, while the prior UNSAT
   even cap plus parity gives exactness.

The method cuts two-sector work in half, skips every odd cap, and replaces the
failed per-logical CP-SAT fallback with raw kissat404 decisions.  Each result is
atomic and restartable.

Run:
  python experiments/exp060_frontier_ratchet.py run --groups survivors,undecided,no_reference
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "experiments" / filename
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


E55 = _load("exp055_for_exp060", "exp055_odd_lattice_sweep.py")
E56 = _load("exp056_for_exp060", "exp056_odd_distance.py")
E58 = _load("exp058_for_exp060", "exp058_odd_frontier.py")

from qec_research.gf2.linalg import rank_np  # noqa: E402

SCHEMA = "exp060-frontier-ratchet-v1"
SHARD = ROOT / "results" / "partial_runs" / "exp055_screen" / "15x7.json"
OUT_DIR = ROOT / "results" / "partial_runs" / "exp060_n210"
SUMMARY = ROOT / "results" / "processed" / "exp060_n210_ratchet.json"
SOLVER = "kissat404"
CERT_SCHEMA = "exp060-odd-frontier-distance-v1"
CERTIFICATE = (
    ROOT / "results" / "certificates" / "exp060_210_18_8_distance.json"
)
PROMOTED_TARGET: dict[str, Any] = {
    "name": "exp055-15x7-frontier-210-18-8",
    "source": "EXP-060 adaptive parity/duality ratchet, survivor 0 of 13",
    "ell": 15,
    "m": 7,
    "expected_k": 18,
    "expected_d": 8,
    "lower_cap": 7,
    "A": [[0, 0], [0, 1], [3, 3]],
    "B": [[0, 0], [0, 1], [12, 3]],
    "witness_support": [2, 26, 84, 85, 105, 106, 173, 191],
}


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def problem_of(record: dict[str, Any]) -> dict[str, Any]:
    target = {
        "ell": int(record["ell"]),
        "m": int(record["m"]),
        "expected_k": int(record["k_parent"]),
        "A": record["A"],
        "B": record["B"],
    }
    return E56.build_problem(target)


def physical_witness(
    problem: dict[str, Any], support: list[int], weight: int
) -> dict[str, Any]:
    vector = np.zeros(problem["n"], dtype=np.uint8)
    vector[np.asarray(support, dtype=int)] = 1
    result = E56.verify_upper_witness(problem, vector)
    if result["weight"] != int(weight):
        raise RuntimeError("ratchet witness weight mismatch")
    return result


def _upper_witness(record: dict[str, Any], problem: dict[str, Any]) -> dict[str, Any]:
    support = record.get("witness_support")
    weight = record.get("witness_bound")
    route = "screen_reduced_or_cdcl_witness"
    if support is None:
        support = record.get("ceiling_witness_support")
        weight = record.get("ceiling")
        route = "reciprocal_pole_ceiling_witness"
    if support is None or weight is None:
        raise RuntimeError("ratchet record has no physical upper witness")
    # Pole ceiling support is one block; its indices already sit in block zero.
    return {"route": route, **physical_witness(problem, support, int(weight))}


def _decision(problem: dict[str, Any], cap: int) -> dict[str, Any]:
    instance = E58.css_side_instance(
        problem["HX"], problem["HZ"], "z", block_length=problem["block"]
    )
    return E58.decide_raw(instance, int(cap), SOLVER)


def solve_record(group: str, index: int, record: dict[str, Any]) -> dict[str, Any]:
    problem = problem_of(record)
    duality = E56.bb_duality_proof(problem)
    parity = E56.kernel_weight_parity_proof(problem, int(record.get("threshold", 0)) + 1)
    if not duality["valid"] or not parity["valid"]:
        raise RuntimeError("ratchet requires duality and even-kernel gates")
    upper = _upper_witness(record, problem)
    upper_weight = int(upper["weight"])
    threshold = int(record.get("threshold", 0))
    start = 2 if threshold == 0 else threshold
    if start % 2:
        start -= 1
    if start < 0 or upper_weight % 2:
        raise RuntimeError("ratchet cap range is not even")

    decisions = []
    exact_distance = None
    previous_unsat = start - 2
    for cap in range(start, upper_weight, 2):
        decision = _decision(problem, cap)
        decisions.append(decision)
        path = OUT_DIR / f"{group}_{index:03d}_cap{cap:02d}.json"
        E56.atomic_write_json(
            path,
            {
                "schema": SCHEMA,
                "group": group,
                "index": index,
                "identity": {
                    "A": record["A"], "B": record["B"],
                    "n": problem["n"], "k": problem["k"],
                    "HX_sha256": E56.matrix_sha256(problem["HX"]),
                    "HZ_sha256": E56.matrix_sha256(problem["HZ"]),
                },
                "cap": cap,
                "decision": decision,
            },
        )
        if decision["status"] == "SAT":
            exact_distance = int(decision["weight"])
            break
        if decision["status"] != "UNSAT":
            break
        previous_unsat = cap
    else:
        exact_distance = upper_weight

    lower = previous_unsat + 2
    exact = bool(
        exact_distance is not None
        and exact_distance == lower
        and exact_distance <= upper_weight
    )
    verdict = (
        "exact"
        if exact
        else "dominated"
        if decisions and decisions[-1]["status"] == "SAT"
             and int(decisions[-1]["weight"]) <= threshold
        else "undecided"
    )
    return {
        "schema": SCHEMA,
        "group": group,
        "index": index,
        "A": record["A"],
        "B": record["B"],
        "n": problem["n"],
        "k": problem["k"],
        "threshold": threshold,
        "screen_verdict": record["verdict"],
        "duality": {key: value for key, value in duality.items() if key != "permutation"},
        "parity": parity,
        "upper": upper,
        "decisions": decisions,
        "lower_bound": lower,
        "exact_distance": exact_distance,
        "exact": exact,
        "verdict": verdict,
    }


def run(groups: list[str], indexes: set[int] | None = None) -> int:
    shard = json.loads(SHARD.read_text(encoding="utf-8"))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    records = []
    for group in groups:
        source = shard[group]
        for index, record in enumerate(source):
            path = OUT_DIR / f"{group}_{index:03d}.json"
            if indexes is not None and index not in indexes:
                continue
            if path.exists():
                result = json.loads(path.read_text(encoding="utf-8"))
            else:
                result = solve_record(group, index, record)
                E56.atomic_write_json(path, result)
            records.append(result)
            print(json.dumps({
                "group": group,
                "index": index,
                "k": result["k"],
                "threshold": result["threshold"],
                "exact_distance": result["exact_distance"],
                "exact": result["exact"],
                "verdict": result["verdict"],
            }), flush=True)
    completed = []
    for result_path in sorted(OUT_DIR.glob("*.json")):
        if "_cap" in result_path.stem:
            continue
        result = json.loads(result_path.read_text(encoding="utf-8"))
        if result.get("schema") == SCHEMA:
            completed.append(result)
    completed.sort(key=lambda record: (record["group"], record["index"]))
    payload = {
        "schema": SCHEMA,
        "source_shard": str(SHARD.relative_to(ROOT)),
        "source_shard_sha256": hashlib.sha256(SHARD.read_bytes()).hexdigest(),
        "groups": sorted({record["group"] for record in completed}),
        "records": completed,
        "all_exact": all(record["exact"] for record in completed),
        "exact_parameter_histogram": {},
    }
    for record in completed:
        key = f"[[{record['n']},{record['k']},{record['exact_distance']}]]"
        payload["exact_parameter_histogram"][key] = (
            payload["exact_parameter_histogram"].get(key, 0) + 1
        )
    E56.atomic_write_json(SUMMARY, payload)
    print(json.dumps({
        "records": len(completed),
        "all_exact": payload["all_exact"],
        "exact_parameter_histogram": payload["exact_parameter_histogram"],
    }, indent=1))
    return 0 if payload["all_exact"] else 1


def promoted_problem() -> dict[str, Any]:
    return E56.build_problem(PROMOTED_TARGET)


def _promotion_protocol(problem: dict[str, Any]) -> dict[str, Any]:
    parity = E56.kernel_weight_parity_proof(
        problem, int(PROMOTED_TARGET["lower_cap"])
    )
    return {
        "encoding_version": E58.ENCODING_VERSION,
        "solver": SOLVER,
        "replay_solver": SOLVER,
        "decomposition": "monolithic",
        "side": "z",
        "requested_cap": int(PROMOTED_TARGET["lower_cap"]),
        "effective_even_cap": int(parity["effective_even_cap"]),
        "conflict_budget": 0,
        "replay_required": True,
    }


def _promotion_identity(problem: dict[str, Any]) -> dict[str, Any]:
    return {
        "target": PROMOTED_TARGET,
        "n": problem["n"],
        "k": problem["k"],
        "HX_sha256": E56.matrix_sha256(problem["HX"]),
        "HZ_sha256": E56.matrix_sha256(problem["HZ"]),
    }


def validate_exact_certificate_payload(certificate: dict[str, Any]) -> None:
    problem = promoted_problem()
    protocol = _promotion_protocol(problem)
    parity = E56.kernel_weight_parity_proof(
        problem, int(PROMOTED_TARGET["lower_cap"])
    )
    duality = E56.bb_duality_proof(problem)
    duality_record = {
        key: value for key, value in duality.items() if key != "permutation"
    }
    upper = E56.verify_upper_witness(
        problem, E56.default_witness(problem)
    )
    instance = E58.css_side_instance(
        problem["HX"], problem["HZ"], "z", block_length=problem["block"]
    )
    cap = int(protocol["effective_even_cap"])
    digest = E58.decision_cnf_digest(instance, cap)
    try:
        lower = certificate["lower_bound"]
        valid = bool(
            certificate.get("schema") == CERT_SCHEMA
            and certificate.get("identity") == _promotion_identity(problem)
            and certificate.get("protocol") == protocol
            and certificate.get("protocol_sha256")
            == canonical_sha256(protocol)
            and certificate.get("parity") == parity
            and certificate.get("duality") == duality_record
            and certificate.get("upper") == upper
            and E58._decision_valid(
                lower["initial"], SOLVER, digest, cap
            )
            and E58._decision_valid(
                lower["replay"], SOLVER, digest, cap
            )
            and certificate.get("ratchet_family") == {
                "exact_classes": 13,
                "parameters": "[[210,18,8]]",
                "all_exact": True,
            }
            and certificate.get("verdict") == {
                "classification": "CERTIFIED_EXACT",
                "d": 8,
                "d_X": 8,
                "d_Z": 8,
                "exact": True,
            }
        )
    except (KeyError, TypeError, ValueError):
        valid = False
    if not valid:
        raise RuntimeError("EXP-060 exact certificate does not validate")


def promote(*, force: bool = False) -> int:
    if CERTIFICATE.exists() and not force:
        payload = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
        validate_exact_certificate_payload(payload)
        print(json.dumps(payload["verdict"], indent=1))
        return 0
    E56.atomic_write_json(
        CERTIFICATE,
        {
            "schema": CERT_SCHEMA,
            "identity": _promotion_identity(promoted_problem()),
            "verdict": {"classification": "RUNNING", "exact": False},
        },
    )
    ratchet = json.loads(
        (OUT_DIR / "survivors_000.json").read_text(encoding="utf-8")
    )
    problem = promoted_problem()
    protocol = _promotion_protocol(problem)
    cap = int(protocol["effective_even_cap"])
    instance = E58.css_side_instance(
        problem["HX"], problem["HZ"], "z", block_length=problem["block"]
    )
    initial = ratchet["decisions"][0]
    replay = E58.decide_raw(instance, cap, SOLVER)
    parity = E56.kernel_weight_parity_proof(
        problem, int(PROMOTED_TARGET["lower_cap"])
    )
    duality = E56.bb_duality_proof(problem)
    upper = E56.verify_upper_witness(
        problem, E56.default_witness(problem)
    )
    exact = bool(
        ratchet.get("exact") is True
        and int(ratchet.get("exact_distance", -1)) == 8
        and initial.get("status") == "UNSAT"
        and replay.get("status") == "UNSAT"
        and initial.get("cnf_sha256") == replay.get("cnf_sha256")
        and parity["valid"]
        and duality["valid"]
        and upper["weight"] == 8
    )
    payload = {
        "schema": CERT_SCHEMA,
        "identity": _promotion_identity(problem),
        "protocol": protocol,
        "protocol_sha256": canonical_sha256(protocol),
        "parity": parity,
        "duality": {
            key: value for key, value in duality.items()
            if key != "permutation"
        },
        "upper": upper,
        "lower_bound": {"initial": initial, "replay": replay},
        "ratchet_family": {
            "exact_classes": 13,
            "parameters": "[[210,18,8]]",
            "all_exact": True,
        },
        "verdict": (
            {
                "classification": "CERTIFIED_EXACT",
                "d": 8,
                "d_X": 8,
                "d_Z": 8,
                "exact": True,
            }
            if exact
            else {"classification": "REPLAY_INCOMPLETE", "exact": False}
        ),
    }
    E56.atomic_write_json(CERTIFICATE, payload)
    if not exact:
        return 1
    validate_exact_certificate_payload(payload)
    print(json.dumps(payload["verdict"], indent=1))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("run", "promote", "validate"))
    parser.add_argument(
        "--groups", default="survivors,undecided,no_reference"
    )
    parser.add_argument("--indexes", default="")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.command == "run":
        groups = [token for token in args.groups.split(",") if token]
        indexes = (
            {int(token) for token in args.indexes.split(",") if token}
            if args.indexes else None
        )
        return run(groups, indexes)
    if args.command == "promote":
        return promote(force=args.force)
    payload = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    validate_exact_certificate_payload(payload)
    print(json.dumps(payload["verdict"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
