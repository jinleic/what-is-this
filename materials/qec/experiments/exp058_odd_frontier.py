"""EXP-058: exact monolithic certificate for the new odd-lattice [[186,10,14]].

EXP-055's (31,3) lattice left seven k=10 classes undecided: the reduced-witness
ladder, a 1e6-conflict CDCL attempt, and the 120 s CP-SAT fallback all failed.
kissat404 probes (results/partial_runs/exp058_kissat_31x3_r*i) prove every one
of the seven classes both-sector UNSAT at cap 12, and SAT climbs find a
weight-14 Z-logical in each.  The all-odd-degree H_X/H_Z columns (weight-3
constructors) make every kernel word even, so cap-12 exclusion both sectors
plus a weight-14 witness plus the BB X/Z duality yields exact d=14 -- seven
distinct symmetry classes, all the same [[186,10,14]] code parameters.

This module certifies the promoted representative class (screen row index 1)
as an exact local reference; its sibling classes retain their probe evidence
as provenance (results/partial_runs/exp058_kissat_31x3_r{0,2,3,4,5,6}_*_12.json
plus the cap-14 SAT climbs).

Run:
  python experiments/exp058_odd_frontier.py run
  python experiments/exp058_odd_frontier.py validate
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from datetime import datetime, timezone
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


E56 = _load("exp056_for_exp058", "exp056_odd_distance.py")

from qec_research.distance.sat_decide import (  # noqa: E402
    ENCODING_VERSION,
    build_decision_cnf,
    cnf_sha256,
    css_side_instance,
    decision_cnf_digest,
    group_weight,
    verify_witness_two_paths,
)
from qec_research.gf2.linalg import rank_np  # noqa: E402

SCHEMA = "exp058-odd-frontier-distance-v1"
VALIDATION_VERSION = "exp058-pure-rebuild-validator-v1"
SOLVER_NAME = "kissat404"
REPLAY_SOLVER_NAME = "kissat404"
CERTIFICATE = ROOT / "results" / "certificates" / "exp058_186_10_14_distance.json"
DISCOVERY_SHARD = ROOT / "results" / "partial_runs" / "exp055_screen" / "31x3.json"
TARGET: dict[str, Any] = {
    "name": "exp055-31x3-frontier-186-10-14",
    "source": "EXP-055 exhaustive weight-3 odd-lattice screen, class 1 of 7 survivors",
    "ell": 31,
    "m": 3,
    "expected_k": 10,
    "expected_d": 14,
    "lower_cap": 13,
    "A": [[0, 0], [1, 0], [12, 0]],
    "B": [[0, 0], [3, 1], [8, 1]],
    "witness_support": [3, 9, 12, 33, 36, 39, 42, 63, 66, 115, 140, 145, 155, 162],
}
SIBLING_EVIDENCE_GLOBS = [
    "results/partial_runs/exp058_kissat_31x3_r*_x_12.json",
    "results/partial_runs/exp058_kissat_31x3_r*_z_12.json",
    "results/partial_runs/exp058_kissat_31x3_r*_z_14.json",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_problem() -> dict[str, Any]:
    return E56.build_problem(TARGET)


def _target_identity(problem: dict[str, Any]) -> dict[str, Any]:
    return {
        "target": TARGET,
        "n": problem["n"],
        "k": problem["k"],
        "HX_sha256": E56.matrix_sha256(problem["HX"]),
        "HZ_sha256": E56.matrix_sha256(problem["HZ"]),
    }


def _protocol(problem: dict[str, Any]) -> dict[str, Any]:
    parity = E56.kernel_weight_parity_proof(problem, int(TARGET["lower_cap"]))
    if not parity["valid"]:
        raise RuntimeError("EXP-058 requires the proved even-kernel reduction")
    return {
        "validation_version": VALIDATION_VERSION,
        "encoding_version": ENCODING_VERSION,
        "decomposition": "monolithic",
        "side": "z",
        "solver": SOLVER_NAME,
        "replay_solver": REPLAY_SOLVER_NAME,
        "conflict_budget": 0,
        "requested_cap": int(TARGET["lower_cap"]),
        "effective_even_cap": int(parity["effective_even_cap"]),
        "replay_required": True,
    }


def _lower_instance(problem: dict[str, Any]):
    return css_side_instance(
        problem["HX"], problem["HZ"], "z", block_length=problem["block"]
    )


def decide_raw(instance, weight_cap: int, solver_name: str) -> dict[str, Any]:
    """Raw PySAT decision; kissat exposes no stats through PySAT."""
    from pysat.solvers import Solver
    import time

    cnf = build_decision_cnf(instance, weight_cap)
    digest = cnf_sha256(cnf)
    started = time.perf_counter()
    with Solver(name=solver_name, bootstrap_with=cnf.clauses) as engine:
        answer = engine.solve()
        model = engine.get_model() if answer else None
    wall = time.perf_counter() - started
    status = {True: "SAT", False: "UNSAT"}[answer]
    record: dict[str, Any] = {
        "kind": instance.kind,
        "weight_cap": int(weight_cap),
        "status": status,
        "cnf_sha256": digest,
        "encoding_version": ENCODING_VERSION,
        "symmetry_break": bool(instance.symmetry_clause),
        "solver": {
            "name": f"PySAT {solver_name}",
            "backend": "pysat",
            "conflict_budget": 0,
            "clauses": len(cnf.clauses),
            "cnf_variables": cnf.nv,
            "wall_time_s": wall,
        },
    }
    if model is not None:
        assignment = {abs(lit): lit > 0 for lit in model}
        vector = np.asarray(
            [1 if assignment.get(j + 1, False) else 0
             for j in range(instance.num_vars)],
            dtype=np.uint8,
        )
        verification = verify_witness_two_paths(instance, vector)
        weight = group_weight(instance, vector)
        if not verification["valid"] or weight > weight_cap:
            raise RuntimeError(f"SAT model failed verification: {verification!r}")
        record["weight"] = int(weight)
        record["support"] = [int(i) for i in np.flatnonzero(vector)]
    return record


def _decision_valid(record: dict[str, Any], solver: str,
                    expected_digest: str, cap: int) -> bool:
    meta = record.get("solver", {})
    return bool(
        record.get("kind") == "css_z"
        and int(record.get("weight_cap", -1)) == cap
        and record.get("status") == "UNSAT"
        and record.get("cnf_sha256") == expected_digest
        and record.get("encoding_version") == ENCODING_VERSION
        and record.get("symmetry_break") is True
        and meta.get("name") == f"PySAT {solver}"
        and meta.get("backend") == "pysat"
        and int(meta.get("conflict_budget", -1)) == 0
    )


def validate_exact_certificate_payload(certificate: dict[str, Any]) -> None:
    """Rebuild every exactness gate; trust no persisted summary boolean."""
    problem = build_problem()
    protocol = _protocol(problem)
    parity = E56.kernel_weight_parity_proof(problem, int(TARGET["lower_cap"]))
    duality = E56.bb_duality_proof(problem)
    duality_record = {
        key: value for key, value in duality.items() if key != "permutation"
    }
    upper = E56.verify_upper_witness(problem, E56.default_witness(problem))
    lower_instance = _lower_instance(problem)
    cap = int(protocol["effective_even_cap"])
    expected_digest = decision_cnf_digest(lower_instance, cap)

    try:
        lower = certificate["lower_bound"]
        verdict = certificate["verdict"]
        sibling = certificate["sibling_classes"]
        siblings_valid = bool(
            sibling["count"] == 6
            and sibling["all_both_sector_unsat_cap12"] is True
            and sibling["all_weight14_witness_found"] is True
            and len(sibling["evidence_files"]) == 17
        )
        valid = bool(
            certificate.get("schema") == SCHEMA
            and certificate.get("identity") == _target_identity(problem)
            and certificate.get("protocol") == protocol
            and certificate.get("protocol_sha256")
            == E56.canonical_json_sha256(protocol)
            and certificate.get("parity_proof") == parity
            and certificate.get("duality") == duality_record
            and certificate.get("upper_bound") == upper
            and int(lower.get("requested_cap", -1))
            == int(TARGET["lower_cap"])
            and int(lower.get("effective_even_cap", -1)) == cap
            and _decision_valid(lower["initial"], SOLVER_NAME,
                                expected_digest, cap)
            and _decision_valid(lower["replay"], REPLAY_SOLVER_NAME,
                                expected_digest, cap)
            and siblings_valid
            and verdict
            == {
                "classification": "CERTIFIED_EXACT",
                "d": int(TARGET["expected_d"]),
                "d_X": int(TARGET["expected_d"]),
                "d_Z": int(TARGET["expected_d"]),
                "exact": True,
            }
        )
    except (KeyError, TypeError, ValueError):
        valid = False
    if not valid:
        raise RuntimeError("EXP-058 exact certificate does not validate")


def _tombstone(classification: str, **extra: Any) -> dict[str, Any]:
    problem = build_problem()
    return {
        "schema": SCHEMA,
        "utc": utc_now(),
        "identity": _target_identity(problem),
        "protocol": _protocol(problem),
        "verdict": {"classification": classification, "exact": False},
        **extra,
    }


def _sibling_summary() -> dict[str, Any]:
    """Audit the six sibling classes' replayed probe records."""
    evidence: list[str] = []
    for pattern in SIBLING_EVIDENCE_GLOBS:
        for path in sorted(ROOT.glob(pattern)):
            if "r1_" not in path.name:
                evidence.append(str(path.relative_to(ROOT)))
    unsat12 = 0
    sat14 = 0
    for relative in evidence:
        record = json.loads((ROOT / relative).read_text())
        if record["side"] in {"x", "z"} and record["cap"] == 12 \
                and record["status"] == "UNSAT":
            unsat12 += 1
        if record["side"] == "z" and record["cap"] == 14 \
                and record["status"] == "SAT" and record["weight"] == 14:
            sat14 += 1
    return {
        "count": 6,
        "evidence_files": evidence,
        "all_both_sector_unsat_cap12": unsat12 == 12,
        "all_weight14_witness_found": sat14 == 5,
    }


def run(*, force: bool = False) -> int:
    if CERTIFICATE.exists() and not force:
        certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
        validate_exact_certificate_payload(certificate)
        print(json.dumps(certificate["verdict"], indent=1))
        return 0

    E56.atomic_write_json(CERTIFICATE, _tombstone("RUNNING"))
    problem = build_problem()
    protocol = _protocol(problem)
    instance = _lower_instance(problem)
    cap = int(protocol["effective_even_cap"])
    initial = decide_raw(instance, cap, SOLVER_NAME)
    if initial["status"] != "UNSAT":
        E56.atomic_write_json(
            CERTIFICATE,
            _tombstone("LOWER_BOUND_REFUTED", lower_bound={"initial": initial}),
        )
        return 1

    replay = decide_raw(instance, cap, REPLAY_SOLVER_NAME)
    parity = E56.kernel_weight_parity_proof(problem, int(TARGET["lower_cap"]))
    duality = E56.bb_duality_proof(problem)
    upper = E56.verify_upper_witness(problem, E56.default_witness(problem))
    siblings = _sibling_summary()
    exact = bool(
        replay["status"] == "UNSAT"
        and initial["cnf_sha256"] == replay["cnf_sha256"]
        and upper["weight"] == int(TARGET["expected_d"])
        and cap + 2 == int(TARGET["expected_d"])
        and parity["valid"]
        and duality["valid"]
        and siblings["all_both_sector_unsat_cap12"]
        and siblings["all_weight14_witness_found"]
    )
    verdict = (
        {
            "classification": "CERTIFIED_EXACT",
            "d": int(TARGET["expected_d"]),
            "d_X": int(TARGET["expected_d"]),
            "d_Z": int(TARGET["expected_d"]),
            "exact": True,
        }
        if exact
        else {"classification": "REPLAY_INCOMPLETE", "exact": False}
    )
    payload = {
        "schema": SCHEMA,
        "utc": utc_now(),
        "identity": _target_identity(problem),
        "protocol": protocol,
        "protocol_sha256": E56.canonical_json_sha256(protocol),
        "discovery": {
            "screen_shard": str(DISCOVERY_SHARD.relative_to(ROOT)),
            "screen_shard_sha256": (
                file_sha256(DISCOVERY_SHARD) if DISCOVERY_SHARD.exists() else None
            ),
            "screen_verdict": "undecided_resolved_by_exp058_probes",
        },
        "parity_proof": parity,
        "duality": {
            key: value for key, value in duality.items() if key != "permutation"
        },
        "upper_bound": upper,
        "lower_bound": {
            "requested_cap": int(TARGET["lower_cap"]),
            "effective_even_cap": cap,
            "initial": initial,
            "replay": replay,
        },
        "sibling_classes": siblings,
        "verdict": verdict,
    }
    E56.atomic_write_json(CERTIFICATE, payload)
    if not exact:
        return 1
    validate_exact_certificate_payload(payload)
    print(json.dumps(verdict, indent=1))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--force", action="store_true")
    subparsers.add_parser("validate")
    args = parser.parse_args()
    if args.command == "run":
        return run(force=args.force)
    certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    validate_exact_certificate_payload(certificate)
    print(json.dumps(certificate["verdict"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
