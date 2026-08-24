"""EXP-057: exact monolithic certificate for the new odd-lattice [[170,16,10]].

EXP-055's first screen above n=162 leaves three k=16 classes on the 17x5
lattice without an admissible reference.  Two have immediate weight-2/4
logicals.  The remaining constructor is exactified here.  Odd check-column
degree reduces exclusion through weight 9 to one weight-8 Z-logical query;
the BB block-swap/inversion isometry then gives d_X=d_Z.  The decisive UNSAT
query is replayed and hash-bound to the current monolithic SAT encoding.

Run:
  python experiments/exp057_odd_frontier.py run
  python experiments/exp057_odd_frontier.py validate
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


E56 = _load("exp056_for_exp057", "exp056_odd_distance.py")

from qec_research.distance.sat_decide import (  # noqa: E402
    ENCODING_VERSION,
    css_side_instance,
    decide_weight_bounded,
    decision_cnf_digest,
)
from qec_research.gf2.linalg import rank_np  # noqa: E402

SCHEMA = "exp057-odd-frontier-distance-v1"
VALIDATION_VERSION = "exp057-pure-rebuild-validator-v1"
SOLVER_NAME = "cadical195"
CERTIFICATE = ROOT / "results" / "certificates" / "exp057_170_16_10_distance.json"
DISCOVERY_SHARD = ROOT / "results" / "partial_runs" / "exp055_screen" / "17x5.json"
TARGET: dict[str, Any] = {
    "name": "exp055-17x5-frontier-170-16-10",
    "source": "EXP-055 exhaustive weight-3 odd-lattice screen",
    "ell": 17,
    "m": 5,
    "expected_k": 16,
    "expected_d": 10,
    "lower_cap": 9,
    "A": [[0, 0], [1, 0], [9, 1]],
    "B": [[0, 0], [2, 1], [15, 1]],
    "witness_support": [2, 4, 7, 9, 14, 43, 46, 53, 84, 85],
}


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
        raise RuntimeError("EXP-057 requires the proved even-kernel reduction")
    return {
        "validation_version": VALIDATION_VERSION,
        "encoding_version": ENCODING_VERSION,
        "decomposition": "monolithic",
        "side": "z",
        "solver": SOLVER_NAME,
        "conflict_budget": 0,
        "requested_cap": int(TARGET["lower_cap"]),
        "effective_even_cap": int(parity["effective_even_cap"]),
        "replay_required": True,
    }


def _lower_instance(problem: dict[str, Any]):
    return css_side_instance(
        problem["HX"], problem["HZ"], "z", block_length=problem["block"]
    )


def _decision_valid(record: dict[str, Any], expected_digest: str, cap: int) -> bool:
    solver = record.get("solver", {})
    return bool(
        record.get("kind") == "css_z"
        and int(record.get("weight_cap", -1)) == cap
        and record.get("status") == "UNSAT"
        and record.get("cnf_sha256") == expected_digest
        and record.get("encoding_version") == ENCODING_VERSION
        and record.get("symmetry_break") is True
        and solver.get("name") == f"PySAT {SOLVER_NAME}"
        and solver.get("backend") == "pysat"
        and int(solver.get("conflict_budget", -1)) == 0
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
            and _decision_valid(lower["initial"], expected_digest, cap)
            and _decision_valid(lower["replay"], expected_digest, cap)
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
        raise RuntimeError("EXP-057 exact certificate does not validate")


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


def run(*, force: bool = False) -> int:
    if CERTIFICATE.exists() and not force:
        certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
        validate_exact_certificate_payload(certificate)
        print(json.dumps(certificate["verdict"], indent=1))
        return 0

    # Revoke any older exact claim before starting a fresh protocol execution.
    E56.atomic_write_json(CERTIFICATE, _tombstone("RUNNING"))
    problem = build_problem()
    protocol = _protocol(problem)
    instance = _lower_instance(problem)
    cap = int(protocol["effective_even_cap"])
    initial = decide_weight_bounded(
        instance, cap, solver_name=SOLVER_NAME, conflict_budget=0
    )
    if initial["status"] != "UNSAT":
        E56.atomic_write_json(
            CERTIFICATE,
            _tombstone("LOWER_BOUND_REFUTED", lower_bound={"initial": initial}),
        )
        return 1

    replay = decide_weight_bounded(
        instance, cap, solver_name=SOLVER_NAME, conflict_budget=0
    )
    parity = E56.kernel_weight_parity_proof(problem, int(TARGET["lower_cap"]))
    duality = E56.bb_duality_proof(problem)
    upper = E56.verify_upper_witness(problem, E56.default_witness(problem))
    exact = bool(
        replay["status"] == "UNSAT"
        and initial["cnf_sha256"] == replay["cnf_sha256"]
        and upper["weight"] == int(TARGET["expected_d"])
        and cap + 2 == int(TARGET["expected_d"])
        and parity["valid"]
        and duality["valid"]
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
            "screen_verdict": "no_reference",
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
