"""EXP-068: persist physical witnesses for EXP-055's 57 fallback dominations.

The original exact CP-SAT fallback retained only ``d_found``.  That value is not
replayable domination evidence.  This repair finds one explicit Z-logical at or
below each threshold on the 29 independent records, transports the 14 n=210
witnesses to the two CRT-equivalent presentations, and makes the shared shard
validator reject every domination without a physical proof.

Run with at most two CP-SAT cores to keep workstation CPU modest:
  python experiments/exp068_screen_proof_repair.py run --workers 2
  python experiments/exp068_screen_proof_repair.py apply
  python experiments/exp068_screen_proof_repair.py validate
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.distance.exact import min_weight_with_parity  # noqa: E402
from qec_research.distance.sat_decide import css_logical_bases  # noqa: E402
from qec_research.gf2.linalg import rank_np  # noqa: E402


def _load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "experiments" / filename
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


E55 = _load("exp055_for_exp068", "exp055_odd_lattice_sweep.py")
E56 = _load("exp056_for_exp068", "exp056_odd_distance.py")
E59 = _load("exp059_for_exp068", "exp059_coprime_transport.py")

SCHEMA = "exp068-screen-fallback-witness-v1"
STATE_DIR = ROOT / "results" / "partial_runs" / "exp068_screen_witnesses"
SOURCE_LATTICES = {(15, 5), (9, 9), (27, 3), (31, 3), (33, 3), (15, 7)}
CRT_TARGETS = ((21, 5), (35, 3))


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record_key(record: dict[str, Any]) -> str:
    return json.dumps(
        {
            "ell": int(record["ell"]),
            "m": int(record["m"]),
            "A": record["A"],
            "B": record["B"],
            "k_parent": int(record["k_parent"]),
        },
        sort_keys=True,
        separators=(",", ":"),
    )


def evidence_path(record: dict[str, Any]) -> Path:
    digest = hashlib.sha256(record_key(record).encode()).hexdigest()[:16]
    return STATE_DIR / f"{int(record['ell'])}x{int(record['m'])}_{digest}.json"


def source_records() -> list[dict[str, Any]]:
    out = []
    for ell, m in sorted(SOURCE_LATTICES):
        shard = json.loads(
            (E55.SCREEN_DIR / f"{ell}x{m}.json").read_text(encoding="utf-8")
        )
        for index, record in enumerate(shard["records"]):
            if record["verdict"] == "dominated":
                out.append({"index": index, "record": record})
    if len(out) != 29:
        raise RuntimeError(f"expected 29 independent fallback records, got {len(out)}")
    return out


def _verify_support(record: dict[str, Any], support: list[int]) -> dict[str, Any]:
    HX, HZ = E55.E53.bb_from_terms(
        int(record["ell"]), int(record["m"]), record["A"], record["B"]
    )
    vector = np.zeros(int(record["n"]), dtype=np.uint8)
    vector[np.asarray(support, dtype=int)] = 1
    return {
        "commutes_with_HX": bool(not np.any(HX @ vector % 2)),
        "outside_Z_stabilizer": bool(
            rank_np(np.vstack([HZ, vector])) == rank_np(HZ) + 1
        ),
        "weight": int(vector.sum()),
        "vector_sha256": E56.matrix_sha256(vector[None, :]),
    }


def validate_evidence(payload: dict[str, Any], record: dict[str, Any]) -> None:
    verification = payload.get("verification", {})
    valid = bool(
        payload.get("schema") == SCHEMA
        and payload.get("identity") == E55._screen_record_identity(record)
        and isinstance(payload.get("support"), list)
        and int(payload.get("weight", -1)) == len(payload["support"])
        and int(payload["weight"]) <= int(record["threshold"])
        and verification == _verify_support(record, payload["support"])
        and verification.get("commutes_with_HX") is True
        and verification.get("outside_Z_stabilizer") is True
        and int(verification.get("weight", -1)) == int(payload["weight"])
    )
    if not valid:
        raise RuntimeError("EXP-068 fallback witness evidence failed validation")


def recover_one(
    item: dict[str, Any], time_limit_s: float, solver_workers: int
) -> dict[str, Any]:
    record = item["record"]
    path = evidence_path(record)
    if path.exists():
        payload = json.loads(path.read_text(encoding="utf-8"))
        validate_evidence(payload, record)
        return payload
    HX, HZ = E55.E53.bb_from_terms(
        int(record["ell"]), int(record["m"]), record["A"], record["B"]
    )
    started = time.perf_counter()
    I = E55.intersect(
        E55.ann_basis(record["A"], int(record["ell"]), int(record["m"])),
        E55.ann_basis(record["B"], int(record["ell"]), int(record["m"])),
    )
    J = I[:, E55.bar_permutation(int(record["ell"]), int(record["m"]))]
    seed = int(hashlib.sha256(record_key(record).encode()).hexdigest()[:16], 16)
    reduced = E55.reduced_witness_bound(
        J,
        HZ,
        int(record["ell"]),
        int(record["m"]),
        tries=2_000,
        keep=64,
        seed=seed,
        target=int(record["threshold"]),
    )
    if (
        reduced.get("bound") is not None
        and int(reduced["bound"]) <= int(record["threshold"])
    ):
        support = [int(index) for index in reduced["witness_support"]]
        payload = {
            "schema": SCHEMA,
            "identity": E55._screen_record_identity(record),
            "route": "deterministic deep information-set coset reduction",
            "search_side": "Z",
            "sector": None,
            "status": "WITNESS",
            "weight": len(support),
            "support": support,
            "attempts": [
                {
                    "route": "information_set",
                    "tries": 2_000,
                    "classes_tested": int(reduced["classes_tested"]),
                    "seed": seed,
                }
            ],
            "protocol": {
                "tries": 2_000,
                "keep": 64,
                "seed": seed,
                "upper_bound": int(record["threshold"]),
            },
            "verification": _verify_support(record, support),
            "wall_time_s": time.perf_counter() - started,
        }
        validate_evidence(payload, record)
        E56.atomic_write_json(path, payload)
        return payload
    LX, LZ = css_logical_bases(HX, HZ)
    n = int(record["n"])
    groups = [[index] for index in range(n)]
    problem = E56.build_problem(
        {
            "ell": int(record["ell"]),
            "m": int(record["m"]),
            "A": record["A"],
            "B": record["B"],
            "expected_k": int(record["k_parent"]),
        }
    )
    duality = E56.bb_duality_proof(problem)
    permutation = np.eye(n, dtype=np.uint8)[duality["permutation"]]
    # CP-SAT is a fail-closed fallback when the deterministic reduction misses.
    started = time.perf_counter()
    attempts = []
    searches = (("X", HZ, LZ), ("Z", HX, LX))
    for search_side, checks, detectors in searches:
        for sector, detector in enumerate(detectors):
            value, lower, status, support = min_weight_with_parity(
                checks,
                groups,
                detector,
                n,
                time_limit_s=float(time_limit_s),
                workers=int(solver_workers),
                upper_bound=int(record["threshold"]),
            )
            attempts.append(
                {
                    "side": search_side,
                    "sector": sector,
                    "status": status,
                    "value": value,
                    "lower_bound": lower,
                }
            )
            if support is None:
                continue
            vector = np.zeros(n, dtype=np.uint8)
            vector[np.asarray(support, dtype=int)] = 1
            if search_side == "X":
                vector = (vector @ permutation.T) % 2
            support = [int(index) for index in np.flatnonzero(vector)]
            route = (
                "single-sector CP-SAT X witness plus exact BB duality"
                if search_side == "X"
                else "single-sector CP-SAT Z witness feasibility"
            )
            payload = {
                "schema": SCHEMA,
                "identity": E55._screen_record_identity(record),
                "route": route,
                "search_side": search_side,
                "sector": sector,
                "status": status,
                "weight": len(support),
                "support": support,
                "attempts": attempts,
                "protocol": {
                    "time_limit_s_per_sector": float(time_limit_s),
                    "workers": int(solver_workers),
                    "upper_bound": int(record["threshold"]),
                },
                "verification": _verify_support(record, support),
                "wall_time_s": time.perf_counter() - started,
            }
            validate_evidence(payload, record)
            E56.atomic_write_json(path, payload)
            return payload
    raise RuntimeError(
        f"no fallback witness recovered for {record_key(record)}: {attempts}"
    )


def run(workers: int, time_limit_s: float) -> int:
    if workers < 1 or workers > 2:
        raise ValueError("EXP-068 permits one or two total CP-SAT workers")
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    items = source_records()
    # The original fallback used two workers inside each model. Preserve that
    # search behaviour while limiting the whole repair to at most two cores.
    with ThreadPoolExecutor(max_workers=1) as pool:
        results = list(
            pool.map(
                lambda item: recover_one(item, time_limit_s, workers), items
            )
        )
    print(
        json.dumps(
            {
                "records": len(results),
                "max_total_cpsat_workers": workers,
                "all_physical": all(
                    result["verification"]["commutes_with_HX"]
                    and result["verification"]["outside_Z_stabilizer"]
                    for result in results
                ),
            },
            indent=1,
        )
    )
    return 0


def _bind_record(record: dict[str, Any], payload: dict[str, Any]) -> dict[str, Any]:
    validate_evidence(payload, record)
    if E55._fallback_witness_evidence_valid(record):
        return json.loads(json.dumps(record))
    output = json.loads(json.dumps(record))
    output["pre_exp068_reduced_witness"] = {
        "bound": record.get("witness_bound"),
        "support": record.get("witness_support"),
    }
    output["witness_bound"] = int(payload["weight"])
    output["witness_support"] = payload["support"]
    path = evidence_path(output)
    output["fallback_witness"] = {
        "schema": SCHEMA,
        "path": str(path.relative_to(ROOT)),
        "sha256": file_sha256(path),
        "route": payload["route"],
        "sector": payload["sector"],
    }
    return output


def _rewrite_shard(
    ell: int,
    m: int,
    replacements: dict[str, dict],
    metadata: dict[str, Any] | None = None,
) -> dict:
    path = E55.SCREEN_DIR / f"{ell}x{m}.json"
    shard = json.loads(path.read_text(encoding="utf-8"))
    records = [
        replacements.get(record_key(record), record) for record in shard["records"]
    ]
    verdicts: dict[str, int] = {}
    for record in records:
        verdicts[record["verdict"]] = verdicts.get(record["verdict"], 0) + 1
    rebound = {
        **shard,
        **(metadata or {}),
        "records": records,
        "verdicts": verdicts,
        "survivors": [r for r in records if r["verdict"] == "survivor"],
        "no_reference": [r for r in records if r["verdict"] == "no_reference"],
        "undecided": [r for r in records if r["verdict"] == "undecided"],
        "exp068_proof_repair": {
            "schema": SCHEMA,
            "records_rebound": len(replacements),
            "all_domination_witnesses_physical": True,
        },
    }
    E55._validate_screen_shard_aggregates(rebound)
    E55._validate_screen_shard_records(rebound)
    E56.atomic_write_json(path, rebound)
    return rebound


def apply() -> int:
    source_replacements: dict[tuple[int, int], dict[str, dict]] = {}
    for item in source_records():
        record = item["record"]
        payload = json.loads(evidence_path(record).read_text(encoding="utf-8"))
        bound = _bind_record(record, payload)
        source_replacements.setdefault(
            (int(record["ell"]), int(record["m"])), {}
        )[record_key(record)] = bound
    rewritten = {
        lattice: _rewrite_shard(*lattice, replacements)
        for lattice, replacements in source_replacements.items()
    }

    source_lattice = (15, 7)
    source_shard = rewritten[source_lattice]
    for target in CRT_TARGETS:
        mapping = E59.coordinate_transport(*source_lattice, *target)
        target_path = E55.SCREEN_DIR / f"{target[0]}x{target[1]}.json"
        target_shard = json.loads(target_path.read_text(encoding="utf-8"))
        target_by_key = {
            record_key(record): record for record in target_shard["records"]
        }
        replacements = {}
        for source_record in source_shard["records"]:
            if source_record["verdict"] != "dominated":
                continue
            transported = E59._transport_record(
                source_record, source_lattice, target, mapping
            )
            target_record = target_by_key.get(record_key(transported))
            if target_record is None:
                raise RuntimeError("CRT fallback witness target identity missing")
            support = transported["witness_support"]
            payload = {
                "schema": SCHEMA,
                "identity": E55._screen_record_identity(target_record),
                "route": "exact CRT transport of CP-SAT witness",
                "sector": source_record["fallback_witness"]["sector"],
                "status": "TRANSPORTED",
                "weight": len(support),
                "support": support,
                "source": source_record["fallback_witness"],
                "transport": transported["transported_from"],
                "verification": _verify_support(target_record, support),
                "wall_time_s": 0.0,
            }
            path = evidence_path(target_record)
            validate_evidence(payload, target_record)
            E56.atomic_write_json(path, payload)
            replacements[record_key(target_record)] = _bind_record(
                target_record, payload
            )
        source_path = E55.SCREEN_DIR / f"{source_lattice[0]}x{source_lattice[1]}.json"
        transport_metadata = {
            **target_shard["transport"],
            "source_shard_sha256": file_sha256(source_path),
            "mapping_sha256": E59.matrix_sha256(mapping[None, :]),
            "records_transported": len(source_shard["records"]),
            "physical_witnesses_rechecked": sum(
                record["verdict"].startswith("dominated")
                for record in source_shard["records"]
            ),
            "all_physical_witnesses_valid": True,
            "valid": True,
        }
        _rewrite_shard(
            *target, replacements, metadata={"transport": transport_metadata}
        )

    old_screen = json.loads(E55.SCREEN_OUT.read_text(encoding="utf-8"))
    lattices = ",".join(
        f"{int(shard['ell'])}x{int(shard['m'])}"
        for shard in old_screen["lattices"]
    )
    E63 = _load("exp063_for_exp068", "exp063_reference_rebind.py")
    E63.run(lattices.split(","))
    args = argparse.Namespace(
        k_min=8, k_max=24, time_limit=120.0, lattices=lattices, min_n=18
    )
    E55.assemble_screen(args)
    E55.certify_screen_survivors(args)
    E67 = _load("exp067_for_exp068", "exp067_n234_connected_cluster.py")
    E66 = _load("exp066_for_exp068", "exp066_n234_frontier.py")
    E67._refresh_exp066_summary(E55, E66)
    return validate()


def validate() -> int:
    screen = json.loads(E55.SCREEN_OUT.read_text(encoding="utf-8"))
    repaired = 0
    for shard in screen["lattices"]:
        E55._validate_screen_shard_aggregates(shard)
        E55._validate_screen_shard_records(shard)
        for record in shard["records"]:
            if record["verdict"] == "dominated":
                if not E55._fallback_witness_evidence_valid(record):
                    raise RuntimeError("fallback witness binding is stale")
                repaired += 1
    if repaired != 57:
        raise RuntimeError(f"expected 57 repaired records, got {repaired}")
    print(json.dumps({"repaired_records": repaired, "valid": True}, indent=1))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("run", "apply", "validate"))
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--time-limit", type=float, default=120.0)
    args = parser.parse_args()
    if args.command == "run":
        return run(args.workers, args.time_limit)
    if args.command == "apply":
        return apply()
    return validate()


if __name__ == "__main__":
    raise SystemExit(main())
