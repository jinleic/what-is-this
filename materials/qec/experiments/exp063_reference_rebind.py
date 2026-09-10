"""EXP-063: proof-preserving reference-set rebinding for EXP-055 shards.

Adding exact references monotonically raises domination thresholds.  Therefore:

* every prior domination proof remains valid without rerunning a solver;
* a prior survivor/undecided/no-reference row becomes dominated when an already
  verified physical witness is at or below the new threshold;
* otherwise the row stays survivor/undecided and is never silently promoted.

Each rewrite archives the original shard, recomputes every threshold under the
current validator, rechecks every carried physical witness, and records the old
shard hash plus per-row transition counts.  CNF bytes are untouched.

Run:
  python experiments/exp063_reference_rebind.py run --lattices 3x3,...,15x7
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import shutil
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


E55 = _load("exp055_for_exp063", "exp055_odd_lattice_sweep.py")
E56 = _load("exp056_for_exp063", "exp056_odd_distance.py")

SCHEMA = "exp063-monotone-reference-rebind-v1"
SCREEN_DIR = ROOT / "results" / "partial_runs" / "exp055_screen"
ARCHIVE_DIR = ROOT / "results" / "partial_runs" / "exp063_reference_rebind"
RATCHET_DIR = ROOT / "results" / "partial_runs" / "exp060_n210"


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record_key(record: dict[str, Any]) -> str:
    return json.dumps(
        {"ell": record["ell"], "m": record["m"],
         "A": record["A"], "B": record["B"]},
        sort_keys=True, separators=(",", ":"),
    )


def ratchet_records() -> dict[str, dict[str, Any]]:
    out = {}
    if not RATCHET_DIR.exists():
        return out
    for path in sorted(RATCHET_DIR.glob("*.json")):
        if "_cap" in path.stem and "_probe_" not in path.stem:
            continue
        record = json.loads(path.read_text(encoding="utf-8"))
        if all(key in record for key in ("A", "B", "n", "k")):
            key = json.dumps(
                {"ell": record["n"] // (2 * 7), "m": 7,
                 "A": record["A"], "B": record["B"]},
                sort_keys=True, separators=(",", ":"),
            )
            out[key] = record
    return out


def _ratchet_witness(result: dict[str, Any]) -> tuple[int, list[int]] | None:
    if isinstance(result.get("decision"), dict):
        decision = result["decision"]
        if decision.get("status") == "SAT":
            return (
                int(decision["weight"]),
                [int(index) for index in decision["support"]],
            )
    distance = result.get("exact_distance")
    if not result.get("exact") or distance is None:
        return None
    for decision in result.get("decisions", []):
        if decision.get("status") == "SAT" and int(decision.get("weight", -1)) == int(distance):
            return int(distance), [int(index) for index in decision["support"]]
    upper = result.get("upper", {})
    if int(upper.get("weight", -1)) == int(distance):
        return int(distance), [int(index) for index in upper["support"]]
    return None


def _physical_witness_valid(record: dict[str, Any]) -> bool:
    support = record.get("witness_support")
    if support is None:
        return False
    HX, HZ = E55.E53.bb_from_terms(
        int(record["ell"]), int(record["m"]), record["A"], record["B"]
    )
    vector = np.zeros(int(record["n"]), dtype=np.uint8)
    vector[np.asarray(support, dtype=int)] = 1
    return bool(
        int(vector.sum()) == int(record["witness_bound"])
        and not np.any(HX @ vector % 2)
        and E55.rank_np(np.vstack([HZ, vector])) == E55.rank_np(HZ) + 1
    )

def _ceiling_witness_valid(record: dict[str, Any]) -> bool:
    support = record.get("ceiling_witness_support")
    if support is None:
        return False
    HX, HZ = E55.E53.bb_from_terms(
        int(record["ell"]), int(record["m"]), record["A"], record["B"]
    )
    vector = np.zeros(int(record["n"]), dtype=np.uint8)
    vector[np.asarray(support, dtype=int)] = 1
    return bool(
        int(vector.sum()) == int(record["ceiling"])
        and not np.any(HX @ vector % 2)
        and E55.rank_np(np.vstack([HZ, vector])) == E55.rank_np(HZ) + 1
    )
def _validate_bound_archive(payload: dict[str, Any]) -> None:
    binding = payload.get("reference_rebind")
    if not isinstance(binding, dict):
        return
    relative = binding.get("archive")
    if not isinstance(relative, str):
        raise RuntimeError("reference-rebind archive path is missing")
    archive = ROOT / relative
    if (
        not archive.is_file()
        or file_sha256(archive) != binding.get("old_shard_sha256")
    ):
        raise RuntimeError("reference-rebind archive hash is stale")
    archived = json.loads(archive.read_text(encoding="utf-8"))
    E55._validate_screen_shard_aggregates(archived)
    E55._validate_screen_shard_records(archived, allow_stale_thresholds=True)




def rebind_shard(path: Path, ratchets: dict[str, dict[str, Any]]) -> dict[str, Any]:
    original_bytes = path.read_bytes()
    old_shard_sha256 = hashlib.sha256(original_bytes).hexdigest()
    payload = json.loads(original_bytes)
    E55._validate_screen_shard_aggregates(payload)
    # The original shard's thresholds are pinned to the previous battery by
    # construction; repairing exactly that staleness is this tool's purpose.
    # Every identity, witness, and verdict-structure check still runs here,
    # and the rebound shard below is validated strictly against the live
    # battery before it is written.
    E55._validate_screen_shard_records(payload, allow_stale_thresholds=True)
    _validate_bound_archive(payload)
    old_protocol = payload["protocol"]
    k_min, k_max = old_protocol["k_range"]
    new_protocol = E55._screen_protocol(
        old_protocol["census_sha256"], E55._reference_fingerprint(),
        int(k_min), int(k_max), float(old_protocol["time_limit_s"]),
    )
    if old_protocol == new_protocol:
        E55._validate_screen_shard_records(payload)
        return payload

    transitions: dict[str, int] = {}
    records = []
    for source in payload["records"]:
        record = json.loads(json.dumps(source))
        old_threshold = int(record.get("threshold", 0))
        new_threshold, new_source = E55.domination_threshold(
            int(record["n"]), int(record["k_parent"])
        )
        if new_threshold < old_threshold:
            raise RuntimeError("reference rebinding is not monotone")
        record["threshold"] = int(new_threshold)
        record["threshold_source"] = new_source
        old_verdict = record["verdict"]
        new_verdict = old_verdict

        if old_verdict.startswith("dominated"):
            pass
        else:
            ceiling = record.get("ceiling")
            witness = record.get("witness_bound")
            result = ratchets.get(record_key(record))
            exact_witness = _ratchet_witness(result) if result else None
            if ceiling is not None and int(ceiling) <= new_threshold:
                if not _ceiling_witness_valid(record):
                    raise RuntimeError("rebound ceiling witness failed verification")
                new_verdict = "dominated_by_ceiling"
            elif witness is not None and int(witness) <= new_threshold:
                if not _physical_witness_valid(record):
                    raise RuntimeError("rebound physical witness failed verification")
                new_verdict = "dominated_by_witness"
            elif exact_witness and exact_witness[0] <= new_threshold:
                record["witness_bound"] = exact_witness[0]
                record["witness_support"] = exact_witness[1]
                record["reference_rebind_ratchet"] = {
                    "schema": result["schema"],
                    "exact_distance": result.get("exact_distance"),
                    "witness_only": result.get("exact_distance") is None,
                    "source_group": result["group"],
                    "source_index": result["index"],
                    "solver": result.get("solver"),
                }
                new_verdict = "dominated_by_cdcl_witness"
            elif new_threshold > old_threshold:
                new_verdict = (
                    "no_reference" if new_threshold == 0 else "undecided"
                )
        record["verdict"] = new_verdict
        transition = f"{old_verdict}->{new_verdict}"
        transitions[transition] = transitions.get(transition, 0) + 1
        if new_verdict == "dominated_by_ceiling":
            if not _ceiling_witness_valid(record):
                raise RuntimeError("rebound ceiling witness failed verification")
        elif new_verdict in {"dominated_by_witness", "dominated_by_cdcl_witness"}:
            if not _physical_witness_valid(record):
                raise RuntimeError("rebound physical witness failed verification")
        records.append(record)
    verdicts: dict[str, int] = {}
    for record in records:
        verdicts[record["verdict"]] = verdicts.get(record["verdict"], 0) + 1
    version = E55.REFERENCE_VALIDATION_VERSION
    archive = ARCHIVE_DIR / version / old_shard_sha256 / path.name
    archive.parent.mkdir(parents=True, exist_ok=True)
    if archive.exists():
        archived_bytes = archive.read_bytes()
        archived = json.loads(archived_bytes)
        E55._validate_screen_shard_aggregates(archived)
        E55._validate_screen_shard_records(archived, allow_stale_thresholds=True)
        if archived_bytes != original_bytes:
            raise RuntimeError("reference-rebind archive disagrees with source shard")
    else:
        archive.write_bytes(original_bytes)
    rebound = {
        **payload,
        "protocol": new_protocol,
        "records": records,
        "verdicts": verdicts,
        "survivors": [record for record in records if record["verdict"] == "survivor"],
        "no_reference": [record for record in records if record["verdict"] == "no_reference"],
        "undecided": [record for record in records if record["verdict"] == "undecided"],
        "reference_rebind": {
            "schema": SCHEMA,
            "old_shard_sha256": old_shard_sha256,
            "old_reference_sha256": old_protocol["reference_sha256"],
            "new_reference_sha256": new_protocol["reference_sha256"],
            "old_validation_version": old_protocol["reference_validation_version"],
            "new_validation_version": new_protocol["reference_validation_version"],
            "transitions": transitions,
            "all_thresholds_monotone": True,
            "all_carried_witnesses_rechecked": True,
            "archive": str(archive.relative_to(ROOT)),
        },
    }
    if isinstance(payload.get("transport"), dict):
        source = payload["transport"]["source_lattice"]
        source_path = E55.SCREEN_DIR / f"{source[0]}x{source[1]}.json"
        rebound["transport"] = {
            **payload["transport"],
            "source_shard_sha256": file_sha256(source_path),
        }
    E55._validate_screen_shard_aggregates(rebound)
    E55._validate_screen_shard_records(rebound)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(rebound, indent=1, sort_keys=True) + "\n")
    temporary.replace(path)
    return rebound


def run(lattices: list[str]) -> int:
    ratchets = ratchet_records()
    outputs = []
    for lattice in lattices:
        path = SCREEN_DIR / f"{lattice}.json"
        result = rebind_shard(path, ratchets)
        outputs.append({
            "lattice": lattice,
            "version": result["protocol"]["reference_validation_version"],
            "verdicts": result["verdicts"],
        })
        print(json.dumps(outputs[-1]), flush=True)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("run",))
    parser.add_argument("--lattices", required=True)
    args = parser.parse_args()
    return run([token for token in args.lattices.split(",") if token])


if __name__ == "__main__":
    raise SystemExit(main())
