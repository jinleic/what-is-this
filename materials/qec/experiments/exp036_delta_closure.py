"""EXP-036: SAT closure of the delta>0 parent-domination gap.

EXP-027 left 87 of 155 delta>0 catalogue rows UNDECIDED because parent (and
PBB) distances were uncertifiable within CP-SAT budgets.  EXP-035 proved that
CDCL on equisatisfiable CNFs decides exactly this instance class orders of
magnitude faster (Theorem F campaign: four UNSAT sector proofs in minutes
after days of CP-SAT UNKNOWN).  This experiment turns that method on the
undecided rows.

Frozen comparison protocol (2026-08-15, before any EXP-036 record existed)
---------------------------------------------------------------------------
For a row with parent CSS code P ([[n, k+delta]]) and PBB member B ([[n, k]]):

* ``decide(code, cap)``: SAT => verified witness => d <= weight; UNSAT =>
  d > cap; UNDECIDED_BUDGET => nothing.  Parent decisions are per CSS side;
  the parent is UNSAT at cap iff BOTH sides are UNSAT at cap.
* Maintain verified witness weights U_p (parent), U_b (PBB), initialised
  from untrusted hints by climbing until SAT (an UNSAT at the catalogue
  value is itself a finding: the claimed distance was too low).
* Loop:
  - U_p > U_b: parent UNSAT at U_b - 1  => DOMINATION  (d_p >= U_b >= d_b);
    a SAT descends U_p and the loop continues.
  - U_b > U_p: PBB UNSAT at U_p        => REVERSAL    (d_b > U_p >= d_p);
    a SAT descends U_b and the loop continues.  The cap is U_p, NOT
    U_p - 1: a tie is weak domination, never a reversal, because the parent
    carries more logical qubits at the same length.
  - U_p == U_b = U: parent UNSAT at U - 1 => d_p = U exactly and
    d_b <= U = d_p => DOMINATION; a SAT descends U_p and the loop continues.
* Weak domination (d_p >= d_b) counts as DOMINATION_PROVED, matching
  EXP-027's convention.  Every witness is re-verified through two
  independent GF(2) paths before any bound moves; every UNSAT carries its
  solver statistics and budget.

Records are per-row JSON under ``results/partial_runs/exp036/``, atomic and
resumable; the assembled artifact never mutates the EXP-027 evidence.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.codes.bicycle import build_bb  # noqa: E402
from qec_research.distance.sat_decide import (  # noqa: E402
    css_side_instance,
    decide_by_sectors,
    decide_weight_bounded,
    symplectic_instance,
)
from qec_research.gf2.linalg import rank_np  # noqa: E402

_EXP027_SPEC = importlib.util.spec_from_file_location(
    "exp027_delta_audit", ROOT / "experiments" / "exp027_delta_audit.py"
)
assert _EXP027_SPEC and _EXP027_SPEC.loader
E27 = importlib.util.module_from_spec(_EXP027_SPEC)
_EXP027_SPEC.loader.exec_module(E27)

EXP027_ARTIFACT = ROOT / "results" / "processed" / "exp027_delta_audit.json"
STATE_DIR = ROOT / "results" / "partial_runs" / "exp036"
PARTIAL = ROOT / "results" / "partial_runs" / "exp036_delta_closure.json"
PROCESSED = ROOT / "results" / "processed" / "exp036_delta_closure.json"

SOLVER_NAME = "cadical195"
SCHEMA = "exp036-delta-closure-row-v1"
MAX_CLIMB = 40  # absolute ceiling on witness-search climbing


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # Process-unique temp name: concurrent shards must never share a .tmp
    # path, or interleaved writes can publish a corrupt artifact.
    tmp = path.with_suffix(f"{path.suffix}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(path)


def row_identity(row: dict[str, Any], index: int) -> dict[str, Any]:
    parent_spec, HX, HZ = E27.parent_matrices(row)
    pbb_spec, code = E27.pbb_code(row)
    return {
        "schema": SCHEMA,
        "catalogue_index": index,
        "label": E27.catalogue_label(row, index),
        "code_id": str(row.get("code_id", "")),
        "n": 2 * int(row["ell"]) * int(row["m"]),
        "delta": int(row["delta"]) if "delta" in row else None,
        "parent_fingerprint": E27.matrix_fingerprint(HX, HZ),
        "pbb_fingerprint": E27.matrix_fingerprint(code.H),
        "catalogue_d_untrusted": int(row["d"]),
    }


def row_state_path(index: int) -> Path:
    return STATE_DIR / f"row_{index:04d}.json"


def load_row_state(index: int, identity: dict[str, Any]) -> dict[str, Any]:
    path = row_state_path(index)
    if path.exists():
        state = json.loads(path.read_text(encoding="utf-8"))
        for key in identity:
            if state.get(key) != identity[key]:
                raise RuntimeError(
                    f"row {index}: persisted state identity mismatch on {key}: "
                    f"stored {state.get(key)!r} != rebuilt {identity[key]!r} - "
                    "nothing may resume or assemble until this is investigated"
                )
        return state
    state = dict(identity)
    state["calls"] = []
    state["verdict"] = "PENDING"
    state["witnesses"] = {}
    return state


def verified_row_state(index: int, row: dict[str, Any]) -> dict[str, Any] | None:
    """Load a persisted row state ONLY if it matches the rebuilt identity.

    Returns None when no state exists; fail-stops on any identity mismatch.
    """

    if not row_state_path(index).exists():
        return None
    return load_row_state(index, row_identity(row, index))


def record_call(
    state: dict[str, Any],
    *,
    target: str,
    record: dict[str, Any],
) -> None:
    entry = {
        "target": target,
        "kind": record["kind"],
        "weight_cap": record["weight_cap"],
        "status": record["status"],
        "cnf_sha256": record.get("cnf_sha256"),
        "encoding_version": record.get("encoding_version"),
        "decomposition": record.get("decomposition", "monolithic"),
        "solver": record["solver"],
        "utc": utc_now(),
    }
    if record.get("vector") is not None:
        entry["weight"] = record["weight"]
        entry["verification"] = record["verification"]
        state["witnesses"][f"{target}_w{record['weight']}"] = {
            "vector": record["vector"],
            "weight": record["weight"],
            "verification": record["verification"],
        }
    state["calls"].append(entry)


def run_decision(
    instance,
    cap: int,
    *,
    conflict_budget: int,
    decomposition: str,
):
    """Dispatch one decision under an explicitly recorded encoding mode.

    ``monolithic`` keeps the original single-CNF disjunction encoding;
    ``sectors`` uses the exact per-functional decomposition.  The mode is
    stored on every call so a replay rebuilds the identical CNF(s) and hash.
    """

    if decomposition == "sectors":
        return decide_by_sectors(
            instance,
            cap,
            solver_name=SOLVER_NAME,
            conflict_budget=conflict_budget,
        )
    if decomposition == "monolithic":
        return decide_weight_bounded(
            instance,
            cap,
            solver_name=SOLVER_NAME,
            conflict_budget=conflict_budget,
        )
    raise ValueError(f"unknown decomposition {decomposition!r}")


class BudgetExhausted(Exception):
    pass


def build_instances(row: dict[str, Any]) -> dict[str, Any]:
    """Instances with the sound BB translation symmetry break enabled.

    Every catalogue row is bivariate-bicycle derived with block length
    lm = ell*m, so the translation group acts transitively on each block and
    the single-clause symmetry break is satisfiability-preserving.
    """

    _, HX, HZ = E27.parent_matrices(row)
    _, pbb = E27.pbb_code(row)
    block = int(row["ell"]) * int(row["m"])
    return {
        "parent_x": css_side_instance(HX, HZ, "x", block_length=block),
        "parent_z": css_side_instance(HX, HZ, "z", block_length=block),
        "pbb": symplectic_instance(pbb, pbb.logical_basis(), block_length=block),
    }


def validated_bounds(
    instances: dict[str, Any], state: dict[str, Any]
) -> dict[str, Any]:
    """Re-derive witness uppers and UNSAT lowers from re-verified evidence.

    Nothing stored is trusted: every witness vector is re-checked through two
    independent GF(2) paths against freshly rebuilt instances, and UNSAT
    lower bounds are collected from call statuses (their proof-grade replay
    happens in ``replay_decisive_unsats``).
    """

    from qec_research.distance.sat_decide import (
        group_weight,
        verify_witness_two_paths,
    )

    upper: dict[str, int | None] = {"parent": None, "pbb": None}
    target_witness_min: dict[str, int] = {}
    for key, witness in state.get("witnesses", {}).items():
        target = key.rsplit("_w", 1)[0]
        instance = instances[target]
        vector = np.asarray(witness["vector"], dtype=np.uint8)
        verification = verify_witness_two_paths(instance, vector)
        weight = group_weight(instance, vector)
        if not verification["valid"] or weight != int(witness["weight"]):
            raise RuntimeError(
                f"stored witness {key} failed reverification "
                f"(recomputed weight {weight}, valid={verification['valid']})"
            )
        side = "parent" if target.startswith("parent") else "pbb"
        if upper[side] is None or weight < upper[side]:
            upper[side] = weight
        if target not in target_witness_min or weight < target_witness_min[target]:
            target_witness_min[target] = weight
    unsat_caps: dict[str, list[int]] = {"parent_x": [], "parent_z": [], "pbb": []}
    for entry in state.get("calls", []):
        if entry["status"] == "UNSAT":
            target = entry["target"]
            cap = int(entry["weight_cap"])
            if target in target_witness_min and cap >= target_witness_min[target]:
                raise RuntimeError(
                    f"internal contradiction on {target}: a verified witness "
                    f"of weight {target_witness_min[target]} coexists with a "
                    f"stored UNSAT at cap {cap} - the call log is defective "
                    "or tampered"
                )
            unsat_caps[target].append(cap)
    parent_lower_gt = None
    both = set(unsat_caps["parent_x"]) & set(unsat_caps["parent_z"])
    if both:
        parent_lower_gt = max(both)
    pbb_lower_gt = max(unsat_caps["pbb"]) if unsat_caps["pbb"] else None
    return {
        "U_p": upper["parent"],
        "U_b": upper["pbb"],
        "parent_lower_gt": parent_lower_gt,
        "pbb_lower_gt": pbb_lower_gt,
    }


def derive_verdict(bounds: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    """Frozen decision rules applied to re-derived bounds only.

    DOMINATION: parent UNSAT at cap >= U_b - 1  =>  d_p >= U_b >= d_b.
    REVERSAL:   pbb    UNSAT at cap >= U_p      =>  d_b >  U_p >= d_p.
    A tie can only ever route to DOMINATION (weak, parent has more k).
    """

    U_p, U_b = bounds["U_p"], bounds["U_b"]
    parent_lower_gt = bounds["parent_lower_gt"]
    pbb_lower_gt = bounds["pbb_lower_gt"]
    if U_b is not None and parent_lower_gt is not None and parent_lower_gt >= U_b - 1:
        return "DOMINATION_PROVED", {
            "U_p": U_p,
            "U_b": U_b,
            "parent_lower_gt": parent_lower_gt,
        }
    if U_p is not None and pbb_lower_gt is not None and pbb_lower_gt >= U_p:
        return "CERTIFIED_REVERSAL", {
            "U_p": U_p,
            "U_b": U_b,
            "pbb_lower_gt": pbb_lower_gt,
        }
    return "PENDING", {"U_p": U_p, "U_b": U_b}


def decisive_unsat_calls(
    verdict: str, bounds: dict[str, Any], state: dict[str, Any]
) -> list[dict[str, Any]]:
    """The UNSAT calls the verdict rests on; these must replay."""

    decisive: list[dict[str, Any]] = []
    if verdict == "DOMINATION_PROVED":
        cap = bounds["parent_lower_gt"]
        for target in ("parent_x", "parent_z"):
            entry = next(
                e
                for e in state["calls"]
                if e["target"] == target
                and e["status"] == "UNSAT"
                and int(e["weight_cap"]) == cap
            )
            decisive.append(entry)
    elif verdict == "CERTIFIED_REVERSAL":
        cap = bounds["pbb_lower_gt"]
        entry = next(
            e
            for e in state["calls"]
            if e["target"] == "pbb"
            and e["status"] == "UNSAT"
            and int(e["weight_cap"]) == cap
        )
        decisive.append(entry)
    return decisive


def replay_decisive_unsats(
    instances: dict[str, Any],
    verdict: str,
    bounds: dict[str, Any],
    state: dict[str, Any],
) -> list[dict[str, Any]]:
    """Independently re-prove every decisive UNSAT from rebuilt matrices.

    The replay never trusts the stored call: the CNF is rebuilt from the
    fresh instance, its canonical hash is recorded, and the solver must
    return UNSAT again.  Any SAT or budget outcome fail-stops.
    """

    stamps: list[dict[str, Any]] = []
    for entry in decisive_unsat_calls(verdict, bounds, state):
        record = run_decision(
            instances[entry["target"]],
            int(entry["weight_cap"]),
            conflict_budget=0,
            decomposition=entry.get("decomposition", "monolithic"),
        )
        if record["status"] != "UNSAT":
            raise RuntimeError(
                f"replay refuted stored UNSAT: {entry['target']} at cap "
                f"{entry['weight_cap']} returned {record['status']} - the "
                "stored record is defective or tampered"
            )
        stamps.append(
            {
                "target": entry["target"],
                "weight_cap": int(entry["weight_cap"]),
                "status": "UNSAT",
                "cnf_sha256": record["cnf_sha256"],
                "encoding_version": record["encoding_version"],
                "decomposition": record.get("decomposition", "monolithic"),
                "solver": record["solver"],
                "utc": utc_now(),
            }
        )
    return stamps


def reverify_row(index: int, row: dict[str, Any]) -> dict[str, Any]:
    """Certificate-grade reverification of one decided row.

    Rebuilds identity and instances, re-verifies every stored witness,
    re-derives the verdict from validated evidence only, replays every
    decisive UNSAT, and stamps the row state.  Fail-stops on any mismatch
    with the stored verdict.
    """

    state = verified_row_state(index, row)
    if state is None:
        raise RuntimeError(f"row {index}: no state to reverify")
    instances = build_instances(row)
    bounds = validated_bounds(instances, state)
    verdict, derived_bounds = derive_verdict(bounds)
    if verdict != state.get("verdict"):
        raise RuntimeError(
            f"row {index}: stored verdict {state.get('verdict')!r} does not "
            f"re-derive from validated evidence (got {verdict!r})"
        )
    if verdict in {"DOMINATION_PROVED", "CERTIFIED_REVERSAL"}:
        state["replay"] = {
            "verdict": verdict,
            "bounds": derived_bounds,
            "stamps": replay_decisive_unsats(instances, verdict, bounds, state),
            "utc": utc_now(),
        }
        atomic_write_json(row_state_path(index), state)
    return state


class RowSolver:
    """All SAT interaction for one catalogue row, with persistent state."""

    def __init__(
        self,
        index: int,
        row: dict[str, Any],
        *,
        conflict_budget: int,
        decomposition: str = "monolithic",
    ) -> None:
        self.index = index
        self.row = row
        self.conflict_budget = conflict_budget
        self.decomposition = decomposition
        parent_spec, HX, HZ = E27.parent_matrices(row)
        self.HX, self.HZ = HX, HZ
        _, self.pbb = E27.pbb_code(row)
        self.identity = row_identity(row, index)
        self.state = load_row_state(index, self.identity)
        # Single source of truth for encodings, so solving, replay and
        # validation cannot drift apart.
        self.instances = build_instances(row)

    # -- state-aware primitive ------------------------------------------------
    def cached_status(self, target: str, cap: int) -> dict[str, Any] | None:
        for entry in self.state["calls"]:
            if (
                entry["target"] == target
                and entry["weight_cap"] == cap
                and entry["status"] != "UNDECIDED_BUDGET"
            ):
                return entry
        return None

    def decide(self, target: str, cap: int) -> dict[str, Any]:
        cached = self.cached_status(target, cap)
        if cached is not None:
            return cached
        record = run_decision(
            self.instances[target],
            cap,
            conflict_budget=self.conflict_budget,
            decomposition=self.decomposition,
        )
        record_call(self.state, target=target, record=record)
        self.persist()
        return self.state["calls"][-1]

    def persist(self) -> None:
        atomic_write_json(row_state_path(self.index), self.state)

    # -- composite questions --------------------------------------------------
    def parent_unsat_at(self, cap: int) -> tuple[bool, int | None]:
        """(both_sides_unsat, sat_witness_weight_if_any)."""

        best_sat: int | None = None
        for target in ("parent_x", "parent_z"):
            entry = self.decide(target, cap)
            if entry["status"] == "SAT":
                weight = entry["weight"]
                best_sat = weight if best_sat is None else min(best_sat, weight)
            elif entry["status"] == "UNDECIDED_BUDGET":
                raise BudgetExhausted(f"{target} at cap {cap}")
        if best_sat is not None:
            return False, best_sat
        return True, None

    def pbb_unsat_at(self, cap: int) -> tuple[bool, int | None]:
        entry = self.decide("pbb", cap)
        if entry["status"] == "UNSAT":
            return True, None
        if entry["status"] == "SAT":
            return False, entry["weight"]
        raise BudgetExhausted(f"pbb at cap {cap}")

    def find_witness(self, which: str, hint: int) -> int:
        """Climb from the untrusted hint until SAT; return verified weight."""

        cap = max(1, int(hint))
        while cap <= MAX_CLIMB:
            if which == "parent":
                unsat, sat_weight = self.parent_unsat_at(cap)
            else:
                unsat, sat_weight = self.pbb_unsat_at(cap)
            if not unsat and sat_weight is not None:
                return sat_weight
            cap += 1
        raise RuntimeError(f"no {which} witness found up to weight {MAX_CLIMB}")

    # -- the frozen comparison loop -------------------------------------------
    def compare(self) -> dict[str, Any]:
        row = self.row
        # Untrusted starting caps; a SAT result is a verified witness and an
        # UNSAT during the climb certifies "d > cap" for that code.
        hint = int(row["d"])
        U_p = self.state.get("U_p") or self.find_witness("parent", hint)
        self.state["U_p"] = U_p
        self.persist()
        U_b = self.state.get("U_b") or self.find_witness("pbb", hint)
        self.state["U_b"] = U_b
        self.persist()

        while True:
            if U_p > U_b:
                unsat, sat_weight = self.parent_unsat_at(U_b - 1)
                if unsat:
                    return self.finish(
                        "DOMINATION_PROVED",
                        U_p=U_p,
                        U_b=U_b,
                        parent_lower_gt=U_b - 1,
                    )
                U_p = sat_weight
                self.state["U_p"] = U_p
                self.persist()
            elif U_b > U_p:
                unsat, sat_weight = self.pbb_unsat_at(U_p)
                if unsat:
                    return self.finish(
                        "CERTIFIED_REVERSAL",
                        U_p=U_p,
                        U_b=U_b,
                        pbb_lower_gt=U_p,
                    )
                U_b = sat_weight
                self.state["U_b"] = U_b
                self.persist()
            else:
                unsat, sat_weight = self.parent_unsat_at(U_p - 1)
                if unsat:
                    return self.finish(
                        "DOMINATION_PROVED",
                        U_p=U_p,
                        U_b=U_b,
                        parent_lower_gt=U_p - 1,
                        tie_exact_parent=True,
                    )
                U_p = sat_weight
                self.state["U_p"] = U_p
                self.persist()

    def finish(self, verdict: str, **bounds: Any) -> dict[str, Any]:
        self.state["verdict"] = verdict
        self.state["bounds"] = {
            key: value for key, value in bounds.items() if value is not None
        }
        self.state["completed_utc"] = utc_now()
        self.persist()
        return self.state


def undecided_rows() -> list[tuple[int, dict[str, Any]]]:
    artifact = json.loads(EXP027_ARTIFACT.read_text(encoding="utf-8"))
    verdicts = {
        item["catalogue_index"]: item["verdict"]
        for item in artifact["per_item_results"]
    }
    rows = E27.load_catalogue()
    out: list[tuple[int, dict[str, Any]]] = []
    for index, row in enumerate(rows):
        verdict = verdicts.get(index)
        if verdict in {"UNDECIDED", "UNDECIDED_PARENT_BELOW_PBB_UPPER_BOUND"}:
            out.append((index, row))
    return out

def replay_stamps_valid(
    instances: dict[str, Any],
    verdict: str,
    bounds: dict[str, Any],
    state: dict[str, Any],
) -> bool:
    """Structural validity of replay stamps against rebuilt evidence.

    Every decisive UNSAT derived fresh from the validated call log must have
    a stamp with matching target/cap, status UNSAT, the current encoding
    version, and a CNF hash equal to the hash of the CNF rebuilt now from
    fresh matrices.  This binds stamps to the exact instances; proof-grade
    re-solving lives in the ``verify`` subcommand.
    """

    from qec_research.distance.sat_decide import (
        ENCODING_VERSION,
        decision_cnf_digest,
    )

    replay = state.get("replay") or {}
    if replay.get("verdict") != verdict:
        return False
    stamps = {
        (stamp.get("target"), int(stamp.get("weight_cap", -1))): stamp
        for stamp in replay.get("stamps", [])
    }
    try:
        expected = decisive_unsat_calls(verdict, bounds, state)
    except StopIteration:
        return False
    if not expected:
        return False
    for entry in expected:
        key = (entry["target"], int(entry["weight_cap"]))
        stamp = stamps.get(key)
        if stamp is None or stamp.get("status") != "UNSAT":
            return False
        if stamp.get("encoding_version") != ENCODING_VERSION:
            return False
        # The stamp must hash-bind under the SAME encoding mode the call was
        # proved in; a monolithic hash can never validate a sector proof.
        mode = entry.get("decomposition", "monolithic")
        if stamp.get("decomposition", "monolithic") != mode:
            return False
        expected_digest = decision_cnf_digest(
            instances[entry["target"]],
            int(entry["weight_cap"]),
            decomposition=mode,
        )
        if stamp.get("cnf_sha256") != expected_digest:
            return False
    return True



def assemble() -> dict[str, Any]:
    items = []
    counts: dict[str, int] = {}
    for index, row in undecided_rows():
        state = verified_row_state(index, row)
        if state is None:
            counts["PENDING"] = counts.get("PENDING", 0) + 1
            continue
        stored = state.get("verdict", "PENDING")
        replay_ok = False
        if stored in {"DOMINATION_PROVED", "CERTIFIED_REVERSAL"}:
            # A decided verdict counts as certified only when it re-derives
            # from validated evidence AND its replay stamps hash-bind to the
            # decisive CNFs rebuilt now from fresh matrices.
            instances = build_instances(row)
            bounds = validated_bounds(instances, state)
            derived, _ = derive_verdict(bounds)
            if derived != stored:
                raise RuntimeError(
                    f"row {index}: stored verdict {stored!r} does not "
                    f"re-derive from validated evidence (got {derived!r})"
                )
            replay_ok = replay_stamps_valid(instances, stored, bounds, state)
            if not replay_ok:
                stored = f"{stored}_REPLAY_PENDING"
        verdict = stored
        counts[verdict] = counts.get(verdict, 0) + 1
        items.append(
            {
                key: state.get(key)
                for key in (
                    "catalogue_index",
                    "label",
                    "code_id",
                    "n",
                    "delta",
                    "U_p",
                    "U_b",
                    "bounds",
                    "parent_fingerprint",
                    "pbb_fingerprint",
                    "catalogue_d_untrusted",
                )
            }
            | {
                "verdict": verdict,
                "replay_verified": replay_ok,
                "num_calls": len(state.get("calls", [])),
            }
        )
    payload = {
        "experiment": "EXP-036",
        "schema": "exp036-delta-closure-v2",
        "protocol": __doc__,
        "exp027_artifact_sha256": sha256_bytes(EXP027_ARTIFACT.read_bytes()),
        "solver": SOLVER_NAME,
        "generated_utc": utc_now(),
        "verdict_counts": counts,
        "items": sorted(items, key=lambda item: item["catalogue_index"]),
        "environment": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
        },
    }
    destination = PARTIAL
    if counts and all(
        verdict in {"DOMINATION_PROVED", "CERTIFIED_REVERSAL"}
        for verdict in counts
    ):
        destination = PROCESSED
    atomic_write_json(destination, payload)
    payload["_written_to"] = str(destination.relative_to(ROOT))
    return payload


def run(args: argparse.Namespace) -> int:
    targets = undecided_rows()
    if args.ns:
        wanted = set(args.ns)
        targets = [
            (index, row)
            for index, row in targets
            if 2 * int(row["ell"]) * int(row["m"]) in wanted
        ]
    if args.indexes:
        wanted_indexes = set(args.indexes)
        targets = [
            (index, row) for index, row in targets if index in wanted_indexes
        ]
    print(f"rows selected: {len(targets)}")
    decided = 0
    for index, row in targets:
        label = E27.catalogue_label(row, index)
        existing = verified_row_state(index, row)
        if existing is not None and existing.get("verdict") in {
            "DOMINATION_PROVED",
            "CERTIFIED_REVERSAL",
        }:
            print(f"row {index:04d} {label}: resume hit {existing['verdict']}")
            decided += 1
            continue
        solver = RowSolver(
            index,
            row,
            conflict_budget=args.conflict_budget,
            decomposition=args.decomposition,
        )
        started = time.perf_counter()
        try:
            state = solver.compare()
            wall = time.perf_counter() - started
            print(
                f"row {index:04d} {label}: {state['verdict']} "
                f"U_p={state['U_p']} U_b={state['U_b']} in {wall:.1f}s"
            )
            decided += 1
        except BudgetExhausted as blocked:
            solver.state["verdict"] = "UNDECIDED_BUDGET"
            solver.state["blocked_on"] = str(blocked)
            solver.persist()
            wall = time.perf_counter() - started
            print(f"row {index:04d} {label}: UNDECIDED_BUDGET ({blocked}) {wall:.1f}s")
    summary = assemble()
    print(json.dumps(summary["verdict_counts"], indent=2, sort_keys=True))
    print(f"wrote {summary['_written_to']}")
    return 0 if decided == len(targets) else 2


def verify(args: argparse.Namespace) -> int:
    """Proof-grade reverification: re-solve every decisive UNSAT."""

    verified = 0
    failures = 0
    for index, row in undecided_rows():
        state = verified_row_state(index, row)
        if state is None or state.get("verdict") not in {
            "DOMINATION_PROVED",
            "CERTIFIED_REVERSAL",
        }:
            continue
        if args.ns and 2 * int(row["ell"]) * int(row["m"]) not in set(args.ns):
            continue
        if args.indexes and index not in set(args.indexes):
            continue
        if not args.force:
            instances = build_instances(row)
            bounds = validated_bounds(instances, state)
            if replay_stamps_valid(instances, state["verdict"], bounds, state):
                verified += 1
                continue
        label = E27.catalogue_label(row, index)
        started = time.perf_counter()
        try:
            state = reverify_row(index, row)
            wall = time.perf_counter() - started
            stamps = state["replay"]["stamps"]
            print(
                f"row {index:04d} {label}: {state['verdict']} replay OK "
                f"({len(stamps)} decisive UNSAT re-proved) in {wall:.1f}s"
            )
            verified += 1
        except RuntimeError as defect:
            failures += 1
            print(f"row {index:04d} {label}: REPLAY FAILURE - {defect}")
    summary = assemble()
    print(json.dumps(summary["verdict_counts"], indent=2, sort_keys=True))
    print(f"wrote {summary['_written_to']}")
    return 0 if failures == 0 else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    runner = sub.add_parser("run")
    runner.add_argument("--ns", type=int, nargs="*", help="restrict to code lengths n")
    runner.add_argument("--indexes", type=int, nargs="*", help="restrict to catalogue indexes")
    runner.add_argument(
        "--conflict-budget",
        type=int,
        default=0,
        help="CaDiCaL conflict budget per call; 0 = unbounded",
    )
    runner.add_argument(
        "--decomposition",
        choices=("sectors", "monolithic"),
        default="monolithic",
        help="nontriviality encoding: single-CNF disjunction (default; sector "
        "decomposition measured 8x slower on these instances, FR-022) or the "
        "exact per-functional sector decomposition",
    )
    sub.add_parser("assemble")
    verifier = sub.add_parser("verify")
    verifier.add_argument("--ns", type=int, nargs="*", help="restrict to code lengths n")
    verifier.add_argument("--indexes", type=int, nargs="*", help="restrict to catalogue indexes")
    verifier.add_argument("--force", action="store_true", help="re-replay even stamped rows")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.command == "run":
        return run(args)
    if args.command == "verify":
        return verify(args)
    summary = assemble()
    print(json.dumps(summary["verdict_counts"], indent=2, sort_keys=True))
    print(f"wrote {summary['_written_to']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
