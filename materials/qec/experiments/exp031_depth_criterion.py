"""EXP-031: test a general one-ancilla syndrome-depth criterion.

Lemma C1 in ``proofs/pbb_structure.md`` is a criterion on a *schedule*, not a
closed-form predicate on a check matrix: for every check pair (a,b), the number
of anticommuting shared qubits where a precedes b must be even.  Here the pure
function :func:`parity_depth_criterion` turns that statement into an independent
PySAT decision problem at depth ``w_max``.  It uses only check supports and
Pauli labels.  The prediction under test is

    minimum depth = w_max + int(depth-w_max parity-colouring is infeasible).

This is falsifiable: even if depth w_max is obstructed, depth w_max+1 need not
be feasible a priori.  OR-Tools ``cpsat_schedule`` is the independent ground
truth.  We also measure (but do not adopt) the weaker closed-form candidate
"there exists an anticommuting-overlap pair".
"""

from __future__ import annotations

# Pre-registered protocol constants.  Keep these above all experiment logic.
EXPERIMENT = "EXP-031"
MASTER_SEED = 20260812
RANDOM_COMMUTING_SEED = 20260812
N_RANDOM_COMMUTING = 150
RANDOM_N_RANGE = (8, 14)
RANDOM_CHECK_RANGE = (4, 8)
RANDOM_WEIGHT_RANGE = (3, 6)
CP_SAT_TIME_LIMIT_S = 60.0
CP_SAT_WORKERS = 8
CP_SAT_MAX_EXTRA_DEPTH = 4
PYSAT_SOLVER = "cadical195"
PBB_TARGET = (144, 12, 12)
EXPECTED_BB_INSTANCES = 47
EXPECTED_PBB_INSTANCES = 14

SMOKE_RANDOM_INSTANCES = 1
import contextlib
import argparse
import hashlib
import io
import itertools
import json
import random
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from pysat.solvers import Solver

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

import exp021_verify_equations as exp021  # noqa: E402
from qec_research.circuits.mixed_stabilizer import (  # noqa: E402
    PauliSupport,
    generator_supports,
)
from qec_research.circuits.bicycle_schedule import (  # noqa: E402
    bb_supports_and_orbits,
    pbb_supports_and_orbits,
)
from qec_research.circuits.scheduling import (  # noqa: E402
    cpsat_schedule,
    depth_lower_bound,
    verify_schedule,
)
from qec_research.codes.bicycle import (  # noqa: E402
    BRAVYI_BB,
    BBSpec,
    PBBSpec,
    bb_stabilizer,
    build_pbb,
)
from qec_research.symplectic.core import lambda_swap  # noqa: E402

CATALOGUE = (
    ROOT
    / "third_party"
    / "qcode-discovery"
    / "results"
    / "campaign7_publication_merged.jsonl"
)
OUTPUT = ROOT / "results" / "processed" / "exp031_depth_criterion.json"


@dataclass(frozen=True)
class CriterionResult:
    """Exact result of the independent depth-w_max colouring/parity CSP."""

    w_max: int
    max_qubit_degree: int
    feasible_at_w_max: bool
    obstruction_present: bool
    witness_slot: dict[tuple[int, int], int] | None
    num_variables: int
    num_clauses: int
    anticommuting_pair_count: int
    anticommuting_incidence_count: int


class _VarPool:
    def __init__(self) -> None:
        self.top = 0

    def new(self) -> int:
        self.top += 1
        return self.top


def _normalise_supports(supports: Sequence[PauliSupport]) -> list[PauliSupport]:
    """Validate and copy a check-support description into canonical order."""
    ordered = sorted(supports, key=lambda support: support.index)
    if [support.index for support in ordered] != list(range(len(ordered))):
        raise ValueError("check indices must be contiguous 0..m-1")
    out: list[PauliSupport] = []
    for support in ordered:
        paulis = dict(sorted(support.paulis.items()))
        if not paulis:
            raise ValueError("zero-weight checks are outside this experiment")
        if any(not isinstance(q, int) or q < 0 for q in paulis):
            raise ValueError("qubit indices must be nonnegative integers")
        if any(pauli not in {"X", "Y", "Z"} for pauli in paulis.values()):
            raise ValueError("Pauli labels must be X, Y, or Z")
        out.append(PauliSupport(support.index, paulis))
    return out


def _anticommuting_overlaps(
    supports: Sequence[PauliSupport],
) -> dict[tuple[int, int], tuple[int, ...]]:
    """Compute Lemma-C1 sets J(a,b), independently of scheduling.py."""
    overlaps: dict[tuple[int, int], tuple[int, ...]] = {}
    for left, right in itertools.combinations(supports, 2):
        shared = sorted(set(left.paulis) & set(right.paulis))
        js = tuple(q for q in shared if left.paulis[q] != right.paulis[q])
        if js:
            overlaps[(left.index, right.index)] = js
    return overlaps


def weak_anticommuting_overlap_predicate(
    supports: Sequence[PauliSupport],
) -> bool:
    """Weak closed-form candidate: at least one nonempty J(a,b)."""
    return bool(_anticommuting_overlaps(_normalise_supports(supports)))


def _add_at_most_one(clauses: list[list[int]], variables: Iterable[int]) -> None:
    variables = list(variables)
    for left, right in itertools.combinations(variables, 2):
        clauses.append([-left, -right])


def _depth_colouring_cnf(
    supports: Sequence[PauliSupport], depth: int
) -> tuple[list[list[int]], dict[tuple[int, int, int], int], int, dict]:
    """Encode incidence colouring plus the exact Lemma-C1 parity equations."""
    supports = _normalise_supports(supports)
    if depth < 1:
        raise ValueError("depth must be positive")
    overlaps = _anticommuting_overlaps(supports)
    if any(len(js) % 2 for js in overlaps.values()):
        raise ValueError("input checks do not commute: an anticommuting overlap is odd")

    pool = _VarPool()
    clauses: list[list[int]] = []
    edges = [(support.index, q) for support in supports for q in support.paulis]
    colour = {
        (check, qubit, slot): pool.new()
        for check, qubit in edges
        for slot in range(depth)
    }

    # Each incidence gets exactly one colour (time slot).
    for check, qubit in edges:
        variables = [colour[(check, qubit, slot)] for slot in range(depth)]
        clauses.append(variables)
        _add_at_most_one(clauses, variables)

    # Proper bipartite edge colouring: no check ancilla or data qubit is reused.
    by_check: dict[int, list[int]] = {}
    by_qubit: dict[int, list[int]] = {}
    for check, qubit in edges:
        by_check.setdefault(check, []).append(qubit)
        by_qubit.setdefault(qubit, []).append(check)
    for slot in range(depth):
        for check, qubits in by_check.items():
            _add_at_most_one(
                clauses, (colour[(check, qubit, slot)] for qubit in qubits)
            )
        for qubit, checks in by_qubit.items():
            _add_at_most_one(
                clauses, (colour[(check, qubit, slot)] for check in checks)
            )

    # o_(a,b,j) means check a acts before check b at shared data qubit j.
    # One-hot colours make the implications below an exact definition of o.
    order: dict[tuple[int, int, int], int] = {}
    for (left, right), shared in overlaps.items():
        order_variables: list[int] = []
        for qubit in shared:
            before = pool.new()
            order[(left, right, qubit)] = before
            order_variables.append(before)
            for left_slot in range(depth):
                for right_slot in range(depth):
                    if left_slot == right_slot:
                        continue  # already forbidden by the data-qubit constraint
                    consequence = before if left_slot < right_slot else -before
                    clauses.append(
                        [
                            -colour[(left, qubit, left_slot)],
                            -colour[(right, qubit, right_slot)],
                            consequence,
                        ]
                    )

        # Lemma C1: exclude every odd assignment of the order variables.
        for bits in itertools.product((0, 1), repeat=len(order_variables)):
            if sum(bits) % 2:
                clauses.append(
                    [
                        -variable if bit else variable
                        for variable, bit in zip(order_variables, bits)
                    ]
                )

    degrees: dict[int, int] = {}
    for _, qubit in edges:
        degrees[qubit] = degrees.get(qubit, 0) + 1
    metadata = {
        "edges": edges,
        "overlaps": overlaps,
        "max_qubit_degree": max(degrees.values(), default=0),
    }
    return clauses, colour, pool.top, metadata


def parity_depth_criterion(
    supports: Sequence[PauliSupport],
    *,
    symmetry_orbits: dict[tuple[int, int], int] | None = None,
) -> CriterionResult:
    """Pure support/Pauli predicate for a depth-w_max parity obstruction.

    The implementation has no dependence on OR-Tools, lattice labels, or code
    parameters.  It exactly decides whether a proper w_max-colouring satisfying
    every Lemma-C1 parity equation exists, using an independent SAT backend.
    ``symmetry_orbits`` optionally restricts the tested schedules to a declared
    translation-invariant subclass; this is recorded and never silently mixed
    with unrestricted random-instance results.
    """
    supports = _normalise_supports(supports)
    w_max = max(support.weight for support in supports)
    clauses, colour, num_variables, metadata = _depth_colouring_cnf(supports, w_max)
    if symmetry_orbits is not None:
        grouped: dict[int, list[tuple[int, int]]] = {}
        for edge in metadata["edges"]:
            grouped.setdefault(symmetry_orbits[edge], []).append(edge)
        for edges in grouped.values():
            representative = edges[0]
            for edge in edges[1:]:
                for slot in range(w_max):
                    left = colour[(*representative, slot)]
                    right = colour[(*edge, slot)]
                    clauses.extend(([-left, right], [left, -right]))
    with Solver(name=PYSAT_SOLVER, bootstrap_with=clauses) as solver:
        feasible = bool(solver.solve())
        model = solver.get_model() if feasible else None

    witness: dict[tuple[int, int], int] | None = None
    if model is not None:
        positive = {literal for literal in model if literal > 0}
        witness = {}
        for check, qubit in metadata["edges"]:
            chosen = [
                slot
                for slot in range(w_max)
                if colour[(check, qubit, slot)] in positive
            ]
            if len(chosen) != 1:
                raise AssertionError("PySAT model does not assign one slot per incidence")
            witness[(check, qubit)] = chosen[0]
        verification = verify_schedule(supports, witness)
        if not verification["valid"]:
            raise AssertionError(f"invalid PySAT witness: {verification}")

    overlaps = metadata["overlaps"]
    return CriterionResult(
        w_max=w_max,
        max_qubit_degree=metadata["max_qubit_degree"],
        feasible_at_w_max=feasible,
        obstruction_present=not feasible,
        witness_slot=witness,
        num_variables=num_variables,
        num_clauses=len(clauses),
        anticommuting_pair_count=len(overlaps),
        anticommuting_incidence_count=sum(len(js) for js in overlaps.values()),
    )


def _commutes_via_lambda(H: np.ndarray) -> bool:
    product = H.astype(np.int64) @ lambda_swap(H).T.astype(np.int64)
    return not bool((product & 1).any())


def _pauli_rows(H: np.ndarray) -> list[str]:
    n = H.shape[1] // 2
    labels = np.array(["I", "X", "Z", "Y"])
    out = []
    for row in H:
        encoded = row[:n].astype(int) + 2 * row[n:].astype(int)
        out.append("".join(labels[encoded]))
    return out


def _row_weight(row: np.ndarray) -> int:
    n = row.size // 2
    return int((row[:n] | row[n:]).sum())


def _make_css_row(n: int, support: Sequence[int], pauli: str) -> np.ndarray:
    row = np.zeros(2 * n, dtype=np.uint8)
    offset = 0 if pauli == "X" else n
    row[offset + np.asarray(support, dtype=int)] = 1
    return row


def _random_commuting_instance(
    rng: random.Random, instance_index: int
) -> tuple[np.ndarray, dict]:
    """Build CSS bases, then mixed rows as products of commuting bases."""
    for construction_attempt in range(1, 20001):
        n = rng.randint(*RANDOM_N_RANGE)
        num_checks = rng.randint(*RANDOM_CHECK_RANGE)
        num_mixed = 1 if num_checks <= 5 else 2
        num_base = num_checks - num_mixed
        num_x = max(1, num_base // 2)
        num_z = num_base - num_x

        x_rows: list[np.ndarray] = []
        seen: set[bytes] = set()
        for _ in range(num_x):
            for _ in range(200):
                weight = rng.randint(*RANDOM_WEIGHT_RANGE)
                support = sorted(rng.sample(range(n), weight))
                row = _make_css_row(n, support, "X")
                if row.tobytes() not in seen:
                    seen.add(row.tobytes())
                    x_rows.append(row)
                    break
            else:
                break
        if len(x_rows) != num_x:
            continue

        z_rows: list[np.ndarray] = []
        for _ in range(num_z):
            for _ in range(1000):
                weight = rng.randint(*RANDOM_WEIGHT_RANGE)
                support = sorted(rng.sample(range(n), weight))
                row = _make_css_row(n, support, "Z")
                if row.tobytes() in seen:
                    continue
                if all(int((row[n:] @ xrow[:n]) & 1) == 0 for xrow in x_rows):
                    seen.add(row.tobytes())
                    z_rows.append(row)
                    break
            else:
                break
        if len(z_rows) != num_z:
            continue

        bases = x_rows + z_rows
        mixed_rows: list[np.ndarray] = []
        factors: list[list[int]] = []
        # Products of an X base and a Z base are mixed and commute with all bases.
        candidates = [(x, num_x + z) for x in range(num_x) for z in range(num_z)]
        rng.shuffle(candidates)
        for x_index, z_index in candidates:
            row = bases[x_index] ^ bases[z_index]
            weight = _row_weight(row)
            if RANDOM_WEIGHT_RANGE[0] <= weight <= RANDOM_WEIGHT_RANGE[1]:
                key = row.tobytes()
                if key not in seen:
                    seen.add(key)
                    mixed_rows.append(row)
                    factors.append([x_index, z_index])
                    if len(mixed_rows) == num_mixed:
                        break
        if len(mixed_rows) != num_mixed:
            continue

        H = np.vstack(bases + mixed_rows).astype(np.uint8)
        weights = [_row_weight(row) for row in H]
        if not all(RANDOM_WEIGHT_RANGE[0] <= weight <= RANDOM_WEIGHT_RANGE[1] for weight in weights):
            continue
        if not _commutes_via_lambda(H):
            raise AssertionError("by-construction random checks failed lambda_swap verification")
        supports = generator_supports(H)
        w_max = max(weights)
        if depth_lower_bound(supports, n) != w_max:
            # The conjecture is stated relative to w_max, so avoid a separate,
            # trivial data-degree obstruction in the random validation tier.
            continue
        metadata = {
            "generator": "random CSS bases plus mixed products",
            "instance_index": instance_index,
            "construction_attempt": construction_attempt,
            "n": n,
            "num_checks": num_checks,
            "num_x_bases": num_x,
            "num_z_bases": num_z,
            "num_mixed_products": num_mixed,
            "mixed_product_factors": factors,
            "pauli_rows": _pauli_rows(H),
            "lambda_swap_commutation_verified": True,
        }
        return H, metadata
    raise RuntimeError(f"could not construct random commuting instance {instance_index}")


def _exp021_bb_instances() -> tuple[list[dict], str]:
    """Reproduce exp021 Part II's exact 7+40 BB population and RNG state."""
    rng = random.Random(MASTER_SEED)
    old_len = len(exp021.RESULTS)
    with contextlib.redirect_stdout(io.StringIO()):
        # These are exactly the calls before Part II in exp021's __main__ block.
        exp021.check_E1_E3(rng)
        exp021.check_E5(rng)
    del exp021.RESULTS[old_len:]

    items: list[dict] = []
    for label, spec in BRAVYI_BB.items():
        items.append({"label": f"Bravyi {label}", "source": "published", "spec": spec})
    for index in range(40):
        items.append(
            {
                "label": f"exp021 seeded-random BB {index:02d}",
                "source": "exp021_part_ii_seeded_random",
                "spec": exp021.rand_bb(rng),
            }
        )
    if len(items) != EXPECTED_BB_INSTANCES:
        raise AssertionError(f"expected {EXPECTED_BB_INSTANCES} BB instances")
    serialised = [
        {
            "label": item["label"],
            "ell": item["spec"].ell,
            "m": item["spec"].m,
            "A": item["spec"].A,
            "B": item["spec"].B,
        }
        for item in items
    ]
    fingerprint = hashlib.sha256(
        json.dumps(serialised, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return items, fingerprint


def _pbb_instances() -> list[dict]:
    items = []
    for row_number, line in enumerate(CATALOGUE.read_text().splitlines(), start=1):
        row = json.loads(line)
        if (row.get("n"), row.get("k"), row.get("d")) != PBB_TARGET:
            continue
        label = (
            row.get("code_id")
            or row.get("code")
            or row.get("label")
            or f"catalogue-line-{row_number}"
        )
        spec = PBBSpec(
            ell=row["ell"],
            m=row["m"],
            A=[tuple(term) for term in row["A_terms"]],
            B=[tuple(term) for term in row["B_terms"]],
            C=[tuple(term) for term in (row["C_terms"] or [])],
            D=[tuple(term) for term in (row["D_terms"] or [])],
            name=str(label),
        )
        items.append(
            {
                "label": str(label),
                "row_number": row_number,
                "row": row,
                "spec": spec,
            }
        )
    if len(items) != EXPECTED_PBB_INSTANCES:
        raise AssertionError(
            f"expected {EXPECTED_PBB_INSTANCES} headline PBB rows, got {len(items)}"
        )
    return items


def _serialise_slot(slot: dict[tuple[int, int], int] | None) -> list[list[int]] | None:
    if slot is None:
        return None
    return [[check, qubit, depth] for (check, qubit), depth in sorted(slot.items())]


def _certify_minimum_depth(
    supports: list[PauliSupport],
    n: int,
    w_max: int,
    solver_seed: int,
    *,
    symmetry_orbits: dict[tuple[int, int], int] | None = None,
    criterion_witness: dict[tuple[int, int], int] | None = None,
) -> dict:
    """CP-SAT ground truth within the explicitly recorded schedule class."""
    trace = []
    first_feasible = None
    # OR-Tools rejects repeated hints to one orbit variable.  PySAT witnesses
    # contain one assignment per edge, so only hint unrestricted models.
    seed_slot = criterion_witness if symmetry_orbits is None else None
    for depth in range(w_max, w_max + CP_SAT_MAX_EXTRA_DEPTH + 1):
        result = cpsat_schedule(
            supports,
            n,
            T=depth,
            time_limit_s=CP_SAT_TIME_LIMIT_S,
            workers=CP_SAT_WORKERS,
            random_seed=solver_seed,
            deterministic=False,
            seed_slot=seed_slot,
            symmetry_orbits=symmetry_orbits,
        )
        valid = bool(result.verification and result.verification.get("valid"))
        trace.append(
            {
                "depth": depth,
                "status": result.status,
                "wall_s": result.wall_time_s,
                "valid": valid,
            }
        )
        if result.status == "OPTIMAL":
            if result.slot is None or not valid:
                raise AssertionError("CP-SAT returned OPTIMAL without a valid schedule")
            first_feasible = {
                "depth": depth,
                "slot": _serialise_slot(result.slot),
                "verification": result.verification,
            }
            break
        if result.status != "INFEASIBLE":
            return {
                "decided": False,
                "reason": f"solver status {result.status} at depth {depth}",
                "trace": trace,
                "certified_depth": None,
            }
    if first_feasible is None:
        return {
            "decided": False,
            "reason": (
                f"no feasible schedule through w_max+{CP_SAT_MAX_EXTRA_DEPTH}; "
                "minimum remains undecided"
            ),
            "trace": trace,
            "certified_depth": None,
        }
    if any(entry["status"] != "INFEASIBLE" for entry in trace[:-1]):
        raise AssertionError("minimum-depth trace contains an uncertified lower depth")
    return {
        "decided": True,
        "reason": "exact(OPTIMAL) after certified INFEASIBLE lower depths",
        "trace": trace,
        "certified_depth": first_feasible["depth"],
        "slot": first_feasible["slot"],
        "verification": first_feasible["verification"],
    }


def certified_infeasible_at(trace: list[dict], depth: int) -> bool:
    """True iff the solver CERTIFIED INFEASIBLE at exactly this depth.

    Timeouts/UNKNOWN never count: only a real infeasibility certificate can
    refute a prediction (FR-014).
    """
    return any(e["depth"] == depth and e["status"] == "INFEASIBLE" for e in trace)


def classify_prediction(*, decided: bool, certified_depth, predicted_depth: int,
                        w_max: int, trace: list[dict]) -> dict:
    """Three-valued prediction/two-value-law classification (FR-014).

    A certified INFEASIBLE at the predicted depth refutes the prediction even
    when the true minimum remains unknown; certified INFEASIBLE at both w_max
    and w_max+1 refutes the two-value law in the tested schedule class.  The
    original classifier only compared decided minima, silently filing certified
    refutations under 'undecided'.
    """
    pred_refuted = certified_infeasible_at(trace, predicted_depth)
    two_value_refuted = (certified_infeasible_at(trace, w_max)
                         and certified_infeasible_at(trace, w_max + 1))
    if decided:
        matches = bool(certified_depth == predicted_depth)
        two_value = certified_depth in {w_max, w_max + 1}
    else:
        matches = False if pred_refuted else None
        two_value = False if two_value_refuted else None
    return {
        "matches_predicted_depth": matches,
        "two_value_law_holds": two_value,
        "prediction_refuted_in_class": pred_refuted,
        "two_value_refuted_in_class": two_value_refuted,
    }


def _evaluate(
    *,
    instance_id: str,
    family: str,
    H: np.ndarray,
    full_spec: dict,
    solver_seed: int,
    symmetry_orbits: dict[tuple[int, int], int] | None = None,
) -> dict:
    if not _commutes_via_lambda(H):
        raise AssertionError(f"{instance_id} does not commute under lambda_swap")
    n = H.shape[1] // 2
    supports = generator_supports(H)
    ordinary_lb = depth_lower_bound(supports, n)

    criterion_start = time.perf_counter()
    criterion = parity_depth_criterion(
        supports, symmetry_orbits=symmetry_orbits
    )
    criterion_wall = time.perf_counter() - criterion_start
    predicted_depth = criterion.w_max + int(criterion.obstruction_present)
    weak = weak_anticommuting_overlap_predicate(supports)

    certification = _certify_minimum_depth(
        supports,
        n,
        criterion.w_max,
        solver_seed,
        symmetry_orbits=symmetry_orbits,
        criterion_witness=criterion.witness_slot,
    )
    certified_depth = certification["certified_depth"]
    decided = certification["decided"]
    classification = classify_prediction(
        decided=decided,
        certified_depth=certified_depth,
        predicted_depth=predicted_depth,
        w_max=criterion.w_max,
        trace=certification["trace"],
    )
    reduction_match = None
    if certification["trace"]:
        at_w = certification["trace"][0]["status"]
        if at_w in {"OPTIMAL", "INFEASIBLE"}:
            reduction_match = (at_w == "OPTIMAL") == criterion.feasible_at_w_max

    schedule_class = (
        "translation-invariant" if symmetry_orbits is not None else "unrestricted"
    )
    return {
        "instance_id": instance_id,
        "family": family,
        "n": n,
        "num_checks": len(supports),
        "full_spec": full_spec,
        "lambda_swap_commutation_verified": True,
        "schedule_class": schedule_class,
        "w_max": criterion.w_max,
        "max_qubit_degree": criterion.max_qubit_degree,
        "ordinary_incidence_lower_bound": ordinary_lb,
        "criterion": {
            "definition": "depth-w_max incidence-colouring plus all Lemma-C1 even-order equations",
            "backend": f"PySAT {PYSAT_SOLVER}",
            "schedule_class": schedule_class,
            "feasible_at_w_max": criterion.feasible_at_w_max,
            "parity_obstruction_present": criterion.obstruction_present,
            "predicted_depth": predicted_depth,
            "num_variables": criterion.num_variables,
            "num_clauses": criterion.num_clauses,
            "anticommuting_pair_count": criterion.anticommuting_pair_count,
            "anticommuting_incidence_count": criterion.anticommuting_incidence_count,
            "witness_verified": criterion.witness_slot is not None,
            "wall_s": criterion_wall,
        },
        "weak_closed_form_candidate": {
            "definition": "exists a check pair with nonempty anticommuting overlap J(a,b)",
            "predicts_obstruction": weak,
            "agrees_with_exact_criterion": weak == criterion.obstruction_present,
        },
        "ground_truth": certification,
        "solver_seed": solver_seed,
        "decided": decided,
        "matches_predicted_depth": classification["matches_predicted_depth"],
        "prediction_refuted_in_class": classification["prediction_refuted_in_class"],
        "two_value_refuted_in_class": classification["two_value_refuted_in_class"],
        "reduction_law_matches_at_w_max": reduction_match,
        "two_value_law_holds": classification["two_value_law_holds"],
    }


def _mismatch_record(record: dict) -> dict:
    """Retain the full reproducible instance specification in every mismatch."""
    return {
        "instance_id": record["instance_id"],
        "w_max": record["w_max"],
        "schedule_class": record["schedule_class"],
        "predicted_depth": record["criterion"]["predicted_depth"],
        "certified_depth": record["ground_truth"]["certified_depth"],
        "decided": record["decided"],
        "prediction_refuted_in_class": record["prediction_refuted_in_class"],
        "two_value_refuted_in_class": record["two_value_refuted_in_class"],
        "trace_statuses": [(e["depth"], e["status"])
                           for e in record["ground_truth"]["trace"]],
        "family": record["family"],
        "full_spec": record["full_spec"],
        "w_max": record["w_max"],
        "max_qubit_degree": record["max_qubit_degree"],
        "criterion": record["criterion"],
        "ground_truth": record["ground_truth"],
        "reduction_law_matches_at_w_max": record["reduction_law_matches_at_w_max"],
        "two_value_law_holds": record["two_value_law_holds"],
    }


def decide_verdict(n_mismatches: int, n_decided: int, n_undecided: int,
                   decided_floor: int = 150) -> tuple[str, str]:
    """Verdict polarity contract (FR-013): a counterexample FALSIFIES the conjecture.

    NEGATIVE      -> at least one decided instance violates the criterion.
    POSITIVE      -> zero mismatches on a sufficiently large decided set.
    INCONCLUSIVE  -> too few decided instances to support either claim.

    A 'POSITIVE' verdict must never be emitted on the mismatch path: downstream
    integration reads this single field, and both outcomes mapping to POSITIVE
    made the artifact machine-unreadable (advisory 2026-08-12).
    """
    if n_mismatches:
        return "NEGATIVE", (
            f"criterion FALSIFIED: {n_mismatches} instance(s) violate the proposed "
            "depth criterion (certified in-class refutations count even when "
            "the true minimum remains undecided)"
        )
    if n_decided >= decided_floor:
        return "POSITIVE", (
            f"criterion holds on all {n_decided} decided instances; theorem candidate"
        )
    return "INCONCLUSIVE", (
        f"only {n_decided} instances were certified; "
        f"{n_undecided} remained undecided"
    )


def main() -> None:
    start = time.time()
    records: list[dict] = []
    bb_items, bb_fingerprint = _exp021_bb_instances()

    serial = 0
    for item in bb_items:
        spec: BBSpec = item["spec"]
        code, _, orbits = bb_supports_and_orbits(spec)
        full_spec = {
            "source": item["source"],
            "construction": "CSS bivariate-bicycle",
            "ell": spec.ell,
            "m": spec.m,
            "A_terms": [list(term) for term in spec.A],
            "B_terms": [list(term) for term in spec.B],
        }
        seed = (MASTER_SEED + 1009 * serial) % 2_147_483_647
        record = _evaluate(
            instance_id=f"bb-{serial:03d}",
            family="exp021_part_ii_bb",
            H=code.H,
            full_spec={"label": item["label"], **full_spec},
            solver_seed=seed,
            symmetry_orbits=orbits,
        )
        records.append(record)
        serial += 1
        print(
            f"{record['instance_id']} {item['label']}: "
            f"pred={record['criterion']['predicted_depth']} "
            f"cert={record['ground_truth']['certified_depth']} "
            f"status={record['ground_truth']['trace'][-1]['status']}",
            flush=True,
        )

    for pbb_index, item in enumerate(_pbb_instances()):
        spec: PBBSpec = item["spec"]
        code, _, orbits = pbb_supports_and_orbits(spec)
        row = item["row"]
        full_spec = {
            "source": "campaign7_publication_merged.jsonl",
            "catalogue_row_number": item["row_number"],
            "label": item["label"],
            "ell": spec.ell,
            "m": spec.m,
            "n": row["n"],
            "k": row["k"],
            "d": row["d"],
            "A_terms": [list(term) for term in spec.A],
            "B_terms": [list(term) for term in spec.B],
            "C_terms": [list(term) for term in spec.C],
            "D_terms": [list(term) for term in spec.D],
        }
        seed = (MASTER_SEED + 1009 * serial) % 2_147_483_647
        record = _evaluate(
            instance_id=f"pbb144-{pbb_index:02d}",
            family="pbb_144_12_12",
            H=code.H,
            full_spec=full_spec,
            solver_seed=seed,
            symmetry_orbits=orbits,
        )
        records.append(record)
        serial += 1
        print(
            f"{record['instance_id']} {item['label']}: "
            f"pred={record['criterion']['predicted_depth']} "
            f"cert={record['ground_truth']['certified_depth']} "
            f"status={record['ground_truth']['trace'][-1]['status']}",
            flush=True,
        )

    random_rng = random.Random(RANDOM_COMMUTING_SEED)
    for random_index in range(N_RANDOM_COMMUTING):
        H, full_spec = _random_commuting_instance(random_rng, random_index)
        seed = (MASTER_SEED + 1009 * serial) % 2_147_483_647
        record = _evaluate(
            instance_id=f"random-{random_index:03d}",
            family="random_commuting",
            H=H,
            full_spec=full_spec,
            solver_seed=seed,
        )
        records.append(record)
        serial += 1
        print(
            f"{record['instance_id']}: pred={record['criterion']['predicted_depth']} "
            f"cert={record['ground_truth']['certified_depth']} "
            f"status={record['ground_truth']['trace'][-1]['status']}",
            flush=True,
        )

    decided = [record for record in records if record["decided"]]
    undecided = [record for record in records if not record["decided"]]
    decided_matches = sum(1 for r in decided if r["matches_predicted_depth"] is True)
    refuted_undecided = sum(1 for r in undecided
                            if r["prediction_refuted_in_class"] or r["two_value_refuted_in_class"])
    mismatch_records = [
        _mismatch_record(record)
        for record in records
        if record["matches_predicted_depth"] is False
    ]
    reduction_mismatches = [
        _mismatch_record(record)
        for record in decided
        if record["reduction_law_matches_at_w_max"] is False
    ]
    two_value_mismatches = [
        _mismatch_record(record)
        for record in records
        if record["two_value_law_holds"] is False
    ]

    weak_tp = sum(
        record["weak_closed_form_candidate"]["predicts_obstruction"]
        and record["criterion"]["parity_obstruction_present"]
        for record in records
    )
    weak_fp = sum(
        record["weak_closed_form_candidate"]["predicts_obstruction"]
        and not record["criterion"]["parity_obstruction_present"]
        for record in records
    )
    weak_tn = sum(
        not record["weak_closed_form_candidate"]["predicts_obstruction"]
        and not record["criterion"]["parity_obstruction_present"]
        for record in records
    )
    weak_fn = sum(
        not record["weak_closed_form_candidate"]["predicts_obstruction"]
        and record["criterion"]["parity_obstruction_present"]
        for record in records
    )
    weak_mismatches = [
        {
            "instance_id": record["instance_id"],
            "family": record["family"],
            "full_spec": record["full_spec"],
            "weak_predicts_obstruction": record["weak_closed_form_candidate"][
                "predicts_obstruction"
            ],
            "exact_csp_obstruction": record["criterion"][
                "parity_obstruction_present"
            ],
        }
        for record in records
        if not record["weak_closed_form_candidate"]["agrees_with_exact_criterion"]
    ]

    verdict, verdict_reason = decide_verdict(
        len(mismatch_records), len(decided), len(undecided)
    )

    payload = {
        "experiment": EXPERIMENT,
        "protocol_constants": {
            "master_seed": MASTER_SEED,
            "random_commuting_seed": RANDOM_COMMUTING_SEED,
            "n_random_commuting": N_RANDOM_COMMUTING,
            "random_n_range_inclusive": list(RANDOM_N_RANGE),
            "random_check_range_inclusive": list(RANDOM_CHECK_RANGE),
            "random_weight_range_inclusive": list(RANDOM_WEIGHT_RANGE),
            "cp_sat_time_limit_s_per_depth": CP_SAT_TIME_LIMIT_S,
            "cp_sat_workers": CP_SAT_WORKERS,
            "cp_sat_max_extra_depth": CP_SAT_MAX_EXTRA_DEPTH,
            "pysat_solver": PYSAT_SOLVER,
            "pbb_target_n_k_d": list(PBB_TARGET),
            "criterion_reading": (
                "C3/Lemma C1 defines schedule validity: every check-pair count of "
                "anticommuting shared qubits with a-before-b is even. The pure "
                "criterion is exact infeasibility of that incidence-colouring CSP at w_max."
            ),
            "ground_truth": (
                "cpsat_schedule in the explicitly recorded schedule class; exact only "
                "for OPTIMAL schedules after INFEASIBLE certificates at every smaller "
                "tested depth. Structured BB/PBB tiers use translation-invariant "
                "colourings; random sets use unrestricted colourings."
            ),
        },
        "seeds": {
            "master": MASTER_SEED,
            "random_commuting_generator": RANDOM_COMMUTING_SEED,
            "per_instance_cp_sat": {
                record["instance_id"]: record["solver_seed"] for record in records
            },
        },
        "population": {
            "attempted": len(records),
            "exp021_part_ii_bb": len(bb_items),
            "exp021_part_ii_bb_composition": (
                "7 published BRAVYI_BB plus 40 bit-identical seeded-random BB specs"
            ),
            "exp021_part_ii_bb_sha256": bb_fingerprint,
            "pbb_144_12_12": EXPECTED_PBB_INSTANCES,
            "random_commuting": N_RANDOM_COMMUTING,
        },
        "n_tested": len(decided),
        "n_undecided": len(undecided),
        "n_refuted_undecided": refuted_undecided,
        "matches": decided_matches,
        "mismatches": mismatch_records,
        "reduction_law_mismatches": reduction_mismatches,
        "two_value_law_mismatches": two_value_mismatches,
        "undecided": [
            {
                "instance_id": record["instance_id"],
                "family": record["family"],
                "full_spec": record["full_spec"],
                "ground_truth": record["ground_truth"],
            }
            for record in undecided
        ],
        "weak_closed_form_candidate": {
            "definition": "exists (a,b) with nonempty anticommuting overlap J(a,b)",
            "truth_definition": "exact PySAT depth-w_max obstruction",
            "true_positive": weak_tp,
            "false_positive": weak_fp,
            "true_negative": weak_tn,
            "false_negative": weak_fn,
            "false_positive_rate": weak_fp / (weak_fp + weak_tn)
            if weak_fp + weak_tn
            else None,
            "false_negative_rate": weak_fn / (weak_fn + weak_tp)
            if weak_fn + weak_tp
            else None,
            "mismatches": weak_mismatches,
        },
        "per_item_results": records,
        "verdict": verdict,
        "verdict_reason": verdict_reason,
        "wall_s": round(time.time() - start, 3),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"wrote {OUTPUT.relative_to(ROOT)}: decided={len(decided)} "
        f"matches={payload['matches']} mismatches={len(mismatch_records)} "
        f"undecided={len(undecided)} verdict={verdict}",
        flush=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--smoke",
        action="store_true",
        help=(
            "run a one-instance import/criterion smoke check; never writes the "
            "canonical artifact"
        ),
    )
    args = parser.parse_args()
    if args.smoke:
        smoke_rng = random.Random(RANDOM_COMMUTING_SEED)
        smoke_H, _ = _random_commuting_instance(smoke_rng, 0)
        smoke_result = parity_depth_criterion(generator_supports(smoke_H))
        partial_dir = ROOT / "results" / "partial_runs"
        partial_dir.mkdir(parents=True, exist_ok=True)
        partial_path = partial_dir / "exp031_SMOKE_seed20260812_rows1of211.json"
        partial_path.write_text(
            json.dumps(
                {
                    "_NOT_CANONICAL": "clean partial smoke run; canonical artifact untouched",
                    "experiment": EXPERIMENT,
                    "seed": RANDOM_COMMUTING_SEED,
                    "attempted": SMOKE_RANDOM_INSTANCES,
                    "full_scope": 211,
                    "criterion_obstruction": smoke_result.obstruction_present,
                },
                indent=2,
            )
            + "\n"
        )
        print(f"partial smoke -> {partial_path.relative_to(ROOT)}")
    else:
        main()
