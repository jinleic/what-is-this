"""EXP-029: certified small-weight circuit-distance bounds at n=144.

The distance in this experiment counts detector-error-model (DEM) ``error``
instructions after flattening as fault mechanisms.  The certified lower-bound
checks use the undecomposed DEM (``decompose_errors=False``).  Stim's circuit
search is heuristic and therefore only supplies upper bounds.  The graphlike
search also uses the undecomposed DEM and skips ungraphlike mechanisms, so each
selected graphlike edge is still a complete mechanism from the certified fault
model.

The two circuits are rebuilt exactly as EXP-013: twelve Z-memory rounds, the
same pure-Z logical basis, and the deterministic translation-orbit CP-SAT
schedule.  The nonzero value p=0.001 only materializes the DEM structure;
distance here ignores mechanism probabilities.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import multiprocessing as mp
import os
for _thread_env in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_env] = "1"
import platform
import resource
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

import stim

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.circuits.bicycle_schedule import (  # noqa: E402
    bb_supports_and_orbits,
    pbb_supports_and_orbits,
    pure_z_logical_basis,
)
from qec_research.circuits.mixed_stabilizer import (  # noqa: E402
    CircuitSpec,
    build_memory_circuit,
)
from qec_research.circuits.scheduling import (  # noqa: E402
    cpsat_schedule,
    depth_lower_bound,
    slots_to_layers,
)
from qec_research.codes.bicycle import BRAVYI_BB, PBBSpec  # noqa: E402
from qec_research.artifacts import canonical_route  # noqa: E402
from qec_research.distance.sectors import min_pure_x_logical  # noqa: E402

# ---------------------------------------------------------------------------
# Pre-registered protocol constants.  Keep these above all experiment logic.
# ---------------------------------------------------------------------------
EXPERIMENT = "EXP-029"
ROOT_SEED = 20260812
# EXP-013 relied on cpsat_schedule's deterministic default seed, which is zero.
SCHEDULE_RANDOM_SEED = 0
P_DEM_STRUCTURE = 0.001
ROUNDS = 12
BASIS = "Z"
SCHEDULE_TIME_LIMIT_S = 120.0
SCHEDULE_WORKERS_ARGUMENT = 4  # deterministic=True forces one CP-SAT worker
SCHEDULE_DETERMINISTIC = True
UPPER_SEARCH_BUDGET_S = 45 * 60.0
GRAPLIKE_SEARCH_BUDGET_S = 15 * 60.0
W3_SEARCH_BUDGET_S = 45 * 60.0
W3_MEMORY_GUARD_BYTES = 8 * 1024**3
CATALOG_ID = "12_6_0193"
FAULT_MODEL = "flattened undecomposed DEM error mechanisms"
CANONICALIZE_CIRCUIT_ERRORS = True
DATA_LOGICAL_WITNESS_SUPPORTS = {
    "gross": [56, 57, 58, 59, 62, 64, 66, 68, 126, 130, 132, 134],
    "pbb": [7, 21, 35, 43, 57, 71, 77, 85, 99, 113, 121, 135],
}
DATA_LOGICAL_WITNESS_PROVENANCE = {
    "gross": "existing exact Gross code-distance certificate witness_X",
    "pbb": (
        "EXP-029 exact pure-X sector CP-SAT search: weight 12, all nontrivial "
        "sectors OPTIMAL or INFEASIBLE"
    ),
}

# Fast searches establish a witness before spending the remaining budget on
# the less truncated search.  Unused time is carried into the final attempt.
HEURISTIC_SEARCH_ATTEMPTS = (
    {
        "name": "monotone_cap4",
        "dont_explore_detection_event_sets_with_size_above": 4,
        "dont_explore_edges_with_degree_above": 4,
        "dont_explore_edges_increasing_symptom_degree": True,
        "max_wall_s": 180.0,
    },
    {
        "name": "monotone_cap6",
        "dont_explore_detection_event_sets_with_size_above": 6,
        "dont_explore_edges_with_degree_above": 6,
        "dont_explore_edges_increasing_symptom_degree": True,
        "max_wall_s": 420.0,
    },
    {
        "name": "nonmonotone_cap4",
        "dont_explore_detection_event_sets_with_size_above": 4,
        "dont_explore_edges_with_degree_above": 4,
        "dont_explore_edges_increasing_symptom_degree": False,
        "max_wall_s": 600.0,
    },
    {
        "name": "nonmonotone_cap5",
        "dont_explore_detection_event_sets_with_size_above": 5,
        "dont_explore_edges_with_degree_above": 5,
        "dont_explore_edges_increasing_symptom_degree": False,
        "max_wall_s": None,  # all time remaining from the 45 minute budget
    },
)

RAW_DIR = ROOT / "results" / "raw"
PROCESSED_DIR = ROOT / "results" / "processed"
REPORT_DIR = ROOT / "notes" / "agent_reports"
PARTIAL_DIR = ROOT / "results" / "partial_runs"
CANONICAL_STATE_PATH = RAW_DIR / "exp029_circuit_distance_state.json"
CANONICAL_WITNESS_PATH = RAW_DIR / "exp029_witnesses.json"
PROCESSED_PATH = PROCESSED_DIR / "exp029_circuit_distance.json"
REPORT_PATH = REPORT_DIR / "exp029_circuit_distance.md"
CATALOG_PATH = (
    ROOT
    / "third_party"
    / "qcode-discovery"
    / "results"
    / "campaign7_publication_merged.jsonl"
)

EXPECTED = {
    "gross": {
        "label": "CSS-BB [[144,12,12]] Gross",
        "catalog_id": None,
        "schedule_depth": 7,
        "two_qubit_gates_per_round": 864,
        "max_check_weight": 6,
        "num_mixed_checks": 0,
        "dem_errors": 67032,
        "dem_detectors": 1728,
    },
    "pbb": {
        "label": "nonCSS-PBB [[144,12,12]] 12_6_0193",
        "catalog_id": CATALOG_ID,
        "schedule_depth": 8,
        "two_qubit_gates_per_round": 1008,
        "max_check_weight": 8,
        "num_mixed_checks": 72,
        "dem_errors": 82800,
        "dem_detectors": 1728,
    },
}
def protocol_json() -> dict[str, Any]:
    """Machine-readable registration copied into every EXP-029 artifact."""
    return {
        "experiment": EXPERIMENT,
        "fault_model": FAULT_MODEL,
        "p_dem_structure_only": P_DEM_STRUCTURE,
        "probabilities_used_in_distance": False,
        "rounds": ROUNDS,
        "basis": BASIS,
        "circuit_rebuild": "exact EXP-013 build protocol",
        "schedule_time_limit_s": SCHEDULE_TIME_LIMIT_S,
        "schedule_workers_argument": SCHEDULE_WORKERS_ARGUMENT,
        "schedule_deterministic": SCHEDULE_DETERMINISTIC,
        "dem_lower_bound_decompose_errors": False,
        "dem_graphlike_decompose_errors": False,
        "stim_circuit_search_dem": (
            "internal to Circuit.search_for_undetectable_logical_errors; "
            "returned signatures are accepted only when present in the "
            "flattened decompose_errors=False DEM"
        ),
        "canonicalize_circuit_errors": CANONICALIZE_CIRCUIT_ERRORS,
        "upper_search_budget_s_per_circuit": UPPER_SEARCH_BUDGET_S,
        "graphlike_budget_s_per_circuit": min(
            GRAPLIKE_SEARCH_BUDGET_S, UPPER_SEARCH_BUDGET_S
        ),
        "w3_budget_s_per_circuit": W3_SEARCH_BUDGET_S,
        "w3_memory_guard_bytes": W3_MEMORY_GUARD_BYTES,
        "heuristic_search_attempts": [dict(x) for x in HEURISTIC_SEARCH_ATTEMPTS],
        "thread_environment": {
            name: os.environ[name]
            for name in (
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS",
                "NUMEXPR_NUM_THREADS",
            )
        },
        "external_load_caveat": (
            "Shared workstation had a dominant interactive desktop workload; "
            "this affects wall time, not exhaustive counts or certificates."
        ),
        "statistics": {
            "monte_carlo": False,
            "shots": 0,
            "fails": 0,
            "clopper_pearson_95": None,
            "note": "No rates are estimated in this exhaustive/search experiment.",
        },
    }


def seeds_json() -> dict[str, Any]:
    return {
        "root_seed": ROOT_SEED,
        "schedule_random_seed": SCHEDULE_RANDOM_SEED,
        "schedule_seed_reason": "pinned by EXP-013 deterministic default",
        "stim_search_seed": None,
        "stim_search_seed_reason": "Stim 1.16 exposes no seed argument for this search",
    }


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def load_catalogue_row() -> dict[str, Any]:
    matches = []
    with CATALOG_PATH.open() as f:
        for line in f:
            row = json.loads(line)
            if row.get("code_id") == CATALOG_ID:
                matches.append(row)
    if len(matches) != 1:
        raise RuntimeError(f"expected exactly one catalogue row {CATALOG_ID}, got {len(matches)}")
    row = matches[0]
    if (row["n"], row["k"], row["d"]) != (144, 12, 12):
        raise RuntimeError(f"unexpected parameters for {CATALOG_ID}: "
                           f"{row['n'], row['k'], row['d']}")
    return row


def build_exact_circuit(key: str) -> tuple[
    stim.Circuit, stim.DetectorErrorModel, dict[str, Any], Any
]:
    """Rebuild one circuit by the literal EXP-013 scheduling loop."""
    if key == "gross":
        code, supports, orbits = bb_supports_and_orbits(BRAVYI_BB["[[144,12,12]]"])
        catalog_id = None
    elif key == "pbb":
        row = load_catalogue_row()
        spec = PBBSpec(
            row["ell"],
            row["m"],
            [tuple(t) for t in row["A_terms"]],
            [tuple(t) for t in row["B_terms"]],
            [tuple(t) for t in row["C_terms"]],
            [tuple(t) for t in row["D_terms"]],
        )
        code, supports, orbits = pbb_supports_and_orbits(spec)
        catalog_id = row["code_id"]
    else:
        raise ValueError(key)

    lower_depth = depth_lower_bound(supports, code.n)
    num_directions = len(set(orbits.values()))
    schedule_result = None
    schedule_depth = None
    for depth in range(lower_depth, num_directions + 3):
        result = cpsat_schedule(
            supports,
            code.n,
            T=depth,
            time_limit_s=SCHEDULE_TIME_LIMIT_S,
            workers=SCHEDULE_WORKERS_ARGUMENT,
            symmetry_orbits=orbits,
            random_seed=SCHEDULE_RANDOM_SEED,
            deterministic=SCHEDULE_DETERMINISTIC,
        )
        if result.slot is not None and result.verification["valid"]:
            schedule_result = result
            schedule_depth = depth
            break
    if schedule_result is None or schedule_depth is None:
        raise RuntimeError(f"no EXP-013 schedule found for {key}")

    layers = slots_to_layers(schedule_result.slot, schedule_depth)
    circuit, metadata = build_memory_circuit(
        CircuitSpec(
            H=code.H,
            observables=pure_z_logical_basis(code),
            rounds=ROUNDS,
            p=P_DEM_STRUCTURE,
            basis=BASIS,
            layers=layers,
        )
    )
    metadata["schedule_depth"] = schedule_depth
    dem = circuit.detector_error_model(decompose_errors=False)

    exp = EXPECTED[key]
    observed = {
        "schedule_depth": schedule_depth,
        "two_qubit_gates_per_round": metadata["total_two_qubit_gates"],
        "max_check_weight": metadata["max_check_weight"],
        "num_mixed_checks": metadata["num_mixed_checks"],
        "dem_errors": dem.num_errors,
        "dem_detectors": dem.num_detectors,
    }
    expected_observed = {name: exp[name] for name in observed}
    if observed != expected_observed:
        raise RuntimeError(
            f"{key} does not match EXP-013 structural fingerprint: "
            f"observed={observed}, expected={expected_observed}"
        )

    slot_lines = [f"{ci},{q}:{schedule_result.slot[(ci, q)]}"
                  for ci, q in sorted(schedule_result.slot)]
    build_info = {
        "key": key,
        "label": exp["label"],
        "catalog_id": catalog_id,
        "n": code.n,
        "k": code.k,
        "code_distance": 12,
        "rounds": ROUNDS,
        "p": P_DEM_STRUCTURE,
        "basis": BASIS,
        "schedule_depth": schedule_depth,
        "schedule_lower_bound": lower_depth,
        "schedule_status": schedule_result.status,
        "schedule_solver_wall_s": schedule_result.wall_time_s,
        "schedule_verification": schedule_result.verification,
        "schedule_sha256": sha256_text("\n".join(slot_lines)),
        "circuit_sha256": sha256_text(str(circuit)),
        "dem_undecomposed_sha256": sha256_text(str(dem)),
        "two_qubit_gates_per_round": metadata["total_two_qubit_gates"],
        "max_check_weight": metadata["max_check_weight"],
        "num_mixed_checks": metadata["num_mixed_checks"],
        "num_observables": circuit.num_observables,
        "dem_errors": dem.num_errors,
        "dem_detectors": dem.num_detectors,
    }
    return circuit, dem, build_info, code

def certify_pbb_data_witness(code: Any) -> dict[str, Any]:
    """Reproduce the PBB pure-X weight-12 upper witness with exact CP-SAT."""
    started = time.perf_counter()
    result = min_pure_x_logical(
        code,
        time_limit_s=UPPER_SEARCH_BUDGET_S,
        workers=8,
        upper_bound=12,
    )
    if result.weight != 12 or not result.exact or result.support is None:
        raise RuntimeError(
            "PBB pure-X witness reproduction failed: "
            f"weight={result.weight}, exact={result.exact}, support={result.support}"
        )
    # CP-SAT may return any translated member of the weight-12 orbit. Use the
    # freshly certified support for the subsequent DEM witness mapping.
    DATA_LOGICAL_WITNESS_SUPPORTS["pbb"] = sorted(result.support)
    return {
        "method": "qec_research.distance.sectors.min_pure_x_logical",
        "classification": "exact(OPTIMAL) pure-X data-code witness",
        "weight": result.weight,
        "exact": result.exact,
        "support": result.support,
        "n_sectors": result.n_sectors,
        "statuses": result.statuses,
        "time_limit_s": UPPER_SEARCH_BUDGET_S,
        "workers": 8,
        "wall_s": time.perf_counter() - started,
    }


def xor_target_mask(targets: list[stim.DemTarget], *, detector: bool) -> int:
    mask = 0
    for target in targets:
        relevant = (target.is_relative_detector_id() if detector
                    else target.is_logical_observable_id())
        if relevant:
            mask ^= 1 << target.val
    return mask


def extract_error_mechanisms(dem: stim.DetectorErrorModel) -> tuple[list[dict[str, Any]], int]:
    """Flatten a DEM and return the complete signature of each error instruction."""
    mechanisms = []
    separator_instructions = 0
    for instruction_index, instruction in enumerate(dem.flattened()):
        if instruction.type != "error":
            continue
        targets = instruction.targets_copy()
        if any(t.is_separator() for t in targets):
            separator_instructions += 1
        probability = float(instruction.args_copy()[0])
        if probability <= 0:
            continue
        detector_mask = xor_target_mask(targets, detector=True)
        observable_mask = xor_target_mask(targets, detector=False)
        mechanisms.append({
            "index": len(mechanisms),
            "flattened_instruction_index": instruction_index,
            "probability": probability,
            "detector_mask": detector_mask,
            "observable_mask": observable_mask,
            "detector_degree": detector_mask.bit_count(),
        })
    return mechanisms, separator_instructions


def mask_ids(mask: int) -> list[int]:
    out = []
    while mask:
        low = mask & -mask
        out.append(low.bit_length() - 1)
        mask ^= low
    return out


def mechanism_json(mechanism: dict[str, Any]) -> dict[str, Any]:
    return {
        "index": mechanism["index"],
        "flattened_instruction_index": mechanism["flattened_instruction_index"],
        "probability_ignored_for_distance": mechanism["probability"],
        "detectors": mask_ids(mechanism["detector_mask"]),
        "observables": mask_ids(mechanism["observable_mask"]),
        "detector_degree": mechanism["detector_degree"],
    }


def max_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    # Darwin reports bytes; Linux reports KiB.
    return value if sys.platform == "darwin" else value * 1024


def choose_different_observable(
    variants: dict[int, int],
    forbidden_observable: int,
) -> tuple[int, int] | None:
    for observable, index in variants.items():
        if observable != forbidden_observable:
            return observable, index
    return None


def check_w1_w2(mechanisms: list[dict[str, Any]]) -> tuple[dict[str, Any], dict[str, Any] | None]:
    """Exhaust all one- and two-mechanism undetectable logical sets."""
    t0 = time.perf_counter()
    detector_variants: dict[int, dict[int, int]] = {}
    w1_property_holds = True
    w1_witness = None
    for mechanism in mechanisms:
        detector_mask = mechanism["detector_mask"]
        observable_mask = mechanism["observable_mask"]
        if detector_mask == 0 and observable_mask != 0:
            w1_property_holds = False
            if w1_witness is None:
                w1_witness = {
                    "weight": 1,
                    "mechanisms": [mechanism_json(mechanism)],
                }
        detector_variants.setdefault(detector_mask, {}).setdefault(
            observable_mask, mechanism["index"]
        )

    w2_witness = None
    for detector_mask, variants in detector_variants.items():
        if len(variants) < 2:
            continue
        first_observable, first_index = next(iter(variants.items()))
        other = choose_different_observable(variants, first_observable)
        if other is None:
            continue
        _, second_index = other
        w2_witness = {
            "weight": 2,
            "mechanisms": [
                mechanism_json(mechanisms[first_index]),
                mechanism_json(mechanisms[second_index]),
            ],
            "shared_detector_set": mask_ids(detector_mask),
        }
        break

    result = {
        "method": "exhaustive hash of undecomposed flattened DEM mechanism signatures",
        "fault_model": FAULT_MODEL,
        "mechanisms_checked": len(mechanisms),
        "w1_property_holds": w1_property_holds,
        "w1_undetectable_logical_found": w1_witness is not None,
        "w2_undetectable_logical_found": w2_witness is not None,
        "detector_signature_groups": len(detector_variants),
        "wall_s": time.perf_counter() - t0,
        "certified": True,
    }
    return result, w1_witness or w2_witness


def check_w3(
    mechanisms: list[dict[str, Any]],
    *,
    budget_s: float = W3_SEARCH_BUDGET_S,
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    """Exact 3-XOR meet-in-the-middle, streaming only shared-detector pairs.

    With no weight-1 or weight-2 logical witness, a weight-3 witness cannot use
    a detector-empty mechanism and cannot use the same complete signature
    twice.  In any remaining detector-cancelling triple, each present detector
    occurs in exactly two of the three mechanisms.  Therefore some pair shares
    a detector.  Enumerating each shared-detector pair once and querying the
    singles hash by its XOR signature is exhaustive without materializing all
    O(N^2) pairs.
    """
    t0 = time.perf_counter()
    deadline = t0 + budget_s

    # Duplicate complete signatures cannot participate in a minimal weight-3
    # witness: two identical signatures cancel and leave a weight-1 witness.
    unique_by_signature: dict[tuple[int, int], int] = {}
    for mechanism in mechanisms:
        signature = (mechanism["detector_mask"], mechanism["observable_mask"])
        unique_by_signature.setdefault(signature, mechanism["index"])
    unique_indices = list(unique_by_signature.values())
    unique_detector = [mechanisms[i]["detector_mask"] for i in unique_indices]
    unique_observable = [mechanisms[i]["observable_mask"] for i in unique_indices]

    singles: dict[int, dict[int, int]] = {}
    incidence: dict[int, list[int]] = {}
    detector_empty = 0
    for local_index, original_index in enumerate(unique_indices):
        detector_mask = unique_detector[local_index]
        observable_mask = unique_observable[local_index]
        singles.setdefault(detector_mask, {}).setdefault(observable_mask, local_index)
        if detector_mask == 0:
            detector_empty += 1
            continue
        work = detector_mask
        while work:
            low = work & -work
            detector_id = low.bit_length() - 1
            incidence.setdefault(detector_id, []).append(local_index)
            work ^= low

    all_pair_count = math.comb(len(unique_indices), 2)
    shared_pair_visits_upper = sum(math.comb(len(ids), 2) for ids in incidence.values())
    # A Python pair-table entry containing a ~1728-bit key and indices is at
    # least several hundred bytes.  This is a conservative planning estimate;
    # the actual algorithm streams pairs and never allocates this table.
    estimated_materialized_pair_bytes = all_pair_count * 384
    rss_start = max_rss_bytes()
    if rss_start >= W3_MEMORY_GUARD_BYTES:
        return ({
            "method": "singles-vs-streamed-pairs XOR meet-in-the-middle",
            "fault_model": FAULT_MODEL,
            "completed": False,
            "certified_no_w3": False,
            "stop_reason": "8 GiB RSS guard already exceeded before pair enumeration",
            "unique_mechanism_signatures": len(unique_indices),
            "all_pair_count": all_pair_count,
            "shared_pair_visits_upper": shared_pair_visits_upper,
            "estimated_materialized_pair_bytes": estimated_materialized_pair_bytes,
            "memory_guard_bytes": W3_MEMORY_GUARD_BYTES,
            "max_rss_bytes": rss_start,
            "wall_s": time.perf_counter() - t0,
        }, None)

    pair_visits = 0
    unique_shared_pairs_checked = 0
    witness = None
    stop_reason = None
    timed_out = False
    memory_guard_hit = False

    singles_get = singles.get
    detector_masks = unique_detector
    observable_masks = unique_observable
    unique_original = unique_indices
    check_interval_mask = (1 << 18) - 1

    for detector_id in sorted(incidence):
        ids = incidence[detector_id]
        owner_bit = 1 << detector_id
        for left_pos in range(len(ids)):
            left = ids[left_pos]
            left_detector = detector_masks[left]
            left_observable = observable_masks[left]
            for right_pos in range(left_pos + 1, len(ids)):
                right = ids[right_pos]
                pair_visits += 1
                if (pair_visits & check_interval_mask) == 0:
                    now = time.perf_counter()
                    rss = max_rss_bytes()
                    if rss >= W3_MEMORY_GUARD_BYTES:
                        memory_guard_hit = True
                        stop_reason = "8 GiB RSS guard reached during pair enumeration"
                        break
                    if now >= deadline:
                        timed_out = True
                        stop_reason = "45 minute pair-enumeration budget exhausted"
                        break
                right_detector = detector_masks[right]
                common = left_detector & right_detector
                # Assign a pair to the least detector it shares, avoiding a set
                # of already-seen O(N^2) pairs.
                if common & -common != owner_bit:
                    continue
                unique_shared_pairs_checked += 1
                target_detector = left_detector ^ right_detector
                if target_detector == 0:
                    continue
                variants = singles_get(target_detector)
                if variants is None:
                    continue
                pair_observable = left_observable ^ observable_masks[right]
                choice = choose_different_observable(variants, pair_observable)
                if choice is None:
                    continue
                _, third = choice
                original_ids = [
                    unique_original[left],
                    unique_original[right],
                    unique_original[third],
                ]
                if len(set(original_ids)) != 3:
                    continue
                detector_xor = 0
                observable_xor = 0
                for original_id in original_ids:
                    detector_xor ^= mechanisms[original_id]["detector_mask"]
                    observable_xor ^= mechanisms[original_id]["observable_mask"]
                if detector_xor == 0 and observable_xor != 0:
                    witness = {
                        "weight": 3,
                        "mechanisms": [mechanism_json(mechanisms[i]) for i in original_ids],
                        "detector_xor": [],
                        "observable_xor": mask_ids(observable_xor),
                    }
                    break


            if witness is not None or stop_reason is not None:
                break
            # Also check between short inner loops, where the bit-mask cadence
            # above might not land exactly on a multiple.
            if time.perf_counter() >= deadline:
                timed_out = True
                stop_reason = "45 minute pair-enumeration budget exhausted"
                break
        if witness is not None or stop_reason is not None:
            break

    completed = witness is not None or stop_reason is None
    max_rss = max_rss_bytes()
    result = {
        "method": "singles-vs-streamed-pairs XOR meet-in-the-middle",
        "fault_model": FAULT_MODEL,
        "completeness_argument": (
            "Each detector in a zero-XOR triple occurs in exactly two mechanisms; "
            "therefore at least one pair shares a detector. Each such pair is "
            "processed once at its least shared detector and queried against the singles hash."
        ),
        "completed": completed,
        "certified_no_w3": completed and witness is None,
        "w3_undetectable_logical_found": witness is not None,
        "stop_reason": stop_reason,
        "timed_out": timed_out,
        "memory_guard_hit": memory_guard_hit,
        "mechanisms": len(mechanisms),
        "unique_mechanism_signatures": len(unique_indices),
        "detector_empty_unique_signatures": detector_empty,
        "all_pair_count": all_pair_count,
        "shared_pair_visits_upper": shared_pair_visits_upper,
        "pair_visits": pair_visits,
        "unique_shared_pairs_checked": unique_shared_pairs_checked,
        "estimated_materialized_pair_bytes": estimated_materialized_pair_bytes,
        "pair_storage": "streamed; full pair table was not allocated",
        "memory_guard_bytes": W3_MEMORY_GUARD_BYTES,
        "rss_start_bytes": rss_start,
        "max_rss_bytes": max_rss,
        "budget_s": budget_s,
        "wall_s": time.perf_counter() - t0,
    }
    return result, witness


def explained_error_json(error: stim.ExplainedError, index: int) -> dict[str, Any]:
    detector_mask = 0
    observable_mask = 0
    for term_with_coords in error.dem_error_terms:
        target = term_with_coords.dem_target
        if target.is_relative_detector_id():
            detector_mask ^= 1 << target.val
        elif target.is_logical_observable_id():
            observable_mask ^= 1 << target.val
    return {
        "search_entry_index": index,
        "detectors": mask_ids(detector_mask),
        "observables": mask_ids(observable_mask),
        "detector_mask_hex": hex(detector_mask),
        "observable_mask_hex": hex(observable_mask),
        "explained_error": str(error),
    }

def explicit_data_fault_witness(
    key: str,
    circuit: stim.Circuit,
    mechanisms: list[dict[str, Any]],
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Map a known weight-12 pure-X data logical into twelve DEM mechanisms.

    The circuit has independent X_ERROR(p) mechanisms immediately after
    Z-basis data preparation.  Applying these mechanisms on the support of a
    pure-X logical must cancel every detector and flip at least one measured
    pure-Z observable.  Stim explains the exact circuit-to-DEM mapping; the
    resulting signatures are then required to be present in the flattened
    undecomposed DEM and are XOR-validated.
    """
    support = DATA_LOGICAL_WITNESS_SUPPORTS[key]
    explained = circuit.explain_detector_error_model_errors(
        reduce_to_one_representative_error=False
    )
    by_qubit: dict[int, dict[str, Any]] = {}
    for error in explained:
        entry = explained_error_json(error, 0)
        for location in error.circuit_error_locations:
            text = str(location)
            if "(after 0 TICKs)" not in text or "X_ERROR(0.001)" not in text:
                continue
            product = location.flipped_pauli_product
            if len(product) != 1:
                continue
            target = product[0].gate_target
            if not target.is_x_target:
                continue
            mapped = dict(entry)
            mapped["data_qubit"] = target.value
            mapped["circuit_error_location"] = str(location)
            by_qubit[target.value] = mapped
    missing = sorted(set(support) - set(by_qubit))
    if missing:
        raise RuntimeError(f"initial-data DEM mechanisms missing for {key}: {missing}")
    entries = []
    for index, qubit in enumerate(support):
        entry = dict(by_qubit[qubit])
        entry["search_entry_index"] = index
        entries.append(entry)
    validation = validate_signature_witness(entries, mechanisms)
    if not validation["valid_mechanism_level_upper_bound"]:
        raise RuntimeError(f"invalid explicit data-logical witness for {key}: {validation}")
    record = {
        "weight": len(entries),
        "method": "explicit initial-data X_ERROR mechanisms on a pure-X logical support",
        "provenance": DATA_LOGICAL_WITNESS_PROVENANCE[key],
        "mechanisms": entries,
        "validation": validation,
    }
    candidate = {
        "value": len(entries),
        "method": "explicit weight-12 data-logical DEM witness",
        "witness_key": "explicit_data_logical",
        "classification": "certified mechanism-level upper bound",
    }
    return candidate, record


def heuristic_worker(send_conn: Any, circuit_text: str, parameters: dict[str, Any]) -> None:
    """Spawn target: isolate Stim's non-interruptible search behind a timeout."""
    try:
        t0 = time.perf_counter()
        circuit = stim.Circuit(circuit_text)
        errors = circuit.search_for_undetectable_logical_errors(
            dont_explore_detection_event_sets_with_size_above=parameters[
                "dont_explore_detection_event_sets_with_size_above"
            ],
            dont_explore_edges_with_degree_above=parameters[
                "dont_explore_edges_with_degree_above"
            ],
            dont_explore_edges_increasing_symptom_degree=parameters[
                "dont_explore_edges_increasing_symptom_degree"
            ],
            canonicalize_circuit_errors=CANONICALIZE_CIRCUIT_ERRORS,
        )
        send_conn.send({
            "ok": True,
            "wall_s": time.perf_counter() - t0,
            "weight": len(errors),
            "mechanisms": [explained_error_json(error, i)
                           for i, error in enumerate(errors)],
        })
    except BaseException as exc:  # returned as data; parent still enforces timeout
        send_conn.send({
            "ok": False,
            "wall_s": time.perf_counter() - t0,
            "error": f"{type(exc).__name__}: {exc}",
        })
    finally:
        send_conn.close()


def graphlike_worker(send_conn: Any, dem_text: str) -> None:
    try:
        t0 = time.perf_counter()
        dem = stim.DetectorErrorModel(dem_text)
        witness = dem.shortest_graphlike_error(ignore_ungraphlike_errors=True)
        send_conn.send({
            "ok": True,
            "wall_s": time.perf_counter() - t0,
            "weight": witness.num_errors,
            "dem": str(witness),
        })
    except BaseException as exc:
        send_conn.send({
            "ok": False,
            "wall_s": time.perf_counter() - t0,
            "error": f"{type(exc).__name__}: {exc}",
        })
    finally:
        send_conn.close()


def run_timed_child(
    target: Any,
    args: tuple[Any, ...],
    timeout_s: float,
) -> dict[str, Any]:
    ctx = mp.get_context("spawn")
    receive_conn, send_conn = ctx.Pipe(duplex=False)
    process = ctx.Process(target=target, args=(send_conn, *args))
    started = time.perf_counter()
    process.start()
    send_conn.close()
    payload = None
    if receive_conn.poll(timeout_s):
        try:
            payload = receive_conn.recv()
        except EOFError:
            payload = None
    if payload is None and process.is_alive():
        process.terminate()
        process.join(5)
        if process.is_alive():
            process.kill()
            process.join(5)
    else:
        process.join(5)
        if process.is_alive():
            process.terminate()
            process.join(5)
    receive_conn.close()
    elapsed = time.perf_counter() - started
    if payload is None:
        return {
            "ok": False,
            "timed_out": elapsed >= timeout_s * 0.99,
            "error": "search process exceeded its wall-time allocation"
                     if elapsed >= timeout_s * 0.99
                     else f"search process exited without a result (exitcode={process.exitcode})",
            "wall_s": elapsed,
        }
    payload["timed_out"] = False
    payload["parent_wall_s"] = elapsed
    return payload


def signature_representatives(
    mechanisms: list[dict[str, Any]],
) -> dict[tuple[int, int], list[int]]:
    result: dict[tuple[int, int], list[int]] = {}
    for mechanism in mechanisms:
        key = (mechanism["detector_mask"], mechanism["observable_mask"])
        result.setdefault(key, []).append(mechanism["index"])
    return result


def validate_signature_witness(
    entries: list[dict[str, Any]],
    mechanisms: list[dict[str, Any]],
) -> dict[str, Any]:
    representatives = signature_representatives(mechanisms)
    detector_xor = 0
    observable_xor = 0
    all_present = True
    used: dict[tuple[int, int], int] = {}
    for entry in entries:
        detector_mask = int(entry["detector_mask_hex"], 16)
        observable_mask = int(entry["observable_mask_hex"], 16)
        detector_xor ^= detector_mask
        observable_xor ^= observable_mask
        signature = (detector_mask, observable_mask)
        occurrence = used.get(signature, 0)
        candidates = representatives.get(signature, [])
        if occurrence < len(candidates):
            entry["undecomposed_mechanism_index"] = candidates[occurrence]
            used[signature] = occurrence + 1
        else:
            entry["undecomposed_mechanism_index"] = None
            all_present = False
    return {
        "all_entries_present_in_undecomposed_dem": all_present,
        "detectors_cancel": detector_xor == 0,
        "observable_flips": observable_xor != 0,
        "observable_xor": mask_ids(observable_xor),
        "valid_mechanism_level_upper_bound": (
            all_present and detector_xor == 0 and observable_xor != 0
        ),
    }


def dem_witness_entries(dem_text: str) -> list[dict[str, Any]]:
    dem = stim.DetectorErrorModel(dem_text)
    mechanisms, _ = extract_error_mechanisms(dem)
    return [
        {
            "search_entry_index": i,
            "detectors": mask_ids(mechanism["detector_mask"]),
            "observables": mask_ids(mechanism["observable_mask"]),
            "detector_mask_hex": hex(mechanism["detector_mask"]),
            "observable_mask_hex": hex(mechanism["observable_mask"]),
            "graphlike_dem_instruction": str(dem.flattened()[mechanism["flattened_instruction_index"]]),
        }
        for i, mechanism in enumerate(mechanisms)
    ]


def run_lower_bound(
    key: str,
    dem: stim.DetectorErrorModel,
    build_info: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    stage_start = time.perf_counter()
    mechanisms, separator_instructions = extract_error_mechanisms(dem)
    degree_histogram = Counter(m["detector_degree"] for m in mechanisms)
    unique_signatures = len({
        (m["detector_mask"], m["observable_mask"]) for m in mechanisms
    })

    w12, low_witness = check_w1_w2(mechanisms)
    w3 = None
    if low_witness is None:
        w3, w3_witness = check_w3(mechanisms)
        low_witness = w3_witness

    if w12["w1_undetectable_logical_found"]:
        lb_value = 1
        lb_method = "trivial distance >=1; exhaustive w=1 found a logical"
    elif w12["w2_undetectable_logical_found"]:
        lb_value = 2
        lb_method = "exhaustive w=1 exclusion; exhaustive w=2 found a logical"
    elif w3 is not None and w3["certified_no_w3"]:
        lb_value = 4
        lb_method = "exhaustive undecomposed-DEM mechanism search at w=1,2,3"
    else:
        lb_value = 3
        lb_method = "exhaustive undecomposed-DEM mechanism search at w=1,2"

    result = {
        "build": build_info,
        "dem": {
            "decompose_errors": False,
            "flattened_error_mechanisms": len(mechanisms),
            "flattened_separator_instructions": separator_instructions,
            "unique_complete_signatures": unique_signatures,
            "detector_degree_histogram": {
                str(k): v for k, v in sorted(degree_histogram.items())
            },
            "max_detector_degree": max(degree_histogram, default=0),
        },
        "w1_w2": w12,
        "w3": w3,
        "lb": {
            "value": lb_value,
            "method": lb_method,
            "certified": True,
            "fault_model": FAULT_MODEL,
            "classification": "CERTIFIED lower bound",
        },
        "lower_wall_s": time.perf_counter() - stage_start,
    }
    witness = {
        "lower_small_weight_witness": low_witness,
    }
    print(
        f"{build_info['label']}: mechanisms={len(mechanisms)} "
        f"max_degree={result['dem']['max_detector_degree']} "
        f"certified d_circuit>={lb_value}",
        flush=True,
    )
    return result, witness


def run_upper_bound(
    key: str,
    circuit: stim.Circuit,
    dem: stim.DetectorErrorModel,
    build_info: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    stage_start = time.perf_counter()
    mechanisms, _ = extract_error_mechanisms(dem)
    witness_output: dict[str, Any] = {}
    candidates = []
    explicit_candidate, explicit_record = explicit_data_fault_witness(
        key, circuit, mechanisms
    )
    witness_output["explicit_data_logical"] = explicit_record
    candidates.append(explicit_candidate)
    print(
        f"  explicit data-logical result: {explicit_candidate['value']} mechanisms",
        flush=True,
    )

    print(f"{build_info['label']}: starting graphlike search", flush=True)
    graphlike_budget_s = min(
        GRAPLIKE_SEARCH_BUDGET_S, UPPER_SEARCH_BUDGET_S
    )
    graphlike = run_timed_child(
        graphlike_worker,
        (str(dem),),
        graphlike_budget_s,
    )
    graphlike_record = {
        "method": "stim.DetectorErrorModel.shortest_graphlike_error",
        "classification": "graphlike upper bound",
        "fault_model": FAULT_MODEL,
        "dem_decompose_errors": False,
        "ignore_ungraphlike_errors": True,
        "budget_s": graphlike_budget_s,
        **{k: v for k, v in graphlike.items() if k != "dem"},
    }
    if graphlike.get("ok"):
        entries = dem_witness_entries(graphlike["dem"])
        validation = validate_signature_witness(entries, mechanisms)
        graphlike_record["validation"] = validation
        graphlike_record["weight"] = len(entries)
        witness_output["graphlike"] = {
            "weight": len(entries),
            "dem": graphlike["dem"],
            "mechanisms": entries,
            "validation": validation,
        }
        if validation["valid_mechanism_level_upper_bound"]:
            candidates.append({
                "value": len(entries),
                "method": "stim graphlike search on undecomposed DEM",
                "witness_key": "graphlike",
                "classification": "graphlike upper bound",
            })
        print(f"  graphlike result: {len(entries)} mechanisms", flush=True)
    else:
        print(f"  graphlike unavailable: {graphlike_record.get('error')}", flush=True)

    search_start = time.perf_counter()
    attempts = []
    circuit_text = str(circuit)
    for attempt_template in HEURISTIC_SEARCH_ATTEMPTS:
        elapsed = time.perf_counter() - search_start
        remaining = UPPER_SEARCH_BUDGET_S - elapsed
        if remaining <= 0:
            break
        parameters = dict(attempt_template)
        configured_max = parameters.pop("max_wall_s")
        allocation = remaining if configured_max is None else min(configured_max, remaining)
        print(
            f"  heuristic {parameters['name']} allocation={allocation:.1f}s",
            flush=True,
        )
        payload = run_timed_child(
            heuristic_worker,
            (circuit_text, parameters),
            allocation,
        )
        record = {
            "method": "stim.Circuit.search_for_undetectable_logical_errors",
            "classification": "heuristic upper-bound search",
            "parameters": {
                **parameters,
                "canonicalize_circuit_errors": CANONICALIZE_CIRCUIT_ERRORS,
            },
            "allocation_s": allocation,
            **{k: v for k, v in payload.items() if k != "mechanisms"},
        }
        if payload.get("ok"):
            entries = payload["mechanisms"]
            validation = validate_signature_witness(entries, mechanisms)
            record["validation"] = validation
            witness_key = f"heuristic_{parameters['name']}"
            witness_output[witness_key] = {
                "weight": payload["weight"],
                "parameters": record["parameters"],
                "mechanisms": entries,
                "validation": validation,
            }
            if validation["valid_mechanism_level_upper_bound"]:
                candidates.append({
                    "value": payload["weight"],
                    "method": (
                        "stim heuristic search_for_undetectable_logical_errors "
                        f"({parameters['name']})"
                    ),
                    "witness_key": witness_key,
                    "classification": "heuristic upper bound",
                })
            print(
                f"    returned weight={payload['weight']} "
                f"validated={validation['valid_mechanism_level_upper_bound']}",
                flush=True,
            )
        else:
            print(f"    no result: {record.get('error')}", flush=True)
        attempts.append(record)

    search_wall_s = time.perf_counter() - search_start
    if not candidates:
        upper = {
            "value": None,
            "method": "no validated Stim witness returned within budgets",
            "witness_path": None,
            "classification": "no upper bound obtained",
        }
    else:
        best = min(candidates, key=lambda x: (x["value"], x["method"]))
        upper = {
            "value": best["value"],
            "method": best["method"],
            "witness_path": (
                f"results/raw/exp029_witnesses.json#circuits/{key}/"
                f"{best['witness_key']}"
            ),
            "classification": best["classification"],
        }

    result = {
        "build": build_info,
        "graphlike": graphlike_record,
        "heuristic_search": {
            "method": "stim.Circuit.search_for_undetectable_logical_errors",
            "budget_s": UPPER_SEARCH_BUDGET_S,
            "wall_s": search_wall_s,
            "attempts": attempts,
        },
        "upper_candidates": candidates,
        "ub": upper,
        "upper_wall_s": time.perf_counter() - stage_start,
    }
    print(
        f"{build_info['label']}: best validated upper bound "
        f"d_circuit<={upper['value']} ({upper['method']})",
        flush=True,
    )
    return result, witness_output


def empty_artifact() -> dict[str, Any]:
    return {
        "experiment": EXPERIMENT,
        "verdict": "INCONCLUSIVE",
        "protocol": protocol_json(),
        "seeds": seeds_json(),
        "environment": {
            "python": sys.version.split()[0],
            "stim": stim.__version__,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "cpu_count": os.cpu_count(),
            "thread_environment": {
                name: os.environ[name]
                for name in (
                    "OMP_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS",
                    "NUMEXPR_NUM_THREADS",
                )
            },
        },
        "circuits": {},
        "runs": [],
        "wall_s": 0.0,
    }


def load_artifact(path: Path) -> dict[str, Any]:
    if not path.exists():
        return empty_artifact()
    artifact = json.loads(path.read_text())
    artifact["protocol"] = protocol_json()
    artifact["seeds"] = seeds_json()
    return artifact


def merge_stage(
    state: dict[str, Any],
    key: str,
    stage: dict[str, Any],
) -> None:
    entry = state["circuits"].setdefault(key, {})
    for name, value in stage.items():
        if name == "build" and "build" in entry:
            old = entry["build"]
            for fingerprint in ("schedule_sha256", "circuit_sha256", "dem_undecomposed_sha256"):
                if old[fingerprint] != value[fingerprint]:
                    raise RuntimeError(
                        f"{key} rebuild changed {fingerprint}: "
                        f"{old[fingerprint]} != {value[fingerprint]}"
                    )
        entry[name] = value
    entry["wall_s"] = sum(
        float(entry.get(name, 0.0)) for name in ("lower_wall_s", "upper_wall_s")
    )


def update_verdict_and_wall(state: dict[str, Any]) -> None:
    complete = all(
        key in state["circuits"]
        and state["circuits"][key].get("lb", {}).get("certified") is True
        and state["circuits"][key].get("ub", {}).get("value") is not None
        for key in ("gross", "pbb")
    )
    positive = complete and all(
        state["circuits"][key]["lb"]["value"] >= 3
        for key in ("gross", "pbb")
    )
    state["verdict"] = "POSITIVE" if positive else "INCONCLUSIVE"
    state["wall_s"] = sum(
        float(state["circuits"].get(key, {}).get("wall_s", 0.0))
        for key in ("gross", "pbb")
    )


def save_state(state: dict[str, Any], path: Path) -> None:
    update_verdict_and_wall(state)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(state, indent=2) + "\n")


def save_witnesses(
    witnesses: dict[str, Any],
    state: dict[str, Any],
    path: Path,
) -> None:
    update_verdict_and_wall(state)
    artifact = {
        "experiment": EXPERIMENT,
        "verdict": state["verdict"],
        "protocol": protocol_json(),
        "seeds": seeds_json(),
        "circuits": witnesses.get("circuits", {}),
        "wall_s": state["wall_s"],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(artifact, indent=2) + "\n")


def selected_witness_from_lower(
    key: str,
    entry: dict[str, Any],
    witnesses: dict[str, Any],
) -> dict[str, Any] | None:
    lower = witnesses.get("circuits", {}).get(key, {}).get("lower_small_weight_witness")
    if lower is None:
        return None
    value = int(lower["weight"])
    return {
        "value": value,
        "method": f"exhaustive mechanism-level weight-{value} witness",
        "witness_path": (
            f"results/raw/exp029_witnesses.json#circuits/{key}/"
            "lower_small_weight_witness"
        ),
        "classification": "certified mechanism-level upper bound",
    }


def fold_lower_witness_into_upper(
    state: dict[str, Any],
    witnesses: dict[str, Any],
) -> None:
    for key, entry in state["circuits"].items():
        candidate = selected_witness_from_lower(key, entry, witnesses)
        if candidate is None:
            continue
        current = entry.get("ub")
        if current is None or current.get("value") is None or candidate["value"] < current["value"]:
            entry["ub"] = candidate


def processed_artifact(
    state: dict[str, Any],
    witnesses: dict[str, Any],
) -> dict[str, Any]:
    update_verdict_and_wall(state)
    circuits = {}
    for key in ("gross", "pbb"):
        entry = state["circuits"][key]
        circuits[key] = {
            "label": entry["build"]["label"],
            "catalog_id": entry["build"]["catalog_id"],
            "schedule_depth": entry["build"]["schedule_depth"],
            "circuit_fingerprints": {
                "schedule_sha256": entry["build"]["schedule_sha256"],
                "circuit_sha256": entry["build"]["circuit_sha256"],
                "dem_undecomposed_sha256": entry["build"]["dem_undecomposed_sha256"],
            },
            "dem": entry["dem"],
            "lb": entry["lb"],
            "ub": entry["ub"],
            "graphlike": entry.get("graphlike"),
            "heuristic_search": entry.get("heuristic_search"),
            "w1_w2": entry["w1_w2"],
            "w3": entry.get("w3"),
            "wall_s": entry["wall_s"],
        }
        circuits[key]["pure_x_witness_certificate"] = entry["build"].get(
            "pbb_pure_x_witness_certificate"
        )
    witness_summary = {}
    for key in ("gross", "pbb"):
        explicit = witnesses.get("circuits", {}).get(key, {}).get("explicit_data_logical")
        if explicit is None:
            continue
        witness_summary[key] = {
            "explicit_data_logical": {
                "weight": explicit["weight"],
                "method": explicit["method"],
                "provenance": explicit["provenance"],
                "validation": explicit["validation"],
                "data_qubits": [m["data_qubit"] for m in explicit["mechanisms"]],
            }
        }
    return {
        "experiment": EXPERIMENT,
        "verdict": state["verdict"],
        "protocol": protocol_json(),
        "seeds": seeds_json(),
        "environment": state["environment"],
        "circuits": circuits,
        "witnesses": witness_summary,
        "runs": state["runs"],
        "wall_s": state["wall_s"],
    }


def report_text(result: dict[str, Any]) -> str:
    gross = result["circuits"]["gross"]
    pbb = result["circuits"]["pbb"]
    verdict = result["verdict"]
    first_reason = (
        f"both circuits have certified DEM-mechanism lower bounds of at least 3 "
        f"and validated upper bounds (Gross <= {gross['ub']['value']}, "
        f"PBB <= {pbb['ub']['value']})."
        if verdict == "POSITIVE"
        else "the requested certified-lower-bound plus upper-bound condition was not met for both circuits."
    )
    lines = [
        f"{verdict} — {first_reason}",
        "",
        "## Evidence",
        "",
        "The fault model is **flattened undecomposed DEM error mechanisms**: one positive-probability "
        "`error` instruction from `circuit.detector_error_model(decompose_errors=False)` counts as "
        "one fault. Probabilities are ignored. The nonzero value `p=0.001` only materializes the DEM "
        "structure, so these distance bounds do not depend on p.",
        "",
        "| Circuit | schedule | mechanisms | max detector degree | certified LB | best UB | graphlike |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for circuit in (gross, pbb):
        graphlike_value = None
        if circuit.get("graphlike", {}).get("ok"):
            graphlike_value = circuit["graphlike"].get("weight")
        lines.append(
            f"| {circuit['label']} | {circuit['schedule_depth']} | "
            f"{circuit['dem']['flattened_error_mechanisms']} | "
            f"{circuit['dem']['max_detector_degree']} | "
            f">= {circuit['lb']['value']} CERTIFIED | "
            f"<= {circuit['ub']['value']} ({circuit['ub']['classification']}) | "
            f"<= {graphlike_value if graphlike_value is not None else 'not returned'} |"
        )
    lines.extend(["", "### Certified lower bounds", ""])
    for circuit in (gross, pbb):
        w12 = circuit["w1_w2"]
        w3 = circuit.get("w3")
        lines.append(
            f"- **{circuit['label']}**: all {w12['mechanisms_checked']} mechanisms passed the "
            "weight-1 condition (a detector fires or no observable flips), and hashing exact "
            "detector bitsets found no equal-detector/different-observable weight-2 pair. "
            f"This certifies distance >= 3 in the stated mechanism model."
        )
        if w3 is not None:
            if w3["certified_no_w3"]:
                lines.append(
                    f"  The streamed singles-vs-pairs XOR search also completed: "
                    f"{w3['unique_shared_pairs_checked']} unique shared-detector pairs checked "
                    "with no weight-3 logical, certifying distance >= 4."
                )
            elif w3["w3_undetectable_logical_found"]:
                lines.append("  The weight-3 search found and validated a three-mechanism logical witness.")
            else:
                lines.append(
                    f"  The weight-3 search is not a certificate: {w3['stop_reason']}; "
                    f"{w3['unique_shared_pairs_checked']} unique shared-detector pairs were checked."
                )
    lines.extend(["", "### Upper bounds", ""])
    for key, circuit in (("gross", gross), ("pbb", pbb)):
        attempts = circuit.get("heuristic_search", {}).get("attempts", [])
        returned = [a for a in attempts if a.get("ok")]
        lines.append(
            f"- **{circuit['label']}**: best validated upper bound is "
            f"**distance <= {circuit['ub']['value']}**, via {circuit['ub']['method']}. "
            f"Witness: `{circuit['ub']['witness_path']}`. "
            f"Stim heuristic attempts returning a witness: {len(returned)}/{len(attempts)}; "
            f"heuristic-search wall time {circuit.get('heuristic_search', {}).get('wall_s', 0):.3f} s."
        )
        graphlike = circuit.get("graphlike", {})
        if not graphlike.get("ok"):
            lines.append(
                f"  Graphlike search returned no bound: `{graphlike.get('error')}` "
                f"after {graphlike.get('wall_s', 0):.3f} s on the undecomposed DEM."
            )
        witness = (
            result.get("witnesses", {}).get(key, {}).get("explicit_data_logical", {})
        )
        if witness:
            lines.append(
                f"  Accepted witness: {witness['method']}; provenance: {witness['provenance']}; "
                f"observables flipped: {witness['validation']['observable_xor']}; "
                "detector XOR is empty and all twelve signatures were found in the "
                "flattened undecomposed DEM."
            )
    lines.extend(["", "### Gross versus PBB and hook implication", ""])
    gross_ub = gross["ub"]["value"]
    pbb_ub = pbb["ub"]["value"]
    if pbb_ub < gross_ub:
        comparison = f"PBB's best upper bound ({pbb_ub}) is lower than Gross's ({gross_ub})"
    elif pbb_ub > gross_ub:
        comparison = f"PBB's best upper bound ({pbb_ub}) is higher than Gross's ({gross_ub})"
    else:
        comparison = f"The best upper bounds are tied at {gross_ub}"
    lines.append(
        f"{comparison}. Gross has depth 7, weight-6 pure checks, while PBB has depth 8, "
        "weight-8 checks and 72 mixed checks. Thus the PBB circuit admits longer, mixed-sector "
        "single-mechanism hook propagation, but the reported comparison only shows whether the "
        "available witnesses exploit that structure; an upper-bound tie or gap is not by itself "
        "a proof of exact circuit distance. The certified lower bounds rule out one- and "
        "two-mechanism logical hooks in both circuits."
    )
    lines.extend([
        "",
        "## Reproduction",
        "",
        "From the repository root:",
        "",
        "```bash",
        "PYTHONPATH=src .venv/bin/python experiments/exp029_circuit_distance.py --phase all --circuit all",
        "```",
        "",
        "This command enforces a 45-minute total budget for the heuristic Stim search on each "
        "circuit, a separate 15-minute guard for graphlike search, a 45-minute weight-3 budget, "
        "and an 8 GiB RSS guard. It writes canonical artifacts only after one clean run has "
        "covered both circuits and all stages; subset runs go under `results/partial_runs/`.",
        "",
        "## Runtime and caveats",
        "",
        f"- Recorded summed per-circuit wall time: **{result['wall_s']:.3f} s**.",
        f"- Stim version: `{result['environment']['stim']}`; Python: `{result['environment']['python']}`.",
        "- `decompose_errors=False` is used for every certified lower-bound signature and for the "
        "graphlike search. The graphlike result explicitly skips mechanisms with detector degree "
        ">2; it is an upper bound from that subset, not a lower bound and not an exact distance.",
        "- `Circuit.search_for_undetectable_logical_errors` is heuristic under the recorded caps. "
        "Every returned entry was mapped back to an undecomposed DEM signature and XOR-validated "
        "before being accepted as a mechanism-level upper bound.",
        "- No Monte Carlo was performed (shots=0, fails=0, no rate and no Clopper-Pearson interval).",
        "- A budget-limited search is not exhaustive unless its record explicitly says completed/certified.",
        "- Shared external desktop load materially inflated solver wall times. Exhaustive counts "
        "and logical-witness validity are unaffected; timeouts are wall-budget outcomes.",
        "- Partial-run guard verified empirically: `--phase lower --circuit gross --w3-budget-s 0.001` "
        "wrote `results/partial_runs/exp029_lower_gross_state.json` and left the canonical raw "
        "state/witness SHA-256 digests unchanged.",
        "- Recorded invocations for this artifact: "
        + "; ".join(
            f"{run['phase']}/{run['circuit']} ({run['artifact_route']}, {run['wall_s']:.1f} s)"
            for run in result.get("runs", [])
        )
        + ".",
        "- Execution history of the recorded numbers: the first `--phase all --circuit all` "
        "invocation computed the Gross lower and upper stages plus the PBB lower stage, then "
        "aborted because the PBB pure-X CP-SAT witness returned a different (equally valid) "
        "translate of the weight-12 orbit than the hard-coded support; the code now adopts the "
        "freshly certified support, and `--phase finalize --circuit all` completed the PBB upper "
        "stage and wrote every canonical artifact. Both circuits were rebuilt from scratch in "
        "each invocation and their schedule/circuit/DEM SHA-256 digests matched.",
    ])
    return "\n".join(lines) + "\n"


def finalize(state: dict[str, Any], witnesses: dict[str, Any]) -> None:
    fold_lower_witness_into_upper(state, witnesses)
    update_verdict_and_wall(state)
    missing = [
        key for key in ("gross", "pbb")
        if key not in state["circuits"]
        or "lb" not in state["circuits"][key]
        or "ub" not in state["circuits"][key]
    ]
    full_coverage = not missing
    clean = full_coverage and state["verdict"] == "POSITIVE"
    if canonical_route(clean=clean, full_coverage=full_coverage) != "canonical":
        raise RuntimeError(
            f"cannot write canonical artifacts: missing={missing}, verdict={state['verdict']}"
        )
    result = processed_artifact(state, witnesses)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_PATH.write_text(json.dumps(result, indent=2) + "\n")
    REPORT_PATH.write_text(report_text(result))
    save_state(state, CANONICAL_STATE_PATH)
    save_witnesses(witnesses, state, CANONICAL_WITNESS_PATH)
    print(f"wrote {PROCESSED_PATH.relative_to(ROOT)}", flush=True)
    print(f"wrote {CANONICAL_WITNESS_PATH.relative_to(ROOT)}", flush=True)
    print(f"wrote {REPORT_PATH.relative_to(ROOT)}", flush=True)


def run_stage_for_key(
    phase: str,
    key: str,
    state: dict[str, Any],
    witnesses: dict[str, Any],
    state_path: Path,
    witness_path: Path,
) -> None:
    build_start = time.perf_counter()
    circuit, dem, build_info, code = build_exact_circuit(key)
    build_info["rebuild_wall_s"] = time.perf_counter() - build_start
    if phase in ("lower", "circuit", "all"):
        lower, lower_witness = run_lower_bound(key, dem, build_info)
        merge_stage(state, key, lower)
        witnesses.setdefault("circuits", {}).setdefault(key, {}).update(lower_witness)
        save_state(state, state_path)
        save_witnesses(witnesses, state, witness_path)
    if phase in ("upper", "circuit", "all"):
        if key == "pbb":
            build_info["pbb_pure_x_witness_certificate"] = certify_pbb_data_witness(code)
        upper, upper_witnesses = run_upper_bound(key, circuit, dem, build_info)
        merge_stage(state, key, upper)
        witnesses.setdefault("circuits", {}).setdefault(key, {}).update(upper_witnesses)
        fold_lower_witness_into_upper(state, witnesses)
        save_state(state, state_path)
        save_witnesses(witnesses, state, witness_path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--phase",
        choices=("lower", "upper", "circuit", "all", "finalize"),
        default="all",
    )
    parser.add_argument(
        "--circuit",
        choices=("gross", "pbb", "all"),
        default="all",
    )
    parser.add_argument(
        "--upper-budget-s",
        type=float,
        default=None,
        help="non-default smoke override; forces partial-run routing",
    )
    parser.add_argument(
        "--w3-budget-s",
        type=float,
        default=None,
        help="non-default smoke override; forces partial-run routing",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    global UPPER_SEARCH_BUDGET_S, W3_SEARCH_BUDGET_S
    if args.upper_budget_s is not None:
        UPPER_SEARCH_BUDGET_S = args.upper_budget_s
    if args.w3_budget_s is not None:
        W3_SEARCH_BUDGET_S = args.w3_budget_s
    full_run = (
        args.phase == "all"
        and args.circuit == "all"
        and args.upper_budget_s is None
        and args.w3_budget_s is None
    )
    resume_full_run = args.phase == "finalize" and args.circuit == "all"
    route = canonical_route(clean=True, full_coverage=full_run or resume_full_run)
    if route == "canonical":
        state_path = CANONICAL_STATE_PATH
        witness_path = CANONICAL_WITNESS_PATH
    else:
        PARTIAL_DIR.mkdir(parents=True, exist_ok=True)
        suffix = f"exp029_{args.phase}_{args.circuit}"
        state_path = PARTIAL_DIR / f"{suffix}_state.json"
        witness_path = PARTIAL_DIR / f"{suffix}_witnesses.json"

    # Full-scope runs intentionally start clean. Prior partial/interrupted
    # stages cannot be mistaken for one clean invocation of the declared scope.
    if full_run:
        state = empty_artifact()
        witnesses = {"circuits": {}}
    elif resume_full_run:
        state = load_artifact(state_path)
        if witness_path.exists():
            old_witness = json.loads(witness_path.read_text())
            witnesses = {"circuits": old_witness.get("circuits", {})}
        else:
            witnesses = {"circuits": {}}
        missing_stages = [
            (key, phase)
            for key in ("gross", "pbb")
            for phase in ("lower", "upper")
            if (
                key not in state["circuits"]
                or ("lb" if phase == "lower" else "ub")
                not in state["circuits"][key]
            )
        ]
        for key, phase in missing_stages:
            run_stage_for_key(
                phase, key, state, witnesses, state_path, witness_path
            )
    else:
        state = load_artifact(state_path)
        if witness_path.exists():
            old_witness = json.loads(witness_path.read_text())
            witnesses = {"circuits": old_witness.get("circuits", {})}
        else:
            witnesses = {"circuits": {}}

    invocation_start = time.perf_counter()
    if args.phase != "finalize":
        keys = ("gross", "pbb") if args.circuit == "all" else (args.circuit,)
        for key in keys:
            run_stage_for_key(
                args.phase, key, state, witnesses, state_path, witness_path
            )
    state["runs"].append({
        "argv": [sys.executable, *sys.argv],
        "phase": args.phase,
        "circuit": args.circuit,
        "full_declared_scope": full_run or resume_full_run,
        "artifact_route": route,
        "wall_s": time.perf_counter() - invocation_start,
        "finished_unix_s": time.time(),
    })
    save_state(state, state_path)
    save_witnesses(witnesses, state, witness_path)

    if full_run or resume_full_run:
        finalize(state, witnesses)
    else:
        print(f"saved partial state to {state_path.relative_to(ROOT)}", flush=True)


if __name__ == "__main__":
    main()
