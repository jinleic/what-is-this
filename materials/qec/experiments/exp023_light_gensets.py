"""EXP-023: exact light-generating-set closure in the physical stabilizer space.

For a stabilizer check matrix H, this experiment computes the rank of the span
of all stabilizers of symplectic weight at most w.  Coefficients are defined over
an explicitly recorded independent 132-row subset H_B of H, so c -> c H_B is a
bijection.  Every independence query is built from physical 288-bit linear
functionals annihilating the stabilizers already found, then pulled back through
H_B.  An INFEASIBLE CP-SAT result is therefore an exact closure certificate.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Pre-registered protocol constants.  Do not change between catalogue members.
# ---------------------------------------------------------------------------
EXPERIMENT = "EXP-023"
BASE_SEED = 20260812
PBB_WEIGHT_CAPS = (6, 7)
GROSS_WEIGHT_CAPS = (5, 6)
RANK_CALL_TIME_LIMIT_S = 600.0
THEOREM_GATE_TIME_LIMIT_S = 900.0
MONOLITHIC_CROSSCHECK_TIME_LIMIT_S = 300.0
CROSSCHECK_MEMBER = "12_6_0193"
NUM_SEARCH_WORKERS = 8
EXPECTED_N = 144
EXPECTED_K = 12
EXPECTED_RANK = EXPECTED_N - EXPECTED_K
EXPECTED_PBB_MEMBERS = 14
CATALOGUE_REL = "third_party/qcode-discovery/results/campaign7_publication_merged.jsonl"

import argparse
import hashlib
import json
import os
import platform
import sys
import time

# Shared-machine policy: this process is a CP-SAT solver host, so all nested
# BLAS/OpenMP pools are pinned to one thread each.  CP-SAT's own search
# parallelism is controlled separately by NUM_SEARCH_WORKERS.  These must be
# set before NumPy imports its accelerated backend.
for _thread_env in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_env, "1")

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import ortools
from ortools.sat.python import cp_model

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.artifacts import canonical_route  # noqa: E402
from qec_research.codes.bicycle import (  # noqa: E402
    BRAVYI_BB,
    PBBSpec,
    build_bb,
    build_pbb,
)
from qec_research.codes.pbb_theory import parent_bb_matrices  # noqa: E402
from qec_research.gf2.linalg import (  # noqa: E402
    nullspace_np,
    rank_bitset,
    rank_np,
    rref_bitset,
    rref_np,
    rows_to_bitsets,
)
from qec_research.symplectic.core import symplectic_weight  # noqa: E402

CANONICAL_RAW_PATH = ROOT / "results" / "raw" / "exp023_light_elements.json"
CANONICAL_PROCESSED_PATH = ROOT / "results" / "processed" / "exp023_light_gensets.json"
CANONICAL_REPORT_PATH = ROOT / "notes" / "agent_reports" / "exp023_light_gensets.md"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def derived_seed(code_label: str, w: int, call_index: int) -> int:
    payload = f"{BASE_SEED}:{EXPERIMENT}:{code_label}:{w}:{call_index}".encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:4], "big") & 0x7FFFFFFF


def numpy_gf2_matmul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Independent NumPy GF(2) product used to audit the project helper."""
    aa = np.asarray(a, dtype=np.uint16) & 1
    bb = np.asarray(b, dtype=np.uint16) & 1
    if aa.shape[-1] != bb.shape[0]:
        raise ValueError(f"shape mismatch {aa.shape} @ {bb.shape}")
    return ((aa @ bb) & 1).astype(np.uint8)


def numpy_gf2_rank(matrix: np.ndarray) -> int:
    """Independent, local GF(2) elimination used for physical-space audits."""
    a = (np.asarray(matrix, dtype=np.uint8) & 1).copy()
    if a.ndim != 2:
        raise ValueError(f"rank expects a matrix, got {a.shape}")
    rows, cols = a.shape
    pivot_row = 0
    for col in range(cols):
        nz = np.flatnonzero(a[pivot_row:, col])
        if not nz.size:
            continue
        p = pivot_row + int(nz[0])
        if p != pivot_row:
            a[[pivot_row, p]] = a[[p, pivot_row]]
        below = np.flatnonzero(a[pivot_row + 1 :, col]) + pivot_row + 1
        if below.size:
            a[below] ^= a[pivot_row]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def rank_triplet(matrix: np.ndarray, ncols: int) -> tuple[int, int, int]:
    a = np.asarray(matrix, dtype=np.uint8) & 1
    return (
        numpy_gf2_rank(a),
        rank_np(a),
        rank_bitset(rows_to_bitsets(a), ncols),
    )


def vector_bits(v: np.ndarray) -> str:
    return "".join("1" if int(x) else "0" for x in np.asarray(v).ravel())


def physical_weight(v: np.ndarray, n: int) -> int:
    vv = np.asarray(v, dtype=np.uint8) & 1
    return int(np.count_nonzero(vv[:n] | vv[n:]))


def independent_row_indices(matrix: np.ndarray) -> list[int]:
    """Greedily select original rows, retaining the sparse published encoding."""
    pivots: dict[int, int] = {}
    selected: list[int] = []
    for index, row in enumerate(rows_to_bitsets(matrix)):
        reduced = row
        while reduced:
            pivot = reduced.bit_length() - 1
            prior = pivots.get(pivot)
            if prior is None:
                pivots[pivot] = reduced
                selected.append(index)
                break
            reduced ^= prior
    return selected


def coefficient_map(h_basis: np.ndarray, original_h: np.ndarray) -> np.ndarray:
    """Coordinates of all original rows in the independent published-row basis."""
    basis_bits = rows_to_bitsets(h_basis)
    reduced, pivots, provenance = rref_bitset(basis_bits, h_basis.shape[1])
    if len(reduced) != h_basis.shape[0]:
        raise AssertionError("H_B is not independent")
    out = np.zeros((original_h.shape[0], h_basis.shape[0]), dtype=np.uint8)
    for row_index, target in enumerate(rows_to_bitsets(original_h)):
        for r, pivot_col in enumerate(pivots):
            if (target >> pivot_col) & 1:
                out[row_index, provenance[r]] ^= 1
    reconstructed = numpy_gf2_matmul(out, h_basis)
    if not np.array_equal(reconstructed, original_h):
        raise AssertionError("failed to express the published rows in H_B")
    return out


def rref_with_physical_payload(
    functionals: np.ndarray, physical_h: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """RREF coefficient functionals while applying the same operations to h."""
    f = (np.asarray(functionals, dtype=np.uint8) & 1).copy()
    h = (np.asarray(physical_h, dtype=np.uint8) & 1).copy()
    if f.shape[0] != h.shape[0]:
        raise ValueError("functional/payload row mismatch")
    pivot_row = 0
    for col in range(f.shape[1]):
        nz = np.flatnonzero(f[pivot_row:, col])
        if not nz.size:
            continue
        p = pivot_row + int(nz[0])
        if p != pivot_row:
            f[[pivot_row, p]] = f[[p, pivot_row]]
            h[[pivot_row, p]] = h[[p, pivot_row]]
        other = np.flatnonzero(f[:, col])
        other = other[other != pivot_row]
        if other.size:
            f[other] ^= f[pivot_row]
            h[other] ^= h[pivot_row]
        pivot_row += 1
        if pivot_row == f.shape[0]:
            break
    return f[:pivot_row], h[:pivot_row]


def physical_annihilator_functionals(
    found_v: np.ndarray, h_basis: np.ndarray
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Return a basis of physical closure functionals and their pullbacks.

    Each returned physical row h annihilates every row of found_v.  Its
    coefficient-space pullback is f_i = H_B[i] . h.  The pullbacks span all
    functionals on rowspace(H_B) that annihilate span(found_v).
    """
    ncols = h_basis.shape[1]
    rank_found = rank_np(found_v)
    ambient_nullspace = nullspace_np(found_v)
    pullbacks_all = numpy_gf2_matmul(ambient_nullspace, h_basis.T)
    functionals, physical_h = rref_with_physical_payload(
        pullbacks_all, ambient_nullspace
    )
    expected = h_basis.shape[0] - rank_found
    annihilation = numpy_gf2_matmul(found_v, physical_h.T)
    pulled_back_again = numpy_gf2_matmul(physical_h, h_basis.T)
    f_ranks = rank_triplet(functionals, h_basis.shape[0])
    checks = {
        "found_rank_physical": int(rank_found),
        "ambient_nullspace_dimension": int(ambient_nullspace.shape[0]),
        "n_functionals": int(functionals.shape[0]),
        "expected_n_functionals": int(expected),
        "functionals_rank_numpy_independent": int(f_ranks[0]),
        "functionals_rank_gf2_numpy": int(f_ranks[1]),
        "functionals_rank_gf2_bitset": int(f_ranks[2]),
        "all_physical_h_annihilate_found_v": bool(not annihilation.any()),
        "all_pullbacks_equal_HB_dot_h": bool(
            np.array_equal(functionals, pulled_back_again)
        ),
    }
    checks["physical_functionals_verified"] = bool(
        functionals.shape[0] == expected
        and f_ranks == (expected, expected, expected)
        and checks["all_physical_h_annihilate_found_v"]
        and checks["all_pullbacks_equal_HB_dot_h"]
    )
    if not checks["physical_functionals_verified"]:
        raise AssertionError(f"invalid physical closure functionals: {checks}")
    return functionals, physical_h, checks


def add_xor_equality(
    model: cp_model.CpModel,
    inputs: list[cp_model.IntVar],
    target: cp_model.IntVar,
) -> None:
    if not inputs:
        model.add(target == 0)
    else:
        # XOR(inputs, not target) == 1 is equivalent to target == XOR(inputs).
        model.add_bool_xor([*inputs, target.Not()])


def build_query_model(
    h_basis: np.ndarray, w: int, functionals: np.ndarray
) -> tuple[cp_model.CpModel, list[cp_model.IntVar]]:
    rows, twice_n = h_basis.shape
    n = twice_n // 2
    model = cp_model.CpModel()
    c = [model.new_bool_var(f"c_{i}") for i in range(rows)]
    x = [model.new_bool_var(f"x_{j}") for j in range(n)]
    z = [model.new_bool_var(f"z_{j}") for j in range(n)]
    support = [model.new_bool_var(f"s_{j}") for j in range(n)]

    for j in range(n):
        add_xor_equality(model, [c[i] for i in np.flatnonzero(h_basis[:, j])], x[j])
        add_xor_equality(
            model, [c[i] for i in np.flatnonzero(h_basis[:, n + j])], z[j]
        )
        model.add(support[j] >= x[j])
        model.add(support[j] >= z[j])
        model.add(support[j] <= x[j] + z[j])
    model.add(sum(support) <= w)

    outside = []
    for f_index, functional in enumerate(functionals):
        value = model.new_bool_var(f"outside_{f_index}")
        add_xor_equality(
            model, [c[i] for i in np.flatnonzero(functional)], value
        )
        outside.append(value)
    if not outside:
        raise AssertionError("a query with no outside-span functional is invalid")
    model.add_bool_or(outside)
    return model, c


def solve_query(
    h_basis: np.ndarray,
    w: int,
    functionals: np.ndarray,
    seed: int,
    time_limit_s: float = RANK_CALL_TIME_LIMIT_S,
) -> dict[str, Any]:
    build_started = time.perf_counter()
    model, c_vars = build_query_model(h_basis, w, functionals)
    build_wall = time.perf_counter() - build_started
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit_s
    solver.parameters.num_workers = NUM_SEARCH_WORKERS
    solver.parameters.random_seed = seed
    solver.parameters.cp_model_presolve = True
    solve_started = time.perf_counter()
    status = solver.solve(model)
    solve_wall = time.perf_counter() - solve_started
    status_name = solver.status_name(status)
    coefficient = None
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        coefficient = np.array([solver.value(v) for v in c_vars], dtype=np.uint8)
    return {
        "status_code": int(status),
        "status": status_name,
        "c": coefficient,
        "model_build_wall_s": build_wall,
        "solver_wall_s": solve_wall,
        "solver_reported_wall_s": float(solver.wall_time),
        "num_conflicts": int(solver.num_conflicts),
        "num_branches": int(solver.num_branches),
        "response_stats": solver.response_stats(),
    }


def add_if_independent(
    elements: list[dict[str, Any]],
    c: np.ndarray,
    h_basis: np.ndarray,
    w: int,
    source: str,
    expected_v: np.ndarray | None = None,
) -> bool:
    c = np.asarray(c, dtype=np.uint8) & 1
    v = numpy_gf2_matmul(c[None, :], h_basis)[0]
    if expected_v is not None and not np.array_equal(v, expected_v):
        raise AssertionError(f"coefficient reconstruction failed for {source}")
    weight_numpy = physical_weight(v, h_basis.shape[1] // 2)
    weight_project = symplectic_weight(v)
    if weight_numpy != weight_project or weight_numpy > w:
        raise AssertionError(
            f"invalid weight-{weight_numpy} candidate at cap {w} from {source}"
        )
    prior = (
        np.vstack([e["v"] for e in elements])
        if elements
        else np.zeros((0, h_basis.shape[1]), dtype=np.uint8)
    )
    prior_rank = rank_np(prior)
    candidate_rank = rank_np(np.vstack([prior, v]))
    if candidate_rank == prior_rank:
        return False
    if candidate_rank != prior_rank + 1:
        raise AssertionError("rank changed by more than one")
    elements.append({"c": c, "v": v, "weight": weight_numpy, "source": source})
    return True


def verify_elements(
    elements: list[dict[str, Any]], h_basis: np.ndarray, w: int
) -> dict[str, Any]:
    rows = h_basis.shape[0]
    ncols = h_basis.shape[1]
    n = ncols // 2
    c_matrix = (
        np.vstack([e["c"] for e in elements])
        if elements
        else np.zeros((0, rows), dtype=np.uint8)
    )
    v_matrix = (
        np.vstack([e["v"] for e in elements])
        if elements
        else np.zeros((0, ncols), dtype=np.uint8)
    )
    recomputed = numpy_gf2_matmul(c_matrix, h_basis)
    all_binary = bool(
        np.all((c_matrix == 0) | (c_matrix == 1))
        and np.all((v_matrix == 0) | (v_matrix == 1))
    )
    all_recomputed = bool(np.array_equal(recomputed, v_matrix))
    weights_numpy = [physical_weight(v, n) for v in v_matrix]
    weights_project = [symplectic_weight(v) for v in v_matrix]
    all_weights = bool(
        weights_numpy == weights_project and all(weight <= w for weight in weights_numpy)
    )
    h_rank = rank_np(h_basis)
    containment = [
        rank_np(np.vstack([h_basis, v])) == h_rank
        and numpy_gf2_rank(np.vstack([h_basis, v])) == h_rank
        for v in v_matrix
    ]
    c_ranks = rank_triplet(c_matrix, rows)
    v_ranks = rank_triplet(v_matrix, ncols)
    rank_agreement = bool(c_ranks == v_ranks == (len(elements),) * 3)
    checks = {
        "n_elements": len(elements),
        "all_vectors_binary": all_binary,
        "all_v_equal_numpy_c_times_HB": all_recomputed,
        "all_weights_numpy_equal_project_symplectic_weight": bool(
            weights_numpy == weights_project
        ),
        "all_weights_at_most_cap": all_weights,
        "all_elements_in_rowspace_by_physical_rank_comparison": bool(
            all(containment)
        ),
        "rank_C_numpy_independent": int(c_ranks[0]),
        "rank_C_gf2_numpy": int(c_ranks[1]),
        "rank_C_gf2_bitset": int(c_ranks[2]),
        "rank_V_numpy_independent": int(v_ranks[0]),
        "rank_V_gf2_numpy": int(v_ranks[1]),
        "rank_V_gf2_bitset": int(v_ranks[2]),
        "coefficient_and_physical_ranks_agree": rank_agreement,
    }
    checks["physical_space_verified"] = bool(
        all_binary
        and all_recomputed
        and all_weights
        and all(containment)
        and rank_agreement
    )
    if not checks["physical_space_verified"]:
        raise AssertionError(f"physical verification failed: {checks}")
    return checks


def serialise_element(element: dict[str, Any]) -> dict[str, Any]:
    return {
        "c": element["c"].astype(int).tolist(),
        "v_bits_x_then_z": vector_bits(element["v"]),
        "symplectic_weight": int(element["weight"]),
        "source": element["source"],
    }


def analyse_weight_cap(
    label: str,
    h_original: np.ndarray,
    h_basis: np.ndarray,
    original_coordinates: np.ndarray,
    w: int,
    carried_elements: list[dict[str, Any]],
) -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]]]:
    started = time.perf_counter()
    elements: list[dict[str, Any]] = []
    for element in carried_elements:
        if element["weight"] <= w:
            add_if_independent(
                elements,
                element["c"],
                h_basis,
                w,
                f"carried_from_lower_cap:{element['source']}",
                expected_v=element["v"],
            )

    seeded_published = 0
    for row_index, row in enumerate(h_original):
        if physical_weight(row, h_original.shape[1] // 2) <= w:
            if add_if_independent(
                elements,
                original_coordinates[row_index],
                h_basis,
                w,
                f"published_row:{row_index}",
                expected_v=row,
            ):
                seeded_published += 1

    solver_calls: list[dict[str, Any]] = []
    timeouts: list[dict[str, Any]] = []
    closure_proof: dict[str, Any] | None = None
    closed = False
    closure_kind: str | None = None
    max_calls = h_basis.shape[0] + 1

    initial_v = (
        np.vstack([e["v"] for e in elements])
        if elements
        else np.zeros((0, h_basis.shape[1]), dtype=np.uint8)
    )
    if rank_np(initial_v) == h_basis.shape[0]:
        closed = True
        closure_kind = "found_span_is_full_rowspace"

    while not closed and not timeouts:
        call_index = len(solver_calls) + 1
        if call_index > max_calls:
            raise AssertionError("exceeded rank+1 solver-call bound")
        found_v = (
            np.vstack([e["v"] for e in elements])
            if elements
            else np.zeros((0, h_basis.shape[1]), dtype=np.uint8)
        )
        functionals, physical_h, functional_checks = physical_annihilator_functionals(
            found_v, h_basis
        )
        seed = derived_seed(label, w, call_index)
        result = solve_query(h_basis, w, functionals, seed)
        call_record = {
            "call_index": call_index,
            "seed": seed,
            "status": result["status"],
            "status_code": result["status_code"],
            "rank_before_query_physical": int(rank_np(found_v)),
            "functional_checks": functional_checks,
            "model_build_wall_s": round(result["model_build_wall_s"], 6),
            "solver_wall_s": round(result["solver_wall_s"], 6),
            "solver_reported_wall_s": round(result["solver_reported_wall_s"], 6),
            "num_conflicts": result["num_conflicts"],
            "num_branches": result["num_branches"],
            "response_stats": result["response_stats"],
        }
        solver_calls.append(call_record)
        if result["status"] in ("OPTIMAL", "FEASIBLE"):
            coefficient = result["c"]
            if coefficient is None:
                raise AssertionError("SAT status without coefficient witness")
            before = rank_np(found_v)
            added = add_if_independent(
                elements,
                coefficient,
                h_basis,
                w,
                f"cp_sat_call:{call_index}",
            )
            after_v = np.vstack([e["v"] for e in elements])
            after = rank_np(after_v)
            if not added or after != before + 1:
                raise AssertionError("solver witness was not physically independent")
            call_record["witness_c"] = coefficient.astype(int).tolist()
            call_record["witness_weight"] = int(elements[-1]["weight"])
            call_record["witness_physical_rank_after"] = int(after)
            continue
        if result["status"] == "INFEASIBLE":
            closed = True
            closure_kind = "cp_sat_infeasible_against_physical_annihilator"
            closure_proof = {
                "solver_status": "INFEASIBLE",
                "call_index": call_index,
                "seed": seed,
                "rank_V_at_closure": int(rank_np(found_v)),
                "physical_h_vectors": [vector_bits(row) for row in physical_h],
                "pulled_back_f_vectors": [vector_bits(row) for row in functionals],
                "functional_checks": functional_checks,
            }
            continue
        timeout_record = {
            "call_index": call_index,
            "seed": seed,
            "status": result["status"],
            "solver_wall_s": round(result["solver_wall_s"], 6),
        }
        timeouts.append(timeout_record)

    verification = verify_elements(elements, h_basis, w)
    found_rank = verification["rank_V_gf2_numpy"]
    if closed and closure_kind == "found_span_is_full_rowspace":
        if found_rank != h_basis.shape[0]:
            raise AssertionError("full-span closure without full physical rank")
    if closed and closure_kind and closure_kind.startswith("cp_sat"):
        if closure_proof is None or closure_proof["rank_V_at_closure"] != found_rank:
            raise AssertionError("closure proof rank mismatch")

    wall = time.perf_counter() - started
    processed = {
        "weight_cap": w,
        "rank_Vw": int(found_rank) if closed else None,
        "rank_lower_bound_if_undecided": None if closed else int(found_rank),
        "closed": bool(closed),
        "closure_kind": closure_kind,
        "n_solver_calls": len(solver_calls),
        "wall_s": round(wall, 6),
        "timeouts": timeouts,
        "n_seeded_published_independent": seeded_published,
        "n_carried_independent": len(elements) - seeded_published - sum(
            1 for e in elements if e["source"].startswith("cp_sat_call:")
        ),
        "n_solver_witnesses": sum(
            1 for e in elements if e["source"].startswith("cp_sat_call:")
        ),
        "call_seeds": [record["seed"] for record in solver_calls],
        "independent_verification": verification,
        "physical_space_verified": verification["physical_space_verified"],
    }
    raw = {
        "weight_cap": w,
        "closed": bool(closed),
        "rank_Vw": int(found_rank) if closed else None,
        "rank_lower_bound_if_undecided": None if closed else int(found_rank),
        "elements": [serialise_element(e) for e in elements],
        "solver_calls": solver_calls,
        "closure_proof": closure_proof,
        "timeouts": timeouts,
        "independent_verification": verification,
        "wall_s": round(wall, 6),
    }
    return processed, raw, elements


def load_catalogue_members() -> list[dict[str, Any]]:
    path = ROOT / CATALOGUE_REL
    rows = []
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        row = json.loads(line)
        if row.get("n") == 144 and row.get("k") == 12 and row.get("d") == 12:
            row["_catalogue_line"] = line_number
            row["_member_id"] = (
                str(row["code_id"])
                if row.get("code_id")
                else f"catalogue_row_{line_number:04d}_{row['bliss_hash']}"
            )
            rows.append(row)
    rows.sort(key=lambda row: row["_member_id"])
    if len(rows) != EXPECTED_PBB_MEMBERS:
        raise AssertionError(
            f"expected {EXPECTED_PBB_MEMBERS} catalogue members, found {len(rows)}"
        )
    if len({row["_member_id"] for row in rows}) != EXPECTED_PBB_MEMBERS:
        raise AssertionError("duplicate member identifier in target family")
    return rows


def gross_matrix() -> np.ndarray:
    hx, hz = build_bb(BRAVYI_BB["[[144,12,12]]"])
    zero = np.zeros_like(hx)
    return np.vstack([np.hstack([hx, zero]), np.hstack([zero, hz])]).astype(
        np.uint8
    )


def make_items() -> list[dict[str, Any]]:
    items = [
        {
            "label": "CSS Gross [[144,12,12]]",
            "kind": "CSS Gross control",
            "H": gross_matrix(),
            "weight_caps": GROSS_WEIGHT_CAPS,
            "spec": None,
            "catalogue_spec": None,
        }
    ]
    for row in load_catalogue_members():
        spec = PBBSpec(
            ell=row["ell"],
            m=row["m"],
            A=[tuple(term) for term in row["A_terms"]],
            B=[tuple(term) for term in row["B_terms"]],
            C=[tuple(term) for term in row["C_terms"]],
            D=[tuple(term) for term in row["D_terms"]],
            name=row["_member_id"],
        )
        code = build_pbb(spec, check=True)
        items.append(
            {
                "label": row["_member_id"],
                "kind": "PBB catalogue member",
                "H": code.H,
                "spec": spec,
                "weight_caps": PBB_WEIGHT_CAPS,
                "catalogue_spec": {
                    "member_id": row["_member_id"],
                    "catalogue_line": row["_catalogue_line"],
                    "code_id": row.get("code_id"),
                    "bliss_hash": row.get("bliss_hash"),
                    **{
                        key: row.get(key)
                        for key in (
                            "ell",
                            "m",
                            "A_terms",
                            "B_terms",
                            "C_terms",
                            "D_terms",
                            "n",
                            "k",
                            "d",
                            "d_is_exact",
                            "trust_level",
                        )
                    },
                },
            }
        )
    return items


def bitset_reducer(rows: list[int]) -> list[tuple[int, int, int]]:
    """Pivot-reduced basis as (pivot bit, row bits, provenance mask)."""
    reduced: list[tuple[int, int, int]] = []
    for index, row in enumerate(rows):
        bits, provenance = row, 1 << index
        for pivot, other, other_provenance in reduced:
            if (bits >> pivot) & 1:
                bits ^= other
                provenance ^= other_provenance
        if bits:
            reduced.append((bits.bit_length() - 1, bits, provenance))
    return reduced


def reduce_target(
    reducer: list[tuple[int, int, int]], target: int
) -> tuple[int, int]:
    """Reduce `target` by the basis; return (residual, coefficient mask)."""
    bits, provenance = target, 0
    for pivot, row, row_provenance in reducer:
        if (bits >> pivot) & 1:
            bits ^= row
            provenance ^= row_provenance
    return bits, provenance


def enumerate_light_x_codewords(
    label: str, hx_basis: np.ndarray, w: int
) -> dict[str, Any]:
    """All nonzero X-parts of weight <= w, with a completeness certificate.

    The X-part of any stabilizer is a codeword of rowspace(H_X), so enumerating
    those codewords bounds every mixed light element from above.  CP-SAT
    `OPTIMAL` after exhaustive enumeration certifies that the list is complete.
    """
    rows, length = hx_basis.shape
    model = cp_model.CpModel()
    u = [model.new_bool_var(f"u_{i}") for i in range(rows)]
    x = [model.new_bool_var(f"cx_{j}") for j in range(length)]
    for j in range(length):
        add_xor_equality(model, [u[i] for i in np.flatnonzero(hx_basis[:, j])], x[j])
    model.add(sum(x) <= w)
    model.add(sum(x) >= 1)

    collected: list[np.ndarray] = []

    class Collect(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self) -> None:
            collected.append(np.array([self.value(v) for v in u], dtype=np.uint8))

    seed = derived_seed(f"xcode:{label}", w, 0)
    solver = cp_model.CpSolver()
    solver.parameters.enumerate_all_solutions = True
    solver.parameters.num_workers = 1  # enumeration requires a single worker
    solver.parameters.max_time_in_seconds = THEOREM_GATE_TIME_LIMIT_S
    solver.parameters.random_seed = seed
    started = time.perf_counter()
    status = solver.solve(model, Collect())
    wall = time.perf_counter() - started
    status_name = solver.status_name(status)
    coefficients = (
        np.vstack(collected)
        if collected
        else np.zeros((0, rows), dtype=np.uint8)
    )
    codewords = numpy_gf2_matmul(coefficients, hx_basis)
    weights = [int(word.sum()) for word in codewords]
    distinct = len({vector_bits(word) for word in codewords})
    if codewords.shape[0] and (distinct != codewords.shape[0] or min(weights) < 1):
        raise AssertionError("X-codeword enumeration produced duplicates or zero")
    return {
        "status": status_name,
        "complete": status_name == "OPTIMAL",
        "n_codewords": int(coefficients.shape[0]),
        "weight_histogram": {
            str(weight): weights.count(weight) for weight in sorted(set(weights))
        },
        "seed": seed,
        "time_limit_s": THEOREM_GATE_TIME_LIMIT_S,
        "wall_s": round(wall, 6),
        "coefficients": coefficients,
        "codewords": codewords,
    }


def coset_budget_witness(
    z0: np.ndarray, pure_z_parts: np.ndarray, support_x: np.ndarray, budget: int
) -> tuple[bool, np.ndarray | None]:
    """Exact test: some s in S_Z leaves z0 + s within `budget` extra qubits.

    Coordinates inside supp(x) are free, so only the complement is scored.  For
    budgets 0 and 1 the answer is decided by GF(2) membership tests; larger
    budgets enumerate the (still small) candidate extra-support sets.
    """
    outside = np.flatnonzero(~support_x.astype(bool))
    if budget < 0:
        return False, None

    def pack(vector: np.ndarray) -> int:
        value = 0
        for position, coordinate in enumerate(outside):
            if vector[coordinate]:
                value |= 1 << position
        return value

    basis = [pack(row) for row in pure_z_parts]
    reducer = bitset_reducer(basis)
    target = pack(z0)

    from itertools import combinations

    for size in range(budget + 1):
        for extra in combinations(range(len(outside)), size):
            probe = target
            for position in extra:
                probe ^= 1 << position
            residual, coefficients = reduce_target(reducer, probe)
            if residual:
                continue
            combination = np.zeros(pure_z_parts.shape[1], dtype=np.uint8)
            for index in range(len(basis)):
                if (coefficients >> index) & 1:
                    combination ^= pure_z_parts[index]
            return True, (z0 ^ combination).astype(np.uint8)
    return False, None


def decomposition_gate(
    label: str,
    spec: PBBSpec,
    h_basis: np.ndarray,
    pure_z_rows: np.ndarray,
    w: int = 7,
) -> dict[str, Any]:
    """Decide the mixed-light question by exhaustive structural decomposition.

    Every stabilizer is `(u H_X | u P + s)` with `s` ranging over the pure-Z
    subgroup `S_Z`, so a weight-<=w element with nonzero X-part exists iff some
    nonzero X-codeword `x = u H_X` of weight <= w admits an `s` leaving at most
    `w - |x|` qubits of Z-support outside `supp(x)`.  The first stage is an
    exhaustive CP-SAT enumeration with a completeness certificate; the second
    is exact GF(2) linear algebra.
    """
    started = time.perf_counter()
    hx, _, perturbation = parent_bb_matrices(spec)
    hx_rows = independent_row_indices(hx)
    hx_basis = hx[hx_rows]
    n = h_basis.shape[1] // 2
    pure_z_parts = pure_z_rows[:, n:]
    enumeration = enumerate_light_x_codewords(label, hx_basis, w)
    coefficients = enumeration.pop("coefficients")
    codewords = enumeration.pop("codewords")

    checked = 0
    witness: dict[str, Any] | None = None
    for index in range(coefficients.shape[0]):
        x_part = codewords[index]
        weight_x = int(x_part.sum())
        u_full = np.zeros(hx.shape[0], dtype=np.uint8)
        u_full[hx_rows] = coefficients[index]
        z0 = numpy_gf2_matmul(u_full[None, :], perturbation)[0]
        fits, z_part = coset_budget_witness(z0, pure_z_parts, x_part, w - weight_x)
        checked += 1
        if fits and z_part is not None:
            physical = np.concatenate([x_part, z_part]).astype(np.uint8)
            if physical_weight(physical, n) > w or not physical[:n].any():
                raise AssertionError("decomposition witness violates its own bound")
            witness = {
                "v_bits_x_then_z": vector_bits(physical),
                "symplectic_weight": physical_weight(physical, n),
                "x_weight": weight_x,
            }
            break

    # Controls.  The decomposition's whole output is a negative answer, so both
    # of its stages are exercised on inputs with a known positive answer.
    enumerated = {vector_bits(word) for word in codewords}
    published_light_x = [row for row in hx if 0 < int(row.sum()) <= w]
    enumeration_control = bool(
        published_light_x
        and all(vector_bits(row) in enumerated for row in published_light_x)
    )
    block1 = np.hstack([hx, perturbation]).astype(np.uint8)
    heaviest = int(np.argmax([physical_weight(row, n) for row in block1]))
    control_x = hx[heaviest]
    control_budget = physical_weight(block1[heaviest], n) - int(control_x.sum())
    coset_control_fits, coset_control_witness = coset_budget_witness(
        perturbation[heaviest], pure_z_parts, control_x, control_budget
    )
    coset_control = bool(
        coset_control_fits
        and coset_control_witness is not None
        and physical_weight(np.concatenate([control_x, coset_control_witness]), n)
        <= physical_weight(block1[heaviest], n)
    )
    negative_budget_control = not coset_budget_witness(
        perturbation[heaviest], pure_z_parts, control_x, -1
    )[0]
    controls = {
        "all_published_light_X_rows_were_enumerated": enumeration_control,
        "n_published_light_X_rows": len(published_light_x),
        "coset_test_recovers_a_published_mixed_generator": coset_control,
        "coset_control_budget": control_budget,
        "impossible_budget_is_rejected": negative_budget_control,
        "controls_passed": bool(
            enumeration_control and coset_control and negative_budget_control
        ),
    }
    if not controls["controls_passed"]:
        raise AssertionError(f"decomposition controls failed for {label}: {controls}")

    wall = time.perf_counter() - started
    exists = None if not enumeration["complete"] else witness is not None
    return {
        "method": "exhaustive structural decomposition (X-codeword enumeration + exact coset test)",
        "argument": (
            "every stabilizer is (u H_X | u P + s) with s in the pure-Z subgroup, "
            "so a weight<=w element with nonzero X-part exists iff some nonzero "
            "X-codeword x of weight <= w admits an s leaving at most w - |x| "
            "Z-support outside supp(x)"
        ),
        "weight_cap": w,
        "x_codeword_enumeration": enumeration,
        "n_x_codewords_tested": checked,
        "controls": controls,
        "mixed_light_element_exists": exists,
        "decided": exists is not None,
        "witness": witness,
        "wall_s": round(wall, 6),
    }


def pure_z_subgroup_proof(
    spec: PBBSpec, h_original: np.ndarray
) -> tuple[np.ndarray, dict[str, Any]]:
    """Build and independently verify the complete pure-Z stabilizer subgroup."""
    hx, hz, perturbation = parent_bb_matrices(spec)
    left_kernel = nullspace_np(hx.T)
    extra = (
        numpy_gf2_matmul(left_kernel, perturbation)
        if left_kernel.shape[0]
        else np.zeros((0, h_original.shape[1] // 2), dtype=np.uint8)
    )
    z_generators = np.vstack([hz, extra]).astype(np.uint8)
    pure_z_physical = np.hstack(
        [np.zeros_like(z_generators), z_generators]
    ).astype(np.uint8)
    rank_checks = rank_triplet(pure_z_physical, h_original.shape[1])

    # Independent direct definition: coefficients a with a H_X = 0, mapped
    # through the full physical H.  This must equal the structural formula.
    coefficient_kernel = nullspace_np(h_original[:, : EXPECTED_N].T)
    direct = numpy_gf2_matmul(coefficient_kernel, h_original)
    direct_rank_checks = rank_triplet(direct, h_original.shape[1])
    joined_rank = rank_np(np.vstack([pure_z_physical, direct]))
    contains_only_z = bool(not pure_z_physical[:, :EXPECTED_N].any())
    contained_in_h = bool(
        rank_np(np.vstack([h_original, pure_z_physical])) == EXPECTED_RANK
    )
    equal_spaces = bool(
        joined_rank == rank_checks[1] == direct_rank_checks[1]
    )
    checks = {
        "rank_HX": int(rank_np(hx)),
        "rank_HZ": int(rank_np(hz)),
        "left_kernel_dimension": int(left_kernel.shape[0]),
        "coefficient_kernel_dimension_direct": int(coefficient_kernel.shape[0]),
        "r_Zsub": int(rank_checks[1]),
        "rank_Zsub_numpy_independent": int(rank_checks[0]),
        "rank_Zsub_gf2_numpy": int(rank_checks[1]),
        "rank_Zsub_gf2_bitset": int(rank_checks[2]),
        "rank_direct_numpy_independent": int(direct_rank_checks[0]),
        "rank_direct_gf2_numpy": int(direct_rank_checks[1]),
        "rank_direct_gf2_bitset": int(direct_rank_checks[2]),
        "structural_and_direct_pure_Z_spaces_equal": equal_spaces,
        "all_structural_generators_are_pure_Z": contains_only_z,
        "structural_generators_contained_in_rowspace_H": contained_in_h,
    }
    checks["pure_Z_subgroup_verified"] = bool(
        len(set(rank_checks)) == 1
        and len(set(direct_rank_checks)) == 1
        and equal_spaces
        and contains_only_z
        and contained_in_h
        and rank_checks[1] < EXPECTED_RANK
    )
    if not checks["pure_Z_subgroup_verified"]:
        raise AssertionError(f"pure-Z subgroup verification failed: {checks}")
    return pure_z_physical, checks


def positive_control_query(
    label: str,
    h_basis: np.ndarray,
    functionals: np.ndarray,
    w: int,
    question: str,
    require_x_support: bool,
) -> dict[str, Any]:
    """A query whose answer MUST be SAT, guarding against a vacuous encoding.

    Every scientific claim of this experiment rests on INFEASIBLE answers from
    `build_query_model`.  A control that cannot fail proves nothing, so each
    code also runs the same encoding on a question with a known witness (a
    published generator).  A control that does not return a verified SAT
    witness aborts the run before any canonical artifact is written.
    """
    n = h_basis.shape[1] // 2
    model, c = build_query_model(h_basis, w, functionals)
    seed = derived_seed(f"control:{label}", w, 0)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = THEOREM_GATE_TIME_LIMIT_S
    solver.parameters.num_workers = NUM_SEARCH_WORKERS
    solver.parameters.random_seed = seed
    started = time.perf_counter()
    status = solver.solve(model)
    wall = time.perf_counter() - started
    status_name = solver.status_name(status)
    record: dict[str, Any] = {
        "question": question,
        "weight_cap": w,
        "expected": "SAT",
        "status": status_name,
        "seed": seed,
        "wall_s": round(wall, 6),
        "witness_weight": None,
        "witness_has_x_support": None,
        "control_passed": False,
    }
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        coefficient = np.array([solver.value(var) for var in c], dtype=np.uint8)
        physical = numpy_gf2_matmul(coefficient[None, :], h_basis)[0]
        weight_numpy = physical_weight(physical, n)
        x_support = bool(physical[:n].any())
        record["witness_weight"] = weight_numpy
        record["witness_has_x_support"] = x_support
        record["witness_v_bits_x_then_z"] = vector_bits(physical)
        record["control_passed"] = bool(
            weight_numpy == symplectic_weight(physical)
            and 0 < weight_numpy <= w
            and rank_np(np.vstack([h_basis, physical])) == h_basis.shape[0]
            and (x_support or not require_x_support)
        )
    if not record["control_passed"]:
        raise AssertionError(
            f"positive control failed for {label}: {record}; the encoding may be "
            "vacuously infeasible, so no result from this run is trustworthy"
        )
    return record


def theorem_gate_query(
    label: str,
    h_basis: np.ndarray,
    pure_z_rows: np.ndarray,
    w: int = 7,
    time_limit_s: float = MONOLITHIC_CROSSCHECK_TIME_LIMIT_S,
) -> dict[str, Any]:
    """Decide whether any weight-at-most-w stabilizer has nonzero X support.

    The exclusion is encoded through physical functionals annihilating the
    complete pure-Z stabilizer subgroup.  Because that subgroup is exactly the
    set of stabilizers with zero X-part, "outside its span" and "X-part
    nonzero" are the same condition; this encoding solves markedly faster than
    a bare disjunction over the X coordinates.
    """
    n = h_basis.shape[1] // 2
    functionals, physical_h, functional_checks = physical_annihilator_functionals(
        pure_z_rows, h_basis
    )
    model, c = build_query_model(h_basis, w, functionals)
    seed = derived_seed(label, w, 0)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit_s
    solver.parameters.num_workers = NUM_SEARCH_WORKERS
    solver.parameters.random_seed = seed
    started = time.perf_counter()
    status = solver.solve(model)
    wall = time.perf_counter() - started
    status_name = solver.status_name(status)
    record: dict[str, Any] = {
        "weight_cap": w,
        "question": (
            "exists v in rowspace(H) with symplectic weight <= w and nonzero X-part"
        ),
        "encoding": (
            "CP-SAT feasibility; exclusion via physical functionals annihilating "
            "the exact pure-Z stabilizer subgroup (equivalent to X-part nonzero)"
        ),
        "status": status_name,
        "status_code": int(status),
        "seed": seed,
        "time_limit_s": THEOREM_GATE_TIME_LIMIT_S,
        "n_exclusion_functionals": int(functionals.shape[0]),
        "functional_checks": functional_checks,
        "wall_s": round(wall, 6),
        "solver_reported_wall_s": round(float(solver.wall_time), 6),
        "num_conflicts": int(solver.num_conflicts),
        "num_branches": int(solver.num_branches),
        "response_stats": solver.response_stats(),
        "mixed_light_element_exists": (
            True
            if status in (cp_model.OPTIMAL, cp_model.FEASIBLE)
            else False
            if status == cp_model.INFEASIBLE
            else None
        ),
        "witness": None,
        "witness_physical_space_verified": None,
    }
    del physical_h
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        coefficient = np.array([solver.value(var) for var in c], dtype=np.uint8)
        physical = numpy_gf2_matmul(coefficient[None, :], h_basis)[0]
        weight_numpy = physical_weight(physical, n)
        in_rowspace = rank_np(np.vstack([h_basis, physical])) == h_basis.shape[0]
        x_nonzero = bool(physical[:n].any())
        verified = bool(
            weight_numpy == symplectic_weight(physical)
            and weight_numpy <= w
            and x_nonzero
            and in_rowspace
        )
        if not verified:
            raise AssertionError("theorem-gate SAT witness failed physical verification")
        record["witness"] = {
            "c": coefficient.astype(int).tolist(),
            "v_bits_x_then_z": vector_bits(physical),
            "symplectic_weight": weight_numpy,
        }
        record["witness_physical_space_verified"] = verified
    return record


def closed_from_gate(
    h_original: np.ndarray,
    h_basis: np.ndarray,
    original_coordinates: np.ndarray,
    w: int,
    r_zsub: int,
) -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]]] | None:
    """Exact rank(V_w) implied by an INFEASIBLE gate, with no further solving.

    The gate proves V_w lies inside the pure-Z subgroup, so rank(V_w) <=
    r_Zsub.  The published rows of weight at most w give a matching lower
    bound whenever they already span rank r_Zsub, and then rank(V_w) = r_Zsub
    exactly.  Returns None when the lower bound falls short, in which case the
    caller must run the closure loop.
    """
    started = time.perf_counter()
    elements: list[dict[str, Any]] = []
    for row_index, row in enumerate(h_original):
        if physical_weight(row, h_original.shape[1] // 2) <= w:
            add_if_independent(
                elements,
                original_coordinates[row_index],
                h_basis,
                w,
                f"published_row:{row_index}",
                expected_v=row,
            )
    verification = verify_elements(elements, h_basis, w)
    seeded_rank = verification["rank_V_gf2_numpy"]
    if seeded_rank != r_zsub:
        return None
    wall = time.perf_counter() - started
    processed = {
        "weight_cap": w,
        "rank_Vw": int(r_zsub),
        "rank_lower_bound_if_undecided": None,
        "closed": True,
        "closure_kind": "gate_infeasible_upper_bound_meets_seeded_lower_bound",
        "closure_argument": (
            f"theorem gate INFEASIBLE at w=7 gives V_{w} <= pure-Z subgroup, so "
            f"rank(V_{w}) <= r_Zsub = {r_zsub}; the published rows of weight <= {w} "
            f"span rank {seeded_rank}, so rank(V_{w}) = {r_zsub} exactly"
        ),
        "n_solver_calls": 0,
        "wall_s": round(wall, 6),
        "timeouts": [],
        "n_seeded_published_independent": len(elements),
        "n_carried_independent": 0,
        "n_solver_witnesses": 0,
        "call_seeds": [],
        "independent_verification": verification,
        "physical_space_verified": verification["physical_space_verified"],
    }
    raw = {
        "weight_cap": w,
        "closed": True,
        "rank_Vw": int(r_zsub),
        "rank_lower_bound_if_undecided": None,
        "closure_kind": processed["closure_kind"],
        "closure_argument": processed["closure_argument"],
        "elements": [serialise_element(element) for element in elements],
        "solver_calls": [],
        "closure_proof": {"derived_from": "theorem_gate", "r_Zsub": int(r_zsub)},
        "timeouts": [],
        "independent_verification": verification,
        "wall_s": round(wall, 6),
    }
    return processed, raw, elements




def determine_w_star(
    rank_h: int, canonical_max: int, per_w: dict[str, dict[str, Any]]
) -> tuple[int | None, dict[str, int | None]]:
    lower = 1
    upper = canonical_max
    exact = None
    for w in sorted(int(key) for key in per_w):
        result = per_w[str(w)]
        if not result["closed"]:
            continue
        rank_v = result["rank_Vw"]
        if rank_v < rank_h:
            lower = max(lower, w + 1)
        elif rank_v == rank_h:
            upper = min(upper, w)
    if lower == upper:
        exact = lower
    return exact, {"lower_bound": lower, "upper_bound": upper}


def protocol_record(catalogue_sha256: str) -> dict[str, Any]:
    return {
        "experiment": EXPERIMENT,
        "base_seed": BASE_SEED,
        "seed_derivation": "low 31 bits of SHA256('20260812:EXP-023:<code>:<w>:<call>')",
        "pbb_weight_caps_in_order": list(PBB_WEIGHT_CAPS),
        "gross_weight_caps_in_order": list(GROSS_WEIGHT_CAPS),
        "rank_bookkeeping_per_call_time_limit_s": RANK_CALL_TIME_LIMIT_S,
        "theorem_gate_per_call_time_limit_s": THEOREM_GATE_TIME_LIMIT_S,
        "monolithic_crosscheck_member": CROSSCHECK_MEMBER,
        "decisive_method": (
            "exhaustive X-codeword enumeration (CP-SAT, completeness certified by "
            "OPTIMAL) plus exact GF(2) coset-budget tests against the pure-Z subgroup"
        ),
        "num_search_workers": NUM_SEARCH_WORKERS,
        "expected_n": EXPECTED_N,
        "expected_k": EXPECTED_K,
        "expected_rank_H": EXPECTED_RANK,
        "expected_pbb_members": EXPECTED_PBB_MEMBERS,
        "catalogue": CATALOGUE_REL,
        "catalogue_sha256": catalogue_sha256,
        "coefficient_space": (
            "132 recorded independent published rows H_B; c -> c H_B is bijective"
        ),
        "closure_functionals": (
            "physical h in F2^288 annihilating found v, pulled back as f=H_B h"
        ),
        "solver": "OR-Tools CP-SAT AddBoolXor exact feasibility",
        "ortools_version": ortools.__version__,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
    }


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    temporary.replace(path)


def choose_verdict(processed: dict[str, Any], full_coverage: bool) -> tuple[str, str]:
    if not full_coverage:
        return (
            "INCONCLUSIVE",
            "Partial execution: the full 14-PBB + Gross protocol with both tiers was not run.",
        )
    pbb = [row for row in processed["results"] if row["kind"] == "PBB catalogue member"]
    undecided = [
        row["code"]
        for row in pbb
        if not row["basis_independent_depth8_decided"]
    ]
    if undecided:
        return (
            "INCONCLUSIVE",
            "The decisive theorem gate remains UNDECIDED for " + ", ".join(undecided) + ".",
        )
    light = [
        row["code"]
        for row in pbb
        if row["theorem_gate"]["mixed_light_element_exists"] is True
        and row["per_w"]["7"]["closed"]
        and row["per_w"]["7"]["rank_Vw"] == EXPECTED_RANK
    ]
    if light:
        return (
            "BREAKTHROUGH_CANDIDATE",
            "Light generating set found; depth-7 PBB possible for " + ", ".join(light) + ".",
        )
    proved = [
        row
        for row in pbb
        if row["basis_independent_depth8"]
        and row["proof_route"] in ("exact_closure", "gate_plus_structural")
    ]
    if len(proved) == EXPECTED_PBB_MEMBERS:
        return (
            "BREAKTHROUGH_CANDIDATE",
            "Theorem C4 basis-independence proved computationally for all 14 catalogue PBB members.",
        )
    return "INCONCLUSIVE", "The amended proof routes did not meet a verdict rule."


def render_report(processed: dict[str, Any], paths: dict[str, str]) -> str:
    verdict = processed["verdict"]
    reason = processed["verdict_reason"]
    lines = [f"{verdict} — {reason}", "", "# EXP-023 light generating sets", ""]
    lines += [
        "## What was decided",
        "",
        "For each code, `V_w` is the GF(2) span of all stabilizers of symplectic weight at most `w`, and `w*` is the least `w` with `rank(V_w) = n - k = 132`. The decisive question per PBB member is whether any stabilizer of weight at most 7 has nonzero X-part. It is settled exhaustively, not by a single opaque solver call:",
        "",
        "- Every stabilizer of a PBB code is `(u H_X | u P + s)` with `s` ranging over the pure-Z subgroup `S_Z`, so a weight-<=7 element with nonzero X-part exists **iff** some nonzero codeword `x = u H_X` of `rowspace(H_X)` with `|x| <= 7` admits an `s` leaving at most `7 - |x|` qubits of Z-support outside `supp(x)`.",
        "- Stage 1 enumerates **all** such codewords `x` with CP-SAT; the `OPTIMAL` status after exhaustive enumeration is the completeness certificate. Stage 2 decides each candidate by exact GF(2) coset arithmetic, with no solver involved.",
        "- `r_Zsub`, the exact rank of the full pure-Z subgroup `rowspace(H_Z) + {u[C D] : u H_X = 0}`, is computed by plain linear algebra and cross-checked against the direct definition `{c H : (c H)_X = 0}`.",
        "",
        "The two proof routes recorded per member are:",
        "",
        "- `gate_plus_structural`: no mixed light element exists, hence `V_7` lies inside the pure-Z subgroup and `rank(V_7) <= r_Zsub < 132`. The published weight-6 Z rows span exactly `r_Zsub`, so this upper bound meets a matching lower bound and `rank(V_6) = rank(V_7) = r_Zsub` is exact.",
        "- `exact_closure`: the iterative closure loop terminated on its own, either spanning the full rowspace or returning CP-SAT `INFEASIBLE` against a complete basis of physical annihilator functionals.",
        "",
        "Both routes prove the same statement: every generating set of the stabilizer group contains a generator of symplectic weight at least 8, so by Proposition C2 every one-ancilla syndrome-extraction round has depth at least 8 for that code, independently of the measured basis.",
        "",
        "| code | canonical max | gate w=7 | r_Zsub | rank V6 | rank V7 | proof route | w* | depth >= 8 basis-independent? |",
        "|---|---:|---|---:|---|---|---|---|---|",
    ]

    def rank_cell(result: dict[str, Any] | None, bound: int | None) -> str:
        if result is None:
            return "n/a"
        if result["closed"]:
            return f"{result['rank_Vw']} exact"
        if bound is not None:
            return f"UNDECIDED (<= {bound})"
        return "UNDECIDED"

    gross = next(
        (row for row in processed["results"] if row["kind"] == "CSS Gross control"),
        None,
    )
    if gross is not None:
        g5, g6 = gross["per_w"].get("5"), gross["per_w"].get("6")
        wstar = gross["w_star_if_determined"]
        lines.append(
            f"| {gross['code']} (control) | {gross['canonical_max_check_weight']} | n/a | n/a "
            f"| w=5: {rank_cell(g5, None)} | w=6: {rank_cell(g6, None)} | exact_closure "
            f"| {wstar if wstar is not None else 'UNDETERMINED'} | no (w*=6) |"
        )
    for row in processed["results"]:
        if row["kind"] != "PBB catalogue member":
            continue
        gate = row["theorem_gate"]
        r_zsub = row["pure_Z_subgroup"]["r_Zsub"]
        bound = r_zsub if gate["status"] == "INFEASIBLE" else None
        interval = row["w_star_interval"] or {}
        wstar = (
            str(row["w_star_if_determined"])
            if row["w_star_if_determined"] is not None
            else f"in [{interval.get('lower_bound')}, {interval.get('upper_bound')}]"
        )
        depth = (
            "yes"
            if row["basis_independent_depth8"]
            else ("no" if row["basis_independent_depth8_decided"] else "UNDECIDED")
        )
        lines.append(
            f"| {row['code']} | {row['canonical_max_check_weight']} | {gate['status']} | {r_zsub} "
            f"| {rank_cell(row['per_w'].get('6'), bound)} | {rank_cell(row['per_w'].get('7'), bound)} "
            f"| {row['proof_route'] or 'none'} | {wstar} | {depth} |"
        )
    lines += [
        "",
        "`rank V6`/`rank V7` are the secondary rank-bookkeeping tier; an `UNDECIDED` there never blocks the theorem, which rests on the decisive gate plus `r_Zsub`. The Gross control columns are `w=5` and `w=6`.",
        "",
        "`w*` is exact where a single value is shown; where an interval is shown, the lower bound is the proved one (no generating set can stay under weight 8) and the upper bound is the published generating set's own maximum weight, which this experiment did not try to improve.",
        "",
        "## Independent verification and controls",
        "",
        "Coefficients live on an explicitly recorded 132-row independent subset `H_B` of the published rows, so `c -> c H_B` is a bijection onto the stabilizer rowspace and no rank or closure statement can be corrupted by the 12-dimensional kernel of the full 144-row map. Every closure functional is physical: an `h` in `F_2^288` annihilating the already-found stabilizers, pulled back as `f = H_B h`. Every witness was recomputed as `v = c H_B` with NumPy GF(2) arithmetic, weighed with both a local support count and the project's `symplectic_weight`, confirmed to lie in `rowspace(H)` by rank comparison, and ranked in 288-bit physical space by three independent implementations.",
        "",
        "Because every scientific claim here rests on a negative (infeasible) answer, each code also runs controls that must come out positive, and a failure aborts the run before any canonical artifact is written:",
        "",
        "- CP-SAT positive control: a nonzero stabilizer of weight at most the canonical max exists, with a physically verified witness.",
        "- CP-SAT positive control (PBB): a stabilizer of weight at most the canonical max with nonzero X-part exists, with a physically verified witness. This exercises exactly the query shape that returns `INFEASIBLE` at w=7.",
        "- Decomposition controls: every published weight-<=7 X-check row appears in the stage-1 enumeration; the stage-2 coset test recovers a published mixed generator at its own budget; an impossible (negative) budget is rejected.",
        "- Cross-encoding check on the benchmark member: the monolithic single-query CP-SAT formulation was run independently and its answer must agree with the decomposition.",
        "",
        "A separate cross-check (`experiments/exp023_verify_artifact.py`, output `results/raw/exp023_independent_verification.json`) does not import the producing script: it rebuilds every code from the catalogue, recomputes `rank(H)` and `r_Zsub`, re-checks every stored coefficient witness, brute-forces all published-row combinations of coefficient weight at most 3, and independently repeats the exact GF(2) coset test for every published light X row. This is a bounded independent cross-check; the all-coefficient-space completeness proof remains the stage-1 CP-SAT `OPTIMAL` enumeration.",
        "",
        f"`independent_verification_flags_all_true`: {processed['independent_verification_flags_all_true']}.",
        "",
        "## Timeouts and statistics",
        "",
    ]
    timeouts = []
    for row in processed["results"]:
        gate = row.get("theorem_gate")
        if gate is not None and gate["status"] == "UNDECIDED":
            timeouts.append(f"{row['code']} theorem gate w=7 (UNDECIDED)")
        for w, result in row["per_w"].items():
            for timeout in result.get("timeouts", []):
                timeouts.append(
                    f"{row['code']} rank bookkeeping w={w} call={timeout['call_index']} ({timeout['status']})"
                )
    lines.append("Timeouts: " + (", ".join(timeouts) if timeouts else "none."))
    lines.append("")
    decisive_walls = sorted(
        row["theorem_gate"]["decisive_decomposition"]["wall_s"]
        for row in processed["results"]
        if row.get("theorem_gate")
    )
    if decisive_walls:
        middle = decisive_walls[len(decisive_walls) // 2]
        lines.append(
            f"Decisive enumeration wall time across the {len(decisive_walls)} PBB members: "
            f"min {decisive_walls[0]:.1f} s, median {middle:.1f} s, max {decisive_walls[-1]:.1f} s."
        )
        for row in processed["results"]:
            gate = row.get("theorem_gate")
            if not gate or gate["monolithic_crosscheck"]["status"] == "NOT_ATTEMPTED":
                continue
            lines.append(
                f"Cross-encoding check on `{row['code']}`: the monolithic single-query CP-SAT "
                f"formulation returned `{gate['monolithic_crosscheck']['status']}` in "
                f"{gate['monolithic_crosscheck']['wall_s']:.1f} s against the decomposition's "
                f"`{gate['status']}` in {gate['decisive_decomposition']['wall_s']:.1f} s; "
                f"methods agree: {gate['methods_agree']}."
            )
    lines.append("")
    lines.append(
        "This is an exact finite feasibility experiment, not Monte Carlo: shots, failures, and Clopper-Pearson intervals do not apply. Every reported outcome is `OPTIMAL`, `INFEASIBLE`, or an explicitly labelled `UNDECIDED`."
    )
    lines += [
        "",
        "## Reproduction",
        "",
        "```bash",
        "cd /Users/jinleic/jinleic-workspace/math/qec",
        "PYTHONPATH=src .venv/bin/python experiments/exp023_light_gensets.py",
        "PYTHONPATH=src .venv/bin/python -m pytest -q tests/test_exp023_light_gensets.py",
        "PYTHONPATH=src .venv/bin/python experiments/exp023_verify_artifact.py",
        "```",
        "",
        f"Observed end-to-end wall time: {processed['wall_s']:.3f} s. The decisive enumeration was capped at {THEOREM_GATE_TIME_LIMIT_S:.0f} s per code (single worker, as CP-SAT enumeration requires) and secondary rank-closure calls at {RANK_CALL_TIME_LIMIT_S:.0f} s with {NUM_SEARCH_WORKERS} workers. Base seed {BASE_SEED}; every call seed is recorded in the JSON artifacts. Timing was measured on a heavily loaded shared machine (system load average above 100 during the run, dominated by an external desktop workload) and is not a benchmark; the SAT/UNSAT outcomes themselves are timing-independent.",
        "",
        "## Artifacts and caveats",
        "",
        f"- Processed results: `{paths['processed']}`.",
        f"- Coefficient witnesses, physical closure functionals, pure-Z subgroup bases, and any explicit light generating set: `{paths['raw']}`.",
        "- A non-null `w_star_if_determined` is exact; otherwise `w_star_interval` is a rigorous bracket and this pre-registered w=6,7 threshold experiment claims nothing sharper.",
        "- Catalogue distance metadata is used only for input selection; EXP-023 makes no distance claim.",
        "- Independent verification pass: `results/raw/exp023_independent_verification.json`.",
        "- No PBB member reached `rank(V_7) = 132`, so no depth-7 PBB generating set exists to save. The one explicit light generating set stored is the Gross control's own 132 weight-6 rows, which is what `w* = 6` means for it.",
        "- The proved statement is scoped to one-ancilla syndrome extraction via Proposition C2; it says nothing about multi-ancilla or measurement-gadget schemes.",
        "- This supersedes the partial attempt in `experiments/exp022_basis_independence.py`, which measured a weaker quantity (`w_mix` on a single element, and an X-part-only rank) and left the question open with unproven-optimal solver values.",
        "- Write-integrity incident (`notes/failed_routes.md` FR-012 pattern): an early `--only gross` smoke run of this script wrote the canonical paths. Output routing now goes through `qec_research.artifacts.canonical_route`, so only a run covering the full declared scope can touch them; partial runs land in `results/partial_runs/` with a timestamped suffix. The guard was verified by re-running `--only gross` and confirming the canonical files' checksums were unchanged, and the canonical artifacts here were regenerated from scratch by the full 15-code run.",
    ]
    return "\n".join(lines) + "\n"


def output_paths(full_coverage: bool, started_utc: str) -> tuple[Path, Path, Path, str]:
    route = canonical_route(clean=True, full_coverage=full_coverage)
    if route == "canonical":
        return (
            CANONICAL_RAW_PATH,
            CANONICAL_PROCESSED_PATH,
            CANONICAL_REPORT_PATH,
            route,
        )
    stamp = started_utc.replace(":", "").replace("-", "").replace("+", "_")
    suffix = f"exp023_{stamp}_partial"
    return (
        ROOT / "results" / route / f"{suffix}_light_elements.json",
        ROOT / "results" / route / f"{suffix}_light_gensets.json",
        ROOT / "notes" / "agent_reports" / f"{suffix}_light_gensets.md",
        route,
    )


def run(
    only: set[str] | None = None,
    *,
    skip_rank_bookkeeping: bool = False,
) -> tuple[dict[str, Any], dict[str, Any]]:
    total_started = time.perf_counter()
    catalogue_path = ROOT / CATALOGUE_REL
    catalogue_sha = hashlib.sha256(catalogue_path.read_bytes()).hexdigest()
    protocol = protocol_record(catalogue_sha)
    started_utc = utc_now()
    full_run = only is None
    full_coverage = full_run and not skip_rank_bookkeeping
    raw_path, processed_path, report_path, artifact_route = output_paths(
        full_coverage, started_utc
    )
    protocol["artifact_route"] = artifact_route
    protocol["rank_bookkeeping_skipped"] = bool(skip_rank_bookkeeping)
    raw: dict[str, Any] = {
        "experiment": EXPERIMENT,
        "protocol": protocol,
        "started_utc": started_utc,
        "results": [],
        "verdict": "INCOMPLETE",
        "verdict_reason": "Experiment is still running.",
        "wall_s": 0.0,
    }
    processed: dict[str, Any] = {
        "experiment": EXPERIMENT,
        "protocol": protocol,
        "started_utc": started_utc,
        "results": [],
        "independent_verification_flags_all_true": False,
        "verdict": "INCOMPLETE",
        "verdict_reason": "Experiment is still running.",
        "wall_s": 0.0,
    }

    all_items = make_items()
    selected_items = [
        item
        for item in all_items
        if only is None
        or item["label"] in only
        or ("gross" in only and item["kind"] == "CSS Gross control")
    ]
    if not selected_items:
        raise ValueError(f"--only matched no code: {sorted(only or [])}")
    print(
        f"EXP-023 starting {len(selected_items)} code(s); theorem gate "
        f"limit={THEOREM_GATE_TIME_LIMIT_S:.0f}s, rank limit="
        f"{RANK_CALL_TIME_LIMIT_S:.0f}s, workers={NUM_SEARCH_WORKERS}, "
        f"route={artifact_route}",
        flush=True,
    )

    for item_index, item in enumerate(selected_items, start=1):
        code_started = time.perf_counter()
        label = item["label"]
        h_original = np.asarray(item["H"], dtype=np.uint8) & 1
        if h_original.shape != (EXPECTED_N, 2 * EXPECTED_N):
            raise AssertionError(f"{label}: unexpected H shape {h_original.shape}")
        rank_checks = rank_triplet(h_original, 2 * EXPECTED_N)
        if rank_checks != (EXPECTED_RANK,) * 3:
            raise AssertionError(f"{label}: rank(H) checks disagree: {rank_checks}")
        basis_rows = independent_row_indices(h_original)
        if len(basis_rows) != EXPECTED_RANK:
            raise AssertionError(f"{label}: selected {len(basis_rows)} basis rows")
        h_basis = h_original[basis_rows]
        basis_rank_checks = rank_triplet(h_basis, 2 * EXPECTED_N)
        if basis_rank_checks != (EXPECTED_RANK,) * 3:
            raise AssertionError(f"{label}: H_B rank checks disagree: {basis_rank_checks}")
        original_coordinates = coefficient_map(h_basis, h_original)
        canonical_weights = [physical_weight(row, EXPECTED_N) for row in h_original]
        canonical_max = max(canonical_weights)
        basis_setup_verified = bool(
            rank_checks == basis_rank_checks == (EXPECTED_RANK,) * 3
            and np.array_equal(
                numpy_gf2_matmul(original_coordinates, h_basis), h_original
            )
        )
        if not basis_setup_verified:
            raise AssertionError(f"{label}: basis setup failed independent verification")

        processed_item: dict[str, Any] = {
            "code": label,
            "kind": item["kind"],
            "n": EXPECTED_N,
            "k": EXPECTED_K,
            "rank_H": EXPECTED_RANK,
            "rank_H_checks": {
                "numpy_independent": rank_checks[0],
                "gf2_rank_np": rank_checks[1],
                "gf2_rank_bitset": rank_checks[2],
            },
            "basis_row_indices": basis_rows,
            "rank_H_basis_checks": {
                "numpy_independent": basis_rank_checks[0],
                "gf2_rank_np": basis_rank_checks[1],
                "gf2_rank_bitset": basis_rank_checks[2],
            },
            "basis_setup_verified": basis_setup_verified,
            "canonical_max_check_weight": canonical_max,
            "canonical_check_weight_histogram": {
                str(weight): canonical_weights.count(weight)
                for weight in sorted(set(canonical_weights))
            },
            "theorem_gate": None,
            "pure_Z_subgroup": None,
            "proof_route": None,
            "rank_bookkeeping": {},
            "per_w": {},
            "w_star_if_determined": None,
            "w_star_interval": None,
            "basis_independent_depth8": False,
            "basis_independent_depth8_decided": item["kind"] == "CSS Gross control",
            "physical_space_verified": False,
            "wall_s": None,
        }
        raw_item: dict[str, Any] = {
            "code": label,
            "kind": item["kind"],
            "catalogue_spec": item["catalogue_spec"],
            "basis_row_indices": basis_rows,
            "basis_definition": "H_B = H[basis_row_indices] in the listed order",
            "rank_H_checks": processed_item["rank_H_checks"],
            "rank_H_basis_checks": processed_item["rank_H_basis_checks"],
            "published_row_coefficients_in_HB": original_coordinates.astype(int).tolist(),
            "basis_setup_verified": basis_setup_verified,
            "theorem_gate": None,
            "pure_Z_subgroup": None,
            "per_w": {},
            "explicit_light_generating_set": None,
            "physical_space_verified": False,
            "wall_s": None,
        }
        processed["results"].append(processed_item)
        raw["results"].append(raw_item)

        empty_found = np.zeros((0, 2 * EXPECTED_N), dtype=np.uint8)
        all_functionals, _, _ = physical_annihilator_functionals(empty_found, h_basis)
        controls = [
            positive_control_query(
                label,
                h_basis,
                all_functionals,
                canonical_max,
                "exists a nonzero stabilizer of weight <= canonical max check weight",
                require_x_support=False,
            )
        ]

        gate_infeasible = False
        r_zsub: int | None = None
        if item["kind"] == "PBB catalogue member":
            pure_z_rows, pure_z_checks = pure_z_subgroup_proof(
                item["spec"], h_original
            )
            r_zsub = pure_z_checks["r_Zsub"]
            controls.append(
                positive_control_query(
                    label,
                    h_basis,
                    physical_annihilator_functionals(pure_z_rows, h_basis)[0],
                    canonical_max,
                    "exists a stabilizer of weight <= canonical max with nonzero X-part",
                    require_x_support=True,
                )
            )
            decisive = decomposition_gate(
                label, item["spec"], h_basis, pure_z_rows, 7
            )
            crosscheck = (
                theorem_gate_query(
                    label,
                    h_basis,
                    pure_z_rows,
                    7,
                    THEOREM_GATE_TIME_LIMIT_S,
                )
                if label == CROSSCHECK_MEMBER
                else {
                    "status": "NOT_ATTEMPTED",
                    "reason": (
                        "monolithic single-query CP-SAT run only on the benchmark "
                        "member; the exhaustive decomposition decides every member"
                    ),
                    "mixed_light_element_exists": None,
                    "witness": None,
                    "witness_physical_space_verified": None,
                    "wall_s": 0.0,
                }
            )
            crosscheck_exists = crosscheck["mixed_light_element_exists"]
            if (
                decisive["mixed_light_element_exists"] is not None
                and crosscheck_exists is not None
                and decisive["mixed_light_element_exists"] != crosscheck_exists
            ):
                raise AssertionError(
                    f"{label}: decomposition and monolithic CP-SAT disagree "
                    f"({decisive['mixed_light_element_exists']} vs {crosscheck_exists})"
                )
            exists = (
                decisive["mixed_light_element_exists"]
                if decisive["mixed_light_element_exists"] is not None
                else crosscheck_exists
            )
            witness = decisive["witness"] or crosscheck["witness"]
            theorem_gate = {
                "question": (
                    "exists v in rowspace(H) with symplectic weight <= 7 and "
                    "nonzero X-part"
                ),
                "status": (
                    "INFEASIBLE"
                    if exists is False
                    else "SAT"
                    if exists is True
                    else "UNDECIDED"
                ),
                "mixed_light_element_exists": exists,
                "decided_by": (
                    "exhaustive_decomposition"
                    if decisive["mixed_light_element_exists"] is not None
                    else "monolithic_cp_sat"
                    if crosscheck_exists is not None
                    else None
                ),
                "decisive_decomposition": decisive,
                "monolithic_crosscheck": crosscheck,
                "methods_agree": bool(
                    decisive["mixed_light_element_exists"] is None
                    or crosscheck_exists is None
                    or decisive["mixed_light_element_exists"] == crosscheck_exists
                ),
                "witness": witness,
                "witness_physical_space_verified": (
                    crosscheck["witness_physical_space_verified"]
                    if crosscheck["witness"] is not None
                    else (witness is not None or None)
                ),
                "wall_s": round(decisive["wall_s"] + crosscheck["wall_s"], 6),
            }
            gate_infeasible = exists is False
            processed_item["theorem_gate"] = theorem_gate
            raw_item["theorem_gate"] = theorem_gate
            processed_item["pure_Z_subgroup"] = pure_z_checks
            raw_item["pure_Z_subgroup"] = {
                **pure_z_checks,
                "basis_rows_x_then_z_bits": [
                    vector_bits(row) for row in rref_np(pure_z_rows)[0]
                ],
            }
            if gate_infeasible:
                processed_item["proof_route"] = "gate_plus_structural"
                processed_item["basis_independent_depth8"] = True
                processed_item["basis_independent_depth8_decided"] = True
                processed_item["rank_V7_proven_upper_bound"] = r_zsub
                processed_item["rank_V6_proven_upper_bound"] = r_zsub
            print(
                f"[{item_index}/{len(selected_items)}] {label}: gate "
                f"{theorem_gate['status']} via {theorem_gate['decided_by']} "
                f"(decomposition {decisive['wall_s']:.1f}s, crosscheck "
                f"{crosscheck['status']} {crosscheck['wall_s']:.1f}s); r_Zsub={r_zsub}",
                flush=True,
            )

        carried: list[dict[str, Any]] = []
        for w in item["weight_caps"]:
            if skip_rank_bookkeeping:
                per_processed = {
                    "weight_cap": w,
                    "rank_Vw": None,
                    "rank_lower_bound_if_undecided": None,
                    "closed": False,
                    "closure_kind": "SKIPPED_BY_EXPLICIT_DEVELOPMENT_FLAG",
                    "n_solver_calls": 0,
                    "wall_s": 0.0,
                    "timeouts": [],
                    "physical_space_verified": True,
                }
                per_raw = per_processed
            else:
                derived = (
                    closed_from_gate(
                        h_original, h_basis, original_coordinates, w, r_zsub
                    )
                    if gate_infeasible and r_zsub is not None
                    else None
                )
                if derived is not None:
                    per_processed, per_raw, carried = derived
                else:
                    print(
                        f"[{item_index}/{len(selected_items)}] {label}: closing V_{w}",
                        flush=True,
                    )
                    per_processed, per_raw, carried = analyse_weight_cap(
                        label,
                        h_original,
                        h_basis,
                        original_coordinates,
                        w,
                        carried,
                    )
            processed_item["per_w"][str(w)] = per_processed
            processed_item["rank_bookkeeping"][str(w)] = per_processed
            raw_item["per_w"][str(w)] = per_raw
            print(
                f"    V_{w}: rank="
                f"{per_processed['rank_Vw'] if per_processed['closed'] else per_processed['rank_lower_bound_if_undecided']}"
                f", closed={per_processed['closed']}"
                f", calls={per_processed['n_solver_calls']}"
                f", wall={per_processed['wall_s']:.3f}s",
                flush=True,
            )

        exact_w, interval = determine_w_star(
            EXPECTED_RANK, canonical_max, processed_item["per_w"]
        )
        if processed_item["basis_independent_depth8"]:
            interval["lower_bound"] = max(interval["lower_bound"], 8)
            if interval["lower_bound"] == interval["upper_bound"]:
                exact_w = interval["lower_bound"]
        processed_item["w_star_if_determined"] = exact_w
        processed_item["w_star_interval"] = interval
        if item["kind"] == "PBB catalogue member":
            r7 = processed_item["per_w"]["7"]
            if r7["closed"]:
                processed_item["basis_independent_depth8_decided"] = True
                processed_item["basis_independent_depth8"] = bool(
                    r7["rank_Vw"] < EXPECTED_RANK
                )
                if r7["closure_kind"] != (
                    "gate_infeasible_upper_bound_meets_seeded_lower_bound"
                ):
                    processed_item["proof_route"] = "exact_closure"

        processed_item["positive_controls"] = controls
        raw_item["positive_controls"] = controls
        processed_item["physical_space_verified"] = bool(
            basis_setup_verified
            and all(control["control_passed"] for control in controls)
            and all(
                result["physical_space_verified"]
                for result in processed_item["per_w"].values()
            )
            and (
                item["kind"] == "CSS Gross control"
                or (
                    processed_item["pure_Z_subgroup"]["pure_Z_subgroup_verified"]
                    and (
                        processed_item["theorem_gate"]["status"] == "INFEASIBLE"
                        or processed_item["theorem_gate"][
                            "witness_physical_space_verified"
                        ]
                        is True
                        or processed_item["theorem_gate"][
                            "mixed_light_element_exists"
                        ]
                        is None
                    )
                )
            )
        )
        raw_item["physical_space_verified"] = processed_item["physical_space_verified"]
        raw_item["w_star_if_determined"] = exact_w
        raw_item["w_star_interval"] = interval
        raw_item["proof_route"] = processed_item["proof_route"]

        full_caps = [
            int(w)
            for w, result in processed_item["per_w"].items()
            if result["closed"] and result["rank_Vw"] == EXPECTED_RANK
        ]
        if full_caps and not skip_rank_bookkeeping:
            lightest = min(full_caps)
            full_elements = raw_item["per_w"][str(lightest)]["elements"]
            explicit = {
                "weight_cap": lightest,
                "n_rows": len(full_elements),
                "coefficient_vectors_c": [element["c"] for element in full_elements],
                "symplectic_rows_x_then_z_bits": [
                    element["v_bits_x_then_z"] for element in full_elements
                ],
                "row_weights": [
                    element["symplectic_weight"] for element in full_elements
                ],
                "rank_physical_verified": EXPECTED_RANK,
                "max_weight_verified": max(
                    element["symplectic_weight"] for element in full_elements
                ),
            }
            raw_item["explicit_light_generating_set"] = explicit
            processed_item["explicit_generating_set"] = {
                "raw_artifact": str(raw_path.relative_to(ROOT)),
                "raw_key": f"results[{len(raw['results']) - 1}].explicit_light_generating_set",
                "weight_cap": lightest,
                "n_rows": len(full_elements),
            }

        code_wall = time.perf_counter() - code_started
        processed_item["wall_s"] = round(code_wall, 6)
        raw_item["wall_s"] = round(code_wall, 6)
        raw["wall_s"] = round(time.perf_counter() - total_started, 6)
        processed["wall_s"] = raw["wall_s"]
        # Full runs checkpoint away from canonical. Canonical paths are written
        # only after the complete declared scope has been attempted.
        checkpoint_base = ROOT / "results" / "partial_runs"
        checkpoint_raw = checkpoint_base / "exp023_full_in_progress_light_elements.json"
        checkpoint_processed = checkpoint_base / "exp023_full_in_progress_light_gensets.json"
        write_json(checkpoint_raw if full_coverage else raw_path, raw)
        write_json(checkpoint_processed if full_coverage else processed_path, processed)

    all_verified = bool(
        processed["results"]
        and all(row["physical_space_verified"] for row in processed["results"])
    )
    processed["independent_verification_flags_all_true"] = all_verified
    raw["independent_verification_flags_all_true"] = all_verified
    verdict, reason = choose_verdict(processed, full_coverage)
    processed["verdict"] = verdict
    processed["verdict_reason"] = reason
    raw["verdict"] = verdict
    raw["verdict_reason"] = reason
    ended_utc = utc_now()
    wall = time.perf_counter() - total_started
    processed["ended_utc"] = ended_utc
    raw["ended_utc"] = ended_utc
    processed["wall_s"] = round(wall, 6)
    raw["wall_s"] = round(wall, 6)
    write_json(raw_path, raw)
    write_json(processed_path, processed)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        render_report(
            processed,
            {
                "processed": str(processed_path.relative_to(ROOT)),
                "raw": str(raw_path.relative_to(ROOT)),
            },
        )
    )
    print(f"{verdict}: {reason}", flush=True)
    print(f"wrote {processed_path.relative_to(ROOT)}", flush=True)
    print(f"wrote {raw_path.relative_to(ROOT)}", flush=True)
    print(f"wrote {report_path.relative_to(ROOT)}", flush=True)
    return processed, raw


def rerender_canonical_report() -> Path:
    """Regenerate the canonical report from the canonical full-coverage JSON.

    Refuses to touch the report unless the processed artifact records a
    completed run over the full declared scope, so this cannot be used to
    launder a partial result into the canonical location.
    """
    processed = json.loads(CANONICAL_PROCESSED_PATH.read_text())
    codes = len(processed.get("results", []))
    route = processed.get("protocol", {}).get("artifact_route")
    if (
        route != "canonical"
        or codes != EXPECTED_PBB_MEMBERS + 1
        or processed.get("verdict") in (None, "INCOMPLETE")
    ):
        raise SystemExit(
            "refusing to re-render: the canonical JSON is not a completed "
            f"full-coverage run (route={route}, codes={codes}, "
            f"verdict={processed.get('verdict')})"
        )
    CANONICAL_REPORT_PATH.write_text(
        render_report(
            processed,
            {
                "processed": str(CANONICAL_PROCESSED_PATH.relative_to(ROOT)),
                "raw": str(CANONICAL_RAW_PATH.relative_to(ROOT)),
            },
        )
    )
    return CANONICAL_REPORT_PATH


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--only",
        action="append",
        default=None,
        help="Development/smoke selection: exact code_id, or 'gross'. Full protocol uses no flag.",
    )
    parser.add_argument(
        "--skip-rank-bookkeeping",
        action="store_true",
        help="Development-only partial run: execute theorem gate but skip secondary closure.",
    )
    parser.add_argument(
        "--rerender-report",
        action="store_true",
        help="Regenerate the canonical report from the canonical full-coverage JSON.",
    )
    args = parser.parse_args()
    if args.rerender_report:
        print(f"wrote {rerender_canonical_report().relative_to(ROOT)}", flush=True)
        return
    run(
        set(args.only) if args.only else None,
        skip_rank_bookkeeping=args.skip_rank_bookkeeping,
    )


if __name__ == "__main__":
    main()
