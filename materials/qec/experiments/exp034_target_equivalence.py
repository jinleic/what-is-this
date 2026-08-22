"""EXP-034: exact CSS-equivalence audit of catalogue target ``12_6_0193``.

This target-only experiment rebuilds the published stabilizer matrix from the
literal catalogue terms.  It decides CSS equivalence under stabilizer row
operations, qubit permutations, and arbitrary local Hadamards.  It also decides
connectedness of the supplied stabilizer incidence representation.  The larger
arbitrary local-Clifford group is persisted as an exact finite CP-SAT/SAT
formulation but remains unresolved unless a complete witness or infeasibility
certificate is supplied by a later bounded solver run.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CATALOGUE_PATH = (
    ROOT / "third_party" / "qcode-discovery" / "results" / "campaign7_publication_merged.jsonl"
)
OUTPUT_PATH = ROOT / "results" / "processed" / "exp034_target_equivalence.json"
PARTIAL_LC_PATH = ROOT / "results" / "partial_runs" / "exp034_target_full_lc_attempt.json"
TARGET_ID = "12_6_0193"
EXPECTED_LINE = 38
EXPECTED_TERMS = {
    "ell": 12,
    "m": 6,
    "A_terms": [[1, 2], [10, 3], [10, 4]],
    "B_terms": [[0, 0], [1, 5], [11, 4]],
    "C_terms": [[1, 3], [10, 3]],
    "D_terms": [[1, 5], [10, 5]],
}
EXPECTED_CATALOGUE_SHA256 = "c0197f0cb5e03bfc94215555ed0fbfefaa72d1c0e3fb89a31fc5045c5f7e25ba"
EXPECTED_RAW_LINE_SHA256 = "0494a1c2f33e27638b9ded66f8648ab7a842d55450e4369bba5b11fca51f261f"
EXPECTED_H_SHA256 = "a1c5d317ef0218bc42fe1ca02347162ef09822f9096358a17e2ea00529b740db"

sys.path.insert(0, str(ROOT / "src"))

from qec_research.codes.bicycle import PBBSpec, build_pbb  # noqa: E402
from qec_research.equivalence.css import (  # noqa: E402
    css_rank_invariants,
    full_lc_cpsat_formulation,
    full_lc_linear_certificate,
    full_lc_parity_equations,
    local_h_affine_system,
    matrix_sha256_payload,
    permute_qubits,
    projection_invariant,
    pure_intersection_basis,
    solve_affine_system,
    stabilizer_direct_sum_certificate,
    stabilizer_incidence_components,
)
from qec_research.gf2.linalg import rank_bitset, rank_np, rows_to_bitsets  # noqa: E402


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def non_vacuousness_control() -> dict[str, Any]:
    """Persist the [[5,1,3]] control proving the linear route is not vacuous.

    The perfect code is LC-inequivalent to CSS (exhaustive over all 6^5
    per-qubit assignments) yet its linear relaxation is feasible, so a
    linear XOR contradiction is informative rather than automatic.
    """

    from itertools import product

    rows = ("XZZXI", "IXZZX", "XIXZZ", "ZXIXZ")
    n = len(rows[0])
    H = np.zeros((len(rows), 2 * n), dtype=np.uint8)
    for row_index, row in enumerate(rows):
        for qubit, pauli in enumerate(row):
            if pauli in "XY":
                H[row_index, qubit] = 1
            if pauli in "ZY":
                H[row_index, n + qubit] = 1
    linear = full_lc_linear_certificate(H)
    witness = None
    for choices in product(range(6), repeat=n):
        if projection_invariant(H, choices, full_lc=True):
            witness = list(choices)
            break
    control = {
        "code": list(rows),
        "n": n,
        "exhaustive_assignments_checked": 6 ** n,
        "lc_equivalent_assignment": witness,
        "linear_relaxation_feasible": linear["linear_relaxation_feasible"],
        "purpose": (
            "linear-FEASIBLE proves nothing, so this LC-inequivalent, "
            "linear-feasible control certifies that the target's linear "
            "infeasibility is informative rather than automatic"
        ),
    }
    if witness is not None or not linear["linear_relaxation_feasible"]:
        raise AssertionError(f"non-vacuousness control violated: {control!r}")
    return control


def load_target_catalogue_row() -> dict[str, Any]:
    """Load the unique literal target row and reject source drift."""

    raw_lines = CATALOGUE_PATH.read_bytes().splitlines(keepends=True)
    matches: list[tuple[int, dict[str, Any], bytes]] = []
    for line_number, raw_line in enumerate(raw_lines, start=1):
        if not raw_line.strip():
            continue
        row = json.loads(raw_line)
        if row.get("code_id") == TARGET_ID:
            matches.append((line_number, row, raw_line))
    if len(matches) != 1:
        raise AssertionError(f"expected one {TARGET_ID} row, found {len(matches)}")
    line_number, row, raw_line = matches[0]
    if line_number != EXPECTED_LINE:
        raise AssertionError(f"target moved from line {EXPECTED_LINE} to {line_number}")
    for field, expected in EXPECTED_TERMS.items():
        if row[field] != expected:
            raise AssertionError(f"target {field} drifted: {row[field]} != {expected}")
    if sha256_bytes(CATALOGUE_PATH.read_bytes()) != EXPECTED_CATALOGUE_SHA256:
        raise AssertionError("catalogue source hash drifted")
    if sha256_bytes(raw_line) != EXPECTED_RAW_LINE_SHA256:
        raise AssertionError("raw target-row hash drifted")
    return row


def build_target_code(row: dict[str, Any] | None = None):
    """Rebuild exact H from the literal PBB terms."""

    source = load_target_catalogue_row() if row is None else row
    spec = PBBSpec(
        ell=int(source["ell"]),
        m=int(source["m"]),
        A=[tuple(term) for term in source["A_terms"]],
        B=[tuple(term) for term in source["B_terms"]],
        C=[tuple(term) for term in source["C_terms"]],
        D=[tuple(term) for term in source["D_terms"]],
        name=TARGET_ID,
    )
    code = build_pbb(spec)
    if sha256_bytes(matrix_sha256_payload(code.H)) != EXPECTED_H_SHA256:
        raise AssertionError("rebuilt H hash drifted")
    return code


def _independent_incidence_components(H: np.ndarray) -> list[list[int]]:
    """Independent check/check-qubit BFS, separate from union-find helper."""

    n = H.shape[1] // 2
    support = H[:, :n] | H[:, n:]
    qubit_to_rows = [np.flatnonzero(support[:, qubit]).tolist() for qubit in range(n)]
    unseen = set(range(n))
    components: list[list[int]] = []
    while unseen:
        seed = min(unseen)
        component: set[int] = set()
        qubit_frontier = [seed]
        seen_rows: set[int] = set()
        while qubit_frontier:
            qubit = qubit_frontier.pop()
            if qubit in component:
                continue
            component.add(qubit)
            unseen.discard(qubit)
            for row in qubit_to_rows[qubit]:
                if row in seen_rows:
                    continue
                seen_rows.add(row)
                qubit_frontier.extend(
                    int(neighbor) for neighbor in np.flatnonzero(support[row])
                    if int(neighbor) not in component
                )
        components.append(sorted(component))
    return sorted(components, key=lambda component: (-len(component), component))


def _rank_crosscheck(H: np.ndarray) -> dict[str, int | bool]:
    n = H.shape[1] // 2
    numpy_rank = rank_np(H)
    bitset_rank = rank_bitset(rows_to_bitsets(H), 2 * n)
    return {
        "numpy": numpy_rank,
        "bitset": bitset_rank,
        "agree": numpy_rank == bitset_rank,
    }


def _local_h_result(H: np.ndarray) -> dict[str, Any]:
    equations, nvars = local_h_affine_system(H)
    solved = solve_affine_system(equations, nvars)
    contradiction = None
    if solved.contradiction is not None:
        xor_mask = 0
        xor_rhs = 0
        encoded = []
        for mask, rhs in solved.contradiction:
            xor_mask ^= mask
            xor_rhs ^= rhs
            encoded.append({"coefficient_mask_hex": hex(mask), "rhs": rhs})
        if (xor_mask, xor_rhs) != (0, 1):
            raise AssertionError("local-H contradiction certificate is invalid")
        contradiction = {
            "equations": encoded,
            "xor_coefficient_mask_hex": hex(xor_mask),
            "xor_rhs": xor_rhs,
        }
    witness_valid = None
    if solved.assignment is not None:
        witness_valid = projection_invariant(H, solved.assignment, full_lc=False)
        if not witness_valid:
            raise AssertionError("local-H witness failed direct row-space verification")
    return {
        "group": "independent H or I on each of the 144 labelled qubits, plus arbitrary stabilizer row operations",
        "status": "INEQUIVALENT" if not solved.feasible else "EQUIVALENT",
        "complete": True,
        "method": "exact affine GF(2) projection-invariance feasibility",
        "num_variables": nvars,
        "num_distinct_equations": len(equations),
        "coefficient_rank": solved.rank,
        "augmented_rank": solved.augmented_rank,
        "witness": list(solved.assignment) if solved.assignment is not None else None,
        "witness_verified": witness_valid,
        "infeasibility_certificate": contradiction,
        "criterion": (
            "after local H choices, CSS iff S is invariant under the conjugated X-axis projection; "
            "basis/nullspace containment is an affine GF(2) system"
        ),
    }


def solve_full_lc_cpsat(
    H: np.ndarray,
    *,
    time_limit_s: float,
    memory_limit_mb: int,
    workers: int = 1,
    seed: int = 20260813,
) -> dict[str, Any]:
    """Run the complete finite full-LC feasibility model with explicit bounds."""

    if time_limit_s <= 0:
        raise ValueError("time_limit_s must be positive")
    if memory_limit_mb <= 0:
        raise ValueError("memory_limit_mb must be positive")
    if workers <= 0:
        raise ValueError("workers must be positive")

    from ortools.sat.python import cp_model

    matrix = np.asarray(H, dtype=np.uint8) & 1
    equations, n = full_lc_parity_equations(matrix)
    model = cp_model.CpModel()
    choices = [
        [model.NewBoolVar(f"q{qubit}_g{gate}") for gate in range(6)]
        for qubit in range(n)
    ]
    for qubit_choices in choices:
        model.AddExactlyOne(qubit_choices)

    # CP-SAT AddBoolXOr enforces odd parity.  Appending a fixed-true literal
    # turns it into the required even-parity equation.
    parity_one = model.NewBoolVar("parity_one")
    model.Add(parity_one == 1)
    for mask in equations:
        literals = []
        remaining = mask
        while remaining:
            least_bit = remaining & -remaining
            index = least_bit.bit_length() - 1
            literals.append(choices[index // 6][index % 6])
            remaining ^= least_bit
        model.AddBoolXOr(literals + [parity_one])

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(time_limit_s)
    solver.parameters.max_memory_in_mb = int(memory_limit_mb)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.random_seed = int(seed)
    status_code = solver.Solve(model)
    status_name = solver.StatusName(status_code)

    assignment = None
    witness_verified = None
    transformed_css = None
    if status_code in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        assignment = [
            next(gate for gate in range(6) if solver.Value(choices[qubit][gate]))
            for qubit in range(n)
        ]
        witness_verified = projection_invariant(matrix, assignment, full_lc=True)
        if not witness_verified:
            raise AssertionError("CP-SAT full-LC witness failed direct verification")

        from qec_research.equivalence.css import LOCAL_CLIFFORD_ACTIONS

        transformed = np.zeros_like(matrix)
        for qubit, gate in enumerate(assignment):
            columns = matrix[:, [qubit, n + qubit]]
            transformed[:, [qubit, n + qubit]] = (
                columns.astype(np.int64)
                @ LOCAL_CLIFFORD_ACTIONS[gate].astype(np.int64)
                % 2
            ).astype(np.uint8)
        transformed_css = css_rank_invariants(transformed)
        transformed_css["h_shape"] = list(transformed.shape)
        transformed_css["h_sha256"] = sha256_bytes(matrix_sha256_payload(transformed))
        if not transformed_css["css_after_row_operations"]:
            raise AssertionError("rebuilt local-Clifford transform is not CSS")

    linear = full_lc_linear_certificate(matrix)
    if status_code in (cp_model.OPTIMAL, cp_model.FEASIBLE) and not linear[
        "linear_relaxation_feasible"
    ]:
        raise AssertionError(
            "verified full-LC witness contradicts the linear relaxation"
        )

    if status_code == cp_model.INFEASIBLE:
        classification = "INEQUIVALENT"
        complete = True
    elif status_code in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        classification = "EQUIVALENT"
        complete = True
    else:
        classification = "UNRESOLVED"
        complete = False

    return {
        "group": (
            "independent choice of all six binary symplectic single-qubit Clifford "
            "actions on each labelled qubit, plus arbitrary stabilizer row operations"
        ),
        "status": classification,
        "complete": complete,
        "solver_status": status_name,
        "status_code": int(status_code),
        "limits": {
            "time_limit_s": float(time_limit_s),
            "memory_limit_mb": int(memory_limit_mb),
            "workers": int(workers),
            "seed": int(seed),
        },
        "model": {
            "variables": 6 * n,
            "one_hot_constraints": n,
            "deduplicated_nonempty_parity_constraints": len(equations),
            "even_parity_encoding": "AddBoolXOr(equation_literals + [fixed_true])",
        },
        "trace": {
            "wall_time_s": solver.WallTime(),
            "num_conflicts": solver.NumConflicts(),
            "num_branches": solver.NumBranches(),
            "response_stats": solver.ResponseStats(),
        },
        "witness": assignment,
        "witness_verified": witness_verified,
        "transformed_css_rank_certificate": transformed_css,
        "linear_relaxation": linear,
        "infeasibility_certificate": (
            {
                "solver": "CP-SAT complete INFEASIBLE proof",
                "independent_linear_xor_contradiction": linear[
                    "contradiction_certificate"
                ],
            }
            if status_code == cp_model.INFEASIBLE
            else None
        ),
    }


def write_full_lc_attempt(attempt: dict[str, Any]) -> Path:
    """Persist an LC attempt without allowing UNKNOWN to overwrite canonical evidence."""

    if attempt["complete"]:
        path = ROOT / "results" / "certificates" / "exp034_target_full_lc_decision.json"
    else:
        path = PARTIAL_LC_PATH
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(attempt, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)
    return path




def build_payload() -> dict[str, Any]:
    row = load_target_catalogue_row()
    code = build_target_code(row)
    H = code.H
    n = code.n

    ranks = css_rank_invariants(H)
    sx = pure_intersection_basis(H, "X")
    sz = pure_intersection_basis(H, "Z")
    components = stabilizer_incidence_components(H)
    independent_components = _independent_incidence_components(H)
    if components != independent_components:
        raise AssertionError("union-find and bipartite BFS incidence components disagree")
    direct_sum = stabilizer_direct_sum_certificate(H)
    if direct_sum["component_sizes"] != [n]:
        raise AssertionError(f"unexpected target direct-sum components: {direct_sum['component_sizes']}")

    # A deterministic nontrivial permutation verifies invariance of the exact
    # ranks; the proof is bijectivity of the coordinate permutation, not this
    # sample alone.
    # The canonical algebra-only run never imports or invokes CP-SAT.
    # A separate explicit CLI attempt may later promote a complete decision.
    permutation = list(range(n - 1, -1, -1))
    permuted_ranks = css_rank_invariants(permute_qubits(H, permutation))
    compared_rank_fields = (
        "rank_s",
        "rank_s_intersect_x_only",
        "rank_s_intersect_z_only",
        "rank_pure_span",
    )
    permutation_check = all(ranks[field] == permuted_ranks[field] for field in compared_rank_fields)
    if not permutation_check:
        raise AssertionError("qubit permutation changed a CSS rank invariant")

    validation = code.validate()
    if not validation["commutes_numpy"] or not validation["commutes_bitset"]:
        raise AssertionError("rebuilt H is not a valid stabilizer matrix")
    if validation["rank_numpy"] != 132 or validation["k"] != 12:
        raise AssertionError("unexpected target parameters")

    full_lc = full_lc_cpsat_formulation(H)
    linear = full_lc_linear_certificate(H)
    full_lc.update({
        "group": (
            "independent choice of all six binary symplectic single-qubit Clifford actions "
            "on each labelled qubit, plus arbitrary stabilizer row operations"
        ),
        "search_space_size": "6^144",
        "linear_relaxation": linear,
    })
    if not linear["linear_relaxation_feasible"]:
        # Every one-hot assignment satisfies the linear consequences, so the
        # XOR contradiction is a complete solver-free refutation.
        full_lc.update({
            "status": "INEQUIVALENT",
            "complete": True,
            "solver_status": "NOT_REQUIRED_LINEAR_XOR_REFUTATION",
            "witness": None,
            "infeasibility_certificate": {
                "independent_linear_xor_contradiction": linear[
                    "contradiction_certificate"
                ],
            },
            "proof": (
                "The listed parity masks and one-hot parity rows XOR to 0=1 over "
                "GF(2). Any per-qubit local Clifford choice satisfies every listed "
                "row, so no independent single-qubit Clifford composed with "
                "stabilizer row operations maps S to CSS. CSS-ness is invariant "
                "under qubit permutations and a permutation conjugate of a local "
                "Clifford is again a local Clifford, so the refutation covers the "
                "full row-operation x permutation x local-Clifford group."
            ),
            "non_vacuousness_control": non_vacuousness_control(),
        })
        decision_path = (
            ROOT / "results" / "certificates" / "exp034_target_full_lc_decision.json"
        )
        if decision_path.exists():
            decision = json.loads(decision_path.read_text())
            if decision.get("status") == "INEQUIVALENT":
                full_lc["solver_corroboration"] = {
                    "solver": "OR-Tools CP-SAT",
                    "solver_status": decision["solver_status"],
                    "limits": decision["limits"],
                    "trace_wall_time_s": decision["trace"]["wall_time_s"],
                    "artifact": str(decision_path.relative_to(ROOT)),
                }
    else:
        full_lc.update({
            "status": "UNRESOLVED",
            "complete": False,
            "solver_status": "NOT_RUN_HOST_LOAD_GUARD",
            "witness": None,
            "infeasibility_certificate": None,
            "residual_uncertainty": (
                "The target may or may not be CSS-equivalent under arbitrary local Clifford transformations. "
                "Local-H inequivalence is not evidence for this larger group."
            ),
            "resume": (
                "Encode one-hot y[j,g] and the persisted projection/nullspace parity constraints in "
                "CP-SAT or XOR-SAT; accept EQUIVALENT only with a directly verified assignment and "
                "INEQUIVALENT only after a complete solver proof."
            ),
        })

    return {
        "experiment": "EXP-034",
        "target": TARGET_ID,
        "scope": "target-only; no catalogue-wide classification",
        "sources": {
            "catalogue_path": str(CATALOGUE_PATH.relative_to(ROOT)),
            "catalogue_sha256": EXPECTED_CATALOGUE_SHA256,
            "catalogue_line": EXPECTED_LINE,
            "raw_target_line_sha256": EXPECTED_RAW_LINE_SHA256,
            "catalogue_bliss_hash": row["bliss_hash"],
            "terms": {field: row[field] for field in EXPECTED_TERMS},
            "construction": "H = [[A,B | C,D], [0,0 | B^T,A^T]]",
            "h_shape": list(H.shape),
            "h_sha256": EXPECTED_H_SHA256,
            "h_hash_encoding": "ASCII '<rows>x<cols>:little-packbits\\n' followed by row-major np.packbits bitorder=little",
            "pure_x_basis_sha256": sha256_bytes(matrix_sha256_payload(sx)),
            "pure_z_basis_sha256": sha256_bytes(matrix_sha256_payload(sz)),
            "context_hashes": {
                "experiments/exp032_catalogue_dedup.py": "2d36c942cef7d2650edd2848a58adddbdc30f8b928db75031670045a149d6110",
                "results/processed/exp032_catalogue_dedup.json": "e4c692cf8e59b7a73f53d897a026079b65762ae5852fedb93ba8823af8a9ba54",
                "third_party/qcode-discovery/evaluation/clifford_equivalence.py": "278d20b02def9500ceb7e1a8f0f56f6517b849f49a4a2ba80ffb844c7b922eeb",
                "third_party/qcode-discovery/CLAUDE.md": "edf76ecd27970c8a7ec9d90075f258148b263ce348991ed2b6c4a61e64b95e49",
                "sources/manifest.yaml": "f62b5017869894a89a8da51f2d7418fddf45a0f118f372f906e716ba2d18d1a2",
                "sources/bibliography.bib": "1ebe0032c4352a1d1edfa6f48f63746f649f49e3764d8d7a619a4dbcebc47623",
            },
        },
        "validation": {
            "rank": _rank_crosscheck(H),
            "commutes_numpy": validation["commutes_numpy"],
            "commutes_bitset": validation["commutes_bitset"],
            "n": validation["n"],
            "k": validation["k"],
            "num_published_rows": validation["num_checks"],
        },
        "css_after_row_operations": {
            "group": "arbitrary invertible stabilizer generator row operations; labelled qubits fixed",
            "status": "INEQUIVALENT",
            "complete": True,
            **ranks,
            "strict_deficit": int(ranks["rank_s"] - ranks["rank_pure_span"]),
        },
        "css_after_qubit_permutation": {
            "group": "arbitrary permutation of the 144 qubits, applied identically to X and Z coordinates, plus row operations",
            "status": "INEQUIVALENT",
            "complete": True,
            "proof": (
                "A qubit permutation is a bijection preserving the X-only and Z-only coordinate subspaces, "
                "so it preserves rank(S), rank(S intersect X-only), and rank(S intersect Z-only)."
            ),
            "deterministic_reverse_permutation_crosscheck": permutation_check,
        },
        "stabilizer_direct_sum_decomposition": {
            "group_decided": (
                "stabilizer row-space decomposition under arbitrary invertible generator row "
                "operations, arbitrary qubit permutations, and arbitrary independent local "
                "Clifford transformations"
            ),
            "status": "INDECOMPOSABLE",
            "complete": True,
            **direct_sum,
            "proof": (
                "The grouped X/Z column-matroid has one circuit-connected qubit component. "
                "For a partition A|B, rank(H_A)+rank(H_B)=rank(H) iff the stabilizer "
                "row space splits into independently supported subspaces on A and B. "
                "Each local Clifford is an invertible change of basis within one qubit's "
                "two-column block, so it preserves every partition rank in this criterion."
            ),
        },
        "incidence_decomposition": {
            "group_decided": (
                "connected components of the supplied bipartite stabilizer-row/qubit incidence representation; "
                "row operations, local Cliffords, and alternative generating sets are excluded"
            ),
            "status": "CONNECTED_NOT_DIRECT_SUM_IN_THIS_REPRESENTATION",
            "complete": True,
            "num_components": len(components),
            "component_sizes": [len(component) for component in components],
            "component_qubits": components,
            "incidences": int(np.count_nonzero(H[:, :n] | H[:, n:])),
            "independent_union_find_bfs_agree": True,
        },
        "css_after_local_hadamards": _local_h_result(H),
        "css_after_full_local_clifford": full_lc,
        "canonical_claims": [
            "row-operation CSS inequivalence",
            "qubit-permutation-plus-row-operation CSS inequivalence",
            "row-space indecomposability under row operations, qubit permutations, and local Cliffords",
            "connectedness of the supplied stabilizer incidence representation",
            "arbitrary local-H-plus-row-operation CSS inequivalence",
        ] + (
            ["arbitrary full local-Clifford CSS inequivalence"]
            if full_lc["complete"] and full_lc["status"] == "INEQUIVALENT"
            else []
        ),
        "noncanonical_unresolved": (
            []
            if full_lc["complete"]
            else ["arbitrary full local-Clifford CSS equivalence"]
        ),
    }


def write_payload(payload: dict[str, Any]) -> None:
    """Publish atomically without downgrading an existing complete LC decision."""

    if OUTPUT_PATH.exists():
        existing = json.loads(OUTPUT_PATH.read_text())
        existing_lc = existing.get("css_after_full_local_clifford", {})
        candidate_lc = payload.get("css_after_full_local_clifford", {})
        if existing_lc.get("complete") and not candidate_lc.get("complete"):
            payload["css_after_full_local_clifford"] = existing_lc
            payload["canonical_claims"] = existing["canonical_claims"]
            payload["noncanonical_unresolved"] = existing["noncanonical_unresolved"]
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    temporary = OUTPUT_PATH.with_suffix(OUTPUT_PATH.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(OUTPUT_PATH)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solve-full-lc", action="store_true")
    parser.add_argument("--time-limit-s", type=float, default=300.0)
    parser.add_argument("--memory-limit-mb", type=int, default=4096)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--seed", type=int, default=20260813)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    if args.solve_full_lc:
        code = build_target_code()
        attempt = solve_full_lc_cpsat(
            code.H,
            time_limit_s=args.time_limit_s,
            memory_limit_mb=args.memory_limit_mb,
            workers=args.workers,
            seed=args.seed,
        )
        path = write_full_lc_attempt(attempt)
        # Only a complete decision may update canonical evidence.  The fresh
        # payload may already carry a complete linear-certificate decision;
        # the solver attempt is only allowed to strengthen, never to replace
        # a linear refutation with a bare solver status.
        if attempt["complete"]:
            payload = build_payload()
            fresh_lc = payload["css_after_full_local_clifford"]
            if fresh_lc["complete"]:
                if fresh_lc["status"] != attempt["status"]:
                    raise AssertionError(
                        "solver decision contradicts the linear certificate"
                    )
            else:
                payload["css_after_full_local_clifford"] = attempt
            claim = (
                "arbitrary full local-Clifford CSS equivalence"
                if attempt["status"] == "EQUIVALENT"
                else "arbitrary full local-Clifford CSS inequivalence"
            )
            if claim not in payload["canonical_claims"]:
                payload["canonical_claims"].append(claim)
            payload["noncanonical_unresolved"] = []
            write_payload(payload)
        print(json.dumps({
            "attempt": str(path.relative_to(ROOT)),
            "status": attempt["status"],
            "solver_status": attempt["solver_status"],
            "canonical_updated": attempt["complete"],
        }, sort_keys=True))
        return

    payload = build_payload()
    write_payload(payload)
    print(json.dumps({
        "output": str(OUTPUT_PATH.relative_to(ROOT)),
        "ranks": payload["css_after_row_operations"]["independent_checks"]["numpy"],
        "components": payload["incidence_decomposition"]["component_sizes"],
        "local_h": payload["css_after_local_hadamards"]["status"],
        "full_lc": payload["css_after_full_local_clifford"]["status"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
