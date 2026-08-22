"""Focused exact checks for EXP-034 target equivalence evidence."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from itertools import product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location(
    "exp034_target_equivalence", ROOT / "experiments" / "exp034_target_equivalence.py"
)
_EXP = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _EXP
assert _SPEC.loader is not None
_SPEC.loader.exec_module(_EXP)

from qec_research.equivalence.css import (  # noqa: E402
    LOCAL_CLIFFORD_PROJECTIONS,
    css_rank_invariants,
    full_lc_parity_equations,
    full_lc_linear_certificate,
    one_hot_parity_row,
    full_lc_parity_sha256,
    full_lc_cpsat_formulation,
    local_h_affine_system,
    matrix_sha256_payload,
    permute_qubits,
    projection_invariant,
    solve_affine_system,
    stabilizer_direct_sum_certificate,
    stabilizer_direct_sum_components,
    stabilizer_incidence_components,
)


def _paulis(*rows: str) -> np.ndarray:
    n = len(rows[0])
    H = np.zeros((len(rows), 2 * n), dtype=np.uint8)
    for row_index, row in enumerate(rows):
        assert len(row) == n
        for qubit, pauli in enumerate(row):
            if pauli in "XY":
                H[row_index, qubit] = 1
            if pauli in "ZY":
                H[row_index, n + qubit] = 1
    return H


def _full_lc_exhaustive(H: np.ndarray):
    n = H.shape[1] // 2
    for choices in product(range(6), repeat=n):
        if projection_invariant(H, choices, full_lc=True):
            return choices
    return None


def test_rebuilds_literal_target_H_and_hashes():
    row = _EXP.load_target_catalogue_row()
    assert row["code_id"] == "12_6_0193"
    assert row["A_terms"] == [[1, 2], [10, 3], [10, 4]]
    assert row["B_terms"] == [[0, 0], [1, 5], [11, 4]]
    assert row["C_terms"] == [[1, 3], [10, 3]]
    assert row["D_terms"] == [[1, 5], [10, 5]]
    H = _EXP.build_target_code(row).H
    assert H.shape == (144, 288)
    assert hashlib.sha256(matrix_sha256_payload(H)).hexdigest() == _EXP.EXPECTED_H_SHA256


def test_exact_row_operation_css_ranks_and_permutation_invariance():
    H = _EXP.build_target_code().H
    ranks = css_rank_invariants(H)
    assert ranks["independent_checks"]["agree"]
    assert ranks["rank_s"] == 132
    assert ranks["rank_s_intersect_x_only"] == 20
    assert ranks["rank_s_intersect_z_only"] == 66
    assert ranks["rank_pure_span"] == 86
    assert not ranks["css_after_row_operations"]

    permuted = permute_qubits(H, list(range(1, 144)) + [0])
    permuted_ranks = css_rank_invariants(permuted)
    assert [permuted_ranks[key] for key in (
        "rank_s", "rank_s_intersect_x_only", "rank_s_intersect_z_only"
    )] == [132, 20, 66]


def test_target_incidence_is_one_component_under_declared_representation():
    H = _EXP.build_target_code().H
    components = stabilizer_incidence_components(H)
    assert components == [list(range(144))]
    assert _EXP._independent_incidence_components(H) == components


def test_target_rowspace_is_exactly_indivisible_under_row_ops_permutations_and_lc():
    H = _EXP.build_target_code().H
    certificate = stabilizer_direct_sum_certificate(H)
    assert certificate["components"] == [list(range(144))]
    assert certificate["component_sizes"] == [144]
    assert certificate["num_components"] == 1
    assert not certificate["is_direct_sum"]
    assert certificate["rank_s_numpy"] == certificate["rank_s_bitset"] == 132
    assert certificate["sum_component_ranks"] == 132
    assert certificate["reconstructed_rowspace_rank_numpy"] == 132
    assert certificate["reconstructed_rowspace_rank_bitset"] == 132

    permuted = permute_qubits(H, list(range(143, -1, -1)))
    assert stabilizer_direct_sum_components(permuted) == [list(range(144))]

    local_h = H.copy()
    for qubit in range(0, 144, 2):
        local_h[:, [qubit, 144 + qubit]] = local_h[
            :, [144 + qubit, qubit]
        ]
    assert stabilizer_direct_sum_components(local_h) == [list(range(144))]
    local_lc = local_h.copy()
    for qubit in range(1, 144, 3):
        local_lc[:, 144 + qubit] ^= local_lc[:, qubit]
    assert stabilizer_direct_sum_components(local_lc) == [list(range(144))]
    declared = _EXP.build_payload()["stabilizer_direct_sum_decomposition"]
    assert "local Clifford transformations" in declared["group_decided"]


def test_direct_sum_controls_include_misleading_connected_generators():
    direct = _paulis("XXII", "ZZII", "IIXX", "IIZZ")
    direct_certificate = stabilizer_direct_sum_certificate(direct)
    assert direct_certificate["component_sizes"] == [2, 2]
    assert direct_certificate["is_direct_sum"]
    assert direct_certificate["sum_component_ranks"] == 4

    # These supplied rows connect both qubits through XX, but their row space
    # is span{XI, IX}; arbitrary row reduction reveals the exact direct sum.
    misleading = _paulis("XX", "IX")
    assert stabilizer_incidence_components(misleading) == [[0, 1]]
    misleading_certificate = stabilizer_direct_sum_certificate(misleading)
    assert misleading_certificate["components"] == [[0], [1]]
    assert misleading_certificate["component_sizes"] == [1, 1]
    assert misleading_certificate["sum_component_ranks"] == 2

    row_reduced = misleading.copy()
    row_reduced[0] ^= row_reduced[1]
    assert stabilizer_incidence_components(row_reduced) == [[0], [1]]
    assert stabilizer_direct_sum_components(row_reduced) == [[0], [1]]


def test_local_h_positive_and_negative_toy_controls():
    # Positive but not native CSS: H on qubit 1 maps XX/ZZ Bell stabilizers to XZ/ZX.
    positive = _paulis("XZ", "ZX")
    assert not css_rank_invariants(positive)["css_after_row_operations"]
    equations, nvars = local_h_affine_system(positive)
    solved = solve_affine_system(equations, nvars)
    assert solved.feasible
    assert solved.assignment is not None
    assert projection_invariant(positive, solved.assignment, full_lc=False)

    # Negative: a single Y stabilizer is invariant under H and cannot become an X/Z axis.
    negative = _paulis("Y")
    equations, nvars = local_h_affine_system(negative)
    solved = solve_affine_system(equations, nvars)
    assert not solved.feasible
    assert solved.contradiction is not None
    xor_mask = 0
    xor_rhs = 0
    for mask, rhs in solved.contradiction:
        xor_mask ^= mask
        xor_rhs ^= rhs
    assert (xor_mask, xor_rhs) == (0, 1)


def test_full_lc_formulation_matches_exhaustive_toy_controls():
    assert len(LOCAL_CLIFFORD_PROJECTIONS) == 6
    positive = _paulis("Y")
    witness = _full_lc_exhaustive(positive)
    assert witness is not None
    assert projection_invariant(positive, witness, full_lc=True)

    # [[5,1,3]] perfect-code stabilizers are not LC-equivalent to CSS.  Six
    # choices on five qubits is only 7,776 assignments, so this is an exact
    # negative toy control for the same projection-invariance formulation.
    negative = _paulis("XZZXI", "IXZZX", "XIXZZ", "ZXIXZ")
    assert _full_lc_exhaustive(negative) is None
    formulation = full_lc_cpsat_formulation(negative)
    assert formulation["variables"] == 30
    assert formulation["one_hot_constraints"] == 5


def test_full_lc_cpsat_matches_exhaustive_toys_and_routes_unknown(monkeypatch):
    positive = _paulis("Y")
    exhaustive_positive = _full_lc_exhaustive(positive)
    assert exhaustive_positive is not None
    solved_positive = _EXP.solve_full_lc_cpsat(
        positive, time_limit_s=5.0, memory_limit_mb=256
    )
    assert solved_positive["status"] == "EQUIVALENT"
    assert solved_positive["complete"]
    assert solved_positive["witness_verified"]
    transformed = solved_positive["transformed_css_rank_certificate"]
    assert transformed["css_after_row_operations"]
    assert transformed["rank_pure_span"] == transformed["rank_s"]

    negative = _paulis("XZZXI", "IXZZX", "XIXZZ", "ZXIXZ")
    assert _full_lc_exhaustive(negative) is None
    solved_negative = _EXP.solve_full_lc_cpsat(
        negative, time_limit_s=5.0, memory_limit_mb=256
    )
    assert solved_negative["status"] == "INEQUIVALENT"
    assert solved_negative["complete"]
    assert solved_negative["solver_status"] == "INFEASIBLE"

    equations, n = full_lc_parity_equations(negative)
    formulation = full_lc_cpsat_formulation(negative)
    assert formulation["deduplicated_nonempty_parity_constraints"] == len(equations)
    assert n == 5
    assert len(equations) == len(set(equations))
    assert formulation["parity_equations_sha256"] == full_lc_parity_sha256(equations, n)
    assert "parity_equation_masks_hex" not in formulation

    class UnknownSolver:
        def __init__(self):
            class Parameters:
                pass
            self.parameters = Parameters()

        def Solve(self, _model):
            from ortools.sat.python import cp_model
            return cp_model.UNKNOWN

        def StatusName(self, _status):
            return "UNKNOWN"

        def WallTime(self):
            return 0.0

        def NumConflicts(self):
            return 0

        def NumBranches(self):
            return 0

        def ResponseStats(self):
            return "forced timeout control"

    from ortools.sat.python import cp_model
    monkeypatch.setattr(cp_model, "CpSolver", UnknownSolver)
    unknown = _EXP.solve_full_lc_cpsat(
        positive, time_limit_s=0.001, memory_limit_mb=64
    )
    assert unknown["status"] == "UNRESOLVED"
    assert not unknown["complete"]
    assert unknown["witness"] is None
    assert unknown["infeasibility_certificate"] is None


def test_target_local_h_certificate_and_full_lc_status_are_honest():
    payload = _EXP.build_payload()
    local_h = payload["css_after_local_hadamards"]
    assert local_h["complete"]
    assert local_h["status"] == "INEQUIVALENT"
    assert local_h["coefficient_rank"] < local_h["augmented_rank"]
    certificate = local_h["infeasibility_certificate"]
    assert certificate["xor_coefficient_mask_hex"] == "0x0"
    assert certificate["xor_rhs"] == 1

    full_lc = payload["css_after_full_local_clifford"]
    assert full_lc["status"] == "INEQUIVALENT"
    assert full_lc["complete"]
    assert full_lc["witness"] is None
    assert full_lc["solver_status"] == "NOT_REQUIRED_LINEAR_XOR_REFUTATION"
    assert "residual_uncertainty" not in full_lc
    assert (
        "arbitrary full local-Clifford CSS inequivalence"
        in payload["canonical_claims"]
    )
    assert payload["noncanonical_unresolved"] == []


def test_target_full_lc_linear_certificate_is_machine_checked():
    H = _EXP.build_target_code().H
    linear = full_lc_linear_certificate(H)
    assert not linear["linear_relaxation_feasible"]
    assert linear["coefficient_rank"] + 1 == linear["augmented_rank"]
    stored = linear["contradiction_certificate"]
    masks, n = full_lc_parity_equations(H)
    mask_set = set(masks)
    xor_mask = 0
    xor_rhs = 0
    one_hot_rows = 0
    for row in stored["equations"]:
        mask = int(row["mask_hex"], 16)
        xor_mask ^= mask
        xor_rhs ^= row["rhs"]
        if row["rhs"] == 0:
            assert mask in mask_set
        else:
            one_hot_rows += 1
            assert mask == one_hot_parity_row((mask.bit_length() - 1) // 6)
    assert (xor_mask, xor_rhs) == (0, 1)
    assert one_hot_rows == stored["num_one_hot_parity_rows"] == 1
    assert stored["num_parity_rows"] == len(stored["equations"]) - 1

    # Controls: a locally fixable code must stay linear-feasible, and the
    # [[5,1,3]] perfect code shows the relaxation is not vacuous - it stays
    # linear-feasible although full LC equivalence is exhaustively refuted.
    assert full_lc_linear_certificate(_paulis("Y"))["linear_relaxation_feasible"]
    negative = _paulis("XZZXI", "IXZZX", "XIXZZ", "ZXIXZ")
    assert full_lc_linear_certificate(negative)["linear_relaxation_feasible"]

    persisted = json.loads(
        (ROOT / "results" / "processed" / "exp034_target_equivalence.json").read_text()
    )
    persisted_lc = persisted["css_after_full_local_clifford"]
    assert persisted_lc["status"] == "INEQUIVALENT"
    assert (
        persisted_lc["infeasibility_certificate"][
            "independent_linear_xor_contradiction"
        ]
        == stored
    )


def test_persisted_certificate_recomputes_exact_payload():
    output = ROOT / "results" / "processed" / "exp034_target_equivalence.json"
    persisted = json.loads(output.read_text())
    rebuilt = _EXP.build_payload()
    assert persisted == rebuilt
