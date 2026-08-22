"""Exact CSS and local-Clifford equivalence helpers for stabilizer row spaces.

A check matrix uses row vectors ``(x | z)``.  All ranks and containments in
this module are over GF(2), and every CSS decision is a property of the full
stabilizer row space rather than of the supplied generating rows.
"""

from __future__ import annotations
import hashlib

from dataclasses import dataclass
from itertools import product
from typing import Iterable, Sequence

import numpy as np

from ..gf2.linalg import (
    nullspace_np,
    rank_bitset,
    rank_np,
    rows_to_bitsets,
    rref_np,
)

# The six rank-one projections C P_X C^{-1}, one for each one-qubit
# symplectic action C in Sp(2, 2) = GL(2, 2).  Distinct Clifford phases have
# no action on the binary symplectic representation.
LOCAL_CLIFFORD_PROJECTIONS: tuple[np.ndarray, ...] = tuple(
    matrix
    for entries in product((0, 1), repeat=4)
    if (matrix := np.asarray(entries, dtype=np.uint8).reshape(2, 2)).tolist()
    and rank_np(matrix) == 1
    and np.array_equal((matrix @ matrix) & 1, matrix)
)

_P_X = np.asarray([[1, 0], [0, 0]], dtype=np.uint8)
_I2 = np.eye(2, dtype=np.uint8)
_GL2 = tuple(
    matrix
    for entries in product((0, 1), repeat=4)
    if rank_np(matrix := np.asarray(entries, dtype=np.uint8).reshape(2, 2)) == 2
)


def _inverse_2x2(matrix: np.ndarray) -> np.ndarray:
    a, b, c, d = (int(value) for value in matrix.flat)
    inverse = np.asarray([[d, b], [c, a]], dtype=np.uint8)
    if not np.array_equal((matrix @ inverse) & 1, _I2):
        raise AssertionError("matrix is not invertible over GF(2)")
    return inverse


LOCAL_CLIFFORD_ACTIONS: tuple[np.ndarray, ...] = tuple(
    next(
        action
        for action in _GL2
        if np.array_equal(
            (action @ _P_X @ _inverse_2x2(action)) & 1,
            projection,
        )
    )
    for projection in LOCAL_CLIFFORD_PROJECTIONS
)

LOCAL_H_PROJECTIONS: tuple[np.ndarray, ...] = (
    np.asarray([[1, 0], [0, 0]], dtype=np.uint8),
    np.asarray([[0, 0], [0, 1]], dtype=np.uint8),
)


def _check_matrix(H: np.ndarray) -> np.ndarray:
    matrix = np.asarray(H, dtype=np.uint8) & 1
    if matrix.ndim != 2 or matrix.shape[1] % 2:
        raise ValueError(f"expected a two-dimensional symplectic matrix, got {matrix.shape}")
    return matrix


def matrix_sha256_payload(H: np.ndarray) -> bytes:
    """Return the unambiguous payload used to hash a binary matrix certificate."""

    matrix = _check_matrix(H)
    header = f"{matrix.shape[0]}x{matrix.shape[1]}:little-packbits\n".encode()
    return header + np.packbits(matrix, axis=None, bitorder="little").tobytes()


def pure_intersection_basis(H: np.ndarray, kind: str) -> np.ndarray:
    """Return an RREF basis for ``rowspan(H)`` intersected with one Pauli axis.

    ``kind='X'`` imposes a zero Z half and ``kind='Z'`` imposes a zero X
    half.  Coefficients are obtained from the left kernel of the forbidden
    half; a final RREF removes dependencies caused by redundant input rows.
    """

    matrix = _check_matrix(H)
    n = matrix.shape[1] // 2
    if kind == "X":
        forbidden = matrix[:, n:]
    elif kind == "Z":
        forbidden = matrix[:, :n]
    else:
        raise ValueError("kind must be 'X' or 'Z'")
    coefficients = nullspace_np(forbidden.T)
    physical = (coefficients.astype(np.int64) @ matrix.astype(np.int64) % 2).astype(np.uint8)
    return rref_np(physical)[0]


def css_rank_invariants(H: np.ndarray) -> dict[str, object]:
    """Compute independent exact ranks and the row-operation CSS criterion.

    Let ``S=rowspan(H)``, ``S_X=S intersect (F_2^n|0)`` and
    ``S_Z=S intersect (0|F_2^n)``.  Their intersection is zero.  Therefore a
    pure-X/pure-Z basis of S exists iff ``S=S_X direct_sum S_Z``, equivalently
    iff ``dim(S_X)+dim(S_Z)=dim(S)``.  Projection rank-nullity also gives
    ``dim(S_X)=dim(S)-rank(H_Z)`` and the analogous Z formula.
    """

    matrix = _check_matrix(H)
    n = matrix.shape[1] // 2
    sx = pure_intersection_basis(matrix, "X")
    sz = pure_intersection_basis(matrix, "Z")

    rank_s_np = rank_np(matrix)
    rank_x_np = rank_np(sx)
    rank_z_np = rank_np(sz)
    rank_pure_span_np = rank_np(np.vstack([sx, sz]))

    rank_s_bitset = rank_bitset(rows_to_bitsets(matrix), 2 * n)
    rank_x_projection_bitset = rank_bitset(rows_to_bitsets(matrix[:, :n]), n)
    rank_z_projection_bitset = rank_bitset(rows_to_bitsets(matrix[:, n:]), n)
    rank_x_bitset = rank_s_bitset - rank_z_projection_bitset
    rank_z_bitset = rank_s_bitset - rank_x_projection_bitset

    if (rank_s_np, rank_x_np, rank_z_np) != (
        rank_s_bitset,
        rank_x_bitset,
        rank_z_bitset,
    ):
        raise AssertionError("independent NumPy and bitset GF(2) ranks disagree")
    if rank_pure_span_np != rank_x_np + rank_z_np:
        raise AssertionError("pure-X and pure-Z intersections must have trivial intersection")

    return {
        "rank_s": rank_s_np,
        "rank_s_intersect_x_only": rank_x_np,
        "rank_s_intersect_z_only": rank_z_np,
        "rank_pure_span": rank_pure_span_np,
        "rank_x_projection": rank_x_projection_bitset,
        "rank_z_projection": rank_z_projection_bitset,
        "css_after_row_operations": rank_x_np + rank_z_np == rank_s_np,
        "independent_checks": {
            "numpy": [rank_s_np, rank_x_np, rank_z_np],
            "bitset_projection_rank_nullity": [
                rank_s_bitset,
                rank_x_bitset,
                rank_z_bitset,
            ],
            "agree": True,
        },
        "proof": (
            "S_X and S_Z intersect only in zero. A pure-X/pure-Z stabilizer "
            "basis exists iff S is their direct sum, iff dim(S_X)+dim(S_Z)=dim(S)."
        ),
    }


def stabilizer_incidence_components(H: np.ndarray) -> list[list[int]]:
    """Partition qubits in the bipartite check/qubit incidence graph.

    Row operations are deliberately not included in this transformation
    group: this decides decomposability of the supplied stabilizer incidence
    representation.  Isolated qubits appear as singleton components.
    """

    matrix = _check_matrix(H)
    n = matrix.shape[1] // 2
    support = matrix[:, :n] | matrix[:, n:]
    parent = list(range(n))

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(left: int, right: int) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[left_root] = right_root

    for row in support:
        qubits = np.flatnonzero(row)
        if qubits.size:
            first = int(qubits[0])
            for qubit in qubits[1:]:
                union(first, int(qubit))

    groups: dict[int, list[int]] = {}
    for qubit in range(n):
        groups.setdefault(find(qubit), []).append(qubit)
    return sorted(groups.values(), key=lambda group: (-len(group), group))


def stabilizer_direct_sum_components(H: np.ndarray) -> list[list[int]]:
    """Return the finest qubit direct-sum partition of the stabilizer row space.

    For a qubit partition ``E = A disjoint_union B``, the stabilizer row
    space splits as ``S_A direct_sum S_B`` exactly when
    ``rank(H_A) + rank(H_B) = rank(H)``.  The finest additive partition of a
    represented linear matroid is its circuit-connectivity partition.  We
    construct it from all fundamental circuits relative to one column basis,
    then union the X and Z columns of each physical qubit because those two
    coordinates are not allowed to separate.

    Left multiplication by arbitrary invertible row operations preserves all
    column dependencies; a qubit permutation merely permutes the paired
    elements.  Thus this partition is exact under row operations and qubit
    permutations, unlike connectivity of the supplied check-incidence graph.
    """

    matrix = _check_matrix(H)
    n = matrix.shape[1] // 2
    _, basis_columns = rref_np(matrix)
    basis_set = set(basis_columns)
    parent = list(range(2 * n))

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(left: int, right: int) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[left_root] = right_root

    # Every nonbasis element has one unique fundamental circuit with respect
    # to the selected basis.  These circuits generate matroid connectivity.
    for column in range(2 * n):
        if column in basis_set:
            continue
        circuit_columns = basis_columns + [column]
        dependencies = nullspace_np(matrix[:, circuit_columns])
        if dependencies.shape != (1, len(circuit_columns)):
            raise AssertionError("a fundamental circuit must have a unique dependency")
        dependency = dependencies[0]
        if not dependency[-1]:
            raise AssertionError("fundamental dependency omitted its nonbasis column")
        support = [
            circuit_columns[index]
            for index in np.flatnonzero(dependency)
        ]
        for connected_column in support[1:]:
            union(support[0], connected_column)

    # A decomposition is over qubits, never over individual X/Z coordinates.
    for qubit in range(n):
        union(qubit, n + qubit)

    groups: dict[int, list[int]] = {}
    for qubit in range(n):
        groups.setdefault(find(qubit), []).append(qubit)
    return sorted(groups.values(), key=lambda group: (-len(group), group))


def stabilizer_direct_sum_certificate(H: np.ndarray) -> dict[str, object]:
    """Certify the finest row-space direct-sum partition with independent ranks."""

    matrix = _check_matrix(H)
    n = matrix.shape[1] // 2
    components = stabilizer_direct_sum_components(matrix)
    full_rank_np = rank_np(matrix)
    full_rank_bitset = rank_bitset(rows_to_bitsets(matrix), 2 * n)
    if full_rank_np != full_rank_bitset:
        raise AssertionError("full NumPy and bitset ranks disagree")

    component_records: list[dict[str, object]] = []
    reconstructed_bases: list[np.ndarray] = []
    component_rank_sum = 0
    all_columns = set(range(2 * n))
    for component in components:
        columns = sorted(component + [n + qubit for qubit in component])
        complement = sorted(all_columns.difference(columns))
        restricted = matrix[:, columns]
        restricted_rank_np = rank_np(restricted)
        restricted_rank_bitset = rank_bitset(
            rows_to_bitsets(restricted), len(columns)
        )
        if restricted_rank_np != restricted_rank_bitset:
            raise AssertionError("component NumPy and bitset ranks disagree")

        # Coefficients annihilating the complement reconstruct exactly the
        # subspace of stabilizers supported inside this component.
        coefficients = nullspace_np(matrix[:, complement].T)
        supported = (
            coefficients.astype(np.int64) @ matrix.astype(np.int64) % 2
        ).astype(np.uint8)
        supported_basis = rref_np(supported)[0]
        if complement and supported_basis[:, complement].any():
            raise AssertionError("reconstructed component basis leaks outside its block")
        supported_rank_np = rank_np(supported_basis[:, columns])
        supported_rank_bitset = rank_bitset(
            rows_to_bitsets(supported_basis[:, columns]), len(columns)
        )
        if supported_rank_np != supported_rank_bitset:
            raise AssertionError("reconstructed component ranks disagree")
        if supported_rank_np != restricted_rank_np:
            raise AssertionError("component projection and supported-subspace ranks disagree")

        reconstructed_bases.append(supported_basis)
        component_rank_sum += restricted_rank_np
        component_records.append({
            "qubits": component,
            "num_qubits": len(component),
            "restricted_rank_numpy": restricted_rank_np,
            "restricted_rank_bitset": restricted_rank_bitset,
            "reconstructed_supported_rank_numpy": supported_rank_np,
            "reconstructed_supported_rank_bitset": supported_rank_bitset,
        })

    reconstructed = (
        np.vstack(reconstructed_bases)
        if reconstructed_bases
        else np.zeros((0, 2 * n), dtype=np.uint8)
    )
    reconstructed_rank_np = rank_np(reconstructed)
    reconstructed_rank_bitset = rank_bitset(
        rows_to_bitsets(reconstructed), 2 * n
    )
    joint_rank_np = rank_np(np.vstack([matrix, reconstructed]))
    if not (
        component_rank_sum
        == reconstructed_rank_np
        == reconstructed_rank_bitset
        == joint_rank_np
        == full_rank_np
    ):
        raise AssertionError("component bases do not reconstruct the stabilizer row space")

    return {
        "components": components,
        "component_sizes": [len(component) for component in components],
        "num_components": len(components),
        "is_direct_sum": len(components) > 1,
        "rank_s_numpy": full_rank_np,
        "rank_s_bitset": full_rank_bitset,
        "sum_component_ranks": component_rank_sum,
        "reconstructed_rowspace_rank_numpy": reconstructed_rank_np,
        "reconstructed_rowspace_rank_bitset": reconstructed_rank_bitset,
        "joint_original_reconstructed_rank_numpy": joint_rank_np,
        "component_certificates": component_records,
        "criterion": (
            "For every qubit partition A|B, S splits iff "
            "rank(H_A)+rank(H_B)=rank(H). Grouped column-matroid circuit "
            "components give the unique finest such partition."
        ),
    }


def permute_qubits(H: np.ndarray, new_to_old: Sequence[int]) -> np.ndarray:
    """Apply one qubit permutation to both symplectic halves."""

    matrix = _check_matrix(H)
    n = matrix.shape[1] // 2
    permutation = np.asarray(new_to_old, dtype=int)
    if permutation.shape != (n,) or sorted(permutation.tolist()) != list(range(n)):
        raise ValueError("new_to_old must be a permutation of all qubits")
    return np.hstack([matrix[:, :n][:, permutation], matrix[:, n:][:, permutation]])


def projection_invariant(H: np.ndarray, choices: Sequence[int], *, full_lc: bool) -> bool:
    """Verify a local-H or local-Clifford CSS witness by row-space containment."""

    matrix = _check_matrix(H)
    n = matrix.shape[1] // 2
    projections = LOCAL_CLIFFORD_PROJECTIONS if full_lc else LOCAL_H_PROJECTIONS
    if len(choices) != n or any(choice < 0 or choice >= len(projections) for choice in choices):
        raise ValueError("invalid local projection assignment")
    projected = np.zeros_like(matrix)
    for qubit, choice in enumerate(choices):
        columns = matrix[:, [qubit, n + qubit]]
        projected[:, [qubit, n + qubit]] = (
            columns.astype(np.int64) @ projections[choice].astype(np.int64) % 2
        ).astype(np.uint8)
    rank_s = rank_np(matrix)
    return rank_np(np.vstack([matrix, projected])) == rank_s


def local_h_affine_system(H: np.ndarray) -> tuple[list[tuple[int, int]], int]:
    """Return the exact arbitrary-per-qubit local-H CSS feasibility system.

    A CSS row space is invariant under the global X-axis projection.  After
    local Hadamards, conjugating that projection back gives, independently on
    qubit j, P_X (choice 0) or P_Z (choice 1).  Thus local-H-to-CSS is
    equivalent to ``S P(s) subseteq S``.  Testing a basis of S against a basis
    of its ordinary orthogonal complement gives affine equations in the n
    Boolean choices.  Each equation is returned as ``(coefficient_bitmask,
    rhs)``.  This is finite, complete at stabilizer-group level, and contains
    no generator-level assumption.
    """

    matrix = _check_matrix(H)
    n = matrix.shape[1] // 2
    stabilizer_basis = rref_np(matrix)[0]
    orthogonal_basis = nullspace_np(matrix)
    equations: set[tuple[int, int]] = set()
    for stabilizer in stabilizer_basis:
        for orthogonal in orthogonal_basis:
            mask = 0
            rhs = 0
            for qubit in range(n):
                p_x = int(stabilizer[qubit] & orthogonal[qubit])
                p_z = int(stabilizer[n + qubit] & orthogonal[n + qubit])
                rhs ^= p_x
                if p_x ^ p_z:
                    mask |= 1 << qubit
            if mask or rhs:
                equations.add((mask, rhs))
    return sorted(equations), n


@dataclass(frozen=True)
class AffineSolveResult:
    feasible: bool
    rank: int
    augmented_rank: int
    assignment: tuple[int, ...] | None
    contradiction: tuple[tuple[int, int], ...] | None


def solve_affine_system(equations: Iterable[tuple[int, int]], nvars: int) -> AffineSolveResult:
    """Solve affine GF(2) equations and retain an XOR contradiction certificate."""

    equation_list = list(equations)
    coefficient_rank = rank_bitset([mask for mask, _ in equation_list], nvars)
    augmented_rank = rank_bitset(
        [mask | (rhs << nvars) for mask, rhs in equation_list],
        nvars + 1,
    )
    mask_all = (1 << nvars) - 1
    basis: dict[int, int] = {}
    provenance: dict[int, set[int]] = {}
    contradiction: set[int] | None = None

    for index, (mask, rhs) in enumerate(equation_list):
        if mask & ~mask_all or rhs not in (0, 1):
            raise ValueError("invalid affine equation")
        current = mask | (rhs << nvars)
        source = {index}
        while True:
            variables = current & mask_all
            if not variables:
                if (current >> nvars) & 1:
                    contradiction = source
                break
            pivot = (variables & -variables).bit_length() - 1
            if pivot in basis:
                current ^= basis[pivot]
                source ^= provenance[pivot]
            else:
                basis[pivot] = current
                provenance[pivot] = source
                break
        if contradiction is not None:
            break

    if contradiction is not None:
        certificate = tuple(equation_list[index] for index in sorted(contradiction))
        if _xor_equations(certificate) != (0, 1):
            raise AssertionError("invalid affine contradiction provenance")
        return AffineSolveResult(False, coefficient_rank, augmented_rank, None, certificate)

    assignment_bits = 0
    for pivot in sorted(basis, reverse=True):
        row = basis[pivot]
        rhs = (row >> nvars) & 1
        other = (row & mask_all) & ~(1 << pivot)
        if rhs ^ ((other & assignment_bits).bit_count() & 1):
            assignment_bits |= 1 << pivot
    assignment = tuple((assignment_bits >> index) & 1 for index in range(nvars))
    if coefficient_rank != augmented_rank:
        raise AssertionError("rank test found an affine contradiction without provenance")
    return AffineSolveResult(True, coefficient_rank, augmented_rank, assignment, None)


def _xor_equations(equations: Iterable[tuple[int, int]]) -> tuple[int, int]:
    mask = 0
    rhs = 0
    for equation_mask, equation_rhs in equations:
        mask ^= equation_mask
        rhs ^= equation_rhs
    return mask, rhs


def full_lc_parity_equations(H: np.ndarray) -> tuple[list[int], int]:
    """Return deduplicated homogeneous XOR equations for full local LC.

    Bit ``6*j+g`` selects local projection ``g`` on qubit ``j``.  Together
    with one-hot constraints, a mask requires the XOR of its selected literals
    to equal zero.  Empty equations are tautologies and are omitted.
    """

    matrix = _check_matrix(H)
    n = matrix.shape[1] // 2
    stabilizer_basis = rref_np(matrix)[0]
    orthogonal_basis = nullspace_np(matrix)
    equations: set[int] = set()
    for stabilizer in stabilizer_basis:
        for orthogonal in orthogonal_basis:
            mask = 0
            for qubit in range(n):
                stabilizer_local = stabilizer[[qubit, n + qubit]]
                orthogonal_local = orthogonal[[qubit, n + qubit]]
                if not stabilizer_local.any() or not orthogonal_local.any():
                    continue
                for gate, projection in enumerate(LOCAL_CLIFFORD_PROJECTIONS):
                    parity = int(
                        stabilizer_local.astype(np.int64)
                        @ projection.astype(np.int64)
                        @ orthogonal_local.astype(np.int64)
                        % 2
                    )
                    if parity:
                        mask |= 1 << (6 * qubit + gate)
            if mask:
                equations.add(mask)
    return sorted(equations), n


def full_lc_parity_sha256(equations: Sequence[int], nqubits: int) -> str:
    """Hash sorted XOR masks with an unambiguous fixed-width encoding."""

    width = (6 * nqubits + 7) // 8
    digest = hashlib.sha256()
    digest.update(f"full-lc-xor-v1;n={nqubits};width={width};count={len(equations)}\\n".encode())
    for mask in equations:
        digest.update(int(mask).to_bytes(width, "little"))
    return digest.hexdigest()


def one_hot_parity_row(qubit: int) -> int:
    """Linear consequence of exactly-one: the six choices of a qubit XOR to 1."""

    return 0x3F << (6 * qubit)


def full_lc_linear_certificate(H: np.ndarray) -> dict[str, object]:
    """Exact GF(2) feasibility of the full-LC parity system's linear relaxation.

    Every one-hot assignment satisfies both the homogeneous parity masks and
    the per-qubit parity ``XOR_g y[j,g] = 1``.  An XOR contradiction over
    these equations therefore refutes CSS equivalence under arbitrary
    independent single-qubit Cliffords outright, with a certificate
    checkable by re-XOR-ing the listed rows.
    """

    matrix = _check_matrix(H)
    masks, n = full_lc_parity_equations(matrix)
    nvars = 6 * n
    equations = [(mask, 0) for mask in masks]
    equations += [(one_hot_parity_row(qubit), 1) for qubit in range(n)]
    solved = solve_affine_system(equations, nvars)
    contradiction = None
    if solved.contradiction is not None:
        xor_mask, xor_rhs = _xor_equations(solved.contradiction)
        if (xor_mask, xor_rhs) != (0, 1):
            raise AssertionError("stored contradiction does not XOR to 0=1")
        mask_set = set(masks)
        one_hot_rows = {one_hot_parity_row(qubit) for qubit in range(n)}
        for mask, rhs in solved.contradiction:
            if rhs == 0 and mask not in mask_set:
                raise AssertionError("contradiction uses an unknown parity mask")
            if rhs == 1 and mask not in one_hot_rows:
                raise AssertionError("contradiction uses an invalid one-hot row")
        contradiction = {
            "equations": [
                {"mask_hex": hex(mask), "rhs": rhs}
                for mask, rhs in solved.contradiction
            ],
            "num_parity_rows": sum(
                1 for _, rhs in solved.contradiction if rhs == 0
            ),
            "num_one_hot_parity_rows": sum(
                1 for _, rhs in solved.contradiction if rhs == 1
            ),
            "xor_mask_hex": hex(xor_mask),
            "xor_rhs": xor_rhs,
        }
    return {
        "method": (
            "exact affine GF(2) feasibility of the deduplicated parity masks "
            "plus the per-qubit one-hot parity consequences"
        ),
        "variables": nvars,
        "num_equations": len(equations),
        "coefficient_rank": solved.rank,
        "augmented_rank": solved.augmented_rank,
        "linear_relaxation_feasible": solved.feasible,
        "contradiction_certificate": contradiction,
        "implication": (
            "every one-hot assignment satisfies these linear consequences, so "
            "an XOR contradiction proves no independent single-qubit Clifford "
            "maps S to CSS under any stabilizer row operations"
        ),
    }


def full_lc_cpsat_formulation(H: np.ndarray) -> dict[str, object]:
    """Describe a complete finite CP-SAT formulation for arbitrary local LC.

    Variables ``y[j,g]`` select exactly one of the six local rank-one
    projections.  For every stabilizer-basis row h and every ordinary-nullspace
    row q, invariance imposes
    ``XOR_{j,g: h_j P_g q_j^T=1} y[j,g] = 0``.  Feasibility is exactly
    equivalence to CSS under independent single-qubit Cliffords; infeasibility
    needs a completed SAT/CP-SAT proof, not merely this formulation.
    """

    matrix = _check_matrix(H)
    equations, n = full_lc_parity_equations(matrix)
    return {
        "variables": 6 * n,
        "one_hot_constraints": n,
        "raw_parity_constraints": rank_np(matrix) * nullspace_np(matrix).shape[0],
        "deduplicated_nonempty_parity_constraints": len(equations),
        "parity_equations_sha256": full_lc_parity_sha256(equations, n),
        "parity_hash_encoding": (
            "SHA-256 of ASCII full-lc-xor-v1 header followed by sorted masks, "
            "each fixed-width little-endian ceil(6*n/8) bytes"
        ),
        "projection_count_per_qubit": 6,
        "projection_matrices": [projection.tolist() for projection in LOCAL_CLIFFORD_PROJECTIONS],
        "symplectic_action_matrices": [action.tolist() for action in LOCAL_CLIFFORD_ACTIONS],
        "exact_equivalence": (
            "exists one-hot y such that S is invariant under the selected block-diagonal "
            "rank-one projection; this is iff an independent single-qubit Clifford maps S to CSS"
        ),
    }
