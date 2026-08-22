"""EXP-035: independent exact-distance certificate for PBB ``12_6_0193``.

The distance of a general stabilizer with check space S is

    min wt_s(v),  v in S^perp outside S.

For a complete quotient basis L_0,...,L_(2k-1), nonmembership is the exact
disjunction ``<v,L_j>_s = 1 for some j``.  Lattice translations preserve S,
symplectic weight, and pairing.  It is therefore enough to prove infeasibility
for orbit representatives whose translated quotient coordinates span all 2k
dual functionals: any nonzero logical pairs with at least one member of that
spanning union.  This experiment machine-checks a four-representative orbit
cover of rank 24 through independent NumPy and Python-integer GF(2) paths.

Each representative solve asks CP-SAT whether ``wt_s(v) <= 11``, ``v in
S^perp``, and ``<v,f>_s = 1``.  Records are atomic, seeded, capped, and
resumable.  Only four matching INFEASIBLE records whose transported dual rank
is 24, plus an independently checked weight-12 witness, can route canonical.
UNKNOWN/timeouts remain partial.  Catalogue distance fields are untrusted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from ortools.sat.python import cp_model

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.codes.bicycle import PBBSpec, build_pbb  # noqa: E402
from qec_research.gf2.linalg import (  # noqa: E402
    rank_bitset,
    rank_np,
    rref_bitset,
    rref_np,
    row_space_contains,
    rows_to_bitsets,
    solve_bitset,
)
from qec_research.symplectic.core import (  # noqa: E402
    StabilizerCode,
    lambda_swap,
    symplectic_product_matrix,
    symplectic_weight,
)

CODE_ID = "12_6_0193"
CATALOGUE_REL = "third_party/qcode-discovery/results/campaign7_publication_merged.jsonl"
CATALOGUE = ROOT / CATALOGUE_REL
CATALOGUE_LINE = 38
EXPECTED_TERMS = {
    "ell": 12,
    "m": 6,
    "A_terms": [[1, 2], [10, 3], [10, 4]],
    "B_terms": [[0, 0], [1, 5], [11, 4]],
    "C_terms": [[1, 3], [10, 3]],
    "D_terms": [[1, 5], [10, 5]],
}
EXPECTED_N = 144
EXPECTED_RANK = 132
EXPECTED_K = 12
LOWER_EXCLUSION_WEIGHT = 11
EXPECTED_DISTANCE = LOWER_EXCLUSION_WEIGHT + 1
DEFAULT_SEED = 2_026_081_934
DEFAULT_TIME_LIMIT_S = 7_200.0
DEFAULT_WORKERS = 2
DEFAULT_MEMORY_MB = 4_096
# Four cyclic submodules are required by this persisted cover.  Representative
# r uses the logical quotient vector whose basis-coefficient mask is entry r.
ORBIT_REPRESENTATIVE_MASKS = (1 << 2, 1 << 6, 1 << 16, 1 << 22)
STATE_DIR = ROOT / "results" / "partial_runs" / "pbb_12_6_0193_distance_sectors"
CANONICAL = ROOT / "results" / "certificates" / "pbb_12_6_0193_distance.json"
PARTIAL_SUMMARY = ROOT / "results" / "partial_runs" / "pbb_12_6_0193_distance.json"
QUARANTINE = ROOT / "results" / "quarantine" / "pbb_12_6_0193_distance.json"
COUNTEREXAMPLE = (
    ROOT / "results" / "certificates" / "pbb_12_6_0193_distance_counterexample.json"
)
REPORT = ROOT / "results" / "processed" / "exp035_pbb_exact_distance.json"
PARTIAL_REPORT = ROOT / "results" / "partial_runs" / "exp035_pbb_exact_distance.json"
# Independently found by EXP-029's pure-X sector optimizer; this experiment
# treats it solely as an input vector and rechecks every algebraic claim.
DEFAULT_WITNESS_X_SUPPORT = [9, 23, 31, 45, 59, 67, 73, 87, 101, 109, 123, 137]
_SMALL_WEIGHT_PROOF_CACHE: dict[str, dict[str, Any]] = {}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def matrix_sha256(matrix: np.ndarray) -> str:
    array = np.ascontiguousarray(np.asarray(matrix, dtype=np.uint8) & 1)
    return sha256_bytes(array.tobytes())


def canonical_json_sha256(value: Any) -> str:
    return sha256_bytes(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    )


def atomic_write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(value, indent=2, sort_keys=True) + "\n"
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def load_target_catalogue_row() -> tuple[dict[str, Any], str]:
    """Load exactly one target record and return it with the source line text."""
    matches: list[tuple[int, dict[str, Any], str]] = []
    with CATALOGUE.open(encoding="utf-8") as stream:
        for line_number, raw in enumerate(stream, start=1):
            row = json.loads(raw)
            if row.get("code_id") == CODE_ID:
                matches.append((line_number, row, raw.rstrip("\n")))
    if len(matches) != 1:
        raise RuntimeError(f"expected one {CODE_ID} row, found {len(matches)}")
    line_number, row, raw = matches[0]
    if line_number != CATALOGUE_LINE:
        raise RuntimeError(f"target moved from line {CATALOGUE_LINE} to {line_number}")
    actual = {name: row[name] for name in EXPECTED_TERMS}
    if actual != EXPECTED_TERMS:
        raise RuntimeError(f"target source terms changed: {actual!r}")
    return row, raw


def direct_monomial_matrix(ell: int, m: int, a: int, b: int) -> np.ndarray:
    """Rebuild one ring monomial without using the catalogue's code builder."""
    size = ell * m
    matrix = np.zeros((size, size), dtype=np.uint8)
    for i in range(ell):
        for j in range(m):
            source = i * m + j
            target = ((i + a) % ell) * m + (j + b) % m
            matrix[source, target] = 1
    return matrix


def direct_poly_matrix(ell: int, m: int, terms: Iterable[Iterable[int]]) -> np.ndarray:
    matrix = np.zeros((ell * m, ell * m), dtype=np.uint8)
    for a, b in terms:
        matrix ^= direct_monomial_matrix(ell, m, int(a), int(b))
    return matrix


def rebuild_target_independently(row: dict[str, Any]) -> StabilizerCode:
    """Build H directly from exact A/B/C/D terms in the stored convention."""
    ell, m = int(row["ell"]), int(row["m"])
    size = ell * m
    A = direct_poly_matrix(ell, m, row["A_terms"])
    B = direct_poly_matrix(ell, m, row["B_terms"])
    C = direct_poly_matrix(ell, m, row["C_terms"] or [])
    D = direct_poly_matrix(ell, m, row["D_terms"] or [])
    zero = np.zeros((size, size), dtype=np.uint8)
    # Stored convention: (x|z), H = (A B | C D ; 0 0 | B^T A^T).
    H = np.vstack(
        [np.hstack([A, B, C, D]), np.hstack([zero, zero, B.T, A.T])]
    ).astype(np.uint8)
    return StabilizerCode(H, name=CODE_ID)


def rebuild_metadata(row: dict[str, Any], raw_line: str) -> tuple[StabilizerCode, dict[str, Any]]:
    code = rebuild_target_independently(row)
    reference = build_pbb(
        PBBSpec(
            ell=row["ell"],
            m=row["m"],
            A=[tuple(term) for term in row["A_terms"]],
            B=[tuple(term) for term in row["B_terms"]],
            C=[tuple(term) for term in row["C_terms"]],
            D=[tuple(term) for term in row["D_terms"]],
            name=CODE_ID,
        )
    )
    validation = code.validate()
    metadata = {
        "convention": {
            "ring": "GF(2)[x,y]/(x^ell-1,y^m-1)",
            "flattening": "index(i,j)=i*m+j",
            "monomial_action": "P[index(i,j),index(i+a mod ell,j+b mod m)]=1",
            "symplectic_columns": "(x|z)",
            "check_matrix": "H=(A B | C D ; 0 0 | B^T A^T)",
        },
        "source": {
            "path": CATALOGUE_REL,
            "line": CATALOGUE_LINE,
            "code_id": CODE_ID,
            "raw_line_sha256": sha256_bytes((raw_line + "\n").encode("utf-8")),
            "catalogue_file_sha256": sha256_bytes(CATALOGUE.read_bytes()),
            "terms_sha256": canonical_json_sha256(EXPECTED_TERMS),
            "terms": EXPECTED_TERMS,
            "untrusted_catalogue_context": {
                "d": row.get("d"),
                "d_is_exact": row.get("d_is_exact"),
                "d_method": row.get("d_method"),
            },
        },
        "rebuild": {
            "n": code.n,
            "num_stored_checks": code.num_checks,
            "rank_numpy": validation["rank_numpy"],
            "rank_bitset": validation["rank_bitset"],
            "rank_paths_agree": validation["rank_agree"],
            "k": validation["k"],
            "commutes_numpy": validation["commutes_numpy"],
            "commutes_bitset": validation["commutes_bitset"],
            "non_css_as_stored": not validation["is_css"],
            "matrix_sha256": matrix_sha256(code.H),
            "independent_builder_matrix_matches_existing_builder": bool(
                np.array_equal(code.H, reference.H)
            ),
        },
    }
    expected = (
        code.n == EXPECTED_N
        and validation["rank_numpy"] == EXPECTED_RANK
        and validation["rank_bitset"] == EXPECTED_RANK
        and validation["k"] == EXPECTED_K
        and validation["commutes_numpy"]
        and validation["commutes_bitset"]
        and not validation["is_css"]
        and metadata["rebuild"]["independent_builder_matrix_matches_existing_builder"]
    )
    if not expected:
        raise RuntimeError(f"target reconstruction failed: {metadata['rebuild']!r}")
    return code, metadata


def logical_basis_metadata(code: StabilizerCode) -> tuple[np.ndarray, dict[str, Any]]:
    logicals = code.logical_basis()
    expected_rows = 2 * code.k
    if logicals.shape != (expected_rows, 2 * code.n):
        raise RuntimeError(f"bad logical basis shape {logicals.shape}")
    centralizer = not symplectic_product_matrix(logicals, code.H).any()
    combined_rank_np = rank_np(np.vstack([code.H, logicals]))
    combined_rank_bitset = rank_bitset(
        rows_to_bitsets(np.vstack([code.H, logicals])), 2 * code.n
    )
    expected_rank = code.rank() + expected_rows
    meta = {
        "basis_rows": expected_rows,
        "basis_sha256": matrix_sha256(logicals),
        "all_basis_rows_in_centralizer": bool(centralizer),
        "rank_H_plus_basis_numpy": combined_rank_np,
        "rank_H_plus_basis_bitset": combined_rank_bitset,
        "expected_rank_H_plus_basis": expected_rank,
        "basis_is_complete_quotient_basis": bool(
            centralizer
            and combined_rank_np == expected_rank
            and combined_rank_bitset == expected_rank
        ),
    }
    if not meta["basis_is_complete_quotient_basis"]:
        raise RuntimeError(f"incomplete logical quotient basis: {meta!r}")
    return logicals, meta


def coefficient_support(mask: int, width: int) -> list[int]:
    return [index for index in range(width) if (mask >> index) & 1]


def representative_logical(logicals: np.ndarray, representative: int) -> np.ndarray:
    if not 0 <= representative < len(ORBIT_REPRESENTATIVE_MASKS):
        raise ValueError(
            f"representative must be in [0,{len(ORBIT_REPRESENTATIVE_MASKS) - 1}]"
        )
    mask = ORBIT_REPRESENTATIVE_MASKS[representative]
    rows = coefficient_support(mask, logicals.shape[0])
    return np.bitwise_xor.reduce(logicals[rows], axis=0)


def translation_column_permutation(
    ell: int, m: int, delta_i: int, delta_j: int
) -> np.ndarray:
    """Map each old (x|z) column to its translated new column.

    Qubits are ``block * (ell*m) + i*m+j`` for blocks 0 and 1.  The same
    permutation is applied independently to x and z columns.
    """
    dim = ell * m
    base = [
        ((i + delta_i) % ell) * m + (j + delta_j) % m
        for i in range(ell)
        for j in range(m)
    ]
    qubits = np.array(base + [index + dim for index in base], dtype=np.int64)
    return np.concatenate([qubits, qubits + 2 * dim])


def translate_rows(rows: np.ndarray, permutation: np.ndarray) -> np.ndarray:
    source = np.asarray(rows, dtype=np.uint8) & 1
    one_dimensional = source.ndim == 1
    source_2d = source[None, :] if one_dimensional else source
    translated = np.zeros_like(source_2d)
    translated[:, permutation] = source_2d
    return translated[0] if one_dimensional else translated


def _reduce_with_rref_bitset(rows: list[int], pivots: list[int], target: int) -> int:
    reduced = target
    for row, pivot in zip(rows, pivots, strict=True):
        if (reduced >> pivot) & 1:
            reduced ^= row
    return reduced


def _gf2_inverse(matrix: np.ndarray) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=np.uint8) & 1
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("GF(2) inverse requires a square matrix")
    size = matrix.shape[0]
    reduced, pivots = rref_np(
        np.hstack([matrix, np.eye(size, dtype=np.uint8)])
    )
    if reduced.shape[0] != size or pivots != list(range(size)):
        raise ValueError("matrix is singular over GF(2)")
    if not np.array_equal(reduced[:, :size], np.eye(size, dtype=np.uint8)):
        raise RuntimeError("GF(2) inversion did not reduce to identity")
    return reduced[:, size:]


def _symplectic_pair_bits(left: int, right: int, n: int) -> int:
    mask = (1 << n) - 1
    left_x, left_z = left & mask, left >> n
    right_x, right_z = right & mask, right >> n
    return ((left_x & right_z).bit_count() + (left_z & right_x).bit_count()) & 1


def small_weight_centralizer_exclusion(code: StabilizerCode) -> dict[str, Any]:
    """Exhaust all Pauli errors of weight at most five by syndrome MITM.

    There are 3n single-qubit Pauli atoms.  Weight 2--5 zero syndromes are
    respectively atom-pair, atom-pair/atom, pair/pair, and triple/pair XOR
    collisions, with all participating qubits required to be distinct.
    """
    cache_key = matrix_sha256(code.H)
    cached = _SMALL_WEIGHT_PROOF_CACHE.get(cache_key)
    if cached is not None:
        return cached
    n = code.n
    checks = rows_to_bitsets(code.H)
    atom_rows: list[tuple[int, int, int]] = []
    for qubit in range(n):
        for pauli in range(3):  # 0=X, 1=Z, 2=Y
            operator = 0
            if pauli != 1:
                operator |= 1 << qubit
            if pauli != 0:
                operator |= 1 << (n + qubit)
            syndrome = sum(
                _symplectic_pair_bits(operator, check, n) << check_index
                for check_index, check in enumerate(checks)
            )
            atom_rows.append((syndrome, qubit, pauli))

    atoms_by_syndrome: dict[int, list[tuple[int, int]]] = {}
    for syndrome, qubit, pauli in atom_rows:
        atoms_by_syndrome.setdefault(syndrome, []).append((qubit, pauli))
    zero_weight_one = [
        [qubit, pauli]
        for syndrome, qubit, pauli in atom_rows
        if syndrome == 0
    ]

    pairs_by_syndrome: dict[int, list[tuple[int, int, int, int]]] = {}
    pair_count = 0
    weight_four_witness: list[int] | None = None
    for first, (left, left_qubit, left_pauli) in enumerate(atom_rows):
        for right, right_qubit, right_pauli in atom_rows[first + 1 :]:
            if left_qubit == right_qubit:
                continue
            pair_count += 1
            syndrome = left ^ right
            prior = pairs_by_syndrome.setdefault(syndrome, [])
            if weight_four_witness is None:
                for (
                    prior_left_qubit,
                    prior_left_pauli,
                    prior_right_qubit,
                    prior_right_pauli,
                ) in prior:
                    if len(
                        {
                            left_qubit,
                            right_qubit,
                            prior_left_qubit,
                            prior_right_qubit,
                        }
                    ) == 4:
                        weight_four_witness = [
                            left_qubit,
                            left_pauli,
                            right_qubit,
                            right_pauli,
                            prior_left_qubit,
                            prior_left_pauli,
                            prior_right_qubit,
                            prior_right_pauli,
                        ]
                        break
            prior.append((left_qubit, left_pauli, right_qubit, right_pauli))

    weight_two_witness = (
        list(pairs_by_syndrome[0][0]) if 0 in pairs_by_syndrome else None
    )
    weight_three_witness: list[int] | None = None
    for syndrome, pairs in pairs_by_syndrome.items():
        for atom_qubit, atom_pauli in atoms_by_syndrome.get(syndrome, []):
            for left_qubit, left_pauli, right_qubit, right_pauli in pairs:
                if atom_qubit not in (left_qubit, right_qubit):
                    weight_three_witness = [
                        left_qubit,
                        left_pauli,
                        right_qubit,
                        right_pauli,
                        atom_qubit,
                        atom_pauli,
                    ]
                    break
            if weight_three_witness is not None:
                break
        if weight_three_witness is not None:
            break

    weight_five_witness: list[int] | None = None
    triple_count = 0
    for first, (left, left_qubit, left_pauli) in enumerate(atom_rows):
        for second in range(first + 1, len(atom_rows)):
            middle, middle_qubit, middle_pauli = atom_rows[second]
            if left_qubit == middle_qubit:
                continue
            for right, right_qubit, right_pauli in atom_rows[second + 1 :]:
                if right_qubit in (left_qubit, middle_qubit):
                    continue
                triple_count += 1
                syndrome = left ^ middle ^ right
                for (
                    pair_left_qubit,
                    pair_left_pauli,
                    pair_right_qubit,
                    pair_right_pauli,
                ) in pairs_by_syndrome.get(syndrome, []):
                    if len(
                        {
                            left_qubit,
                            middle_qubit,
                            right_qubit,
                            pair_left_qubit,
                            pair_right_qubit,
                        }
                    ) == 5:
                        weight_five_witness = [
                            left_qubit,
                            left_pauli,
                            middle_qubit,
                            middle_pauli,
                            right_qubit,
                            right_pauli,
                            pair_left_qubit,
                            pair_left_pauli,
                            pair_right_qubit,
                            pair_right_pauli,
                        ]
                        break
                if weight_five_witness is not None:
                    break
            if weight_five_witness is not None:
                break
        if weight_five_witness is not None:
            break

    witnesses = {
        "weight_1": zero_weight_one[0] if zero_weight_one else None,
        "weight_2": weight_two_witness,
        "weight_3": weight_three_witness,
        "weight_4": weight_four_witness,
        "weight_5": weight_five_witness,
    }
    complete = all(witness is None for witness in witnesses.values())
    atom_bytes = b"".join(
        syndrome.to_bytes((code.num_checks + 7) // 8, "little")
        for syndrome, _, _ in atom_rows
    )
    pair_bytes = b"".join(
        syndrome.to_bytes((code.num_checks + 7) // 8, "little")
        for syndrome, pairs in sorted(pairs_by_syndrome.items())
        for _ in pairs
    )
    proof = {
        "method": (
            "exact Python-integer GF(2) syndrome meet-in-the-middle over all "
            "3n single-qubit X/Z/Y atoms; distinct-qubit filters enforce "
            "symplectic weight"
        ),
        "checked_through_weight": 5,
        "pauli_type_encoding": {"0": "X", "1": "Z", "2": "Y"},
        "atom_count": len(atom_rows),
        "distinct_qubit_pair_count": pair_count,
        "distinct_qubit_triple_count": triple_count,
        "atom_syndromes_sha256": sha256_bytes(atom_bytes),
        "pair_syndromes_multiset_sha256": sha256_bytes(pair_bytes),
        "zero_syndrome_witnesses": witnesses,
        "no_centralizer_of_weight_at_most_5": complete,
        "certified_distance_lower_bound": 6 if complete else 1,
    }
    if code.n == EXPECTED_N and not complete:
        raise RuntimeError(f"unexpected small centralizer: {witnesses!r}")
    _SMALL_WEIGHT_PROOF_CACHE[cache_key] = proof
    return proof


_WEIGHT_SIX_CACHE: dict[str, dict[str, Any]] = {}


def weight_six_zero_syndrome_classification(code: StabilizerCode) -> dict[str, Any]:
    """Enumerate every weight-6 zero-syndrome Pauli and classify it against S.

    Completeness: a weight-6 centralizer element supports six distinct qubits
    carrying one Pauli atom each, and every 3+3 split of that support is a
    pair of qubit-disjoint distinct-qubit atom triples with equal syndromes.
    Scanning all distinct-qubit triples and matching equal-syndrome,
    qubit-disjoint pairs therefore finds every such element (ten splits per
    element).  Combined with the exact weight-<=5 exclusion, "every weight-6
    zero-syndrome vector lies in the stabilizer row space" is exactly "no
    logical of symplectic weight <= 6", i.e. d >= 7; any nonmember is a
    verified weight-6 logical, i.e. d = 6 exactly.
    """
    cache_key = matrix_sha256(code.H)
    cached = _WEIGHT_SIX_CACHE.get(cache_key)
    if cached is not None:
        return cached
    n = code.n
    checks = rows_to_bitsets(code.H)
    atoms: list[tuple[int, int, int]] = []
    for qubit in range(n):
        for pauli in range(3):  # 0=X, 1=Z, 2=Y
            operator = 0
            if pauli != 1:
                operator |= 1 << qubit
            if pauli != 0:
                operator |= 1 << (n + qubit)
            syndrome = sum(
                _symplectic_pair_bits(operator, check, n) << check_index
                for check_index, check in enumerate(checks)
            )
            atoms.append((syndrome, 1 << qubit, operator))

    atom_count = len(atoms)
    triple_count = 0
    disjoint_equal_syndrome_pairs = 0
    triples_by_syndrome: dict[int, list[tuple[int, int]]] = {}
    distinct: set[int] = set()
    for first in range(atom_count):
        syndrome_first, qubit_first, operator_first = atoms[first]
        for second in range(first + 1, atom_count):
            syndrome_second, qubit_second, operator_second = atoms[second]
            if qubit_second == qubit_first:
                continue
            syndrome_pair = syndrome_first ^ syndrome_second
            qubit_pair = qubit_first | qubit_second
            operator_pair = operator_first | operator_second
            for third in range(second + 1, atom_count):
                syndrome_third, qubit_third, operator_third = atoms[third]
                if qubit_third & qubit_pair:
                    continue
                triple_count += 1
                syndrome = syndrome_pair ^ syndrome_third
                qubits = qubit_pair | qubit_third
                operator = operator_pair | operator_third
                bucket = triples_by_syndrome.get(syndrome)
                if bucket is None:
                    triples_by_syndrome[syndrome] = [(qubits, operator)]
                    continue
                for prior_qubits, prior_operator in bucket:
                    if not (prior_qubits & qubits):
                        disjoint_equal_syndrome_pairs += 1
                        distinct.add(prior_operator | operator)
                bucket.append((qubits, operator))

    qubit_mask = (1 << n) - 1
    H_matrix = np.asarray(code.H, dtype=np.uint8) & 1
    H_rank = rank_np(H_matrix)
    members = 0
    nonmembers: list[dict[str, Any]] = []
    width = 2 * n
    for operator in sorted(distinct):
        x_bits = operator & qubit_mask
        z_bits = operator >> n
        support = x_bits | z_bits
        if support.bit_count() != 6:
            raise AssertionError("collision produced a non-weight-6 vector")
        for check in checks:
            if _symplectic_pair_bits(operator, check, n):
                raise AssertionError("collision produced a nonzero syndrome")
        in_row_space_bitset = row_space_contains(checks, width, operator)
        vector = np.array(
            [(operator >> index) & 1 for index in range(width)], dtype=np.uint8
        )
        in_row_space_numpy = rank_np(np.vstack([H_matrix, vector])) == H_rank
        if in_row_space_bitset != in_row_space_numpy:
            raise AssertionError("bitset and numpy membership disagree")
        if in_row_space_bitset:
            members += 1
        else:
            nonmembers.append(witness_with_crosschecks(code, vector))

    pure_z_weight6_rows = {
        check
        for check in checks
        if (check & qubit_mask) == 0 and (check >> n).bit_count() == 6
    }
    vector_bytes = b"".join(
        operator.to_bytes((width + 7) // 8, "little")
        for operator in sorted(distinct)
    )
    proof = {
        "method": (
            "exact Python-integer GF(2) triple-triple syndrome meet-in-the-"
            "middle over all 3n single-qubit atoms; qubit-disjoint filters "
            "enforce symplectic weight exactly 6; every distinct collision "
            "vector re-verified to have zero syndrome and classified against "
            "the stabilizer row space by independent bitset and NumPy paths"
        ),
        "atom_count": atom_count,
        "distinct_qubit_triple_count": triple_count,
        "qubit_disjoint_equal_syndrome_pairs": disjoint_equal_syndrome_pairs,
        "distinct_weight_six_zero_syndrome_vectors": len(distinct),
        "vectors_sha256": sha256_bytes(vector_bytes),
        "stabilizer_members": members,
        "nonmember_logicals": nonmembers,
        "no_weight6_logical": not nonmembers,
        "stored_pure_z_weight6_rows": len(pure_z_weight6_rows),
        "matches_stored_pure_z_weight6_rows": distinct == pure_z_weight6_rows,
        "complete": True,
        "certified_lower_bound_with_weight5_exclusion": (
            7 if not nonmembers else 6
        ),
    }
    # A target nonmember is a decisive d=6 counterexample; it is returned in
    # full (never raised away) so the router can persist it atomically.
    _WEIGHT_SIX_CACHE[cache_key] = proof
    return proof


def translation_orbit_proof(
    code: StabilizerCode, logicals: np.ndarray
) -> dict[str, Any]:
    """Certify the four representative translation orbits span the dual."""
    ell, m = EXPECTED_TERMS["ell"], EXPECTED_TERMS["m"]
    width = logicals.shape[0]
    if width != 2 * code.k:
        raise ValueError("logical quotient basis has the wrong dimension")

    check_bits = rows_to_bitsets(code.H)
    check_rref, check_pivots, _ = rref_bitset(check_bits, 2 * code.n)
    logical_bits = rows_to_bitsets(logicals)
    pairing = symplectic_product_matrix(logicals, logicals)
    pairing_inverse = _gf2_inverse(pairing)
    pairing_rows = rows_to_bitsets(pairing)
    all_coordinate_rows: list[list[int]] = []
    per_representative: list[dict[str, Any]] = []
    shifts: list[dict[str, Any]] = []

    translated_cache: list[tuple[int, int, np.ndarray, np.ndarray]] = []
    for delta_i in range(ell):
        for delta_j in range(m):
            permutation = translation_column_permutation(
                ell, m, delta_i, delta_j
            )
            translated_checks = translate_rows(code.H, permutation)
            translated_check_bits = rows_to_bitsets(translated_checks)
            rowspace_ok = all(
                _reduce_with_rref_bitset(check_rref, check_pivots, row) == 0
                for row in translated_check_bits
            )
            translated_logicals = translate_rows(logicals, permutation)
            pairing_ok = np.array_equal(
                symplectic_product_matrix(translated_logicals, translated_logicals),
                pairing,
            )
            shifts.append(
                {
                    "shift": [delta_i, delta_j],
                    "column_permutation_sha256": canonical_json_sha256(
                        permutation.astype(int).tolist()
                    ),
                    "permutation_is_bijection": bool(
                        np.array_equal(np.sort(permutation), np.arange(2 * code.n))
                    ),
                    "maps_stabilizer_rowspace_to_itself_bitset": rowspace_ok,
                    "logical_pairing_matrix_preserved_numpy": bool(pairing_ok),
                }
            )
            translated_cache.append(
                (delta_i, delta_j, permutation, translated_logicals)
            )

    for representative, coefficient_mask in enumerate(
        ORBIT_REPRESENTATIVE_MASKS
    ):
        logical = representative_logical(logicals, representative)
        logical_bit = rows_to_bitsets(logical[None, :])[0]
        orbit_coordinates: list[int] = []
        translations: list[dict[str, Any]] = []
        for delta_i, delta_j, permutation, translated_logicals in translated_cache:
            translated = translate_rows(logical, permutation)
            translated_bit = rows_to_bitsets(translated[None, :])[0]
            pairing_numpy = symplectic_product_matrix(
                translated[None, :], logicals
            )[0]
            coordinate_numpy = (
                pairing_numpy.astype(np.int64)
                @ pairing_inverse.astype(np.int64)
                % 2
            ).astype(np.uint8)
            pairing_mask = sum(
                _symplectic_pair_bits(translated_bit, basis_bit, code.n) << index
                for index, basis_bit in enumerate(logical_bits)
            )
            coordinate_bitset = solve_bitset(pairing_rows, width, pairing_mask)
            if coordinate_bitset is None:
                raise RuntimeError("nondegenerate quotient pairing solve failed")
            coordinate_numpy_mask = sum(
                int(bit) << index for index, bit in enumerate(coordinate_numpy)
            )
            reconstructed = 0
            for index in coefficient_support(coordinate_bitset, width):
                reconstructed ^= logical_bits[index]
            residual_in_stabilizer = (
                _reduce_with_rref_bitset(
                    check_rref,
                    check_pivots,
                    translated_bit ^ reconstructed,
                )
                == 0
            )
            centralizer_bitset = all(
                _symplectic_pair_bits(translated_bit, check, code.n) == 0
                for check in check_bits
            )
            pairing_transport = np.array_equal(
                symplectic_product_matrix(
                    translated[None, :], translated_logicals
                ),
                symplectic_product_matrix(logical[None, :], logicals),
            )
            coordinate_ok = coordinate_numpy_mask == coordinate_bitset
            if not (
                coordinate_ok
                and residual_in_stabilizer
                and centralizer_bitset
                and pairing_transport
            ):
                raise RuntimeError(
                    f"translation proof failed at representative {representative}, "
                    f"shift {(delta_i, delta_j)}"
                )
            orbit_coordinates.append(coordinate_bitset)
            translations.append(
                {
                    "shift": [delta_i, delta_j],
                    "quotient_coefficient_mask_hex": f"0x{coordinate_bitset:06x}",
                    "quotient_coefficient_support": coefficient_support(
                        coordinate_bitset, width
                    ),
                    "numpy_and_bitset_coordinates_agree": coordinate_ok,
                    "residual_is_in_stabilizer_rowspace_bitset": (
                        residual_in_stabilizer
                    ),
                    "translated_logical_is_centralizer_bitset": centralizer_bitset,
                    "symplectic_pairing_transport_identity_numpy": bool(
                        pairing_transport
                    ),
                    "symplectic_weight_preserved": (
                        symplectic_weight(translated) == symplectic_weight(logical)
                    ),
                }
            )
        coordinate_rows = [
            [(coordinate >> index) & 1 for index in range(width)]
            for coordinate in orbit_coordinates
        ]
        orbit_rank_numpy = rank_np(np.asarray(coordinate_rows, dtype=np.uint8))
        orbit_rank_bitset = rank_bitset(orbit_coordinates, width)
        all_coordinate_rows.append(orbit_coordinates)
        per_representative.append(
            {
                "representative": representative,
                "original_basis_coefficient_mask_hex": (
                    f"0x{coefficient_mask:06x}"
                ),
                "original_basis_coefficients": [
                    (coefficient_mask >> index) & 1 for index in range(width)
                ],
                "original_basis_coefficient_support": coefficient_support(
                    coefficient_mask, width
                ),
                "logical_vector_sha256": matrix_sha256(logical[None, :]),
                "detector_vector_sha256": matrix_sha256(
                    lambda_swap(logical[None, :])
                ),
                "orbit_size_with_multiplicity": ell * m,
                "orbit_unique_coordinate_count": len(set(orbit_coordinates)),
                "orbit_quotient_rank_numpy": orbit_rank_numpy,
                "orbit_quotient_rank_bitset": orbit_rank_bitset,
                "translations": translations,
            }
        )

    flat_coordinates = [
        coordinate for orbit in all_coordinate_rows for coordinate in orbit
    ]
    flat_numpy = np.asarray(
        [
            [(coordinate >> index) & 1 for index in range(width)]
            for coordinate in flat_coordinates
        ],
        dtype=np.uint8,
    )
    union_rank_numpy = rank_np(flat_numpy)
    union_rank_bitset = rank_bitset(flat_coordinates, width)
    proof = {
        "schema": "pbb-translation-orbit-cover-v1",
        "translation_group": "Z_12 x Z_6",
        "translation_count": ell * m,
        "permutation_convention": (
            "old qubit (block,i,j) maps to "
            "(block,i+delta_i mod 12,j+delta_j mod 6); same on x and z"
        ),
        "stabilizer_translation_checks": shifts,
        "every_translation_is_bijection": all(
            item["permutation_is_bijection"] for item in shifts
        ),
        "every_translation_preserves_stabilizer_rowspace_bitset": all(
            item["maps_stabilizer_rowspace_to_itself_bitset"] for item in shifts
        ),
        "every_translation_preserves_logical_pairing_numpy": all(
            item["logical_pairing_matrix_preserved_numpy"] for item in shifts
        ),
        "logical_pairing_matrix_sha256": matrix_sha256(pairing),
        "logical_pairing_rank_numpy": rank_np(pairing),
        "logical_pairing_rank_bitset": rank_bitset(pairing_rows, width),
        "representatives": per_representative,
        "selected_orbit_union_rank_numpy": union_rank_numpy,
        "selected_orbit_union_rank_bitset": union_rank_bitset,
        "full_dual_rank": width,
        "selected_orbits_span_full_dual": bool(
            union_rank_numpy == width and union_rank_bitset == width
        ),
        "coverage_argument": (
            "If a nonzero quotient class v existed below the cap, "
            "nondegeneracy gives a covered functional f with <v,f>_s=1. "
            "Translation transports representative infeasibility to every f "
            "in its orbit while preserving centralizer membership and weight."
        ),
    }
    if not (
        proof["every_translation_is_bijection"]
        and proof["every_translation_preserves_stabilizer_rowspace_bitset"]
        and proof["every_translation_preserves_logical_pairing_numpy"]
        and proof["selected_orbits_span_full_dual"]
    ):
        raise RuntimeError("translation-orbit coverage proof failed")
    return proof


def detector_for_sector(logicals: np.ndarray, sector: int) -> np.ndarray:
    logical = representative_logical(logicals, sector)
    return lambda_swap(logical[None, :])[0]


def sector_seed(base_seed: int, sector: int) -> int:
    digest = hashlib.sha256(
        f"{CODE_ID}:orbit-representative:{base_seed}:{sector}".encode()
    ).digest()
    return int.from_bytes(digest[:4], "little") & 0x7FFF_FFFF


def sector_identity(
    code: StabilizerCode,
    logicals: np.ndarray,
    sector: int,
    base_seed: int,
    weight_cap: int,
) -> dict[str, Any]:
    detector = detector_for_sector(logicals, sector)
    coefficient_mask = ORBIT_REPRESENTATIVE_MASKS[sector]
    return {
        "schema": "pbb-exact-distance-orbit-sector-v2",
        "code_id": CODE_ID,
        "sector": sector,
        "orbit_representative": sector,
        "logical_basis_coefficient_mask_hex": f"0x{coefficient_mask:06x}",
        "logical_basis_coefficient_support": coefficient_support(
            coefficient_mask, logicals.shape[0]
        ),
        "seed": sector_seed(base_seed, sector),
        "weight_cap": weight_cap,
        "matrix_sha256": matrix_sha256(code.H),
        "logical_basis_sha256": matrix_sha256(logicals),
        "detector_sha256": matrix_sha256(detector[None, :]),
    }


def add_even_parity(model: cp_model.CpModel, variables: list[Any], row: np.ndarray) -> None:
    selected = [variables[j] for j in np.flatnonzero(row)]
    if selected:
        # AddBoolXOr is odd parity; appending constant true enforces even parity.
        model.add_bool_xor(selected + [model.new_constant(1)])


def solve_lower_sector(
    code: StabilizerCode,
    logicals: np.ndarray,
    sector: int,
    *,
    base_seed: int,
    weight_cap: int,
    time_limit_s: float,
    workers: int,
    memory_mb: int,
    log_path: Path,
) -> dict[str, Any]:
    """Decide one bounded sector; INFEASIBLE excludes all errors in it."""
    identity = sector_identity(code, logicals, sector, base_seed, weight_cap)
    n = code.n
    model = cp_model.CpModel()
    variables = [model.new_bool_var(f"v{j}") for j in range(2 * n)]
    weights = [model.new_bool_var(f"w{j}") for j in range(n)]
    for row in lambda_swap(code.H):
        add_even_parity(model, variables, row)
    detector = detector_for_sector(logicals, sector)
    model.add_bool_xor([variables[j] for j in np.flatnonzero(detector)])
    for qubit in range(n):
        model.add_max_equality(weights[qubit], [variables[qubit], variables[n + qubit]])
    model.add(sum(weights) <= weight_cap)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit_s
    solver.parameters.num_workers = workers
    solver.parameters.random_seed = identity["seed"]
    solver.parameters.max_memory_in_mb = memory_mb
    solver.parameters.log_search_progress = True
    solver.parameters.log_to_stdout = False
    solver.parameters.log_to_response = True
    started_utc = datetime.now(timezone.utc).isoformat()
    started = time.perf_counter()
    status_code = solver.solve(model)
    wall = time.perf_counter() - started
    status = solver.status_name(status_code)
    solution: dict[str, Any] | None = None
    if status in {"OPTIMAL", "FEASIBLE"}:
        vector = [int(solver.value(variable)) for variable in variables]
        solution = {
            "vector": vector,
            "x": vector[:n],
            "z": vector[n:],
            "support": [j for j in range(n) if vector[j] or vector[n + j]],
            "symplectic_weight": int(sum(solver.value(weight) for weight in weights)),
        }
    trace = solver.response_stats()
    response = solver.response_proto
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(trace + "\n", encoding="utf-8")
    return {
        **identity,
        "formulation": {
            "variables": "v=(x|z) in GF(2)^(2n)",
            "centralizer": "Lambda(H) v^T = 0",
            "nontriviality_sector": (
                "v . Lambda(sum_j coefficient_j L_j) = 1"
            ),
            "weight": "sum_q OR(x_q,z_q) <= weight_cap",
            "coverage": (
                "translation orbit of this representative; selected orbit "
                "union must have independently checked quotient rank 2k"
            ),
        },
        "solver": {
            "name": "OR-Tools CP-SAT",
            "version": _ortools_version(),
            "status": status,
            "status_code": int(status_code),
            "proof_complete_for_sector": status == "INFEASIBLE",
            "time_limit_s": time_limit_s,
            "workers": workers,
            "random_seed": identity["seed"],
            "max_memory_mb": memory_mb,
            "wall_time_s": wall,
            "started_utc": started_utc,
            "num_conflicts": int(response.num_conflicts),
            "num_branches": int(response.num_branches),
            "deterministic_time": float(response.deterministic_time),
            "trace_path": str(log_path.relative_to(ROOT)),
            "trace_sha256": sha256_bytes((trace + "\n").encode("utf-8")),
        },
        "solution": solution,
    }


def _ortools_version() -> str:
    import ortools

    return ortools.__version__


def pysat_sector_path(sector: int) -> Path:
    return STATE_DIR / f"sector_{sector:02d}_pysat.json"


def pysat_sector_log_path(sector: int) -> Path:
    return STATE_DIR / "traces" / f"sector_{sector:02d}_pysat.txt"


def build_sector_cnf(
    code: StabilizerCode, detector: np.ndarray, weight_cap: int
):
    """Build the exact CNF for one bounded sector.

    Variables 1..2n are the (x|z) bits; auxiliary variables encode Tseitin
    XOR chains (parity constraints) and the sequential-counter cardinality
    bound.  The encoding is equisatisfiable with the CP-SAT model: v in the
    centralizer, <v, detector>_s = 1, symplectic weight <= weight_cap.
    """

    from pysat.card import CardEnc, EncType
    from pysat.formula import CNF, IDPool

    n = code.n
    cnf = CNF()
    vpool = IDPool(start_from=3 * n + 1)

    def xor_constraint(literals: list[int], rhs: int) -> None:
        accumulator: int | None = None
        for literal in literals:
            if accumulator is None:
                accumulator = literal
                continue
            auxiliary = vpool.id()
            cnf.extend(
                [
                    [-accumulator, -literal, -auxiliary],
                    [accumulator, literal, -auxiliary],
                    [accumulator, -literal, auxiliary],
                    [-accumulator, literal, auxiliary],
                ]
            )
            accumulator = auxiliary
        if accumulator is None:
            if rhs:
                cnf.append([])
            return
        cnf.append([accumulator] if rhs else [-accumulator])

    for row in lambda_swap(code.H):
        literals = [int(j) + 1 for j in np.flatnonzero(row)]
        xor_constraint(literals, 0)
    xor_constraint([int(j) + 1 for j in np.flatnonzero(detector)], 1)
    weight_literals = list(range(2 * n + 1, 3 * n + 1))
    for qubit in range(n):
        x_literal = qubit + 1
        z_literal = n + qubit + 1
        w_literal = 2 * n + qubit + 1
        cnf.extend(
            [
                [-x_literal, w_literal],
                [-z_literal, w_literal],
                [x_literal, z_literal, -w_literal],
            ]
        )
    cnf.extend(
        CardEnc.atmost(
            lits=weight_literals,
            bound=weight_cap,
            vpool=vpool,
            encoding=EncType.seqcounter,
        ).clauses
    )
    return cnf


def solve_lower_sector_pysat(
    code: StabilizerCode,
    logicals: np.ndarray,
    sector: int,
    *,
    base_seed: int,
    weight_cap: int,
    conflict_budget: int,
    log_path: Path,
    solver_name: str = "cadical195",
) -> dict[str, Any]:
    """Decide one bounded sector with an independent CDCL SAT engine.

    Frozen amendment (2026-08-15, before any completed sector record
    existed): a representative counts as covered iff EITHER the CP-SAT
    backend returns a proof-complete INFEASIBLE OR this backend returns a
    complete UNSAT on the equisatisfiable CNF.  A SAT answer is accepted
    only after the decoded vector passes the same independent witness
    crosschecks as every other claim.
    """

    import pysat
    from pysat.solvers import Solver

    identity = sector_identity(code, logicals, sector, base_seed, weight_cap)
    n = code.n
    detector = detector_for_sector(logicals, sector)
    cnf = build_sector_cnf(code, detector, weight_cap)
    started_utc = datetime.now(timezone.utc).isoformat()
    started = time.perf_counter()
    with Solver(name=solver_name, bootstrap_with=cnf.clauses) as engine:
        if conflict_budget > 0:
            engine.conf_budget(conflict_budget)
            answer = engine.solve_limited(expect_interrupt=False)
        else:
            answer = engine.solve()
        stats = dict(engine.accum_stats())
        model = engine.get_model() if answer else None
    wall = time.perf_counter() - started
    if answer is True:
        status = "SAT"
    elif answer is False:
        status = "UNSAT"
    else:
        status = "UNDECIDED_BUDGET"
    solution: dict[str, Any] | None = None
    if model is not None:
        assignment = {abs(literal): literal > 0 for literal in model}
        vector = [1 if assignment.get(index + 1, False) else 0 for index in range(2 * n)]
        solution = {
            "vector": vector,
            "x": vector[:n],
            "z": vector[n:],
            "support": [j for j in range(n) if vector[j] or vector[n + j]],
            "symplectic_weight": sum(
                1 for j in range(n) if vector[j] or vector[n + j]
            ),
        }
    trace = json.dumps(
        {"answer": status, "stats": stats, "wall_s": wall}, sort_keys=True
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(trace + "\n", encoding="utf-8")
    return {
        **identity,
        "formulation": {
            "variables": "v=(x|z) in GF(2)^(2n) as CNF literals 1..2n",
            "centralizer": "Tseitin XOR chains for Lambda(H) v^T = 0",
            "nontriviality_sector": "Tseitin XOR chain = 1 for the detector",
            "weight": "w_q <-> OR(x_q,z_q); seqcounter atmost weight_cap",
            "coverage": (
                "translation orbit of this representative; selected orbit "
                "union must have independently checked quotient rank 2k"
            ),
        },
        "solver": {
            "name": f"PySAT {solver_name}",
            "version": pysat.__version__,
            "backend": "pysat",
            "status": status,
            "proof_complete_for_sector": status == "UNSAT",
            "conflict_budget": conflict_budget,
            "clauses": len(cnf.clauses),
            "cnf_variables": cnf.nv,
            "wall_time_s": wall,
            "started_utc": started_utc,
            "accum_stats": stats,
            "trace_path": str(log_path.relative_to(ROOT)),
            "trace_sha256": sha256_bytes((trace + "\n").encode("utf-8")),
        },
        "solution": solution,
    }


def verify_witness_bitset(code: StabilizerCode, vector: np.ndarray) -> dict[str, Any]:
    """Independent pure-Python integer checks of centralizer and row-space claims."""
    vector = np.asarray(vector, dtype=np.uint8) & 1
    if vector.shape != (2 * code.n,):
        raise ValueError(f"bad witness shape {vector.shape}")
    n = code.n
    x_bits = sum(int(vector[j]) << j for j in range(n))
    z_bits = sum(int(vector[n + j]) << j for j in range(n))
    witness_bits = x_bits | (z_bits << n)
    stabilizers = rows_to_bitsets(code.H)
    centralizer = True
    for stabilizer in stabilizers:
        sx = stabilizer & ((1 << n) - 1)
        sz = stabilizer >> n
        if ((x_bits & sz).bit_count() + (z_bits & sx).bit_count()) & 1:
            centralizer = False
            break
    in_stabilizer = row_space_contains(stabilizers, 2 * n, witness_bits)
    support = [j for j in range(n) if (x_bits >> j) & 1 or (z_bits >> j) & 1]
    return {
        "x": vector[:n].astype(int).tolist(),
        "z": vector[n:].astype(int).tolist(),
        "support": support,
        "symplectic_weight": len(support),
        "commutes_with_all_stabilizers_bitset": centralizer,
        "in_stabilizer_row_space_bitset": in_stabilizer,
        "nontrivial_logical_bitset": centralizer and not in_stabilizer,
        "vector_sha256": matrix_sha256(vector[None, :]),
    }


def default_witness(code: StabilizerCode) -> np.ndarray:
    vector = np.zeros(2 * code.n, dtype=np.uint8)
    vector[DEFAULT_WITNESS_X_SUPPORT] = 1
    return vector


def witness_with_crosschecks(code: StabilizerCode, vector: np.ndarray) -> dict[str, Any]:
    bitset = verify_witness_bitset(code, vector)
    centralizer_numpy = not symplectic_product_matrix(vector[None, :], code.H).any()
    in_stabilizer_numpy = rank_np(np.vstack([code.H, vector])) == rank_np(code.H)
    return {
        **bitset,
        "commutes_with_all_stabilizers_numpy": bool(centralizer_numpy),
        "in_stabilizer_row_space_numpy": bool(in_stabilizer_numpy),
        "nontrivial_logical_numpy": bool(centralizer_numpy and not in_stabilizer_numpy),
        "symplectic_weight_numpy": symplectic_weight(vector),
        "matches_expected_distance": bitset["symplectic_weight"] == EXPECTED_DISTANCE,
        "independent_paths_agree": bool(
            centralizer_numpy == bitset["commutes_with_all_stabilizers_bitset"]
            and in_stabilizer_numpy == bitset["in_stabilizer_row_space_bitset"]
            and symplectic_weight(vector) == bitset["symplectic_weight"]
        ),
    }


def sector_path(sector: int) -> Path:
    return STATE_DIR / f"sector_{sector:02d}.json"


def sector_log_path(sector: int) -> Path:
    return STATE_DIR / "traces" / f"sector_{sector:02d}.txt"


def read_valid_sector_record(
    path: Path, expected_identity: dict[str, Any]
) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return None
    for key, value in expected_identity.items():
        if record.get(key) != value:
            return None
    trace_path = record.get("solver", {}).get("trace_path")
    trace_hash = record.get("solver", {}).get("trace_sha256")
    if not trace_path or not trace_hash:
        return None
    absolute_trace = ROOT / trace_path
    if not absolute_trace.exists() or sha256_bytes(absolute_trace.read_bytes()) != trace_hash:
        return None
    return record


def validate_sector_solution(
    code: StabilizerCode,
    logicals: np.ndarray,
    sector: int,
    record: dict[str, Any],
) -> dict[str, Any]:
    """Fail-stop unless a claimed sector solution re-verifies completely.

    A valid solution must have the right shape, symplectic weight at most the
    record's cap, centralizer membership and nontriviality on both GF(2)
    paths with agreement, and pairing exactly 1 with this sector's detector.
    Anything else is a corrupted or defective record and nothing may route
    until it is investigated.
    """

    solution = record["solution"]
    vector = np.asarray(solution["vector"], dtype=np.uint8) & 1
    if vector.shape != (2 * code.n,):
        raise RuntimeError(
            f"sector {sector}: claimed solution has shape {vector.shape}"
        )
    checked = witness_with_crosschecks(code, vector)
    detector = detector_for_sector(logicals, sector)
    pairing = int(vector @ detector % 2)
    if not (
        checked["commutes_with_all_stabilizers_bitset"]
        and checked["commutes_with_all_stabilizers_numpy"]
        and checked["nontrivial_logical_bitset"]
        and checked["nontrivial_logical_numpy"]
        and checked["independent_paths_agree"]
        and checked["symplectic_weight"] <= int(record["weight_cap"])
        and pairing == 1
    ):
        raise RuntimeError(
            f"sector {sector}: claimed solution failed reverification "
            f"(weight={checked['symplectic_weight']}, pairing={pairing}, "
            f"nontrivial_bitset={checked['nontrivial_logical_bitset']}, "
            f"nontrivial_numpy={checked['nontrivial_logical_numpy']})"
        )
    return checked


def record_proof_complete(sector: int, record: dict[str, Any]) -> bool:
    """Accept only the frozen backend/status proof pairs; fail-stop otherwise.

    CP-SAT records may claim completeness only with status INFEASIBLE; PySAT
    records only with status UNSAT.  Any record asserting completeness with a
    different backend/status combination is corrupted or tampered and stops
    the campaign.
    """

    solver = record.get("solver", {})
    claimed = bool(solver.get("proof_complete_for_sector", False))
    if not claimed:
        return False
    name = str(solver.get("name", ""))
    status = str(solver.get("status", ""))
    if name.startswith("OR-Tools CP-SAT") and status == "INFEASIBLE":
        return True
    if name.startswith("PySAT") and status == "UNSAT":
        return True
    raise RuntimeError(
        f"sector {sector}: record claims proof completeness with "
        f"solver={name!r} status={status!r}, outside the frozen backend/"
        "status pairs - nothing may route until this is investigated"
    )


def collect_sector_records(
    code: StabilizerCode,
    logicals: np.ndarray,
    *,
    base_seed: int,
    weight_cap: int,
) -> tuple[list[dict[str, Any]], list[int], list[int]]:
    """Collect the strongest valid record per sector from either backend.

    Amendment (2026-08-15, frozen before any completed sector record): a
    sector is covered iff CP-SAT proves INFEASIBLE or the independent PySAT
    backend proves UNSAT on the equisatisfiable CNF; both records carry the
    same identity binding and hashed traces.
    """

    records: list[dict[str, Any]] = []
    missing: list[int] = []
    unresolved: list[int] = []
    for sector in range(len(ORBIT_REPRESENTATIVE_MASKS)):
        identity = sector_identity(code, logicals, sector, base_seed, weight_cap)
        cpsat_record = read_valid_sector_record(sector_path(sector), identity)
        pysat_record = read_valid_sector_record(
            pysat_sector_path(sector), identity
        )
        proof_records = [
            record
            for record in (cpsat_record, pysat_record)
            if record is not None and record_proof_complete(sector, record)
        ]
        solution_records = [
            record
            for record in (cpsat_record, pysat_record)
            if record is not None and record.get("solution") is not None
        ]
        if proof_records and solution_records:
            checked = witness_with_crosschecks(
                code,
                np.asarray(
                    solution_records[0]["solution"]["vector"], dtype=np.uint8
                ),
            )
            raise RuntimeError(
                f"backend contradiction on sector {sector}: a proof-complete "
                "exclusion coexists with a claimed solution "
                f"(reverification: weight={checked['symplectic_weight']}, "
                f"nontrivial_bitset={checked['nontrivial_logical_bitset']}, "
                f"nontrivial_numpy={checked['nontrivial_logical_numpy']}); "
                "one backend or record is defective - nothing may route until "
                "this is resolved"
            )
        for record in solution_records:
            # Fail-stop on any malformed/tampered claimed solution before it
            # can influence ranking or assembly.
            validate_sector_solution(code, logicals, sector, record)
        if proof_records:
            records.append(proof_records[0])
        elif solution_records:
            # Every validated witness reaches assembly so the upper bound
            # tightens to the minimum recomputed weight.
            records.extend(solution_records)
            unresolved.append(sector)
        elif cpsat_record is not None or pysat_record is not None:
            records.append(cpsat_record if cpsat_record is not None else pysat_record)
            unresolved.append(sector)
        else:
            missing.append(sector)
    return records, missing, unresolved


def next_command(
    unresolved: list[int], *, base_seed: int, time_limit_s: float, workers: int, memory_mb: int
) -> str | None:
    if not unresolved:
        return None
    sector_args = " ".join(str(sector) for sector in unresolved)
    return (
        "uv run python experiments/exp035_pbb_exact_distance.py run "
        f"--sectors {sector_args} --seed {base_seed} --time-limit-s {time_limit_s:g} "
        f"--workers {workers} --memory-mb {memory_mb}"
    )


def assemble_result(
    code: StabilizerCode,
    metadata: dict[str, Any],
    logicals: np.ndarray,
    logical_meta: dict[str, Any],
    records: list[dict[str, Any]],
    missing: list[int],
    unresolved: list[int],
    *,
    base_seed: int,
    time_limit_s: float,
    workers: int,
    memory_mb: int,
) -> dict[str, Any]:
    witness = witness_with_crosschecks(code, default_witness(code))
    orbit_proof = translation_orbit_proof(code, logicals)
    small_weight_proof = small_weight_centralizer_exclusion(code)
    weight_six_proof = weight_six_zero_syndrome_classification(code)
    covered = sorted(
        {
            record["sector"]
            for record in records
            if record_proof_complete(record["sector"], record)
        }
    )
    pending = sorted(set(missing + unresolved))
    covered_coordinates: list[int] = []
    for representative in covered:
        for translation in orbit_proof["representatives"][representative][
            "translations"
        ]:
            covered_coordinates.append(
                int(translation["quotient_coefficient_mask_hex"], 16)
            )
    covered_rows = np.asarray(
        [
            [
                (coordinate >> index) & 1
                for index in range(logicals.shape[0])
            ]
            for coordinate in covered_coordinates
        ],
        dtype=np.uint8,
    )
    if not covered_coordinates:
        covered_rows = np.zeros((0, logicals.shape[0]), dtype=np.uint8)
    covered_rank_numpy = rank_np(covered_rows)
    covered_rank_bitset = rank_bitset(
        covered_coordinates, logicals.shape[0]
    )
    lower_complete = bool(
        orbit_proof["selected_orbits_span_full_dual"]
        and covered_rank_numpy == logicals.shape[0]
        and covered_rank_bitset == logicals.shape[0]
    )
    upper_complete = bool(
        witness["symplectic_weight"] == EXPECTED_DISTANCE
        and witness["nontrivial_logical_bitset"]
        and witness["nontrivial_logical_numpy"]
        and witness["independent_paths_agree"]
    )
    weight6_clean = bool(
        weight_six_proof["complete"]
        and small_weight_proof["no_centralizer_of_weight_at_most_5"]
    )
    candidate_logicals: list[dict[str, Any]] = []
    if weight6_clean and not weight_six_proof["no_weight6_logical"]:
        candidate_logicals.extend(weight_six_proof["nonmember_logicals"])
    for record in records:
        solution = record.get("solution")
        if solution is None:
            continue
        checked = witness_with_crosschecks(
            code, np.asarray(solution["vector"], dtype=np.uint8)
        )
        candidate_logicals.append(checked)
    verified_low_weight = [
        candidate
        for candidate in candidate_logicals
        if candidate["nontrivial_logical_bitset"]
        and candidate["nontrivial_logical_numpy"]
        and candidate["independent_paths_agree"]
        and candidate["symplectic_weight"] <= LOWER_EXCLUSION_WEIGHT
    ]
    counterexample = (
        min(verified_low_weight, key=lambda c: c["symplectic_weight"])
        if verified_low_weight
        else None
    )
    if counterexample is not None and lower_complete:
        raise RuntimeError(
            "inconsistent evidence: complete sector exclusion coexists with a "
            f"verified weight-{counterexample['symplectic_weight']} logical"
        )
    lower_bound = (
        EXPECTED_DISTANCE
        if lower_complete
        else (
            weight_six_proof["certified_lower_bound_with_weight5_exclusion"]
            if weight6_clean and weight_six_proof["no_weight6_logical"]
            else small_weight_proof["certified_distance_lower_bound"]
        )
    )
    if counterexample is not None:
        # A verified nontrivial logical of weight w <= 11 caps the distance
        # at w, contradicting the untrusted catalogue d=12.  Exactness holds
        # when the exclusion floor meets the witness weight.
        upper_bound = counterexample["symplectic_weight"]
        if lower_bound > upper_bound:
            raise RuntimeError(
                "inconsistent evidence: exclusion floor exceeds a verified "
                "witness weight"
            )
        exact = lower_bound == upper_bound
    else:
        exact = lower_complete and upper_complete
        upper_bound = EXPECTED_DISTANCE if upper_complete else None
    statuses: dict[str, str] = {}
    for record in records:
        statuses.setdefault(str(record["sector"]), record["solver"]["status"])
    result = {
        "schema": "pbb-exact-distance-certificate-v3",
        "code_id": CODE_ID,
        **metadata,
        "logical_quotient": logical_meta,
        "translation_orbit_reduction": orbit_proof,
        "small_weight_exclusion": small_weight_proof,
        "weight_six_classification": weight_six_proof,
        "distance_problem": {
            "definition": "min symplectic weight over S^perp \\ S",
            "nontriviality": (
                "For v in S^perp, v is outside S iff its pairing functional "
                "on S^perp/S is nonzero. A nonzero functional evaluates to 1 "
                "on some member of any spanning family; the selected translated "
                "representatives span the full 24-dimensional quotient."
            ),
            "lower_exclusion_weight": LOWER_EXCLUSION_WEIGHT,
            "full_logical_dual_dimension": logicals.shape[0],
            "orbit_representative_count": len(ORBIT_REPRESENTATIVE_MASKS),
            "covered_infeasible_representatives": covered,
            "covered_dual_rank_numpy": covered_rank_numpy,
            "covered_dual_rank_bitset": covered_rank_bitset,
            "missing_representatives": missing,
            "unresolved_representatives": unresolved,
            "representative_statuses": statuses,
            "lower_bound_coverage_complete": lower_complete,
            "timeout_policy": (
                "UNKNOWN/MODEL_INVALID/timeout remains partial and never exact"
            ),
        },
        "solver_protocol": {
            "solver": "OR-Tools CP-SAT exact Boolean/XOR feasibility",
            "version": _ortools_version(),
            "independent_backends": [
                "OR-Tools CP-SAT (proof-complete status INFEASIBLE)",
                "PySAT cadical195 on the equisatisfiable Tseitin-XOR + "
                "seqcounter CNF (proof-complete status UNSAT); frozen "
                "amendment 2026-08-15, before any completed sector record "
                "existed; SAT models re-verified through both GF(2) paths; "
                "cross-backend contradictions fail-stop the campaign",
            ],
            "base_seed": base_seed,
            "sector_seed_derivation": (
                "low 31 bits of "
                "SHA256(code_id:orbit-representative:base_seed:sector)"
            ),
            "time_limit_s_per_sector": time_limit_s,
            "workers_per_sector": workers,
            "max_memory_mb_per_sector": memory_mb,
            "state_directory": str(STATE_DIR.relative_to(ROOT, walk_up=True)),
            "atomic_sector_records": True,
            "trace_per_sector": True,
            "canonical_gate": (
                "matching proof-complete representative records (CP-SAT "
                "INFEASIBLE or PySAT UNSAT on the equisatisfiable CNF) whose "
                "transported dual rank is 24 in both GF(2) paths, plus "
                "independently checked weight-12 witness"
            ),
        },
        "witness": witness,
        "unexpected_low_weight_logical": counterexample,
        "bounds": {
            "certified_lower_bound": lower_bound,
            "certified_upper_bound": upper_bound,
            "distance": lower_bound if exact else None,
            "CERTIFIED_EXACT": exact,
        },
        "resume_command": next_command(
            pending,
            base_seed=base_seed,
            time_limit_s=time_limit_s,
            workers=workers,
            memory_mb=memory_mb,
        ),
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
    }
    return result



def reverify_low_weight_counterexample(candidate: dict[str, Any]) -> bool:
    """Re-derive every claim of a sub-12 counterexample against fresh algebra.

    The router never trusts upstream fields: the target is rebuilt from the
    catalogue and the candidate vector is re-crosschecked through both GF(2)
    paths.  Only a genuinely nontrivial logical whose recomputed weight
    equals the claimed weight and is at most the exclusion cap passes.
    """

    try:
        claimed_weight = int(candidate["symplectic_weight"])
        row, raw = load_target_catalogue_row()
        code, _ = rebuild_metadata(row, raw)
        vector = np.zeros(2 * code.n, dtype=np.uint8)
        for coordinate, value in enumerate(
            list(candidate["x"]) + list(candidate["z"])
        ):
            vector[coordinate] = int(value) & 1
        checked = witness_with_crosschecks(code, vector)
    except Exception:
        return False
    return bool(
        checked["symplectic_weight"] == claimed_weight
        and checked["symplectic_weight_numpy"] == claimed_weight
        and checked["symplectic_weight"] <= LOWER_EXCLUSION_WEIGHT
        and checked["nontrivial_logical_bitset"]
        and checked["nontrivial_logical_numpy"]
        and checked["independent_paths_agree"]
        and checked["support"] == list(candidate["support"])
    )

def route_result(result: dict[str, Any]) -> tuple[str, Path]:
    rebuild = result["rebuild"]
    witness = result["witness"]
    orbit = result["translation_orbit_reduction"]
    small_weight = result["small_weight_exclusion"]
    weight_six = result["weight_six_classification"]
    clean = bool(
        rebuild["n"] == EXPECTED_N
        and rebuild["rank_numpy"] == EXPECTED_RANK
        and rebuild["rank_bitset"] == EXPECTED_RANK
        and rebuild["k"] == EXPECTED_K
        and rebuild["commutes_numpy"]
        and rebuild["commutes_bitset"]
        and witness["independent_paths_agree"]
        and orbit["every_translation_preserves_stabilizer_rowspace_bitset"]
        and orbit["every_translation_preserves_logical_pairing_numpy"]
        and orbit["selected_orbits_span_full_dual"]
        and small_weight["no_centralizer_of_weight_at_most_5"]
        and weight_six["complete"]
    )
    candidate = result.get("unexpected_low_weight_logical")
    if candidate is not None:
        # Decisive counterexample to the untrusted catalogue distance; keep it
        # out of the pre-registered d=12 canonical slot but persist it as its
        # own first-class certificate for protocol amendment.  The claim is
        # re-derived from fresh algebra; anything unverifiable is quarantined.
        if clean and reverify_low_weight_counterexample(candidate):
            return "counterexample", COUNTEREXAMPLE
        return "quarantine", QUARANTINE
    if result["bounds"]["CERTIFIED_EXACT"] and clean:
        return "canonical", CANONICAL
    if clean:
        return "partial", PARTIAL_SUMMARY
    return "quarantine", QUARANTINE


def write_concise_report(
    result: dict[str, Any], records: list[dict[str, Any]]
) -> None:
    exact = result["bounds"]["CERTIFIED_EXACT"]
    bounds = result["bounds"]
    report = {
        "experiment": "EXP-035",
        "code_id": CODE_ID,
        "verdict": (
            f"EXACT d={bounds['distance']}"
            if exact
            else (
                f"PARTIAL {bounds['certified_lower_bound']} <= d <= "
                f"{bounds['certified_upper_bound']}"
            )
        ),
        "artifact_route": result["artifact_route"],
        "artifact": result["artifact_path"],
        "bounds": bounds,
        "representative_statuses": result["distance_problem"][
            "representative_statuses"
        ],
        "missing_representatives": result["distance_problem"][
            "missing_representatives"
        ],
        "unresolved_representatives": result["distance_problem"][
            "unresolved_representatives"
        ],
        "covered_dual_rank_numpy": result["distance_problem"][
            "covered_dual_rank_numpy"
        ],
        "covered_dual_rank_bitset": result["distance_problem"][
            "covered_dual_rank_bitset"
        ],
        "solver_wall_time_s": sum(
            record.get("solver", {}).get("wall_time_s", 0.0)
            for record in records
        ),
        "witness_support": result["witness"]["support"],
        "witness_verified_nontrivial": bool(
            result["witness"]["nontrivial_logical_bitset"]
            and result["witness"]["nontrivial_logical_numpy"]
        ),
        "translation_orbit_union_rank": result["translation_orbit_reduction"][
            "selected_orbit_union_rank_bitset"
        ],
        "resume_command": result["resume_command"],
    }
    atomic_write_json(REPORT if exact else PARTIAL_REPORT, report)


def run_selected_sectors(args: argparse.Namespace) -> int:
    row, raw = load_target_catalogue_row()
    code, metadata = rebuild_metadata(row, raw)
    logicals, logical_meta = logical_basis_metadata(code)
    selected = (
        args.sectors
        if args.sectors is not None
        else list(range(len(ORBIT_REPRESENTATIVE_MASKS)))
    )
    if len(set(selected)) != len(selected):
        raise ValueError("duplicate sector arguments")
    for sector in selected:
        identity = sector_identity(code, logicals, sector, args.seed, LOWER_EXCLUSION_WEIGHT)
        existing = read_valid_sector_record(sector_path(sector), identity)
        if (
            existing is not None
            and existing["solver"]["status"] == "INFEASIBLE"
            and record_proof_complete(sector, existing)
        ):
            print(f"sector {sector:02d}: resume hit INFEASIBLE")
            continue
        record = solve_lower_sector(
            code,
            logicals,
            sector,
            base_seed=args.seed,
            weight_cap=LOWER_EXCLUSION_WEIGHT,
            time_limit_s=args.time_limit_s,
            workers=args.workers,
            memory_mb=args.memory_mb,
            log_path=sector_log_path(sector),
        )
        atomic_write_json(sector_path(sector), record)
        print(
            f"sector {sector:02d}: {record['solver']['status']} "
            f"in {record['solver']['wall_time_s']:.3f}s"
        )
    records, missing, unresolved = collect_sector_records(
        code,
        logicals,
        base_seed=args.seed,
        weight_cap=LOWER_EXCLUSION_WEIGHT,
    )
    result = assemble_result(
        code,
        metadata,
        logicals,
        logical_meta,
        records,
        missing,
        unresolved,
        base_seed=args.seed,
        time_limit_s=args.time_limit_s,
        workers=args.workers,
        memory_mb=args.memory_mb,
    )
    route, path = route_result(result)
    result["artifact_route"] = route
    result["artifact_path"] = str(path.relative_to(ROOT))
    atomic_write_json(path, result)
    write_concise_report(result, records)
    print(json.dumps({"route": route, "path": str(path), "bounds": result["bounds"], "unresolved": result["distance_problem"]["unresolved_representatives"], "missing": missing, "resume_command": result["resume_command"]}, indent=2))
    return 0 if route == "canonical" else 2


def campaign_caps() -> tuple[float, int, int]:
    """Inherit solver caps from the live campaign artifact.

    The partial summary (or, once routed, the canonical artifact) is the
    single source of truth for the running campaign's per-sector caps; the
    module defaults apply only before any artifact exists.
    """

    for artifact_path in (CANONICAL, PARTIAL_SUMMARY):
        if not artifact_path.exists():
            continue
        protocol = json.loads(artifact_path.read_text())["solver_protocol"]
        time_limit_s = float(protocol["time_limit_s_per_sector"])
        workers = int(protocol["workers_per_sector"])
        memory_mb = int(protocol["max_memory_mb_per_sector"])
        if time_limit_s <= 0 or workers <= 0 or memory_mb <= 0:
            raise RuntimeError(
                f"invalid solver caps in {artifact_path.name}: "
                f"{time_limit_s=} {workers=} {memory_mb=}"
            )
        return time_limit_s, workers, memory_mb
    return DEFAULT_TIME_LIMIT_S, DEFAULT_WORKERS, DEFAULT_MEMORY_MB


def run_selected_sectors_pysat(args: argparse.Namespace) -> int:
    """Independent-backend sector runs; identical identity/routing discipline."""

    row, raw = load_target_catalogue_row()
    code, metadata = rebuild_metadata(row, raw)
    logicals, logical_meta = logical_basis_metadata(code)
    selected = (
        args.sectors
        if args.sectors is not None
        else list(range(len(ORBIT_REPRESENTATIVE_MASKS)))
    )
    if len(set(selected)) != len(selected):
        raise ValueError("duplicate sector arguments")
    for sector in selected:
        identity = sector_identity(
            code, logicals, sector, args.seed, LOWER_EXCLUSION_WEIGHT
        )
        for existing_path in (sector_path(sector), pysat_sector_path(sector)):
            existing = read_valid_sector_record(existing_path, identity)
            if existing is not None and record_proof_complete(
                sector, existing
            ):
                print(f"sector {sector:02d}: resume hit proof-complete record")
                break
        else:
            record = solve_lower_sector_pysat(
                code,
                logicals,
                sector,
                base_seed=args.seed,
                weight_cap=LOWER_EXCLUSION_WEIGHT,
                conflict_budget=args.conflict_budget,
                log_path=pysat_sector_log_path(sector),
                solver_name=args.solver,
            )
            solution = record.get("solution")
            if solution is not None:
                record["solution_crosschecks"] = validate_sector_solution(
                    code, logicals, sector, record
                )
            atomic_write_json(pysat_sector_path(sector), record)
            print(
                f"sector {sector:02d}: pysat {record['solver']['status']} "
                f"in {record['solver']['wall_time_s']:.3f}s"
            )
    inherited_time_limit_s, inherited_workers, inherited_memory_mb = (
        campaign_caps()
    )
    records, missing, unresolved = collect_sector_records(
        code,
        logicals,
        base_seed=args.seed,
        weight_cap=LOWER_EXCLUSION_WEIGHT,
    )
    result = assemble_result(
        code,
        metadata,
        logicals,
        logical_meta,
        records,
        missing,
        unresolved,
        base_seed=args.seed,
        time_limit_s=inherited_time_limit_s,
        workers=inherited_workers,
        memory_mb=inherited_memory_mb,
    )
    route, path = route_result(result)
    result["artifact_route"] = route
    result["artifact_path"] = str(path.relative_to(ROOT))
    atomic_write_json(path, result)
    write_concise_report(result, records)
    print(json.dumps({
        "route": route,
        "path": str(path),
        "bounds": result["bounds"],
        "missing": missing,
        "unresolved": unresolved,
    }, indent=2))
    return 0 if route == "canonical" else 2


def assemble_only(args: argparse.Namespace) -> int:
    row, raw = load_target_catalogue_row()
    code, metadata = rebuild_metadata(row, raw)
    logicals, logical_meta = logical_basis_metadata(code)
    records, missing, unresolved = collect_sector_records(
        code, logicals, base_seed=args.seed, weight_cap=LOWER_EXCLUSION_WEIGHT
    )
    result = assemble_result(
        code,
        metadata,
        logicals,
        logical_meta,
        records,
        missing,
        unresolved,
        base_seed=args.seed,
        time_limit_s=args.time_limit_s,
        workers=args.workers,
        memory_mb=args.memory_mb,
    )
    route, path = route_result(result)
    result["artifact_route"] = route
    result["artifact_path"] = str(path.relative_to(ROOT))
    atomic_write_json(path, result)
    write_concise_report(result, records)
    print(json.dumps({"route": route, "path": str(path), "bounds": result["bounds"], "missing": missing, "unresolved": unresolved, "resume_command": result["resume_command"]}, indent=2))
    return 0 if route == "canonical" else 2


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("run", "assemble"):
        sub = subparsers.add_parser(command)
        sub.add_argument("--seed", type=int, default=DEFAULT_SEED)
        sub.add_argument("--time-limit-s", type=float, default=DEFAULT_TIME_LIMIT_S)
        sub.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
        sub.add_argument("--memory-mb", type=int, default=DEFAULT_MEMORY_MB)
        if command == "run":
            sub.add_argument("--sectors", type=int, nargs="*")
    pysat_sub = subparsers.add_parser("run-pysat")
    pysat_sub.add_argument("--seed", type=int, default=DEFAULT_SEED)
    pysat_sub.add_argument("--sectors", type=int, nargs="*")
    pysat_sub.add_argument("--conflict-budget", type=int, default=0,
                           help="0 = run to completion")
    pysat_sub.add_argument("--solver", type=str, default="cadical195")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.command == "run-pysat":
        if args.conflict_budget < 0:
            raise ValueError("conflict budget must be nonnegative")
        return run_selected_sectors_pysat(args)
    if args.time_limit_s <= 0 or args.workers <= 0 or args.memory_mb <= 0:
        raise ValueError("solver caps must be positive")
    if args.command == "run":
        return run_selected_sectors(args)
    return assemble_only(args)


if __name__ == "__main__":
    raise SystemExit(main())
