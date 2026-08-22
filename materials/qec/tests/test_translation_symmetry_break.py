"""Validate the premise behind the BB translation symmetry break.

The break asserts: for a bivariate-bicycle-derived code with block length
lm, the translation group {x^a y^b} acts on qubits as permutations that

  (1) preserve the stabilizer row space (hence the centralizer, hence
      nontriviality of a logical) and preserve weight, and
  (2) act transitively on the lm positions inside each block.

Given (1) and (2), any nontrivial logical of weight <= c can be translated so
that qubit 0 of block A or qubit 0 of block B carries support, so adding the
single clause "qubit 0 of block A supported OR qubit 0 of block B supported"
cannot turn a satisfiable instance unsatisfiable.  Both facts are checked
here algebraically, plus a differential SAT/UNSAT agreement test.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

from qec_research.codes.bicycle import (  # noqa: E402
    BRAVYI_BB,
    build_bb,
    monomial_matrix,
)
from qec_research.distance.sat_decide import (  # noqa: E402
    build_decision_cnf,
    css_side_instance,
    symplectic_instance,
)
from qec_research.gf2.linalg import rank_bitset, rows_to_bitsets  # noqa: E402

_SPEC = importlib.util.spec_from_file_location(
    "exp027_delta_audit", ROOT / "experiments" / "exp027_delta_audit.py"
)
assert _SPEC and _SPEC.loader
E27 = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(E27)


def _row_space_equal(A: np.ndarray, B: np.ndarray) -> bool:
    """Compare GF(2) row spaces of two equal-width matrices."""
    a = rows_to_bitsets(A)
    b = rows_to_bitsets(B)
    width = A.shape[1]
    ra, rb = rank_bitset(a, width), rank_bitset(b, width)
    return ra == rb == rank_bitset(a + b, width)


def _qubit_permutation(ell: int, m: int, a: int, b: int) -> np.ndarray:
    """Permutation of the 2*lm qubits induced by x^a y^b on both blocks."""
    P = monomial_matrix(ell, m, a, b)
    lm = ell * m
    perm = np.empty(2 * lm, dtype=int)
    for src in range(lm):
        dst = int(np.flatnonzero(P[src])[0])
        perm[src] = dst
        perm[lm + src] = lm + dst
    return perm


def test_translation_is_code_automorphism_and_transitive() -> None:
    spec = BRAVYI_BB["[[144,12,12]]"]
    ell, m = spec.ell, spec.m
    HX, HZ = build_bb(spec)
    lm = ell * m
    for a, b in ((1, 0), (0, 1), (3, 2)):
        perm = _qubit_permutation(ell, m, a, b)
        # (1) row spaces preserved on both CSS sides
        assert _row_space_equal(HX, HX[:, perm])
        assert _row_space_equal(HZ, HZ[:, perm])
    # (2) transitivity inside a block: the orbit of position 0 is everything
    reached = {0}
    frontier = [0]
    gens = [_qubit_permutation(ell, m, 1, 0), _qubit_permutation(ell, m, 0, 1)]
    while frontier:
        current = frontier.pop()
        for g in gens:
            nxt = int(g[current])
            if nxt not in reached:
                reached.add(nxt)
                frontier.append(nxt)
    assert reached == set(range(lm))


def test_pbb_translation_preserves_stabilizer_row_space() -> None:
    rows = E27.load_catalogue()
    checked = 0
    for index, row in enumerate(rows):
        if 2 * int(row["ell"]) * int(row["m"]) != 72:
            continue
        _, code = E27.pbb_code(row)
        ell, m = int(row["ell"]), int(row["m"])
        lm = ell * m
        perm = _qubit_permutation(ell, m, 1, 1)
        symplectic_perm = np.concatenate([perm, code.n + perm])
        assert _row_space_equal(code.H, code.H[:, symplectic_perm])
        checked += 1
        if checked == 3:
            break
    assert checked == 3


def _status(instance, cap: int) -> str:
    from pysat.solvers import Solver

    cnf = build_decision_cnf(instance, cap)
    with Solver(name="cadical195", bootstrap_with=cnf.clauses) as engine:
        return "SAT" if engine.solve() else "UNSAT"


def test_symmetry_break_never_changes_satisfiability() -> None:
    """Differential check: symmetrized and plain encodings must agree."""
    rows = E27.load_catalogue()
    targets = [
        (index, row)
        for index, row in enumerate(rows)
        if 2 * int(row["ell"]) * int(row["m"]) == 72
    ][:2]
    assert targets
    for index, row in targets:
        _, HX, HZ = E27.parent_matrices(row)
        _, pbb = E27.pbb_code(row)
        block = int(row["ell"]) * int(row["m"])
        cases = [
            (css_side_instance(HX, HZ, "x"), css_side_instance(HX, HZ, "x", block_length=block)),
            (css_side_instance(HX, HZ, "z"), css_side_instance(HX, HZ, "z", block_length=block)),
            (
                symplectic_instance(pbb, pbb.logical_basis()),
                symplectic_instance(pbb, pbb.logical_basis(), block_length=block),
            ),
        ]
        for plain, symmetric in cases:
            assert symmetric.symmetry_clause
            for cap in (3, 4, 5, 6, 7):
                assert _status(plain, cap) == _status(symmetric, cap), (
                    f"row {index} {plain.kind} cap {cap}: symmetry break changed status"
                )
