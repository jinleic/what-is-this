"""Control tests: our independent construction vs the qcode-discovery convention.

qcode-discovery builds monomial matrices as ``bb_code.eval(poly).lift().T``,
i.e. the transpose of the regular representation.  Ours uses the regular
representation directly.  These differ by the ring involution x^a y^b ->
x^-a y^-b, which is realised on qubits by the permutation (i,j) -> (-i,-j).

This module proves the two conventions give permutation-equivalent codes, so
any distance / k / weight statement transfers.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from qec_research.codes.bicycle import (  # noqa: E402
    BRAVYI_BB, PBBSpec, build_pbb, monomial_matrix, poly_matrix,
    poly_transpose_terms, bb_stabilizer,
)
from qec_research.gf2.linalg import rank_np  # noqa: E402
from qec_research.symplectic.core import StabilizerCode  # noqa: E402


def transposed_convention_pbb(spec: PBBSpec) -> np.ndarray:
    """Build H exactly as qcode-discovery does: every block matrix transposed."""
    ell, m = spec.ell, spec.m
    dim = ell * m
    A = poly_matrix(ell, m, spec.A).T
    B = poly_matrix(ell, m, spec.B).T
    C = poly_matrix(ell, m, spec.C).T if spec.C else np.zeros((dim, dim), np.uint8)
    D = poly_matrix(ell, m, spec.D).T if spec.D else np.zeros((dim, dim), np.uint8)
    zero = np.zeros((dim, dim), dtype=np.uint8)
    top = np.hstack([A, B, C, D])
    bot = np.hstack([zero, zero, B.T, A.T])
    return np.vstack([top, bot]).astype(np.uint8)


def inversion_permutation(ell: int, m: int) -> np.ndarray:
    """Index map for (i,j) -> (-i,-j) on a single ell x m block."""
    p = np.empty(ell * m, dtype=int)
    for i in range(ell):
        for j in range(m):
            p[i * m + j] = ((-i) % ell) * m + ((-j) % m)
    return p


@pytest.mark.parametrize("terms", [[(3, 0), (0, 1), (0, 2)], [(2, 1), (3, 1), (4, 4)]])
def test_transpose_equals_negated_terms(terms):
    ell, m = 12, 6
    lhs = poly_matrix(ell, m, terms).T
    rhs = poly_matrix(ell, m, poly_transpose_terms(ell, m, terms))
    assert np.array_equal(lhs, rhs)


def test_conventions_are_permutation_equivalent():
    """Our H, with qubits relabelled by (i,j)->(-i,-j), spans their row space."""
    spec = PBBSpec(12, 6, [(2, 1), (3, 1), (4, 4)], [(0, 0), (5, 1), (4, 5)],
                   [(3, 4), (4, 1)], [(2, 5), (3, 2)])
    mine = build_pbb(spec).H
    theirs = transposed_convention_pbb(spec)
    ell, m, dim = spec.ell, spec.m, spec.ell * spec.m
    p = inversion_permutation(ell, m)
    full = np.concatenate([p, p + dim])           # qubit permutation on 2*dim qubits
    n = 2 * dim
    permuted = np.hstack([mine[:, :n][:, full], mine[:, n:][:, full]])
    r1, r2 = rank_np(permuted), rank_np(theirs)
    rj = rank_np(np.vstack([permuted, theirs]))
    assert r1 == r2 == rj, f"row spaces differ: {r1} {r2} {rj}"


def test_bb_reference_parameters():
    expect = {"[[72,12,6]]": (72, 12), "[[90,8,10]]": (90, 8), "[[108,8,10]]": (108, 8),
              "[[144,12,12]]": (144, 12), "[[288,12,18]]": (288, 12),
              "[[360,12,<=24]]": (360, 12)}
    for name, (n, k) in expect.items():
        sc = bb_stabilizer(BRAVYI_BB[name])
        v = sc.validate()
        assert v["n"] == n and v["k"] == k, (name, v["n"], v["k"])
        assert v["commutes_numpy"] and v["commutes_bitset"]
        assert v["rank_agree"]
        assert v["max_check_weight"] == 6
        assert not sc.is_direct_sum()


def test_gf2_paths_agree_on_random_matrices():
    from qec_research.gf2.linalg import rank_bitset, rows_to_bitsets, nullspace_np, matmul
    rng = np.random.default_rng(0xC0DE)
    for _ in range(40):
        r, c = rng.integers(1, 40, size=2)
        M = rng.integers(0, 2, size=(int(r), int(c))).astype(np.uint8)
        assert rank_np(M) == rank_bitset(rows_to_bitsets(M), int(c))
        N = nullspace_np(M)
        if N.shape[0]:
            assert not matmul(M, N.T).any()
        assert rank_np(M) + N.shape[0] == int(c)


def test_distance_is_not_declared_exact_without_a_witness():
    """Regression (advisory): the solver used to seed `best = upper_bound` and
    search only strictly below it, so a run in which every sector was
    infeasible returned the seed as `exact=True` with no witness.  That proves
    d >= D, never d = D.  Here we hand it a cap BELOW the true distance, so no
    witness can exist, and require that it refuses to claim exactness."""
    from qec_research.codes.bicycle import BRAVYI_BB, build_bb
    from qec_research.distance.exact import exact_distance_css

    HX, HZ = build_bb(BRAVYI_BB["[[72,12,6]]"])      # true distance 6
    res = exact_distance_css(HX, HZ, time_limit_s=120, workers=4, upper_bound=4)
    assert res["d_X"] is None and res["d_Z"] is None, "no operator of weight <= 4 exists"
    assert res["d_exact"] is False, "must not claim exactness without a witness"
    assert res["d_X_lower_bound"] >= 5 and res["d_Z_lower_bound"] >= 5
    assert res["d_X_all_sectors_decided"] and res["d_Z_all_sectors_decided"]

    # with an adequate cap it must produce a witness of exactly the right weight
    ok = exact_distance_css(HX, HZ, time_limit_s=120, workers=4, upper_bound=6)
    assert ok["d"] == 6 and ok["d_exact"] is True
    for tag in ("X", "Z"):
        w = ok[f"d_{tag}_witness"]
        assert w is not None and len(w) == ok[f"d_{tag}"]


def test_sector_solver_also_refuses_exactness_without_a_witness():
    """Same regression as above, for the sector-restricted solver that feeds
    the Theorem-1/Theorem-2 checks."""
    from qec_research.codes.bicycle import BRAVYI_BB, bb_stabilizer
    from qec_research.distance.sectors import min_pure_z_logical

    code = bb_stabilizer(BRAVYI_BB["[[72,12,6]]"])   # true pure-Z minimum is 6
    r = min_pure_z_logical(code, time_limit_s=120, workers=4, upper_bound=4)
    assert r.weight is None, "no pure-Z logical of weight <= 4 exists"
    assert r.exact is False, "must not claim exactness without a witness"

    ok = min_pure_z_logical(code, time_limit_s=120, workers=4, upper_bound=6)
    assert ok.weight == 6 and ok.exact is True
    assert ok.support is not None and len(ok.support) == 6
    v = np.zeros(2 * code.n, dtype=np.uint8)
    for j in ok.support:
        v[code.n + j] = 1
    assert code.is_logical(v)


def test_parity_forwarding_beats_the_qubit_degree_bound():
    """Guards the scope of Theorem C4 (see proofs/pbb_structure.md, open Q4).

    A counting argument would make the depth no-go unconditional via
    T >= max_j deg(j), derived from "each check-qubit incidence needs its own
    data gate".  Ancilla-ancilla CNOTs refute that PREMISE: this circuit measures
    Z0Z1 and Z0Z2 with ONE gate at the shared qubit q0, where deg(q0) = 2.
    The inequality itself is NOT refuted (the circuit uses 3 layers >= 2) and
    remains open for general multi-ancilla circuits."""
    import itertools
    import stim

    c = stim.Circuit()
    c.append("R", [3, 4])
    c.append("CX", [0, 3])      # a1 ^= q0     <- the only gate touching q0
    c.append("CX", [3, 4])      # a2 ^= a1     <- parity forwarding
    c.append("CX", [1, 3])      # a1 = q0+q1
    c.append("CX", [2, 4])      # a2 = q0+q2
    c.append("M", [3, 4])

    for bits in itertools.product([0, 1], repeat=3):
        pre = stim.Circuit()
        for i, b in enumerate(bits):
            if b:
                pre.append("X", [i])
        sim = stim.TableauSimulator()
        sim.do(pre + c)
        m1, m2 = (int(x) for x in sim.current_measurement_record()[-2:])
        assert m1 == (bits[0] ^ bits[1]), "a1 must measure Z0Z1"
        assert m2 == (bits[0] ^ bits[2]), "a2 must measure Z0Z2"

    touching_q0 = sum(1 for inst in c
                      if inst.name == "CX" and 0 in [t.value for t in inst.targets_copy()])
    assert touching_q0 == 1, (
        "one gate at q0 while deg(q0)=2 -- this refutes the ONE-GATE-PER-INCIDENCE "
        "premise used to derive T >= max_j deg(j)")

    # It does NOT refute the inequality itself: the circuit still uses 3 two-qubit
    # layers against max degree 2.  Whether T >= max_j deg(j) holds for general
    # multi-ancilla circuits is OPEN.
    two_qubit_layers = 3          # q0->a1 | a1->a2 | (q1->a1, q2->a2)
    max_degree = 2                # q0 lies in both checks
    assert two_qubit_layers >= max_degree, (
        "the forwarding circuit must still satisfy T >= max deg -- it refutes the "
        "derivation, not the bound")


def test_feasibility_query_agrees_with_the_reference_symplectic_solver():
    """Regression (advisory): EXP-019's feasibility query originally used
    `constraints = code.H` instead of `lambda_swap(code.H)`.  Ordinary-dot
    constraints are NOT symplectic centralizer membership, and on mixed
    (non-CSS) checks they admit vectors that are not logicals, which would make
    a `DOMINATED` verdict unsound.  Cross-check the query against the reference
    solver on a genuinely non-CSS code."""
    import importlib.util
    from pathlib import Path

    from qec_research.codes.bicycle import PBBSpec, build_pbb
    from qec_research.distance.exact import exact_distance_symplectic

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "e19", str(root / "experiments" / "exp019_heldout_domination.py"))
    e19 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(e19)

    # a small genuinely non-CSS PBB code: Bravyi's [[72,12,6]] parent with a
    # commutation-valid weight-2 perturbation (from the exact null space)
    Q = build_pbb(PBBSpec(ell=6, m=6, A=[(3, 0), (0, 1), (0, 2)],
                          B=[(0, 3), (1, 0), (2, 0)],
                          C=[(0, 0), (0, 3)], D=[]))
    # a stabiliser is MIXED when the row has both X support and Z support --
    # not when a single qubit carries Y
    mixed = Q.H[:, :Q.n].any(axis=1) & Q.H[:, Q.n:].any(axis=1)
    assert mixed.any(), "the test code must actually have mixed stabilisers"
    ref = exact_distance_symplectic(Q, time_limit_s=180, workers=4)
    assert ref.exact and ref.value is not None
    d = ref.value

    found, w, decided = e19.exists_logical_within(Q, d, tl=180, workers=4)
    assert found and w is not None and w <= d, "a logical of weight <= d must exist"

    found2, _, decided2 = e19.exists_logical_within(Q, d - 1, tl=180, workers=4)
    assert not found2 and decided2, "nothing below the true distance may be found"


def test_artifact_routing():
    """Canonical writes require clean AND full coverage (FR-012 regression).

    A passing smoke run must never overwrite the canonical artifact, and a
    failing run must be quarantined, not filed with legitimate partial runs.
    """
    from qec_research.artifacts import canonical_route
    assert canonical_route(True, True) == "canonical"
    assert canonical_route(True, False) == "partial_runs"   # passing 60-row smoke
    assert canonical_route(False, True) == "quarantine"     # failing full run
    assert canonical_route(False, False) == "quarantine"    # failing smoke
