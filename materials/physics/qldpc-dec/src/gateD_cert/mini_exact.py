"""Surface-code (code-capacity) exact-ML generators for the mini sanity check.

Rotated surface code d x d data qubits; Z-error decoding from the X-check
syndrome (paper Sec 7 preamble: code capacity = independent X-flip rate p
decoded from the Z-syndrome). Conventions:
  H: (n_X checks) x (n data) binary check matrix (X-stabilizer supports).
  L: (1 x n) logical Z equivalent-class representative; the degenerate-ML
     logical observable = L @ e for the Z-error vector e.
  The coset sum per class (Eq 7): Z_lam = sum_u w^{|e_lam xor u S|} over the
  trivial-kernel generator basis S.
"""
from __future__ import annotations

import numpy as np


def rotated_surface_code(d: int):
    """Returns H (n_xchecks, n), L (1, n) for the rotated surface code.

    Qubit coords (r, c) with r, c in [0, d-1], qubit id = r*d + c.
    X-type plaquette checks on faces with even (r+c); measure X on faces
    (Z-error syndrome). Boundary handling per standard rotated layout:
    a face is an X-check if (r + c) % 2 == 1 and it touches >= 2 qubits.
    """
    n = d * d
    checks = []
    for r in range(d - 1):
        for c in range(d - 1):
            if (r + c) % 2 == 1:
                face = [ (r, c), (r + 1, c), (r, c + 1), (r + 1, c + 1) ]
                row = np.zeros(n, dtype=np.uint8)
                for (rr, cc) in face:
                    row[rr * d + cc] = 1
                if row.sum() >= 2:
                    checks.append(row)
    # add the two boundary X-checks (top face r=0 and bottom r=d-1 segments)
    for c in range(0, d - 1):
        row = np.zeros(n, dtype=np.uint8)
        row[c], row[c + 1] = 1, 1
        # keep only boundaries with correct parity class for X checks
        if c % 2 == 1:
            checks.append(row)
        row2 = np.zeros(n, dtype=np.uint8)
        row2[(d - 1) * d + c] = 1
        row2[(d - 1) * d + c + 1] = 1
        if c % 2 == 1:
            checks.append(row2)
    H = np.array(checks, dtype=np.uint8)
    # remove duplicate rows (rank reduction from the naive construction)
    from .gf2 import kernel_basis
    _, ind = np.unique(H, axis=0, return_index=True)
    H = H[np.sort(ind)]
    # logical Z operator: a full column of data qubits (c = 0): commutes
    # with every X-check? X-checks on c=0,c=1 faces (r,c faces with r+c odd,
    # c=0 -> r odd): overlap 1 -> anticommute: wrong; use a row of data
    # qubits (r = 0): X-faces with r=0 need c odd: overlap with row r=0 at
    # (0, c), (0, c+1): 2 -> commute. Z-logical = row r=0 (d qubits).
    L = np.zeros((1, n), dtype=np.uint8)
    L[0, :d] = 1
    # verify: H @ L.T == 0
    assert not (H @ L.T % 2).any(), "logical Z fails to commute with X-checks"
    return H, L


def sample_shots(H: np.ndarray, L: np.ndarray, p: float, shots: int, seed: int):
    """IID X flips at rate p; returns (syndromes, errors, obs) arrays."""
    rng = np.random.default_rng(seed)
    n = H.shape[1]
    errs = (rng.random((shots, n)) < p).astype(np.uint8)
    syn = (errs @ H.T) % 2
    obs = (errs @ L.T) % 2
    return syn.astype(np.uint8), errs, obs.astype(np.uint8)


def exact_coset_logZ(H: np.ndarray, L: np.ndarray, lam: float,
                     syndrome: np.ndarray, lam_mech: np.ndarray | None = None):
    """Exact class partition functions for the code-capacity problem by
    direct coset enumeration. lam = ln((1-p)/p) uniform per data qubit.

    Generator basis S = rowspace(H) (reduced to independent rows):
    candidate errors per class: e_lam XOR u^T S. 2^rank terms per class;
    feasible to rank <= 22.

    Returns (n_classes, 2^rank-class logZ array, class reps) with classes =
    1 + k logical shifts (k = 1 for the surface code: classes 0/1).
    """
    from .gf2 import solve_gf2
    n = H.shape[1]
    # independent generator rows
    Hr = _independent_rows(H)
    G = Hr.shape[0]
    # logical shifts: solve [H; L] l = (0, 1)
    Lb = L.reshape(1, n).astype(np.uint8)
    shift = solve_gf2(np.vstack([H, Lb]),
                      np.concatenate([np.zeros(H.shape[0], dtype=np.uint8),
                                      np.ones(1, dtype=np.uint8)]))
    assert shift is not None
    # syndrome solution e0: H e0 = s  (columns of H are the unknowns)
    e0 = solve_gf2(H.astype(np.uint8), syndrome.astype(np.uint8))
    assert e0 is not None, "syndrome not satisfiable"
    reps = np.vstack([e0, e0 ^ shift])               # 2 classes
    Us = ((np.arange(1 << G)[:, None] >> np.arange(G)[None, :]) & 1).astype(np.uint8)
    coset_par = (Us @ Hr) % 2                        # (2^G, n)
    logZs = np.empty((2, ), dtype=np.float64)
    for ci in range(2):
        x = reps[ci][None, :] ^ coset_par
        ex = x.sum(axis=1) * lam
        logZs[ci] = np.logaddexp.reduce(-ex)         # ln sum exp(-E)
    return logZs, reps, Hr


def _independent_rows(H: np.ndarray) -> np.ndarray:
    from .gf2 import gf2_rank
    rows = []
    cur = np.zeros((0, H.shape[1]), dtype=np.uint8)
    for i in range(H.shape[0]):
        cand = np.vstack([cur, H[i:i + 1]])
        if gf2_rank(cand) > gf2_rank(cur):
            cur = cand
            rows.append(i)
    return H[rows]


def _solve_syndrome(H: np.ndarray, s: np.ndarray) -> np.ndarray:
    """Particular syndrome solution e0 with H e0 = s (columns unknowns)."""
    from gf2 import solve_gf2
    e0 = solve_gf2(H.astype(np.uint8), s.astype(np.uint8))
    if e0 is None:
        raise RuntimeError("syndrome not satisfiable: sampling mismatch")
    return e0
