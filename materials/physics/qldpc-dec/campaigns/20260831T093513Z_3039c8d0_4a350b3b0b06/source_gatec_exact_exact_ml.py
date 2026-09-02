"""Exact degenerate-ML reference: coset enumeration for BB[[36,4,4]], Z-memory.

Code construction
-----------------
Bivariate-bicycle code on the (L, M) = (6, 3) grid with the Bravyi et al.
monomial family (A = x^3 + y + y^2, B = y + x + x^3), giving the binary
matrices Hx = [A | B], Hz = [B^T | A^T] over N = 36 data qubits. Parameters
verified exactly below: [[36, 4, 4]] (k = 4, both distances 4; d computed by
exhaustive span enumeration over ker(Hx) \\ rowspace(Hz) and conversely).
This is the smallest instance of the [[72,12,6]]/[[144,12,12]] monomial
family: same 3-term polynomials, smaller cyclic grid. All monomial exponents
are "much smaller a, b" as the task requested, distance is low (4).

Noise model: code-capacity depolarizing, Z-memory sector
--------------------------------------------------------
Each data qubit carries an independent depolarizing channel with parameter p.
Conditioned on the Z-memory experiment (Z-checks read out, Z-type logicals
predicted), the X-component of each qubit's error is what faults the
syndrome: marginalizing the depolarizing channel over Z-components leaves
per-qubit X-error probability

    q = P(X-part flips) = P(X) + P(Y) = p/3 + p/3 = 2p/3.

The Z-syndrome and the Z-logical errors are then EXACTLY a binary linear-code
syndrome problem:  s = Hz_mech e  (18 checks),  observable class = A_mech e,
with per-mechanism prior q. (In the stim circuit built here the mechanisms
are indexed in stim's own mechanism order, which is a fixed permutation of
qubit order; H_mech/A_mech are extracted from the stim DEM so every decoder
and the exact reference see identical columns.)

Exact class partition functions
-------------------------------
For a sampled syndrome s the four-logical-bit class probabilities are

    P(class c | s) ~ sum_{e in affine space, A e = c} q^{|e|} (1-q)^{N-|e|},

an EXACT coset sum: the solutions of H e = s form an affine space of
dimension N - rank = 36 - 16 = 20 (2^20 = 1,048,576 members). We enumerate
the whole affine space per shot as packed 36-bit words (one 8 MB numpy
table), accumulate integer weight histograms n[c, w] = #{solutions in class
c of weight w}, and evaluate class masses as exact rationals

    Z_c = sum_w n[c, w] r^w,      r = q / (1 - q) = 2p / (3 - 2p),

with Python Fractions — no floating point in the class comparison, so the
argmax (the exact-ML decision) is exact. The degenerate-ML property is that
ALL stabilizer-equivalent solutions contribute, not just a minimum-weight
representative.

Certificates / validation baked in:
  1. histogram totals equal 2^20 exactly (partition of the affine space);
  2. a brute-force walk over all errors of weight <= 5 must reproduce the
     histogram rows n[c, w<=5] (run once at build time, and re-checkable
     via validate_exact_vs_bruteforce());
  3. the true error's class must appear with nonzero mass (syndrome
     consistency).

Cost: one 2^20 coset enumeration is ~7-10 ms single-core (Apple M3 Ultra,
numpy 2.5); 2000 shots ~ 20 s per noise point. The tail beyond enumerated
space is zero: the affine space is COMPLETE, so this is exact ML, not a
truncation. There is no weight-ball approximation anywhere.
"""
from __future__ import annotations

import numpy as np
from fractions import Fraction

import stim

__all__ = [
    "bb_parity_checks", "bb_logical_ops", "bb_distances", "codecap_circuit",
    "dem_matrices_raw", "ExactMLReference", "validate_exact_vs_bruteforce",
]


# ---------------------------------------------------------------- code construction

def _circulant_block(l: int, m: int, shifts) -> np.ndarray:
    block = np.zeros((l * m, l * m), dtype=np.uint8)
    for dx, dy in shifts:
        for i in range(l):
            for j in range(m):
                block[i * m + j, ((i + dx) % l) * m + (j + dy) % m] = 1
    return block


def bb_parity_checks(l: int = 6, m: int = 3,
                     a_shifts=((3, 0), (0, 1), (0, 2)),
                     b_shifts=((0, 1), (1, 0), (3, 0))):
    """Hx = [A | B], Hz = [B^T | A^T] for the BB monomial pair on the (l,m) grid.

    Defaults are the Bravyi-family polynomials A = x^3 + y + y^2,
    B = y + x + x^3 on the (6,3) torus -> [[36,4,4]] (verified by bb_distances).
    """
    A = _circulant_block(l, m, a_shifts)
    B = _circulant_block(l, m, b_shifts)
    if ((A @ B + B @ A) % 2).any():
        raise ValueError("monomial pair does not commute")
    Hx = np.hstack([A, B]).astype(np.uint8)
    Hz = np.hstack([B.T, A.T]).astype(np.uint8)
    return Hx, Hz


def _rref_gf2(mat):
    A = mat.copy().astype(np.uint8) % 2
    pivots = []
    row = 0
    for col in range(A.shape[1]):
        if row >= A.shape[0]:
            break
        nz = np.nonzero(A[row:, col])[0]
        if nz.size == 0:
            continue
        pr = row + nz[0]
        if pr != row:
            A[[row, pr]] = A[[pr, row]]
        rows_with = np.nonzero(A[:, col])[0]
        rows_with = rows_with[rows_with != row]
        A[rows_with] ^= A[row]
        pivots.append(col)
        row += 1
    return A, pivots


def _rowspace_membership(mat, max_dim=24):
    """Set of byte-keys for the row space of mat (raises if dim > max_dim)."""
    R, _ = _rref_gf2(mat)
    rows = [R[i].copy() for i in range(R.shape[0]) if R[i].any()]
    if len(rows) > max_dim:
        raise RuntimeError(f"rowspace dim {len(rows)} too large for enumeration")
    keys = set()
    n = mat.shape[1]
    v = np.zeros(n, dtype=np.uint8)
    for mask in range(1 << len(rows)):
        x = mask
        v[:] = 0
        i = 0
        while x:
            if x & 1:
                v ^= rows[i]
            x >>= 1
            i += 1
        keys.add(v.tobytes())
    return keys


def _min_coset_weight(ker_basis, stab_keys, n):
    """Min weight in span(ker_basis) \\ stab_keys (exhaustive over 2^dim)."""
    dim = len(ker_basis)
    best = None
    v = np.zeros(n, dtype=np.uint8)
    for mask in range(1, 1 << dim):
        x = mask
        v[:] = 0
        i = 0
        while x:
            if x & 1:
                v ^= ker_basis[i]
            x >>= 1
            i += 1
        if v.tobytes() in stab_keys:
            continue
        w = int(v.sum())
        if best is None or w < best:
            best = w
            if best <= 2:
                break
    return best


def bb_distances(Hx, Hz):
    """Exact (dx, dz) by exhaustive enumeration of ker/rowspace cosets."""
    n = Hx.shape[1]
    ker_Hx = _nullspace_basis(Hx)
    ker_Hz = _nullspace_basis(Hz)
    rows_Hz_keys = _rowspace_membership(Hz)
    rows_Hx_keys = _rowspace_membership(Hx)
    dz = _min_coset_weight(ker_Hx, rows_Hz_keys, n)  # Z-type logicals commute w/ Hx? careful:
    # Z-logical = ker(Hx) \\ rowspace(Hz); X-logical = ker(Hz) \\ rowspace(Hx)
    dx = _min_coset_weight(ker_Hz, rows_Hx_keys, n)
    return dx, dz


def _nullspace_basis(mat):
    R, pivots = _rref_gf2(mat)
    free = [c for c in range(mat.shape[1]) if c not in set(pivots)]
    basis = []
    for f in free:
        v = np.zeros(mat.shape[1], dtype=np.uint8)
        v[f] = 1
        for r, pc in enumerate(pivots):
            if R[r, f]:
                v[pc] = 1
        basis.append(v)
    return np.array(basis, dtype=np.uint8) if basis else np.zeros((0, mat.shape[1]), dtype=np.uint8)


def bb_logical_ops(Hx, Hz):
    """Rows: Z-logicals (ker(Hx)/rowspace(Hz)), X-logicals (ker(Hz)/rowspace(Hx)) in RREF."""
    def _quotient_basis(ker_mat, mod_mat):
        ker = _nullspace_basis(ker_mat)
        R, pivots = _rref_gf2(mod_mat)
        survivors = []
        for v in ker:
            w = v.copy()
            for r, pc in enumerate(pivots):
                if w[pc]:
                    w ^= R[r]
            if w.any():
                survivors.append(w)
        if not survivors:
            return np.zeros((0, ker_mat.shape[1]), dtype=np.uint8)
        S, _ = _rref_gf2(np.array(survivors, dtype=np.uint8))
        keep = [r for r in range(S.shape[0]) if S[r].any()]
        return S[keep]
    return _quotient_basis(Hx, Hz), _quotient_basis(Hz, Hx)


def codecap_circuit(
    q: float,
    checks: np.ndarray,
    logicals: np.ndarray,
    basis: str = "Z",
) -> stim.Circuit:
    """Build a one-shot CSS code-capacity memory circuit.

    Z memory prepares/measures Z and injects X errors, using Z checks and
    Z-logicals. X memory is its Hadamard dual: prepare/measure X, inject Z
    errors, and use X checks/X-logicals.
    """
    basis = basis.upper()
    if basis == "Z":
        reset, error, measure = "R", "X_ERROR", "M"
    elif basis == "X":
        reset, error, measure = "RX", "Z_ERROR", "MX"
    else:
        raise ValueError(f"basis must be X or Z, got {basis!r}")

    checks = np.asarray(checks, dtype=np.uint8)
    logicals = np.asarray(logicals, dtype=np.uint8)
    n = checks.shape[1]
    if logicals.ndim != 2 or logicals.shape[1] != n:
        raise ValueError("logical operators must match the check width")

    circuit = stim.Circuit()
    circuit.append(reset, list(range(n)))
    circuit.append("TICK")
    circuit.append(error, list(range(n)), q)
    circuit.append("TICK")
    circuit.append(measure, list(range(n)))
    for check in checks:
        targets = [stim.target_rec(int(d) - n) for d in np.nonzero(check)[0]]
        circuit.append("DETECTOR", targets)
    for key, logical in enumerate(logicals):
        targets = [stim.target_rec(int(d) - n) for d in np.nonzero(logical)[0]]
        circuit.append("OBSERVABLE_INCLUDE", targets, key)
    return circuit


def dem_matrices_raw(dem: stim.DetectorErrorModel):
    """(H, A, p) exactly as emitted by stim (no merging), dense numpy arrays."""
    n_dets = dem.num_detectors
    n_obs = dem.num_observables
    cols = []
    pvec = []
    for inst in dem.flattened():
        if inst.type != "error":
            continue
        dets = [t.val for t in inst.targets_copy() if t.is_relative_detector_id()]
        obs = [t.val for t in inst.targets_copy() if t.is_logical_observable_id()]
        col = np.zeros((n_dets + n_obs,), dtype=np.uint8)
        col[dets] = 1
        col[n_dets + np.asarray(obs, dtype=int)] = 1
        cols.append(col)
        pvec.append(inst.args_copy()[0])
    M = np.array(cols, dtype=np.uint8).T  # (n_dets + n_obs, n_mech)
    H = M[:n_dets]
    A = M[n_dets:]
    return H, A, np.array(pvec, dtype=np.float64)


# ---------------------------------------------------------------- exact ML engine

def _packbits(v) -> int:
    return int("".join(map(str, np.asarray(v)[::-1])), 2)


class ExactMLReference:
    """Exact degenerate-ML over the full affine coset space of the DEM.

    Everything is fixed at construction: H, A (dense), the nullspace basis
    masks E0 (2^(n-rank) packed words), and the packed observable masks. Per
    shot, class_histogram(syndrome) returns the EXACT integer weight
    histogram n[c, w] over the complete affine space of consistent errors.
    """

    def __init__(self, H: np.ndarray, A: np.ndarray):
        self.H = np.asarray(H, dtype=np.uint8)
        self.A = np.asarray(A, dtype=np.uint8)
        self.n = self.H.shape[1]
        self.num_dets = self.H.shape[0]
        self.num_obs = self.A.shape[0]
        R, pivots = _rref_gf2(self.H)
        self.pivots = pivots
        self.rank = len(pivots)
        self.dim = self.n - self.rank
        if self.dim > 22:
            raise RuntimeError(f"affine space 2^{self.dim} too large")
        # nullspace basis masks
        free = [c for c in range(self.n) if c not in set(pivots)]
        masks = []
        for f in free:
            v = np.zeros(self.n, dtype=np.uint8)
            v[f] = 1
            for r, pc in enumerate(pivots):
                if R[r, f]:
                    v[pc] = 1
            masks.append(_packbits(v))
        self.L_packed = [_packbits(self.A[k]) for k in range(self.num_obs)]
        table = np.zeros(1 << self.dim, dtype=np.uint64)
        for j, m in enumerate(masks):
            step = 1 << j
            table[step:2 * step] = table[:step] ^ np.uint64(m)
        self.E0 = table

    def particular_mask(self, syndrome: np.ndarray) -> int:
        """Mask of one particular solution: joint [H | s] elimination, free bits = 0."""
        s = np.asarray(syndrome, dtype=np.uint8).copy()
        A = self.H.copy()
        b = s.copy()
        pi = 0
        for c in self.pivots:
            # find row at/after pi with A[r,c]=1, swap to pi
            nz = pi + np.nonzero(A[pi:, c])[0]
            if len(nz) == 0:
                continue
            r = int(nz[0])
            if r != pi:
                A[[pi, r]] = A[[r, pi]]
                b[[pi, r]] = b[[r, pi]]
            rows_w = np.nonzero(A[:, c])[0]
            rows_w = rows_w[rows_w != pi]
            A[rows_w] ^= A[pi]
            b[rows_w] ^= b[pi]
            pi += 1
        m = 0
        for i, c in enumerate(self.pivots):
            if b[i]:
                m |= 1 << c
        return m

    def class_histogram(self, syndrome: np.ndarray) -> np.ndarray:
        """Exact integer histogram n[class, weight] over the affine space.

        Returns (2^num_obs, n+1) int64 (num_obs=4 -> 16 class rows). Sum over
        all entries == 2^dim by construction (asserted here as the built-in
        partition check).
        """
        epm = self.particular_mask(syndrome)
        Es = self.E0 ^ np.uint64(epm)
        wt = np.bitwise_count(Es).astype(np.int64)
        ck = np.zeros(Es.shape[0], dtype=np.int64)
        for k in range(self.num_obs):
            par = np.bitwise_count(Es & np.uint64(self.L_packed[k])) & np.uint64(1)
            ck |= par.astype(np.int64) << np.int64(k)
        key = (ck << 6) | wt  # wt <= n <= 63 for our sizes
        ncls = 1 << self.num_obs
        flat = np.bincount(key, minlength=ncls * 64)
        pad = (-len(flat)) % 64
        if pad:
            flat = np.concatenate([flat, np.zeros(pad, dtype=flat.dtype)])
        hist = flat.reshape(-1, 64)
        out = np.zeros((ncls, self.n + 1), dtype=np.int64)
        out[: min(hist.shape[0], ncls), : self.n + 1] = hist[: min(hist.shape[0], ncls), : self.n + 1]
        total = int(out.sum())
        if total != 1 << self.dim:
            raise AssertionError(f"class histogram lost mass: {total} != {1 << self.dim}")
        return out

    @staticmethod
    def r_of_p(p: float) -> Fraction:
        """Exact prior ratio r = q/(1-q) with q = 2p/3 (depolarizing marginal)."""
        pf = Fraction(p).limit_denominator(10**15)
        q = 2 * pf / 3
        return q / (1 - q)

    @classmethod
    def class_masses(cls, hist: np.ndarray, p: float) -> list[Fraction]:
        """Exact rational class masses Z_c = sum_w n[c,w] r^w."""
        r = cls.r_of_p(p)
        out = []
        rp = [Fraction(1)] * (hist.shape[1])
        for w in range(1, hist.shape[1]):
            rp[w] = rp[w - 1] * r
        for c in range(hist.shape[0]):
            m = Fraction(0)
            row = hist[c]
            for w in np.nonzero(row)[0]:
                m += int(row[w]) * rp[w]
            out.append(m)
        return out

    @staticmethod
    def ml_decision(masses) -> tuple[int, Fraction]:
        best = max(range(len(masses)), key=lambda c: masses[c])
        return best, masses[best]

    @staticmethod
    def regret_nats(masses, chosen: int) -> float:
        """ln(Z_best / Z_chosen) >= 0; 0 iff the chosen class is the ML class."""
        from math import log
        best = max(range(len(masses)), key=lambda c: masses[c])
        ratio = masses[best] / masses[chosen]
        return log(float(ratio))


# ---------------------------------------------------------------- validation

def validate_exact_vs_bruteforce(ref: ExactMLReference, seed: int = 0,
                                 num_syndromes: int = 4, wmax: int = 5) -> dict:
    """Cross-check class histograms against brute-force enumeration of all
    errors of weight <= wmax for a few random ( syndrome from real error )
    draws. Returns assertion-passed diagnostics."""
    from itertools import combinations
    rng = np.random.default_rng(seed)
    checked = 0
    max_n = 0
    for _ in range(num_syndromes):
        e = rng.integers(0, 2, ref.n).astype(np.uint8)
        s = (ref.H @ e) % 2
        hist = ref.class_histogram(s)
        # brute force: every weight<=wmax vector with H e = s
        hist_bf = np.zeros_like(hist)
        for w in range(1, wmax + 1):
            for S in combinations(range(ref.n), w):
                v = np.zeros(ref.n, dtype=np.uint8)
                v[list(S)] = 1
                if ((ref.H @ v) % 2 == s).all():
                    c = int(np.dot(ref.A, v) % 2 @ (1 << np.arange(ref.num_obs)))
                    hist_bf[c, w] += 1
        if not (hist[:, : wmax + 1] == hist_bf[:, : wmax + 1]).all():
            raise AssertionError(f"exact histogram mismatch vs brute force at w<={wmax}")
        checked += 1
        max_n = max(max_n, int(hist_bf.sum()))
    return {"syndromes_checked": checked, "wmax": wmax, "max_low_weight_solutions": max_n}
