"""Gate-D model: coset-MRF construction per arXiv:2608.25545 Sec 7.7.

Input: a stim DEM factored by qldpc_dec.dem_matrices.dem_to_matrices(merge=False)
into H (detectors x mechanisms), A (observables x mechanisms), lam (mechanism
weights ln((1-p)/p) = -ln w, exact Eq-8 match).

Components (paper Sec 7.7 devices i-iii):
  1. Trivial-kernel generator basis S (kernel of [H; L] over GF(2)), then
     greedy pairwise sparsification to mean generator weight ~4.
  2. Candidate classes 1 + k single-logical shifts (spacetime convention,
     Sec 8.4): particular solutions l_lambda of H l = 0, L l = lambda.
  3. Lightening (device ii): reduce each representative to a local weighted
     minimum within its coset by XORing trivial-kernel generators.

AIS state semantics: u in F2^G over the SPARSIFIED basis; the candidate error
for class lambda is x = r_lambda + u^T K (mod 2); energy E(x) = sum_v x_v lam_v.
"""
from __future__ import annotations

import numpy as np

from .gf2 import gf2_rank, kernel_basis, solve_gf2

POPCOUNT_LUT = np.array([bin(i).count("1") for i in range(256)], dtype=np.int16)


class GateDModel:
    """Coset-MRF data for one DEM (one-time build, reused across syndromes)."""

    def __init__(self, H: np.ndarray, L: np.ndarray, lam: np.ndarray,
                 sparse_basis: np.ndarray | None = None,
                 log_shifts: np.ndarray | None = None):
        self.H = np.ascontiguousarray(H, dtype=np.uint8)
        self.L = np.ascontiguousarray(L, dtype=np.uint8)
        self.lam = np.ascontiguousarray(lam, dtype=np.float64)
        self.num_dets, self.n = self.H.shape
        self.num_log = self.L.shape[0]
        import time as _t
        t0 = _t.time()
        if sparse_basis is not None:
            self.K = np.ascontiguousarray(sparse_basis, dtype=np.uint8)
            self.G = self.K.shape[0]
        else:
            K_raw = kernel_basis(np.vstack([self.H, self.L]))
            self.G = K_raw.shape[0]
            self.K, _ = sparsify(K_raw, self.n)
        self.K_weights = POPCOUNT_LUT[
            np.packbits(self.K, axis=1, bitorder="little")
        ].sum(axis=1).astype(np.int64)
        # per-mechanism generator lists (flat CSR over mechanisms)
        nz_m, nz_g = np.nonzero(self.K.T)          # (edges,): mech, gen
        order = np.argsort(nz_m, kind="stable")
        self.col_ptr = np.concatenate([
            np.zeros(1, dtype=np.int64),
            np.cumsum(np.bincount(nz_m, minlength=self.n))
        ]).astype(np.int64)
        self.col_gens = nz_g[order].astype(np.int64)   # generator idx grouped by mech
        # per-generator mechanism lists (flat CSR over generators)
        nz_g2, nz_c2 = np.nonzero(self.K)
        o2 = np.argsort(nz_g2, kind="stable")
        self.row_ptr = np.concatenate([
            np.zeros(1, dtype=np.int64),
            np.cumsum(np.bincount(nz_g2, minlength=self.G))
        ]).astype(np.int64)
        self.row_mechs = nz_c2[o2].astype(np.int64)
        self.basis_seconds = _t.time() - t0
        if log_shifts is not None:
            self.log_shifts = np.ascontiguousarray(log_shifts, dtype=np.uint8)
        else:
            self.log_shifts = self._solve_log_shifts()

    def _solve_log_shifts(self) -> np.ndarray:
        """Particular solutions l_lambda of H l = 0, L l = e_i, i = 1..k."""
        aug = np.vstack([self.H, self.L]).astype(np.uint8)
        I = np.eye(self.num_log, dtype=np.uint8)
        zeros = np.zeros(self.num_dets, dtype=np.uint8)
        rows = []
        for i in range(self.num_log):
            b = np.concatenate([zeros, I[i]]).astype(np.uint8)
            x = solve_gf2(aug, b)
            if x is None:
                raise RuntimeError(f"log shift {i}: inconsistent system")
            rows.append(x.astype(np.uint8))
        return np.array(rows, dtype=np.uint8)

    # ---------- per-syndrome ----------
    def representatives(self, e0: np.ndarray, lighten: bool = True,
                        max_sweeps: int = 8) -> np.ndarray:
        """(1+k, n) class representatives around the seed correction e0:
        row 0 = e0 (lambda=0 class), row 1+i = e0 ^ l_i.  When lighten,
        each is greedily reduced to a local weighted minimum of its coset."""
        reps = np.vstack([
            e0.astype(np.uint8),
            e0.astype(np.uint8)[None, :] ^ self.log_shifts,
        ])
        if not lighten:
            return reps
        w = self.lam
        for ri in range(reps.shape[0]):
            r = reps[ri].copy()
            e_r = float(np.dot(r, w))
            for _ in range(max_sweeps):
                improved = False
                for c in np.argsort(self.K_weights):
                    mechs = self.row_mechs[self.row_ptr[c]:self.row_ptr[c + 1]]
                    if mechs.size == 0:
                        continue
                    # energy delta of toggling r over 'mechs'
                    d = float(w[mechs].sum()) - 2.0 * float(np.dot(r[mechs], w[mechs]))
                    if d < -1e-9:
                        r[mechs] ^= 1
                        e_r += d
                        improved = True
                if not improved:
                    break
            reps[ri] = r
        return reps

    def residuals(self, syndrome: np.ndarray, reps: np.ndarray) -> np.ndarray:
        """delta_v = H[:,v] r + s_v (mod 2) for each representative row:
        the defect parity the generator perturbations must cancel."""
        s = syndrome.astype(np.uint8)
        Hr = (self.H @ reps.T.astype(np.uint64)) % 2           # (dets, 1+k)
        return (Hr.T ^ s[None, :]).astype(np.uint8)            # (1+k, dets)

    def syndrome_of(self, e: np.ndarray) -> np.ndarray:
        return (self.H @ e.astype(np.uint64)) % 2

    def logical_of(self, e: np.ndarray) -> np.ndarray:
        return (self.L @ e.astype(np.uint64)) % 2


# ---------------- sparsifier ----------------
def sparsify(K: np.ndarray, n_bits: int, max_passes: int = 15, target: float = 4.0,
             verbose: bool = False):
    """Greedy pairwise reduction (paper Sec 7.7 device i): repeatedly
    K_i <- K_i xor K_j whenever it strictly lowers |K_i|. Bit-packed rows,
    occupancy-index candidate generation. Span is preserved (unimodular)."""
    import time as _t
    t0 = _t.time()
    m = K.shape[0]
    Kp = np.packbits(K.astype(np.uint8), axis=1, bitorder="little")
    w = POPCOUNT_LUT[Kp].sum(axis=1).astype(np.int64)
    rows_of_bit: dict[int, list[int]] = {}
    for rr in range(m):
        for b in np.nonzero(K[rr])[0]:
            rows_of_bit.setdefault(int(b), []).append(rr)
    bit_rows = {k: np.asarray(v, dtype=np.int64) for k, v in rows_of_bit.items()}
    for ps in range(max_passes):
        improved = 0
        order = np.argsort(-w)
        for ri in order:
            i = int(ri)
            if w[i] <= max(3, int(target) - 1):
                continue
            Ki_bits = np.unpackbits(Kp[i], bitorder="little")[:n_bits]
            nzb = np.nonzero(Ki_bits)[0]
            lists = [bit_rows[int(b)] for b in nzb if int(b) in bit_rows]
            if not lists:
                continue
            cand = np.unique(np.concatenate(lists))
            cand = cand[cand != i]
            if cand.size == 0:
                continue
            ands = POPCOUNT_LUT[(Kp[cand] & Kp[i])].sum(axis=1).astype(np.int64)
            new_wt = w[i] + w[cand] - 2 * ands
            jbest = int(np.argmin(new_wt))
            if new_wt[jbest] < w[i]:
                jj = int(cand[jbest])
                Kp[i] ^= Kp[jj]
                w[i] = int(new_wt[jbest])
                improved += 1
        if verbose:
            print(f"  sparsify pass {ps+1}: improved {improved}, mean wt {w.mean():.3f}, "
                  f"t={_t.time()-t0:.1f}s")
        if improved == 0:
            break
        if w.mean() <= target:
            break
    Kout = np.unpackbits(Kp, axis=1, bitorder="little")[:, :n_bits].astype(np.uint8)
    return Kout, w
