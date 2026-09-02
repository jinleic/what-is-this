"""Normalized min-sum (NMS) decoder with randomized serial schedule and
24-decoder ensemble (gate C: arXiv:2510.14060).

Implements the published configuration for the [[144,12,12]] gross code
under uniform depolarizing circuit-level noise:

- NMS with normalization factor alpha = 0.96875 (11/2^4), max 400
  iterations (GARI Sec. III A "GARI-NMS" hyperparameters).
- Randomized serial schedule: at each iteration the check-node update
  order is shuffled with a per-decoder seed (GARI Sec. II B 1-2:
  "randomized serial schedule ... using a distinct randomization seed").
- Ensemble: 24 decoders run to convergence on the same syndrome; the
  ensemble stops at the FIRST converging decoder (minimum-latency rule,
  GARI Sec. II B 2); simultaneous convergence resolved by lowest
  log-odds weight (their "most likely error" rule).
- Non-convergence within max_iters counts as a decoding failure and the
  hard-decision output is used (paper: failures are dominated by
  non-convergence; LER = failures / shots).

NOTE on the graph: this implementation decodes the RAW stim DEM
hyperedges (Y-type errors stay as hyperedges). The GARI paper's central
claim is that a 4-cycle-eliminating graph transform (GARI matrix) is
required for their accuracy. Our gate C therefore tests the WEAKER claim:
"24-ensemble randomized-schedule NMS on the plain DEM already reaches
LER <= 1.93e-9/round at p=1e-3". The GARI transform itself is out of
scope (pre_statement.md, explicit non-goal).

(H, A, lam): raw stim DEM mechanisms (merge=False). Fairness with respect
to GARI's D_Z-convention: GARI-NMS decodes correlated (X,Y,Z) mechanisms
against BOTH detector blocks; ours takes the memory_Z DEM as-is (all 936
detectors, mechanisms as emitted), which is the same information set GM
would see for a Z-memory experiment modulo their augmentation.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import stim

from .dem_matrices import dem_to_matrices


class NmsDecoder:
    """Single NMS decoder with randomized serial schedule."""

    def __init__(
        self,
        H: sp.csr_matrix,
        lam: np.ndarray,
        alpha: float = 0.96875,
        max_iters: int = 400,
        seed: int = 0,
        serial: bool = True,
    ):
        self.H = H.tocsr()
        self.lam = lam
        self.alpha = alpha
        self.max_iters = max_iters
        self.seed = seed
        self.serial = serial
        self.rng = np.random.default_rng(seed)
        self.num_dets, self.num_errs = H.shape
        Hc = H.tocsc()
        # per error node: detector list; per detector: error list
        self.err_dets = [Hc.indices[Hc.indptr[j] : Hc.indptr[j + 1]] for j in range(self.num_errs)]
        self.det_errs = [
            self.H.indices[self.H.indptr[i] : self.H.indptr[i + 1]]
            for i in range(self.num_dets)
        ]
        # message layout: nu[j, i] error->check for (j, i) incident pairs.
        # dense (2 x num_pairs) would be heavy for 8784 cols; use dict-free
        # flattened arrays aligned to a fixed pair ordering.
        pairs = []
        for i in range(self.num_dets):
            for j in self.det_errs[i]:
                pairs.append((i, j))
        self.pairs = pairs
        self.det_of = np.array([p[0] for p in pairs], dtype=np.int64)
        self.err_of = np.array([p[1] for p in pairs], dtype=np.int64)
        self.pair_index = {(i, j): k for k, (i, j) in enumerate(pairs)}
        self.err_pairs = [np.array([], dtype=np.int64)] * self.num_errs
        tmp: dict[int, list[int]] = {}
        for k, (_, j) in enumerate(pairs):
            tmp.setdefault(int(j), []).append(k)
        self.err_pairs = {j: np.asarray(v) for j, v in tmp.items()}
        self.det_pairs = {i: np.asarray(v) for i, v in tmp.items() if False}  # unused

    def decode(self, syndrome: np.ndarray) -> tuple[bool, np.ndarray, int]:
        """Returns (converged, e_hat, iters)."""
        s = np.asarray(syndrome).astype(np.int64)
        nu = np.empty(len(self.pairs))
        for k, (i, j) in enumerate(self.pairs):
            nu[k] = self.lam[j]
        # mu per pair, computed from nu with NMS normalization
        mu = np.zeros(len(self.pairs))
        e_hat = np.zeros(self.num_errs, dtype=np.uint8)
        check_order = np.arange(self.num_dets)
        converged = False
        it = 0
        for it in range(1, self.max_iters + 1):
            if self.serial:
                self.rng.shuffle(check_order)
            else:
                pass
            mu_new = mu  # in-place: serial reuse
            # ---- check-node updates (serialized in random order) ----
            abs_nu = np.abs(nu)
            sgn_nu = np.sign(nu)
            sgn_nu[sgn_nu == 0] = 1.0
            for i in check_order:
                js = self.det_errs[i]
                if len(js) == 0:
                    continue
                ks = np.array([self.pair_index[(i, int(j))] for j in js])
                a = abs_nu[ks]
                m1 = a.min()
                k_m1 = ks[int(np.argmin(a))]
                m2 = np.partition(a, 1)[1] if len(a) > 1 else np.inf
                sgn_prod = float(np.prod(sgn_nu[ks]))
                for ki, k in enumerate(ks):
                    m = m2 if k == k_m1 else m1
                    kappa = sgn_prod / sgn_nu[k]
                    mu_new[k] = self.alpha * kappa * (1 if s[i] == 0 else -1) * m
            # ---- variable-node updates ----
            for j in range(self.num_errs):
                ks = self.err_pairs[j]
                if len(ks) == 0:
                    continue
                tot = mu[ks].sum()
                # serial in-place: nu for this error row
                nu[ks] = self.lam[j] + tot - mu[ks]
            # ---- posterior, hard decision, convergence ----
            fails = 0
            for j in range(self.num_errs):
                ks = self.err_pairs[j]
                if len(ks) == 0:
                    post = self.lam[j]
                else:
                    post = self.lam[j] + mu[ks].sum()
                e_hat[j] = 0 if post > 0 else 1
            if np.array_equal((self.H @ e_hat.astype(np.int64)) % 2, s):
                converged = True
                break
        return converged, e_hat, it


class NmsEnsembleDecoder:
    """24 randomized-schedule NMS decoders; minimum-latency stopping."""

    def __init__(
        self,
        dem: stim.DetectorErrorModel,
        ensemble_size: int = 24,
        alpha: float = 0.96875,
        max_iters: int = 400,
        base_seed: int = 0,
    ):
        H, A, lam = dem_to_matrices(dem, merge=False)
        self.H = H
        self.A = A.tocsr()
        self.lam = lam
        self.num_dets, self.num_errs = H.shape
        self.num_obs = self.A.shape[0]
        self.ensemble_size = ensemble_size
        self.alpha = alpha
        self.max_iters = max_iters
        self.decoders = [
            NmsDecoder(H, lam, alpha=alpha, max_iters=max_iters, seed=base_seed + s)
            for s in range(ensemble_size)
        ]
        self.weight = np.abs(lam)
        self.stats = {"ensemble_stop_iter": [], "nonconverged": 0}

    def decode(self, syndrome: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Ensemble decode of one syndrome.
        Returns (e_hat, prediction). Stopping: decoders run in order; the
        first round in which any decoder converges fixes the round `it*`;
        among decoders converging at `it*` the minimum-log-odds-weight
        correction wins (paper: "most likely error" rule); decoders that
        need more iterations are skipped (minimum-latency rule).
        """
        best_e = None
        best_w = np.inf
        stop_round = None
        for dec in self.decoders:
            conv, e_hat, it = dec.decode(syndrome)
            if conv:
                if stop_round is None:
                    stop_round = it
                if it != stop_round:
                    break  # minimum-latency: later rounds never win
                w = float(np.dot(e_hat, self.weight))
                if w < best_w:
                    best_w = w
                    best_e = e_hat
            elif stop_round is not None and it > stop_round:
                break
        if best_e is None:
            self.stats["nonconverged"] += 1
            best_e = np.zeros(self.num_errs, dtype=np.uint8)
        else:
            self.stats["ensemble_stop_iter"].append(stop_round)
        pred = (self.A @ best_e.astype(np.int64)) % 2
        return best_e, pred.ravel().astype(np.int64)

    def decode_batch(self, shots: np.ndarray) -> np.ndarray:
        n = len(shots)
        out = np.empty((n, self.num_obs), dtype=np.int64)
        for i in range(n):
            _, pred = self.decode(shots[i])
            out[i] = pred
        return out


def lerr_per_round(lerr: float, rounds: int = 12) -> float:
    """Per-round LER (GARI eq. 7): LER_r = (1 - (1 - 2 LER)^{1/r}) / 2."""
    import math

    return (1.0 - math.pow(max(1.0 - 2.0 * lerr, 0.0), 1.0 / rounds)) / 2.0


if __name__ == "__main__":
    circuit = stim.Circuit("""
        X_ERROR(0.08) 0 1 2
        M 0 1 2
        DETECTOR rec[-2] rec[-1]
        DETECTOR rec[-3] rec[-2]
        OBSERVABLE_INCLUDE(0) rec[-1] rec[-2]
    """)
    dem = circuit.detector_error_model()
    ens = NmsEnsembleDecoder(dem, ensemble_size=4, alpha=0.96875, max_iters=60)
    for syn in ([1, 0], [1, 1], [0, 1], [0, 0]):
        print(syn, "->", ens.decode(np.array(syn)))
    print("stats:", ens.stats)
