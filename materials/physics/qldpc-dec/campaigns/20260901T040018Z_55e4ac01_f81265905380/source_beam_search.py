"""Beam-search decoder, reimplementation from arXiv:2512.07057 (gate B).

Implements Algorithm 3 (beam search decoder) using the min-sum BP of
Appendix A (eqs. A1-A4), branching on the least reliable error node
(smallest |sum of posterior LLRs|), scoring paths by summed reliability
normalized by iterations, pruning to beam width, warm-starting from the
last round's error-to-detector messages. Pure numpy -- no C++ extension,
so runtimes are NOT comparable to the paper's Table II/III; only the
logical-error-rate ratios (gate B primary) are reproduced. Timing claims
from the paper are tagged hardware-conditional in pre_statement.md.

(H, A, lam): raw stim DEM mechanisms (dem_matrices.dem_to_matrices,
merge=False), matching the paper's Appendix A definition. H has 8784
columns for the gross-code Z-memory DEM -- identical to the count of DEM
error mechanisms.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import stim

from .dem_matrices import dem_to_matrices


class _MaskedBP:
    """Masked min-sum BP (paper Algorithm 2) over active (unmasked) nodes."""

    def __init__(self, H: sp.csr_matrix, lam: np.ndarray, masked: set[int],
                 syndrome: np.ndarray):
        self.num_dets, self.num_errs = H.shape
        self.lam = lam
        Hc = H.tocsc()
        # effective syndrome: masked-1 nodes have their incident detector
        # bits flipped (paper: "recast as all-masked-zero problem")
        s = syndrome.astype(np.int64).copy()
        for j, v in masked.items():
            if v:
                dets_j = Hc.indices[Hc.indptr[j] : Hc.indptr[j + 1]]
                s[dets_j] ^= 1
        self.s = s
        self.masked = masked
        active = np.ones(self.num_errs, dtype=bool)
        if masked:
            active[list(masked)] = False
        self.act = np.where(active)[0]
        self.act_set = set(self.act.tolist())
        # detectors per active error node
        self.err_dets = {}
        for j in self.act:
            self.err_dets[j] = Hc.indices[Hc.indptr[j] : Hc.indptr[j + 1]]
        # active neighbor sets per detector
        Hr = H.tocsr()
        self.det_errs = []
        for i in range(self.num_dets):
            self.det_errs.append(
                [j for j in Hr.indices[Hr.indptr[i] : Hr.indptr[i + 1]] if active[j]]
            )
        # active-submatrix product for convergence check
        self.H_active = H[:, self.act].tocsr()
        self.s_active = ()
        # NOTE: masked-1 flips make the reduced syndrome s the syndrome of
        # the all-masked-zero problem; convergence == parity over ACTIVE
        # columns equals s.
        self._s_active = s  # full-length; slice when comparing

    def run(self, edge_msgs: np.ndarray, max_iters: int):
        """Returns (converged, e_hat, sum_llr, iters_run, nu_per_active_err).

        e_hat/sum_llr are full-length vectors with masked entries = 0.
        edge_msgs: initial error->detector LLR per active node (E_j->i(0)).
        """
        num_pairs = sum(len(v) for v in self.det_errs)
        det_of_pair = np.empty(num_pairs, dtype=np.int64)
        err_of_pair = np.empty(num_pairs, dtype=np.int64)
        k = 0
        for i, js in enumerate(self.det_errs):
            for j in js:
                det_of_pair[k] = i
                err_of_pair[k] = j
                k += 1
        # per-active-error node, pair indices and position within detector rows
        nu = np.empty(num_pairs)
        # initialize: nu_j->i(0) = edge_msgs[j] (paper: E_j->i(0) = edge msgs)
        for idx in range(num_pairs):
            nu[idx] = edge_msgs[err_of_pair[idx]]
        err_pair_idx: dict[int, list[int]] = {}
        for idx in range(num_pairs):
            err_pair_idx.setdefault(int(err_of_pair[idx]), []).append(idx)
        err_pair_idx = {j: np.asarray(v) for j, v in err_pair_idx.items()}

        sum_llr = np.zeros(self.num_errs)
        e_hat = np.zeros(self.num_errs, dtype=np.uint8)
        converged = False
        iters_run = 0
        for t in range(1, max_iters + 1):
            iters_run = t
            # ----- check -> error (A1) -----
            absnu = np.abs(nu)
            # per-detector: two smallest |nu| and total sign
            min1 = np.full(self.num_dets, np.inf)
            min2 = np.full(self.num_dets, np.inf)
            argmin = np.full(self.num_dets, -1, dtype=np.int64)
            sgn_prod = np.ones(self.num_dets)
            deg = np.zeros(self.num_dets, dtype=np.int64)
            for i, js in enumerate(self.det_errs):
                if not js:
                    continue
                ks = err_pair_idx_of(i, det_of_pair)
                a = absnu[ks]
                order = np.argsort(a, kind="stable")
                min1[i] = a[order[0]]
                argmin[i] = ks[order[0]]
                if len(ks) > 1:
                    min2[i] = a[order[1]]
                sg = np.sign(nu[ks])
                if np.any(sg == 0):
                    sg = np.where(sg == 0, 1.0, sg)
                sgn_prod[i] = float(np.prod(sg))
                deg[i] = len(ks)
            mu = np.empty(num_pairs)
            for idx in range(num_pairs):
                i = det_of_pair[idx]
                if deg[i] == 0:
                    mu[idx] = 0.0
                    continue
                m = min2[i] if idx == argmin[i] else min1[i]
                # kappa: product of signs of the OTHER neighbors = sgn_prod / own sign
                own = nu[idx]
                own_sign = 1.0 if own >= 0 else -1.0
                kappa = sgn_prod[i] / own_sign
                mu[idx] = kappa * (1 if self.s[i] == 0 else -1) * m
            # ----- error -> check (A2) -----
            for j in self.act:
                ks = err_pair_idx[j]
                tot = mu[ks].sum()
                nu[ks] = self.lam[j] + tot - mu[ks]
            # ----- posteriors (A3), hard decision (A4) -----
            for j in self.act:
                ks = err_pair_idx[j]
                post = self.lam[j] + mu[ks].sum()
                sum_llr[j] += post
                e_hat[j] = 0 if post > 0 else 1
            # ----- convergence on the reduced problem -----
            reduced = (self.H_active @ e_hat[self.act].astype(np.int64)) % 2
            if np.array_equal(reduced.astype(np.int64), self.s):
                converged = True
                break
        return converged, e_hat, sum_llr, iters_run, nu, det_of_pair, err_of_pair, err_pair_idx


def err_pair_idx_of(det, det_of_pair):
    return np.where(det_of_pair == det)[0]


class BeamSearchDecoder:
    """Beam search decoder (arXiv:2512.07057 Algorithm 3)."""

    def __init__(
        self,
        dem: stim.DetectorErrorModel,
        beam_width: int = 8,
        initial_iters: int = 30,
        iters_per_round: int = 20,
        max_rounds: int = 10,
        num_results: int = 1,
    ):
        H, A, lam = dem_to_matrices(dem, merge=False)
        self.H = H
        self.A = A
        self.lam = lam
        self.num_dets, self.num_errs = H.shape
        self.num_obs = A.shape[0]
        self.beam_width = beam_width
        self.initial_iters = initial_iters
        self.iters_per_round = iters_per_round
        self.max_rounds = max_rounds
        self.num_results = num_results
        self.weight = np.abs(lam)

    def decode(self, syndrome: np.ndarray) -> np.ndarray:
        syndrome = np.asarray(syndrome).astype(np.int64)
        results: list[tuple[float, np.ndarray]] = []

        # ---------- seed: standard BP (no masked nodes) ----------
        eng = _MaskedBP(self.H, self.lam, {}, syndrome)
        conv, e_hat, sum_llr, iters, nu, dop, eop, epi = eng.run(None or np.zeros(self.num_errs), self.initial_iters)
        if conv:
            return e_hat
        seed_msgs = np.zeros(self.num_errs)
        for j in eng.act:
            ks = epi[j]
            seed_msgs[j] = nu[ks].mean()
        nxt = int(np.argmin(np.abs(sum_llr)))
        paths = [{"msgs": seed_msgs, "masked": {}, "next": nxt, "score": 0.0}]

        for rnd in range(self.max_rounds):
            next_set = []
            for path in paths:
                for val in (0, 1):
                    masked = dict(path["masked"])
                    masked[path["next"]] = val
                    eng2 = _MaskedBP(self.H, self.lam, masked, syndrome)
                    conv, e_h, s_lr, it2, nu2, _, _, epi2 = eng2.run(
                        path["msgs"], self.iters_per_round
                    )
                    if conv:
                        e_full = e_h.copy()
                        for j, v in masked.items():
                            e_full[j] = v
                        if np.array_equal((self.H @ e_full.astype(np.int64)) % 2, syndrome):
                            wt = float(np.dot(e_full, self.weight))
                            results.append((wt, e_full))
                            if len(results) >= self.num_results:
                                return min(results, key=lambda r: r[0])[1]
                    unmasked = [j for j in range(self.num_errs) if j not in masked]
                    if not unmasked:
                        continue
                    score = float(np.sum(np.abs(s_lr[unmasked])) / max(it2, 1))
                    nextj = int(min(unmasked, key=lambda j: abs(s_lr[j])))
                    # warm-start messages for children
                    child_msgs = np.zeros(self.num_errs)
                    for j in eng2.act:
                        ks = epi2[j]
                        child_msgs[j] = nu2[ks].mean()
                    next_set.append(
                        {"msgs": child_msgs, "masked": masked, "next": nextj, "score": score}
                    )
            if not next_set:
                break
            next_set.sort(key=lambda p: p["score"], reverse=True)
            paths = next_set[: self.beam_width]
            if all(p["next"] < 0 for p in paths):
                break
        if results:
            return min(results, key=lambda r: r[0])[1]
        # no valid solution: best-effort (top path's partial decode)
        e_full = np.zeros(self.num_errs, dtype=np.uint8)
        if paths:
            for j, v in paths[0]["masked"].items():
                e_full[j] = v
        return e_full


class BeamSearchBatchDecoder:
    """Batch adapter: predict observable flips from decoded corrections."""

    def __init__(self, dem: stim.DetectorErrorModel, **kwargs):
        self.dec = BeamSearchDecoder(dem, **kwargs)
        self.A = self.dec.A.tocsr()

    def decode_batch(self, shots: np.ndarray) -> np.ndarray:
        n = len(shots)
        out = np.empty((n, self.A.shape[0]), dtype=np.int64)
        for i in range(n):
            e_hat = self.dec.decode(shots[i])
            out[i] = (self.A @ e_hat.astype(np.int64)) % 2
        return out


if __name__ == "__main__":
    circuit = stim.Circuit("""
        X_ERROR(0.15) 0 1 2
        M 0 1 2
        DETECTOR rec[-2] rec[-1]
        DETECTOR rec[-3] rec[-2]
        OBSERVABLE_INCLUDE(0) rec[-1] rec[-2]
    """)
    dem = circuit.detector_error_model()
    dec = BeamSearchDecoder(dem, beam_width=4, initial_iters=5, iters_per_round=5, max_rounds=4)
    for syn in ([1, 0], [1, 1], [0, 1]):
        print(syn, "->", dec.decode(np.array(syn)))
