"""CRN-AIS engine over the coset MRF (arXiv:2608.25545 Sec 7.3 + 7.7).

Path: the concentrated base is part of the bridge (spacetime device iii):
  gamma_beta(u) prop Bern(q0)(u)^{1-beta} * exp(-beta E(u)),
where E(u) = sum_v lam_v x_v with x = r_lambda XOR parity(u over
generators touching v), and the base is the iid Bernoulli(q0) product over
the G generator coordinates (exactly normalized at beta=0: Z_0 = 1).
beta runs 0 -> 1 on a power-concave schedule; at beta=1 the target is the
true Boltzmann exp(-E) of the class MRF.

Chain kernel: one single-coordinate Metropolis sweep per stage, targeting
gamma_beta. Flipping generator c from 0->1 (s=+1) or 1->0 (s=-1):
  dlog gamma = (1-beta) * s * ln(q0/(1-q0)) - beta * dE,
  dE = sum_{v in supp(c)} lam_v (1 - 2 x_v),
accepted with prob min(1, exp(dlog gamma)). Rows = (class, chain); each
class runs its OWN chain state while the per-step uniforms are shared
across classes (CRN) and chain inits are identical (device: common random
numbers), so the class comparison is a paired one on one stream.

Per-chain AIS weight (log domain; Neal 2001): at each bridge step, evaluate
the unnormalized-density ratio at the PRE-transition state, then apply a
Metropolis kernel invariant for the new beta:
  logW += -(beta_t-beta_{t-1}) * [E(u_{t-1}) + logpdf_base(u_{t-1})].
The weight is a path product; intermediate Markov transitions do not telescope
into an endpoint-only target/base ratio.
"""
from __future__ import annotations

import numpy as np


def beta_schedule(T: int, kind: str = "power2") -> np.ndarray:
    t = np.arange(T + 1, dtype=np.float64) / T
    return t * t if kind == "power2" else t


def mechanism_parity(u: np.ndarray, row_mechs: np.ndarray, row_ptr: np.ndarray,
                     n: int) -> np.ndarray:
    """u^T K mod 2: (Kch, G) uint8 -> (Kch, n) uint8."""
    Kch = u.shape[0]
    out = np.zeros((Kch, n), dtype=np.uint8)
    for c in range(u.shape[1]):
        mechs = row_mechs[row_ptr[c]:row_ptr[c + 1]]
        if mechs.size == 0:
            continue
        sel = np.nonzero(u[:, c])[0]
        if sel.size:
            out[sel[:, None], mechs[None, :]] ^= 1
    return out


def metropolis_sweep(
    u: np.ndarray,
    x: np.ndarray,
    energy: np.ndarray,
    ones: np.ndarray,
    beta: float,
    lam: np.ndarray,
    row_mechs: np.ndarray,
    row_ptr: np.ndarray,
    uniforms: np.ndarray,
    lk: float,
) -> None:
    """One Metropolis scan targeting gamma_beta, updating cached state values.

    ``u`` is (C, Kch, G), ``x`` is (C, Kch, n), and ``uniforms`` is
    (Kch, G) shared across classes. ``energy`` and ``ones`` cache E(x) and
    the Hamming weight of u so AIS path weights need no full-state rescans.
    """
    logU = np.log(uniforms)
    for c in range(u.shape[2]):
        mechs = row_mechs[row_ptr[c]:row_ptr[c + 1]]
        if mechs.size == 0:
            continue
        ucol = u[:, :, c]
        s = np.where(ucol > 0, -1.0, 1.0)
        xs = x[:, :, mechs]
        dE = (1 - 2.0 * xs) @ lam[mechs]
        dlog = (1.0 - beta) * s * lk - beta * dE
        accept = logU[:, c][None, :] < dlog
        if not accept.any():
            continue
        rows_c, rows_k = np.nonzero(accept)
        energy[rows_c, rows_k] += dE[rows_c, rows_k]
        ones[rows_c, rows_k] += s[rows_c, rows_k].astype(np.int64)
        x[rows_c[:, None], rows_k[:, None], mechs[None, :]] ^= 1
        u[rows_c, rows_k, c] ^= 1


def base_logpdf(u_flat: np.ndarray, G: int, q0: float) -> np.ndarray:
    """log Bernoulli(q0) product pdf over G coordinates: (..., G) -> (...)."""
    a = u_flat.sum(axis=-1).astype(np.float64)
    return a * np.log(q0) + (G - a) * np.log1p(-q0)


class AISEngine:
    """All (class, chain) rows of one syndrome, vectorized."""

    def __init__(self, row_mechs: np.ndarray, row_ptr: np.ndarray,
                 lam: np.ndarray, T: int = 64, K: int = 64, q0: float = 0.02,
                 seed: int = 0, beta_kind: str = "power2"):
        self.row_mechs = row_mechs
        self.row_ptr = row_ptr
        self.G = len(row_ptr) - 1
        self.n = len(lam)
        self.lam = np.ascontiguousarray(lam, dtype=np.float64)
        self.T = T
        self.Kch = K
        self.q0 = q0
        self.betas = beta_schedule(T, beta_kind)
        self.seed = seed
        self.lk = float(np.log(q0 / (1.0 - q0)))

    def run(self, reps: np.ndarray, seed_offset: int = 0) -> dict:
        """Estimate each class partition function with standard forward AIS.

        ``reps`` has shape (C, n). Chains start exactly from the normalized
        Bernoulli(q0) base. For bridge t, the importance increment is evaluated
        at u_{t-1}; a Metropolis sweep invariant for beta_t then produces the
        next state. The final transition is omitted because no later weight
        consumes it.
        """
        C = reps.shape[0]
        Kch, G, T = self.Kch, self.G, self.T
        rng = np.random.default_rng(self.seed + seed_offset)
        lam = self.lam

        u_init = (rng.random((Kch, G)) < self.q0).astype(np.uint8)
        par = mechanism_parity(
            u_init, self.row_mechs, self.row_ptr, self.n
        )
        x = reps[:, None, :] ^ par[None, :, :]
        u = np.repeat(u_init[None, :, :], C, axis=0)
        initial_ones = u_init.sum(axis=1, dtype=np.int64)
        ones = np.repeat(initial_ones[None, :], C, axis=0)
        energy = (x * lam).sum(axis=2)
        logW = np.zeros((C, Kch), dtype=np.float64)
        log_q0 = float(np.log(self.q0))
        log_1mq0 = float(np.log1p(-self.q0))

        for t in range(1, T + 1):
            delta_beta = self.betas[t] - self.betas[t - 1]
            log_base = ones * log_q0 + (G - ones) * log_1mq0
            logW -= delta_beta * (energy + log_base)
            if t < T:
                uniforms = rng.random((Kch, G))
                metropolis_sweep(
                    u,
                    x,
                    energy,
                    ones,
                    self.betas[t],
                    lam,
                    self.row_mechs,
                    self.row_ptr,
                    uniforms,
                    self.lk,
                )

        m = logW.max(axis=1)
        logZ = (
            m
            + np.log(np.exp(logW - m[:, None]).sum(axis=1))
            - np.log(Kch)
        )
        return {"logZ": logZ, "logW": logW}
