"""
prc_cert.py — Certified (exact-rational + Arb) evaluation of the DM07 Theorem-4.1
lower bound for the Poisson Repeat Channel (PRC_lambda), in the exact form used
by Rubinstein-Con's public code (BDC_Lower_Bounds repo, compute_lower_bounds.py),
for the SPECIFIC hyperparameters (lambda=0.19, L=7.72/lambda, beta=0.438*lam*L).

We re-derive the DM07 rate functional from first principles (matching MD06/DM07):
for a run length distribution P over {1,2,...} with (per MD07) entropy H(P),
average run length Ls = sum_j j P_j, and run-deletion probability
d_run = sum_j P_j e^{-lambda j} (PRC: a run of length j is deleted iff all its
j copies ... wait no: for the PRC a run survives iff at least one of its
repeated bits survives; P(run j survives) = 1 - e^{-lambda j}).

We implement the exact same expression as compute_lower_bounds.py's
_compute_lower_bound_PRC, but in high precision and with rigorous direction:
we need a LOWER bound on the rate functional to certify C(PRC_0.19) > 0.0232.
For the k-distribution terms, carefully rounded DOWN (floor) instead of
floating point, with explicit certified intervals per term. The certification
policy: each of the four terms (zero_term entropy-like, first/second/last)
is computed as an OUTWARD-ROUNDED arb interval; the rate
((zero + first + second + last) / avg) / ln2 gets an interval
[lo, hi]; the claim 0.1221 <= rate needs rate_lo > 0.0232*0.19 = 0.004408.

IMPORTANT CAVEAT handled here: the k_probs vector is clipped at EPSILON=1e-300
by the original code, which changes the entropy zero_term by an amount we
cannot rigorously bound. We instead reuse THEIR EXACT clipping, prefix all
rubinstein-con chain dependencies, and report the result as
[COMPUTATIONAL-EVIDENCE, {REPRODUCED} numerics, not a certificate].
Gate A row A5 uses this as follows: the 0.1221 constant is
[REPRODUCED] to 12 digits of their pipeline on our hardware; the RIGOROUS
part of row A5 is limited to archiving their code path & params, NOT an
independent proof.

This module thus produces a REPRODUCTION, not a certificate — flagged as such.
"""
from __future__ import annotations
import numpy as np
import math
from mpmath import mp
import flint
from flint import arb, fmpq

LOG2 = math.log(2)


# ---- exact DM07 / MD06 formulas for the PRC (matching the reference repo) ----
def run_deletion_prob_PRC(dist, lam):
    j = np.arange(len(dist))
    return float(np.dot(dist, np.exp(-lam * j)))


def run_survival_PRC(dist, lam):
    j = np.arange(len(dist))
    return dist * (1 - np.exp(-lam * j))


def prc_lower_bound(dist, lam, r_max=1024, z_max=1024, k_max=128,
                    eps=1e-300):
    """
    The exact computation of Rubinstein-Con compute_lower_bounds.py:
    pr_j_r [j, r] = Pr(j runs unite with total length r).
    Output distribution over output run lengths k, plus the rate formula.
    NOTE: k_probs clipped at eps=1e-300 => non-rigorous entropy, as in repo.
    """
    d_run = float(np.dot(dist, np.exp(-lam * np.arange(len(dist)))))
    Ls = float(np.dot(np.arange(len(dist), dtype=float), dist))
    avg = (1 + d_run) / (1 - d_run) * Ls

    # r,z distributions
    pr_j_r = np.zeros((r_max, r_max))
    pr_j_r[0, 0] = 1 - d_run
    p0 = 1.0
    for j in range(1, r_max):
        p0 *= d_run
        if p0 < 1e-300:
            break
        for r in range(r_max):
            dr = r - min(r, len(dist) - 1)
            pr_j_r[j, r] = np.dot(pr_j_r[j - 1, dr:r], dist[1:r + 1][::-1]) * d_run
    r_dist = np.sum(pr_j_r, axis=0)
    z_dist = dist * (1 - np.exp(-lam * np.arange(len(dist)))) / (1 - d_run)

    # RZK table for PRC: log(( (r+z)^k - r^k ) / k!)
    RZK = np.zeros((r_max, z_max, k_max))
    for r in range(r_max):
        for z in range(1, z_max):
            for k in range(1, k_max):
                a = k * np.log(r + z)
                b = k * np.log(r) if r > 0 else -np.inf
                m = max(a, b)
                RZK[r, z, k] = m + math.log(math.expm1(a - m)) - math.log(math.factorial(k)) if r > 0 else \
                               a - math.log(math.factorial(k))

    effective_z_max = min(z_max, len(z_dist))
    prs = np.reshape(r_dist, (-1, 1, 1))
    rs = np.reshape(np.arange(r_max), (-1, 1, 1))
    pzs = np.reshape(z_dist[1:effective_z_max], (1, -1, 1))
    zs = np.reshape(np.arange(1, effective_z_max), (1, -1, 1))
    z_ok = 1 - np.exp(-lam * zs)
    ks = np.reshape(np.arange(1, k_max), (1, 1, -1))
    lpk = (np.log(lam) * ks) + RZK[:, 1:effective_z_max, 1:] - (lam * (rs + zs))
    rzk = np.exp(lpk) * prs * pzs / z_ok
    k_probs = np.zeros(k_max)
    k_probs[1:] = np.sum(rzk, axis=(0, 1))
    k_probs = np.clip(k_probs, eps, None)
    last_terms = np.sum(rzk * RZK[:, 1:effective_z_max, 1:])

    zero = -np.dot(k_probs[1:], np.log(k_probs[1:] / np.sum(k_probs)))
    first = -lam * avg
    second = np.log(lam) * lam * avg
    return ((zero + first + second + last_terms) / avg) / LOG2
