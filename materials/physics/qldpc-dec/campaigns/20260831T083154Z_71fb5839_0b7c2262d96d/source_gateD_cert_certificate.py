"""Paired-bootstrap certificate (arXiv:2608.25545 Sec 7.4).

Given per-chain log-importance-weights logW (C classes x K chains) from one
CRN-AIS run, the decision is lambda_hat = argmax_c log Z_hat_c with
log Z_hat_c = logsumexp(logW_c) - ln K.

Certificate: resample the K chains WITH replacement using the SAME index
draw for every class (CRN pairing preserved); per resample recompute
log Z_hat_b[lambda_hat] - log Z_hat_b[lambda] for every competitor lambda;
certify iff the lower (delta/C_comp) percentile of every difference exceeds
0, delta = 0.05, C_comp = C - 1 competitors (Bonferroni; contract pins
0.05/12 for the 13-class gross set). B >= 2500 resamples [DERIVED].
"""
from __future__ import annotations

import numpy as np

from scipy.special import logsumexp


def logZ_estimates(logW: np.ndarray) -> np.ndarray:
    """Per-class point estimates: logsumexp per class minus ln K. (C,)"""
    Kch = logW.shape[1]
    m = logW.max(axis=1)
    return m + np.log(np.exp(logW - m[:, None]).sum(axis=1)) - np.log(Kch)


def paired_bootstrap(logW: np.ndarray, B: int = 2500, delta: float = 0.05,
                     rng: np.random.Generator | None = None) -> dict:
    """Returns decision diagnostics. logW: (C, K).

    certified: bool -- lower Bonferroni percentile > 0 for ALL competitors.
    margin_nats: min over competitors of (log Z_hat[best] - log Z_hat[comp]).
    p_fail_worst: the worst (largest) competitor's non-certification
    probability = 1 - (bootstrap probability that the difference > 0).
    """
    rng = rng or np.random.default_rng(0)
    logW = np.asarray(logW, dtype=np.float64)
    C, Kch = logW.shape
    logZ = logZ_estimates(logW)
    best = int(np.argmax(logZ))
    competitors = [c for c in range(C) if c != best]
    # margin: gap to the runner-up in the point estimate
    sorted_lz = np.sort(logZ)
    margin = float(sorted_lz[-1] - sorted_lz[-2])

    # paired resampling: one index draw shared by all classes
    idx = rng.integers(0, Kch, size=(B, Kch))
    Wb = logW[:, idx]                                   # (C, B, Kch)
    m = Wb.max(axis=2)
    lz_b = m + np.log(np.exp(Wb - m[:, :, None]).sum(axis=2)) - np.log(Kch)
    diff = lz_b[best][None, :] - lz_b                   # (C, B) vs best

    level = delta / max(1, len(competitors))
    q = (100.0 * level)
    lo_q = np.percentile(diff[competitors], q, axis=1)  # lower percentile
    certified = bool(np.all(lo_q > 0))
    p_fail = float(np.mean(np.min(diff[competitors], axis=0) <= 0))
    return {
        "best": best,
        "margin_nats": margin,
        "certified": certified,
        "worst_lo_q": float(lo_q.min()),
        "bonf_level": level,
        "boot_p_fail": p_fail,
        "logZ": logZ,
    }
