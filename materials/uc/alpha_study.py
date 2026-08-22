"""Alpha-margin study: is alpha = 0.0356069 the best choice for certification?

NUMERICAL GUIDANCE ONLY (float64 sampling + the exact obstruction family).
The certified margin at level t is min Lambda over the feasible 2-pair-orbit
family; a larger margin makes interval certification easier.  This script maps
margin(alpha) at t = psi+1e-4 and t = c*-2e-4 to check whether the inherited
alpha is near-optimal.  No universal claims: minima are over samples plus the
exact 1-D obstruction family.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from entropy import PSI  # noqa: E402

LOG2 = np.log(2.0)
CSTAR = 0.38234553336670272115


def hf(x):
    x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
    o = np.zeros_like(x)
    m = (x > 0) & (x < 1)
    xm = x[m]
    o[m] = -(xm * np.log(xm) + (1 - xm) * np.log(1 - xm)) / LOG2
    return o


def rhf(z):
    z = np.asarray(z, dtype=float)
    return np.clip(2 * (1 - z) * hf(z) - hf((1 - z) ** 2), 0.0, None)


def sstar(a, b):
    return np.minimum(np.maximum(0.5, np.maximum(a, b)), np.minimum(a + b, 1.0))


def lam_vec(p1, q1, p2, q2, w, alpha):
    M = w * (p1 + q1) / 2 + (1 - w) * (p2 + q2) / 2
    L = (w * (hf(p1) + hf(q1)) / 2 + (1 - w) * (hf(p2) + hf(q2)) / 2)
    C = w * hf(sstar(p1, q1)) + (1 - w) * hf(sstar(p2, q2))
    S = (w * (np.sqrt(rhf(p1)) + np.sqrt(rhf(q1))) / 2
         + (1 - w) * (np.sqrt(rhf(p2)) + np.sqrt(rhf(q2))) / 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = 2 * (1 - alpha) * (1 - M) - 1 + alpha * C / L - (1 - alpha) * S * S / L
    return M, L, lam


def polish(x0, alpha, t, iters=3000, seed=1):
    """Crude coordinate random descent from x0, staying feasible."""
    rng = np.random.default_rng(seed)
    x = np.array(x0, dtype=float)

    def val(x):
        M, L, lam = lam_vec(*[np.array([v]) for v in x], alpha)
        if M[0] > t or L[0] <= 1e-14 or not np.isfinite(lam[0]):
            return np.inf
        return lam[0]

    v = val(x)
    step = 0.05
    for i in range(iters):
        if i % 600 == 599:
            step *= 0.3
        cand = np.clip(x + rng.normal(0, step, 5), 0, 1)
        cv = val(cand)
        if cv < v:
            x, v = cand, cv
    return v, x


def min_lambda(alpha, t, n=600_000, seed=0):
    rng = np.random.default_rng(seed)
    p1 = rng.uniform(0, 1, n); q1 = rng.uniform(0, 1, n)
    p2 = rng.uniform(0, 1, n); q2 = rng.uniform(0, 1, n)
    w = rng.uniform(0, 1, n)
    m = n // 3                       # focused cloud near the obstruction
    b0 = 0.3294547385
    p1[:m] = np.clip(b0 + rng.normal(0, 0.05, m), 0, 1)
    q1[:m] = np.clip(b0 + rng.normal(0, 0.05, m), 0, 1)
    p2[:m] = np.clip(b0 + rng.normal(0, 0.05, m), 0, 1)
    q2[:m] = np.clip(1.0 - np.abs(rng.normal(0, 0.05, m)), 0, 1)
    w[:m] = np.clip(0.84 + rng.normal(0, 0.08, m), 0, 1)
    M, L, lam = lam_vec(p1, q1, p2, q2, w, alpha)
    ok = (M <= t) & (L > 1e-12) & np.isfinite(lam)
    i = np.argmin(np.where(ok, lam, np.inf))
    best_rand = lam[i]
    # polish the best random hit
    v_pol, _ = polish([p1[i], q1[i], p2[i], q2[i], w[i]], alpha, t)
    # exact obstruction family (b,b,b,1,w), mean pinned at t
    b = np.linspace(0.02, t - 1e-6, 6000)
    a = (t - b) / (1 - b)
    wv = 1 - 2 * a
    g = (wv > 0) & (wv < 1)
    _, _, lo = lam_vec(b[g], b[g], b[g], np.ones_like(b[g]), wv[g], alpha)
    best_obs = np.nanmin(lo)
    kappa = 2 * (1 - alpha) * (1 - t) - 1        # sink-limit Lambda
    return min(best_rand, v_pol), best_obs, kappa


for lbl, t in (("psi+1e-4", float(PSI) + 1e-4), ("c*-2e-4", CSTAR - 2e-4)):
    print(f"\nt = {lbl} = {t:.12f}")
    print(f"{'alpha':>9} {'min sampled':>12} {'min obs-fam':>12} "
          f"{'sink kappa':>12} {'overall':>12}")
    rows = []
    for alpha in (0.0144, 0.02, 0.0356069, 0.05, 0.07, 0.09, 0.12):
        br, bo, ka = min_lambda(alpha, t)
        overall = min(br, bo, ka)
        rows.append((alpha, overall, br, bo, ka))
        print(f"{alpha:>9.5f} {br:>12.6f} {bo:>12.6f} {ka:>12.6f} {overall:>12.6f}")
    best = max(rows, key=lambda r: r[1])
    cur = [r for r in rows if abs(r[0] - 0.0356069) < 1e-9][0]
    print(f"best-in-grid alpha {best[0]}: margin {best[1]:.6f}; "
          f"current alpha margin {cur[1]:.6f}")
    assert best[1] >= cur[1] - 1e-12
print("\nSampled minima only -- guidance for choosing alpha, not a bound.")
