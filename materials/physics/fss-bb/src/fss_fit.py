"""Finite-size-scaling (FSS) fitting per the exact ansatz of arXiv:2603.19062v3 Sec 2.5:

    WER(p, N) ~ f[(p - p_inf) * N^(1/nu)]

with the scaling function f approximated by a *degree-3 polynomial* in the
scaling variable x = (p - p_inf) * N^(1/nu). (p_inf, nu) are found by
minimizing the residual sum of squares of the polynomial fit.
Primary window: |p - p*(N)| < 0.06, where p*(N) is the pseudo-threshold
(WER = 0.10) of that size; sensitivity windows |p - p*| < 0.04 and < 0.08.
Uncertainty: bootstrap resampling of the WER data (paper: 500 datasets) —
parametric binomial resampling of per-point logical-error counts.

Independent companion estimator (arXiv:2603.19062v3 Sec 3.2 / Fig 2):
linearized FSS  p*(N) ~= p_inf + c * N^(-1/nu)  fit by least squares with
p*(N) values and their bootstrap CIs; intercept = p_inf.

Testable on synthetic data with known exponents (see validate_fss_module.py
in the campaign bundle): generate WER = f_syn((p-p0) N^{1/nu0}) + noise for a
fixed polynomial f_syn and recover (p0, nu0).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

WER_TARGET = 0.10  # pseudo-threshold definition: erasure rate where WER = 0.10 (paper Sec 2.4)


# --------------- FSS ansatz ---------------

def scaling_variable(p: np.ndarray, N: np.ndarray, p_inf: float, nu: float) -> np.ndarray:
    """x = (p - p_inf) * N^(1/nu), the argument of the universal scaling function."""
    p = np.asarray(p, float)
    N = np.asarray(N, float)
    return (p - p_inf) * N ** (1.0 / nu)


def fss_residuals(theta: np.ndarray, p, N, wer, ncoef: int = 4) -> np.ndarray:
    """Residuals of the deg-(ncoef-1) polynomial f(x) fit, over theta=(p_inf, nu, coefs).

    WER ~ c0 + c1 x + c2 x^2 + c3 x^3 with x = (p - p_inf) N^(1/nu).
    """
    p_inf, nu = theta[0], theta[1]
    coefs = theta[2:]
    x = scaling_variable(p, N, p_inf, nu)
    model = np.polynomial.polynomial.polyval(x, coefs)
    return model - wer



def _lstsq_rss(x: np.ndarray, y: np.ndarray, ncoef: int) -> tuple[float, np.ndarray]:
    V = np.vander(x, ncoef, increasing=True)
    c, *_ = np.linalg.lstsq(V, y, rcond=None)
    r = V @ c - y
    return float(r @ r), c


def _multistart_init(p, N, wer, pstar, window, ncoef, bounds):
    """Deterministic coarse grid scan over (p_inf, nu); best nested lstsq wins."""
    m = np.abs(p - pstar) < window
    pm, Nm, ym = p[m], N[m], wer[m]
    p_lo = max(bounds[0][0], float(pm.min() - 0.05))
    p_hi = min(bounds[0][1], float(pm.max() + 0.09))
    best = None
    for nu in np.linspace(bounds[1][0] + 0.1, 3.0, 9):
        for p_inf in np.linspace(p_lo, p_hi, 15):
            x = (pm - p_inf) * Nm ** (1.0 / nu)
            rss, c = _lstsq_rss(x, ym, ncoef)
            if best is None or rss < best[0]:
                best = (rss, p_inf, nu, c)
    return np.r_[best[1], best[2], best[3]]




def _nudged_starts(theta0, lb, ub):
    """Extra least_squares starts with nu nudged both ways (escape local minima)."""
    starts = []
    for d in (-0.35, 0.35):
        t = theta0.copy()
        t[1] = float(np.clip(theta0[1] + d, lb[1] + 1e-9, ub[1] - 1e-9))
        starts.append(t)
    return starts


def fit_fss_win(
    p,
    N,
    wer,
    pstar,
    window: float = 0.06,
    ncoef: int = 4,
    p_inf0: float = 0.5,
    nu0: float = 1.2,
    bounds=((0.0, 0.6), (0.3, 5.0)),
    multistart: bool = True,
):
    """Fit the paper's FSS ansatz on data with |p - p*(N)| < window.

    p, N, wer: per-point arrays (wer = logical error rate, not counts).
    pstar: per-point array of that size's pseudo-threshold (for windowing).
    Init: deterministic coarse (p_inf, nu) grid scan with nested polynomial
    lstsq; least_squares polish from the grid best plus two nu-nudged starts.
    Returns dict with p_inf, nu, coefs, rss, n_pts, and at_bound flags.
    """
    p = np.asarray(p, float)
    N = np.asarray(N, float)
    wer = np.asarray(wer, float)
    pstar = np.asarray(pstar, float)
    m = np.abs(p - pstar) < window
    if m.sum() < ncoef + 4:
        raise ValueError(f"window {window} selects only {m.sum()} points; need >= {ncoef + 4}")

    lb = np.r_[bounds[0][0], bounds[1][0], np.full(ncoef, -np.inf)]
    ub = np.r_[bounds[0][1], bounds[1][1], np.full(ncoef, np.inf)]

    def resid(t):
        return fss_residuals(t, p[m], N[m], wer[m], ncoef)

    if multistart:
        theta0 = _multistart_init(p, N, wer, pstar, window, ncoef, bounds)
    else:
        x0 = scaling_variable(p[m], N[m], p_inf0, nu0)
        V = np.vander(x0, ncoef, increasing=True)
        c0, *_ = np.linalg.lstsq(V, wer[m], rcond=None)
        theta0 = np.r_[p_inf0, nu0, c0]
    theta0 = np.clip(theta0, lb + 1e-12, ub - 1e-12)
    best = None
    for t_start in [theta0] + (_nudged_starts(theta0, lb, ub) if multistart else []):
        sol = least_squares(resid, t_start, bounds=(lb, ub), xtol=1e-12, ftol=1e-12, gtol=1e-12)
        r = float(np.sum(sol.fun**2))
        if best is None or r < best[0]:
            best = (r, sol)
    sol, rss = best[1], best[0]
    at_bound = [
        abs(sol.x[0]) < 1e-6 or abs(sol.x[0] - ub[0]) < 1e-6,
        abs(sol.x[1] - lb[1]) < 1e-6 or abs(sol.x[1] - ub[1]) < 1e-6,
    ]
    return {
        "p_inf": float(sol.x[0]),
        "nu": float(sol.x[1]),
        "coefs": np.asarray(sol.x[2:]),
        "rss": rss,
        "n_pts": int(m.sum()),
        "window": window,
        "win_mask": m,
        "at_bound": at_bound,
    }


def fit_fss_bootstrap(
    counts,
    shots,
    p,
    N,
    pstar,
    window: float = 0.06,
    ncoef: int = 4,
    n_boot: int = 500,
    seed: int = 12345,
    bounds=((0.0, 0.65), (0.3, 5.0)),
) -> dict:
    """Paper-style bootstrap CI: parametric binomial resampling of per-point
    logical-error counts (arXiv:2603.19062v3 Sec 2.5: 500 resampled datasets;
    Sec 2.4 uses 5000 for the pseudo-threshold CIs).

    counts, shots: integer arrays per point.
    Returns dict with point estimates and bootstrap CI percentiles for p_inf, nu.
    """
    counts = np.asarray(counts, float)
    shots = np.asarray(shots, float)
    wer = counts / shots
    base = fit_fss_win(p, N, wer, pstar, window=window, ncoef=ncoef, bounds=bounds)
    rng = np.random.default_rng(seed)
    pinf_b, nu_b = [], []
    ok = 0
    for _ in range(n_boot):
        k = rng.binomial(shots.astype(np.int64), np.clip(wer, 0.0, 1.0))
        try:
            r = fit_fss_win(p, N, k / shots, pstar, window=window, ncoef=ncoef, bounds=bounds)
        except Exception:
            continue
        pinf_b.append(r["p_inf"])
        nu_b.append(r["nu"])
        ok += 1
    out = dict(base)
    out["n_boot_ok"] = ok
    if ok:
        q = np.percentile
        out["p_inf_ci"] = (float(q(pinf_b, 2.5)), float(q(pinf_b, 97.5)))
        out["nu_ci"] = (float(q(nu_b, 2.5)), float(q(nu_b, 97.5)))
        out["p_inf_boot"] = np.asarray(pinf_b)
        out["nu_boot"] = np.asarray(nu_b)
    return out


# --------------- linearized companion estimator ---------------

def fit_linearized_fss(pstar, N, pstar_halfwidth=None, nu: float = 1.18):
    """p*(N) = p_inf + c * N^(-1/nu) via weighted least squares.

    Weights from the pseudo-threshold 95% CI half-widths (1.96 sigma),
    uniform when not provided. Returns dict with p_inf, slope c, and the
    intercept standard error.
    """
    x = np.asarray(N, float) ** (-1.0 / nu)
    y = np.asarray(pstar, float)
    if pstar_halfwidth is None:
        sigma = np.ones_like(y)
    else:
        sigma = np.maximum(np.asarray(pstar_halfwidth, float) / 1.96, 1e-9)
    w = 1.0 / sigma**2
    X = np.vstack([np.ones_like(x), x]).T
    XtWX = (X * w[:, None]).T @ X
    XtWX_inv = np.linalg.inv(XtWX)
    beta = XtWX_inv @ ((X * w[:, None]).T @ y)
    resid = y - X @ beta
    dof = max(len(y) - 2, 1)
    s2 = float((w * resid**2).sum() / dof)
    cov = XtWX_inv * s2
    return {
        "p_inf": float(beta[0]),
        "c": float(beta[1]),
        "p_inf_se": float(np.sqrt(cov[0, 0])),
        "c_se": float(np.sqrt(cov[1, 1])),
        "nu_used": nu,
    }


# --------------- pseudo-threshold search (paper Sec 2.4) ---------------

@dataclass
class BracketPoint:
    p: float
    counts: int = 0
    shots: int = 0

    @property
    def wer(self) -> float:
        return self.counts / self.shots if self.shots else float("nan")


def regula_falsi_illinois(fn, x0: float, x1: float, tol: float = 5e-4, max_iter: int = 10):
    """False-position with the Illinois modification (paper Sec 2.4 step 2).

    fn: WER(p) - 0.10, monotone increasing in p for the codes here.
    Returns (p*, (lo, hi)) where (lo, hi) is the final bracket.
    """
    f0, f1 = fn(x0), fn(x1)
    if f0 * f1 > 0:
        raise ValueError("bracket does not straddle the root")
    side = 0
    for _ in range(max_iter):
        # secant step
        p = x1 - f1 * (x1 - x0) / (f1 - f0)
        fp = fn(p)
        if abs(x1 - x0) < tol:
            return p, (min(x0, x1), max(x0, x1))
        if fp * f1 < 0:
            x0, f0 = x1, f1
            side = 0
        elif side == -1:
            f0 *= 0.5
        elif side == 1:
            f1 *= 0.5
        side = 1 if fp * f0 < 0 else side
        x1, f1 = p, fp
        if abs(x1 - x0) < tol:
            return p, (min(x0, x1), max(x0, x1))
    return p, (min(x0, x1), max(x0, x1))


def adaptive_pseudo_threshold(
    wer_fn, shots_per_eval: int = 200_000, start: float = 0.38, step: float = 0.04,
    bracket_tol: float = 5e-4, max_evals: int = 10,
):
    """Paper's 6-step adaptive search (Sec 2.4): bracket upward from p=0.38 in
    steps of 0.04 until WER > 0.10, then Illinois regula falsi to bracket width
    < 5e-4 (or 10 evaluations); report the final linear interpolant.

    wer_fn: callable p -> (counts, shots). Deterministic use expects the caller
    to bind a seed per evaluation point.
    Returns dict(p=..., bracket=..., n_evals=..., points=[BracketPoint...]).
    """
    points: list[BracketPoint] = []
    cache: dict[float, tuple[int, int]] = {}

    def ev(p):
        if p not in cache:
            cache[p] = wer_fn(p)
        return cache[p]

    def fn(p):
        c, s = ev(p)
        points.append(BracketPoint(p, c, s))
        return c / s - WER_TARGET

    # bracket from start upward; ev() is the single cache — never re-sample p
    p_lo = p_hi = None
    p = start
    n_evals = 0
    c, s = ev(p)
    if c / s <= WER_TARGET:
        p_lo = p
        while n_evals < max_evals:
            p += step
            c, s = ev(p)
            n_evals += 1
            if c / s > WER_TARGET:
                p_hi = p
                break
            p_lo = p
        if p_hi is None:
            raise RuntimeError("failed to bracket the pseudo-threshold with max_evals")
    else:
        # start already above 0.10: step backward (defensive; not hit in this study)
        p_hi = p
        while n_evals < max_evals:
            p -= step
            c, s = ev(p)
            n_evals += 1
            if c / s <= WER_TARGET:
                p_lo = p
                break
            p_hi = p
        if p_lo is None:
            raise RuntimeError("failed to bracket the pseudo-threshold with max_evals")
    lo, hi = min(p_lo, p_hi), max(p_lo, p_hi)
    # Illinois refinement; cached evals make fn() free on already-sampled p
    p_est, (blo, bhi) = regula_falsi_illinois(fn, lo, hi, tol=bracket_tol, max_iter=max_evals - n_evals)
    return {
        "p": float(p_est),
        "bracket": (blo, bhi),
        "n_evals": n_evals + len(points) - len(set(pt.p for pt in points)),
        "points": points,
    }


def pseudo_threshold_bootstrap(conf: dict, n_boot: int = 5000, seed: int = 12345):
    """Parametric bootstrap of the fitted pseudo-threshold (paper Sec 2.4:
    5000-iteration resampling from the binomial at the bracketing points).

    Returns (p_lo, p_hi) 95% CI by refitting the secant slope on resampled WERs.
    """
    pts = [pt for pt in conf["points"] if pt.shots > 0]
    if len(pts) < 2:
        raise ValueError("need >=2 evaluated points")
    ps = np.array([pt.p for pt in pts])
    ws = np.array([pt.wer for pt in pts])
    ns = np.array([pt.shots for pt in pts])
    # keep only the two tightest bracketing points (final interpolant)
    order = np.argsort(ps)
    ps, ws, ns = ps[order], ws[order], ns[order]
    i = np.searchsorted(ws, WER_TARGET)
    if i == 0:
        i = 1
    lo, hi = i - 1, i
    rng = np.random.default_rng(seed)
    est = []
    for _ in range(n_boot):
        kl = rng.binomial(ns[lo], ws[lo])
        kh = rng.binomial(ns[hi], ws[hi])
        wl, wh = kl / ns[lo], kh / ns[hi]
        if wh - wl < 1e-9 and (wh > WER_TARGET) == (wl > WER_TARGET):
            continue
        p = ps[lo] + (WER_TARGET - wl) * (ps[hi] - ps[lo]) / (wh - wl)
        est.append(p)
    est = np.asarray(est)
    return (float(np.percentile(est, 2.5)), float(np.percentile(est, 97.5))), est
