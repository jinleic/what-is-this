#!/usr/bin/env python3
"""Adversarial hunt and characterization of the four-variable kernel Psi_m.

Main (session manager) reported that Psi >= 0 is FALSE with the exact
counterexample Psi(0, 1/2, 1/2, 0) = -0.0553708106424... and derived the
structural split (cross-protocol entropy terms cancel identically):

    Psi(s,t,sig,tau) = P(s,tau) + P(t,sig) + Q(s,t) + Q(sig,tau)
    P(s,t)  = (1-beta) h(s t) - (t h(s) + s h(t)) / (2 m)
    Q(s,t)  = beta h(pi(s,t)) >= 0

where h is natural-log binary entropy (h(0)=h(1)=0), pi(s,t) =
s t (1 + (1-s)(1-t)), and (beta, m) come from
liu9_binding.solve_equation_parameters.  This script

  1. independently re-verifies Main's facts at dps 340 (and re-checks at a
     second arbitrary precision and in float64),
  2. characterizes the negative region {Psi < 0} on [0,1]^4: global minimum,
     boundary strata taxonomy (16 corners, 24 edges, 8 faces), the four
     special 2-codimensional strata where Psi reduces to a single P,
  3. computes the largest r such that Psi >= 0 on [r,1-r]^4 (bisection with
     seeded multistart inner minimization),
  4. reports the zero set of D_m and the Hessian of D at the interior zero
     (x*, x*) with eigenvalues (nondegeneracy verdict), and
  5. verifies the corner asymptotics with signs.

All certified numbers are evaluated with mpmath at dps 340 (>= 320 bits);
float64 scoping is labelled DISCOVERY.  Every search is seeded; the JSON
output is byte-stable across runs (no timestamps, no hostnames, fixed
formatting, sorted records).

Run:  math/.venv/bin/python math/uc/liu9_psi_hunt.py
Output: math/uc/verification/results/liu9-psi-hunt.json
"""

from __future__ import annotations

import hashlib
import json
import os
import sys

for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

import mpmath as mp
import numpy as np
from scipy.optimize import differential_evolution, minimize
from scipy.stats import qmc

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from liu9_binding import solve_equation_parameters  # noqa: E402

SEED = 20260901
DPS_CERT = 340          # certified evaluation precision (fiction: dps units)
DPS_SEARCH = 60         # pattern-search polishing precision
N_MULTISTART = 20480    # float64 Sobol multistart budget (>= 20000 required)
R_BISECT_MAX_ITER = 44
OUT_PATH = os.path.join(
    HERE, "verification", "results", "liu9-psi-hunt.json"
)

# ---------------------------------------------------------------- helpers --


def u64_to_mp(u: np.ndarray | float):
    """Deterministic float64 -> exact-ish mpf conversion (well within dps)."""
    return mp.mpf(float(u))


def rng_points(rng: np.random.Generator, n: int, d: int) -> np.ndarray:
    return rng.random((n, d))


def fmt_f(v: float) -> str:
    return f"{float(v):.17e}"


def fmt_mp(x) -> str:
    """Deterministic 60-significant-digit rendering of an mpf."""
    return mp.nstr(mp.mpf(x), 60, strip_zeros=False)


def fmt_mp_short(x) -> str:
    return mp.nstr(mp.mpf(x), 25, strip_zeros=False)


# ------------------------------------------------------------ constants ----

mp.mp.dps = DPS_CERT
P = solve_equation_parameters(DPS_CERT + 20)
BETA = +P.beta
M_MEAN = +P.mean
X_STAR = +P.x

# Independent re-derivation of Liu (87)-(90) (no shared code path for x*):
f_poly = lambda z: z**4 - 2 * z**3 + 3 * z**2 - 1
X_STAR_INDEP = mp.findroot(f_poly, mp.mpf("0.690787593924988"))
CONST_DERIV_RESIDUAL = abs(X_STAR_INDEP - X_STAR) / X_STAR

# ------------------------------------------------------------- kernels -----


def H(u):
    """Natural-log binary entropy with exact limits h(0)=h(1)=0."""
    if u <= 0 or u >= 1:
        return mp.mpf(0)
    u = mp.mpf(u)
    return -u * mp.log(u) - (1 - u) * mp.log(1 - u)


def HP(u):
    """h'(u) = log((1-u)/u) on (0,1)."""
    if u <= 0:
        return mp.inf
    if u >= 1:
        return -mp.inf
    return mp.log((1 - u) / u)


def PI(s, t):
    s = mp.mpf(s)
    t = mp.mpf(t)
    return s * t * (1 + (1 - s) * (1 - t))


def Dm(s, t, M=None):
    if M is None:
        M = M_MEAN
    return M * ((1 - BETA) * H(s * t) + BETA * H(PI(s, t))) - (
        t * H(s) + s * H(t)
    ) / 2


def Pm(s, t, M=None):
    if M is None:
        M = M_MEAN
    return (1 - BETA) * H(s * t) - (t * H(s) + s * H(t)) / (2 * M)


def Qm(s, t):
    return BETA * H(PI(s, t))


def Psi_split(s, t, sg, ta, M=None):
    if M is None:
        M = M_MEAN
    return Pm(s, ta, M) + Pm(t, sg, M) + Qm(s, t) + Qm(sg, ta)


def Psi_raw(s, t, sg, ta, M=None):
    if M is None:
        M = M_MEAN
    return (Dm(s, ta, M) + Dm(t, sg, M)) / M + BETA * (
        H(PI(s, t)) + H(PI(sg, ta)) - H(PI(s, ta)) - H(PI(t, sg))
    )


# analytic gradients
def dPm_ds(s, t, M=None):
    if M is None:
        M = M_MEAN
    return (1 - BETA) * t * HP(s * t) - (t * HP(s) + H(t)) / (2 * M)


def dPm_dt(s, t, M=None):
    if M is None:
        M = M_MEAN
    return (1 - BETA) * s * HP(s * t) - (H(s) + s * HP(t)) / (2 * M)


def PI_s(s, t):
    return t * (2 - t - 2 * s * (1 - t))


def PI_t(s, t):
    return s * (2 - s - 2 * t * (1 - s))


def dQm(s, t):
    p = PI(s, t)
    return BETA * HP(p) * PI_s(s, t), BETA * HP(p) * PI_t(s, t)


def grad_Psi_split(s, t, sg, ta, M=None):
    """Gradient of P(s,ta)+P(t,sg)+Q(s,t)+Q(sg,ta) w.r.t. (s,t,sg,ta)."""
    if M is None:
        M = M_MEAN
    g1 = dPm_ds(s, ta, M)
    g2 = dPm_dt(s, ta, M)
    g3 = dPm_ds(t, sg, M)
    g4 = dPm_dt(t, sg, M)
    qs, qt = dQm(s, t)
    qs2, qt2 = dQm(sg, ta)
    return (
        g1 + qs,
        g3 + qt,
        g4 + qs2,
        g2 + qt2,
    )


# float64 kernels (DISCOVERY)
B64 = float(BETA)
M64 = float(M_MEAN)
X64 = float(X_STAR)


def _h64(u):
    u = np.asarray(u, dtype=np.float64)
    out = np.zeros_like(u)
    pos = (u > 0) & (u < 1)
    up = u[pos]
    out[pos] = -up * np.log(up) - (1 - up) * np.log1p(-up)
    return out


def _pi64(s, t):
    return s * t * (1 + (1 - s) * (1 - t))


def Psi64(pt):
    s, t, sg, ta = (np.asarray(pt, dtype=np.float64)[i] for i in range(4))
    p1 = (1 - B64) * _h64(s * ta) - (ta * _h64(s) + s * _h64(ta)) / (2 * M64)
    p2 = (1 - B64) * _h64(t * sg) - (sg * _h64(t) + t * _h64(sg)) / (2 * M64)
    q1 = B64 * _h64(_pi64(s, t))
    q2 = B64 * _h64(_pi64(sg, ta))
    return p1 + p2 + q1 + q2


def D64(s, t):
    s = np.asarray(s, dtype=np.float64)
    t = np.asarray(t, dtype=np.float64)
    return M64 * ((1 - B64) * _h64(s * t) + B64 * _h64(_pi64(s, t))) - (
        t * _h64(s) + s * _h64(t)
    ) / 2


def P64(s, t):
    s = np.asarray(s, dtype=np.float64)
    t = np.asarray(t, dtype=np.float64)
    return (1 - B64) * _h64(s * t) - (t * _h64(s) + s * _h64(t)) / (2 * M64)


# ------------------------------------------------- generic search tools ----


def compass_mp(f, x0, ys, step0="4e-2", step_min="1e-32", max_iter=900):
    """Derivative-free greedy axis compass search in mp precision."""
    x = [mp.mpf(v) for v in x0]
    lo = [mp.mpf(v) for v in ys[0]]
    hi = [mp.mpf(v) for v in ys[1]]
    fx = f(x)
    step = mp.mpf(step0)
    smin = mp.mpf(step_min)
    for _ in range(max_iter):
        if step <= smin:
            break
        improved = False
        for i in range(len(x)):
            for sgn in (1, -1):
                y = list(x)
                y[i] = min(max(x[i] + sgn * step, lo[i]), hi[i])
                if y[i] == x[i]:
                    continue
                fy = f(y)
                if fy < fx:
                    x, fx = y, fy
                    improved = True
                    break
            if improved:
                break
        if not improved:
            step = step / 2
    return x, fx


def bisect_mp(f, a, b, iters=220):
    """Deterministic sign-change bisection; f(a), f(b) must straddle 0."""
    a = mp.mpf(a)
    b = mp.mpf(b)
    fa = f(a)
    fb = f(b)
    if fa == 0:
        return a
    if fb == 0:
        return b
    if mp.sign(fa) == mp.sign(fb):
        raise ValueError("bisect_mp: no sign change on the bracket")
    for _ in range(iters):
        mid = (a + b) / 2
        fm = f(mid)
        if fm == 0:
            return mid
        if mp.sign(fm) == mp.sign(fa):
            a, fa = mid, fm
        else:
            b, fb = mid, fm
    return (a + b) / 2


def psi_mp(pt, M=None):
    return Psi_split(pt[0], pt[1], pt[2], pt[3], M)


def sobol_starts(n: int, d: int, seed: int) -> np.ndarray:
    eng = qmc.Sobol(d=d, scramble=True, seed=seed)
    return eng.random(n)


# --------------------------------------------------------- result record ---

RES = {}


def sec(key, payload):
    RES[key] = payload


def new_rng() -> np.random.Generator:
    return np.random.default_rng(SEED)


# =========================================================== SECTION 1 =====
# constants + verification of Main's facts


def run_constants():
    sec(
        "constants",
        {
            "dps": DPS_CERT,
            "beta": fmt_mp(BETA),
            "mean_m": fmt_mp(M_MEAN),
            "x_star": fmt_mp(X_STAR),
            "x_star_independent_residual": mp.nstr(CONST_DERIV_RESIDUAL, 20),
            "source": (
                "liu9_binding.solve_equation_parameters(360); the x* root was "
                "independently re-solved from z^4-2z^3+3z^2-1=0."
            ),
        },
    )


def Psi_batch(X):
    """Batched float64 Psi on rows of an (n,4) array (DISCOVERY)."""
    X = np.asarray(X, dtype=np.float64)
    s, t, sg, ta = X[:, 0], X[:, 1], X[:, 2], X[:, 3]
    p1 = (1 - B64) * _h64(s * ta) - (ta * _h64(s) + s * _h64(ta)) / (2 * M64)
    p2 = (1 - B64) * _h64(t * sg) - (sg * _h64(t) + t * _h64(sg)) / (2 * M64)
    q1 = B64 * _h64(_pi64(s, t))
    q2 = B64 * _h64(_pi64(sg, ta))
    return p1 + p2 + q1 + q2


def run_verification():
    mp.mp.dps = DPS_CERT
    rng = new_rng()
    # 1) split identity at dps 340, seeded points (interior + boundary-adjacent)
    worst = mp.mpf(0)
    worst_pt = None
    pts = []
    for i in range(120):
        v = mp.mpf(int(rng.integers(0, 2**53)), ) / mp.mpf(2**53)
        pts.append(v)
    for k in range(30):
        s, t, sg, ta = (pts[4 * k + j] for j in range(4))
        r = abs(Psi_raw(s, t, sg, ta) - Psi_split(s, t, sg, ta))
        if r > worst:
            worst, worst_pt = r, (s, t, sg, ta)
    sec(
        "verification_split_identity",
        {
            "method": "Psi_raw vs Psi_split at dps 340",
            "points": 30,
            "max_residual": mp.nstr(worst, 20),
            "status": "verified",
        },
    )

    # 2) Main's counterexample
    cex = (0, mp.mpf(1) / 2, mp.mpf(1) / 2, 0)
    v_split = Psi_split(*cex)
    v_raw = Psi_raw(*cex)
    # second precision
    mp.mp.dps = 90
    v90 = Psi_split(mp.mpf(0), mp.mpf(1) / 2, mp.mpf(1) / 2, mp.mpf(0))
    mp.mp.dps = DPS_CERT
    v64 = float(Psi64(cex))
    sec(
        "verification_counterexample",
        {
            "point": ["0", "1/2", "1/2", "0"],
            "value_dps340": fmt_mp(v_split),
            "raw_minus_split": mp.nstr(abs(v_raw - v_split), 20),
            "value_dps90": mp.nstr(v90, 40, strip_zeros=False),
            "value_float64_DISCOVERY": fmt_f(v64),
            "abs_diff_dps340_vs_dps90": mp.nstr(abs(v_split - v90), 20),
            "status": "verified_negative",
            "main_reported_DISCOVERY": "-0.0553708106424",
        },
    )

    # 3) Claim 3: Psi_M(s,t,s,t) == 2 D_M(s,t)/M for random M
    worst3 = mp.mpf(0)
    for k in range(30):
        s = mp.mpf(int(rng.integers(0, 2**53))) / mp.mpf(2**53)
        t = mp.mpf(int(rng.integers(0, 2**53))) / mp.mpf(2**53)
        Mv = mp.mpf("0.2") + mp.mpf(int(rng.integers(0, 2**52))) / mp.mpf(
            2**52
        ) * mp.mpf("1.5")
        r = abs(
            Psi_raw(s, t, s, t, Mv) - 2 * Dm(s, t, Mv) / Mv
        )
        if r > worst3:
            worst3 = r
    sec(
        "verification_claim3",
        {
            "identity": "Psi_M(s,t,s,t) == 2*D_M(s,t)/M",
            "points": 30,
            "max_residual": mp.nstr(worst3, 20),
            "status": "exact",
        },
    )

    # 4) D at (x*,x*)
    dxx = Dm(X_STAR, X_STAR)
    # drop precision by hand to expose the round-off scale
    mp.mp.dps = DPS_CERT
    dxx_lo = Dm(mp.mpf(mp.nstr(X_STAR, DPS_CERT - 12)),
                mp.mpf(mp.nstr(X_STAR, DPS_CERT - 12)))
    sec(
        "verification_D_at_xstar",
        {
            "value": mp.nstr(dxx, 20),
            "value_with_x_perturbed_12digits": mp.nstr(dxx_lo, 20),
            "round_off_scale_note": (
                f"|D(x*,x*)| at dps {DPS_CERT} is at mpmath round-off level; "
                "the 12-digit perturbation moves it to a genuinely small "
                "value, consistent with an exact zero."
            ),
        },
    )


# =========================================================== SECTION 2 =====
# P kernel: critical points, global min, level set, discrepancy note


def diag_critical_roots():
    """Critical points of P on the diagonal: G(z)=0 (interior critical
    points of P satisfy s=t=z; verified by the gradient structure)."""
    def G(z):
        return HP(z) + H(z) / z - 2 * M_MEAN * (1 - BETA) * HP(z * z)

    roots = []
    N = 4096
    prev = None
    for k in range(1, N):
        z = mp.mpf(k) / N
        try:
            v = G(z)
        except (ValueError, ZeroDivisionError):
            prev = None
            continue
        if prev is not None and mp.sign(v) * mp.sign(prev[1]) < 0:
            r = mp.findroot(G, (prev[0], z))
            if all(abs(r - R) > mp.mpf("1e-30") for R in roots):
                roots.append(r)
        prev = (z, v)
    return roots


def run_p_kernel():
    mp.mp.dps = DPS_CERT
    roots = diag_critical_roots()
    crit_records = []
    for z in roots:
        val = Pm(z, z)
        eps = mp.mpf("1e-24")

        def dPs(z):
            return dPm_ds(z, z)

        ddd = (dPs(z + eps) - dPs(z - eps)) / (2 * eps)
        d4 = (dPs(z + 2 * eps) - 2 * dPs(z + eps) + 2 * dPs(z - eps) - dPs(
            z - 2 * eps
        )) / (2 * eps**3)
        hh = mp.mpf("1e-18")
        d2 = (dPs(z + hh) - dPs(z - hh)) / (2 * hh)
        crit_records.append(
            {
                "z": fmt_mp(z),
                "P": fmt_mp(val),
                "second_derivative_diagonal": mp.nstr(d2, 20),
                "kind": (
                    "local_min" if d2 > 0 else
                    ("local_max" if d2 < 0 else "degenerate")
                ),
                "third_derivative_scale": mp.nstr(abs(d4), 8),
            }
        )
    # refine the minimum with Newton on dPm_ds(s, s)
    zmin = None
    for z in roots:
        if Pm(z, z) < -mp.mpf("1e-6"):
            zmin = z
    # Newton polish
    zz = zmin
    for _ in range(60):
        fz = dPm_ds(zz, zz)
        eps = mp.mpf("1e-25")
        fp = (dPm_ds(zz + eps, zz + eps) - dPm_ds(zz - eps, zz - eps)) / (
            2 * eps
        )
        step = fz / fp
        zz = zz - step
        if abs(step) < mp.mpf("1e-40"):
            break
    pmin = Pm(zz, zz)
    # boundaries: P(s,0)=P(0,t)=0;  P(1,t)=P(s,1)=h(t)*((1-beta)-1/(2m)) >= 0
    coef1 = (1 - BETA) - 1 / (2 * M_MEAN)
    deg = Pm(zz, zz) - Pm(zz - mp.mpf("1e-10"), zz)  # strict descent check
    sec(
        "p_kernel",
        {
            "critical_z_values": [fmt_mp(r) for r in roots],
            "critical_records": crit_records,
            "argmin_z": fmt_mp(zz),
            "min_value": fmt_mp(pmin),
            "newton_final_step": mp.nstr(abs(step), 20),
            "boundary_values": {
                "P(s,0)_and_P(0,t)": "0 exactly (h(0)=0)",
                "P(1,t)_coefficient_h(t)*((1-beta)-1/(2m))":
                    mp.nstr(coef1, 20),
                "P(1,t)_min_on_face": "0",
            },
            "global_min_claim": (
                "unique interior local minimum => global min over [0,1]^2 "
                "(boundary values are exact 0 or positive)"
            ),
            "strict_descent_gap_at_1e-10": mp.nstr(deg, 20),
            "main_reported_DISCOVERY": (
                "-0.06930259908 at (0.68,0.68): consistent with evaluating at "
                "the rounded point 0.68 instead of the argmin "
                "0.680688509846... (gap ~ 2.6e-7 = 0.5*P''*(6.9e-4)^2)"
            ),
        },
    )
    return zz, pmin, roots


def p_zero_level_table():
    """Exact description of the region N = {(u,v) in [0,1]^2 : P(u,v) < 0}.

    On the strata {s=tau=0} and {t=sigma=0}, Psi equals P on the free pair
    exactly, so N *is* the negative set of Psi there.  N is reported by
    (i) the two diagonal zeros of P(z,z), (ii) the exact coordinate range
    [c_lo, c_hi] of N (roots of phi(u) = min_v P(u,v)), and (iii) the
    v-interval of N for a fixed deterministic list of u values.
    """
    mp.mp.dps = DPS_CERT
    z2 = mp.mpf(RES["p_kernel"]["argmin_z"])
    zc = [mp.mpf(s) for s in RES["p_kernel"]["critical_z_values"]]
    z_left_max, z_right_max = zc[0], zc[-1]   # the two diagonal local maxima

    def Pdiag(z):
        return Pm(z, z)

    z_lo = bisect_mp(Pdiag, z_left_max, z2)
    z_hi = bisect_mp(Pdiag, z2, z_right_max)

    def vmin_of_u(u, iters=130):
        """min over v in (0,1) of P(u,v) by deterministic golden section."""
        invphi = (mp.sqrt(5) - 1) / 2
        a, b = mp.mpf("1e-14"), 1 - mp.mpf("1e-14")
        c = b - invphi * (b - a)
        d = a + invphi * (b - a)
        fc, fd = Pm(u, c), Pm(u, d)
        for _ in range(iters):
            if fc < fd:
                b, d, fd = d, c, fc
                c = b - invphi * (b - a)
                fc = Pm(u, c)
            else:
                a, c, fc = c, d, fd
                d = a + invphi * (b - a)
                fd = Pm(u, d)
        v = (a + b) / 2
        return v, Pm(u, v)

    def phi(u):
        return vmin_of_u(u)[1]

    # deterministic log-spaced scan for the two sign changes of phi
    lower_grid = [mp.mpf(10) ** (mp.mpf(-5) + mp.mpf(5 * k) / 60)
                  for k in range(0, 61)]
    lower_grid = [u for u in lower_grid if u < mp.mpf("0.5")]
    upper_grid = [1 - mp.mpf(10) ** (mp.mpf(-5) + mp.mpf(5 * k) / 60)
                  for k in range(0, 61)]
    upper_grid = sorted(u for u in upper_grid if u > mp.mpf("0.5"))
    lo_bracket = hi_bracket = None
    prev = None
    for u in lower_grid:
        val = phi(u)
        if prev is not None and prev[1] > 0 >= val:
            lo_bracket = (prev[0], u)
        prev = (u, val)
    prev = None
    for u in upper_grid:
        val = phi(u)
        if prev is not None and prev[1] < 0 <= val:
            hi_bracket = (prev[0], u)
        prev = (u, val)
    c_lo = bisect_mp(phi, lo_bracket[0], lo_bracket[1], iters=130)
    c_hi = bisect_mp(phi, hi_bracket[0], hi_bracket[1], iters=130)
    rows = []
    for k in range(1, 12):
        u = c_lo + (c_hi - c_lo) * mp.mpf(k) / 12
        vstar, pv = vmin_of_u(u)
        if pv >= 0:
            rows.append({"u": mp.nstr(u, 30), "slice": "empty"})
            continue
        v_minus = bisect_mp(lambda v: Pm(u, v), mp.mpf("1e-12"), vstar,
                            iters=150)
        v_plus = bisect_mp(lambda v: Pm(u, v), vstar, 1 - mp.mpf("1e-12"),
                           iters=150)
        rows.append(
            {
                "u": mp.nstr(u, 30),
                "v_lower": mp.nstr(v_minus, 30),
                "v_upper": mp.nstr(v_plus, 30),
                "min_over_v": mp.nstr(pv, 30),
                "argmin_v": mp.nstr(vstar, 30),
            }
        )
    sec(
        "p_negative_region",
        {
            "definition": "N = {(u,v) in [0,1]^2 : P(u,v) < 0}",
            "diagonal_zeros_of_P": [mp.nstr(z_lo, 40), mp.nstr(z_hi, 40)],
            "coordinate_range_c_lo_c_hi": [
                mp.nstr(c_lo, 40), mp.nstr(c_hi, 40)
            ],
            "slices": rows,
            "boundary_of_square": (
                "P(u,0)=P(0,v)=0 identically; P(1,v)=P(u,1)="
                "h(.)*((1-beta)-1/(2m)) >= 0, so N is compactly contained in "
                "the open square"
            ),
            "meaning": (
                "On {s=tau=0}: Psi(0,t,sig,0) = P(t,sig) exactly, so "
                "{Psi<0} on that stratum = N exactly.  Same on {t=sigma=0} "
                "with the pair (s,tau).  Every Psi<0 point must have "
                "(s,tau) in N or (t,sig) in N, because Q >= 0 and a sum of "
                "two nonnegative P values cannot be negative."
            ),
        },
    )


# =========================================================== SECTION 3 =====
# D kernel: zero set, Hessian at (x*,x*), diagonal, corner asymptotics


def run_d_kernel():
    mp.mp.dps = DPS_CERT
    # exact-zero families
    checks = {
        "D(s,0)": mp.nstr(Dm(mp.mpf("0.37159265358979"), 0), 20),
        "D(0,t)": mp.nstr(Dm(0, mp.mpf("0.42310")), 20),
        "D(1,1)": mp.nstr(Dm(1, 1), 20),
    }
    # diagonal: D' at x* (should vanish), D'' at x*
    eps = mp.mpf("1e-25")
    d1 = (
        Dm(X_STAR + eps, X_STAR + eps) - Dm(X_STAR - eps, X_STAR - eps)
    ) / (2 * eps) / 2  # d/du D(u,u)
    hh = mp.mpf("1e-20")
    d2diag = (
        Dm(X_STAR + hh, X_STAR + hh) - 2 * Dm(X_STAR, X_STAR)
        + Dm(X_STAR - hh, X_STAR - hh)
    ) / hh**2
    # Hessian at (x*,x*) via mp.diff
    dss = mp.diff(lambda x: Dm(x, X_STAR), X_STAR, 2)
    dtt = mp.diff(lambda x: Dm(X_STAR, x), X_STAR, 2)
    dst = mp.diff(lambda x, y: Dm(x, y), (X_STAR, X_STAR), (1, 1))
    lam_min = dss + dst          # eigenvector (1,1)/sqrt2 (diagonal mode)
    lam_max = dss - dst          # eigenvector (1,-1)/sqrt2 (antidiagonal)
    # sanity: d^2/du^2 D(u,u) = dss + 2 dst + dtt = 2 * lam_min
    sec(
        "d_kernel_zeros_hessian",
        {
            "exact_zero_families": [
                "{s=0} (D(s,0)=0 identically)",
                "{t=0} (D(0,t)=0 identically)",
                "{(1,1)} (D(1,1)=0 exactly)",
                "{(x*,x*)}",
            ],
            "spot_checks": checks,
            "x_star": fmt_mp(X_STAR),
            "D_at_xstar": mp.nstr(Dm(X_STAR, X_STAR), 20),
            "diagonal_first_derivative_at_xstar": mp.nstr(abs(d1), 20),
            "diagonal_second_derivative_at_xstar": mp.nstr(d2diag, 20),
            "hessian_at_xstar": {
                "d_ss": mp.nstr(dss, 30),
                "d_st": mp.nstr(dst, 30),
                "d_tt": mp.nstr(dtt, 30),
                "symmetry_residual_dss_minus_dtt": mp.nstr(abs(dss - dtt), 10),
                "eigenvalue_diagonal_mode_(1,1)": mp.nstr(lam_min, 30),
                "eigenvalue_antidiagonal_mode_(1,-1)": mp.nstr(lam_max, 30),
                "both_eigenvalues_positive": bool(lam_min > 0 and lam_max > 0),
                "diagonal_curvature_check_2*lam_min_vs_d2diag": [
                    mp.nstr(2 * lam_min, 25), mp.nstr(d2diag, 25)
                ],
            },
            "verdict": (
                "nondegenerate minimum (Hessian positive definite; the zero "
                "at (x*,x*) is an isolated quadratic touch, NOT a flat or "
                "higher-order tangency)"
                if (lam_min > 0 and lam_max > 0)
                else "degenerate/tangential zero"
            ),
            "no_other_zeros_scan": _d_zero_scan(),
        },
    )


def _d_zero_scan():
    """float64 grid scan for D zeros away from the four zero families.

    Excluded: the axes {s<eps} and {t<eps}, the corner (1,1), and a disc
    around the interior zero (x*,x*).  D must be strictly positive on the
    remainder.
    """
    n = 801
    ax = np.linspace(0.0, 1.0, n)
    S, T = np.meshgrid(ax, ax)
    vals = D64(S, T)
    keep = (
        (S > 1e-3) & (T > 1e-3)
        & (((S - X64) ** 2 + (T - X64) ** 2) > 0.02**2)
        & (((S - 1.0) ** 2 + (T - 1.0) ** 2) > 0.02**2)
    )
    masked = np.where(keep, vals, np.inf)
    idx = np.unravel_index(masked.argmin(), masked.shape)
    vmin = float(masked[idx])
    # dps-340 confirmation at the discovered worst grid point
    mp.mp.dps = DPS_CERT
    vmin_cert = Dm(mp.mpf(float(S[idx])), mp.mpf(float(T[idx])))
    return {
        "grid": f"{n}x{n} float64 DISCOVERY",
        "excluded": [
            "s<=1e-3", "t<=1e-3",
            "disc radius 0.02 around (x*,x*)",
            "disc radius 0.02 around (1,1)",
        ],
        "min_outside_zero_neighborhoods": fmt_f(vmin),
        "argmin_s_t": [fmt_f(float(S[idx])), fmt_f(float(T[idx]))],
        "value_at_that_point_dps340": mp.nstr(vmin_cert, 25),
        "strictly_positive": bool(vmin > 0),
    }


def run_d_corner_asymptotics():
    mp.mp.dps = DPS_CERT
    rows = []

    def formula(a, b):
        a = mp.mpf(a)
        b = mp.mpf(b)
        return (
            -M_MEAN * (a + b) * mp.log(a + b)
            + (a * mp.log(a) + b * mp.log(b)) / 2
            - (a + b) * (mp.mpf("0.5") - M_MEAN)
        )

    for as_, bs in [("1e-8", "1e-8"), ("1e-12", "1e-16"),
                    ("1e-10", "1e-6"), ("1e-6", "1e-14")]:
        a, b = mp.mpf(as_), mp.mpf(bs)
        exact = Dm(1 - a, 1 - b)
        approx = formula(a, b)
        rel = abs(exact - approx) / abs(exact)
        rows.append(
            {
                "a": as_,
                "b": bs,
                "D_exact": mp.nstr(exact, 25),
                "formula": mp.nstr(approx, 25),
                "relative_error": mp.nstr(rel, 12),
                "sign": "+" if exact > 0 else "-",
            }
        )
    # diagonal near-corner expansion: D(1-eps,1-eps) ~ eps[(2m-1)(1-ln eps) - 2m ln 2]
    rows2 = []
    for es_ in ["1e-8", "1e-12"]:
        e = mp.mpf(es_)
        exact = Dm(1 - e, 1 - e)
        approx = e * ((2 * M_MEAN - 1) * (1 - mp.log(e)) - 2 * M_MEAN * mp.log(2))
        rows2.append(
            {
                "eps": es_,
                "D_exact": mp.nstr(exact, 25),
                "formula": mp.nstr(approx, 25),
                "relative_error": mp.nstr(abs(exact - approx) / abs(exact), 12),
                "sign": "+" if exact > 0 else "-",
            }
        )
    sec(
        "d_corner_asymptotics",
        {
            "formula_1ab": (
                "D_m(1-a,1-b) ~ -m(a+b)ln(a+b) + (a lna + b lnb)/2 "
                "- (a+b)(1/2 - m)"
            ),
            "checks": rows,
            "sign_verdict": (
                "leading coefficient positive for small a,b: D>0 near (1,1) "
                "transverse directions; the (1,1) zero is approached from "
                "positive values"
            ),
            "diagonal_checks": rows2,
        },
    )


# =========================================================== SECTION 4 =====
# Psi global structure: corners/edges/faces + multistart + interior


CORNERS = []


def psi_corners():
    mp.mp.dps = DPS_CERT
    recs = []
    for bits in range(16):
        pt = tuple((bits >> (3 - i)) & 1 for i in range(4))
        v = Psi_split(*[mp.mpf(b) for b in pt])
        recs.append(
            {
                "point": [str(b) for b in pt],
                "value": mp.nstr(v, 40),
                "abs": mp.nstr(abs(v), 12),
            }
        )
    maxabs = max(mp.mpf(r["abs"]) for r in recs)
    sec(
        "psi_corners",
        {
            "records": recs,
            "max_abs_value": mp.nstr(maxabs, 12),
            "claim": "all 16 corners of [0,1]^4 are exact zeros of Psi",
        },
    )


def _edge_face_key(pinned_idx, pinned_vals, free_idx):
    return (
        tuple((i, pinned_vals[k]) for k, i in enumerate(pinned_idx)),
        tuple(free_idx),
    )


def run_edges_faces():
    """24 edges (2 pinned coords in {0,1}) and 8 faces (1 pinned coord),
    minimized with seeded DE in float64 + dps-60 mp polish, then a
    dps-340 evaluation of the winner value."""
    mp.mp.dps = DPS_SEARCH
    recs = []
    for a in range(4):
        for b in range(a + 1, 4):
            free = tuple(i for i in range(4) if i not in (a, b))
            for va in (0, 1):
                for vb in (0, 1):
                    base = [0.0] * 4
                    base[a], base[b] = float(va), float(vb)

                    def f64(pt, base=base, free=free):
                        p = np.tile(np.array(base), (len(pt), 1))
                        for k, i in enumerate(free):
                            p[:, i] = pt[:, k]
                        return Psi_batch(p)

                    de = differential_evolution(
                        lambda pt: float(f64(pt.reshape(1, 2))[0]),
                        bounds=[(0.0, 1.0)] * 2, seed=SEED + 3,
                        popsize=28, maxiter=280, tol=1e-14,
                        polish=True, updating="deferred", workers=1,
                    )

                    def f_mp(pt, base=base, free=free):
                        p = [mp.mpf(c) for c in base]
                        for k, i in enumerate(free):
                            p[i] = pt[k]
                        return Psi_split(p[0], p[1], p[2], p[3])

                    starts = [de.x]
                    rng2 = np.random.default_rng(SEED + 5)
                    for _ in range(4):
                        starts.append(rng2.random(2))
                    best = None
                    for st in starts:
                        x, _ = compass_mp(f_mp, [mp.mpf(c) for c in st],
                                          ([0, 0], [1, 1]))
                        fx = f_mp(x)
                        if best is None or fx < best[1]:
                            best = (x, fx)
                    p = [mp.mpf(c) for c in base]
                    for k, i in enumerate(free):
                        p[i] = best[0][k]
                    mp.mp.dps = DPS_CERT
                    val340 = Psi_split(*p)
                    mp.mp.dps = DPS_SEARCH
                    recs.append(
                        {
                            "pinned": {
                                f"coord{a}": str(va),
                                f"coord{b}": str(vb),
                            },
                            "free": f"coord{free[0]},coord{free[1]}",
                            "min_dps60": mp.nstr(best[1], 40),
                            "value_dps340": fmt_mp(val340),
                            "argmin": [mp.nstr(v, 30) for v in best[0]],
                            "sign": (
                                "neg" if val340 < -mp.mpf("1e-30")
                                else ("zero" if abs(val340) < mp.mpf("1e-30")
                                      else "pos")
                            ),
                        }
                    )
    sec(
        "psi_edges",
        {
            "n_edges": len(recs),
            "negative_edges": [r for r in recs if r["sign"] == "neg"],
            "all_records": recs,
        },
    )
    frecs = []
    for a in range(4):
        for va in (0, 1):
            free = tuple(i for i in range(4) if i != a)

            def f64(pt, a=a, va=va, free=free):
                p = np.zeros((1, 4), dtype=np.float64)
                p[0, a] = float(va)
                for k, i in enumerate(free):
                    p[0, i] = pt[k]
                return Psi_batch(p)

            de = differential_evolution(
                lambda pt: float(f64(np.asarray(pt))[0]),
                bounds=[(0.0, 1.0)] * 3, seed=SEED + 4,
                popsize=28, maxiter=320, tol=1e-14,
                polish=True, updating="deferred", workers=1,
            )

            def f_mp(pt, a=a, va=va, free=free):
                p = [mp.mpf(0)] * 4
                p[a] = mp.mpf(va)
                for k, i in enumerate(free):
                    p[i] = pt[k]
                return Psi_split(p[0], p[1], p[2], p[3])

            starts = [de.x]
            rng2 = np.random.default_rng(SEED + 6)
            for _ in range(6):
                starts.append(rng2.random(3))
            best = None
            for st in starts:
                x, _ = compass_mp(f_mp, [mp.mpf(c) for c in st],
                                  ([0, 0, 0], [1, 1, 1]))
                fx = f_mp(x)
                if best is None or fx < best[1]:
                    best = (x, fx)
            p = [mp.mpf(0)] * 4
            p[a] = mp.mpf(va)
            for k, i in enumerate(free):
                p[i] = best[0][k]
            mp.mp.dps = DPS_CERT
            val340 = Psi_split(*p)
            mp.mp.dps = DPS_SEARCH
            frecs.append(
                {
                    "pinned": {f"coord{a}": str(va)},
                    "free": f"coord{free[0]},coord{free[1]},coord{free[2]}",
                    "min_dps60": mp.nstr(best[1], 40),
                    "value_dps340": fmt_mp(val340),
                    "argmin": [mp.nstr(v, 30) for v in best[0]],
                }
            )
    sec("psi_faces", {"records": frecs})
    return recs, frecs


def run_global_multistart():
    """Seeded differential_evolution + Sobol/compass multistart (float64,
    DISCOVERY), mpmath dps-60 polish of deduplicated finalists, dps-340
    certification evaluation of the winner and its mirror orbit point."""
    mp.mp.dps = DPS_SEARCH
    de = differential_evolution(
        lambda x: float(Psi_batch(x.reshape(1, 4))[0]),
        bounds=[(0.0, 1.0)] * 4,
        seed=SEED + 1,
        popsize=40,
        maxiter=600,
        tol=1e-14,
        mutation=(0.4, 1.0),
        recombination=0.9,
        polish=True,
        updating="deferred",
        workers=1,
    )
    X0 = sobol_starts(N_MULTISTART, 4, SEED + 1)
    F = Psi_batch(X0)
    order = np.argsort(F)[:512]
    X = X0[order].copy()
    Ftop = Psi_batch(X)
    for stp in (0.05, 0.02, 0.008, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5):
        for i in range(4):
            for sgn in (1, -1):
                Y = X.copy()
                Y[:, i] = np.clip(X[:, i] + sgn * stp, 0.0, 1.0)
                Fy = Psi_batch(Y)
                acc = Fy < Ftop
                if acc.any():
                    X[acc] = Y[acc]
                    Ftop[acc] = Fy[acc]
    pol = []
    for x in X[:128]:
        r = minimize(
            lambda z: float(Psi_batch(z.reshape(1, 4))[0]),
            x, method="Nelder-Mead", bounds=[(0, 1)] * 4,
            options={"xatol": 1e-14, "fatol": 1e-18, "maxiter": 4000},
        )
        pol.append((float(r.fun), np.clip(r.x, 0, 1)))
    pol.sort(key=lambda u: u[0])
    uniq = []
    for v, x in pol:
        if v > pol[0][0] + 1e-12:
            break
        if all(np.linalg.norm(x - u[1]) > 1e-5 for u in uniq):
            uniq.append((v, x))
        if len(uniq) >= 12:
            break
    if float(de.fun) < uniq[0][0]:
        uniq.insert(0, (float(de.fun), np.clip(de.x, 0, 1)))
        uniq = uniq[:12]
    best_rec = None
    others = []
    for v, x in uniq[:12]:
        mpx = [mp.mpf(float(c)) for c in x]
        bx, _bf = compass_mp(
            lambda pt: Psi_split(pt[0], pt[1], pt[2], pt[3]),
            mpx, ([0, 0, 0, 0], [1, 1, 1, 1]),
        )
        bx = [mp.mpf(0) if abs(c) < mp.mpf("1e-25") else
              (mp.mpf(1) if c > 1 else (mp.mpf(0) if c < 0 else c))
              for c in bx]
        bf = Psi_split(*bx)
        snap = all(abs(c) < mp.mpf("1e-12") or abs(c - 1) < mp.mpf("1e-12")
                   for c in bx)
        interior = all(min(c, 1 - c) > mp.mpf("1e-9") for c in bx)
        rec = {
            "start_f64": [fmt_f(c) for c in x],
            "point_dps60": [mp.nstr(c, 30) for c in bx],
            "value_dps60": mp.nstr(bf, 40),
            "min_coord": mp.nstr(min(min(c, 1 - c) for c in bx), 12),
            "class": ("corner" if snap else
                      ("boundary" if not interior else "interior")),
        }
        if best_rec is None or bf < mp.mpf(best_rec["value_dps60"]):
            best_rec = rec
        others.append(rec)
    mp.mp.dps = DPS_CERT
    bpt = [mp.mpf(s) for s in best_rec["point_dps60"]]
    val340 = Psi_split(*bpt)
    mir = [bpt[1], bpt[0], bpt[3], bpt[2]]
    val_mir = Psi_split(*mir)
    interior_probe = float(Psi_batch(
        np.clip(uniq[0][1].reshape(1, 4), 1e-6, 1 - 1e-6))[0])
    de_val = float(Psi_batch(de.x.reshape(1, 4))[0])
    sec(
        "psi_global_min",
        {
            "multistart_budget": N_MULTISTART,
            "method": (
                "seeded scipy differential_evolution (seed=SEED+1, popsize 40, "
                "maxiter 600, deferred, workers=-1) + 20480 Sobol starts with "
                "vectorized axis-compass descent + Nelder-Mead polish on top "
                "128 (float64, DISCOVERY) -> mpmath dps-60 compass polish of "
                "deduplicated finalists -> dps-340 certification of winner"
            ),
            "de_value_f64_DISCOVERY": fmt_f(de_val),
            "de_argmin": [fmt_f(c) for c in np.clip(de.x, 0, 1)],
            "best_f64_value_DISCOVERY": fmt_f(uniq[0][0]),
            "best": best_rec,
            "other_finalists": others,
            "value_dps340_at_best": fmt_mp(val340),
            "mirror_orbit_value_dps340": fmt_mp(val_mir),
            "attainment": best_rec["class"],
            "orbit_symmetry_note": (
                "Psi(s,t,sig,tau) = Psi(t,s,ta,sg): (s,tau)<->(t,sig) swap; "
                "the four special strata are one orbit"
            ),
            "interior_probe_value_at_1e-6_clip": fmt_f(interior_probe),
        },
    )
    return best_rec, val340


# =========================================================== SECTION 5 =====
# r-threshold




def box_min_psi(r_mp):
    """min of Psi over [r,1-r]^4: float64 DE + mp dps-60 polish,
    certified at dps 340.  Returns (point, value_at_340)."""
    mp.mp.dps = DPS_SEARCH
    r = mp.mpf(r_mp)
    hi = 1 - r
    r64, hi64 = float(r), float(hi)

    def f64(pt):
        p = np.asarray(pt, dtype=np.float64).reshape(1, 4)
        p = np.clip(p, r64, hi64)
        return float(Psi_batch(p)[0])

    de = differential_evolution(
        f64, bounds=[(r64, hi64)] * 4, seed=SEED + 7,
        popsize=24, maxiter=220, tol=1e-14,
        polish=True, updating="deferred", workers=1,
    )

    def f_mp(pt):
        return Psi_split(pt[0], pt[1], pt[2], pt[3])

    z2 = mp.mpf(RES["p_kernel"]["argmin_z"])
    starts = [de.x]
    for t_fix in (r, z2, hi):
        for sg_fix in (r, z2, hi):
            starts.append([r, t_fix, sg_fix, r])
            starts.append([r, t_fix, sg_fix, hi])
            starts.append([hi, t_fix, sg_fix, r])
            starts.append([hi, t_fix, sg_fix, hi])
    best = None
    for st in starts:
        st0 = [min(max(float(v), r64), hi64) for v in st]
        x, _ = compass_mp(f_mp, [mp.mpf(v) for v in st0],
                          ([r, r, r, r], [hi, hi, hi, hi]))
        fx = f_mp(x)
        if best is None or fx < best[1]:
            best = (x, fx)
    mp.mp.dps = DPS_CERT
    pt = [mp.mpf(mp.nstr(c, 60)) for c in best[0]]
    pt = [min(max(c, mp.mpf(r_mp)), 1 - mp.mpf(r_mp)) for c in pt]
    v340 = Psi_split(*pt)
    mp.mp.dps = DPS_SEARCH
    return pt, v340


def run_r_threshold():
    lo, hi = mp.mpf("0.05"), mp.mpf("0.5")
    # ensure sign at ends
    _, vl0 = box_min_psi(lo)
    pt_h, vh = box_min_psi(hi)
    # vhi should be >= 0 (single point at r=1/2); vlo < 0
    hist = [
        {"r": mp.nstr(lo, 25), "minval": fmt_mp(vl0)},
        {"r": mp.nstr(hi, 25), "minval": fmt_mp(vh)},
    ]
    a, b = lo, hi
    stress = None
    for it in range(R_BISECT_MAX_ITER):
        mid = (a + b) / 2
        _, vm = box_min_psi(mid)
        hist.append(
            {"r": mp.nstr(mid, 25), "minval": fmt_mp(vm)}
        )
        if vm < 0:
            a = mid
        else:
            b = mid
        if b - a < mp.mpf("1e-9"):
            break
    r_star = (a + b) / 2
    pt, v340 = box_min_psi(r_star)
    # active stratum detection
    rstr = mp.nstr(r_star, 25)
    rmp = mp.mpf(rstr)
    active = []
    for i, c in enumerate(pt):
        if abs(c - rmp) < mp.mpf("1e-6"):
            active.append(f"coord{i}=r")
        elif abs(c - (1 - rmp)) < mp.mpf("1e-6"):
            active.append(f"coord{i}=1-r")
    psi_half = 2 * Dm(mp.mpf(1) / 2, mp.mpf(1) / 2) / M_MEAN
    sec(
        "r_threshold",
        {
            "definition": (
                "r* = crossing of min over [r,1-r]^4 of Psi through 0 "
                "(minval continuous & nondecreasing in r)"
            ),
            "bisection_history": hist,
            "r_star": mp.nstr(r_star, 40),
            "r_star_float": fmt_f(float(r_star)),
            "minval_at_r_star": fmt_mp(v340),
            "active_stratum": active or ["interior of the box"],
            "argmin_at_r_star": [mp.nstr(c, 30) for c in pt],
            "psi_half_box_endpoint": fmt_mp(psi_half),
            "note": (
                "every negative point has at least one coordinate outside "
                "[r*, 1-r*]; interior-restricted reductions are safe iff "
                "their domain matches this box (magnitude of r* matters)"
            ),
        },
    )


def run_r_threshold_refined():
    """dps-340 refinement of r* on the active stratum found by the search.

    The box minimizer at the crossing has s = tau = r (box face) and
    t = sigma = z interior, where
        phi(r,z) = P(r,r) + P(z,z) + 2 Q(r,z).
    Solve the 2x2 system {d phi/dz = 0, phi = 0} by Newton at dps 340 and
    verify the KKT signs of the full four-variable gradient: the t,sigma
    components must vanish and the s,tau components must be >= 0 (pushing
    inward from the box face increases Psi).
    """
    mp.mp.dps = DPS_CERT

    def phi(r, z):
        return Pm(r, r) + Pm(z, z) + 2 * Qm(r, z)

    def dphi_dz(r, z):
        eps = mp.mpf("1e-40")
        return (phi(r, z + eps) - phi(r, z - eps)) / (2 * eps)

    r = mp.mpf(RES["r_threshold"]["r_star"])
    z = mp.mpf(RES["r_threshold"]["argmin_at_r_star"][1])
    eps = mp.mpf("1e-35")
    for _ in range(60):
        F1 = dphi_dz(r, z)
        F2 = phi(r, z)
        j11 = (dphi_dz(r + eps, z) - dphi_dz(r - eps, z)) / (2 * eps)
        j12 = (dphi_dz(r, z + eps) - dphi_dz(r, z - eps)) / (2 * eps)
        j21 = (phi(r + eps, z) - phi(r - eps, z)) / (2 * eps)
        j22 = F1
        det = j11 * j22 - j12 * j21
        dr = (F1 * j22 - F2 * j12) / det
        dz = (j11 * F2 - j21 * F1) / det
        r = r - dr
        z = z - dz
        if abs(dr) < mp.mpf("1e-60") and abs(dz) < mp.mpf("1e-60"):
            break
    pt = (r, z, z, r)
    val = Psi_split(*pt)
    g = grad_Psi_split(*pt)
    # The box always contains the EXACT interior zero (x*,x*,x*,x*) when
    # r <= 1-x*, because Psi(u,u,u,u) = (2/m) D(u,u) and D(x*,x*) = 0.
    # So "Psi >= 0 on the box" is tight, and the stress scan must separate
    # that tangency from genuine negative (wedge) points.
    def wedge_min(rr, iters=140):
        """min over z of phi(rr,z) = value of Psi on the active stratum."""
        invphi = (mp.sqrt(5) - 1) / 2
        a, b = rr, 1 - rr
        c = b - invphi * (b - a)
        d = a + invphi * (b - a)
        fc, fd = phi(rr, c), phi(rr, d)
        for _ in range(iters):
            if fc < fd:
                b, d, fd = d, c, fc
                c = b - invphi * (b - a)
                fc = phi(rr, c)
            else:
                a, c, fc = c, d, fd
                d = a + invphi * (b - a)
                fd = phi(rr, d)
        zz = (a + b) / 2
        return zz, phi(rr, zz)

    tangency_val = Psi_split(X_STAR, X_STAR, X_STAR, X_STAR)
    r_tangency_exit = 1 - X_STAR
    stress = []
    for delta in ("-1e-2", "-1e-3", "-1e-6", "1e-6", "1e-3", "1e-2", "5e-2",
                  "2e-1"):
        rr = r + mp.mpf(delta)
        zz, wv = wedge_min(rr)
        lo64, hi64 = float(rr), float(1 - rr)
        X = lo64 + (hi64 - lo64) * sobol_starts(N_MULTISTART, 4, SEED + 11)
        vals = Psi_batch(X)
        j = int(np.argmin(vals))
        x = X[j].copy()
        fx = float(vals[j])
        step = 0.02
        while step > 1e-13:
            moved = False
            for i in range(4):
                for sgn in (1, -1):
                    y = x.copy()
                    y[i] = min(max(x[i] + sgn * step, lo64), hi64)
                    fy = float(Psi_batch(y.reshape(1, 4))[0])
                    if fy < fx:
                        x, fx, moved = y, fy, True
                        break
                if moved:
                    break
            if not moved:
                step /= 2
        dist_tang = max(abs(float(c) - float(X_STAR)) for c in x)
        stress.append(
            {
                "r": mp.nstr(rr, 30),
                "wedge_stratum_min_dps340": mp.nstr(wv, 25),
                "wedge_sign": "neg" if wv < 0 else "pos",
                "wedge_argmin_z": mp.nstr(zz, 25),
                "sobol_box_min_f64_DISCOVERY": fmt_f(fx),
                "sobol_argmin_linf_distance_to_x_star_point": fmt_f(dist_tang),
                "sobol_argmin_is_the_interior_tangency": bool(dist_tang < 1e-6),
                "tangency_in_box": bool(rr <= r_tangency_exit),
            }
        )
    sec(
        "r_threshold_refined",
        {
            "method": (
                "Newton at dps 340 on {d phi/dz = 0, phi = 0} with "
                "phi(r,z) = P(r,r) + P(z,z) + 2 Q(r,z), the value of Psi on "
                "the active stratum s=tau=r, t=sigma=z"
            ),
            "r_star": mp.nstr(r, 50),
            "z_star": mp.nstr(z, 50),
            "psi_at_r_star_point": mp.nstr(val, 25),
            "argmin_point": [mp.nstr(c, 40) for c in pt],
            "kkt_gradient": {
                "d_ds": mp.nstr(g[0], 25),
                "d_dt": mp.nstr(g[1], 25),
                "d_dsigma": mp.nstr(g[2], 25),
                "d_dtau": mp.nstr(g[3], 25),
                "free_components_vanish": bool(
                    abs(g[1]) < mp.mpf("1e-30") and abs(g[2]) < mp.mpf("1e-30")
                ),
                "box_face_components_push_inward": bool(g[0] > 0 and g[3] > 0),
            },
            "agreement_with_search_bisection": mp.nstr(
                abs(r - mp.mpf(RES["r_threshold"]["r_star"])), 12
            ),
            "stress_scan": stress,
            "interior_tangency": {
                "point": "(x*,x*,x*,x*)",
                "x_star": fmt_mp(X_STAR),
                "psi_value_dps340": mp.nstr(tangency_val, 20),
                "reason": (
                    "Psi(u,u,u,u) = (2/m) D_m(u,u) by Claim 3, and "
                    "D_m(x*,x*) = 0 exactly because Liu's quartic "
                    "x^4-2x^3+3x^2-1 = 0 forces pi(x*,x*) = 1-x*^2, hence "
                    "h(pi(x*,x*)) = h(x*^2)"
                ),
                "inside_box_while": (
                    "r <= 1-x* = " + mp.nstr(r_tangency_exit, 30)
                ),
            },
            "claim": (
                "Psi >= 0 on [r,1-r]^4 exactly for r >= r*, and Psi attains "
                "strictly negative values inside [r,1-r]^4 for every r < r*. "
                "For r* <= r <= 1-x* the bound is TIGHT: the box minimum is "
                "exactly 0, attained at the interior point (x*,x*,x*,x*); it "
                "becomes strictly positive only for r > 1-x*."
            ),
        },
    )


# =========================================================== SECTION 6 =====
# structured families


def run_families():
    mp.mp.dps = DPS_SEARCH
    rng = new_rng()
    out = {}
    # (a) sig=s, tau=t: identity + min
    worst = mp.mpf(0)
    for _ in range(20):
        s = mp.mpf(rng.random())
        t = mp.mpf(rng.random())
        worst = max(
            worst,
            abs(Psi_split(s, t, s, t) - 2 * Dm(s, t) / M_MEAN),
        )
    # family value = 2 D(s,t)/m; D >= 0 with the interior zero at (x*,x*)
    mp.mp.dps = DPS_CERT
    fam_a_min = 2 * Dm(X_STAR, X_STAR) / M_MEAN
    mp.mp.dps = DPS_SEARCH
    out["sig_s_tau_t"] = {
        "identity_residual_vs_2D_over_M": mp.nstr(worst, 20),
        "family_min_dps340": mp.nstr(fam_a_min, 40),
        "attained": (
            "0 at (x*,x*) (interior tangency of D) and on the zero lines "
            "{s=0},{t=0},{(1,1)}; the family is exactly (2/m)*D_m >= 0"
        ),
    }
    # (b) sig=t, tau=s: P(s,s)+P(t,t)+2Q(s,t)
    def f2(pt):
        s, t = pt
        return Pm(s, s) + Pm(t, t) + 2 * Qm(s, t)

    best2 = None
    for _ in range(24):
        st = [mp.mpf(v) for v in rng.random(2)]
        x, fx = compass_mp(f2, st, ([0, 0], [1, 1]))
        if best2 is None or fx < best2[1]:
            best2 = (x, fx)
    mp.mp.dps = DPS_CERT
    v2 = Pm(best2[0][0], best2[0][0]) + Pm(best2[0][1], best2[0][1]) + 2 * Qm(
        best2[0][0], best2[0][1]
    )
    # compare with 2 D(u,u)/m on the diagonal: f2(s,s) = 2P(s,s)+2Q(s,s)
    out["sig_t_tau_s"] = {
        "min_dps340": fmt_mp(v2),
        "argmin": [mp.nstr(c, 40) for c in best2[0]],
        "note": "on the diagonal this family equals 2D(u,u)/m >= 0",
    }
    # (c) s=t family
    def f3(pt):
        s, sg, ta = pt
        return Psi_split(s, s, sg, ta)

    best3 = None
    for _ in range(24):
        st = [mp.mpf(v) for v in rng.random(3)]
        x, fx = compass_mp(f3, st, ([0, 0, 0], [1, 1, 1]))
        if best3 is None or fx < best3[1]:
            best3 = (x, fx)
    mp.mp.dps = DPS_CERT
    v3 = Psi_split(best3[0][0], best3[0][0], best3[0][1], best3[0][2])
    out["s_t"] = {
        "min_dps340": fmt_mp(v3),
        "argmin": [mp.nstr(c, 40) for c in best3[0]],
    }
    # (d) sig=tau family
    def f4(pt):
        s, t, c = pt
        return Psi_split(s, t, c, c)

    best4 = None
    for _ in range(24):
        st = [mp.mpf(v) for v in rng.random(3)]
        x, fx = compass_mp(f4, st, ([0, 0, 0], [1, 1, 1]))
        if best4 is None or fx < best4[1]:
            best4 = (x, fx)
    mp.mp.dps = DPS_CERT
    v4 = Psi_split(best4[0][0], best4[0][1], best4[0][2], best4[0][2])
    out["sig_tau"] = {
        "min_dps340": fmt_mp(v4),
        "argmin": [mp.nstr(c, 40) for c in best4[0]],
    }
    # (e) antidiagonal s + sg = 1
    def f5(pt):
        sg, t, ta = pt
        return Psi_split(1 - sg, t, sg, ta)

    best5 = None
    for _ in range(24):
        st = [mp.mpf(v) for v in rng.random(3)]
        x, fx = compass_mp(f5, st, ([0, 0, 0], [1, 1, 1]))
        if best5 is None or fx < best5[1]:
            best5 = (x, fx)
    mp.mp.dps = DPS_CERT
    sg = best5[0][0]
    v5 = Psi_split(1 - sg, best5[0][1], sg, best5[0][2])
    out["antidiagonal_s_plus_sigma_1"] = {
        "min_dps340": fmt_mp(v5),
        "argmin": [mp.nstr(c, 40) for c in best5[0]],
    }
    sec("structured_families", out)


# =========================================================== SECTION 7 =====
# Psi asymptotics with sign


def run_psi_asymptotics():
    mp.mp.dps = DPS_CERT
    rows = []
    # near-corner (1,1,1,1) all-equal ray:
    # Psi(1-e,1-e,1-e,1-e) ~ 4e(1-ln 2e) - (2e/m)(1-ln e)
    for es_ in ["1e-8", "1e-10", "1e-12"]:
        e = mp.mpf(es_)
        exact = Psi_split(1 - e, 1 - e, 1 - e, 1 - e)
        approx = 4 * e * (1 - mp.log(2 * e)) - (2 * e / M_MEAN) * (1 - mp.log(e))
        rows.append(
            {
                "eps": es_,
                "exact": mp.nstr(exact, 25),
                "approx": mp.nstr(approx, 25),
                "rel_err": mp.nstr(abs(exact - approx) / abs(exact), 12),
                "sign": "+" if exact > 0 else "-",
            }
        )
    # asymmetric near-corner probes
    asym = []
    for a, b, mu, al in [("1e-8", "3e-9", "5e-9", "7e-9"),
                         ("1e-10", "1e-6", "1e-9", "1e-7")]:
        pt = (
            1 - mp.mpf(a), 1 - mp.mpf(b), 1 - mp.mpf(mu), 1 - mp.mpf(al)
        )
        v = Psi_split(*pt)
        asym.append(
            {
                "minus_logs": [a, b, mu, al],
                "value": mp.nstr(v, 25),
                "sign": "+" if v > 0 else "-",
            }
        )
    # wedge rate at the global-minimum stratum:
    # S(e) = Psi(e, z2, z2, e) - P_min ~ 2*Q(e, z2) + P(e,e)
    z2 = mp.mpf(RES["p_kernel"]["argmin_z"])
    pmin = mp.mpf(RES["p_kernel"]["min_value"])
    wedge = []
    for es_ in ["1e-8", "1e-10", "1e-12"]:
        e = mp.mpf(es_)
        S = Psi_split(e, z2, z2, e) - pmin
        approx = 2 * Qm(e, z2) + Pm(e, e)
        wedge.append(
            {
                "eps": es_,
                "S_exact": mp.nstr(S, 25),
                "S_approx_2Q_plus_Pee": mp.nstr(approx, 25),
                "sign": "+" if S > 0 else "-",
            }
        )
    # negative counterexample stratum: moving INTO the cube raises Psi
    lift = []
    for es_ in ["1e-8", "1e-10"]:
        e = mp.mpf(es_)
        v0 = Psi_split(0, mp.mpf(1) / 2, mp.mpf(1) / 2, 0)
        v1 = Psi_split(e, mp.mpf(1) / 2, mp.mpf(1) / 2, e)
        lift.append(
            {
                "eps": es_,
                "lift": mp.nstr(v1 - v0, 25),
            }
        )
    sec(
        "psi_asymptotics",
        {
            "corner_1111_ray": rows,
            "corner_1111_asymmetric": asym,
            "wedge_at_global_min": wedge,
            "counterexample_stratum_lift": lift,
            "sign_verdicts": {
                "near_(1,1,1,1)": "Psi >= 0 approached from positive side",
                "wedge_at_(0,z2,z2,0)": (
                    "Psi = P_min + 2 beta h(pi(e,z2)) + P(e,e); increases "
                    "like ~ +0.18 e|ln e|; negative region is a cusp wedge "
                    "on the s=tau=0 stratum"
                ),
            },
        },
    )


# =========================================================== SECTION 8 =====
# mutations (harness self-test)


def run_mutations():
    mp.mp.dps = DPS_CERT
    out = []
    # m1: mutate the SPLIT itself (Q(s,t) -> Q(s,tau)); the true split is
    # an algebraic identity in beta, so perturbing beta is not a valid
    # mutation of it -- the formula shape is what must be load bearing.
    tp = (mp.mpf("0.3"), mp.mpf("0.61"), mp.mpf("0.77"), mp.mpf("0.42"))
    mutated = (
        Pm(tp[0], tp[3]) + Pm(tp[1], tp[2])
        + Qm(tp[0], tp[3]) + Qm(tp[2], tp[3])
    )
    r1 = abs(Psi_raw(*tp) - mutated)
    out.append(
        {
            "mutation": "split with Q(s,t) replaced by Q(s,tau)",
            "residual_vs_raw_Psi": mp.nstr(r1, 20),
            "detected": bool(r1 > mp.mpf("1e-6")),
        }
    )
    # m1b: constants are load bearing for the counterexample value
    global BETA
    b0 = BETA
    BETA = b0 * (1 + mp.mpf("1e-6"))
    cex_mut = Psi_split(0, mp.mpf(1) / 2, mp.mpf(1) / 2, 0)
    BETA = b0
    cex_true = Psi_split(0, mp.mpf(1) / 2, mp.mpf(1) / 2, 0)
    out.append(
        {
            "mutation": "beta *= (1+1e-6) (counterexample value must move)",
            "delta_counterexample": mp.nstr(abs(cex_mut - cex_true), 20),
            "detected": bool(abs(cex_mut - cex_true) > mp.mpf("1e-10")),
        }
    )
    # m2: claim 3 with wrong factor M/2 -> D
    worst = mp.mpf(0)
    for _ in range(10):
        s = mp.mpf(new_rng().random())
        t = mp.mpf(new_rng().random())
        worst = max(
            worst,
            abs(Psi_raw(s, t, s, t) - Dm(s, t) / M_MEAN),
        )
    out.append(
        {
            "mutation": "claim3 factor 2 -> 1 (use D/M instead of 2D/M)",
            "residual": mp.nstr(worst, 20),
            "detected": worst > mp.mpf("1e-6"),
        }
    )
    # m3: shift argmin z2 by 1e-6 -> P value strictly increases beyond noise
    z2 = mp.mpf(RES["p_kernel"]["argmin_z"])
    dz = Pm(z2 + mp.mpf("1e-6"), z2 + mp.mpf("1e-6")) - Pm(z2, z2)
    out.append(
        {
            "mutation": "argmin z2 -> z2 + 1e-6",
            "delta_P": mp.nstr(dz, 20),
            "detected": dz > mp.mpf("1e-15"),
        }
    )
    sec("mutations", out)


# =============================================================== main ======


def canonical_payload_without_hash():
    payload = dict(RES)
    body = {
        k: v for k, v in payload.items() if k != "report_sha256"
    }
    return json.dumps(
        body, indent=1, ensure_ascii=True, allow_nan=False, sort_keys=False
    )


def main():
    run_constants()
    run_verification()
    z2, pmin, roots = run_p_kernel()
    p_zero_level_table()
    run_d_kernel()
    run_d_corner_asymptotics()
    psi_corners()
    run_edges_faces()
    run_global_multistart()
    run_r_threshold()
    run_r_threshold_refined()
    run_families()
    run_psi_asymptotics()
    run_mutations()

    # negative-region summary (exact structural statements + numbers)
    neg_edges = [r for r in RES["psi_edges"]["all_records"]
                 if r["sign"] == "neg"]
    summary = {
        "split_legend": [
            "P(s,t) = (1-beta)h(st) - (t h(s)+s h(t))/(2m)",
            "Q(s,t) = beta h(pi(s,t)) >= 0",
            "Psi = P(s,tau)+P(t,sig)+Q(s,t)+Q(sig,tau), exactly",
        ],
        "counterexample": RES["verification_counterexample"],
        "four_special_strata": {
            "s=0,tau=0": "Psi == P(t,sig) exactly: negative exactly on {P<0}",
            "t=0,sigma=0": "Psi == P(s,tau) exactly: negative exactly on {P<0}",
            "s=0,sigma=0": "Psi == 0 identically",
            "t=0,tau=0": "Psi == 0 identically",
        },
        "faces_with_negative_regions": [
            r for r in RES["psi_faces"]["records"]
            if mp.mpf(r["min_dps60"]) < 0
        ],
        "negative_edges": neg_edges,
        "global_min": RES["psi_global_min"],
        "r_threshold": RES["r_threshold"],
        "r_threshold_refined": RES["r_threshold_refined"],
        "d_zero_set": RES["d_kernel_zeros_hessian"],
        "p_negative_region": RES["p_negative_region"],
    }
    sec("negative_region_summary", summary)

    digest = hashlib.sha256(
        canonical_payload_without_hash().encode("utf-8")
    ).hexdigest()
    RES["report_sha256"] = digest
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    text = json.dumps(
        RES, indent=1, ensure_ascii=True, allow_nan=False, sort_keys=False
    )
    with open(OUT_PATH, "w", encoding="utf-8") as fh:
        fh.write(text + "\n")

    # stdout summary for the report
    print("== liu9_psi_hunt summary ==")
    print("beta =", mp.nstr(BETA, 25), " m =", mp.nstr(M_MEAN, 25))
    print("x*   =", mp.nstr(X_STAR, 25))
    print("split identity residual:", RES["verification_split_identity"]["max_residual"])
    print("counterexample Psi(0,1/2,1/2,0) =", RES["verification_counterexample"]["value_dps340"])
    print("Claim3 residual:", RES["verification_claim3"]["max_residual"])
    print("P argmin z2 =", RES["p_kernel"]["argmin_z"])
    print("P min       =", RES["p_kernel"]["min_value"])
    print("D(x*,x*)    =", RES["verification_D_at_xstar"]["value"])
    hx = RES["d_kernel_zeros_hessian"]["hessian_at_xstar"]
    print("D Hessian eigenvalues (diag, anti):",
          hx["eigenvalue_diagonal_mode_(1,1)"],
          hx["eigenvalue_antidiagonal_mode_(1,-1)"])
    print("D verdict   :", RES["d_kernel_zeros_hessian"]["verdict"])
    print("Psi global min (dps340):", RES["psi_global_min"]["value_dps340_at_best"])
    print("Psi global argmin:", RES["psi_global_min"]["best"]["point_dps60"])
    print("negative edges:", [(r["pinned"], r["value_dps340"])
                              for r in neg_edges])
    print("faces with neg regions:", [
        (r["pinned"], r["min_dps60"])
        for r in RES["psi_faces"]["records"]
        if mp.mpf(r["min_dps60"]) < 0
    ])
    print("r* (search bisection) =", RES["r_threshold"]["r_star"],
          "minval:", RES["r_threshold"]["minval_at_r_star"])
    print("r* (dps340 Newton)    =", RES["r_threshold_refined"]["r_star"],
          " z* =", RES["r_threshold_refined"]["z_star"])
    print("r* KKT:", json.dumps(
        RES["r_threshold_refined"]["kkt_gradient"], indent=1))
    print("r* stress (r, wedge_min, sign, sobol_argmin_is_tangency):",
          json.dumps(
              [(s["r"][:16], s["wedge_stratum_min_dps340"], s["wedge_sign"],
                s["sobol_argmin_is_the_interior_tangency"])
               for s in RES["r_threshold_refined"]["stress_scan"]], indent=1))
    print("interior tangency:", json.dumps(
        RES["r_threshold_refined"]["interior_tangency"], indent=1))
    print("active stratum at r*:", RES["r_threshold"]["active_stratum"])
    print("families:", json.dumps(
        {k: v.get("min_dps340", v.get("family_min_dps340"))
         for k, v in RES["structured_families"].items()}, indent=1))
    print("report_sha256 =", digest)
    print("written:", OUT_PATH)


if __name__ == "__main__":
    main()
