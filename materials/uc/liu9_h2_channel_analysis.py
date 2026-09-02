#!/usr/bin/env python3
"""Deterministic analytic/discovery report for the Liu-H2 channel routes.

The program writes nothing: its only output is byte-stable JSON on stdout.
Arb-labelled statements are enclosed independently at 320 and 512 bits;
all NumPy/mpmath searches are explicitly labelled DISCOVERY.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

for _key in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS",
):
    os.environ[_key] = "1"

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import mpmath as mp
import numpy as np
import sympy as sp
from flint import arb, ctx
from scipy.optimize import minimize_scalar
from liu9_binding import certify_equation_parameters, solve_equation_parameters

ctx.prec = max(ctx.prec, 320)
DPS = 100
GRID_N = 801


def fs(x, n=17):
    return format(float(x), f".{n}g")


def ms(x, n=70):
    return mp.nstr(x, n)


def br(x):
    return {
        "ball": str(x),
        "lower": x.lower().str(45, more=True),
        "upper": x.upper().str(45, more=True),
    }


def h_mp(x):
    if x == 0 or x == 1:
        return mp.mpf(0)
    return -x * mp.log(x) - (1 - x) * mp.log1p(-x)


def hp_mp(x):
    return mp.log((1 - x) / x)


def pi_mp(s, t):
    return s * t * (1 + (1 - s) * (1 - t))


def h_arb(x):
    if x == 0 or x == 1:
        return arb(0)
    return -(x * x.log() + (arb(1) - x) * (arb(1) - x).log())


def hp_arb(x):
    return ((arb(1) - x) / x).log()


def hpp_arb(x):
    return -arb(1) / (x * (arb(1) - x))


def pi_arb(s, t):
    return s * t * (1 + (1 - s) * (1 - t))


def h_np(x):
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(x)
    mask = (x > 0) & (x < 1)
    z = x[mask]
    out[mask] = -z * np.log(z) - (1 - z) * np.log1p(-z)
    return float(out) if out.ndim == 0 else out


def pi_np(s, t):
    return s * t * (1 + (1 - s) * (1 - t))


def symbolic_facts():
    s, t = sp.symbols("s t", real=True)
    pi = lambda u, v: sp.expand(u * v * (1 + (1 - u) * (1 - v)))
    q = lambda u, v: sp.expand(1 - pi(u, v))
    complement = sp.factor(q(s, t) ** 2 - q(s, s) * q(t, t))
    complement_expected = (s - t) ** 2 * (
        (s + t - 1) ** 2 + 1 - s**2 * t**2
    )
    assert sp.expand(complement - complement_expected) == 0
    gram = sp.factor(pi(s, s) * pi(t, t) - pi(s, t) ** 2)
    assert sp.expand(gram - s**2 * t**2 * (s - t) ** 2) == 0

    Axy, Ayx, Bxx, Byy = sp.symbols("Axy Ayx Bxx Byy")
    px, py, qx, qy = sp.symbols("px py qx qy")
    m_entry = (
        Axy - px * qy + Ayx - qx * py
        + Bxx + px * qx + Byy + py * qy
    )
    n_entry = Axy + Ayx + Bxx + Byy
    assert sp.expand(m_entry - n_entry - (px - py) * (qx - qy)) == 0
    return {
        "status": "PROVED exact symbolic algebra",
        "protocol_deficit": (
            "pi(s,s)pi(t,t)-pi(s,t)^2=s^2t^2(s-t)^2"
        ),
        "log_complement_displacement": (
            "[1-pi(s,t)]^2-[1-pi(s,s)][1-pi(t,t)]="
            "(s-t)^2((s+t-1)^2+1-s^2t^2)>=0"
        ),
        "rank_one_decomposition_residual": "0 in free symbols",
    }


def arb_evidence():
    rows = []
    mutants = []
    for bits in (320, 512):
        ctx.prec = bits
        p = certify_equation_parameters(solve_equation_parameters(110))
        s, t = arb(3) / 8, arb(7) / 8
        kss = h_arb(pi_arb(s, s))
        ktt = h_arb(pi_arb(t, t))
        kst = h_arb(pi_arb(s, t))
        energy = kss + ktt - 2 * kst
        c = p.beta * energy
        flipped = -c
        no_cross = p.beta * (kss + ktt)
        a = 1 - arb(1) / (2 * p.mean)
        kappa = 2 * p.beta / a

        x = p.x
        z = pi_arb(x, x)
        zs = x * (1 + (1 - 2 * x) * (1 - x))
        zst = 2 - 4 * x + 4 * x * x
        k12 = hpp_arb(z) * zs * zs + hp_arb(z) * zst
        prod = x * x
        p11 = ((1 - p.beta) * hpp_arb(prod) * x * x
               - x * hpp_arb(x) / (2 * p.mean))
        p12 = ((1 - p.beta) * (hpp_arb(prod) * x * x + hp_arb(prod))
               - hp_arb(x) / p.mean)
        k11 = hpp_arb(z) * zs * zs + hp_arb(z) * (-2 * x * (1 - x))
        r11 = p11 + p.beta * k11
        r12 = p12 + p.beta * k12
        local = -4 * p.beta * k12 / (r11 - r12)
        rows.append({
            "bits": bits,
            "cpd_witness_K_energy": br(energy),
            "cpd_witness_c": br(c),
            "cpd_witness_negative": bool(energy.upper() < 0 and c.upper() < 0),
            "negative_direction": (
                "upper(c)<0: exact c is <= the Arb upper endpoint"
            ),
            "kappa_boundary": br(kappa),
            "kappa_between_1_and_2": bool(
                kappa.lower() > 1 and kappa.upper() < 2
            ),
            "kappa_direction": (
                "upper(kappa)<2 proves the required upper inequality; "
                "lower(kappa)>1 proves the other side"
            ),
            "K12_at_xstar": br(k12),
            "xstar_local_ratio": br(local),
            "local_below_boundary": bool(local.upper() < kappa.lower()),
        })
        for name, value in (
            ("flip_channel_sign", flipped),
            ("drop_minus_2_cross_term", no_cross),
        ):
            mutants.append({
                "bits": bits,
                "name": name,
                "sound_upper": c.upper().str(30, more=True),
                "mutant_upper": value.upper().str(30, more=True),
                "sound_negative_guard": bool(c.upper() < 0),
                "mutant_negative_guard": bool(value.upper() < 0),
                "fired": bool(c.upper() < 0 and not (value.upper() < 0)),
            })
    if not all(r["cpd_witness_negative"] for r in rows):
        raise AssertionError("CPD witness not certified")
    if not all(r["kappa_between_1_and_2"] for r in rows):
        raise AssertionError("kappa enclosure failed")
    if not all(r["fired"] for r in mutants):
        raise AssertionError("mandatory mutant survived")
    return {
        "label": "CERTIFIED Arb",
        "signed_measure": "d=delta_(3/8)-delta_(7/8), sum(d)=0",
        "precisions": rows,
        "mutations": mutants,
        "mutations_all_fail": True,
    }


def pair_discovery(p):
    with mp.workdps(DPS):
        beta, mean = p.beta, p.mean
        pd = lambda s: pi_mp(s, s)
        pdp = lambda s: 2 * s * (1 + (1 - s) ** 2) - 2 * s * s * (1 - s)
        deriv = lambda s: hp_mp(pd(s)) * pdp(s) - 2 * hp_mp(s)
        smin = mp.findroot(deriv, (mp.mpf("0.3"), mp.mpf("0.4")))
        defect = h_mp(pd(smin)) - 2 * h_mp(smin)

        s0, t0 = mp.mpf(3) / 8, mp.mpf(7) / 8
        wss = mp.log1p(-pi_mp(s0, s0))
        wtt = mp.log1p(-pi_mp(t0, t0))
        wst = mp.log1p(-pi_mp(s0, t0))
        mid = (wss + wtt) / 2
        g = lambda w: h_mp(1 - mp.exp(w))
        mid_def = g(wss) + g(wtt) - 2 * g(mid)
        shift_def = 2 * (g(mid) - g(wst))
        pc = mp.findroot(
            lambda z: mp.log(z / (1 - z)) - 1 / z,
            (mp.mpf("0.75"), mp.mpf("0.85")), solver="anderson",
        )

        a = 1 - 1 / (2 * mean)
        kappa = 2 * beta / a
        approach = []
        for power in (2, 4, 8, 16):
            s = mp.power(10, -power)
            c = beta * (h_mp(pd(s)) - 2 * h_mp(s))
            ratio = -c / (a * h_mp(s))
            approach.append({
                "s": f"1e-{power}",
                "minus_c_over_R": ms(ratio, 55),
                "kappa_minus_ratio": ms(kappa - ratio, 18),
            })

        phihalf2 = beta * h_mp(pd(mp.mpf(1) / 2))
        psmall = mp.mpf("1e-12")
        family_ratio = (
            kappa - psmall * phihalf2 / (a * h_mp(mp.mpf(1) / 2))
        )
        return {
            "label": "DISCOVERY mpmath 100 dps",
            "continuum_pair_minimum_candidate": {
                "supports": [ms(smin, 80), "1"],
                "stationarity_residual": ms(deriv(smin), 8),
                "K_energy": ms(defect, 80),
                "beta_scaled_c": ms(beta * defect, 80),
                "warning": "not a global Arb minimum certificate",
            },
            "log_complement_test_at_3_8_7_8": {
                "w_ss": ms(wss, 45),
                "w_tt": ms(wtt, 45),
                "w_st": ms(wst, 45),
                "midpoint": ms(mid, 45),
                "w_st_minus_midpoint": ms(wst - mid, 45),
                "midpoint_convexity_defect": ms(mid_def, 45),
                "cross_log_displacement_defect": ms(shift_def, 45),
                "sum": ms(mid_def + shift_def, 45),
                "g_curvature_change_probability": ms(pc, 55),
                "g_curvature_change_w": ms(mp.log(1 - pc), 55),
            },
            "kappa_formula": "4*m*beta/(2*m-1)",
            "kappa_value": ms(kappa, 80),
            "boundary_approach": approach,
            "endpoint_mixture_near_sharp": {
                "measures": (
                    "mu=(1-p)delta_0+p delta_(1/2),nu=delta_1,p=1e-12"
                ),
                "ratio": ms(family_ratio, 65),
                "exact_formula": (
                    "kappa-p*beta*K(s,s)/(a*h(s)), a=1-1/(2m)"
                ),
            },
            "naive_tensorization_obstruction": {
                "measure": "lambda=(delta_0+delta_(1/2))/2",
                "self_B_minus_variance_phi": ms(-phihalf2 / 4, 55),
                "meaning": (
                    "negative, so integrating point-mass comparisons does "
                    "not prove the general-measure comparison"
                ),
            },
        }


def matrix_discovery(p):
    beta, mean, xstar = float(p.beta), float(p.mean), float(p.x)
    K = lambda s, t: h_np(pi_np(s, t))
    P = lambda s, t: ((1 - beta) * h_np(s * t)
                      - (t * h_np(s) + s * h_np(t)) / (2 * mean))
    grid = np.linspace(0, 1, GRID_N)
    col = grid[:, None]
    kmat = K(col, col.T)
    pmat = P(col, col.T)
    rmat = pmat + beta * kmat
    phi = np.sqrt(np.maximum(0, beta * np.diag(kmat)))
    amat = pmat + phi[:, None] * phi[None, :]
    bmat = beta * kmat - phi[:, None] * phi[None, :]

    row = kmat.mean(1)
    centered = kmat - row[:, None] - row[None, :] + kmat.mean()
    eig, vec = np.linalg.eigh(centered)
    d = vec[:, 0].copy()
    d -= d.mean()
    d /= np.linalg.norm(d)
    if d @ grid < 0:
        d = -d
    pos, neg = np.maximum(d, 0), np.maximum(-d, 0)
    pos /= pos.sum()
    neg /= neg.sum()
    signed = pos - neg

    def structure(w):
        top = np.argsort(w)[-8:][::-1]
        active = np.flatnonzero(w >= 0.01 * w.max())
        return {
            "mean": fs(w @ grid),
            "effective_atoms": fs(1 / (w @ w)),
            "one_percent_peak_range": [fs(grid[active[0]]), fs(grid[active[-1]])],
            "top_eight": [[fs(grid[i]), fs(w[i])] for i in top],
        }

    off = ~np.eye(GRID_N, dtype=bool)
    boff = np.where(off, bmat, np.inf)
    ib = np.unravel_index(np.argmin(boff), boff.shape)
    interior = amat[1:-1, 1:-1]
    ia0 = np.unravel_index(np.argmin(interior), interior.shape)
    ia = (ia0[0] + 1, ia0[1] + 1)

    cpair = beta * (
        np.diag(kmat)[:, None] + np.diag(kmat)[None, :] - 2 * kmat
    )
    ratio = np.full_like(rmat, -np.inf)
    mask = rmat > 1e-13
    ratio[mask] = -cpair[mask] / rmat[mask]
    ir = np.unravel_index(np.argmax(ratio), ratio.shape)
    kappa = 4 * mean * beta / (2 * mean - 1)
    margin = cpair + kappa * rmat

    rng = np.random.default_rng(20260901)
    max_resid = 0.0
    min_A, min_B, min_N = np.inf, np.inf, np.inf
    max_ratio = -np.inf
    for _ in range(4096):
        xs, ys = rng.random(3), rng.random(3)
        mu = rng.dirichlet(np.full(3, 0.7))
        nu = rng.dirichlet(np.full(3, 0.7))
        kxx = K(xs[:, None], xs[None, :])
        kyy = K(ys[:, None], ys[None, :])
        kxy = K(xs[:, None], ys[None, :])
        pxy = P(xs[:, None], ys[None, :])
        qx = np.sqrt(beta * np.diag(kxx))
        qy = np.sqrt(beta * np.diag(kyy))
        Axy = pxy + qx[:, None] * qy[None, :]
        Bxx = beta * kxx - qx[:, None] * qx[None, :]
        Byy = beta * kyy - qy[:, None] * qy[None, :]
        direct = beta * mu @ kxx @ mu + 2 * mu @ pxy @ nu + beta * nu @ kyy @ nu
        split = (mu @ Bxx @ mu + nu @ Byy @ nu + 2 * mu @ Axy @ nu
                 + (mu @ qx - nu @ qy) ** 2)
        max_resid = max(max_resid, abs(direct - split))
        min_A = min(min_A, float(Axy.min()))
        min_B = min(min_B, float(Bxx.min()), float(Byy.min()))
        N = Axy + Axy.T + Bxx + Byy
        v = qx - qy
        M = pxy + pxy.T + beta * kxx + beta * kyy
        max_resid = max(max_resid, float(np.max(np.abs(M - N - np.outer(v, v)))))
        min_N = min(min_N, float(N.min()))
        rcross = mu @ (pxy + beta * kxy) @ nu
        chan = beta * (mu @ kxx @ mu + nu @ kyy @ nu - 2 * mu @ kxy @ nu)
        if rcross > 1e-13 and chan < 0:
            max_ratio = max(max_ratio, float(-chan / rcross))

    half = 0.5
    dcross = beta * h_np(half)
    ddiag = beta * K(half, half)
    badS = np.array([
        [0, 0, -dcross], [0, 0, -dcross], [-dcross, -dcross, ddiag]
    ])
    bind_d = beta * K(xstar, xstar)

    qhalf = np.sqrt(beta * K(half, half))
    def lower(t):
        qt = np.sqrt(beta * K(t, t))
        return -1e100 if qt == 0 else -P(half, t) / (qhalf * qt)
    ts = np.linspace(1e-8, 1 - 1e-8, 20001)
    vals = np.array([lower(t) for t in ts])
    j = int(np.argmax(vals))
    fit = minimize_scalar(
        lambda t: -lower(float(t)),
        bounds=(ts[max(0, j - 3)], ts[min(len(ts) - 1, j + 3)]),
        method="bounded", options={"xatol": 1e-15},
    )

    return {
        "label": "DISCOVERY deterministic float64",
        "grid": {"points": GRID_N, "spacing": "1/800"},
        "centered_K": {
            "lambda_min_raw": fs(eig[0]),
            "lambda_min_divided_by_N": fs(eig[0] / GRID_N),
            "lambda_max_raw": fs(eig[-1]),
            "negative_count_below_minus_1e-10": int(np.sum(eig < -1e-10)),
            "positive_count_above_1e-10": int(np.sum(eig > 1e-10)),
            "sign_changes": int(np.sum(np.signbit(d[1:]) != np.signbit(d[:-1]))),
            "positive_part": structure(pos),
            "negative_part": structure(neg),
            "K_energy_unit_positive_negative_mass": fs(signed @ kmat @ signed),
            "beta_scaled_c": fs(beta * signed @ kmat @ signed),
            "minus_c_over_R": fs(-(beta * signed @ kmat @ signed) / (pos @ rmat @ neg)),
        },
        "sandwich": {
            "A_min_grid": fs(amat.min()),
            "A_min_interior": {"value": fs(amat[ia]), "point": [fs(grid[ia[0]]), fs(grid[ia[1]])]},
            "B_min_grid": fs(bmat.min()),
            "B_min_off_diagonal": {"value": fs(bmat[ib]), "point": [fs(grid[ib[0]]), fs(grid[ib[1]])]},
            "warning": "roundoff-scale values are not certificates",
        },
        "decomposition_checks": {
            "draws": 4096,
            "max_abs_residual": fs(max_resid),
            "min_A": fs(min_A),
            "min_B": fs(min_B),
            "min_N_entry": fs(min_N),
        },
        "phi_slack": {
            "s": "1/2",
            "canonical_phi": fs(qhalf),
            "off_diagonal_minimum_scale": fs(-fit.fun),
            "minimum_isolated_scale_including_diagonal": fs(max(
                -fit.fun, np.sqrt(max(0, -P(half, half) / qhalf**2))
            )),
            "maximizing_t": fs(fit.x),
            "example_scale": "0.95",
            "diagonal_lower_scale": fs(np.sqrt(max(0, -P(half, half) / qhalf**2))),
            "upper_scale": "1",
        },
        "failed_natural_splits": {
            "channel_matrix_geometry": "((0,1),(0,1),(1/2,0))",
            "channel_matrix_lambda_min": fs(np.linalg.eigvalsh(badS)[0]),
            "sign_clip_geometry": "((xstar,0),(xstar,0),(0,xstar))",
            "sign_clip_lambda_formula": "beta*K(xstar,xstar)*(1-sqrt(2))",
            "sign_clip_lambda": fs(bind_d * (1 - np.sqrt(2))),
        },
        "route4": {
            "max_pair_grid_ratio": fs(ratio[ir]),
            "pair_grid_point": [fs(grid[ir[0]]), fs(grid[ir[1]])],
            "min_c_plus_kappa_R": fs(margin.min()),
            "max_4096_random_measure_ratio": fs(max_ratio),
        },
        "block_psd": {
            "R_grid_lambda_min": fs(np.linalg.eigvalsh(beta * kmat + pmat)[0]),
            "twoQ_minus_R_grid_lambda_min": fs(np.linalg.eigvalsh(beta * kmat - pmat)[0]),
        },
    }


def bernstein(p):
    with mp.workdps(DPS):
        def ek(n): return h_mp(1 / mp.sqrt(n))
        def ep(n):
            return ((1 - p.beta) * ek(n) + h_mp(1 / (2 * mp.sqrt(n))) / p.mean
                    + mp.log(2) / (2 * p.mean * mp.sqrt(n)))
        def er(n): return ep(n) + p.beta * ek(n)
        def et(n): return 2 * p.beta * ek(n) + 2 * ep(n)
        values = [{
            "n": n, "K_tail": ms(ek(n), 25), "P_tail": ms(ep(n), 25),
            "R_tail": ms(er(n), 25), "T_tail": ms(et(n), 25),
            "label": "DISCOVERY evaluation of rigorous formula",
        } for n in (64, 256, 1024, 4096, 65536, 10**6)]
        orders = []
        for target in map(mp.mpf, ("1e-1", "1e-2", "1e-3", "1e-6")):
            lo, hi = 4, 8
            while et(hi) > target:
                lo, hi = hi, 2 * hi
            while lo + 1 < hi:
                mid = (lo + hi) // 2
                if et(mid) <= target: hi = mid
                else: lo = mid
            orders.append({"target": ms(target, 4), "first_n": hi, "tail": ms(et(hi), 25)})

    ctx.prec = 320
    pa = certify_equation_parameters(solve_equation_parameters(110))
    arb_rows = []
    for n in (64, 256, 1024, 4096):
        root = arb(n).sqrt()
        e = (2 * h_arb(arb(1) / root)
             + 2 * h_arb(arb(1) / (2 * root)) / pa.mean
             + arb(2).log() / (pa.mean * root))
        arb_rows.append({
            "n": n, "T_tail_expression": br(e),
            "direction": (
                "error<=exact expression<=Arb upper endpoint, so the upper endpoint is the rigorous side"
            ),
        })
    return {
        "basis": "tensor-product Bernstein degree n",
        "rigorous_formulas": {
            "K": "h(1/sqrt(n))",
            "P": "(1-beta)h(1/sqrt(n))+h(1/(2sqrt(n)))/m+log(2)/(2m sqrt(n))",
            "R": "h(1/sqrt(n))+h(1/(2sqrt(n)))/m+log(2)/(2m sqrt(n))",
            "T": "2h(1/sqrt(n))+2h(1/(2sqrt(n)))/m+log(2)/(m sqrt(n))",
            "proof": (
                "|pi_s|,|pi_t|<=1; E|Bin(n,x)/n-x|<=1/(2sqrt(n)); "
                "|h(u)-h(v)|<=h(|u-v|); Jensen"
            ),
        },
        "values": values,
        "orders_needed": orders,
        "arb_320_upper_bounds": arb_rows,
        "obstruction": (
            "inf T=0, so no finite additive uniform tail can prove T>=0 without a separate zero-local factorization"
        ),
    }


def main():
    mp.mp.dps = DPS
    p = solve_equation_parameters(DPS)
    exact = symbolic_facts()
    cert = arb_evidence()
    pair = pair_discovery(p)
    matrix = matrix_discovery(p)
    tails = bernstein(p)
    report = {
        "claim": "Liu H2 channel-theory routes",
        "parameters": {"beta": ms(p.beta, 80), "m": ms(p.mean, 80), "xstar": ms(p.x, 80)},
        "exact_algebra": exact,
        "route1": {
            "verdict": "REFUTED",
            "certified_witness": cert,
            "continuum_discovery": pair["continuum_pair_minimum_candidate"],
            "log_complement": pair["log_complement_test_at_3_8_7_8"],
            "centered_spectrum": matrix["centered_K"],
        },
        "route2": {
            "verdict": "PROVED-CONDITIONAL-ON-THE-SANDWICH",
            "phi": "sqrt(Q2(s,s))",
            "A": "P2(s,t)+phi(s)phi(t)",
            "B": "Q2(s,t)-phi(s)phi(t)",
            "sandwich": "A>=0 and B>=0, equivalently 0<=B<=R",
            "measure_identity": (
                "T=<mu,B mu>+<nu,B nu>+2<mu,A nu>+(<mu,phi>-<nu,phi>)^2"
            ),
            "matrix_decomposition": (
                "M_ij=N_ij+v_i v_j; v_i=phi(x_i)-phi(y_i); "
                "N_ij=A(x_i,y_j)+A(x_j,y_i)+B(x_i,x_j)+B(y_i,y_j)"
            ),
            "numeric_evidence": {"sandwich": matrix["sandwich"], "checks": matrix["decomposition_checks"]},
            "phi_forcedness": {
                "diagonal_interval": "max(0,-P2(s,s))<=phi(s)^2<=Q2(s,s)",
                "forced_at": "s=0,xstar,1; at xstar -P2=Q2",
                "not_forced_elsewhere": matrix["phi_slack"],
                "consequence": "slack away from xstar cannot remove the forced zero at (xstar,xstar)",
            },
            "failed_obvious_splits": matrix["failed_natural_splits"],
            "closure": "dimension-free once the separately owned sandwich certificates are complete",
        },
        "route3": {
            "verdict": "OPEN-WITH-OBSTRUCTION",
            "bernstein_tail": tails,
            "psd_failure": matrix["block_psd"],
        },
        "route4": {
            "verdict": "OPEN-WITH-OBSTRUCTION",
            "kappa_candidate_formula": pair["kappa_formula"],
            "kappa_candidate": pair["kappa_value"],
            "arb_kappa_and_local_ratio": cert["precisions"],
            "boundary_approach": pair["boundary_approach"],
            "pair_and_measure_search": matrix["route4"],
            "endpoint_measure_family": pair["endpoint_mixture_near_sharp"],
            "tensorization_obstruction": pair["naive_tensorization_obstruction"],
            "status": (
                "necessary lower bound and kappa<2 certified; global pair upper bound and general-measure extension remain open"
            ),
        },
        "mutations": cert["mutations"],
        "mutations_all_fail": cert["mutations_all_fail"],
        "limitations": [
            "float64/mpmath extrema are DISCOVERY only",
            "this theory file does not duplicate either scalar Arb box certificate",
            "the continuum pair minimum and global kappa upper bound are not certified here",
        ],
    }
    body = json.dumps(report, sort_keys=True, separators=(",", ":"))
    report["payload_sha256"] = hashlib.sha256(body.encode()).hexdigest()
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
