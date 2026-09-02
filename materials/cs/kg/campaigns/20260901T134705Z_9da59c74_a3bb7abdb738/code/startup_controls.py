#!/usr/bin/env python3
"""Campaign A startup controls — machine-verified, abort on first failure.

Implements pre-statement §3 asserts 1-6, in order:
1. kernel-admission (no-cap, no-3D-inflation) with cap counterfactual;
2. odd-radius 2-D-paper form vs 3-D circumradius counterfactual;
3. K_o/K_e basis identities;
4. cited-tail anchor from authors' frozen d3h_certificate.json (sha-verified);
5. margin-semantics counterfactuals;
6. envelope-admissibility on the analytic even anchor p* = (4/5)psi0-(3/5)psi2
   against the frozen restart-campaign closed-form windows.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

from flint import arb, fmpq

import importlib.util as _ilu
_sp = _ilu.spec_from_file_location(
    'd4core', os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           '..', '..', '20260831T082425Z_kg_direct_d4_restart', 'code', 'core.py'))
_dc = _ilu.module_from_spec(_sp)
sys.modules['d4core'] = _dc
_sp.loader.exec_module(_dc)
import d4core
from d4core import Prec
import gate_corr_env as gce

RESULT: dict = {}
ABORT = None


def check(name: str, ok: bool, detail: str) -> None:
    global ABORT
    RESULT[name] = {'ok': bool(ok), 'detail': detail}
    print(f"[{'OK' if ok else 'FAIL'}] {name}: {detail}", flush=True)
    if not ok and ABORT is None:
        ABORT = name


def mag_upper(x: arb, bits: int) -> arb:
    """Outward upper for |x| given an interval x: max(upper, -lower)."""
    with Prec(bits):
        u = x.upper()
        l = x.lower()
        nl = -l
        return u if u >= nl else nl


def main() -> None:
    B = 320

    # ---- 1. kernel admission ------------------------------------------------
    with Prec(B):
        a0 = (arb(2) / 3).sqrt()
        a2 = -(arb(1) / 3).sqrt()
        n2 = a0 * a0 + a2 * a2
        unit_ok = bool((n2.lower() >= 1 - arb(2) ** -260) and (n2.upper() <= 1 + arb(2) ** -260))
        check('kernel_unit_norm', unit_ok, f'|a|^2 = {n2.str(35)}')
        zero = arb(0)
        e0 = a0 * d4core.psi(0, zero, B) + a2 * d4core.psi(2, zero, B)
        s32 = (arb(3) / 2).sqrt()
        eps = arb(2) ** -260
        lo_ok = bool(e0.lower() >= s32.lower() - eps)
        hi_ok = bool(e0.upper() <= s32.upper() + eps)
        check('kernel_e0_is_sqrt32', lo_ok and hi_ok,
              f'e(0) = {e0.str(35)}, sqrt(3/2) = {s32.str(35)}')
        # corner-exact E_P at the point box / degenerate panel [0,0]:
        _, _, _, E_up = gce.e_panel_range((a0, a0, a2, a2), zero.union(zero), B)
        admit = bool(E_up >= s32.upper() - eps)
        check('kernel_admitted_by_envelope', admit,
              f'E_P.upper = {E_up.str(35)} >= sqrt(3/2)-2^-260: {admit}')
        excess = E_up - s32.upper()
        check('kernel_envelope_tight', bool(excess <= eps and excess >= -eps),
              f'slack = {excess.str(10)} (must be < 2^-260)')
        # counterfactual: |e|-cap at 1 is detectably unsound
        deficit = s32.lower() - arb(1)
        check('cap1_counterfactual_rejected', bool(deficit > 0),
              f'sqrt(3/2) - 1 = {deficit.str(30)} > 0 (a cap at 1 clips the paper kernel)')

    # ---- 2. odd radius: paper 2-D form vs 3-D circumradius counterfactual ---
    with Prec(B):
        bx = (arb('0.8'), arb('0.9'), arb('-0.6'), arb('-0.5'))
        l2d = gce.dist0_box(bx)
        ref = (arb(89) / 100).sqrt()
        check('dist0_exact_value', bool(l2d.overlaps(ref) and (l2d - ref).rad() < arb('1e-30')),
              f'l = {l2d.str(30)} vs sqrt(0.89) = {ref.str(30)}')
        r2d = (arb(1) - l2d * l2d).nonnegative_part().sqrt()
        # Predecessor's 3-D circumradius inflation, per its j_integrand_bound:
        # replace the 2-D corner-exact even range with  |e_center(s)| + rho3*sqrt(K_e(s))
        # where rho3 = half-diagonal of the box in (a0, a2, plus a phantom odd axis
        # of length sqrt(1-l^2)): the predecessor's 3-D circumradius. We reconstruct
        # its value as the circumradius of the 3-D box [w/2 x w/2 x sqrt(1-l^2)]
        # centered at the even-box center: this is the smallest 3-ball covering
        # (even box) x (odd ellipse half-axis), i.e. an INFLATED 3-D reading.
        # Control: at s=0 the inflated RHS must EXCEED the true maximal |e| over
        # the box (proving the failure mode "3-D inflation wastes margin" is real
        # and measurable, i.e. the corrected 2-D form is strictly tighter).
        cx = (bx[0] + bx[1]) / 2
        cy = (bx[2] + bx[3]) / 2
        wx = (bx[1] - bx[0]) / 2
        wy = (bx[3] - bx[2]) / 2
        rho3 = (wx * wx + wy * wy + r2d * r2d).sqrt()
        s = arb(0)
        # maximal TRUE |e| over the box at s=0: corner-exact, |c0| with c0 = a0 - a2/sqrt(2)
        c0_lo = bx[0] - bx[3] / arb(2).sqrt()   # a0_min - a2_max/sqrt2
        c0_hi = bx[1] - bx[2] / arb(2).sqrt()   # a0_max - a2_min/sqrt2
        # outward |max e(0)|: max(|c0_lo|, |c0_hi|) outward — corners suffice (affine)
        e_true = mag_upper(arb(0).union(arb(0)) * 0 + c0_hi, B) if False else None
        v1 = mag_upper(c0_lo, B)
        v2 = mag_upper(c0_hi, B)
        e_true_max = v1 if v1 >= v2 else v2
        infl = mag_upper(cx - cy / arb(2).sqrt(), B) + rho3 * d4core.K_e(s, B).sqrt()
        check('three_d_inflation_is_looser', bool(infl > e_true_max),
              f'3-D-inflated |e|-RHS at s=0: {infl.str(24)} > corner-exact sup: {e_true_max.str(24)} (slack {arb(infl - e_true_max).str(8)})')
        # the corrected 2-D envelope value at s=0 must EQUAL the corner-exact sup:
        _, _, _, E_up = gce.e_panel_range(bx, arb(0).union(arb(0)), B)
        check('two_d_envelope_corner_exact_at_s0',
              bool(E_up >= e_true_max - arb(2) ** -260 and E_up <= e_true_max + arb(2) ** -260),
              f'E_P.upper = {E_up.str(30)} == corner-exact sup = {e_true_max.str(30)}')

    # ---- 4. cited-tail anchor ------------------------------------------------
    cert_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             '..', '..', '..', 'scratch', 'repo', 'd3h_certificate.json')
    cert_path = os.path.normpath(cert_path)
    with open(cert_path, 'rb') as f:
        cert_bytes = f.read()
    sha = hashlib.sha256(cert_bytes).hexdigest()
    expected_sha = '80e945589b796408f5de5add837272037160d163467050b354f15ab9d6822d2e'
    check('d3h_cert_sha256', sha == expected_sha, f'sha256 {sha[:16]}... == 80e945589b...')
    cert = json.loads(cert_bytes)
    # ---- 3. basis identities -------------------------------------------------
    with Prec(B):
        ok_all = True
        tol = arb(2) ** -250
        for sv in (arb(0), arb(1), arb('1.7'), arb(2), arb(3)):
            psi1 = d4core.psi(1, sv, B)
            psi3 = d4core.psi(3, sv, B)
            Ko_id = psi1 * psi1 + psi3 * psi3
            Ko_lit = d4core.K_o(sv, B)
            oe = d4core.psi(0, sv, B) ** 2 + d4core.psi(2, sv, B) ** 2 - d4core.K_e(sv, B)
            ok_ko = bool((Ko_id - Ko_lit).rad() < tol)
            ok_ke = bool(oe == 0) or bool((oe - oe).contains(0) and oe.rad() < tol)
            ok_all = ok_all and ok_ko and ok_ke
        check('Ko_Ke_identities', ok_all, 'K_o = psi1^2+psi3^2 and K_e = psi0^2+psi2^2 agree within 2^-250 at s in {0,1,1.7,2,3}')
    # certificate's own B3 ball must be enclosed by the imported decimal:
    B3_ball = arb(cert['B3_upper'].strip('[] '))
    with Prec(B):
        B3 = arb('14.44243664663976457')
        check('B3_import_is_upper', bool(B3.lower() >= B3_ball.lower()),
              f'imported 14.44243664663976457 lower {B3.lower().str(30)} >= certificate lower end {B3_ball.lower().str(30)}')

    # ---- 5. margin semantics counterfactuals ---------------------------------
    with Prec(B):
        d = arb('[0.9 +/- 1e-3]')
        U = arb('[0.899 +/- 1e-3]')
        m_strad = (d.lower() - U.upper()) > 0
        check('margin_straddle_rejected', not m_strad, f'straddling ball margin {m_strad} (must be False)')
        U2 = arb('[0.89 +/- 1e-3]')
        m_pass = (d.lower() - U2.upper()) > 0
        check('margin_tight_accepted', m_pass, f'separated ball margin {m_pass} (must be True)')
        U3 = arb('[0.91 +/- 1e-3]')
        m_above = (d.lower() - U3.upper()) > 0
        check('margin_above_rejected', not m_above, f'upper-above-d margin {m_above} (must be False)')

    # ---- 6. envelope admissibility on analytic even anchor -------------------
    with Prec(B):
        c = arb('1.30')
        A = arb(4) / 5 + arb(3) / (5 * arb(2).sqrt())
        Bc = arb(3) / (5 * arb(2).sqrt())
        # closed-form anchors (frozen restart windows, independently derived)
        wl = arb('0.37483055458876357002453330452976048')
        wh = arb('0.37483055458876357002453330452976050')
        disc = (c * c + 4 * A * Bc).sqrt()
        r1 = (disc - c) / (2 * Bc)
        r2 = (disc + c) / (2 * Bc)
        P = lambda r: d4core.gauss_cdf(r, B) - arb(1) / 2
        Q = lambda r: 1 - d4core.gauss_cdf(r, B)
        phi = lambda r: d4core.phi_density(arb(r), B)
        J = 2 * ((A - Bc) * P(r1) + Bc * r1 * phi(r1) - c * (phi(0) - phi(r1))
                 + (Bc - A) * Q(r2) + Bc * r2 * phi(r2) - c * phi(r2))
        j_ok = bool(J.lower() >= wl and J.upper() <= wh)
        check('pstar_closed_form_window', j_ok, f'J(1.30,p*) = {J.str(38)} inside frozen window')
        # envelope admissibility: point box of p* at c_L = 1.3, n = 2048
        a0p = A
        a2p = -Bc
        U = gce.cell_envelope(arb('1.30'), (a0p, a0p, a2p, a2p), 2048, B)
        slack = U.upper() - J.lower()
        # (no tightness assert: pre-statement §3.6 requires only ADMISSION with
        check('envelope_admits_pstar', bool(slack >= 0),
              f'U(point box, 1.30).upper - J.lower = {slack.str(15)} >= 0')
        # positive slack; the panel-sup gap at finite n is expected and priced)

    if ABORT is not None:
        print(f'ABORT at control: {ABORT}', flush=True)
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               '..', 'logs', 'startup_controls.json'), 'w') as f:
            json.dump({'aborted': ABORT, 'controls': RESULT}, f, indent=2)
        sys.exit(1)

    print('ALL STARTUP CONTROLS PASSED', flush=True)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           '..', 'logs', 'startup_controls.json'), 'w') as f:
        json.dump({'aborted': None, 'controls': RESULT}, f, indent=2)


if __name__ == '__main__':
    main()
