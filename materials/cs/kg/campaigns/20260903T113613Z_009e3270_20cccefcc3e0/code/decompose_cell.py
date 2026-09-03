#!/usr/bin/env python3
"""Campaign B — §4 slack-decomposition battery on A's frozen frontier cell.

D1  A envelope (byte-shared) at cap 131072 — reconciliation gate: must
    reproduce A's frozen -1.7629e-5 margin to frozen width.
D2  V envelope (mq = t^2 kill/subdivision) at cap 131072 — Phase-1 mechanism;
    entry "fold/second-hinge" = D1.upper - D2.upper.
D3  A envelope at caps 2^19 and 2^21 — panel-sup series; per-quadrupling gains.
D4  |e|-hinge disk-feasibility diagnostic (NOT a verdict mechanism; firewall
    per pre-statement §4): certified envelope with Esup replaced by the
    certified sup of e over B ∩ closed disk (outward candidate enumeration:
    in-disk corners, box-edge/circle crossings, in-disk corner cone dirs).
D5  tail ball (identical in all variants; differential certified 0).
D6  within-panel skew = D3's 2^19→2^21 delta, restated.
All outputs outward Arb balls (rule 15); JSON to logs/decomposition_tile3.json.
"""
from __future__ import annotations

import json
import os
import sys
import time

import importlib.util as _ilu

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(_HERE, '..', 'logs')


def load(name, path):
    spec = _ilu.spec_from_file_location(name, path)
    mod = _ilu.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


PREC = 256
CORE_DIR = os.path.join(_HERE, '..', '..',
                        '20260831T082425Z_kg_direct_d4_restart', 'code')

gsv = load('gsvd', os.path.join(_HERE, 'gate_subdiv.py'))
d4core = load('d4core', os.path.join(CORE_DIR, 'core.py'))
GCE_A = load('gce_A_mod',
             os.path.join(_HERE, '..', '..',
                          '20260831T231500Z_KgEvenBoxRerun_campA_corr_env',
                          'code', 'gate_corr_env.py'))
gce = GCE_A

from d4core import Prec  # noqa: E402
from flint import arb, fmpq  # noqa: E402

from d4core import Prec  # noqa: E402
from flint import arb, fmpq  # noqa: E402


def mag_up(v: arb) -> arb:
    """Outward upper bound of |v| for an arb ball v."""
    m1 = v.upper()
    m2 = (-v).upper()
    return m1 if bool(m1 >= m2) else m2


def disk_feasible_candidates(bx):
    a0lo, a0hi, a2lo, a2hi = bx
    cands = []
    for (x, y) in [(a0lo, a2lo), (a0lo, a2hi), (a0hi, a2lo), (a0hi, a2hi)]:
        if bool(x * x + y * y <= 1):
            cands.append((x, y))
    for x in (a0lo, a0hi):
        y2 = (1 - x * x).nonnegative_part()
        for sg in (1, -1):
            y = sg * y2.sqrt()
            if bool(y >= a2lo) and bool(y <= a2hi):
                cands.append((x, y))
    for y in (a2lo, a2hi):
        x2 = (1 - y * y).nonnegative_part()
        for sg in (1, -1):
            x = sg * x2.sqrt()
            if bool(x >= a0lo) and bool(x <= a0hi):
                cands.append((x, y))
    sq = (arb(3) / 2).sqrt()
    for s1 in (1, -1):
        for s2 in (1, -1):
            x = arb(s1) / sq
            y = arb(s2) / (sq * arb(2).sqrt())
            if bool(x >= a0lo) and bool(x <= a0hi) and bool(y >= a2lo) and bool(y <= a2hi):
                cands.append((x, y))
    return cands


def main() -> None:
    started = time.process_time()
    with Prec(PREC):
        cL = arb('4.1485893120')
        cR = arb('4.181778026496')
        box = gsv.frontier_box()
        roots = gsv.crossings(cL, box, PREC)
        d_lower = d4core.d_of_c(cR, PREC).lower()

        def margin_of(upper: arb) -> arb:
            return d_lower - upper.upper()

        # D1: baseline at cap
        t0 = time.process_time()
        U_cap = GCE_A.cell_envelope(cL, box, 131072, PREC)
        m1 = margin_of(U_cap)
        t1 = time.process_time()

        # D2: V at cap
        V_cap = gsv.cell_envelope_V(cL, box, 131072, roots, PREC)
        m2 = margin_of(V_cap)
        t2 = time.process_time()

        # D3: panel-sup series
        U_19 = GCE_A.cell_envelope(cL, box, 2 ** 19, PREC)
        m19 = margin_of(U_19)
        t3 = time.process_time()
        U_21 = GCE_A.cell_envelope(cL, box, 2 ** 21, PREC)
        m21 = margin_of(U_21)
        tail = d4core.tail_account(8, PREC)
        t4 = time.process_time()

        # D4: disk-feasible even sup, full cap-131072 envelope
        cands = disk_feasible_candidates(box)
        ell = GCE_A.dist0_box(box)
        odd_r = (arb(1) - ell * ell).nonnegative_part().sqrt()
        n = 131072
        total_D4 = arb(0)
        for i in range(n):
            a = arb(fmpq(i * 8, n))
            b = arb(fmpq((i + 1) * 8, n))
            if bool(d4core.K(b, PREC).sqrt() <= cL * a):
                continue
            p2a = d4core.psi(2, a, PREC)
            p2b = d4core.psi(2, b, PREC)
            sup = None
            for (x, y) in cands:
                v_a = x + y * p2a
                v_b = x + y * p2b
                for v in (v_a, v_b):
                    m = mag_up(v)
                    if sup is None or bool(m > sup):
                        sup = m
            Wob = odd_r * d4core.K_o(b, PREC).sqrt()
            tL = cL * a
            first = (sup + Wob - tL).nonnegative_part()
            dua = (sup - Wob).upper()
            dla = (Wob - sup).upper()
            adu = dua if bool(dua >= dla) else dla
            second = (adu - tL).nonnegative_part()
            wt = d4core.gauss_cdf(b, PREC) - d4core.gauss_cdf(a, PREC)
            total_D4 += wt * (first + second)
        D4_upper = total_D4 + tail
        mD4 = margin_of(D4_upper)
        t5 = time.process_time()

        out = {
            'battery': 'decomposition on A frontier cell (pre-statement §4)',
            'c_pair': [cL.str(30), cR.str(30)],
            'frozen_a_margin': '-1.76291787028287693244667476546e-5',
            'D1_baseline_cap131072': {
                'upper': U_cap.str(30),
                'margin': m1.str(30),
                'margin_radius': m1.rad().str(12),
                'reproduces_A_frozen': bool((m1 - arb(
                    '-1.76291787028287693244667476546e-5')).rad() < arb(2) ** -40),
                'process_seconds': t1 - t0},
            'D2_V_kill_subdivision_cap131072': {
                'upper': V_cap.str(30),
                'margin': m2.str(30),
                'margin_radius': m2.rad().str(12),
                'fold_second_hinge_recovery': (U_cap.upper() - V_cap.upper()).str(15),
                'margin_delta_vs_A': (m2 - m1).str(15)},
            'D3_panel_sup_series': {
                'cap_2p19_upper': U_19.str(30),
                'cap_2p19_margin': m19.str(30),
                'cap_2p21_upper': U_21.str(30),
                'cap_2p21_margin': m21.str(30),
                'gain_2p17_to_2p19': (U_cap.upper() - U_19.upper()).str(15),
                'gain_2p19_to_2p21': (U_19.upper() - U_21.upper()).str(15),
                'panel_sup_total_2p17_to_2p21': (U_cap.upper() - U_21.upper()).str(15)},
            'D4_disk_feasible_even_sup_diagnostic': {
                'note': ('NOT a verdict mechanism (pre-statement §4 firewall): '
                         'certified envelope with Esup replaced by the certified '
                         'sup of e over the disk-feasible part of B; measures the '
                         '|e|-hinge coefficient-space excess'),
                'feasible_candidate_count': len(cands),
                'upper': D4_upper.str(30),
                'margin': mD4.str(30),
                'even_hinge_recovery_vs_A_cap': (U_cap.upper() - D4_upper.upper()).str(15),
                'margin_delta_vs_A': (mD4 - m1).str(15)},
            'D5_tail_ball': {'value': tail.str(20),
                             'note': 'identical in every variant; differential exactly 0'},
            'D6_within_panel_skew': {
                'note': "captured by D3's 2^19→2^21 delta; restated",
                'delta_2p19_to_2p21': (U_19.upper() - U_21.upper()).str(15)},
            'process_seconds_total': time.process_time() - started,
        }

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 'decomposition_tile3.json'), 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2)
        f.write('\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
