#!/usr/bin/env python3
"""Mandatory analytic and direct-envelope anchors for the restart campaign."""
from __future__ import annotations

import json
import os

import mpmath as mp
from flint import arb

import core
import gate_direct as gate

ANALYTIC_BITS = 300


def inside(ball: arb, low: str, high: str) -> bool:
    return ball.lower() >= arb(low) and ball.upper() <= arb(high)


def analytic_pstar(c_value: str, bits: int = ANALYTIC_BITS):
    """Closed-form J(c,p*) for p*=(4/5)psi0-(3/5)psi2=A-Bs^2."""
    with core.Prec(bits):
        c = arb(c_value)
        sqrt2 = arb(2).sqrt()
        A = arb(4) / 5 + arb(3) / (5 * sqrt2)
        B = arb(3) / (5 * sqrt2)
        discriminant = (c * c + 4 * A * B).sqrt()
        r1 = (discriminant - c) / (2 * B)
        r2 = (discriminant + c) / (2 * B)
        phi0 = core.phi_density(arb(0), bits)
        phi1 = core.phi_density(r1, bits)
        phi2 = core.phi_density(r2, bits)
        P1 = core.gauss_cdf(r1, bits) - arb(1) / 2
        Q2 = arb(1) - core.gauss_cdf(r2, bits)
        inner = ((A - B) * P1 + B * r1 * phi1
                 - c * (phi0 - phi1))
        outer = ((B - A) * Q2 + B * r2 * phi2 - c * phi2)
        J = 2 * (inner + outer)
        return {
            'c': c,
            'A': A,
            'B': B,
            'r1': r1,
            'r2': r2,
            'inner_halfline': inner,
            'outer_halfline': outer,
            'J': J,
        }


def mpmath_quadrature() -> dict:
    """Non-certifying independent quadrature of the two active components."""
    mp.mp.dps = 80
    sqrt2 = mp.sqrt(2)
    A = mp.mpf(4) / 5 + mp.mpf(3) / (5 * sqrt2)
    B = mp.mpf(3) / (5 * sqrt2)
    c = mp.mpf('1.30')
    disc = mp.sqrt(c * c + 4 * A * B)
    r1 = (disc - c) / (2 * B)
    r2 = (disc + c) / (2 * B)
    normal = 1 / mp.sqrt(2 * mp.pi)

    def density(s):
        return normal * mp.exp(-s * s / 2)

    inner = mp.quad(lambda s: (A - B * s * s - c * s) * density(s), [0, r1])
    outer = mp.quad(lambda s: (B * s * s - A - c * s) * density(s), [r2, mp.inf])
    return {
        'digits': 80,
        'r1': mp.nstr(r1, 75),
        'r2': mp.nstr(r2, 75),
        'inner_halfline': mp.nstr(inner, 75),
        'outer_halfline': mp.nstr(outer, 75),
        'J': mp.nstr(2 * (inner + outer), 75),
        'evidence': 'COMPUTATIONAL-EVIDENCE only',
    }


def main() -> None:
    analytic = analytic_pstar('1.30')
    with core.Prec(ANALYTIC_BITS):
        A = analytic['A']
        B = analytic['B']
        c = analytic['c']
        r1 = analytic['r1']
        r2 = analytic['r2']
        J = analytic['J']

        # Exact coefficient/basis convention and active-region semantics.
        norm_residual = (arb(4) / 5) ** 2 + (arb(3) / 5) ** 2 - 1
        assert norm_residual.contains(0)
        assert (A - B - arb(4) / 5).contains(0)
        zero_of_p = (A / B).sqrt()
        q_inner_r1 = A - B * r1 * r1 - c * r1
        q_outer_r2 = B * r2 * r2 - c * r2 - A
        assert q_inner_r1.contains(0), q_inner_r1
        assert q_outer_r2.contains(0), q_outer_r2
        assert 0 < r1 < zero_of_p < r2
        assert A > 0 and B > 0
        assert -2 * B * r1 - c < 0
        assert 2 * B * r2 - c > 0

        # Pre-registered known-value windows.
        assert inside(r1,
                      '0.7554755377609217537989341027304',
                      '0.7554755377609217537989341027306'), r1
        assert inside(r2,
                      '3.8196049229026276928692596718514',
                      '3.8196049229026276928692596718516'), r2
        assert inside(J,
                      '0.37483055458876357002453330452976048',
                      '0.37483055458876357002453330452976050'), J

        d130 = core.d_of_c(arb('1.30'), ANALYTIC_BITS)
        d145 = core.d_of_c(arb('1.45'), ANALYTIC_BITS)
        equal_margin = d130 - J
        coarse_margin = d145 - J
        assert inside(equal_margin,
                      '0.032236520487306532595117305408006',
                      '0.032236520487306532595117305408008'), equal_margin
        assert inside(coarse_margin,
                      '-0.002148854794213660928402532894742',
                      '-0.002148854794213660928402532894740'), coarse_margin
        assert equal_margin.lower() > 0
        assert coarse_margin.upper() < 0

    with core.Prec(gate.PREC):
        roots = gate.root_boxes(12)
        assert len(roots) == 132
        tiles = gate.c_tiles('1.30', '1.45')
        gate.assert_tile_cover('1.30', '1.45', tiles)

        exact_weight = gate.cdf_weight(arb(0), arb(1))
        right_rectangle = core.phi_density(arb(1))
        left_rectangle = core.phi_density(arb(0))
        assert right_rectangle.upper() < exact_weight.lower()
        assert exact_weight.upper() < left_rectangle.lower()

        # Strict-endpoint decision plants around the actual d(1) threshold.
        d_plant = core.d_of_c(arb(1))
        assert gate.strict_margin_lower(
            arb(1), d_plant.lower() - arb('0.0001')) > 0
        assert gate.strict_margin_lower(
            arb(1), d_plant.upper() + arb('0.0001')) < 0
        assert gate.strict_margin_lower(
            arb(1), d_plant + arb(0, arb('0.001'))) < 0

        eps = arb('0.0001')
        box = (arb('0.8') - eps, arb('0.8') + eps,
               arb('-0.6') - eps, arb('-0.6') + eps)
        direct_upper = gate.cell_envelope(arb('1.30'), box, 512)
        analytic_lower_256 = arb(J.str(80)).lower()
        assert direct_upper.upper() >= analytic_lower_256
        planted_pseudo_upper = analytic_lower_256 - arb('0.000001')
        assert not planted_pseudo_upper.upper() >= analytic_lower_256

        equal_direct_margin = gate.strict_margin_lower(arb('1.30'), direct_upper)
        ratio_direct_margin = gate.strict_margin_lower(arb('1.326'), direct_upper)
        assert equal_direct_margin > 0, equal_direct_margin
        assert ratio_direct_margin > 0, ratio_direct_margin

        direct_record = {
            'cell_center': ['0.8', '-0.6'],
            'cell_radius_each_axis': '0.0001',
            'n_panels': 512,
            'upper': direct_upper.str(50),
            'upper_radius': direct_upper.rad().str(20),
            'equal_c_pair': ['1.30', '1.30'],
            'equal_c_certified_margin_lower': equal_direct_margin.str(50),
            'equal_c_margin_radius': equal_direct_margin.rad().str(20),
            'ratio_1p02_pair': ['1.30', '1.326'],
            'ratio_1p02_certified_margin_lower': ratio_direct_margin.str(50),
            'ratio_1p02_margin_radius': ratio_direct_margin.rad().str(20),
        }

    quadrature = mpmath_quadrature()
    analytic_record = {
        'precision_bits': ANALYTIC_BITS,
        'A': analytic['A'].str(60),
        'B': analytic['B'].str(60),
        'r1': analytic['r1'].str(60),
        'r2': analytic['r2'].str(60),
        'zero_of_p': zero_of_p.str(60),
        'inner_halfline': analytic['inner_halfline'].str(60),
        'outer_halfline': analytic['outer_halfline'].str(60),
        'J_1p30': analytic['J'].str(60),
        'd_1p30': d130.str(60),
        'd_1p45': d145.str(60),
        'equal_c_margin': equal_margin.str(60),
        'coarse_pair_margin': coarse_margin.str(60),
        'evidence': 'MACHINE-VERIFIED Arb closed form',
    }
    result = {
        'status': 'PASS',
        'analytic': analytic_record,
        'mpmath_cross_check': quadrature,
        'root_box_count_n12': 132,
        'c_tile_cover_control': {
            'domain': ['1.30', '1.45'],
            'ratio': '1.02',
            'tiles': len(tiles),
            'adjacency_and_endpoints_asserted': True,
        },
        'direct_envelope': direct_record,
        'controls': {
            'basis_and_active_regions': 'PASS',
            'known_windows': 'PASS',
            'weight_counterfactual': 'PASS',
            'margin_counterfactuals': 'PASS',
            'analytic_containment_and_pseudo_upper_plant': 'PASS',
        },
    }
    out = os.path.join(os.path.dirname(__file__), '..', 'logs', 'startup_controls.json')
    with open(out, 'w', encoding='utf-8') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
