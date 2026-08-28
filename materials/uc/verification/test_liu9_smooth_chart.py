#!/usr/bin/env python3
"""Mutation tests for the Liu H2 smooth-chart certificate.

A control that cannot fail proves nothing.  Each test below rejects one
specific way the chart lemma could be wrong, and several of them reject a
mistake that was actually made while writing it.

Run directly (this file is not importable as a package module):
    math/.venv/bin/python -I -B math/uc/verification/test_liu9_smooth_chart.py
"""

from __future__ import annotations

import os
import sys
import unittest
from fractions import Fraction

for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

HERE = os.path.dirname(os.path.abspath(__file__))
UC = os.path.dirname(HERE)
if UC not in sys.path:
    sys.path.insert(0, UC)

from flint import arb, ctx  # noqa: E402

import liu9_smooth_chart as chart  # noqa: E402
from liu9_binding import (  # noqa: E402
    certify_curvatures,
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_boundary_layer import entropy_at_star  # noqa: E402

ctx.prec = max(ctx.prec, 320)

MP = solve_equation_parameters(90)
PARAMS = certify_equation_parameters(MP)
CEILING = chart.smooth_chart_ceiling(PARAMS)


class SmoothChartControlTest(unittest.TestCase):

    def test_centre_is_a_critical_point(self):
        """The mean-value form is void without grad = 0 at the centre."""
        report = chart.critical_point_check(
            PARAMS, (Fraction(0), Fraction(1, 3), Fraction(1)))
        self.assertLess(report["worst_magnitude"], 1e-60)

    def test_hessian_matches_the_published_curvatures(self):
        """The jet must agree with liu9_binding's independent certificate."""
        report = chart.hessian_cross_check(
            PARAMS, (Fraction(1, 4), Fraction(1, 2)))
        self.assertLess(report["worst_difference"], 1e-40)

    def test_dropping_the_ehpi_component_weight_is_caught(self):
        """EHPI carries (1-q) and q.  Omitting them moves gap(P*) off zero.

        This is not hypothetical: the first draft of the chart omitted exactly
        this factor, and the symptom was gap(P*) = 0.0553 instead of 0.  The
        test rebuilds the wrong functional and asserts the critical-point check
        would have rejected it.
        """
        parameters = PARAMS
        x = parameters.x
        s = chart.Jet2.variable(x, "s")
        d = chart.Jet2.variable(arb(0), "d")
        q = arb(0)
        weights, supports, masses = chart.chart(parameters, q, s, d, arb(0))
        unweighted = chart.Jet2.constant(0)
        for component in range(2):
            base = 3 * component
            for i in range(3):
                for k in range(3):
                    y, z = supports[base + i], supports[base + k]
                    one = chart.Jet2.constant(1)
                    protocol = y * z * (one + (one - y) * (one - z))
                    unweighted = unweighted + (masses[i] * masses[k]
                                               * chart.jet_entropy(protocol))
        correct = chart.gap_jet(parameters, q, s, d, arb(0))
        ehxy = chart.Jet2.constant(0)
        for i in range(6):
            for k in range(6):
                ehxy = ehxy + weights[i] * weights[k] * chart.jet_entropy(
                    supports[i] * supports[k])
        ehx = chart.Jet2.constant(0)
        for i in range(6):
            ehx = ehx + weights[i] * chart.jet_entropy(supports[i])
        beta = parameters.beta
        mutated = (chart.Jet2.constant(1 - beta) * ehxy
                   + chart.Jet2.constant(beta) * unweighted
                   - ehx)
        self.assertTrue(correct.v.contains(0),
                        "the correct gap must vanish at P*")
        self.assertFalse(mutated.v.contains(0),
                         "the unweighted EHPI must be detectably wrong")
        self.assertGreater(float(mutated.v.lower()), 0.05)

    def test_pencil_beats_the_eigenvalue_ratio(self):
        """lambda_min/lambda_max is the wrong comparison and must be rejected.

        At small q(1-q) the gap Hessian's dd entry is O(q(1-q)) while the
        distance Hessian's ss entry is not, so the ratio of separate extreme
        eigenvalues collapses to zero even though the pencil is fine.  The
        generalized comparison is the one that respects the cancellation.
        """
        x = PARAMS.x
        h_star = entropy_at_star(PARAMS)
        curvatures = certify_curvatures(PARAMS)
        q = arb(1) / 1000
        weight = q * (1 - q)
        hss = h_star * curvatures.face
        hdd = h_star * curvatures.split_unit * weight
        dss = 2 * PARAMS.p * x**2
        ddd = 2 * (PARAMS.p**2 + PARAMS.p * x**2) * weight
        zero = arb(0)
        ratio = hdd / dss
        self.assertLess(float(ratio.upper()), float(CEILING.lower()) / 100,
                        "the naive ratio should collapse to a tiny fraction "
                        "of the ceiling at small q")
        self.assertTrue(
            chart._pencil_is_psd(hss, zero, hdd, dss, zero, ddd,
                                 arb(3) / 10),
            "the pencil must still certify kappa=0.3 at q=1/1000")

    def test_kappa_never_exceeds_the_ceiling(self):
        """No kappa above H*min(A/(2px^2), C/(2(p^2+px^2))) can hold."""
        certificate = chart.certify_chart(PARAMS, Fraction(1, 2048),
                                          Fraction(1, 4), q_cells=8)
        self.assertGreater(certificate.kappa, 0)
        self.assertLessEqual(float(certificate.kappa),
                             float(CEILING.upper()))

    def test_refuted_radius_is_load_bearing(self):
        """1/2048 certifies and 1/1024 does not; if both pass, the grid lies."""
        good = chart.certify_chart(PARAMS, Fraction(1, 2048), Fraction(1, 4),
                                   q_cells=8)
        bad = chart.certify_chart(PARAMS, Fraction(1, 1024), Fraction(1, 4),
                                  q_cells=8)
        self.assertTrue(good.certified())
        self.assertFalse(bad.certified(),
                         "radius 1/1024 must not certify, or the naive "
                         "enclosure has silently improved and the recorded "
                         "obstruction constant is stale")

    def test_entropy_endpoint_convention_is_not_applied_to_small_arguments(self):
        """h(0)=0 holds only for the identically constant zero jet.

        Applying it to a merely small argument would discard the entropy
        singularity that the boundary-layer modules exist to handle.
        """
        constant_zero = chart.Jet2.constant(0)
        self.assertTrue(chart.jet_entropy(constant_zero).is_constant_zero())
        tiny = chart.Jet2.variable(arb(10) ** -30, "s")
        value = chart.jet_entropy(tiny)
        self.assertFalse(value.is_constant_zero())
        self.assertGreater(float(value.gs.lower()), 60.0,
                           "h'(1e-30) is about 69, not 0")

    def test_r_is_a_gauge(self):
        """Fixing r=0 must be a quotient, not a restriction."""
        report = chart.gauge_check(PARAMS, Fraction(1, 2048),
                                   (Fraction(0), Fraction(1, 50)))
        self.assertLess(report["worst_difference"], 1e-60)

    def test_enclosure_width_still_inflates_linearly(self):
        """The recorded obstruction constant must stay honest.

        If a future change makes the enclosure tighter, this test fails and
        forces the documented constant and cell estimate to be updated rather
        than left stale.
        """
        coarse = chart.certify_chart(PARAMS, Fraction(1, 512), Fraction(1, 4),
                                     q_cells=8)
        fine = chart.certify_chart(PARAMS, Fraction(1, 2048), Fraction(1, 4),
                                   q_cells=8)
        self.assertGreater(coarse.hessian_width, fine.hessian_width)
        constant = fine.hessian_width * 2048
        self.assertGreater(constant, 400.0)
        self.assertLess(constant, 2000.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
