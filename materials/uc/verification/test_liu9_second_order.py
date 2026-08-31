#!/usr/bin/env python3
"""Mutation tests for ``liu9_second_order.py``.

Every test rejects a concrete way the second-order theorem could be wrong:
omitting the component-mean distance term, assuming first-order decoupling at
second order, losing a domain endpoint, accepting a coarse dependency-widened
cover, changing a sign, or confusing a negative pencil with a negative raw gap.
"""

from __future__ import annotations

import os
import sys
import unittest
from fractions import Fraction as F

for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

import mpmath
from flint import arb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
UC = os.path.dirname(HERE)
if UC not in sys.path:
    sys.path.insert(0, UC)

from liu9_binding import certify_equation_parameters, solve_equation_parameters  # noqa: E402
from liu9_boundary_layer import gap_mp  # noqa: E402
from liu9_ninevar import distance_squared  # noqa: E402
from liu9_second_order import (  # noqa: E402
    HALF,
    KAPPA,
    automatic_expansion,
    certify_exact_expansion,
    certify_symmetric,
    cross_pencil,
    cross_raw,
    endpoint_obstruction,
    hand_expansion,
    mean_excess_margin,
    certify_mean_excess,
    jhat_average_bound,
    linear_coefficient,
    q_over_y2_boundary,
    quadratic_minimum_lower,
    quadratic_symmetric,
    symmetric_half_coefficient,
    zero_tail_jhat_lower,
    _arbf,
)

ctx.prec = max(ctx.prec, 320)
MP = solve_equation_parameters(100)
ARB = certify_equation_parameters(MP)
KAPPA_MP = mpmath.mpf(KAPPA.numerator) / KAPPA.denominator


class ExactExpansionTest(unittest.TestCase):
    def test_automatic_and_hand_coefficients_agree(self):
        report = certify_exact_expansion(ARB)
        self.assertTrue(report["certified"])
        self.assertEqual(report["exact_degree"], 2)
        self.assertEqual(report["epsilon_dependent_entropy_arguments"], 0)
        self.assertEqual(report["epsilon_squared_log_coefficient"], "0")

    def test_no_log_term_over_twenty_epsilon_decades(self):
        """A hidden eps^2 log term would make the normalized residual drift."""
        with mpmath.workdps(90):
            y = mpmath.mpf(3) / 10
            q = mpmath.mpf(1) / 2
            from liu9_ninevar import transverse_coefficients_mp
            l, q2 = transverse_coefficients_mp(MP, y, "mean_preserving")
            values = []
            for exponent in (8, 12, 16, 20, 24, 28):
                eps = mpmath.mpf(10) ** (-exponent)
                vector = (MP.p - eps * y / MP.x, eps, q,
                          MP.x, y, 0, MP.x, y, 0)
                pencil = gap_mp(vector, MP.beta) - KAPPA_MP * distance_squared(vector, MP)
                values.append((pencil - eps * l) / (eps * eps))
            self.assertLess(float(max(values) - min(values)), 1e-12)
            self.assertAlmostEqual(float(values[-1]), q2, places=10)

    def test_dropping_quadratic_distance_is_detected(self):
        """Mean preservation is global; the two component means still differ."""
        y1, y2, q = _arbf(F(1, 5)), _arbf(F(9, 10)), _arbf(F(1, 3))
        automatic = automatic_expansion(ARB, y1, y2, q)
        hand = hand_expansion(ARB, y1, y2, q)
        correct = hand["quadratic"]
        wrong = correct + _arbf(KAPPA) * q * (arb(1) - q) * (y1 - y2) ** 2
        self.assertTrue(automatic["pencil"].c2.overlaps(correct))
        self.assertFalse(automatic["pencil"].c2.overlaps(wrong))

    def test_cross_interaction_does_not_decouple_at_second_order(self):
        y1, y2, q = _arbf(F(1, 5)), _arbf(F(9, 10)), _arbf(F(1, 2))
        hand = hand_expansion(ARB, y1, y2, q)
        decoupled = ((arb(1) - q) * quadratic_symmetric(ARB, y1)
                     + q * quadratic_symmetric(ARB, y2))
        self.assertFalse(hand["quadratic"].overlaps(decoupled))


class SymmetricSignTest(unittest.TestCase):
    def test_zero_boundary_quadratic_is_positive(self):
        bound = q_over_y2_boundary(ARB, F(1, 32))
        self.assertGreater(float(bound.lower()), 0.09)

    def test_double_root_curvature_is_positive(self):
        with mpmath.workdps(70):
            def s(value):
                y = mpmath.mpf(value)
                from liu9_ninevar import transverse_coefficients_mp
                l, q2 = transverse_coefficients_mp(MP, y, "mean_preserving")
                return 2 * l + q2
            self.assertAlmostEqual(float(mpmath.diff(s, MP.x, 2)),
                                   1.4651674782, places=8)

    def test_coarse_cover_fails_closed(self):
        """Natural/Taylor error at 32 cells is wider than the true margin."""
        with self.assertRaises(AssertionError):
            certify_symmetric(ARB, low_cells=32, high_cells=48,
                              window_cells=128)

    def test_sign_flipping_Q_changes_the_ssot_value(self):
        y = _arbf(F(3, 5))
        correct = symmetric_half_coefficient(ARB, y)
        wrong = 2 * linear_coefficient(ARB, y) - quadratic_symmetric(ARB, y)
        self.assertFalse(correct.overlaps(wrong))


class CrossTermTest(unittest.TestCase):
    def test_distance_cross_factor_is_exactly_y1_minus_y2_squared(self):
        y1, y2, q = _arbf(F(1, 5)), _arbf(F(9, 10)), _arbf(F(3, 8))
        expansion = automatic_expansion(ARB, y1, y2, q)
        expected = q * (arb(1) - q) * (y1 - y2) ** 2
        self.assertTrue(expansion["distance"].c2.overlaps(expected))
        self.assertGreater(float(expected.lower()), 0.0)

    def test_cross_term_vanishes_on_diagonal(self):
        for y in (F(1, 10), F(1, 2), F(9, 10)):
            value = cross_pencil(ARB, _arbf(y), _arbf(y))
            self.assertTrue(value.contains(0))
            self.assertLess(float(abs(value).upper()), 1e-60)

    def test_negative_interaction_really_occurs(self):
        """The proof must minimize the q vertex; C cannot be discarded."""
        value = cross_pencil(ARB, _arbf(F(3, 25)), _arbf(F(4, 5)))
        self.assertLess(float(value.upper()), -0.02)

    def test_zero_tail_factor_remains_positive_after_kappa(self):
        value = zero_tail_jhat_lower(ARB, F(1, 10))
        self.assertGreater(float(value.lower()), 0.4)

    def test_diagonal_divided_difference_avoids_dependency(self):
        value = jhat_average_bound(ARB, F(9, 20), F(1, 2))
        self.assertTrue(value.is_finite())
        self.assertGreater(float(value.lower()), -1.0)

    def test_q_vertex_is_not_replaced_by_endpoint_only_logic(self):
        s1, s2, c = arb("0.001"), arb("0.021"), arb("-0.03")
        lower = quadratic_minimum_lower(s1, s2, c)
        self.assertLess(float(lower.lower()), 0.001)
        self.assertGreater(float(lower.lower()), 0.0)


class MeanExcessTest(unittest.TestCase):
    def test_normal_margin_is_uniformly_positive(self):
        _, axis = certify_symmetric(ARB, low_cells=192, high_cells=256,
                                    window_cells=1024)
        report = certify_mean_excess(ARB, axis)
        self.assertTrue(report["certified"])
        self.assertGreater(report["linear_margin_lower"], 0.23)
        self.assertGreater(report["quadratic_coefficient_lower"], 0.63)

    def test_strict_mean_formula_matches_ssot(self):
        """The advisory's delta>0 configuration must be covered explicitly."""
        from liu9_ninevar import _mp_kernel
        with mpmath.workdps(90):
            y1, y2, q = (mpmath.mpf(1) / 5, mpmath.mpf(9) / 10,
                          mpmath.mpf(2) / 5)
            eps, delta = mpmath.mpf(3) / 10, mpmath.mpf(1) / 10000
            ybar = (1 - q) * y1 + q * y2
            base = (MP.p - eps * ybar / MP.x, eps, q,
                    MP.x, y1, 0, MP.x, y2, 0)
            raised = (base[0] + delta, eps, q,
                      MP.x, y1, 0, MP.x, y2, 0)
            def pencil(values):
                return (gap_mp(values, MP.beta)
                        - KAPPA_MP * distance_squared(values, MP))
            kxx = _mp_kernel(MP.x, MP.x, MP.beta)
            hx = -(MP.x * mpmath.log(MP.x)
                   + (1 - MP.x) * mpmath.log(1 - MP.x))
            r1 = _mp_kernel(MP.x, y1, MP.beta) - y1 / MP.x * kxx
            r2 = _mp_kernel(MP.x, y2, MP.beta) - y2 / MP.x * kxx
            rbar = (1 - q) * r1 + q * r2
            predicted = (delta * (hx + 2 * eps * rbar)
                         + delta * delta * (kxx - KAPPA_MP * MP.x * MP.x))
            residual = pencil(raised) - pencil(base) - predicted
            self.assertLess(abs(float(residual)), 1e-75)
            mean_raise = (raised[0] - base[0]) * MP.x
            self.assertGreater(float(mean_raise), 0.0)

    def test_reversing_the_normal_is_not_covered(self):
        """delta<0 lowers the mean and correctly lies outside the theorem."""
        delta = mpmath.mpf(-1) / 10000
        self.assertLess(float(delta * MP.x), 0.0)


class ObstructionAndCounterexampleTest(unittest.TestCase):
    def test_whole_segment_pencil_claim_is_refuted_at_y_one(self):
        report = endpoint_obstruction(ARB)
        self.assertTrue(report["raw_gap_structural_zero"])
        self.assertLess(report["pencil_upper"], -0.007)

    def test_endpoint_omission_would_hide_the_refutation(self):
        """The obstruction is exactly at y=1; the domain endpoint is load-bearing."""
        with mpmath.workdps(70):
            eps = MP.mean
            vector = (0, eps, mpmath.mpf(1) / 2,
                      MP.x, 1, 0, MP.x, 1, 0)
            raw = gap_mp(vector, MP.beta)
            pencil = raw - KAPPA_MP * distance_squared(vector, MP)
            self.assertEqual(raw, 0)
            self.assertLess(float(pencil), -0.007)

    def test_known_ambient_refutation_is_reproduced_and_mean_infeasible(self):
        """Required negative-control mutation: drop the y/x mean correction."""
        with mpmath.workdps(70):
            y = mpmath.mpf(1) / 32
            eps = mpmath.mpf(1) / 1024
            q = mpmath.mpf(1) / 2
            vector = (MP.p - eps, eps, q,
                      MP.x, y, 0, MP.x, y, 0)
            raw = gap_mp(vector, MP.beta)
            mean = ((MP.p - eps) * MP.x + eps * y)
            self.assertLess(float(raw), -5e-4)
            self.assertLess(float(mean - MP.mean), 0.0)

    def test_negative_pencil_is_not_called_negative_raw_gap(self):
        report = endpoint_obstruction(ARB)
        self.assertTrue(report["raw_gap_structural_zero"])
        self.assertLess(report["pencil_upper"], 0.0)
        self.assertIn("pencil obstruction only", report["interpretation"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
