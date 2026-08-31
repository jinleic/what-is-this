#!/usr/bin/env python3
"""Mutation tests for the endpoint support theorem and exact inward-q lift."""

from __future__ import annotations

import math
import os
import sys
import unittest
from fractions import Fraction
from pathlib import Path

for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

import mpmath
from flint import arb

HERE = Path(__file__).resolve().parent
UC = HERE.parent
if str(UC) not in sys.path:
    sys.path.insert(0, str(UC))

from liu9_binding import (  # noqa: E402
    certify_equation_parameters, solve_equation_parameters,
)
from liu9_boundary_layer import gap_mp  # noqa: E402
from liu9_objective import arb_fraction  # noqa: E402
from liu9_qendpoint_lift import (  # noqa: E402
    _bisect_coordinate, _full_bound, _quadratic_remainder_lower,
    _refined_zero_mass_penalty, _single_mass_high_endpoint_lower,
    _swap_components, q_gap_coefficients, qendpoint_shift_lower,
)
from liu9_qone_face import (  # noqa: E402
    _cross_slope_ceiling, best_qone_endpoint_lower, canonical_qone_box,
    contract_qone_mean, qone_feasible_shift_lower, reduce_simplex_face,
    zero_mass_face, zero_mass_total_upper,
)

MP = solve_equation_parameters(100)
PARAMETERS = certify_equation_parameters(MP)


class EndpointSupportTest(unittest.TestCase):
    def setUp(self):
        self.blocker = canonical_qone_box((
            (13 / 32, 7 / 16),
            (0.0, 1 / 32),
            (1 / 32, 1 / 16),
            (13 / 32, 7 / 16),
            (15 / 16, 1.0),
        ))

    def test_exact_blocker_is_certified_after_mean_contraction(self):
        contracted, count = contract_qone_mean(
            self.blocker, PARAMETERS.mean)
        self.assertIsNotNone(contracted)
        self.assertGreaterEqual(count, 3)
        self.assertGreater(contracted[6][0], 0.057)
        self.assertGreater(contracted[8][0], 0.996)
        lower, method = best_qone_endpoint_lower(
            contracted, PARAMETERS.mean, PARAMETERS.beta)
        self.assertGreater(lower, 0)
        self.assertIn("endpoint", method)

    def test_skipping_mean_contraction_is_rejected_by_sign(self):
        lower, _ = best_qone_endpoint_lower(
            self.blocker, PARAMETERS.mean, PARAMETERS.beta)
        self.assertLess(lower, 0)

    def test_contractor_retains_binding_mean_corner(self):
        contracted, _ = contract_qone_mean(self.blocker, PARAMETERS.mean)
        with mpmath.workdps(80):
            a1 = mpmath.mpf(13) / 32
            b1 = mpmath.mpf(1) / 16
            a3 = 1 - a1
            b5 = (MP.mean - a1 * b1) / a3
        self.assertLessEqual(contracted[0][0], float(a1))
        self.assertGreaterEqual(contracted[0][1], float(a1))
        self.assertLessEqual(contracted[6][0], float(b1))
        self.assertGreaterEqual(contracted[6][1], float(b1))
        self.assertLessEqual(contracted[8][0], float(b5))
        self.assertGreaterEqual(contracted[8][1], float(b5))


class ExactQPolynomialTest(unittest.TestCase):
    def setUp(self):
        self.values = (
            mpmath.mpf("0.40625"), mpmath.mpf("0.001"),
            mpmath.mpf("0.97"),
            mpmath.mpf("0.2"), mpmath.mpf("0.4"), mpmath.mpf("0.8"),
            mpmath.mpf("0.0625"), mpmath.mpf("0.42"), mpmath.mpf("0.999"),
        )
        self.box = tuple((float(value), float(value)) for value in self.values)

    def test_quadratic_reconstructs_sanctioned_raw_gap(self):
        g0, g1, g2, _, _ = q_gap_coefficients(
            self.box, PARAMETERS.beta)
        with mpmath.workdps(80):
            exact = gap_mp(self.values, MP.beta)
            radius = 1 - self.values[2]
        polynomial = g0 + arb(float(radius)) * g1 + arb(float(radius**2)) * g2
        self.assertLessEqual(float(polynomial.lower()), float(exact))
        self.assertGreaterEqual(float(polynomial.upper()), float(exact))

    def test_dropping_linear_coefficient_is_detected(self):
        g0, _, g2, _, _ = q_gap_coefficients(self.box, PARAMETERS.beta)
        with mpmath.workdps(80):
            exact = gap_mp(self.values, MP.beta)
            radius = 1 - self.values[2]
        mutated = g0 + arb(float(radius**2)) * g2
        self.assertGreater(abs(float(mutated.mid()) - float(exact)), 1e-3)

    def test_active_mean_bound_handles_m0_greater_than_m1(self):
        values = (
            mpmath.mpf(0), mpmath.mpf(0), mpmath.mpf("0.99"),
            mpmath.mpf("0.5"), mpmath.mpf("0.5"), mpmath.mpf("0.9"),
            mpmath.mpf("0.5"), mpmath.mpf("0.5"), mpmath.mpf("0.62"),
        )
        box = tuple((float(value), float(value)) for value in values)
        g0, g1, g2, mean0, mean1 = q_gap_coefficients(
            box, PARAMETERS.beta)
        self.assertGreater(mean0.lower(), mean1.upper())
        q_lo = values[2]
        active_target = 1 - (1 - PARAMETERS.mean) / arb(float(q_lo))
        active = canonical_qone_box(
            tuple(box[index] for index in (0, 1, 6, 7, 8)))
        lam = Fraction(-1, 32)
        base, _ = qone_feasible_shift_lower(
            active, active_target, PARAMETERS.beta, lam)
        radius = 1 - q_lo
        candidate = base + arb(float(radius)) * g1 + arb(float(radius**2)) * g2
        with mpmath.workdps(80):
            exact = gap_mp(values, MP.beta)
            mean = (1 - q_lo) * values[5] + q_lo * values[8]
        self.assertGreaterEqual(mean, MP.mean)
        self.assertLessEqual(float(candidate.lower()), float(exact))
    def test_common_radius_avoids_decorrelation_mutation(self):
        lower = _quadratic_remainder_lower(
            0.6, 0.8, arb(1), arb(-1))
        # r-r^2 is increasing on [0.2,0.4], so the exact lower is 0.16.
        self.assertGreaterEqual(float(lower), 0.16 - 1e-15)
        decorrelated = (
            arb(0.2).union(arb(0.4)) * arb(1)
            + arb(0.2).union(arb(0.4))**2 * arb(-1)
        ).lower()
        self.assertGreater(float(lower), float(decorrelated))

    def test_convex_q_remainder_checks_interior_vertex(self):
        lower = _quadratic_remainder_lower(
            0.2, 1.0, arb(-1), arb(1))
        self.assertLessEqual(float(lower), -0.25)
        self.assertGreater(float(lower), -0.250000000000001)
    def test_component_swap_is_exact_symmetry(self):
        swapped_box = _swap_components(self.box)
        exact_q = 1 - self.values[2]
        self.assertLessEqual(swapped_box[2][0], float(exact_q))
        self.assertGreaterEqual(swapped_box[2][1], float(exact_q))
        swapped = (
            self.values[0], self.values[1], exact_q,
            *self.values[6:9], *self.values[3:6],
        )
        with mpmath.workdps(80):
            original_gap = gap_mp(self.values, MP.beta)
            swapped_gap = gap_mp(swapped, MP.beta)
        self.assertLess(abs(original_gap - swapped_gap), mpmath.mpf("1e-70"))

    def test_positive_feasibility_multiplier_is_rejected(self):
        with self.assertRaises(ValueError):
            qendpoint_shift_lower(
                self.box, PARAMETERS.mean, PARAMETERS.beta, Fraction(1, 100))

    def test_q_zero_layer_requires_separate_endpoint(self):
        mutated = list(self.box)
        mutated[2] = (0.0, self.box[2][1])
        with self.assertRaises(ValueError):
            qendpoint_shift_lower(
                tuple(mutated), PARAMETERS.mean, PARAMETERS.beta, Fraction(0))


class ZeroMassLayerTest(unittest.TestCase):
    def test_factored_face_closes_previous_deep_blocker(self):
        box = (
            (0.0, 0.00390625), (0.0, 0.0078125), (0.5, 0.5078125),
            (0.0, 0.03125), (0.46875, 0.5), (0.96875, 1.0),
            (0.0078125, 0.015625), (0.5, 0.5078125),
            (0.234375, 0.2421875),
        )
        lower, _, method, _ = _full_bound(box, PARAMETERS)
        self.assertGreater(lower, 0)
        self.assertTrue(method.startswith("mass-zero"))

    def test_omitting_mass_loss_would_false_clear_wide_layer(self):
        box = (
            (0.0, 0.03), (0.0, 0.03), (0.5, 0.5078125),
            (0.0, 0.03125), (0.46875, 0.5), (0.96875, 1.0),
            (0.0078125, 0.015625), (0.5, 0.5078125),
            (0.234375, 0.2421875),
        )
        face = zero_mass_face(box)
        face_lower = _single_mass_high_endpoint_lower(
            face, PARAMETERS.mean, PARAMETERS.beta, Fraction(0))
        penalty = _refined_zero_mass_penalty(
            box, face, PARAMETERS.beta, Fraction(0))
        self.assertGreater(face_lower, 0)
        self.assertLess(face_lower - penalty, 0)


class OutwardRoundingTest(unittest.TestCase):
    def test_simplex_contractor_retains_halfway_exact_point(self):
        tiny = Fraction(3, 2**55)
        exact_a1 = 1 - tiny
        nearest = float(exact_a1)
        a1_box = (
            math.nextafter(nearest, -math.inf),
            math.nextafter(nearest, math.inf),
        )
        a2 = float(tiny)
        box = canonical_qone_box((
            a1_box, (a2, a2), (1.0, 1.0), (1.0, 1.0), (1.0, 1.0),
        ))
        contracted, _ = contract_qone_mean(box, PARAMETERS.mean)
        self.assertIsNotNone(contracted)
        self.assertGreaterEqual(
            Fraction.from_float(contracted[0][1]), exact_a1)

    def test_simplex_face_does_not_round_near_equality_to_equal(self):
        below_half = math.nextafter(0.5, 0.0)
        box = canonical_qone_box((
            (0.5, 0.75), (below_half, 0.5),
            (0.25, 0.5), (0.25, 0.5), (0.5, 1.0),
        ))
        _, reduced = reduce_simplex_face(box)
        self.assertFalse(reduced)

    def test_q_remainder_uses_exact_one_minus_q(self):
        q = Fraction(3, 2**55)
        lower = _quadratic_remainder_lower(
            float(q), float(q), arb(-1), arb(0))
        expected = -arb_fraction(1 - q)
        self.assertLessEqual(lower, expected)

    def test_cross_slope_float_product_is_outward(self):
        y = math.nextafter(1.0, 0.0)
        x = 3.0 / 16.0
        ceiling = _cross_slope_ceiling(x, x, y)
        with mpmath.workdps(80):
            z = mpmath.mpf(y) * mpmath.mpf(x)
            actual = mpmath.mpf(x) * mpmath.log((1 - z) / z)
        self.assertGreaterEqual(float(ceiling.upper()), float(actual))

    def test_zero_mass_halfway_sum_is_not_rounded_down(self):
        upper = math.nextafter(1.0 / 32.0, math.inf)
        box = (
            (0.0, 1.0 / 32.0), (0.0, upper), (1.0, 1.0),
            (0.5, 0.5), (0.5, 0.5), (0.5, 0.5),
            (0.5, 0.5), (0.5, 0.5), (0.5, 0.5),
        )
        self.assertGreater(zero_mass_total_upper(box).upper(), arb(1) / 16)

    def test_unsplittable_preferred_coordinate_falls_back(self):
        lo = 0.6785225868225097
        adjacent = math.nextafter(lo, math.inf)
        box = (
            (lo, adjacent), (0.0, 0.0), (0.75, 1.0),
            (0.0, 0.0), (0.0, 0.0), (0.0, 0.0),
            (0.0, 0.0), (0.0, 0.0), (0.0, 0.0),
        )
        left, right = _bisect_coordinate(box, 0)
        self.assertEqual(left[0], box[0])
        self.assertEqual(right[0], box[0])
        self.assertEqual(left[2][1], right[2][0])
        self.assertLess(left[2][0], left[2][1])
        self.assertLess(right[2][0], right[2][1])
if __name__ == "__main__":
    unittest.main(verbosity=2)
