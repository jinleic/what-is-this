#!/usr/bin/env python3
"""Mutation and sanity tests for the centered Liu H2 chart certificate.

Run from the repository root with

    ./.venv/bin/python -I -B uc/verification/test_liu9_chart_centered.py
"""

from __future__ import annotations

import dataclasses
import json
import math
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

from flint import arb  # noqa: E402

from liu9_binding import (  # noqa: E402
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_chart_centered import (  # noqa: E402
    CellEnclosure,
    DomainCertificate,
    HessianBox,
    Jet3,
    RadiusCertificate,
    _canonical_digest,
    build_cell_enclosure,
    candidate_holds,
    centered_entry,
    certify_radius,
    certify_support_domain,
    compose_unary,
    entropy_third_formula,
    finite_difference_check,
    jet_entropy,
    smooth_chart_ceiling,
)

MP = solve_equation_parameters(130)
ARB = certify_equation_parameters(MP)
RADIUS = Fraction(1, 256)
Q_FLOOR = Fraction(1, 4)
REPORT = os.path.join(HERE, "results", "liu9-chart-centered.json")


def _width(value: arb) -> float:
    return float(value.upper()) - float(value.lower())


class ThirdOrderJetControlTest(unittest.TestCase):

    def test_full_jet_matches_independent_finite_differences(self) -> None:
        report = finite_difference_check(MP, ARB, 120)
        self.assertEqual(report["derivative_checks"], 12)
        self.assertTrue(report["passed"])
        self.assertLess(report["worst_relative_error"], 1e-20)

    def test_missing_faa_di_bruno_cross_term_is_rejected(self) -> None:
        """Omitting h''*u_ss*u_d from the ssd partial must be visible."""
        argument = Jet3(
            arb("0.37"),
            gs=arb("0.2"),
            gd=arb("-0.15"),
            hss=arb("0.11"),
            hsd=arb("0.07"),
            hdd=arb("-0.09"),
            tssd=arb("0.03"),
        )
        composed = jet_entropy(argument)
        u = argument.v
        first = (1 - u).log() - u.log()
        second = -1 / u - 1 / (1 - u)
        third = entropy_third_formula(u)
        expected = (
            third * argument.gs**2 * argument.gd
            + second * (
                argument.hss * argument.gd
                + 2 * argument.hsd * argument.gs
            )
            + first * argument.tssd
        )
        mutated = (
            third * argument.gs**2 * argument.gd
            + second * (2 * argument.hsd * argument.gs)
            + first * argument.tssd
        )
        self.assertTrue(composed.tssd.overlaps(expected))
        self.assertFalse(composed.tssd.overlaps(mutated))

    def test_wrong_entropy_third_sign_or_power_is_rejected(self) -> None:
        """An affine argument isolates h''' from every other chain-rule term."""
        u = arb("0.37")
        composed = jet_entropy(Jet3(u, gs=arb(1)))
        correct = (1 - 2 * u) / (u**2 * (1 - u) ** 2)
        wrong_sign = (2 * u - 1) / (u**2 * (1 - u) ** 2)
        wrong_power = (1 - 2 * u) / (u * (1 - u))
        self.assertTrue(composed.tsss.overlaps(correct))
        self.assertFalse(composed.tsss.overlaps(wrong_sign))
        self.assertFalse(composed.tsss.overlaps(wrong_power))

    def test_generic_composition_has_all_mixed_pairings(self) -> None:
        """A polynomial outer map independently checks compose_unary itself."""
        argument = Jet3(
            arb("0.4"),
            gs=arb("0.3"),
            gd=arb("0.2"),
            hss=arb("0.1"),
            hsd=arb("-0.05"),
            hdd=arb("0.07"),
            tssd=arb("0.04"),
            tsdd=arb("-0.03"),
        )
        value = argument.v**3
        composed = compose_unary(
            argument,
            value,
            3 * argument.v**2,
            6 * argument.v,
            arb(6),
        )
        expected_ssd = (
            6 * argument.gs**2 * argument.gd
            + 6 * argument.v * (
                argument.hss * argument.gd
                + 2 * argument.hsd * argument.gs
            )
            + 3 * argument.v**2 * argument.tssd
        )
        expected_sdd = (
            6 * argument.gs * argument.gd**2
            + 6 * argument.v * (
                argument.hdd * argument.gs
                + 2 * argument.hsd * argument.gd
            )
            + 3 * argument.v**2 * argument.tsdd
        )
        self.assertTrue(composed.tssd.overlaps(expected_ssd))
        self.assertTrue(composed.tsdd.overlaps(expected_sdd))

    def test_entropy_endpoint_convention_requires_a_constant_endpoint(self) -> None:
        self.assertEqual(jet_entropy(Jet3.constant(0)).v, 0)
        self.assertEqual(jet_entropy(Jet3.constant(1)).v, 0)
        with self.assertRaises(ArithmeticError):
            jet_entropy(Jet3(arb(0), gs=arb(1)))
        with self.assertRaises(ArithmeticError):
            jet_entropy(Jet3(arb(1), gd=arb(1)))


class CenteredCertificateControlTest(unittest.TestCase):

    def test_omitting_the_radius_factor_is_rejected(self) -> None:
        enclosure = centered_entry(
            arb(1), arb(2), arb(3), Fraction(1, 16)
        )
        honest_width = 2 * (2 + 3) / 16
        mutated_width = 2 * (2 + 3)
        self.assertAlmostEqual(_width(enclosure), honest_width, delta=1e-7)
        self.assertGreater(abs(_width(enclosure) - mutated_width), 9.0)

    def test_kappa_above_the_ceiling_is_rejected_before_the_pencil(self) -> None:
        ceiling = smooth_chart_ceiling(ARB)
        holds, reason, cell = candidate_holds(
            [], Fraction(2, 5), ceiling
        )
        self.assertFalse(holds)
        self.assertEqual(reason, "above_ceiling")
        self.assertIsNone(cell)
        below_holds, _, _ = candidate_holds([], Fraction(3, 8), ceiling)
        self.assertTrue(below_holds)

    def test_support_range_check_is_load_bearing(self) -> None:
        domain = certify_support_domain(ARB, RADIUS, Q_FLOOR)
        self.assertTrue(domain.checked)
        self.assertTrue(domain.safe)
        with self.assertRaises(ValueError):
            build_cell_enclosure(
                ARB,
                RADIUS,
                Q_FLOOR,
                Q_FLOOR + RADIUS,
                None,
            )
        unchecked = dataclasses.replace(domain, checked=False)
        with self.assertRaises(ValueError):
            build_cell_enclosure(
                ARB,
                RADIUS,
                Q_FLOOR,
                Q_FLOOR + RADIUS,
                unchecked,
            )

        claimed_without_check = RadiusCertificate(
            radius=RADIUS,
            q_floor=Q_FLOOR,
            q_cells=1,
            domain=unchecked,
            kappa=Fraction(1, 100),
            kappa_resolution=Fraction(1, 1 << 20),
            ceiling=smooth_chart_ceiling(ARB),
            centered_hss_width=0.1,
            naive_hss_width=1.0,
            centered_inflation=25.6,
            naive_inflation=256.0,
            naive_failure=None,
            third_subdivisions=1,
            center_containment_checks=6,
            full_overlap_checks=6,
            next_failure=None,
            next_failing_q_cell=None,
        )
        self.assertFalse(claimed_without_check.certified())

    def test_centered_hessian_contains_the_independent_naive_center(self) -> None:
        domain = certify_support_domain(ARB, RADIUS, Q_FLOOR)
        cell = build_cell_enclosure(
            ARB,
            RADIUS,
            Q_FLOOR,
            Q_FLOOR + RADIUS,
            domain,
            third_subdivisions=2,
        )
        self.assertEqual(cell.center_containment_checks, 6)
        self.assertEqual(cell.full_overlap_checks, 6)
        self.assertIsNotNone(cell.naive_gap)

    def test_requested_transition_radius_and_ceiling_gate(self) -> None:
        failing = certify_radius(
            ARB,
            Fraction(1, 128),
            bisection_steps=12,
            third_subdivisions=2,
        )
        passing = certify_radius(
            ARB,
            RADIUS,
            bisection_steps=12,
            third_subdivisions=2,
        )
        self.assertFalse(failing.certified())
        self.assertTrue(passing.certified())
        self.assertGreater(passing.kappa, 0)
        self.assertLessEqual(
            float(passing.kappa), float(passing.ceiling.lower())
        )
        self.assertTrue(passing.domain.checked and passing.domain.safe)

    def test_interval_psd_control_rejects_a_large_off_diagonal(self) -> None:
        gap = HessianBox(arb(1), arb(2), arb(1))
        distance = HessianBox(arb(0), arb(0), arb(0))
        dummy = CellEnclosure(
            Q_FLOOR,
            Q_FLOOR + RADIUS,
            gap,
            distance,
            None,
            distance,
            "not needed",
            0,
            0,
        )
        holds, reason, _ = candidate_holds(
            [dummy], Fraction(0), smooth_chart_ceiling(ARB)
        )
        self.assertFalse(holds)
        self.assertEqual(reason, "pencil")

    def test_committed_report_has_a_canonical_digest_and_finite_numbers(self) -> None:
        with open(REPORT, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
        reported = payload.pop("report_sha256")
        self.assertEqual(reported, _canonical_digest(payload))

        def check_finite(value: object) -> None:
            if isinstance(value, float):
                self.assertTrue(math.isfinite(value))
            elif isinstance(value, dict):
                for child in value.values():
                    check_finite(child)
            elif isinstance(value, list):
                for child in value:
                    check_finite(child)

        check_finite(payload)
        largest = payload["largest_positive_certificate"]
        self.assertEqual(largest["radius"], "1/256")
        self.assertGreater(largest["kappa_float"], 0)
        self.assertTrue(payload["sanity_gates"]["all_kappas_below_ceiling"])
        self.assertTrue(
            payload["sanity_gates"]["all_certified_radii_checked_support"]
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
