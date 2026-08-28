#!/usr/bin/env python3
"""Mutation tests for the Liu H2 zero-support boundary-layer estimate.

A control that cannot fail proves nothing.  These tests check that
``uc/liu9_boundary_layer.py`` rejects each way the estimate could be wrong:
a threshold that is too large, a mis-differentiated closed form, a discarded
term that is not actually signed, and an inflated constant.

Run from the repository root with

    math/.venv/bin/python -I -B math/uc/verification/test_liu9_boundary_layer.py
"""

from __future__ import annotations

import os
import sys
import unittest
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
UC = os.path.dirname(HERE)
if UC not in sys.path:
    sys.path.insert(0, UC)

import mpmath  # noqa: E402

from liu9_boundary_layer import (  # noqa: E402
    _hp_mp,
    _pi_first_mp,
    _pi_mp,
    _weights,
    autodiff_partial,
    binding_coefficient_check,
    certify_constants,
    closed_form_factor,
    closed_form_partial,
    constants_are_certified,
    gap_mp,
    sweep_paired_link,
)
from liu9_binding import (  # noqa: E402
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_objective import evaluate_mpmath  # noqa: E402

MP = solve_equation_parameters(120)
ARB = certify_equation_parameters(MP)


def _sample_points() -> list[tuple[tuple[mpmath.mpf, ...], int]]:
    mp = mpmath.mpf
    return [
        ((MP.p, 1 - MP.p, mp(0), MP.x, mp(1) / 64, mp(0), mp(0), mp(0), mp(0)), 1),
        ((mp(4) / 10, mp(35) / 100, mp(3) / 10, mp(9) / 10, mp(1) / 50,
          mp(8) / 10, mp(7) / 10, mp(1) / 1000, mp(1) / 2), 1),
        ((mp(1) / 2, mp(1) / 4, mp(7) / 10, mp(1), mp(1), mp(1) / 100,
          mp(95) / 100, mp(1) / 300, mp(1)), 4),
        ((mp(9) / 10, mp(1) / 10, mp(1) / 1000000, mp(8) / 10, mp(1) / 40,
          mp(0), mp(0), mp(0), mp(1)), 1),
    ]


class BoundaryLayerControlTest(unittest.TestCase):

    def test_threshold_is_load_bearing(self) -> None:
        """1/8 and 1/16 must be rejected; 1/24 and below must be accepted."""
        self.assertFalse(constants_are_certified(
            certify_constants(ARB, Fraction(1, 8))))
        self.assertFalse(constants_are_certified(
            certify_constants(ARB, Fraction(1, 16))))
        for denominator in (24, 32, 64, 128, 1024):
            self.assertTrue(
                constants_are_certified(
                    certify_constants(ARB, Fraction(1, denominator))),
                f"y0=1/{denominator} should be certifiable")

    def test_constants_reject_an_out_of_range_threshold(self) -> None:
        with self.assertRaises(ValueError):
            certify_constants(ARB, Fraction(1, 3))
        with self.assertRaises(ValueError):
            certify_constants(ARB, Fraction(0))

    def test_closed_form_matches_forward_mode(self) -> None:
        with mpmath.workdps(80):
            for values, index in _sample_points():
                closed = closed_form_partial(values, index, MP.beta)
                auto = autodiff_partial(values, index, MP.beta)
                scale = max(abs(closed), abs(auto), mpmath.mpf(1))
                self.assertLess(float(abs(closed - auto) / scale), 1e-60)

    def test_gap_agrees_with_the_published_transcription(self) -> None:
        """_GapOps must only skip the quotient, never change N or D."""
        with mpmath.workdps(60):
            values = _sample_points()[1][0]
            terms = evaluate_mpmath(values, MP.beta, dps=60)
            reference = terms.numerator - terms.ehx
            self.assertLess(
                float(abs(gap_mp(values, MP.beta) - reference)), 1e-45)

    def test_dropping_the_entropy_term_is_caught(self) -> None:
        """Removing -h'(y) shifts the leading coefficient by exactly one."""
        with mpmath.workdps(80):
            for values, index in _sample_points():
                correct = closed_form_factor(values, index, MP.beta)
                y = values[3 + index]
                mutated = correct + _hp_mp(y)
                auto = autodiff_partial(values, index, MP.beta)
                weights = _weights(values, mpmath.mpf(1))
                self.assertGreater(
                    float(abs(weights[index] * mutated - auto)), 1e-6,
                    "dropping the denominator term must break the agreement")

    def test_halving_the_cross_term_is_caught(self) -> None:
        """The factor two comes from the symmetry of the double sum."""
        with mpmath.workdps(80):
            for values, index in _sample_points():
                one = mpmath.mpf(1)
                weights = _weights(values, one)
                support = tuple(values[3:9])
                cross = mpmath.mpf(0)
                for other, point in enumerate(support):
                    if point == 0:
                        continue
                    cross += weights[other] * point * _hp_mp(
                        support[index] * point)
                mutated = (closed_form_factor(values, index, MP.beta)
                           - (1 - MP.beta) * cross)
                auto = autodiff_partial(values, index, MP.beta)
                self.assertGreater(
                    float(abs(weights[index] * mutated - auto)), 1e-6)

    def test_binding_slope_reproduces_the_certified_coefficient(self) -> None:
        report = binding_coefficient_check(MP, ARB, 140)
        self.assertLess(abs(report["slope_rows"][-1]["difference_float"]),
                        1e-30)

    def test_paired_discard_needs_a_quarter(self) -> None:
        """Above y0=1/4 the paired argument pi<=1/2 is simply false."""
        with self.assertRaises(ValueError):
            sweep_paired_link(Fraction(1, 3), 4, 4)
        with mpmath.workdps(40):
            y, z = mpmath.mpf(6) / 10, mpmath.mpf(1)
            self.assertGreater(float(_pi_mp(y, z)), 0.5)
            self.assertLess(
                float(_pi_first_mp(y, z) * _hp_mp(_pi_mp(y, z))), 0.0,
                "the discarded term is negative there, so the threshold is "
                "load-bearing")

    def test_estimate_is_tight_enough_to_be_falsifiable(self) -> None:
        """A tenfold inflated constant must break on a real stratum point."""
        constants = certify_constants(ARB, Fraction(1, 32))
        honest = mpmath.mpf(str(float(constants.m_uniform.lower())))
        with mpmath.workdps(80):
            worst = mpmath.inf
            for values, index in _sample_points():
                zeroed = list(values)
                zeroed[3 + index] = mpmath.mpf(0)
                increment = gap_mp(values, MP.beta) - gap_mp(
                    tuple(zeroed), MP.beta)
                weights = _weights(values, mpmath.mpf(1))
                predicted = honest * weights[index] * values[3 + index]
                if predicted <= 0:
                    continue
                ratio = increment / predicted
                self.assertGreaterEqual(float(ratio), 1.0)
                worst = min(worst, ratio)
            self.assertLess(float(worst), 10.0,
                            "a tenfold inflation of m0 must be refutable")


if __name__ == "__main__":
    unittest.main(verbosity=2)
