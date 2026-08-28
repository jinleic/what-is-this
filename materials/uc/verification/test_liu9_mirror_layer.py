#!/usr/bin/env python3
"""Mutation tests for the Liu H2 mirror boundary-layer estimate.

A control that cannot fail proves nothing.  These tests check that
``uc/liu9_mirror_layer.py`` rejects each way the estimate could be wrong:
a collapsed beta split, a decorative tube hypothesis, an invalid paired-link
threshold, an unsound dual bound, and a quotient/raw-gap units mix-up.

Run from the repository root with

    math/.venv/bin/python -I -B math/uc/verification/test_liu9_mirror_layer.py
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

from liu9_binding import (  # noqa: E402
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_boundary_layer import (  # noqa: E402
    _support,
    _weights,
    closed_form_factor,
    entropy_at_star,
    mean_of,
)
from liu9_mirror_layer import (  # noqa: E402
    MU_GRID,
    SUPPORT_VARIABLE_INDEX,
    asymptotic_coefficient,
    best_constants,
    certify_constants,
    certify_theta,
    radius_ceiling,
    smooth_kappa_ceiling,
    sweep_cross_link,
    sweep_paired_link,
)

MP = solve_equation_parameters(120)
ARB = certify_equation_parameters(MP)
T0 = Fraction(1, 64)
RHO = Fraction(1, 10)


def _arb_midpoint(value: object) -> float:
    return (float(value.lower()) + float(value.upper())) / 2


def _sampled_dual_maximum(mu: int, t0: Fraction, points: int) -> mpmath.mpf:
    """Sample phi(b)-mu*g(b) on an endpoint-inclusive uniform mesh."""
    with mpmath.workdps(80):
        lo = mpmath.mpf(1) / 2
        hi = 1 - mpmath.mpf(t0.numerator) / t0.denominator
        worst = -mpmath.inf
        for index in range(points):
            b = lo + (hi - lo) * index / (points - 1)
            phi = b * mpmath.log(b / (1 - b))
            g = (b * (b - MP.x)) ** 2
            worst = max(worst, phi - mu * g)
        return +worst


def _discriminating_slope() -> mpmath.mpf:
    """Measure the W1=1/2, A1=1 slope without using the report's prediction."""
    with mpmath.workdps(220):
        one = mpmath.mpf(1)
        half = one / 2
        zero = mpmath.mpf(0)
        values_by_exponent = []
        for exponent in (40, 80):
            moving = one - mpmath.mpf(10) ** (-exponent)
            values = (half, half, half, one, moving, zero,
                      MP.x, MP.x, zero)
            values_by_exponent.append(
                (exponent, closed_form_factor(values, 1, MP.beta)))
        (near_exponent, near), (far_exponent, far) = values_by_exponent
        return ((far - near)
                / ((far_exponent - near_exponent) * mpmath.log(10)))


class MirrorLayerControlTest(unittest.TestCase):

    def test_leading_coefficient_matches_the_measured_slope(self) -> None:
        """All six measured slopes must retain the separate W1 and A1 terms."""
        with mpmath.workdps(120):
            report = asymptotic_coefficient(MP, 120)
            rows = report["rows"]
            self.assertEqual(len(rows), 6)
            expected_masses = (
                (1 - MP.p, 1 - MP.p),
                (mpmath.mpf(2) / 5, mpmath.mpf(2) / 5),
                (mpmath.mpf(1) / 2, mpmath.mpf(1) / 2),
                (mpmath.mpf(7) / 10, mpmath.mpf(7) / 10),
                (mpmath.mpf(1) / 2, mpmath.mpf(1)),
                (mpmath.mpf(1), mpmath.mpf(1)),
            )
            for row, (w1, a1) in zip(rows, expected_masses):
                measured = mpmath.mpf(row["measured_slope"])
                predicted = 1 - 2 * (1 - MP.beta) * w1 - 2 * MP.beta * a1
                self.assertLess(abs(measured - predicted),
                                mpmath.mpf("1e-18"))
                self.assertLess(abs(mpmath.mpf(row["difference"])),
                                mpmath.mpf("1e-30"))

            discriminating = rows[4]
            self.assertEqual(discriminating["weight_at_support_one"], "0.5")
            self.assertEqual(discriminating["component_mass_at_support_one"],
                             "1.0")
            self.assertLess(abs(_discriminating_slope() + MP.beta),
                            mpmath.mpf("1e-30"))

    def test_dropping_the_beta_split_is_caught(self) -> None:
        """Substituting global W1 for within-component A1 predicts zero."""
        with mpmath.workdps(120):
            w1 = mpmath.mpf(1) / 2
            a1 = mpmath.mpf(1)
            correct = 1 - 2 * (1 - MP.beta) * w1 - 2 * MP.beta * a1
            mutated = 1 - 2 * (1 - MP.beta) * w1 - 2 * MP.beta * w1
            self.assertLess(abs(correct + MP.beta), mpmath.mpf("1e-100"))
            self.assertLess(abs(mutated), mpmath.mpf("1e-100"))
            self.assertGreater(abs(_discriminating_slope() - mutated),
                               mpmath.mpf("0.09"))

    def test_refuted_threshold_is_load_bearing(self) -> None:
        """At rho=1/10, t0=1/32 fails while the thinner 1/64 layer works."""
        refuted = best_constants(
            ARB, Fraction(1, 32), RHO, pieces=512, mus=MU_GRID)
        certified = best_constants(
            ARB, Fraction(1, 64), RHO, pieces=512, mus=MU_GRID)
        self.assertFalse(refuted.certified())
        self.assertTrue(certified.certified())

    def test_theta_is_a_genuine_upper_bound(self) -> None:
        """The interval lambda bounds samples and cannot be reduced by 10%."""
        for mu in (12, 42, 140, 520):
            dual = certify_theta(ARB.x, T0, RHO, mu, pieces=512)
            certified_upper = float(dual["lambda"].upper())
            sampled_maximum = float(_sampled_dual_maximum(mu, T0, 4097))
            self.assertLessEqual(
                sampled_maximum, certified_upper + 1e-12,
                f"mu={mu} exceeds its certified lambda")
            self.assertGreater(certified_upper, 0.0)
            self.assertGreater(
                sampled_maximum, 0.9 * certified_upper,
                f"mu={mu} has an upper bound too loose to falsify")

    def test_dual_multiplier_choice_is_sound(self) -> None:
        """Each multiplier is valid and best_constants really takes the max."""
        pieces = 256
        candidates = [
            certify_constants(ARB, T0, RHO, mu, pieces=pieces)
            for mu in MU_GRID
        ]
        for candidate in candidates:
            sampled_maximum = float(
                _sampled_dual_maximum(candidate.mu, T0, 513))
            self.assertLessEqual(
                sampled_maximum,
                float(candidate.theta_lambda.upper()) + 1e-12,
                f"mu={candidate.mu} is not individually a valid dual bound")

        best = best_constants(
            ARB, T0, RHO, pieces=pieces, mus=MU_GRID)
        best_midpoint = _arb_midpoint(best.kappa)
        candidate_midpoints = [_arb_midpoint(item.kappa)
                               for item in candidates]
        self.assertIn(best.mu, MU_GRID)
        self.assertLessEqual(best_midpoint, max(candidate_midpoints) + 1e-12)
        certified_candidates = [item for item in candidates
                                if item.certified()]
        self.assertTrue(certified_candidates)
        for candidate in certified_candidates:
            self.assertGreaterEqual(
                best_midpoint + 1e-12, _arb_midpoint(candidate.kappa),
                f"best_constants lost certified mu={candidate.mu}")

    def test_tube_bound_on_W1_is_needed(self) -> None:
        """Moving beyond the reported critical radius makes M negative."""
        row = next(item for item in radius_ceiling(ARB)["rows"]
                   if item["t0"] == str(T0))
        inflated_rho = Fraction(1, 4)
        self.assertLess(row["critical_rho_lower"], float(inflated_rho))
        constants = certify_constants(
            ARB, T0, inflated_rho, mu=42, pieces=64)
        self.assertLessEqual(float(constants.leading.upper()), 0.0)

    def test_reduction_raises_the_mean(self) -> None:
        """Raising every constructed mirror atom to one preserves feasibility."""
        mp = mpmath.mpf
        with mpmath.workdps(80):
            values = (mp("0.2"), mp("0.3"), mp("0.4"),
                      mp("0.99"), mp("0.4"), mp("0.985"),
                      mp("0.3"), mp("0.999"), mp("0.7"))
            support = _support(values)
            weights = _weights(values, mp(1))
            raised = list(values)
            expected_increase = mp(0)
            cutoff = 1 - mp(T0.numerator) / T0.denominator
            for index, point in enumerate(support):
                if point > cutoff:
                    raised[SUPPORT_VARIABLE_INDEX[index]] = mp(1)
                    expected_increase += weights[index] * (1 - point)
            before = mean_of(values, mp(1))
            after = mean_of(tuple(raised), mp(1))
            self.assertGreater(after, before)
            self.assertLess(abs((after - before) - expected_increase),
                            mp("1e-70"))

    def test_paired_and_cross_links_reject_a_wrong_threshold(self) -> None:
        """The paired rewriting is invalid once y=1-t may fall below 1/2."""
        with self.assertRaises(ValueError):
            sweep_paired_link(Fraction(3, 5), 4, 4)
        paired = sweep_paired_link(Fraction(1, 2), 4, 4)
        cross = sweep_cross_link(T0, 3, 4)
        self.assertEqual(paired["boxes"], 16)
        self.assertEqual(cross["boxes"], 15)

    def test_kappa_is_reported_in_raw_gap_units(self) -> None:
        """The raw-gap ceiling is H* times liu9_tube's quotient ceiling."""
        quotient = smooth_kappa_ceiling(ARB, units="phi")
        raw_gap = smooth_kappa_ceiling(ARB, units="gap")
        h_star = entropy_at_star(ARB)
        expected = quotient * h_star
        quotient_midpoint = _arb_midpoint(quotient)
        raw_midpoint = _arb_midpoint(raw_gap)
        expected_midpoint = _arb_midpoint(expected)
        self.assertAlmostEqual(raw_midpoint, expected_midpoint, delta=1e-15)
        self.assertAlmostEqual(quotient_midpoint, 0.70046750915,
                               delta=5e-11)
        self.assertAlmostEqual(raw_midpoint, 0.3871250877692325,
                               delta=5e-13)
        self.assertGreater(abs(quotient_midpoint - raw_midpoint), 0.3)


    def test_w1_bound_uses_the_stratum_infimum_not_the_endpoint(self):
        """g(1-t0), not g(1).  The substitution overclaims the radius.

        `g(b)=[b(b-x)]^2` is smallest at the inner edge of the mirror
        stratum, so that is the divisor the W1 bound must use.  Substituting
        the endpoint `g(1)=[1*(1-x)]^2=0.0956` -- which is the number the
        module docstring quotes as motivation, one line above the real bound
        -- inflates the critical radius from 0.1795 to 0.2061 at t0=1/32.
        That is a 15% overclaim in the unsafe direction, so it must be
        rejected explicitly rather than left to a reader's care.
        """
        from fractions import Fraction as F
        import liu9_mirror_layer as M
        from flint import arb
        x = ARB.x
        for t0 in (F(1, 32), F(1, 64)):
            report = M.certify_face_infimum(x, t0)
            infimum = arb(report["infimum_g_lower"])
            endpoint = arb(report["endpoint_g_lower"])
            self.assertLess(float(infimum), float(endpoint),
                            "g must be strictly smaller at the inner edge")
            factor = 2 * (1 - ARB.beta)
            honest = (infimum * (1 - 2 * ARB.beta) / factor).sqrt()
            inflated = (endpoint * (1 - 2 * ARB.beta) / factor).sqrt()
            self.assertGreater(float(inflated.lower()),
                               float(honest.upper()),
                               "the endpoint substitution must be detectably "
                               "more permissive")
            certified = M.best_constants(ARB, t0, F(1, 10), pieces=512)
            self.assertLessEqual(
                float(certified.w1_ceiling.lower()),
                float((arb(1) / 100 / infimum).upper()) + 1e-12)
            self.assertGreater(
                float(certified.w1_ceiling.upper()),
                float((arb(1) / 100 / endpoint).lower()),
                "the module must be using the larger, honest W1 ceiling")

    def test_face_infimum_rejects_a_stratum_that_reaches_below_x(self):
        """If 1-t0 <= x the infimum is not at the inner edge at all."""
        from fractions import Fraction as F
        import liu9_mirror_layer as M
        with self.assertRaises(AssertionError):
            M.certify_face_infimum(ARB.x, F(2, 5))


    def test_w1_division_stays_conservative_when_the_face_is_wide(self):
        """The safe end must govern even when g's enclosure is not tiny.

        At 320 bits the face cost `g(1-t0)` encloses to about 1e-70, so lower,
        midpoint and upper agree in every printed digit and any substitution
        among them would be invisible.  This test widens `x` deliberately so
        the three ends separate, then checks that the ceiling the module would
        consume is the one built from the SMALLEST face cost -- the largest,
        most permissive-to-refute W1 -- and that a midpoint substitution is
        detectably weaker.
        """
        from fractions import Fraction as F
        import liu9_mirror_layer as M
        from flint import arb
        wide = ARB.x + arb(0).union(arb(1) / 1000)
        t0 = F(1, 64)
        rho2 = M._arbf(F(1, 10)) ** 2
        face = M._g(1 - M._arbf(t0), wide)
        self.assertGreater(float(face.upper()) - float(face.lower()), 1e-6,
                           "the widened face must actually be wide")
        conservative = rho2 / arb(face.lower())
        midpoint = rho2 / arb(face.mid())
        interval = rho2 / face
        self.assertGreater(float(conservative.lower()),
                           float(midpoint.upper()),
                           "the midpoint substitution must be weaker")
        self.assertGreaterEqual(float(interval.upper()),
                                float(conservative.lower()),
                                "interval division must expose the safe end")
        leading_safe = 1 - 2 * ARB.beta - 2 * (1 - ARB.beta) * interval
        leading_mid = 1 - 2 * ARB.beta - 2 * (1 - ARB.beta) * midpoint
        self.assertLess(float(leading_safe.lower()),
                        float(leading_mid.lower()),
                        "the certified leading coefficient must be the "
                        "smaller, honest one")


if __name__ == "__main__":
    unittest.main(verbosity=2)
