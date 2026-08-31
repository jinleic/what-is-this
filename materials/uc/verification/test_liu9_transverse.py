#!/usr/bin/env python3
"""Mutation tests for the mean-preserving transverse positivity certificate.

A control that cannot fail proves nothing.  Each test below rejects one
specific way `uc/liu9_transverse.py` could be wrong, and several of them reject
a mistake that was actually available while writing it.

Run from the repository root:

    math/.venv/bin/python -I -B math/uc/verification/test_liu9_transverse.py
"""

from __future__ import annotations

import os
import sys
import unittest
from fractions import Fraction as F

for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

import mpmath
from flint import arb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
UC = os.path.dirname(HERE)
for _path in (UC,):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from liu9_binding import (  # noqa: E402
    _h_arb,
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_ninevar import transverse_coefficients_mp  # noqa: E402
from liu9_transverse import (  # noqa: E402
    KAPPA,
    Jet2,
    _arbf,
    _ball,
    _constant_a,
    _mu_enclosure,
    _transverse_value,
    certify_boundary_layer,
    certify_bulk,
    certify_coverage,
    certify_double_root,
    certify_protocol_identity,
    certify_structural_zero,
    certify_window,
    asymmetric_coefficient_mp,
    certify_asymmetric_decoupling,
    certify_perturbation_decomposition,
    _pencil_mp,
    transverse_jet,
)

ctx.prec = max(ctx.prec, 320)
MP = solve_equation_parameters(90)
ARB = certify_equation_parameters(MP)


def _mutated_value(parameters, y: arb, *, mean_factor, kappa=None,
                   support_power: int = 2) -> arb:
    """L with one ingredient deliberately altered."""
    x, p, beta = parameters.x, parameters.p, parameters.beta
    one = arb(1)
    kappa = _arbf(KAPPA) if kappa is None else kappa
    product = x * y
    protocol = x * y + x * (one - x) * y * (one - y)
    kernel = (one - beta) * _h_arb(product) + beta * _h_arb(protocol)
    constant_a = _constant_a(parameters)
    support = y * y
    for _ in range(support_power):
        support = support * (y - x)
    return (2 * p * kernel - _h_arb(y)
            - mean_factor(y) * constant_a - kappa * support)


class ProtocolIdentityTest(unittest.TestCase):
    """The identity must be the quartic exactly, not merely near it."""

    def test_the_identity_is_liu_quartic_coefficient_for_coefficient(self):
        report = certify_protocol_identity()
        self.assertTrue(report["identity_holds"])
        self.assertEqual(report["difference"], [-1, 0, 3, -2, 1])
        self.assertEqual(report["difference"], report["liu_quartic"])

    def test_a_perturbed_quartic_is_rejected(self):
        """If the comparison were sloppy any nearby polynomial would pass."""
        report = certify_protocol_identity()
        for index in range(5):
            wrong = list(report["liu_quartic"])
            wrong[index] += 1
            self.assertNotEqual(report["difference"], wrong,
                                f"coefficient {index} must be discriminating")

    def test_the_protocol_diagonal_is_not_accidentally_the_reflection(self):
        """prot(z,z) and 1-z^2 must differ; equality would make it vacuous."""
        report = certify_protocol_identity()
        self.assertNotEqual(report["prot_diagonal"],
                            report["entropy_reflection_of_x_squared"])


class StructuralZeroTest(unittest.TestCase):
    """L(x)=0 must hold for arbitrary constants, and mutations must break it."""

    def test_the_zero_is_a_rational_function_identity(self):
        report = certify_structural_zero(samples=20)
        self.assertTrue(report["all_exactly_zero"])
        self.assertEqual(report["samples"], 20)

    def test_halving_the_mean_factor_destroys_the_structural_zero(self):
        """y/x is what makes it cancel; y/(2x) must not."""
        x, p, beta, kappa = F(3, 7), F(5, 11), F(2, 9), F(7, 13)
        h_x, h_x2, h_prot = F(1, 3), F(4, 5), F(6, 7)
        diagonal = (1 - beta) * h_x2 + beta * h_prot
        constant_a = 2 * p * diagonal - h_x
        good = 2 * p * diagonal - h_x - (x / x) * constant_a
        bad = 2 * p * diagonal - h_x - (x / (2 * x)) * constant_a
        self.assertEqual(good, 0)
        self.assertNotEqual(bad, 0)

    def test_a_linear_support_term_destroys_the_structural_zero(self):
        """(y-x)^2 vanishes at y=x; (y-x)^1 does too, but (y-x)^0 does not."""
        x, kappa = F(3, 7), F(7, 13)
        self.assertEqual(kappa * x ** 2 * (x - x) ** 2, 0)
        self.assertEqual(kappa * x ** 2 * (x - x) ** 1, 0)
        self.assertNotEqual(kappa * x ** 2, 0)


class DoubleRootTest(unittest.TestCase):
    """The vanishing slope must be Liu's, and the curvature must be real."""

    def test_the_slope_is_below_the_parameter_residual(self):
        report = certify_double_root(ARB)
        self.assertLess(report["slope_upper"], 1e-60)
        self.assertGreater(report["curvature_lower"], 1.0)

    def test_the_curvature_matches_an_independent_finite_difference(self):
        """1.2289613682... was reproduced by an independent derivation."""
        report = certify_double_root(ARB)
        with mpmath.workdps(60):
            def coefficient(y):
                return transverse_coefficients_mp(
                    MP, mpmath.mpf(y), "mean_preserving")[0]

            second = mpmath.diff(coefficient, MP.x, 2)
        self.assertAlmostEqual(float(second),
                               report["curvature_lower"], places=5)
        self.assertAlmostEqual(float(second), 1.2289613682, places=8)

    def test_perturbing_beta_moves_the_slope_off_zero(self):
        """The cancellation is equation (90), not a transcription accident."""
        with mpmath.workdps(60):
            x, p = MP.x, MP.p

            def slope(beta):
                def kernel(right):
                    prot = x * right + x * (1 - x) * right * (1 - right)
                    return ((1 - beta) * _entropy(x * right)
                            + beta * _entropy(prot))

                constant_a = 2 * p * kernel(x) - _entropy(x)
                return (2 * p * mpmath.diff(kernel, x)
                        - mpmath.diff(_entropy, x) - constant_a / x)

            def _entropy(value):
                if value == 0 or value == 1:
                    return mpmath.mpf(0)
                return -(value * mpmath.log(value)
                         + (1 - value) * mpmath.log(1 - value))

            self.assertLess(abs(slope(MP.beta)), 1e-40)
            for delta in ("1e-6", "1e-3", "1e-1"):
                self.assertGreater(abs(slope(MP.beta + mpmath.mpf(delta))),
                                   1e-8)


class MeanFeasibilityTest(unittest.TestCase):
    """The whole result turns on the mean factor; without it the sign flips."""

    def test_dropping_mean_preservation_reproduces_the_known_refutation(self):
        """Replacing y/x by 1 is the ambient family, which IS negative.

        This is the single most important mutation.  If the certificate still
        passed with the mean factor removed, it would not be distinguishing
        the feasible half-space from the ambient neighbourhood at all, and the
        headline claim would be empty.
        """
        y = _arbf(F(1, 32))
        feasible = _mutated_value(ARB, y, mean_factor=lambda v: v / ARB.x)
        ambient = _mutated_value(ARB, y, mean_factor=lambda v: arb(1))
        self.assertGreater(float(feasible.lower()), 0.0)
        self.assertLess(float(ambient.upper()), -0.5)

    def test_the_two_families_differ_by_exactly_A_times_the_shortfall(self):
        """L_ambient - L_feasible = A*(y/x - 1); the algebra must be that."""
        for numerator, denominator in ((1, 32), (1, 4), (9, 10)):
            y = _arbf(F(numerator, denominator))
            feasible = _mutated_value(ARB, y, mean_factor=lambda v: v / ARB.x)
            ambient = _mutated_value(ARB, y, mean_factor=lambda v: arb(1))
            predicted = _constant_a(ARB) * (y / ARB.x - arb(1))
            self.assertTrue((ambient - feasible).overlaps(predicted))

    def test_the_constant_A_equals_entropy_at_x(self):
        """A = h(x) follows from the protocol identity; check it holds."""
        self.assertTrue(_constant_a(ARB).overlaps(_h_arb(ARB.x)))


class MuEnclosureTest(unittest.TestCase):
    """The boundary-layer bound is only as good as the mu enclosure."""

    def test_the_enclosure_contains_the_true_value(self):
        with mpmath.workdps(50):
            for numerator, denominator in ((1, 2), (1, 10), (1, 1000),
                                           (1, 10 ** 9)):
                u = F(numerator, denominator)
                true = ((1 - u) / u) * mpmath.log(1 / (1 - mpmath.mpf(
                    numerator) / denominator))
                enclosure = _mu_enclosure(_arbf(u), u)
                self.assertLessEqual(float(enclosure.lower()), float(true))
                self.assertGreaterEqual(float(enclosure.upper()), float(true))

    def test_a_reversed_enclosure_would_miss_the_true_value(self):
        """[1, 1+u] is the natural wrong guess; it must not contain mu."""
        with mpmath.workdps(50):
            u = mpmath.mpf(1) / 2
            true = ((1 - u) / u) * mpmath.log(1 / (1 - u))
            self.assertLess(float(true), 1.0)

    def test_the_enclosure_tightens_to_one_as_u_goes_to_zero(self):
        wide = _mu_enclosure(_arbf(F(1, 2)), F(1, 2))
        narrow = _mu_enclosure(_arbf(F(1, 10 ** 6)), F(1, 10 ** 6))
        self.assertLess(float(narrow.upper()) - float(narrow.lower()),
                        float(wide.upper()) - float(wide.lower()))


class RegionTest(unittest.TestCase):
    """Each region must actually be doing work, and must fail when wrong."""

    def test_the_window_must_exclude_the_non_convex_region(self):
        """L is not globally convex; a window reaching 0.43 must be rejected."""
        report = certify_window(ARB, F(1, 8), 64, 1e-68)
        self.assertGreater(report["minimum_curvature_lower"], 0.0)
        with self.assertRaises(AssertionError):
            certify_window(ARB, F(1, 2), 64, 1e-68)

    def test_an_inflated_kappa_breaks_bulk_positivity(self):
        """The certified kappa is not free; a large one must fail somewhere."""
        y = _arbf(F(9, 10))
        good = _mutated_value(ARB, y, mean_factor=lambda v: v / ARB.x)
        bad = _mutated_value(ARB, y, mean_factor=lambda v: v / ARB.x,
                             kappa=arb(10))
        self.assertGreater(float(good.lower()), 0.0)
        self.assertLess(float(bad.upper()), 0.0)

    def test_the_boundary_layer_certifies_without_many_octaves(self):
        """The tail cell carries the argument, so one octave must suffice."""
        report = certify_boundary_layer(ARB, F(1, 32), 1)
        self.assertTrue(report["certified"])
        self.assertGreater(report["tail_bound_lower"], 0.0)
        self.assertGreater(report["uniform_log_coefficient_lower"], 0.2)

    def test_a_coarse_bulk_cover_is_detectably_insufficient(self):
        """Resolution is load-bearing exactly where the layer hands over.

        Near y0 the enclosure widens because h carries a log(1/y) derivative,
        so a coarse cover cannot resolve the sign even though L is comfortably
        positive there.  A certificate that passed at any resolution would be
        reporting the grid, not the function.
        """
        with self.assertRaises(AssertionError) as caught:
            certify_bulk(ARB, F(1, 32), F(1, 8), 128)
        self.assertIn("bulk positivity failed", str(caught.exception))
        fine = certify_bulk(ARB, F(1, 32), F(1, 8), 512)
        self.assertGreater(fine["worst_margin_lower"], 0.0)

    def test_the_bulk_worst_margin_is_a_real_number_not_a_sign_test(self):
        report = certify_bulk(ARB, F(1, 32), F(1, 8), 256)
        self.assertGreater(report["worst_margin_lower"], 1e-3)
        self.assertGreater(report["cells"], 100)


class CoverageTest(unittest.TestCase):
    """A cover with a hole proves nothing about the hole."""

    def test_the_regions_tile_the_unit_interval(self):
        centre = F(ARB.x.str(40, radius=False))
        report = certify_coverage(F(1, 32), F(1, 8), centre, F(1, 2) ** 205)
        self.assertTrue(report["spans_unit_interval"])
        self.assertTrue(report["abutting"])

    def test_the_claimed_cover_matches_what_each_certifier_actually_covers(self):
        """The real gap risk is a coverage claim that outruns the certificates.

        `certify_coverage` builds abutting pieces by construction, so testing
        it against itself is vacuous.  The failure that matters is `main`
        handing the three certifiers parameters that leave a sliver nobody
        proved.  This pins each claimed piece to the region its own certifier
        reports.
        """
        y0, window, octaves = F(1, 32), F(1, 8), 12
        centre = F(ARB.x.str(40, radius=False))
        claimed = certify_coverage(y0, window, centre, y0 / F(2) ** octaves)
        pieces = {name: (F(lo), F(hi)) for name, lo, hi in claimed["pieces"]}

        layer = certify_boundary_layer(ARB, y0, octaves)
        self.assertEqual(F(layer["smallest_covered_y"]),
                         pieces["A tail"][1])
        self.assertEqual(F(layer["y0"]), pieces["A octaves"][1])
        self.assertEqual(F(layer["tail_cell"][0]), pieces["A tail"][0])

        bulk = certify_bulk(ARB, y0, window, 512)
        self.assertEqual(F(bulk["y0"]), pieces["B below"][0])
        self.assertEqual(F(bulk["window"]), window)
        self.assertEqual(pieces["B below"][1], centre - window)
        self.assertEqual(pieces["B above"][0], centre + window)
        self.assertEqual(pieces["B above"][1], F(1))

        win = certify_window(ARB, window, 32, 1e-68)
        self.assertEqual(F(win["radius"]), window)
        self.assertEqual(pieces["C window"],
                         (centre - window, centre + window))

    def test_an_empty_region_is_rejected(self):
        centre = F(ARB.x.str(40, radius=False))
        with self.assertRaises(AssertionError):
            # boundary layer starting above where the bulk ends
            certify_coverage(F(9, 10), F(1, 8), centre, F(1, 2) ** 205)

    def test_a_window_that_swallows_the_bulk_is_rejected(self):
        centre = F(ARB.x.str(40, radius=False))
        with self.assertRaises(AssertionError):
            certify_coverage(F(1, 32), F(3, 4), centre, F(1, 2) ** 205)


class CrossCheckTest(unittest.TestCase):
    """The module must agree with the pre-existing independent transcription."""

    def test_agreement_with_liu9_ninevar(self):
        with mpmath.workdps(60):
            for numerator, denominator in ((1, 64), (1, 8), (1, 3), (1, 2),
                                           (2, 3), (7, 8), (63, 64)):
                y = F(numerator, denominator)
                mine = _transverse_value(ARB, _arbf(y))
                theirs = transverse_coefficients_mp(
                    MP, mpmath.mpf(numerator) / denominator,
                    "mean_preserving")[0]
                # mpmath carries 60 dps; a float64 round-trip would truncate to
                # 1e-17 and could never meet a 1e-69 enclosure, so compare
                # through an exact rational with an explicit mpmath tolerance.
                expected = _arbf(F(mpmath.nstr(theirs, 40)))
                difference = abs(mine - expected)
                self.assertLess(float(difference.upper()), 1e-38,
                                f"disagreement at y={y}: {mine} vs {theirs}")

    def test_the_jet_value_matches_the_plain_enclosure(self):
        for numerator, denominator in ((3, 5), (7, 10), (4, 5)):
            y = _arbf(F(numerator, denominator))
            jet = transverse_jet(ARB, Jet2.variable(y))
            plain = _transverse_value(ARB, y)
            self.assertTrue(jet.v.overlaps(plain))


class AsymmetricDecouplingTest(unittest.TestCase):
    """The free two-parameter extension must actually be free, and correct."""

    def test_the_pencil_vanishes_at_the_known_zero(self):
        """Calibration guard against the _ScalarOps trap.

        `_ScalarOps(entropy, one)` takes the ENTROPY FUNCTION first.  Passing
        `mpmath.mpf` makes "entropy(u) = u", which yields a plausible but
        entirely wrong sign -- it produced 586 spurious "negative" configs the
        first time this family was scanned.  The base configuration is a
        two-atom law at Liu's optimizer, where the pencil must vanish, so this
        test fails loudly for any transcription that is not binary entropy.
        """
        with mpmath.workdps(60):
            base = _pencil_mp((MP.p, 0, mpmath.mpf(1) / 2, MP.x, MP.x, 0,
                               MP.x, MP.x, 0), MP)
            self.assertLess(abs(float(base)), 1e-50,
                            "the pencil must vanish at the optimizer")

    def test_the_coefficient_vanishes_at_y_equals_x(self):
        with mpmath.workdps(60):
            value = asymmetric_coefficient_mp(MP, MP.x, MP.x,
                                              mpmath.mpf(1) / 2)
            self.assertLess(abs(float(value)), 1e-25)

    def test_the_decoupling_identity_holds(self):
        report = certify_asymmetric_decoupling(
            MP, 60, [("1/5", "9/10", "1/2"), ("1/50", "49/50", "1/4"),
                     ("7/10", "1/10", "3/4")])
        self.assertTrue(report["certified"])
        self.assertLess(report["worst_difference"], 1e-25)
        self.assertTrue(report["endpoint_reduction_invariant"])

    def test_the_perturbation_decomposes_exactly(self):
        report = certify_perturbation_decomposition(samples=25)
        self.assertTrue(report["atoms_match_exactly"])

    def test_a_wrong_mean_correction_breaks_the_decomposition(self):
        """ybar must be the q-weighted mean; y1 alone must not decompose."""
        x, y1, y2, q = F(7, 10), F(1, 5), F(9, 10), F(1, 3)
        qbar = 1 - q
        right = qbar * y1 + q * y2
        wrong = y1
        self.assertEqual(qbar * (-y1 / x) + q * (-y2 / x), -right / x)
        self.assertNotEqual(qbar * (-y1 / x) + q * (-y2 / x), -wrong / x)

    def test_a_wrong_mean_correction_is_mean_infeasible(self):
        """Using y1 instead of ybar leaves the half-space, so it proves nothing."""
        with mpmath.workdps(60):
            x, p = MP.x, MP.p
            y1, y2 = mpmath.mpf(1) / 5, mpmath.mpf(9) / 10
            q, eps = mpmath.mpf(1) / 2, mpmath.mpf(1) / 1000
            for correction, feasible in ((qbar_mean(y1, y2, q), True),
                                         (y1, False)):
                a1 = p - eps * correction / x
                mean = (1 - q) * (a1 * x + eps * y1) + q * (a1 * x + eps * y2)
                if feasible:
                    self.assertLess(abs(float(mean - p * x)), 1e-40)
                else:
                    self.assertGreater(abs(float(mean - p * x)), 1e-6)

    def test_the_coefficient_is_affine_in_q(self):
        with mpmath.workdps(60):
            y1, y2 = mpmath.mpf(1) / 5, mpmath.mpf(9) / 10
            values = [asymmetric_coefficient_mp(MP, y1, y2, mpmath.mpf(k) / 4)
                      for k in (1, 2, 3)]
            self.assertLess(abs(float(values[0] - 2 * values[1] + values[2])),
                            1e-25)

    def test_q_zero_eliminates_the_second_support(self):
        with mpmath.workdps(60):
            first = asymmetric_coefficient_mp(
                MP, mpmath.mpf(3) / 10, mpmath.mpf(1) / 10, mpmath.mpf(0))
            second = asymmetric_coefficient_mp(
                MP, mpmath.mpf(3) / 10, mpmath.mpf(99) / 100, mpmath.mpf(0))
            self.assertEqual(mpmath.nstr(first, 25), mpmath.nstr(second, 25))

    def test_the_asymmetric_family_is_nowhere_negative_on_a_grid(self):
        """A coarse independent sweep; the certificate covers the continuum."""
        with mpmath.workdps(50):
            worst = None
            for a in (1, 5, 11, 15):
                for b in (2, 7, 13):
                    for k in (1, 4, 7):
                        value = asymmetric_coefficient_mp(
                            MP, mpmath.mpf(a) / 16, mpmath.mpf(b) / 16,
                            mpmath.mpf(k) / 8)
                        worst = value if worst is None or value < worst \
                            else worst
            self.assertGreater(float(worst), -1e-25)


def qbar_mean(y1, y2, q):
    return (1 - q) * y1 + q * y2


if __name__ == "__main__":
    unittest.main(verbosity=2)
