#!/usr/bin/env python3
"""Mutation controls for the Liu nine-variable transverse audit.

Run from the math directory with

    ./.venv/bin/python -I -B uc/verification/test_liu9_ninevar.py
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import unittest
from fractions import Fraction
from pathlib import Path

for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

HERE = Path(__file__).resolve().parent
UC = HERE.parent
if str(UC) not in sys.path:
    sys.path.insert(0, str(UC))

from flint import arb  # noqa: E402
import mpmath  # noqa: E402

from liu9_binding import (  # noqa: E402
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_ninevar import (  # noqa: E402
    CHART_POINTS,
    DECISIVE_EPSILON,
    DECISIVE_Q,
    DECISIVE_Y,
    KAPPA,
    REPORT_PATH,
    EpsilonPolynomial,
    _arbf,
    _canonical_digest,
    active_split_embedding,
    active_split_polynomial_values,
    chart_embedding,
    exact_parameter_intervals,
    exact_sign_certificate,
    formula_has_epsilon_log_argument,
    h_arb,
    mean_preserving_embedding,
    polynomial_expansion,
    support_opening_polynomial_values,
    validate_chart_agreement,
    validate_log_term_claim,
    validate_sign_claim,
    validate_y_grid,
    y_grid,
)
from liu9_objective import _ArbOps, _formula  # noqa: E402


MP = solve_equation_parameters(110)
ARB = certify_equation_parameters(MP)
EXACT_PARAMETERS = exact_parameter_intervals(ARB)
EXACT_SIGN = exact_sign_certificate(EXACT_PARAMETERS)


class NineVariableMutationTest(unittest.TestCase):

    def test_explicit_chart_embedding_has_the_documented_order(self) -> None:
        with mpmath.workdps(100):
            q = mpmath.mpf(2) / 5
            s = MP.x + mpmath.mpf(1) / 1024
            d = -mpmath.mpf(1) / 768
            r = mpmath.mpf(1) / 80
            values = chart_embedding(MP, q, s, d, r)
            self.assertEqual(len(values), 9)
            self.assertEqual(values[2], q)
            self.assertEqual(values[4], 0)
            self.assertEqual(values[5], 0)
            self.assertEqual(values[7], 0)
            self.assertEqual(values[8], 0)
            self.assertLess(
                abs(values[0] - MP.mean / s), mpmath.mpf("1e-90")
            )
            mean = (
                (1 - q) * (values[0] * values[3] + values[1] * values[4]
                           + (1 - values[0] - values[1]) * values[5])
                + q * (values[0] * values[6] + values[1] * values[7]
                       + (1 - values[0] - values[1]) * values[8])
            )
            self.assertLess(abs(mean - MP.mean), mpmath.mpf("1e-90"))

    def test_wrong_transverse_sign_is_rejected(self) -> None:
        validate_sign_claim(EXACT_SIGN.leading, "negative")
        with self.assertRaises(AssertionError):
            validate_sign_claim(EXACT_SIGN.leading, "positive")
        self.assertLess(EXACT_SIGN.leading.hi, Fraction(-1, 2))
        self.assertLess(EXACT_SIGN.full_pencil.hi, Fraction(-1, 2048))

    def test_omitted_epsilon_log_term_is_rejected_when_present(self) -> None:
        active_values = active_split_polynomial_values(
            ARB, _arbf(DECISIVE_Y), _arbf(DECISIVE_Q), arb(1)
        )
        self.assertFalse(formula_has_epsilon_log_argument(
            active_values, ARB.beta, h_arb, arb(1)
        ))
        validate_log_term_claim(
            active_values, ARB.beta, h_arb, arb(1), False
        )

        opening_values = support_opening_polynomial_values(ARB)
        self.assertTrue(formula_has_epsilon_log_argument(
            opening_values, ARB.beta, h_arb, arb(1)
        ))
        with self.assertRaises(AssertionError):
            validate_log_term_claim(
                opening_values, ARB.beta, h_arb, arb(1), False
            )
        validate_log_term_claim(
            opening_values, ARB.beta, h_arb, arb(1), True
        )

    def test_y_grid_cannot_drop_either_endpoint(self) -> None:
        grid = y_grid(intervals=16)
        validate_y_grid(grid)
        self.assertEqual(grid[0], Fraction(1, 32))
        self.assertEqual(grid[-1], Fraction(31, 32))
        with self.assertRaises(AssertionError):
            validate_y_grid(grid[1:])
        with self.assertRaises(AssertionError):
            validate_y_grid(grid[:-1])

    def test_chart_values_are_reproduced_and_bad_embedding_is_rejected(self) -> None:
        report = validate_chart_agreement(ARB, points=CHART_POINTS[:3])
        self.assertTrue(report["passed"])
        self.assertEqual(report["checked_points"], 3)

        def bad_embedding(parameters, q, s, d, r):
            values = list(chart_embedding(parameters, q, s, d, r))
            values[6] += _arbf(Fraction(1, 4096))
            return tuple(values)

        with self.assertRaises(AssertionError):
            validate_chart_agreement(
                ARB, points=CHART_POINTS[:3], embedding=bad_embedding
            )

    def test_automatic_formula_expansion_has_negative_linear_pencil(self) -> None:
        expansion = polynomial_expansion(
            ARB, _arbf(DECISIVE_Y), _arbf(DECISIVE_Q), "active_split"
        )
        pencil = expansion["pencil"]
        self.assertTrue(pencil.c0.contains(0))
        self.assertLess(pencil.c1, 0)
        self.assertEqual(expansion["epsilon_dependent_entropy_arguments"], 0)

        epsilon = _arbf(DECISIVE_EPSILON)
        values = active_split_embedding(
            ARB, _arbf(DECISIVE_Y), epsilon, _arbf(DECISIVE_Q)
        )
        terms = _formula(values, ARB.beta, _ArbOps(False))
        direct_gap = terms.numerator - terms.ehx
        direct_distance = __import__("liu9_ninevar").distance_squared(values, ARB)
        direct_pencil = direct_gap - _arbf(KAPPA) * direct_distance
        residual = direct_pencil - pencil.evaluate(epsilon)
        self.assertTrue(residual.contains(0))
        self.assertLess(direct_pencil, 0)

    def test_mean_preserving_control_is_distinct_and_positive(self) -> None:
        epsilon = _arbf(DECISIVE_EPSILON)
        y = _arbf(DECISIVE_Y)
        q = _arbf(DECISIVE_Q)
        values = mean_preserving_embedding(ARB, y, epsilon, q)
        mean = values[0] * ARB.x + values[1] * y
        self.assertTrue((mean - ARB.mean).contains(0))
        expansion = polynomial_expansion(
            ARB, y, q, "mean_preserving"
        )
        self.assertGreater(expansion["pencil"].c1, 0)
        self.assertGreater(EXACT_SIGN.mean_preserving_leading.lo,
                           Fraction(1, 64))

    def test_report_digest_and_units_are_load_bearing(self) -> None:
        payload = json.loads(REPORT_PATH.read_text(encoding="utf-8"))
        digest = payload.pop("report_sha256")
        self.assertEqual(digest, _canonical_digest(payload))
        self.assertEqual(payload["verdict"], "REFUTED")
        self.assertEqual(payload["kappa"]["units"],
                         "raw-gap per distance-squared")
        self.assertFalse(payload["counterexample"]["liu_mean_feasible"])
        self.assertEqual(
            payload["counterexample"]["units_identity"],
            "gap=EHX*(objective-1)",
        )
        canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        self.assertEqual(hashlib.sha256(canonical.encode()).hexdigest(), digest)


if __name__ == "__main__":
    unittest.main(verbosity=2)
