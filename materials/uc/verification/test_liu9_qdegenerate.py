#!/usr/bin/env python3
"""Mutation controls for the q-degenerate endpoint certificate.

Each test names and rejects a specific false replacement for the theorem:
a wrong q power, a missing degeneracy factor, q-only symmetry, a q-regime gap,
an inflated tube radius, an inflated endpoint threshold or coefficient, an
omitted active-mean compensation term, or a tampered report.

Run from the ``math`` directory with

    ./.venv/bin/python -I -B uc/verification/test_liu9_qdegenerate.py
"""

from __future__ import annotations

import json
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

import mpmath  # noqa: E402

from liu9_binding import (  # noqa: E402
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_qdegenerate import (  # noqa: E402
    CORE_KAPPA,
    DEFAULT_INTERIOR_Q_MIN,
    DEFAULT_Q_STAR,
    DEFAULT_SMOOTH_CUTOFF,
    DEFAULT_TUBE_RADIUS,
    ENDPOINT_KAPPA,
    ENDPOINT_RATIO,
    _arbf,
    _canonical_digest,
    _core,
    _entropy,
    cubic_weight,
    certify,
    endpoint_delta2,
    endpoint_raw_first_variation,
    formula_d_jet,
    formula_q_jet,
    pure_distance2_mp,
    pure_gap_mp,
    regime_cover_is_complete,
    seam_geometry_holds,
    symmetry_crosscheck,
)

MP = solve_equation_parameters(150)
ARB = certify_equation_parameters(MP)
RESULT = os.path.join(HERE, "results", "liu9-qdegenerate.json")


class QDegenerateMutationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.certificate = certify(ARB)

    def test_wrong_power_of_q_is_rejected(self) -> None:
        """Wrong lemma: replace each cubic displacement weight by its square."""
        q = Fraction(1, 7)
        correct = q * (1 - q) * (1 - 2 * q)
        self.assertEqual(cubic_weight(q, 3), correct)
        self.assertNotEqual(
            cubic_weight(q, 2), correct,
            "the power-two mutation must not pass for the cubic term",
        )

        with mpmath.workdps(100):
            d = mpmath.mpf("1e-16")
            q1, q2 = mpmath.mpf("1e-3"), mpmath.mpf("1e-6")
            ratio1 = pure_gap_mp(q1, d, MP, 100) / pure_distance2_mp(q1, d, MP)
            ratio2 = pure_gap_mp(q2, d, MP, 100) / pure_distance2_mp(q2, d, MP)
            self.assertLess(abs(ratio1 - ratio2), mpmath.mpf("1e-12"))
            wrong1 = pure_gap_mp(q1, d, MP, 100) / (q1**2 * (1 - q1) * d**2)
            wrong2 = pure_gap_mp(q2, d, MP, 100) / (q2**2 * (1 - q2) * d**2)
            self.assertGreater(abs(wrong2 / wrong1), 100,
                               "an extra q power must visibly diverge")

    def test_omitted_degeneracy_factor_is_rejected(self) -> None:
        """Wrong lemma: claim a positive d curvature at the exact q endpoint."""
        with mpmath.workdps(90):
            d = mpmath.mpf(1) / 16
            endpoint_gap = pure_gap_mp(mpmath.mpf(0), d, MP, 90)
            endpoint_dist = pure_distance2_mp(mpmath.mpf(0), d, MP)
            endpoint_hessian = formula_d_jet(mpmath.mpf(0), d, MP).second
            self.assertLess(abs(endpoint_gap), mpmath.mpf("1e-80"))
            self.assertLess(abs(endpoint_dist), mpmath.mpf("1e-80"))
            self.assertLess(abs(endpoint_hessian), mpmath.mpf("1e-80"))
            self.assertGreater(
                float(formula_d_jet(mpmath.mpf(1) / 4, 0, MP).second), 0.0)
            self.assertFalse(endpoint_gap >= mpmath.mpf(CORE_KAPPA.numerator)
                             / CORE_KAPPA.denominator * d**2)

    def test_q_only_symmetry_is_rejected(self) -> None:
        """Wrong lemma: q<->1-q without simultaneously swapping P0 and P1."""
        report = symmetry_crosscheck(MP, 100)
        self.assertLess(report["full_component_swap_residual"], 1e-80)
        self.assertGreater(
            report["q_only_mutation_difference"], 1e-8,
            "the q-only map must change this asymmetric transcribed point",
        )

    def test_threshold_gap_is_rejected(self) -> None:
        """Wrong lemma: endpoint q range stops before the interior range starts."""
        self.assertTrue(regime_cover_is_complete(
            DEFAULT_Q_STAR, DEFAULT_INTERIOR_Q_MIN))
        self.assertFalse(regime_cover_is_complete(
            Fraction(1, 8192), Fraction(1, 4096)),
            "(1/8192,1/4096) leaves a genuine uncovered q interval",
        )

    def test_inflated_tube_radius_is_rejected(self) -> None:
        """Wrong lemma: double rho while retaining eps_sm and q_* unchanged."""
        self.assertTrue(seam_geometry_holds(
            ARB, DEFAULT_SMOOTH_CUTOFF, DEFAULT_Q_STAR,
            DEFAULT_TUBE_RADIUS))
        self.assertFalse(seam_geometry_holds(
            ARB, DEFAULT_SMOOTH_CUTOFF, DEFAULT_Q_STAR,
            Fraction(1, 1024)),
            "rho=1/1024 violates rho^2<=p^2 eps^2 q_*(1-q_*)",
        )

    def test_inflated_q_star_is_rejected(self) -> None:
        """Wrong lemma: use q_*=1/2048 with the certified q-remainder bound."""
        self.assertLessEqual(
            _arbf(DEFAULT_Q_STAR), self.certificate.q_star_ceiling)
        self.assertGreater(
            _arbf(Fraction(1, 2048)), self.certificate.q_star_ceiling,
            "the doubled endpoint threshold exceeds the directed ceiling",
        )

    def test_inflated_endpoint_ratio_is_rejected(self) -> None:
        """Wrong lemma: strengthen H*D>=delta^2/4 to H*D>=delta^2/3."""
        with mpmath.workdps(100):
            d = -mpmath.mpf(1) / 4
            actual = endpoint_raw_first_variation(d, MP) / endpoint_delta2(d, MP)
            self.assertGreater(actual, mpmath.mpf(ENDPOINT_RATIO.numerator)
                               / ENDPOINT_RATIO.denominator)
            self.assertLess(actual, mpmath.mpf(1) / 3,
                            "the one-third mutation is falsified at d=-1/4")
            self.assertGreater(float(self.certificate.endpoint_kappa_lower),
                               float(ENDPOINT_KAPPA))

    def test_omitted_mean_compensation_is_rejected_by_autodiff(self) -> None:
        """Wrong lemma: delete -H*core*(M(Q)-m) from endpoint D."""
        with mpmath.workdps(100):
            d = mpmath.mpf(1) / 8
            direct = formula_q_jet(mpmath.mpf(0), d, MP).first
            correct = endpoint_raw_first_variation(d, MP)
            h_star = MP.p * _entropy(MP.x)
            mutated = correct + h_star * _core(MP) * MP.p * d
            self.assertLess(abs(direct - correct), mpmath.mpf("1e-80"))
            self.assertGreater(abs(direct - mutated), mpmath.mpf("1e-3"),
                               "the active-mean term is load-bearing")

    def test_tampered_or_cornerless_report_is_rejected(self) -> None:
        """Wrong lemma: alter a structural claim or omit boundary corners."""
        with open(RESULT, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
        recorded = payload.pop("report_sha256")
        self.assertEqual(_canonical_digest(payload), recorded)
        self.assertEqual(
            payload["certificate"]["endpoint_sweep"]["boundary_points"], 4)
        self.assertEqual(
            payload["certificate"]["cubic_sweep"]["corner_points"], 4)
        mutated = json.loads(json.dumps(payload))
        mutated["structural_findings"]["q_factors_cancel"] = False
        mutated["certificate"]["endpoint_sweep"]["boundary_points"] = 0
        self.assertNotEqual(_canonical_digest(mutated), recorded,
                            "canonical hashing must reject the wrong report")


if __name__ == "__main__":
    unittest.main(verbosity=2)
