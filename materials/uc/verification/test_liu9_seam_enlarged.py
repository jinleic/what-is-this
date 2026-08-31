#!/usr/bin/env python3
"""Mutation tests for the enlarged q-degenerate seam certificate."""

from __future__ import annotations

import os
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

HERE = Path(__file__).resolve().parent
UC = HERE.parent
if str(UC) not in sys.path:
    sys.path.insert(0, str(UC))

from liu9_binding import certify_equation_parameters, solve_equation_parameters  # noqa: E402
from liu9_qdegenerate import (  # noqa: E402
    certify,
    regime_cover_is_complete,
    seam_geometry_holds,
)

ARB = certify_equation_parameters(solve_equation_parameters(100))


class EnlargedSeamTest(unittest.TestCase):
    def test_selected_enlarged_seam_is_certified(self):
        certificate = certify(
            ARB, smooth_cutoff=F(1, 32), q_star=F(1, 2254),
            interior_q_min=F(1, 4096), tube_radius=F(1, 1701))
        self.assertTrue(certificate.certified())
        self.assertEqual(certificate.tube_radius, F(1, 1701))

    def test_exact_geometry_inequality_holds(self):
        self.assertTrue(seam_geometry_holds(
            ARB, F(1, 32), F(1, 2254), F(1, 1701)))

    def test_q_star_above_the_certified_ceiling_is_rejected(self):
        with self.assertRaises(AssertionError):
            certify(ARB, smooth_cutoff=F(1, 32), q_star=F(1, 2200),
                    interior_q_min=F(1, 4096), tube_radius=F(1, 1728))

    def test_inflated_tube_radius_is_rejected(self):
        with self.assertRaises(AssertionError):
            certify(ARB, smooth_cutoff=F(1, 32), q_star=F(1, 2254),
                    interior_q_min=F(1, 4096), tube_radius=F(1, 1600))

    def test_endpoint_and_interior_ranges_must_overlap(self):
        self.assertTrue(regime_cover_is_complete(F(1, 2256), F(1, 4096)))
        self.assertFalse(regime_cover_is_complete(F(1, 4096), F(1, 1024)))


if __name__ == "__main__":
    unittest.main(verbosity=2)
