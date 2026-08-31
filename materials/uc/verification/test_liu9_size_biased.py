#!/usr/bin/env python3
"""Mutation tests for the universal size-biased q-one kernel theorem."""

from __future__ import annotations

import hashlib
import json
import os
import sys
import unittest
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

from liu9_binding import certify_equation_parameters, solve_equation_parameters
from liu9_boundary_layer import gap_mp
from liu9_size_biased import (
    certify_optimizer_box, direct_kernel_range, kernel_range,
)

REPORT = HERE / "results/liu9-size-biased-qone-kernel.json"
MP = solve_equation_parameters(100)
PARAMETERS = certify_equation_parameters(MP)


def _digest(payload):
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode()).hexdigest()


class SizeBiasedKernelTest(unittest.TestCase):
    def test_report_is_authenticated_and_complete(self):
        report = json.loads(REPORT.read_text())
        claimed = report.pop("report_sha256")
        self.assertEqual(_digest(report), claimed)
        self.assertEqual(report["claim_status"], "PROVED")
        self.assertEqual(report["bulk_cells_per_axis"], 1024)
        self.assertGreater(report["optimizer"]["determinant_lower"], 0)
        self.assertIn("small_support", report["analytic_lemmas"])
        self.assertIn("near_one_tail", report["analytic_lemmas"])

    def test_exact_size_biased_identity_matches_gap_mp(self):
        with mpmath.workdps(100):
            masses = (mpmath.mpf("0.2"), mpmath.mpf("0.3"),
                      mpmath.mpf("0.5"))
            points = (mpmath.mpf("0.1"), mpmath.mpf("0.6"),
                      mpmath.mpf("0.9"))
            mean = sum(a * x for a, x in zip(masses, points))
            values = (
                masses[0], masses[1], mpmath.mpf(1),
                mpmath.mpf("0.5"), mpmath.mpf("0.5"), mpmath.mpf("0.5"),
                *points,
            )
            exact = gap_mp(values, MP.beta)
            total = mpmath.mpf(0)
            for i, x in enumerate(points):
                for j, y in enumerate(points):
                    interval, _ = kernel_range(
                        arb(str(x)), arb(str(y)), arb(str(mean)),
                        PARAMETERS.beta)
                    total += masses[i] * masses[j] * mpmath.mpf(
                        interval.mid().str(100).split()[0].lstrip("["))
            reconstructed = total / mean
        self.assertLess(abs(exact - reconstructed), mpmath.mpf("1e-60"))

    def test_protocol_term_is_load_bearing(self):
        x = arb("0.616")
        with_protocol = direct_kernel_range(
            x, x, PARAMETERS.mean, PARAMETERS.beta)
        without_protocol = direct_kernel_range(
            x, x, PARAMETERS.mean, arb(0))
        self.assertGreater(with_protocol, 0)
        self.assertLess(without_protocol, 0)

    def test_optimizer_hessian_mutation_guard(self):
        certificate = certify_optimizer_box(32)
        self.assertGreater(certificate["determinant_lower"], 0.05)
        self.assertIn("+/-", certificate["center_value"])
        self.assertIn("+/-", certificate["center_gradient"][0])


if __name__ == "__main__":
    unittest.main(verbosity=2)
