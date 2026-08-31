#!/usr/bin/env python3
"""Fail-closed tests for the composed Liu H2 local tube certificate."""

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

HERE = Path(__file__).resolve().parent
UC = HERE.parent
if str(UC) not in sys.path:
    sys.path.insert(0, str(UC))

from liu9_tube import (  # noqa: E402
    PIECEWISE_REPORTS,
    _canonical_json_digest,
    _load_verified_report,
    certify_piecewise_tube,
    load_global_blocker,
)


class PiecewiseTubeTest(unittest.TestCase):
    def test_all_component_digests_replay(self):
        for name in PIECEWISE_REPORTS:
            report = _load_verified_report(name)
            self.assertEqual(len(report["report_sha256"]), 64)

    def test_composed_radius_is_the_q_seam(self):
        report = certify_piecewise_tube()
        self.assertEqual(report["claim_status"], "PROVED")
        self.assertEqual(report["tube_radius"], "1/1701")
        self.assertEqual(report["bottleneck"],
                         "q-degenerate seam at rho=1/1701")

    def test_every_named_region_is_covered(self):
        report = certify_piecewise_tube()
        self.assertTrue(all(report["coverage"].values()))
        self.assertEqual(set(report["coverage"]), {
            "interior_smooth_chart", "zero_support_boundary",
            "mirror_support_boundary", "q_degenerate_endpoints",
            "simultaneous_second_order_insertions", "strict_mean_half_space",
            "q_seam_complete", "quotient_gauges_removed",
        })

    def test_mass_half_is_not_a_geometric_radius(self):
        """The proof uses atom relabelling; it must not report rho=1/2."""
        report = certify_piecewise_tube()
        self.assertEqual(report["mass_chart_radius"], "1/2")
        self.assertNotEqual(report["tube_radius"], report["mass_chart_radius"])
        self.assertIn("smaller fragment", report["canonical_mass_argument"])

    def test_stored_piecewise_digest_is_canonical(self):
        path = UC / "verification/results/liu9-piecewise-tube.json"
        report = json.loads(path.read_text())
        recorded = report.pop("report_sha256")
        self.assertEqual(recorded, _canonical_json_digest(report))

    def test_global_blocker_report_is_open_and_digest_checked(self):
        report = load_global_blocker()
        self.assertEqual(report["claim_status"], "OPEN")
        self.assertFalse(report["negative_mean_feasible_raw_gap_found"])
        self.assertEqual(
            report["endpoint_support"]["claim_status"], "PROVED")
        self.assertEqual(report["qone_kernel"]["claim_status"], "PROVED")
        self.assertEqual(
            report["block_kernel"]["claim_status"], "OPEN_REDUCTION")

    def test_local_report_does_not_overclaim_global_hypothesis(self):
        report = certify_piecewise_tube()
        self.assertIn("complement", report["not_claimed"])
        self.assertIn("still", report["not_claimed"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
