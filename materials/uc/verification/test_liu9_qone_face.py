#!/usr/bin/env python3
"""Mutation and fail-closed tests for the q=1 active-face certifier."""

from __future__ import annotations

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

from liu9_binding import solve_equation_parameters  # noqa: E402
from liu9_boundary_layer import gap_mp  # noqa: E402
from liu9_pilot import GradientUnavailable, Jet, _JetOps  # noqa: E402
from liu9_qone_face import (  # noqa: E402
    SOURCE_DEFAULT,
    ZERO_LAYER,
    canonical_qone_box,
    feasible,
    initial_active_boxes,
    reduce_simplex_face,
    reduce_zero_layer,
    run,
)

MP = solve_equation_parameters(100)


class QOneProjectionTest(unittest.TestCase):
    def test_source_reduces_to_230_active_boxes(self):
        boxes, meta = initial_active_boxes(SOURCE_DEFAULT)
        self.assertEqual(meta["q_at_one_source_boxes"], 48996)
        self.assertEqual(meta["unique_active_boxes"], 230)
        self.assertEqual(len(boxes), 230)

    def test_inactive_supports_are_gauges_only_at_exact_q_one(self):
        with mpmath.workdps(70):
            active = (MP.p, 0, 1, 0.1, 0.2, 0.3,
                      MP.x, mpmath.mpf(1) / 3, 0)
            changed = list(active)
            changed[3:6] = [0.9, 0.8, 0.7]
            self.assertEqual(gap_mp(active, MP.beta), gap_mp(changed, MP.beta))
            changed[2] = mpmath.mpf(99) / 100
            self.assertNotEqual(gap_mp(active, MP.beta), gap_mp(changed, MP.beta))

    def test_forced_simplex_face_becomes_exact_and_drops_b5(self):
        box = canonical_qone_box((
            (5/16, 11/32), (11/16, 23/32),
            (5/32, 3/16), (1/4, 9/32), (31/32, 1),
        ))
        reduced, applied = reduce_simplex_face(box)
        self.assertTrue(applied)
        self.assertEqual(reduced[0], (5/16, 5/16))
        self.assertEqual(reduced[1], (11/16, 11/16))
        self.assertEqual(reduced[8], (0.5, 0.5))

    def test_nonbinding_simplex_box_is_not_reduced(self):
        box = canonical_qone_box((
            (0, 1/4), (0, 1/4), (0.2, 0.3), (0.3, 0.4), (0.5, 0.6)))
        unchanged, applied = reduce_simplex_face(box)
        self.assertFalse(applied)
        self.assertEqual(unchanged, box)

    def test_zero_layer_reduces_only_wholly_covered_intervals(self):
        base = canonical_qone_box((
            (0, 0.5), (0, 0.5),
            (0, float(ZERO_LAYER)), (0.25, 0.5), (0.5, 1),
        ))
        reduced, count = reduce_zero_layer(base)
        self.assertEqual(count, 1)
        self.assertEqual(reduced[6], (0.0, 0.0))
        wider = list(base)
        wider[6] = (0.0, float(ZERO_LAYER) * 2)
        unchanged, count = reduce_zero_layer(tuple(wider))
        self.assertEqual(count, 0)
        self.assertEqual(unchanged[6], wider[6])


class EndpointJetTest(unittest.TestCase):
    def test_fixed_endpoint_entropy_is_exact_constant(self):
        jet = Jet(arb(0), (arb(0), arb(0)))
        result = _JetOps._entropy(jet)
        self.assertEqual(result.value, 0)
        self.assertTrue(all(value == 0 for value in result.gradient))

    def test_nonzero_endpoint_direction_still_fails(self):
        jet = Jet(arb(0), (arb(1), arb(0)))
        with self.assertRaises(GradientUnavailable):
            _JetOps._entropy(jet)


class QOneRunTest(unittest.TestCase):
    def test_mean_infeasible_cell_is_rejected_before_bounding(self):
        from liu9_binding import certify_equation_parameters
        parameters = certify_equation_parameters(MP)
        box = canonical_qone_box((
            (0, 1/32), (3/4, 25/32),
            (1/32, 1/16), (5/16, 3/8), (9/16, 5/8),
        ))
        self.assertFalse(feasible(box, parameters.mean))


    def test_short_run_has_sane_depth_and_no_false_proof(self):
        report = run(SOURCE_DEFAULT, budget=501, seconds=20, checkpoint=None)
        self.assertEqual(report["source"]["unique_active_boxes"], 230)
        self.assertLessEqual(report["run"]["maximum_depth"], 10)
        self.assertGreater(report["run"]["residual_count"], 0)
        self.assertEqual(report["claim_status"], "NUMERICAL")

    def test_known_ambient_negative_is_mean_infeasible(self):
        with mpmath.workdps(70):
            y = mpmath.mpf(1) / 32
            eps = mpmath.mpf(1) / 1024
            values = (MP.p - eps, eps, mpmath.mpf(1) / 2,
                      MP.x, y, 0, MP.x, y, 0)
            raw = gap_mp(values, MP.beta)
            mean = (MP.p - eps) * MP.x + eps * y
            self.assertLess(float(raw), -5e-4)
            self.assertLess(mean, MP.mean)


if __name__ == "__main__":
    unittest.main(verbosity=2)
