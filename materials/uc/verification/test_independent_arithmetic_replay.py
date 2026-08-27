#!/usr/bin/env python3
"""Behavioral tests for the independent UC interval replay implementation."""

from __future__ import annotations

import ast
from fractions import Fraction
import importlib.util
import os
from pathlib import Path
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
MATH = HERE.parent.parent
CAMPAIGN = (
    MATH
    / "uc/campaigns/"
    / "cert3_20260818T212601Z_425f109c15b64a6198785c6cebbbdaab_2f23a58ebdb8"
)
LOCK = HERE / "campaign-lock.json"
MODULE_PATH = HERE / "independent_arithmetic_replay.py"


def persistent_temp_path(suffix):
    descriptor, name = tempfile.mkstemp(prefix="uc-independent-test-", suffix=suffix)
    os.close(descriptor)
    return Path(name)

def load_module():
    spec = importlib.util.spec_from_file_location("independent_arithmetic_replay", MODULE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load independent verifier")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class IndependentArithmeticReplayTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_module()

    def test_source_imports_no_frozen_certifier_module(self):
        tree = ast.parse(MODULE_PATH.read_text(encoding="utf-8"))
        forbidden = {
            "arbcore",
            "bound_kkt",
            "cert2",
            "cert3",
            "cert3_replay",
            "cert3_par",
            "diag_exhaust",
            "entropy",
        }
        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(alias.name.split(".", 1)[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module.split(".", 1)[0])
        self.assertTrue(forbidden.isdisjoint(imported), imported & forbidden)

    def test_exact_mean_contractor_preserves_decimal_outer_cutoff(self):
        box = ((1.0, 1.0), (1.0, 1.0), (0.0, 0.0), (0.0, 0.0), (0.0, 1.0))
        contracted = self.mod.contract_mean(box, Fraction(3, 8))
        self.assertIsNotNone(contracted)
        self.assertLessEqual(Fraction(repr(contracted[4][1])), Fraction(3, 8))
        self.assertGreaterEqual(
            Fraction(repr(self.mod.next_decimal_up(Fraction(3, 8)))),
            Fraction(3, 8),
        )

    def test_campaign_target_upper_endpoint_matches_frozen_geometry(self):
        observed = self.mod.target_upper_fraction("0.3820660112501052")
        expected = Fraction(115472366449356389981023, 2**78)
        self.assertEqual(observed, expected)

    def test_all_campaign_atom_dyadics_have_exact_decimal_repr(self):
        for depth in range(14):
            denominator = 2**depth
            for numerator in range(denominator + 1):
                value = numerator / denominator
                self.assertEqual(
                    Fraction(repr(value)),
                    Fraction(numerator, denominator),
                )

    def test_interval_ad_encloses_weight_derivative(self):
        box = (
            (0.30, 0.31),
            (0.34, 0.35),
            (0.20, 0.21),
            (0.70, 0.71),
            (0.60, 0.61),
        )
        gradient = self.mod.ad_gradient(box, "0.3820660112501052", "1.0")
        self.assertIsNotNone(gradient)
        self.assertEqual(len(gradient), 5)
        self.assertTrue(all(component.is_finite() for component in gradient))

    def test_touching_sstar_branch_keeps_mixed_bound_positive(self):
        box = (
            (0.5546875, 0.5625),
            (0.5625, 0.5703125),
            (0.0, 0.015625),
            (0.296875, 0.3125),
            (0.5078125, 0.51171875),
        )
        rh_cap, _, _ = self.mod._ensure_global_caps()
        previous = self.mod.ctx.prec
        try:
            self.mod.ctx.prec = self.mod.REPLAY_PREC_BITS
            bound = self.mod._best_centered(
                box,
                self.mod._a("0.3820660112501052"),
                tuple(
                    self.mod._a(value)
                    for value in self.mod.PINNED_CENTER_LAMBDAS
                ),
                (0, 1, 4),
                rh_cap,
            )
        finally:
            self.mod.ctx.prec = previous
        self.assertIsNotNone(bound)
        self.assertGreaterEqual(bound, self.mod.ZERO)

    def test_false_mean_infeasible_event_fails_closed(self):
        root = ((0.0, 1.0),) * 4 + ((0.5, 9.0 / 16.0),)
        trace = persistent_temp_path(".bin")
        trace.write_bytes(bytes((self.mod.EVENT_MEAN_INFEASIBLE,)))
        with self.assertRaises(self.mod.VerificationError):
            self.mod.replay_trace(
                trace,
                root,
                "0.3820660112501052",
                ["0.0", "1.0"],
                face_min_width=0.001,
                face_box_budget=100000,
            )
    
    def test_unknown_and_residual_opcodes_fail_closed(self):
        root = ((0.0, 1.0),) * 4 + ((0.5, 9.0 / 16.0),)
        for event in (15, self.mod.EVENT_RESIDUAL):
            with self.subTest(event=event):
                trace = persistent_temp_path(".bin")
                trace.write_bytes(bytes((event,)))
                with self.assertRaises(self.mod.VerificationError):
                    self.mod.replay_trace(
                        trace,
                        root,
                        "0.3820660112501052",
                        ["0.0", "1.0"],
                        face_min_width=0.001,
                        face_box_budget=100000,
                    )

    def test_real_slice_prefix_recomputes_without_frozen_code(self):
        report = self.mod.audit_campaign_slice(
            LOCK,
            CAMPAIGN,
            0,
            max_events=20_000,
            checkpoint_path=None,
            resume=False,
        )
        self.assertEqual(report["outcome"], "PARTIAL_AUDIT")
        self.assertEqual(report["processed"], 20_000)
        self.assertGreater(report["pending_stack"], 0)

    def test_authenticated_checkpoint_resumes_and_rejects_tampering(self):
        checkpoint = persistent_temp_path(".json")
        first = self.mod.audit_campaign_slice(
            LOCK,
            CAMPAIGN,
            0,
            max_events=1_000,
            checkpoint_path=checkpoint,
            resume=False,
        )
        self.assertEqual(first["processed"], 1_000)
        second = self.mod.audit_campaign_slice(
            LOCK,
            CAMPAIGN,
            0,
            max_events=2_000,
            checkpoint_path=checkpoint,
            resume=True,
        )
        self.assertEqual(second["processed"], 2_000)
        checkpoint.write_bytes(checkpoint.read_bytes() + b"tamper")
        with self.assertRaises(self.mod.VerificationError):
            self.mod.audit_campaign_slice(
                LOCK,
                CAMPAIGN,
                0,
                max_events=3_000,
                checkpoint_path=checkpoint,
                resume=True,
            )


if __name__ == "__main__":
    unittest.main()
