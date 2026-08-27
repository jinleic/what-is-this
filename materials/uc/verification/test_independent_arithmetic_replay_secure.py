#!/usr/bin/env python3
"""Adversarial contracts for the byte-zero UC independent verifier."""

from __future__ import annotations

import ast
from fractions import Fraction
import importlib.util
import json
import math
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
SEAL = HERE / "results/runtime-environment-seal-secure.json"
MODULE_PATH = HERE / "independent_arithmetic_replay_secure.py"


def load_module():
    spec = importlib.util.spec_from_file_location("uc_secure_replay_test", MODULE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load secure independent verifier")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def persistent_temp_path(suffix):
    descriptor, name = tempfile.mkstemp(prefix="uc-secure-test-", suffix=suffix)
    os.close(descriptor)
    return Path(name)


def entropy(value):
    if value <= 0.0 or value >= 1.0:
        return 0.0
    return -value * math.log2(value) - (1.0 - value) * math.log2(1.0 - value)


def sstar(left, right):
    return min(max(0.5, left, right), min(left + right, 1.0))


def sqrt_rh(value):
    rh = 2.0 * (1.0 - value) * entropy(value) - entropy((1.0 - value) ** 2)
    return math.sqrt(max(rh, 0.0))


def phi(point, alpha=0.0356069, target=0.3820660112501052, lam=1.0):
    p1, q1, p2, q2, weight = point
    other = 1.0 - weight
    mean1, mean2 = (p1 + q1) / 2.0, (p2 + q2) / 2.0
    entropy1 = (entropy(p1) + entropy(q1)) / 2.0
    entropy2 = (entropy(p2) + entropy(q2)) / 2.0
    correction1 = entropy(sstar(p1, q1))
    correction2 = entropy(sstar(p2, q2))
    sigma1 = (sqrt_rh(p1) + sqrt_rh(q1)) / 2.0
    sigma2 = (sqrt_rh(p2) + sqrt_rh(q2)) / 2.0
    mean = weight * mean1 + other * mean2
    total_entropy = weight * entropy1 + other * entropy2
    correction = weight * correction1 + other * correction2
    sigma = weight * sigma1 + other * sigma2
    beta = 1.0 - alpha
    return (
        (2.0 * beta * (1.0 - mean) - 1.0) * total_entropy
        + alpha * correction
        - beta * sigma * sigma
        + lam * (mean - target)
    )


class SecureReplayTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_module()

    def test_no_frozen_arithmetic_import(self):
        tree = ast.parse(MODULE_PATH.read_text(encoding="utf-8"))
        forbidden = {
            "arbcore",
            "bound_kkt",
            "cert2",
            "cert3",
            "cert3_par",
            "cert3_replay",
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

    def test_target_geometry_is_exactly_frozen(self):
        self.assertEqual(
            self.mod.target_upper_fraction("0.3820660112501052"),
            Fraction(115472366449356389981023, 2**78),
        )

    def test_ad_contains_independent_finite_differences(self):
        points = (
            (0.18, 0.27, 0.55, 0.72, 0.63),
            (0.22, 0.31, 0.60, 0.72, 0.40),
        )
        step = 1.0e-7
        for point in points:
            gradient = self.mod.ad_gradient(
                tuple((value, value) for value in point),
                "0.3820660112501052",
                "1.0",
            )
            self.assertIsNotNone(gradient)
            for coordinate, enclosure in enumerate(gradient):
                plus = list(point)
                minus = list(point)
                plus[coordinate] += step
                minus[coordinate] -= step
                numerical = (phi(plus) - phi(minus)) / (2.0 * step)
                self.assertLessEqual(float(enclosure.lower()) - 2.0e-6, numerical)
                self.assertGreaterEqual(float(enclosure.upper()) + 2.0e-6, numerical)

    def test_runtime_seal_is_bound_and_validated(self):
        record = self.mod._validate_runtime_seal(SEAL)
        self.assertEqual(record["file_count"], 343)
        forged = json.loads(SEAL.read_text())
        forged["file_count"] += 1
        forged_path = persistent_temp_path(".json")
        forged_path.write_text(json.dumps(forged, separators=(",", ":")) + "\n")
        with self.assertRaises(self.mod.VerificationError):
            self.mod._validate_runtime_seal(forged_path)

    def test_resume_can_never_produce_full_pass(self):
        with self.assertRaisesRegex(
            self.mod.VerificationError,
            "cannot produce a full PASS",
        ):
            self.mod.audit_campaign_slice(
                LOCK,
                CAMPAIGN,
                0,
                runtime_seal_path=SEAL,
                max_events=None,
                checkpoint_path=persistent_temp_path(".json"),
                resume=True,
            )

    def test_byte_zero_partial_report_discloses_provenance(self):
        report = self.mod.audit_campaign_slice(
            LOCK,
            CAMPAIGN,
            0,
            runtime_seal_path=SEAL,
            max_events=20_000,
            checkpoint_path=None,
            resume=False,
        )
        self.assertEqual(report["outcome"], "PARTIAL_AUDIT")
        self.assertEqual(report["initial_trace_offset"], 0)
        self.assertIs(report["resume_used"], False)
        self.assertEqual(report["runtime_seal"]["file_count"], 343)


if __name__ == "__main__":
    unittest.main()
