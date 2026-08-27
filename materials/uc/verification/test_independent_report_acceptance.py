#!/usr/bin/env python3
"""Fail-closed tests for secure replay report acceptance."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import unittest

HERE = Path(__file__).resolve().parent
LOCK_PATH = HERE / "campaign-lock.json"
SECURE_DIR = HERE / "results/independent-arithmetic/secure-full"
PREFLIGHT = HERE / "results/independent-arithmetic/secure-preflight-20k.json"
LAUNCH_MANIFEST = SECURE_DIR / "launch-manifest.json"
LIVE_ATTESTATION = SECURE_DIR / "live-runtime-midrun.json"
REPORT_LOCK = SECURE_DIR / "report-lock.json"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load %s" % path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class IndependentAcceptanceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.aggregate = load(
            "uc_aggregate_test", HERE / "aggregate_independent_reports.py"
        )
        cls.locker = load(
            "uc_report_locker_test", HERE / "lock_independent_reports.py"
        )
        cls.lock = json.loads(LOCK_PATH.read_text())
        cls.lock_raw = LOCK_PATH.read_bytes()

    def test_pinned_source_and_lock_digests_match(self):
        verifier = HERE / "independent_arithmetic_replay_secure.py"
        self.assertEqual(
            hashlib.sha256(verifier.read_bytes()).hexdigest(),
            self.aggregate.EXPECTED_VERIFIER_SHA256,
        )
        self.assertEqual(
            hashlib.sha256(self.lock_raw).hexdigest(),
            self.aggregate.EXPECTED_LOCK_SHA256,
        )

    def test_partial_report_cannot_be_accepted_as_slice_pass(self):
        raw = PREFLIGHT.read_bytes()
        expected = self.lock["slices"][0]["expected_replay_tallies"]
        with self.assertRaisesRegex(
            self.aggregate.VerificationError,
            "not PASS",
        ):
            self.aggregate.validate_report(
                PREFLIGHT,
                {
                    "size": len(raw),
                    "raw_sha256": hashlib.sha256(raw).hexdigest(),
                },
                0,
                self.lock["slices"][0]["run_id"],
                expected,
                self.lock["slices"][0]["trace"]["size"],
                self.aggregate.EXPECTED_CAMPAIGN_ID,
                self.lock["slices"][0]["trace"]["sha256"],
                self.lock["launch"]["raw_sha256"],
                self.lock["slices"][0]["result"]["raw_sha256"],
                self.aggregate.EXPECTED_VERIFIER_SHA256,
                self.aggregate.EXPECTED_LOCK_SHA256,
                {},
            )

    def test_completed_reports_lock_and_aggregate(self):
        generated = self.locker.create_lock(
            SECURE_DIR,
            LAUNCH_MANIFEST,
            LIVE_ATTESTATION,
        )
        self.assertEqual(len(generated["reports"]), 8)
        raw_digest = hashlib.sha256(REPORT_LOCK.read_bytes()).hexdigest()
        aggregate = self.aggregate.aggregate_reports(
            SECURE_DIR,
            LOCK_PATH,
            REPORT_LOCK,
            raw_digest,
        )
        self.assertEqual(aggregate["outcome"], "PASS")
        self.assertEqual(aggregate["total_processed"], 488_465_854)


if __name__ == "__main__":
    unittest.main()
