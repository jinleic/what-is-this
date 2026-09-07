"""Regression tests for fixed-five certificate bundle closure."""

from __future__ import annotations

import copy
import json
import sys
import unittest
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import check_involution_f5_certificate_coverage as coverage  # noqa: E402

_ROOT = Path(__file__).resolve().parents[1]
_W_BATCH = _ROOT / "data" / "involution_f5_w_dfs_certificates.json"


class TestSourceZeroArchiveSchema(unittest.TestCase):
    def setUp(self):
        document = json.loads(_W_BATCH.read_text())
        self.record = document["records"][0]
        self.assertIn("seed_manifest", self.record)

    def test_source_zero_archives_pass_the_closure_checker(self):
        self.assertEqual(coverage._check_w_record(self.record), 0)

    def test_source_zero_corrupt_gzip_digest_is_rejected(self):
        record = copy.deepcopy(self.record)
        record["veripb_proof"]["gzip_sha256"] = "0" * 64
        with self.assertRaisesRegex(
                coverage.CoverageViolation, "proof archive hash mismatch"):
            coverage._check_w_record(record)

    def test_source_zero_archive_path_cannot_escape_data_directory(self):
        record = copy.deepcopy(self.record)
        record["veripb_proof"]["relative_path"] = (
            "r55/data/../../etc/passwd")
        with self.assertRaisesRegex(
                coverage.CoverageViolation, "invalid proof archive path"):
            coverage._check_w_record(record)


class TestSeedSchemaRestriction(unittest.TestCase):
    def setUp(self):
        document = json.loads(_W_BATCH.read_text())
        self.record = document["records"][0]
        self.assertIn("seed_manifest", self.record)

    def test_seed_manifest_on_nonzero_source_is_rejected(self):
        record = copy.deepcopy(self.record)
        record["source_index"] = 5
        with self.assertRaisesRegex(
                coverage.CoverageViolation, "restricted to source zero"):
            coverage._check_w_record(record)

    def test_wrong_seed_manifest_path_is_rejected(self):
        record = copy.deepcopy(self.record)
        record["seed_manifest"] = "r55/data/somewhere_else.json"
        with self.assertRaisesRegex(
                coverage.CoverageViolation, "restricted to source zero"):
            coverage._check_w_record(record)

    def test_noninteger_seed_source_is_rejected(self):
        record = copy.deepcopy(self.record)
        record["source_index"] = "0"
        with self.assertRaisesRegex(
                coverage.CoverageViolation, "invalid W certificate record"):
            coverage._check_w_record(record)

    def test_seed_manifest_in_hierarchical_record_is_rejected(self):
        record = {
            "source_index": 0,
            "status": "VERIPB_CAKEPB_VERIFIED",
            "unchecked_deletion_present": False,
            "seed_manifest": coverage.W_SEED_MANIFEST,
        }
        with self.assertRaisesRegex(
                coverage.CoverageViolation,
                "invalid in hierarchical records"):
            coverage._check_hierarchy_record(record, set(), {0: 1}, {})

    def test_seed_unchecked_deletion_flag_is_rejected(self):
        record = copy.deepcopy(self.record)
        record["unchecked_deletion_present"] = True
        with self.assertRaisesRegex(
                coverage.CoverageViolation, "permits unchecked deletion"):
            coverage._check_w_record(record)


if __name__ == "__main__":
    unittest.main()
