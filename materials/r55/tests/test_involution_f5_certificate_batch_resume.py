"""Resume-safety tests for the fixed-five certificate batch drivers."""

from __future__ import annotations

import copy
import json
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import campaign_runtime  # noqa: E402
import certify_involution_f5_hierarchical_batch as h_batch  # noqa: E402
import certify_involution_f5_w_dfs_batch as w_batch  # noqa: E402


def _w_document(directory: str, records: list[dict], target: int = 3) -> dict:
    output = Path(directory) / "w.json"
    w_batch._write_progress(output, records, target, time.monotonic())
    return json.loads(output.read_text())


def _h_document(directory: str, records: list[dict], target: int = 3) -> dict:
    output = Path(directory) / "h.json"
    h_batch._write_progress(output, records, target, time.monotonic())
    return json.loads(output.read_text())


class TestWPartialLedgerValidation(unittest.TestCase):
    def test_valid_seed_only_ledger_delegates_record_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            records = [w_batch._seed_source_zero()]
            document = _w_document(directory, records)
            with mock.patch.object(
                    w_batch.coverage, "_check_w_record", return_value=0
            ) as checker:
                w_batch._validate_partial_ledger(
                    document, records, {0, 1, 2}, 3)
            checker.assert_called_once_with(records[0])
            self.assertFalse((Path(directory) / "w.json.tmp").exists())

    def test_stale_schema_campaign_input_and_toolchain_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            base = _w_document(directory, [])
            mutations = {
                "schema": lambda document: document.__setitem__(
                    "schema_version", True),
                "campaign": lambda document: document.__setitem__(
                    "campaign_id", "foreign_campaign"),
                "input": lambda document: document["input_dependency"].__setitem__(
                    "sha256", "0" * 64),
                "toolchain": lambda document: document["toolchain"].__setitem__(
                    "veripb_binary_sha256", "0" * 64),
            }
            for label, mutate in mutations.items():
                with self.subTest(label=label):
                    document = copy.deepcopy(base)
                    mutate(document)
                    with self.assertRaisesRegex(
                            w_batch.BatchViolation, "stale"):
                        w_batch._validate_partial_ledger(
                            document, [], {0, 1, 2}, 3)

    def test_duplicate_source_index_is_rejected_before_archive_checks(self):
        with tempfile.TemporaryDirectory() as directory:
            records = [{"source_index": 1}, {"source_index": 1}]
            document = _w_document(directory, records)
            with self.assertRaisesRegex(
                    w_batch.BatchViolation, "duplicate record for source 1"):
                w_batch._validate_partial_ledger(
                    document, records, {0, 1, 2}, 3)

    def test_foreign_source_index_is_rejected_before_archive_checks(self):
        with tempfile.TemporaryDirectory() as directory:
            records = [{"source_index": 9}]
            document = _w_document(directory, records)
            with self.assertRaisesRegex(
                    w_batch.BatchViolation, "foreign record for source 9"):
                w_batch._validate_partial_ledger(
                    document, records, {0, 1, 2}, 3)

    def test_noninteger_source_index_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            records = [{"source_index": "1"}]
            document = _w_document(directory, records)
            with self.assertRaisesRegex(
                    w_batch.BatchViolation, "integer source_index"):
                w_batch._validate_partial_ledger(
                    document, records, {0, 1, 2}, 3)

    def test_unverified_status_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            records = [{"source_index": 1, "status": "PENDING"}]
            document = _w_document(directory, records)
            with self.assertRaisesRegex(
                    w_batch.BatchViolation, "invalid W certificate record"):
                w_batch._validate_partial_ledger(
                    document, records, {0, 1, 2}, 3)

    def test_unchecked_deletion_policy_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            records = [{
                "source_index": 1,
                "status": "VERIPB_CAKEPB_VERIFIED",
                "unchecked_deletion_present": True,
            }]
            document = _w_document(directory, records)
            with self.assertRaisesRegex(
                    w_batch.BatchViolation, "permits unchecked deletion"):
                w_batch._validate_partial_ledger(
                    document, records, {0, 1, 2}, 3)

    def test_inconsistent_coverage_and_disposition_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            base = _w_document(directory, [])
            coverage_document = copy.deepcopy(base)
            coverage_document["coverage"]["certified_support_orbits"] = 1
            with self.assertRaisesRegex(
                    w_batch.BatchViolation, "coverage is inconsistent"):
                w_batch._validate_partial_ledger(
                    coverage_document, [], {0, 1, 2}, 3)
            typed_records = [{"source_index": 1}]
            typed_coverage = _w_document(directory, typed_records)
            typed_coverage["coverage"]["certified_support_orbits"] = True
            with self.assertRaisesRegex(
                    w_batch.BatchViolation, "coverage is inconsistent"):
                w_batch._validate_partial_ledger(
                    typed_coverage, typed_records, {0, 1, 2}, 3)
            disposition_document = copy.deepcopy(base)
            disposition_document["disposition"] = (
                "F5_W_NEGATIVE_550_VERIPB_CAKEPB_CERTIFIED")
            with self.assertRaisesRegex(
                    w_batch.BatchViolation, "disposition is inconsistent"):
                w_batch._validate_partial_ledger(
                    disposition_document, [], {0, 1, 2}, 3)


    def test_fresh_seed_archives_are_validated_before_use(self):
        seed = {"source_index": 0}
        with (mock.patch.object(
                w_batch, "_seed_source_zero", return_value=seed),
              mock.patch.object(
                  w_batch.coverage, "_check_w_record",
                  side_effect=w_batch.coverage.CoverageViolation(
                      "missing proof archive")) as checker):
            with self.assertRaisesRegex(
                    w_batch.BatchViolation,
                    "source-zero certificate seed.*missing proof archive"):
                w_batch._validated_seed_source_zero()
        checker.assert_called_once_with(seed)

    def test_exact_legacy_toolchain_migrates_after_record_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "w.json"
            records = [{"source_index": 1}]
            document = _w_document(directory, records)
            document["toolchain"] = w_batch._legacy_toolchain()
            campaign_runtime.atomic_write_json(output, document)
            order = []
            write_json = campaign_runtime.atomic_write_json

            def check_record(_record):
                order.append("record")
                return 1

            def write_migration(path, migrated):
                order.append("write")
                write_json(path, migrated)

            with (mock.patch.object(
                    w_batch.coverage, "_check_w_record",
                    side_effect=check_record),
                  mock.patch.object(
                      w_batch.campaign_runtime, "atomic_write_json",
                      side_effect=write_migration)):
                migrated = w_batch._load_validated_output(
                    output, {0, 1, 2}, 3)

            self.assertEqual(order, ["record", "write"])
            expected = w_batch._expected_toolchain()
            self.assertEqual(migrated["toolchain"], expected)
            persisted = json.loads(output.read_text())
            self.assertEqual(persisted["toolchain"], expected)
            self.assertEqual([
                profile["heap"]
                for profile in persisted["toolchain"][
                    "cakepb_memory_profiles_mib"]
            ], [65536, 90112, 114688])
            self.assertFalse((Path(directory) / "w.json.tmp").exists())

    def test_legacy_toolchain_is_not_rewritten_if_record_validation_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "w.json"
            records = [{"source_index": 1}]
            document = _w_document(directory, records)
            legacy = w_batch._legacy_toolchain()
            document["toolchain"] = legacy
            campaign_runtime.atomic_write_json(output, document)
            with mock.patch.object(
                    w_batch.coverage, "_check_w_record",
                    side_effect=w_batch.coverage.CoverageViolation(
                        "missing proof archive")):
                with self.assertRaisesRegex(
                        w_batch.BatchViolation, "missing proof archive"):
                    w_batch._load_validated_output(
                        output, {0, 1, 2}, 3)
            persisted = json.loads(output.read_text())
            self.assertEqual(persisted["toolchain"], legacy)

    def test_arbitrary_stale_toolchain_is_not_treated_as_predecessor(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "w.json"
            records = [{"source_index": 1}]
            document = _w_document(directory, records)
            document["toolchain"]["cakepb_memory_profiles_mib"] = [
                {"heap": 65536, "stack": 4096},
                {"heap": 114688, "stack": 4096},
            ]
            campaign_runtime.atomic_write_json(output, document)
            with mock.patch.object(
                    w_batch.coverage, "_check_w_record") as checker:
                with self.assertRaisesRegex(
                        w_batch.BatchViolation, "stale toolchain metadata"):
                    w_batch._load_validated_output(
                        output, {0, 1, 2}, 3)
            checker.assert_not_called()


    def test_fresh_seed_counts_against_batch_allowance(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            targets = [{"source_index": 0}, {"source_index": 1}]
            seed = {"source_index": 0}
            with (mock.patch.object(
                    w_batch, "_load_targets", return_value=targets),
                  mock.patch.object(
                      w_batch, "_validated_seed_source_zero",
                      return_value=seed),
                  mock.patch.object(
                      w_batch.completion, "load_instance") as load_instance):
                document = w_batch._run_locked(
                    root / "w.json", root / "proofs", root / "work",
                    1.0, max_new_records=1)
            self.assertEqual(document["records"], [seed])
            load_instance.assert_not_called()


class TestDriverLockInheritance(unittest.TestCase):
    def test_w_run_passes_held_driver_descriptor(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with (
                mock.patch.object(
                    w_batch.campaign_runtime, "DRIVER_LOCK",
                    root / "driver.lock"),
                mock.patch.object(
                    w_batch, "_run_locked", return_value={}) as run_locked,
            ):
                w_batch.run(
                    root / "w.json", root / "proofs", root / "work", 1.0)
            descriptor = run_locked.call_args.kwargs["driver_lock_fd"]
            self.assertIs(type(descriptor), int)

    def test_hierarchy_run_passes_held_driver_descriptor(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with (
                mock.patch.object(
                    h_batch.campaign_runtime, "DRIVER_LOCK",
                    root / "driver.lock"),
                mock.patch.object(
                    h_batch, "_run_locked", return_value={}) as run_locked,
            ):
                h_batch.run(
                    root / "h.json", root / "proofs", root / "work", 1.0)
            descriptor = run_locked.call_args.kwargs["driver_lock_fd"]
            self.assertIs(type(descriptor), int)


class TestCalibratedCakePbProfiles(unittest.TestCase):
    def test_w_large_formula_uses_calibrated_profile_last(self):
        self.assertEqual(
            w_batch._cakepb_profiles(350_000_001),
            (114688,))
        self.assertEqual(
            w_batch._cakepb_profiles(250_000_001),
            (90112, 114688))
        self.assertEqual(
            w_batch._cakepb_profiles(250_000_000),
            (65536, 90112, 114688))

    def test_hierarchy_large_formula_uses_calibrated_profile_last(self):
        self.assertEqual(
            h_batch._cakepb_profiles(350_000_001, True),
            (114688,))
        self.assertEqual(
            h_batch._cakepb_profiles(250_000_001, True),
            (90112, 114688))
        self.assertEqual(
            h_batch._cakepb_profiles(250_000_000, True),
            (65536, 90112, 114688))
        self.assertEqual(
            h_batch._cakepb_profiles(500_000_000, False),
            (8192, 32768, 65536, 90112, 114688))


class TestHierarchicalPartialLedgerValidation(unittest.TestCase):
    def test_empty_ledger_with_current_metadata_is_valid(self):
        with tempfile.TemporaryDirectory() as directory:
            document = _h_document(directory, [])
            h_batch._validate_partial_ledger(
                document, [], {0, 1, 2}, 3, set(), {}, {})
            self.assertFalse((Path(directory) / "h.json.tmp").exists())

    def test_foreign_campaign_and_stale_toolchain_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            base = _h_document(directory, [])
            campaign = copy.deepcopy(base)
            campaign["campaign_id"] = "involution_f5_w_dfs_certificates"
            with self.assertRaisesRegex(
                    h_batch.HierarchyViolation, "stale campaign metadata"):
                h_batch._validate_partial_ledger(
                    campaign, [], {0, 1, 2}, 3, set(), {}, {})
            toolchain = copy.deepcopy(base)
            toolchain["toolchain"]["veripb_commit"] = "foreign"
            with self.assertRaisesRegex(
                    h_batch.HierarchyViolation, "stale toolchain metadata"):
                h_batch._validate_partial_ledger(
                    toolchain, [], {0, 1, 2}, 3, set(), {}, {})

    def test_duplicate_and_foreign_sources_are_rejected_structurally(self):
        with tempfile.TemporaryDirectory() as directory:
            duplicate = [{"source_index": 1}, {"source_index": 1}]
            duplicate_document = _h_document(directory, duplicate)
            with self.assertRaisesRegex(
                    h_batch.HierarchyViolation,
                    "duplicate record for source 1"):
                h_batch._validate_partial_ledger(
                    duplicate_document, duplicate, {0, 1, 2}, 3,
                    set(), {}, {})
            foreign = [{"source_index": 9}]
            foreign_document = _h_document(directory, foreign)
            with self.assertRaisesRegex(
                    h_batch.HierarchyViolation,
                    "foreign record for source 9"):
                h_batch._validate_partial_ledger(
                    foreign_document, foreign, {0, 1, 2}, 3,
                    set(), {}, {})

    def test_record_validation_is_delegated_to_coverage_checker(self):
        with tempfile.TemporaryDirectory() as directory:
            records = [{"source_index": 1}]
            document = _h_document(directory, records)
            with mock.patch.object(
                    h_batch.coverage, "_check_hierarchy_record",
                    return_value=1
            ) as checker:
                h_batch._validate_partial_ledger(
                    document, records, {0, 1, 2}, 3,
                    set(), {1: 0}, {})
            checker.assert_called_once_with(records[0], set(), {1: 0}, {})

    def test_inconsistent_coverage_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            document = _h_document(directory, [])
            document["coverage"]["remaining_support_orbits"] = 99
            with self.assertRaisesRegex(
                    h_batch.HierarchyViolation, "coverage is inconsistent"):
                h_batch._validate_partial_ledger(
                    document, [], {0, 1, 2}, 3, set(), {}, {})


class TestDriverLock(unittest.TestCase):
    def test_w_and_h_drivers_share_one_nonblocking_lock(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            lock = root / "driver.lock"
            with mock.patch.object(campaign_runtime, "DRIVER_LOCK", lock):
                with campaign_runtime.FileLock(lock):
                    with self.assertRaises(campaign_runtime.LockBusy):
                        w_batch.run(
                            root / "w.json", root / "proofs", root / "work",
                            1.0)
                    with self.assertRaises(campaign_runtime.LockBusy):
                        h_batch.run(
                            root / "h.json", root / "proofs", root / "work",
                            1.0)
            self.assertTrue(lock.is_file())


if __name__ == "__main__":
    unittest.main()
