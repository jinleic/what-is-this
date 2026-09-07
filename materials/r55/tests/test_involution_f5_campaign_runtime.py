"""Contract tests for campaign lock files and durable JSON replacement."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import campaign_runtime  # noqa: E402


class TestConfiguredLockPaths(unittest.TestCase):
    def test_driver_scheduler_and_host_locks_are_distinct_and_exact(self):
        scratch = campaign_runtime.WORKSPACE / "scratch/r55-involution-f5"
        self.assertEqual(
            campaign_runtime.DRIVER_LOCK,
            scratch / ".certificate-driver.lock")
        self.assertEqual(
            campaign_runtime.SCHEDULER_LOCK,
            scratch / ".certificate-scheduler.lock")
        self.assertEqual(
            campaign_runtime.HOST_HEAVY_LOCK,
            campaign_runtime.WORKSPACE / "scratch/.host-heavy-job.lock")
        self.assertEqual(
            len({
                campaign_runtime.DRIVER_LOCK,
                campaign_runtime.SCHEDULER_LOCK,
                campaign_runtime.HOST_HEAVY_LOCK,
            }),
            3)


class TestFileLock(unittest.TestCase):
    def test_second_holder_is_refused_until_release(self):
        with tempfile.TemporaryDirectory() as directory:
            lock_path = Path(directory) / "campaign.lock"
            first = campaign_runtime.FileLock(lock_path)
            first.acquire()
            try:
                with self.assertRaises(campaign_runtime.LockBusy):
                    campaign_runtime.FileLock(lock_path).acquire()
            finally:
                first.release()
            second = campaign_runtime.FileLock(lock_path)
            second.acquire()  # must not raise
            second.release()

    def test_context_manager_releases_and_lock_file_persists(self):
        with tempfile.TemporaryDirectory() as directory:
            lock_path = Path(directory) / "campaign.lock"
            with campaign_runtime.FileLock(lock_path) as lock:
                self.assertTrue(lock_path.is_file())
                self.assertIsInstance(lock.fileno(), int)
            with campaign_runtime.FileLock(lock_path):
                pass
            self.assertTrue(lock_path.is_file())
            with self.assertRaisesRegex(RuntimeError, "lock is not held"):
                lock.fileno()

    def test_blocking_mode_uses_exclusive_lock_without_nonblocking_flag(self):
        with tempfile.TemporaryDirectory() as directory:
            lock_path = Path(directory) / "host-heavy.lock"
            with mock.patch.object(campaign_runtime.fcntl, "flock") as flock:
                lock = campaign_runtime.FileLock(lock_path, blocking=True)
                lock.acquire()
                lock.release()
            self.assertEqual(
                flock.call_args_list[0].args[1], campaign_runtime.fcntl.LOCK_EX)


class TestAtomicWriteJson(unittest.TestCase):
    def test_writes_canonical_json_without_temp_remnants(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ledger.json"
            campaign_runtime.atomic_write_json(path, {"b": 1, "a": [1, 2]})
            self.assertEqual(
                path.read_text(),
                json.dumps({"b": 1, "a": [1, 2]}, indent=2, sort_keys=True)
                + "\n")
            self.assertEqual(
                [entry.name for entry in Path(directory).iterdir()],
                ["ledger.json"])

    def test_rewrite_replaces_the_whole_document(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ledger.json"
            campaign_runtime.atomic_write_json(path, {"stale": True})
            campaign_runtime.atomic_write_json(path, {"fresh": 1})
            self.assertEqual(json.loads(path.read_text()), {"fresh": 1})

    def test_serialization_failure_preserves_the_previous_document(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ledger.json"
            campaign_runtime.atomic_write_json(path, {"keep": True})
            temp = Path(directory) / "ledger.json.tmp"
            temp.write_text("stale")
            with self.assertRaises(TypeError):
                campaign_runtime.atomic_write_json(path, {"bad": object()})
            with self.assertRaises(ValueError):
                campaign_runtime.atomic_write_json(
                    path, {"bad": float("nan")})
            self.assertEqual(json.loads(path.read_text()), {"keep": True})
            self.assertFalse(temp.exists())

    def test_load_json_object_rejects_corrupt_payloads(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ledger.json"
            path.write_text("{not json")
            with self.assertRaisesRegex(ValueError, "cannot load"):
                campaign_runtime.load_json_object(path, ValueError)

    def test_load_json_object_rejects_duplicate_members(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ledger.json"
            path.write_text('{"records": [], "records": [1]}')
            with self.assertRaisesRegex(ValueError, "duplicate JSON member"):
                campaign_runtime.load_json_object(path, ValueError)

    def test_load_json_object_rejects_nonstandard_constants(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ledger.json"
            path.write_text('{"coverage": NaN}')
            with self.assertRaisesRegex(
                    ValueError, "nonstandard JSON constant"):
                campaign_runtime.load_json_object(path, ValueError)


if __name__ == "__main__":
    unittest.main()
