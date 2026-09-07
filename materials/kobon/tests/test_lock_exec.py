#!/usr/bin/env python3
"""Cross-process exclusion tests for the permanent lock executor."""

from __future__ import annotations

import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
DECIDE = HERE.parent / "n12_decide"
LOCK_EXEC = DECIDE / "lock_exec.py"
sys.path.insert(0, str(DECIDE))
import lock_exec as lock_exec_module  # noqa: E402

class LockExecTests(unittest.TestCase):
    def invoke(self, lock: Path, *command: str, nonblocking: bool = False):
        argv = [sys.executable, "-I", "-B", str(LOCK_EXEC)]
        if nonblocking:
            argv.append("--nonblocking")
        argv.extend([str(lock), "--", *command])
        return subprocess.run(argv, capture_output=True, text=True, timeout=5)

    def test_second_process_is_refused_and_lock_file_persists(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            lock = Path(temporary) / "permanent.lock"
            holder = subprocess.Popen([
                sys.executable,
                "-I",
                "-B",
                str(LOCK_EXEC),
                str(lock),
                "--",
                sys.executable,
                "-c",
                "import time; time.sleep(0.5)",
            ])
            try:
                deadline = time.monotonic() + 2
                while not lock.exists() and time.monotonic() < deadline:
                    time.sleep(0.01)
                refused = self.invoke(
                    lock, sys.executable, "-c", "pass", nonblocking=True)
                self.assertEqual(refused.returncode, 75)
                self.assertIn("REFUSED: lock held", refused.stderr)
            finally:
                self.assertEqual(holder.wait(timeout=2), 0)
            accepted = self.invoke(
                lock, sys.executable, "-c", "pass", nonblocking=True)
            self.assertEqual(accepted.returncode, 0, accepted.stderr)
            self.assertTrue(lock.is_file())

    def test_child_keeps_lock_after_wrapper_sigkill(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            lock = root / "inherited.lock"
            ready = root / "child.ready"
            holder = subprocess.Popen([
                sys.executable,
                "-I",
                "-B",
                str(LOCK_EXEC),
                str(lock),
                "--",
                sys.executable,
                "-c",
                (
                    "from pathlib import Path; import sys,time; "
                    "Path(sys.argv[1]).write_text('ready'); time.sleep(0.6)"
                ),
                str(ready),
            ])
            deadline = time.monotonic() + 2
            while not ready.exists() and time.monotonic() < deadline:
                time.sleep(0.01)
            self.assertTrue(ready.is_file())
            holder.kill()
            holder.wait(timeout=2)
            refused = self.invoke(
                lock, sys.executable, "-c", "pass", nonblocking=True)
            self.assertEqual(refused.returncode, 75)
            time.sleep(0.7)
            accepted = self.invoke(
                lock, sys.executable, "-c", "pass", nonblocking=True)
            self.assertEqual(accepted.returncode, 0, accepted.stderr)

    def test_terminate_group_escalates_for_ignoring_child(self) -> None:
        child = subprocess.Popen([
            sys.executable,
            "-c",
            (
                "import signal,time; "
                "signal.signal(signal.SIGTERM, signal.SIG_IGN); "
                "time.sleep(30)"
            ),
        ], start_new_session=True)
        try:
            time.sleep(0.05)
            started = time.monotonic()
            lock_exec_module._terminate_group(child, grace_seconds=0.05)
            self.assertIsNotNone(child.poll())
            self.assertLess(time.monotonic() - started, 1.0)
        finally:
            if child.poll() is None:
                child.kill()
                child.wait()


if __name__ == "__main__":
    unittest.main()
