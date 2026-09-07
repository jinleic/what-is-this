#!/usr/bin/env python3
"""Small real-process regressions for the Ising research resource guard."""
import contextlib
import ctypes
import errno
import fcntl
import importlib.util
import os
import signal
from pathlib import Path
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('research_guard_test', ROOT / 'tools' / 'resource_guard.py')
guard = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = guard
spec.loader.exec_module(guard)


class GuardBoundaries(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.scratch = Path(tempfile.mkdtemp(prefix='ising3d-guard-regression-'))
        print('Retained small regression files:', cls.scratch, flush=True)

    def setUp(self):
        self.cwd = self.scratch / self._testMethodName
        self.cwd.mkdir()
        self.runtime = self.cwd / 'runtime'
        self.runtime_patch = patch.object(guard, 'RUNTIME', self.runtime)
        self.runtime_patch.start()

    def tearDown(self):
        self.runtime_patch.stop()

    def command(self, code):
        script = self.cwd / 'worker.py'
        script.write_text(code)
        return [sys.executable, '-I', '-B', str(script)]

    def test_native_group_enumeration_includes_owned_process(self):
        buffer = (ctypes.c_int * 1024)()
        samples, denied = guard.group_samples(os.getpgrp(), buffer)
        self.assertIn(os.getpid(), {item['pid'] for item in samples})
        self.assertEqual(denied, [])

    def test_setuid_descendant_does_not_abort_an_otherwise_clean_worker(self):
        # /bin/ps is setuid root, so proc_pid_rusage returns EPERM for it. The
        # Claude CLI spawns it at startup; sampling must survive that member.
        command = self.command(
            'import os, signal, subprocess, time\n'
            'children = []\n'
            'for _ in range(5):\n'
            '    child = subprocess.Popen(["/bin/ps", "-A"], stdout=subprocess.DEVNULL)\n'
            '    time.sleep(0.01)\n'
            '    os.kill(child.pid, signal.SIGSTOP)\n'
            '    children.append(child)\n'
            'time.sleep(0.3)\n'
            'for child in children:\n'
            '    os.kill(child.pid, signal.SIGCONT)\n'
            '    child.wait()\n'
        )
        result = guard.supervise(command, self.cwd)
        self.assertEqual((result['guard_exit_code'], result['child_exit_code']), (0, 0))
        self.assertEqual(result['reason'], 'completed')
        self.assertGreaterEqual(result['unmeasurable_processes'], 1,
                                'the unmeasurable-descendant path was never exercised')

    def test_unmeasurable_descendant_cannot_outlive_its_accounting_window(self):
        marker = self.cwd / 'descendant.pid'
        command = self.command(
            'import os, time\nfrom pathlib import Path\n'
            f'marker = Path({str(marker)!r})\n'
            'if os.fork() == 0:\n'
            '    marker.write_text(str(os.getpid()))\n'
            '    time.sleep(20)\n'
            'else:\n'
            '    while not marker.exists(): time.sleep(0.01)\n'
            '    time.sleep(20)\n'
        )
        measure = guard.sample_process

        def deny_descendant(pid):
            if marker.is_file() and pid == int(marker.read_text()):
                raise OSError(errno.EPERM, 'proc_pid_rusage', pid)
            return measure(pid)

        with patch.object(guard, 'sample_process', deny_descendant), \
                patch.object(guard, 'DENIED_GRACE_SECONDS', 0.2):
            result = guard.supervise(command, self.cwd)
        with contextlib.suppress(ProcessLookupError, ValueError):
            os.kill(int(marker.read_text()), signal.SIGKILL)
        self.assertEqual(result['guard_exit_code'], 1)
        # The supervisor classifies on this exact mark (GUARD_UNMEASURABLE_MARK).
        self.assertIn('unmeasurable descendant', result['reason'])
        self.assertIn('accounting window', result['reason'])
        self.assertGreaterEqual(result['unmeasurable_processes'], 1)

    def test_teardown_gives_up_instead_of_spinning_on_an_unkillable_member(self):
        # A member that never reports exit (privileged, or wedged in the kernel)
        # must not hold the guard inside its own finally block forever.
        buffer = (ctypes.c_int * 1024)()
        immortal = [{'pid': 2 ** 22, 'exit_ticks': 0, 'resident_bytes': 0,
                     'footprint_bytes': 0, 'cpu_seconds': 0.0, 'children_cpu_seconds': 0.0,
                     'written_bytes': 0, 'start_ticks': 1}]
        with patch.object(guard, 'group_samples', lambda pid, buf: (immortal, [4242])), \
                patch.object(guard, 'signal_owned_group', lambda pid, signum: None), \
                patch.object(guard, 'CLEANUP_SECONDS', 0.2):
            started = time.monotonic()
            quiet, denied = guard.stop_remaining_owned_group(os.getpgrp(), buffer)
        self.assertFalse(quiet, 'teardown claimed a quiet group while a member was alive')
        self.assertEqual(denied, [4242])
        self.assertLess(time.monotonic() - started, 5.0)

    def test_short_worker_failure_preserves_exit_status(self):
        result = guard.supervise(self.command('raise SystemExit(7)\n'), self.cwd)
        self.assertEqual((result['guard_exit_code'], result['child_exit_code']), (7, 7))
        self.assertEqual(result['reason'], 'completed')

    def test_resource_failure_overrides_an_inner_success(self):
        # Lower only this test's threshold; do not allocate pressure on the Mac.
        with patch.object(guard, 'RSS_CAP', 1):
            result = guard.supervise(self.command('print("INNER PASS", flush=True)\n'), self.cwd)
        self.assertEqual(result['guard_exit_code'], 1)
        self.assertIn('memory ceiling', result['reason'])

    def test_busy_worker_slot_refuses_before_launch(self):
        marker = self.cwd / 'must-not-be-written'
        command = self.command(f'from pathlib import Path\nPath({str(marker)!r}).write_text("bad")\n')
        self.runtime.mkdir()
        with (self.runtime / 'worker.lock').open('a') as holder:
            fcntl.flock(holder.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            with self.assertRaisesRegex(RuntimeError, 'single worker slot'):
                guard.supervise(command, self.cwd)
        self.assertFalse(marker.exists())

    def test_disk_reserve_refuses_before_launch(self):
        marker = self.cwd / 'must-not-be-written'
        command = self.command(f'from pathlib import Path\nPath({str(marker)!r}).write_text("bad")\n')
        with patch.object(guard, 'free_bytes', return_value=guard.FREE_FLOOR - 1):
            with self.assertRaisesRegex(RuntimeError, '50 GiB reserve'):
                guard.supervise(command, self.cwd)
        self.assertFalse(marker.exists())

    def test_cpu_admission_pacing_includes_final_off_time(self):
        command = self.command('import time\nend=time.process_time()+0.15\nwhile time.process_time()<end:\n    pass\n')
        result = guard.supervise(command, self.cwd)
        self.assertEqual(result['guard_exit_code'], 0)
        self.assertGreaterEqual(result['cpu_seconds'], 0.15)
        self.assertLessEqual(result['average_cpu_percent'], 35.05)

    def test_exited_leader_cannot_leave_a_running_descendant(self):
        marker = self.cwd / 'descendant.pid'
        command = self.command(
            'import os, signal, time\nfrom pathlib import Path\n'
            f'marker = Path({str(marker)!r})\n'
            'if os.fork() == 0:\n'
            '    signal.signal(signal.SIGTERM, signal.SIG_IGN)\n'
            '    marker.write_text(str(os.getpid()))\n'
            '    time.sleep(3)\n'
            'else:\n'
            '    while not marker.exists(): time.sleep(0.01)\n'
        )
        result = guard.supervise(command, self.cwd)
        descendant = int(marker.read_text())
        try:
            exited = bool(guard.sample_process(descendant)['exit_ticks'])
        except ProcessLookupError:
            exited = True
        self.assertTrue(exited, 'guard released its slot while its descendant was alive')
        self.assertEqual(result['guard_exit_code'], 1)
        self.assertIn('background descendants', result['reason'])

    def test_short_writer_cannot_bypass_total_write_ceiling(self):
        payload = self.cwd / 'payload'
        command = self.command(
            'import os\n'
            f'with open({str(payload)!r}, "wb") as stream:\n'
            '    stream.write(b"x" * 4096)\n'
            '    stream.flush()\n'
            '    os.fsync(stream.fileno())\n'
        )
        with patch.object(guard, 'WRITE_CAP', 1):
            result = guard.supervise(command, self.cwd)
        self.assertEqual(result['guard_exit_code'], 1)
        self.assertIn('write ceiling', result['reason'])
        self.assertGreater(result['written_bytes'], 0)


if __name__ == '__main__':
    unittest.main()
