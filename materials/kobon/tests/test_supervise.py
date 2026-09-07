"""Focused tests for the n=12 v2 campaign supervisor (math/kobon/n12_decide/supervise.py).

Covers: manifest mismatch refusal, campaign lock exclusion, resume of
exact-manifest terminal rows vs retry of nonterminal rows, unknown/duplicate
terminal row rejection, exact 15-cube all-UNSAT detection, incomplete
nonzero status logic, and process-group cleanup on INT/TERM via injectable
subprocesses.  Run with:

    python3 math/kobon/tests/test_supervise.py
"""

import hashlib
import io
import json
import signal
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from pathlib import Path

N12_DIR = Path(__file__).resolve().parents[1] / "n12_decide"
if str(N12_DIR) not in sys.path:
    sys.path.insert(0, str(N12_DIR))

import supervise  # noqa: E402


FIFTEEN = ["simple", "par_noc", "par_conc"] + [f"c{i}" for i in range(12)]


def manifest_text(cubes):
    doc = {
        "schema": supervise.MANIFEST_SCHEMA,
        "case_count": len(cubes),
        "case_order": list(cubes),
        "claim_scope": "test",
        "solver": {"path": "/usr/bin/true", "sha256": "0" * 64,
                   "version": "test"},
        "instances": [{"cube": c, "path": f"cnf/{c}.cnf"} for c in cubes],
    }
    return json.dumps(doc, indent=2, sort_keys=True) + "\n"


def digest_of(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


class FakeChild:
    """Popen double: settle() models process-group death."""

    def __init__(self, recorder, pgid, rc, settle_on_term):
        self._recorder = recorder
        self.pgid = pgid
        self._rc = rc
        self._settle_on_term = settle_on_term
        self._settled = threading.Event()
        self.returncode = None

    def settle(self):
        if not self._settled.is_set():
            self.returncode = self._rc
            self._settled.set()

    def wait(self, timeout=None):
        if self._settled.wait(timeout):
            return self.returncode
        return None

    def alive(self):
        return self.returncode is None


class SpawnRecorder:
    """Injectable spawn/killpg: records commands and process-group signals.

    By default each child settles (exits) immediately with its scripted rc;
    signal tests pass ``auto_settle=False`` so only killpg settles them.
    """

    def __init__(self, rcs, settle_on_term=True, auto_settle=True):
        self._rcs = list(rcs)
        self._settle_on_term = settle_on_term
        self._auto_settle = auto_settle
        self.spawn_cmds = []
        self.spawn_logs = []
        self.children = []
        self.killpg_calls = []
        self._base = 41000
        self._lock = threading.Lock()

    def spawn(self, cmd, log_fh):
        with self._lock:
            index = len(self.spawn_cmds)
            self.spawn_cmds.append(tuple(cmd))
            self.spawn_logs.append(getattr(log_fh, "name", "?"))
            rc = self._rcs[index]
            if rc == 10:
                log_fh.write(b"s SATISFIABLE\nv 1 0\n")
            elif rc == 20:
                log_fh.write(b"s UNSATISFIABLE\n")
            elif rc >= 0:
                log_fh.write(b"s UNKNOWN\n")
            log_fh.flush()
            child = FakeChild(self, self._base + index, rc,
                              self._settle_on_term)
            self.children.append(child)
            if self._auto_settle:
                child.settle()
            return child

    def killpg(self, pgid, sig):
        self.killpg_calls.append((pgid, sig))
        child = self.children[pgid - self._base]
        if sig == signal.SIGKILL or child._settle_on_term:
            child.settle()


class SupervisorTestCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.tmp = Path(self._tmp.name)

    # -- helpers ----------------------------------------------------------

    def build(self, cubes, recorder, jobs=1, pre_seeded=None, **kw):
        text = manifest_text(cubes)
        (self.tmp / "manifest.json").write_text(text, encoding="utf-8")
        if pre_seeded is not None:
            self.seed_ledger(text, pre_seeded)
        out, err = io.StringIO(), io.StringIO()
        sup = supervise.Supervisor(
            self.tmp, jobs=jobs, manifest_text=lambda: text,
            spawn=recorder.spawn, killpg=recorder.killpg,
            out=out, err=err, **kw)
        self.addCleanup(sup._restore_signal_handlers)
        return sup, out, err, text

    def ledger_row(self, digest, cube, attempt, verdict, exit_status=None,
                   log_path=None, log_sha=None):
        if exit_status is None:
            exit_status = {"SAT": "10", "UNSAT": "20"}.get(verdict, "1")
        if log_path is None:
            log_path = f"logs.v2.test/{cube}.attempt{attempt:03d}.log"
        if log_sha is None:
            log_sha = "0" * 64
        return "\t".join([
            cube, str(attempt), digest, "2026-09-03T00:00:00Z",
            "2026-09-03T00:00:01Z", exit_status, verdict, "1", "?",
            log_path, log_sha,
        ])

    def terminal_row(self, text, cube, attempt, verdict):
        digest = digest_of(text)
        _, logs = supervise.v2_namespace(self.tmp, digest)
        logs.mkdir(exist_ok=True)
        content = (
            b"s SATISFIABLE\nv 1 0\n"
            if verdict == "SAT" else b"s UNSATISFIABLE\n")
        path = logs / f"{cube}.attempt{attempt:03d}.log"
        path.write_bytes(content)
        return self.ledger_row(
            digest, cube, attempt, verdict,
            log_path=f"{logs.name}/{path.name}",
            log_sha=hashlib.sha256(content).hexdigest())

    def seed_ledger(self, text, rows):
        ledger, _ = supervise.v2_namespace(self.tmp, digest_of(text))
        lines = [supervise.LEDGER_HEADER] + list(rows)
        ledger.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return ledger

    def ledger_rows(self, sup):
        lines = sup.ledger_path.read_text(encoding="utf-8").splitlines()
        self.assertEqual(lines[0], supervise.LEDGER_HEADER)
        return [dict(zip(supervise.LEDGER_COLUMNS, line.split("\t")))
                for line in lines[1:] if line]

    def assert_prelaunch(self, sup, digest, resumed):
        self.assertTrue(sup.launches_path.is_file())
        records = [json.loads(line) for line
                   in sup.launches_path.read_text(encoding="utf-8")
                   .splitlines() if line]
        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertEqual(record["binding"], supervise.PRELAUNCH_BINDING)
        self.assertEqual(record["manifest_sha256"], digest)
        self.assertEqual(record["resumed_terminal"], resumed)
        self.assertEqual(record["v2_ledger"], sup.ledger_path.name)
        self.assertEqual(record["evidence_scope"], "DISCOVERY_ONLY")
        self.assertEqual(record["cnf_engine_provenance"], "UNVERIFIED")

    def assert_lock_released(self, sup):
        self.assertIsNotNone(sup._lock)
        self.assertIsNone(sup._lock._fd)
        probe = supervise.CampaignLock(sup.lock_path)
        probe.acquire()
        try:
            self.assertTrue(sup.lock_path.exists())
        finally:
            probe.release()

    # -- refusal paths ------------------------------------------------------

    def test_manifest_mismatch_refuses_before_launch(self):
        text = manifest_text(["c0", "c1"])
        (self.tmp / "manifest.json").write_text(text + "  ",
                                                encoding="utf-8")
        recorder = SpawnRecorder([20, 20])
        sup = supervise.Supervisor(
            self.tmp, manifest_text=lambda: text, spawn=recorder.spawn,
            killpg=recorder.killpg, out=io.StringIO(), err=io.StringIO())
        with self.assertRaises(supervise.Refused):
            sup.run()
        self.assertEqual(recorder.spawn_cmds, [])
        self.assertFalse(sup.launches_path.exists())
        self.assertFalse(
            (self.tmp / f"logs.v2.{digest_of(text)[:16]}").exists())
        # The campaign lock file is permanent: created, never deleted.
        self.assertTrue((self.tmp / supervise.LOCK_NAME).exists())
        self.assert_lock_released(sup)

    def test_campaign_lock_exclusion(self):
        holder_code = (
            "import fcntl, os, sys\n"
            "fd = os.open(sys.argv[1], os.O_RDWR | os.O_CREAT, 0o644)\n"
            "fcntl.flock(fd, fcntl.LOCK_EX)\n"
            "print('LOCKED', flush=True)\n"
            "sys.stdin.read()\n")
        holder = subprocess.Popen(
            [sys.executable, "-c", holder_code,
             str(self.tmp / supervise.LOCK_NAME)],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, text=True)

        def stop_holder():
            if holder.stdin is not None and not holder.stdin.closed:
                holder.stdin.close()
            try:
                holder.wait(timeout=5)
            except subprocess.TimeoutExpired:
                holder.kill()
                holder.wait()
            for stream in (holder.stdout, holder.stderr):
                if stream is not None:
                    stream.close()

        self.addCleanup(stop_holder)
        self.assertEqual(holder.stdout.readline().strip(), "LOCKED")
        recorder = SpawnRecorder([20, 20])
        sup, _, _, text = self.build(["c0", "c1"], recorder)
        with self.assertRaises(supervise.Refused):
            sup.run()
        self.assertEqual(recorder.spawn_cmds, [])
        self.assertFalse(sup.launches_path.exists())
        self.assertFalse(
            (self.tmp / f"ledger.v2.{digest_of(text)[:16]}.tsv").exists())
        self.assertIsNone(holder.poll())

    def test_lock_released_after_normal_completion(self):
        recorder = SpawnRecorder([20])
        sup, _, _, _ = self.build(["c0"], recorder)
        self.assertEqual(sup.run(), supervise.EXIT_OK)
        self.assert_lock_released(sup)

    # -- resume semantics ---------------------------------------------------

    def test_terminal_rows_resume_and_nonterminal_rows_retry(self):
        text = manifest_text(["c0", "c1"])
        digest = digest_of(text)
        seeded = [
            self.terminal_row(text, "c0", 1, "UNSAT"),
            self.ledger_row(digest, "c1", 1, "NO-VERDICT(rc=1)",
                            exit_status="1"),
        ]
        recorder = SpawnRecorder([20])
        sup, out, _, text = self.build(["c0", "c1"], recorder,
                                       pre_seeded=seeded)
        self.assertEqual(sup.run(), supervise.EXIT_OK)
        # Only the nonterminal lane re-ran; the terminal lane resumed.
        self.assertEqual(len(recorder.spawn_cmds), 1)
        rows = self.ledger_rows(sup)
        self.assertEqual(len(rows), 3)
        c0_rows = [r for r in rows if r["cube"] == "c0"]
        c1_rows = [r for r in rows if r["cube"] == "c1"]
        self.assertEqual([r["attempt"] for r in c0_rows], ["1"])
        self.assertEqual([r["attempt"] for r in c1_rows], ["1", "2"])
        self.assertEqual(c1_rows[-1]["verdict"], "UNSAT")
        self.assertEqual(c1_rows[-1]["exit_status"], "20")
        # Attempt log exists and its recorded hash matches the file.
        log_path = sup.logs_dir / "c1.attempt002.log"
        self.assertTrue(log_path.is_file())
        self.assertEqual(c1_rows[-1]["log_sha256"],
                         hashlib.sha256(log_path.read_bytes()).hexdigest())
        self.assertNotEqual(c1_rows[-1]["log_sha256"], "0" * 64)
        self.assert_prelaunch(sup, digest, resumed=1)
        self.assertIn("resumed=1 pending=1", out.getvalue())

    def test_orphan_attempt_log_is_never_rewritten(self):
        recorder = SpawnRecorder([20])
        sup, _, _, text = self.build(["c0"], recorder)
        _, logs = supervise.v2_namespace(self.tmp, digest_of(text))
        logs.mkdir()
        orphan = logs / "c0.attempt001.log"
        orphan.write_bytes(b"partial orphan output")
        self.assertEqual(sup.run(), supervise.EXIT_OK)
        self.assertEqual(orphan.read_bytes(), b"partial orphan output")
        rows = self.ledger_rows(sup)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["attempt"], "2")
        self.assertEqual(rows[0]["log_path"],
                         f"{logs.name}/c0.attempt002.log")


    def test_unknown_cube_terminal_row_is_rejected(self):
        text = manifest_text(["c0", "c1"])
        seeded = [self.ledger_row(digest_of(text), "zzz", 1, "UNSAT")]
        recorder = SpawnRecorder([20, 20])
        sup, _, _, _ = self.build(["c0", "c1"], recorder,
                                  pre_seeded=seeded)
        with self.assertRaises(supervise.Refused):
            sup.run()
        self.assertEqual(recorder.spawn_cmds, [])
        self.assertFalse(sup.launches_path.exists())

    def test_duplicate_terminal_rows_are_rejected(self):
        text = manifest_text(["c0", "c1"])
        digest = digest_of(text)
        seeded = [
            self.ledger_row(digest, "c0", 1, "UNSAT"),
            self.ledger_row(digest, "c0", 2, "UNSAT"),
        ]
        recorder = SpawnRecorder([20])
        sup, _, _, _ = self.build(["c0", "c1"], recorder,
                                  pre_seeded=seeded)
        with self.assertRaises(supervise.Refused):
            sup.run()
        self.assertEqual(recorder.spawn_cmds, [])
        self.assertFalse(sup.launches_path.exists())

    def test_tampered_terminal_log_is_rejected(self):
        text = manifest_text(["c0", "c1"])
        row = self.terminal_row(text, "c0", 1, "UNSAT")
        _, logs = supervise.v2_namespace(self.tmp, digest_of(text))
        (logs / "c0.attempt001.log").write_bytes(b"s SATISFIABLE\nv 1 0\n")
        recorder = SpawnRecorder([20])
        sup, _, _, _ = self.build(["c0", "c1"], recorder,
                                  pre_seeded=[row])
        with self.assertRaises(supervise.Refused):
            sup.run()
        self.assertEqual(recorder.spawn_cmds, [])
        self.assertFalse(sup.launches_path.exists())

    def test_foreign_manifest_digest_and_bad_header_are_rejected(self):
        text = manifest_text(["c0", "c1"])
        foreign = [
            self.ledger_row("b" * 64, "c0", 1, "UNSAT"),
        ]
        recorder = SpawnRecorder([20])
        sup, _, _, _ = self.build(["c0", "c1"], recorder,
                                  pre_seeded=foreign)
        with self.assertRaises(supervise.Refused):
            sup.run()
        self.assertEqual(recorder.spawn_cmds, [])
        self.assert_lock_released(sup)

        bad_header = self.tmp / "manifest.json"
        self.assertEqual(bad_header.read_text(encoding="utf-8"),
                         manifest_text(["c0", "c1"]))
        recorder2 = SpawnRecorder([20])
        ledger, _ = supervise.v2_namespace(self.tmp, digest_of(text))
        ledger.write_text("cube\tverdict\n", encoding="utf-8")
        sup2, _, _, _ = self.build(["c0", "c1"], recorder2)
        with self.assertRaises(supervise.Refused):
            sup2.run()
        self.assertEqual(recorder2.spawn_cmds, [])
        self.assert_lock_released(sup2)

    # -- outcome semantics --------------------------------------------------
    def test_exact_fifteen_cube_all_unsat_exits_zero(self):
        recorder = SpawnRecorder([20] * 15)
        sup, out, _, text = self.build(FIFTEEN, recorder, jobs=4)
        self.assertEqual(sup.run(), supervise.EXIT_OK)
        digest = digest_of(text)
        rows = self.ledger_rows(sup)
        self.assertEqual(len(rows), 15)
        self.assertEqual(sorted(r["cube"] for r in rows), sorted(FIFTEEN))
        for row in rows:
            self.assertEqual(row["verdict"], "UNSAT")
            self.assertEqual(row["attempt"], "1")
            self.assertEqual(row["manifest_sha256"], digest)
        transcript = out.getvalue()
        self.assertIn("solver-UNSAT", transcript)
        self.assertIn("exactly one terminal row per named cube", transcript)
        self.assertIn("discovery evidence ONLY", transcript)
        self.assertIn("does not prove", transcript)
        self.assertIn("byte-regeneration/provenance/cover audit", transcript)
        self.assertIn("DRAT", transcript)
        self.assertIn("CUBES_DONE", transcript)

    def test_incomplete_coverage_exits_nonzero(self):
        rcs = [20] * 15
        rcs[7] = 1  # one lane dies with NO-VERDICT(rc=1)
        recorder = SpawnRecorder(rcs)
        sup, out, _, _ = self.build(FIFTEEN, recorder, jobs=3)
        self.assertEqual(sup.run(), supervise.EXIT_INCOMPLETE)
        rows = self.ledger_rows(sup)
        terminal = [r for r in rows if r["verdict"] in
                    supervise.TERMINAL_VERDICTS]
        self.assertEqual(len(terminal), 14)
        self.assertIn("INCOMPLETE", out.getvalue())
        self.assertIn("No result is claimed", out.getvalue())

    def test_sat_is_terminal_discovery_and_remaining_lanes_are_moot(self):
        recorder = SpawnRecorder([10])
        sup, out, _, _ = self.build(FIFTEEN, recorder, jobs=1)
        self.assertEqual(sup.run(), supervise.EXIT_OK)
        self.assertEqual(len(recorder.spawn_cmds), 1)
        self.assertIn("Boolean candidate ONLY", out.getvalue())
        self.assertIn("moot", out.getvalue())
        rows = self.ledger_rows(sup)
        self.assertEqual([r["verdict"] for r in rows], ["SAT"])

    # -- process-group cleanup ----------------------------------------------

    def run_in_thread(self, sup):
        result = {}
        thread = threading.Thread(
            target=lambda: result.update(rc=sup.run()), daemon=True)
        thread.start()
        deadline = time.monotonic() + 5.0
        while not sup.children and time.monotonic() < deadline:
            time.sleep(0.005)
        self.assertTrue(sup.children, "supervisor never spawned a child")
        return thread, result

    def test_signal_terminates_child_process_group(self):
        recorder = SpawnRecorder([-signal.SIGTERM, 20],
                                 settle_on_term=True, auto_settle=False)
        sup, out, _, _ = self.build(["c0", "c1"], recorder)
        thread, result = self.run_in_thread(sup)
        pgid = recorder.children[0].pgid
        sup._on_signal(signal.SIGTERM, None)
        thread.join(5.0)
        self.assertFalse(thread.is_alive())
        self.assertEqual(result["rc"],
                         supervise.SIGNAL_EXIT_BASE + signal.SIGTERM)
        self.assertEqual(recorder.killpg_calls, [(pgid, signal.SIGTERM)])
        self.assertFalse(recorder.children[0].alive())
        rows = self.ledger_rows(sup)
        self.assertEqual([r["verdict"] for r in rows], ["INTERRUPTED"])
        self.assertEqual(len(recorder.spawn_cmds), 1)  # no further lanes
        self.assertIn("INTERRUPTED", out.getvalue())
        self.assertIn("CUBES_DONE", out.getvalue())

    def test_signal_escalates_to_sigkill_after_grace(self):
        recorder = SpawnRecorder([-9, 20], settle_on_term=False,
                                 auto_settle=False)
        sup, _, _, _ = self.build(["c0", "c1"], recorder, kill_grace=0.05)
        thread, result = self.run_in_thread(sup)
        pgid = recorder.children[0].pgid
        sup._on_signal(signal.SIGTERM, None)
        thread.join(5.0)
        self.assertFalse(thread.is_alive())
        self.assertEqual(result["rc"],
                         supervise.SIGNAL_EXIT_BASE + signal.SIGTERM)
        self.assertEqual(recorder.killpg_calls,
                         [(pgid, signal.SIGTERM), (pgid, signal.SIGKILL)])
        self.assertFalse(recorder.children[0].alive())
        rows = self.ledger_rows(sup)
        self.assertEqual([r["verdict"] for r in rows], ["KILLED(sig=9)"])


if __name__ == "__main__":
    unittest.main()
