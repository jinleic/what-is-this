#!/usr/bin/env python3
"""Supervised, lock-guarded runner for the n=12 target-39 cube campaign (v2).

Safety contract (run_cubes.sh is only an exec wrapper; this file is the
campaign implementation and is hash-pinned by manifest.json):

* Holds the permanent nonblocking fcntl campaign lock
  ``.campaign.lock`` for its full lifetime; a second supervisor refuses
  to launch while one is held (the lock file is never deleted).
* Validates the already-generated ``manifest.json`` byte-for-byte against a
  fresh regeneration over the current 15 CNFs, engine, exec wrapper,
  supervisor implementation, and the exact resolved Kissat binary.  Any
  drift refuses the launch before anything is spawned.  This binds those
  existing bytes together but does NOT prove that the current engine emitted
  the existing CNFs; byte-regeneration/provenance audit remains mandatory
  before certification.
* Appends a prelaunch manifest digest record to ``launches.jsonl``
  (append-only; the historical RETROSPECTIVE_UNBOUND record is never
  rewritten).
* Keeps all new results in a manifest-keyed v2 namespace:
  ``ledger.v2.<manifest[:16]>.tsv`` and ``logs.v2.<manifest[:16]>/``.
  The legacy five-column ``verdicts.tsv`` and ``logs/`` remain untouched
  and are interpreted as discovery-only.
* Puts each Kissat in its own process group (``start_new_session``) and on
  INT/TERM/exit terminates and reaps every child (TERM, then KILL after a
  grace period).
* Bounds concurrency (JOBS, default 7) and the per-lane wall budget
  (TIMEOUT seconds, 0 = unlimited, passed to kissat ``--time``).
* Solver stdout is preserved verbatim per attempt in the v2 log (including
  the ``s`` line and ``v`` model lines for SAT), and every attempt appends
  manifest digest, cube, exit status, verdict, wall time, peak RSS, log
  path, and log hash to the v2 ledger.
* Resume honors only exact-manifest SAT/UNSAT rows; nonterminal rows retry
  automatically; unknown or duplicate terminal rows refuse the run.
  All-UNSAT is declared only when every named cube has exactly one
  terminal row and all are UNSAT.
* Evidence boundary: a solver SAT is a Boolean candidate only; a proofless
  all-UNSAT is discovery evidence only.  DRAT emission/verification is
  deliberately out of scope here.
"""

from __future__ import annotations

import argparse
import atexit
import fcntl
import hashlib
import json
import os
import queue
import re
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

if str(Path(__file__).resolve().parent) not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parent))

import make_manifest

HERE = Path(__file__).resolve().parent
TIME_BIN = Path("/usr/bin/time")

DEFAULT_JOBS = 7
EXIT_OK = 0
EXIT_INCOMPLETE = 1
EXIT_REFUSED = 2
SIGNAL_EXIT_BASE = 128
KEY_LEN = 16
KILL_GRACE_SECONDS = 10.0
WAIT_POLL_SECONDS = 0.2
LOCK_NAME = ".campaign.lock"
MANIFEST_SCHEMA = "kobon-n12-decision-manifest/2"
PRELAUNCH_BINDING = "MANIFEST_PREVALIDATED_V2"
TERMINAL_VERDICTS = frozenset({"SAT", "UNSAT"})
VERDICT_EXIT_CODES = {10: "SAT", 20: "UNSAT"}
LEDGER_COLUMNS = (
    "cube", "attempt", "manifest_sha256", "started_utc", "finished_utc",
    "exit_status", "verdict", "secs", "peak_rss_gb", "log_path",
    "log_sha256",
)
LEDGER_HEADER = "\t".join(LEDGER_COLUMNS)
_RSS_RE = re.compile(rb"([0-9]+)[ \t]+maximum resident set size")
_STATUS_RE = re.compile(
    rb"^s (SATISFIABLE|UNSATISFIABLE)\r?$", re.MULTILINE)
_MODEL_RE = re.compile(rb"^v(?:[ \t]|$)", re.MULTILINE)
_APPEND_LOCK = threading.Lock()


class Refused(RuntimeError):
    """The launch is refused before any solver is spawned."""


def utc_now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def append_line(path: Path, line: str) -> None:
    """Append one fsynced record through O_APPEND without interleaving."""
    payload = (line + "\n").encode("utf-8")
    with _APPEND_LOCK:
        fd = os.open(
            path,
            os.O_WRONLY | os.O_CREAT | os.O_APPEND | os.O_CLOEXEC,
            0o644,
        )
        try:
            view = memoryview(payload)
            while view:
                written = os.write(fd, view)
                if written <= 0:
                    raise OSError("short append write")
                view = view[written:]
            os.fsync(fd)
        finally:
            os.close(fd)


def _inspect_log(log_path: Path) -> dict:
    try:
        data = Path(log_path).read_bytes()
    except OSError:
        return {
            "sha256": "-",
            "peak_rss_gb": "?",
            "status": None,
            "has_model": False,
        }
    rss_matches = _RSS_RE.findall(data)
    rss = ("?" if not rss_matches else
           f"{int(rss_matches[-1]) / (1 << 30):.3f}")
    statuses = _STATUS_RE.findall(data)
    if len(statuses) == 1:
        status = {
            b"SATISFIABLE": "SAT",
            b"UNSATISFIABLE": "UNSAT",
        }[statuses[0]]
    elif statuses:
        status = "AMBIGUOUS"
    else:
        status = None
    return {
        "sha256": hashlib.sha256(data).hexdigest(),
        "peak_rss_gb": rss,
        "status": status,
        "has_model": bool(_MODEL_RE.search(data)),
    }


def _default_manifest_text() -> str:
    """Deterministically regenerate the manifest text (make_manifest.py)."""
    return json.dumps(make_manifest.build(), indent=2, sort_keys=True) + "\n"


def v2_namespace(here: Path, digest: str) -> tuple[Path, Path]:
    """Manifest-keyed v2 ledger path and log directory."""
    key = digest[:KEY_LEN]
    here = Path(here)
    return here / f"ledger.v2.{key}.tsv", here / f"logs.v2.{key}"


class CampaignLock:
    """Permanent nonblocking fcntl campaign lock; the file is never deleted."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self._fd: int | None = None

    def acquire(self) -> None:
        fd = os.open(self.path, os.O_RDWR | os.O_CREAT | os.O_CLOEXEC, 0o644)
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            os.close(fd)
            raise Refused(
                f"campaign lock {self.path} is held by another process: "
                f"{exc}") from exc
        self._fd = fd  # held for the supervisor's full lifetime

    def release(self) -> None:
        if self._fd is not None:
            try:
                fcntl.flock(self._fd, fcntl.LOCK_UN)
            finally:
                os.close(self._fd)
                self._fd = None


class Child:
    """Minimal uniform view over a real Popen for the supervisor."""

    def __init__(self, popen: subprocess.Popen):
        self._popen = popen

    @property
    def pgid(self) -> int:
        return self._popen.pid  # start_new_session: pid == pgid

    def wait(self, timeout: float | None = None):
        try:
            return self._popen.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            return None

    def alive(self) -> bool:
        return self._popen.poll() is None


def default_spawn(cmd: list[str], log_fh) -> Child:
    """Spawn one solver lane in its own session/process group."""
    popen = subprocess.Popen(
        cmd,
        stdout=log_fh,
        stderr=subprocess.STDOUT,
        start_new_session=True,
        close_fds=True,
    )
    return Child(popen)


class Supervisor:
    """Runs the manifest-pinned cube campaign to a resumable v2 ledger."""

    def __init__(
        self,
        here: Path,
        *,
        jobs: int = DEFAULT_JOBS,
        timeout: int = 0,
        manifest_text=None,
        spawn=None,
        killpg=None,
        kill_grace: float = KILL_GRACE_SECONDS,
        now=utc_now,
        out=None,
        err=None,
    ):
        self.here = Path(here).resolve()
        self.lock_path = self.here / LOCK_NAME
        self.manifest_path = self.here / "manifest.json"
        self.launches_path = self.here / "launches.jsonl"
        self.jobs = int(jobs)
        if self.jobs < 1:
            raise Refused("JOBS must be at least 1")
        self.timeout = int(timeout)
        if self.timeout < 0:
            raise Refused("TIMEOUT must be nonnegative")
        self._manifest_text_fn = manifest_text or _default_manifest_text
        self._spawn_fn = spawn or default_spawn
        self._killpg = killpg or os.killpg
        self.kill_grace = float(kill_grace)
        if self.kill_grace < 0:
            raise Refused("kill grace must be nonnegative")
        self._now = now
        self._out = out if out is not None else sys.stdout
        self._err = err if err is not None else sys.stderr

        self.digest: str | None = None
        self.case_order: list[str] = []
        self.cnf_by_cube: dict[str, Path] = {}
        self.solver_path: str | None = None
        self.solver_sha256 = "?"
        self.solver_version = "?"
        self.ledger_path: Path | None = None
        self.logs_dir: Path | None = None
        self.terminal: dict[str, dict] = {}
        self.attempts: dict[str, int] = {}
        self.row_count = 0

        self.stop = threading.Event()
        self.stop_signum: int | None = None
        self.discovered_sat = False
        self._child_lock = threading.Lock()
        self._rows_lock = threading.Lock()
        self.children: set = set()
        self.signalled: set[int] = set()
        self._old_handlers: dict = {}
        self._lock: CampaignLock | None = None

    # -- plumbing ---------------------------------------------------------

    def _say(self, message: str) -> None:
        print(message, file=self._out)

    def _workspace(self) -> Path:
        parents = self.here.parents
        return parents[2] if len(parents) > 2 else self.here

    def _say_banner(self, resumed: int, pending: int) -> None:
        self._say(
            "=== kobon n=12 t=39 cube cover: supervised deciding (v2) ===")
        self._say(
            f"manifest pinned {self.digest} (CNFs x{len(self.case_order)}, "
            f"engine, wrapper, supervisor, kissat {self.solver_version} "
            f"@ {self.solver_path})")
        self._say(
            "binding is byte identity only: it does not prove that the "
            "current engine emitted the existing CNFs")
        self._say(
            f"ledger {self.ledger_path.name}  logs {self.logs_dir.name}/")
        self._say(
            f"jobs={self.jobs} timeout={self.timeout} "
            f"resumed={resumed} pending={pending}")

    # -- launch validation -------------------------------------------------

    def _pin_manifest(self) -> bytes:
        try:
            expected = self._manifest_text_fn().encode("utf-8")
        except (make_manifest.ManifestFailure, OSError,
                subprocess.SubprocessError) as exc:
            raise Refused(f"manifest regeneration failed: {exc}") from exc
        try:
            actual = self.manifest_path.read_bytes()
        except OSError as exc:
            raise Refused(
                f"cannot read manifest {self.manifest_path}: {exc}") from exc
        if actual != expected:
            raise Refused(
                "manifest.json does not match current inputs byte-for-byte "
                "(a CNF, engine, wrapper, supervisor, or the resolved kissat "
                "binary drifted); regenerate with make_manifest.py before "
                "launching")
        self.digest = hashlib.sha256(expected).hexdigest()
        return actual

    def _parse_manifest(self, raw: bytes) -> None:
        try:
            doc = json.loads(raw.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise Refused(f"manifest is not valid UTF-8 JSON: {exc}") from exc
        if doc.get("schema") != MANIFEST_SCHEMA:
            raise Refused(
                f"manifest schema {doc.get('schema')!r} != "
                f"{MANIFEST_SCHEMA!r}")
        case_order = doc.get("case_order")
        if (not isinstance(case_order, list) or not case_order
                or not all(isinstance(cube, str) and cube
                           for cube in case_order)
                or len(set(case_order)) != len(case_order)):
            raise Refused(
                "manifest case_order must be a nonempty duplicate-free list")
        if doc.get("case_count") != len(case_order):
            raise Refused("manifest case_count disagrees with case_order")
        solver = doc.get("solver")
        if not isinstance(solver, dict) or not solver.get("path"):
            raise Refused("manifest lacks a solver record with a path")
        instances = doc.get("instances")
        if not isinstance(instances, list) or len(instances) != len(case_order):
            raise Refused("manifest instances disagree with case_order")
        workspace = self._workspace()
        for instance in instances:
            if not isinstance(instance, dict):
                raise Refused("manifest instance is not an object")
            cube = instance.get("cube")
            path = instance.get("path")
            if cube not in case_order:
                raise Refused(
                    f"manifest instance cube {cube!r} is not in case_order")
            if cube in self.cnf_by_cube or not path:
                raise Refused(f"manifest instance for {cube!r} is malformed")
            cnf = Path(path)
            self.cnf_by_cube[cube] = (
                cnf if cnf.is_absolute() else workspace / cnf)
        self.case_order = [str(cube) for cube in case_order]
        self.solver_path = str(solver["path"])
        self.solver_sha256 = str(solver.get("sha256", "?"))
        self.solver_version = str(solver.get("version", "?"))

    # -- v2 ledger ---------------------------------------------------------

    def _read_ledger(self):
        rows: list[dict] = []
        if not self.ledger_path.exists():
            return rows, {}, {}
        raw = self.ledger_path.read_text(encoding="utf-8")
        if not raw.strip():
            return rows, {}, {}
        lines = raw.splitlines()
        if lines[0] != LEDGER_HEADER:
            raise Refused(
                f"v2 ledger {self.ledger_path.name} has an unexpected "
                "header; refusing to interpret it")
        named = set(self.case_order)
        terminal: dict[str, dict] = {}
        terminal_lines: dict[str, int] = {}
        attempts: dict[str, int] = {}
        for lineno, line in enumerate(lines[1:], start=2):
            if not line.strip():
                raise Refused(f"v2 ledger line {lineno}: blank row")
            fields = line.split("\t")
            if len(fields) != len(LEDGER_COLUMNS):
                raise Refused(
                    f"v2 ledger line {lineno}: expected "
                    f"{len(LEDGER_COLUMNS)} tab-separated fields, found "
                    f"{len(fields)}")
            row = dict(zip(LEDGER_COLUMNS, fields))
            cube = row["cube"]
            if cube not in named:
                raise Refused(
                    f"v2 ledger line {lineno}: unknown cube {cube!r} is not "
                    "in the pinned manifest case list")
            if row["manifest_sha256"] != self.digest:
                raise Refused(
                    f"v2 ledger line {lineno}: manifest digest "
                    f"{row['manifest_sha256']!r} does not match the pinned "
                    f"{self.digest}; a manifest-keyed ledger holds one "
                    "manifest's rows only")
            try:
                attempt = int(row["attempt"])
            except ValueError:
                raise Refused(
                    f"v2 ledger line {lineno}: non-integer attempt") from None
            if attempt < 1:
                raise Refused(
                    f"v2 ledger line {lineno}: attempt must be >= 1")
            if row["verdict"] in TERMINAL_VERDICTS:
                if cube in terminal:
                    raise Refused(
                        f"v2 ledger line {lineno}: duplicate terminal row "
                        f"for cube {cube!r} (previous at line "
                        f"{terminal_lines[cube]}); refusing ambiguous resume")
                terminal[cube] = row
                terminal_lines[cube] = lineno
            attempts[cube] = max(attempts.get(cube, 0), attempt)
            rows.append(row)
        for cube, row in terminal.items():
            self._validate_terminal_row(row, terminal_lines[cube])
        return rows, terminal, attempts

    def _validate_terminal_row(self, row: dict, lineno: int) -> None:
        cube = row["cube"]
        verdict = row["verdict"]
        expected_exit = {"SAT": "10", "UNSAT": "20"}[verdict]
        if row["exit_status"] != expected_exit:
            raise Refused(
                f"v2 ledger line {lineno}: terminal {verdict} for cube "
                f"{cube!r} has exit {row['exit_status']!r}, expected "
                f"{expected_exit}")
        attempt = int(row["attempt"])
        expected_rel = (
            f"{self.logs_dir.name}/{cube}.attempt{attempt:03d}.log")
        if row["log_path"] != expected_rel:
            raise Refused(
                f"v2 ledger line {lineno}: log path {row['log_path']!r} "
                f"does not exactly match {expected_rel!r}")
        info = _inspect_log(self.here / expected_rel)
        if info["sha256"] == "-":
            raise Refused(
                f"v2 ledger line {lineno}: terminal log is missing: "
                f"{expected_rel}")
        if row["log_sha256"] != info["sha256"]:
            raise Refused(
                f"v2 ledger line {lineno}: terminal log hash mismatch for "
                f"cube {cube!r}")
        if info["status"] != verdict:
            raise Refused(
                f"v2 ledger line {lineno}: terminal log status "
                f"{info['status']!r} does not match {verdict}")
        if verdict == "SAT" and not info["has_model"]:
            raise Refused(
                f"v2 ledger line {lineno}: SAT log has no preserved model")

    def _ensure_ledger_file(self) -> None:
        if self.ledger_path.exists() and self.ledger_path.stat().st_size > 0:
            return
        append_line(self.ledger_path, LEDGER_HEADER)

    def _append_prelaunch(self, resumed: int) -> None:
        record = {
            "binding": PRELAUNCH_BINDING,
            "cnf_engine_provenance": "UNVERIFIED",
            "evidence_scope": "DISCOVERY_ONLY",
            "jobs": self.jobs,
            "manifest_sha256": self.digest,
            "process_id": os.getpid(),
            "resumed_terminal": resumed,
            "solver_path": self.solver_path,
            "solver_sha256": self.solver_sha256,
            "started_utc": self._now(),
            "timeout_seconds": self.timeout,
            "v2_ledger": self.ledger_path.name,
            "v2_logs_dir": self.logs_dir.name,
        }
        append_line(
            self.launches_path,
            json.dumps(record, sort_keys=True, separators=(",", ":")))

    def _record_row(self, row: dict) -> None:
        append_line(
            self.ledger_path,
            "\t".join(str(row[column]) for column in LEDGER_COLUMNS))
        with self._rows_lock:
            self.row_count += 1
            if row["verdict"] in TERMINAL_VERDICTS:
                self.terminal[row["cube"]] = row
                if row["verdict"] == "SAT":
                    self.discovered_sat = True
        self._say(
            f"[done] {row['cube']} -> {row['verdict']} (attempt "
            f"{row['attempt']}, exit {row['exit_status']}, {row['secs']}s, "
            f"peak {row['peak_rss_gb']} GB, log {row['log_path']})")

    # -- children ----------------------------------------------------------

    def _register(self, child) -> None:
        with self._child_lock:
            self.children.add(child)

    def _unregister(self, child) -> None:
        with self._child_lock:
            self.children.discard(child)
            self.signalled.discard(child.pgid)

    def _signal(self, child, sig: int) -> None:
        try:
            self._killpg(child.pgid, sig)
        except ProcessLookupError:
            pass
        except OSError as exc:
            print(
                f"[error] cannot signal process group {child.pgid}: {exc}",
                file=self._err)
            return
        with self._child_lock:
            self.signalled.add(child.pgid)

    def _request_terminate(self, child) -> None:
        with self._child_lock:
            already = child.pgid in self.signalled
        if not already:
            self._signal(child, signal.SIGTERM)

    def _on_signal(self, signum, _frame) -> None:
        self.stop_signum = signum
        self.stop.set()
        with self._child_lock:
            children = list(self.children)
        for child in children:
            self._request_terminate(child)

    def _wait_child(self, child):
        grace_deadline = None
        while True:
            rc = child.wait(timeout=WAIT_POLL_SECONDS)
            if rc is not None:
                return rc
            if self.stop.is_set():
                now = time.monotonic()
                if grace_deadline is None:
                    self._request_terminate(child)
                    grace_deadline = now + self.kill_grace
                elif now >= grace_deadline:
                    self._signal(child, signal.SIGKILL)
                    grace_deadline = now + self.kill_grace

    def _terminate_and_wait(self, children) -> None:
        running = [child for child in children if child.alive()]
        for child in running:
            self._request_terminate(child)
        deadline = time.monotonic() + self.kill_grace
        for child in running:
            while child.alive():
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    break
                child.wait(timeout=min(WAIT_POLL_SECONDS, remaining))
        for child in running:
            if child.alive():
                self._signal(child, signal.SIGKILL)
        for child in running:
            while child.alive():
                child.wait(timeout=WAIT_POLL_SECONDS)

    def _ensure_dead(self, child) -> None:
        self._terminate_and_wait([child])

    def _sweep(self) -> None:
        with self._child_lock:
            children = list(self.children)
        self._terminate_and_wait(children)
        for child in children:
            self._unregister(child)

    def _atexit_cleanup(self) -> None:
        self.stop.set()
        self._sweep()

    def _install_signal_handlers(self) -> None:
        for sig in (signal.SIGINT, signal.SIGTERM):
            try:
                self._old_handlers[sig] = signal.getsignal(sig)
                signal.signal(sig, self._on_signal)
            except (ValueError, OSError):
                self._old_handlers.pop(sig, None)

    def _restore_signal_handlers(self) -> None:
        for sig, handler in self._old_handlers.items():
            try:
                signal.signal(sig, handler)
            except (ValueError, OSError):
                pass
        self._old_handlers = {}

    # -- lanes ---------------------------------------------------------------

    def _solver_cmd(self, cnf: Path) -> list[str]:
        # --quiet keeps logs small yet still prints the s line and, for SAT,
        # the v model lines; stdout is preserved verbatim in the lane log.
        core = [self.solver_path, "--quiet"]
        if self.timeout > 0:
            core.append(f"--time={self.timeout}")
        core.append(str(cnf))
        if TIME_BIN.exists():
            # /usr/bin/time -l appends the max-RSS line to the merged log.
            return [str(TIME_BIN), "-l"] + core
        return core

    def _verdict_for(self, exit_status, log_info: dict) -> str:
        verdict = VERDICT_EXIT_CODES.get(exit_status)
        if verdict is not None:
            if log_info["status"] != verdict:
                return f"INVALID-OUTPUT(rc={exit_status})"
            if verdict == "SAT" and not log_info["has_model"]:
                return "INVALID-SAT-MODEL"
            return verdict
        if isinstance(exit_status, int) and exit_status < 0:
            sig = -exit_status
            if self.stop.is_set() and sig == signal.SIGTERM:
                return "INTERRUPTED"
            return f"KILLED(sig={sig})"
        return f"NO-VERDICT(rc={exit_status})"

    def _open_attempt_log(self, cube: str):
        """Reserve a never-before-used attempt log without truncating files."""
        with self._rows_lock:
            attempt = self.attempts.get(cube, 0) + 1
            while True:
                log_path = self.logs_dir / (
                    f"{cube}.attempt{attempt:03d}.log")
                try:
                    log_fh = open(log_path, "xb")
                except FileExistsError:
                    attempt += 1
                    continue
                self.attempts[cube] = attempt
                return attempt, log_path, log_fh

    def _attempt(self, cube: str) -> None:
        attempt, log_path, log_fh = self._open_attempt_log(cube)
        started = self._now()
        wall0 = time.monotonic()
        child = None
        completed = False
        exit_status = "-"
        verdict = "SPAWN-ERROR"
        try:
            with log_fh:
                child = self._spawn_fn(self._solver_cmd(self.cnf_by_cube[cube]),
                                       log_fh)
                self._register(child)
                exit_status = self._wait_child(child)
                completed = True
        except Exception as exc:  # attempts are recorded, never fatal mid-run
            if child is None:
                verdict, exit_status = "SPAWN-ERROR", "-"
            else:
                verdict, exit_status = "IO-ERROR", "-"
            print(f"[error] {cube}: {verdict}: {exc!r}", file=self._err)
        finally:
            if child is not None:
                self._ensure_dead(child)
                self._unregister(child)
        log_info = _inspect_log(log_path)
        if completed:
            verdict = self._verdict_for(exit_status, log_info)
        secs = int(round(time.monotonic() - wall0))
        try:
            log_rel = str(log_path.relative_to(self.here))
        except ValueError:
            log_rel = str(log_path)
        self._record_row({
            "cube": cube,
            "attempt": attempt,
            "manifest_sha256": self.digest,
            "started_utc": started,
            "finished_utc": self._now(),
            "exit_status": exit_status,
            "verdict": verdict,
            "secs": secs,
            "peak_rss_gb": log_info["peak_rss_gb"],
            "log_path": log_rel,
            "log_sha256": log_info["sha256"],
        })

    def _worker(self, work: "queue.Queue[str]") -> None:
        while True:
            try:
                cube = work.get_nowait()
            except queue.Empty:
                return
            try:
                if self.stop.is_set() or self.discovered_sat:
                    continue
                self._attempt(cube)
            except Exception as exc:  # keep lanes alive; row may be missing
                print(f"[error] lane failure on {cube}: {exc!r}",
                      file=self._err)
            finally:
                work.task_done()

    def _run_workers(self) -> None:
        pending = [c for c in self.case_order if c not in self.terminal]
        if (not pending or self.stop.is_set()
                or self.discovered_sat):
            return
        work: "queue.Queue[str]" = queue.Queue()
        for cube in pending:
            work.put(cube)
        lanes = min(self.jobs, len(pending))
        workers = [
            threading.Thread(target=self._worker, args=(work,),
                             name=f"kobon-lane-{index}", daemon=True)
            for index in range(lanes)
        ]
        atexit.register(self._atexit_cleanup)
        try:
            for worker in workers:
                worker.start()
            for worker in workers:
                worker.join()
        finally:
            if any(worker.is_alive() for worker in workers):
                self.stop.set()
            self._sweep()
            for worker in workers:
                if worker.is_alive():
                    worker.join()
            atexit.unregister(self._atexit_cleanup)

    # -- outcome -------------------------------------------------------------

    def _decide(self) -> int:
        if self.stop.is_set():
            self._say(
                f"INTERRUPTED by signal {self.stop_signum}: children were "
                "terminated and reaped; this session's attempts are recorded "
                "as nonterminal rows.  Rerun to resume (exact-manifest "
                "SAT/UNSAT rows are kept).")
            return SIGNAL_EXIT_BASE + (self.stop_signum or 0)
        named = self.case_order
        terminal = self.terminal
        sat = [c for c in named
               if terminal.get(c, {}).get("verdict") == "SAT"]
        unsat = [c for c in named
                 if terminal.get(c, {}).get("verdict") == "UNSAT"]
        total = len(named)
        if sat:
            skipped = total - len(terminal)
            moot = (
                f"  {skipped} lane(s) have no terminal row; further "
                "coverage is moot for this discovery outcome."
                if skipped > 0 else "")
            self._say(
                "DISCOVERY: solver-SAT on cube(s) " + ", ".join(sat) + "."
                + moot)
            self._say(
                "SAT is a Boolean candidate ONLY: extract the model, "
                "straighten it, and pass exact-rational "
                "verify_selection(minimum=39).  SAT alone proves no "
                "39-triangle configuration.")
            return EXIT_OK
        if len(unsat) == total and len(terminal) == total:
            self._say(
                f"DISCOVERY: all {total} relaxed target-39 cubes are "
                "solver-UNSAT with exactly one terminal row per named cube, "
                "supporting K_gen(12) <= 38.")
            self._say(
                "This is discovery evidence ONLY: the run is proof-free.  "
                "The manifest hashes the current engine and existing CNFs "
                "independently; it does not prove that the engine emitted "
                "those CNFs.  Theorem-grade closure still requires a "
                "byte-regeneration/provenance/cover audit plus per-cube "
                "DRAT emission and drat-trim VERIFIED.")
            return EXIT_OK
        pending = [c for c in named if c not in terminal]
        listed = (", ".join(pending) if len(pending) <= 20
                  else f"{len(pending)} lanes")
        self._say(
            f"INCOMPLETE: {len(terminal)}/{total} cubes have a terminal "
            f"verdict (UNSAT={len(unsat)} SAT={len(sat)}); without one: "
            f"{listed}.")
        self._say(
            "No result is claimed.  Rerun to resume; nonterminal lanes "
            "retry automatically.")
        return EXIT_INCOMPLETE

    # -- entry ----------------------------------------------------------------

    def run(self) -> int:
        self._lock = CampaignLock(self.lock_path)
        self._lock.acquire()  # Refused propagates; never unlink lock file
        try:
            return self._run_locked()
        finally:
            # Stop spawning before the final sweep.  Keep the campaign lock
            # held until every registered process group has been terminated
            # and reaped, then deterministically release it for in-process
            # callers as well as the normal exec-wrapper lifecycle.
            self.stop.set()
            try:
                self._sweep()
            finally:
                try:
                    self._restore_signal_handlers()
                finally:
                    self._lock.release()

    def _run_locked(self) -> int:
        raw = self._pin_manifest()
        self._parse_manifest(raw)
        self.ledger_path, self.logs_dir = v2_namespace(self.here, self.digest)
        self.logs_dir.mkdir(parents=True, exist_ok=True)
        rows, terminal, attempts = self._read_ledger()
        self.terminal = terminal
        self.attempts = attempts
        self.row_count = len(rows)
        self.discovered_sat = any(
            row["verdict"] == "SAT" for row in self.terminal.values())
        self._ensure_ledger_file()
        self._append_prelaunch(resumed=len(self.terminal))
        pending = [c for c in self.case_order if c not in self.terminal]
        self._say_banner(resumed=len(self.terminal), pending=len(pending))
        self._install_signal_handlers()
        self._run_workers()
        self._say("")
        self._say(f"=== v2 LEDGER {self.ledger_path.name} ===")
        unsat = sum(1 for r in self.terminal.values()
                    if r["verdict"] == "UNSAT")
        sat = sum(1 for r in self.terminal.values() if r["verdict"] == "SAT")
        self._say(
            f"attempts={self.row_count} terminal={len(self.terminal)} "
            f"(UNSAT={unsat} SAT={sat})")
        code = self._decide()
        self._say("CUBES_DONE")
        return code


def _env_int(name: str, default: int) -> int:
    raw = (os.environ.get(name) or "").strip()
    if not raw:
        return default
    try:
        return int(raw)
    except ValueError:
        raise Refused(
            f"environment {name}={raw!r} is not an integer") from None


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Supervised, manifest-pinned n=12 target-39 cube "
                    "campaign runner (v2).")
    parser.add_argument(
        "--jobs", type=int, default=None,
        help="max concurrent Kissat lanes (env JOBS, default 7)")
    parser.add_argument(
        "--timeout", type=int, default=None,
        help="per-lane wall budget seconds, 0 unlimited (env TIMEOUT)")
    args = parser.parse_args(argv)
    try:
        jobs = (args.jobs if args.jobs is not None
                else _env_int("JOBS", DEFAULT_JOBS))
        timeout = (args.timeout if args.timeout is not None
                   else _env_int("TIMEOUT", 0))
        return Supervisor(HERE, jobs=jobs, timeout=timeout).run()
    except Refused as exc:
        print(f"REFUSED: {exc}", file=sys.stderr)
        return EXIT_REFUSED
    except make_manifest.ManifestFailure as exc:
        print(f"REFUSED: manifest generation failed: {exc}", file=sys.stderr)
        return EXIT_REFUSED


if __name__ == "__main__":
    sys.exit(main())
