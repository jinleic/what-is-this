#!/usr/bin/env python3
"""Produce manifest-bound binary-DRAT certificates for Kobon n=12 cubes.

This is independent of the proof-free discovery run.  It reads only the
current decision manifest and never imports legacy ledgers or marker files.
"""
from __future__ import annotations
import argparse
import atexit
import contextlib
import errno
import fcntl
import hashlib
import json
import math
import os
import re
import resource
import secrets
import shutil
import signal
import stat
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterator, Sequence
if str(Path(__file__).resolve().parent) not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_manifest
WORKSPACE = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
EXPECTED_DECISION_SHA256 = '48671f44f726dc4272438d0797beff936878c39b9b8d26ec110546b1136db7cf'
DECISION_SCHEMA = 'kobon-n12-decision-manifest/2'
CERT_SCHEMA = 'kobon-n12-unsat-certificate-manifest/1'
LEDGER_SCHEMA = 'kobon-n12-unsat-certificate-ledger/1'
ROW_SCHEMA = 'kobon-n12-unsat-certificate-attempt/2'
RECHECK_SCHEMA = 'kobon-n12-unsat-certificate-recheck/1'
STATUS_SCHEMA = 'kobon-n12-unsat-certificate-status/1'
EXIT_OK, EXIT_INCOMPLETE, EXIT_REFUSED = (0, 1, 2)
# Preserve a 64 GiB host reserve while one compressed proof grows.
DEFAULT_DISK_FLOOR_GIB = 64.0
DEFAULT_MEMORY_LIMIT_GIB = 8.0
DEFAULT_DISK_POLL_SECONDS = 5.0
DEFAULT_PROCESS_POLL_SECONDS = 1.0
DEFAULT_KILL_GRACE_SECONDS = 10.0
SUCCESS = frozenset({'CERTIFIED_UNSAT', 'VERIFIED_INPUT_UP_UNSAT'})
OUTCOMES = SUCCESS | frozenset({'CHECKER_EXIT_FAILED', 'CHECKER_REJECTED', 'DISK_FLOOR_ABORT', 'INPUT_UP_PROBE_FAILED', 'INTERRUPTED', 'MEMORY_LIMIT_ABORT', 'ORPHANED_ATTEMPT', 'PROOF_MISSING_OR_EMPTY', 'RESOURCE_WAIT_TIMEOUT', 'ROGUE_CAKEPB_ABORT', 'SAT_NOT_CERTIFIED', 'SOLVER_EVIDENCE_MISMATCH', 'SOLVER_FAILED', 'SOLVER_LAUNCH_FAILED', 'SOLVER_RESOURCE_LIMIT', 'SOLVER_TIMEOUT', 'XZ_INTEGRITY_FAILED', 'XZ_STREAM_FAILED'})
EVIDENCE = {'CERTIFIED_UNSAT': 'PROOF_CONSUMING_BINARY_DRAT_XZ', 'VERIFIED_INPUT_UP_UNSAT': 'INPUT_UNIT_PROPAGATION_UNSAT'}
_SHA_RE = re.compile('[0-9a-f]{64}\\Z')
_ATTEMPT_RE = re.compile('[A-Za-z0-9][A-Za-z0-9._-]{7,127}\\Z')
_MEMORY_ERROR_RE = re.compile(b'out of memory|memory allocation|cannot allocate', re.IGNORECASE)
_HANDLED_SIGNALS = (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)

class Refused(RuntimeError):
    """The launch or persisted evidence is unsafe."""

class LockBusy(Refused):
    """A permanent fcntl lock is held by another process."""

class ResourceWaitTimeout(Refused):
    """The configured shared-resource wait expired."""

class RunInterrupted(RuntimeError):

    def __init__(self, signum: int):
        super().__init__(f'interrupted by signal {signum}')
        self.signum = signum

@dataclass(frozen=True)
class Config:
    workspace: Path
    here: Path
    manifest_path: Path
    checker_path: Path
    driver_lock_path: Path
    host_lock_path: Path
    expected_decision_sha256: str
    expected_cases: tuple[str, ...]
    xz_path: Path | None = None
    pgrep_path: Path | None = None
    ps_path: Path | None = None
    disk_floor_bytes: int = int(DEFAULT_DISK_FLOOR_GIB * (1 << 30))
    memory_limit_bytes: int | None = int(DEFAULT_MEMORY_LIMIT_GIB * (1 << 30))
    disk_poll_seconds: float = DEFAULT_DISK_POLL_SECONDS
    process_poll_seconds: float = DEFAULT_PROCESS_POLL_SECONDS
    solver_timeout_seconds: float = 0.0
    checker_timeout_seconds: float = 0.0
    resource_wait_seconds: float = 0.0
    kill_grace_seconds: float = DEFAULT_KILL_GRACE_SECONDS

def default_config() -> Config:
    return Config(workspace=WORKSPACE, here=HERE, manifest_path=HERE / 'manifest.json', checker_path=WORKSPACE / 'scratch/kobon-audit/tools/drat-trim/drat-trim', driver_lock_path=HERE / '.certify-unsat.lock', host_lock_path=WORKSPACE / 'scratch/.host-heavy-job.lock', expected_decision_sha256=EXPECTED_DECISION_SHA256, expected_cases=tuple(make_manifest.CUBES))

class FileLock:
    """Exclusive nonblocking lock on a permanent file."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self.stream = None

    @property
    def fd(self) -> int | None:
        return None if self.stream is None else self.stream.fileno()

    def acquire(self) -> None:
        if self.stream is not None:
            raise RuntimeError(f'lock already held: {self.path}')
        self.path.parent.mkdir(parents=True, exist_ok=True)
        stream = self.path.open('a', encoding='utf-8')
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as error:
            stream.close()
            if error.errno in {errno.EACCES, errno.EAGAIN, errno.EWOULDBLOCK}:
                raise LockBusy(f'lock held: {self.path}') from error
            raise
        self.stream = stream

    def release(self) -> None:
        stream, self.stream = (self.stream, None)
        if stream is None:
            return
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
        finally:
            stream.close()

    def __enter__(self) -> 'FileLock':
        self.acquire()
        return self

    def __exit__(self, *_exception) -> None:
        self.release()

def canonical_json(document: object) -> bytes:
    try:
        text = json.dumps(document, allow_nan=False, ensure_ascii=True, separators=(',', ':'), sort_keys=True)
    except (TypeError, ValueError) as error:
        raise Refused(f'non-canonical JSON value: {error}') from error
    return (text + '\n').encode('ascii')

def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        while (block := stream.read(1 << 20)):
            digest.update(block)
    return digest.hexdigest()

def _file(path: Path, *, executable: bool=False) -> Path:
    try:
        path = Path(path).resolve(strict=True)
        mode = path.stat().st_mode
    except OSError as error:
        raise Refused(f'missing file {path}: {error}') from error
    if not stat.S_ISREG(mode):
        raise Refused(f'not a regular file: {path}')
    if executable and (not os.access(path, os.X_OK)):
        raise Refused(f'not executable: {path}')
    return path

def _directory(path: Path) -> Path:
    path = Path(path)
    if path.exists() or path.is_symlink():
        if path.is_symlink() or not path.is_dir():
            raise Refused(f'not a plain artifact directory: {path}')
    else:
        path.mkdir(parents=True)
    return path.resolve()

def _portable(path: Path, workspace: Path) -> str:
    path = path.resolve()
    try:
        return str(path.relative_to(workspace.resolve()))
    except ValueError:
        return str(path)

def file_record(path: Path, workspace: Path, *, executable: bool=False) -> dict:
    path = _file(path, executable=executable)
    return {'bytes': path.stat().st_size, 'path': _portable(path, workspace), 'sha256': sha256_file(path)}

def _fsync(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)

def atomic_write(path: Path, data: bytes) -> None:
    """Sibling temp, file fsync, replace, and directory fsync."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(f'.{path.name}.tmp-{os.getpid()}-{secrets.token_hex(8)}')
    descriptor = None
    try:
        descriptor = os.open(temp, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC, 384)
        view = memoryview(data)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise OSError('short write')
            view = view[written:]
        os.fsync(descriptor)
        os.close(descriptor)
        descriptor = None
        os.replace(temp, path)
        _fsync(path.parent)
    except BaseException:
        if descriptor is not None:
            os.close(descriptor)
        # Preserve even an incomplete atomic-write sibling for explicit audit.
        raise

def utc(ns: int) -> str:
    seconds, nanos = divmod(ns, 1000000000)
    prefix = time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime(seconds))
    return f'{prefix}.{nanos:09d}Z'

def inspect_dimacs(path: Path) -> tuple[int, int, int]:
    header = None
    clauses_seen = 0
    maximum = 0
    clause_open = False
    with path.open('rb') as stream:
        for line_number, raw in enumerate(stream, 1):
            fields = raw.strip().split()
            if not fields or fields[0] == b'c':
                continue
            if header is None:
                if len(fields) != 4 or fields[:2] != [b'p', b'cnf']:
                    raise Refused(f'{path}:{line_number}: invalid DIMACS header')
                try:
                    header = (int(fields[2]), int(fields[3]))
                except ValueError as error:
                    raise Refused(f'{path}:{line_number}: invalid DIMACS size') from error
                if min(header) < 0:
                    raise Refused(f'{path}:{line_number}: negative DIMACS size')
                continue
            if fields[:2] == [b'p', b'cnf']:
                raise Refused(f'{path}:{line_number}: duplicate DIMACS header')
            for token in fields:
                try:
                    literal = int(token)
                except ValueError as error:
                    raise Refused(f'{path}:{line_number}: invalid DIMACS token') from error
                if literal:
                    clause_open = True
                    maximum = max(maximum, abs(literal))
                else:
                    clause_open = False
                    clauses_seen += 1
    if header is None:
        raise Refused(f'{path}: missing DIMACS header')
    variables, clauses = header
    if clause_open or clauses_seen != clauses or maximum > variables:
        raise Refused(f'{path}: DIMACS mismatch: declared ({variables}, {clauses}), found max variable {maximum} and {clauses_seen} clauses')
    return (variables, clauses, maximum)

def _status_lines(data: bytes) -> list[bytes]:
    return [line for line in data.replace(b'\r', b'\n').split(b'\n') if line.startswith(b's ')]

def solver_verdict(data: bytes) -> str | None:
    return {(b's SATISFIABLE',): 'SAT', (b's UNSATISFIABLE',): 'UNSAT', (b's UNKNOWN',): 'UNKNOWN'}.get(tuple(_status_lines(data)))

def checker_verdict(data: bytes) -> str | None:
    return {(b's VERIFIED',): 'VERIFIED', (b's NOT VERIFIED',): 'NOT_VERIFIED', (b's DERIVATION',): 'DERIVATION'}.get(tuple(_status_lines(data)))

def _has_line(data: bytes, line: bytes) -> bool:
    return line in data.replace(b'\r', b'\n').split(b'\n')

def _memory_limit(limit: int | None):
    # Darwin maps a process above ordinary finite RLIMIT_AS values before
    # exec, so lowering the limit in preexec_fn fails every launch.  The
    # parent still enforces the same byte limit over aggregate group RSS.
    if (
        limit is None
        or sys.platform == 'darwin'
        or not hasattr(resource, 'RLIMIT_AS')
    ):
        return None

    def apply() -> None:
        _soft, hard = resource.getrlimit(resource.RLIMIT_AS)
        target = limit if hard == resource.RLIM_INFINITY else min(limit, hard)
        resource.setrlimit(resource.RLIMIT_AS, (target, hard))
    return apply

class Certifier:

    def __init__(self, config: Config | None=None, *, manifest_builder: Callable[[], dict]=make_manifest.build, disk_free: Callable[[Path], int] | None=None, out=None):
        self.config = config or default_config()
        self.manifest_builder = manifest_builder
        self.disk_free = disk_free or (lambda path: shutil.disk_usage(path).free)
        self.out = out if out is not None else sys.stdout
        self._validate_config()
        self.workspace = self.config.workspace.resolve()
        self.here = self.config.here.resolve()
        self.driver_lock = FileLock(self.config.driver_lock_path)
        self.host_lock: FileLock | None = None
        self.active: dict[int, subprocess.Popen] = {}
        self.groups: set[int] = set()
        self.stop_signal: int | None = None
        self.handlers: dict[int, object] = {}
        self.atexit_registered = False
        self.decision = None
        self.decision_sha = None
        self.tools: dict[str, dict] = {}
        self.child_env: dict[str, str] = {}
        self.child_env_sha = None
        self.bindings: dict[str, dict] = {}
        self.certificate_sha = None
        self.namespace = self.proof_dir = self.log_dir = None
        self.start_dir = self.recheck_dir = self.ledger_path = None
        self.ledger = None
        self.ledger_bytes = None

    def _validate_config(self) -> None:
        c = self.config
        if not c.expected_cases or len(set(c.expected_cases)) != len(c.expected_cases):
            raise Refused('expected case order is empty or duplicate')
        if not _SHA_RE.fullmatch(c.expected_decision_sha256):
            raise Refused('malformed expected manifest digest')
        if c.disk_floor_bytes < 0:
            raise Refused('disk floor must be nonnegative')
        if c.memory_limit_bytes is not None and c.memory_limit_bytes <= 0:
            raise Refused('memory limit must be positive or unlimited')
        numbers = {'disk poll': c.disk_poll_seconds, 'process poll': c.process_poll_seconds, 'solver timeout': c.solver_timeout_seconds, 'checker timeout': c.checker_timeout_seconds, 'resource wait': c.resource_wait_seconds, 'kill grace': c.kill_grace_seconds}
        for name, value in numbers.items():
            if not math.isfinite(value) or value < 0:
                raise Refused(f'{name} must be finite and nonnegative')
        if not c.disk_poll_seconds or not c.process_poll_seconds:
            raise Refused('poll intervals must be positive')

    def _lock_fds(self) -> tuple[int, ...]:
        descriptors = [self.driver_lock.fd]
        if self.host_lock is not None:
            descriptors.append(self.host_lock.fd)
        return tuple((fd for fd in descriptors if fd is not None))

    def _spawn(self, argv: Sequence[str], **kwargs) -> subprocess.Popen:
        requested_preexec = kwargs.pop('preexec_fn', None)
        previous_mask = signal.pthread_sigmask(signal.SIG_BLOCK, _HANDLED_SIGNALS)

        def prepare_child() -> None:
            signal.pthread_sigmask(signal.SIG_UNBLOCK, _HANDLED_SIGNALS)
            if requested_preexec is not None:
                requested_preexec()
        try:
            self._check_signal()
            process = subprocess.Popen(
                list(argv), start_new_session=True, close_fds=True,
                pass_fds=self._lock_fds(), preexec_fn=prepare_child, **kwargs
            )
            # Register before unmasking: a delivered signal can now terminate
            # both the new process group and every inherited-lock holder.
            self.active[process.pid] = process
            self.groups.add(process.pid)
        finally:
            signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)
        self._check_signal()
        return process

    @staticmethod
    def _group_exists(pgid: int) -> bool:
        try:
            os.killpg(pgid, 0)
        except ProcessLookupError:
            return False
        except PermissionError:
            return True
        return True

    def _finish_group(self, pgid: int) -> None:
        if not self._group_exists(pgid):
            self.groups.discard(pgid)
            return
        try:
            os.killpg(pgid, signal.SIGTERM)
        except ProcessLookupError:
            self.groups.discard(pgid)
            return
        deadline = time.monotonic() + self.config.kill_grace_seconds
        while self._group_exists(pgid) and time.monotonic() < deadline:
            time.sleep(min(0.05, max(0.001, deadline - time.monotonic())))
        if self._group_exists(pgid):
            try:
                os.killpg(pgid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        self.groups.discard(pgid)

    def _complete(self, process: subprocess.Popen) -> None:
        self.active.pop(process.pid, None)
        self._finish_group(process.pid)

    def _signal_groups(self, signum: int) -> None:
        for pgid in list(self.groups):
            try:
                os.killpg(pgid, signum)
            except ProcessLookupError:
                pass

    def _terminate_all(self) -> None:
        groups = list(self.groups)
        processes = list(self.active.values())
        for pgid in groups:
            try:
                os.killpg(pgid, signal.SIGTERM)
            except ProcessLookupError:
                pass
        deadline = time.monotonic() + self.config.kill_grace_seconds
        for process in processes:
            try:
                process.wait(timeout=max(0.0, deadline - time.monotonic()))
            except subprocess.TimeoutExpired:
                pass
        while any((self._group_exists(pgid) for pgid in groups)) and time.monotonic() < deadline:
            time.sleep(min(0.05, max(0.001, deadline - time.monotonic())))
        for pgid in groups:
            if self._group_exists(pgid):
                try:
                    os.killpg(pgid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
        for process in processes:
            try:
                process.wait()
            except BaseException:
                pass
            self.active.pop(process.pid, None)
        self.groups.difference_update(groups)

    def _install_signals(self) -> None:

        def stop(signum, _frame) -> None:
            if self.stop_signal is None:
                self.stop_signal = signum
            self._signal_groups(signal.SIGTERM)
        for signum in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
            self.handlers[signum] = signal.signal(signum, stop)
        atexit.register(self._shutdown)
        self.atexit_registered = True

    def _restore_signals(self) -> None:
        for signum, handler in self.handlers.items():
            signal.signal(signum, handler)
        self.handlers.clear()
        if self.atexit_registered:
            atexit.unregister(self._shutdown)
            self.atexit_registered = False

    def _shutdown(self) -> None:
        self._terminate_all()
        self.driver_lock.release()

    def _check_signal(self) -> None:
        if self.stop_signal is not None:
            raise RunInterrupted(self.stop_signal)

    def _capture(self, argv: Sequence[str], *, env=None) -> tuple[int, bytes, bytes]:
        process = self._spawn(argv, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
        while True:
            try:
                stdout, stderr = process.communicate(timeout=0.25)
                break
            except subprocess.TimeoutExpired:
                if self.stop_signal is not None:
                    self._terminate_all()
        self._complete(process)
        self._check_signal()
        return (process.returncode, stdout, stderr)

    def _resolve_tool(self, name: str, configured: Path | None) -> Path:
        if configured is not None:
            return _file(configured, executable=True)
        found = shutil.which(name)
        if found is None:
            raise Refused(f'{name} not found during certificate construction')
        return _file(Path(found), executable=True)

    def _tool_record(self, path: Path, version_args: Sequence[str] | None) -> dict:
        path = _file(path, executable=True)
        record = file_record(path, self.workspace, executable=True)
        record['path'] = str(path)
        if version_args is None:
            record.update({'version': None, 'version_argv': None, 'version_stdout_sha256': None, 'version_stderr_sha256': None})
            return record
        argv = [str(path), *version_args]
        code, stdout, stderr = self._capture(argv, env=self.child_env)
        if code or not stdout.strip():
            raise Refused(f'version command failed for {path}: {code}')
        record.update({'version': stdout.decode('utf-8', 'strict').strip(), 'version_argv': argv, 'version_stdout_sha256': sha256_bytes(stdout), 'version_stderr_sha256': sha256_bytes(stderr)})
        return record

    def _verify_tool(self, name: str) -> None:
        expected = self.tools[name]
        actual = file_record(Path(expected['path']), self.workspace, executable=True)
        if actual['bytes'] != expected['bytes'] or actual['sha256'] != expected['sha256']:
            raise Refused(f'pinned tool changed: {name}')

    def _verify_tools(self) -> None:
        for name in sorted(self.tools):
            self._verify_tool(name)
        found = shutil.which('xz', path=self.child_env['PATH'])
        if found is None or Path(found).resolve() != Path(self.tools['xz']['path']):
            raise Refused('sanitized PATH does not resolve pinned XZ')

    def _verify_cnf(self, case: str) -> None:
        binding = self.bindings[case]
        path = Path(binding['resolved_path'])
        actual = file_record(path, self.workspace)
        if actual['bytes'] != binding['bytes'] or actual['sha256'] != binding['sha256']:
            raise Refused(f'CNF bytes changed: {case}')
        if inspect_dimacs(path) != (binding['variables'], binding['clauses'], binding['max_variable_seen']):
            raise Refused(f'CNF structure changed: {case}')

    def _load_decision(self) -> None:
        if self.config.manifest_path.is_symlink():
            raise Refused('decision manifest must not be a symlink')
        rebuilt = self.manifest_builder()
        regenerated = (json.dumps(rebuilt, indent=2, sort_keys=True) + '\n').encode()
        try:
            stored = self.config.manifest_path.read_bytes()
        except OSError as error:
            raise Refused(f'cannot read decision manifest: {error}') from error
        digest = sha256_bytes(regenerated)
        if stored != regenerated:
            raise Refused('decision manifest is stale: rebuilt bytes differ')
        if digest != self.config.expected_decision_sha256:
            raise Refused(f'wrong decision manifest digest: {digest}')
        order, instances = (rebuilt.get('case_order'), rebuilt.get('instances'))
        if rebuilt.get('schema') != DECISION_SCHEMA:
            raise Refused('wrong decision manifest schema')
        if order != list(self.config.expected_cases):
            raise Refused('decision manifest case_order is not the required split')
        if rebuilt.get('case_count') != len(order):
            raise Refused('decision manifest case_count mismatch')
        if not isinstance(instances, list) or [item.get('cube') if isinstance(item, dict) else None for item in instances] != order:
            raise Refused('instances do not exactly follow case_order')
        bindings = {}
        for case, item in zip(order, instances):
            if not isinstance(item, dict) or set(item) != {'bytes', 'clauses', 'cube', 'path', 'sha256', 'variables'}:
                raise Refused(f'malformed instance: {case}')
            relative = Path(item['path'])
            if relative.is_absolute() or '..' in relative.parts:
                raise Refused(f'non-relative CNF path: {case}')
            candidate = self.workspace / relative
            if candidate.is_symlink():
                raise Refused(f'CNF path is a symlink: {case}')
            path = _file(candidate)
            try:
                path.relative_to(self.workspace)
            except ValueError as error:
                raise Refused(f'CNF escapes workspace: {case}') from error
            if any((not isinstance(item[key], int) or isinstance(item[key], bool) or item[key] < 0 for key in ('bytes', 'variables', 'clauses'))) or (not isinstance(item['sha256'], str) or not _SHA_RE.fullmatch(item['sha256'])):
                raise Refused(f'malformed CNF metadata: {case}')
            actual = file_record(path, self.workspace)
            variables, clauses, maximum = inspect_dimacs(path)
            if actual['bytes'] != item['bytes'] or actual['sha256'] != item['sha256']:
                raise Refused(f'CNF size/hash mismatch: {case}')
            if (variables, clauses) != (item['variables'], item['clauses']):
                raise Refused(f'CNF header mismatch: {case}')
            bindings[case] = {'bytes': item['bytes'], 'clauses': item['clauses'], 'manifest_path': item['path'], 'max_variable_seen': maximum, 'resolved_path': str(path), 'sha256': item['sha256'], 'variables': item['variables']}
        self.decision, self.decision_sha, self.bindings = (rebuilt, digest, bindings)

    def _build_certificate_manifest(self) -> None:
        assert self.decision is not None and self.decision_sha is not None
        xz_path = self._resolve_tool('xz', self.config.xz_path)
        self.child_env = {'LANG': 'C', 'LC_ALL': 'C', 'PATH': str(xz_path.parent), 'TZ': 'UTC'}
        self.child_env_sha = sha256_bytes(canonical_json(self.child_env))
        found = shutil.which('xz', path=self.child_env['PATH'])
        if found is None or Path(found).resolve() != xz_path:
            raise Refused('sanitized producer PATH does not resolve pinned XZ')
        source_solver = self.decision.get('solver')
        if not isinstance(source_solver, dict) or set(source_solver) != {'bytes', 'path', 'sha256', 'version'}:
            raise Refused('malformed decision-manifest solver')
        solver_path = Path(source_solver['path'])
        if not solver_path.is_absolute():
            solver_path = self.workspace / solver_path
        solver = self._tool_record(solver_path, ['--version'])
        if any((solver[key] != source_solver[key] for key in ('bytes', 'sha256', 'version'))):
            raise Refused('Kissat differs from decision manifest')
        self.tools = {'checker': self._tool_record(self.config.checker_path, None), 'kissat': solver, 'pgrep': self._tool_record(self._resolve_tool('pgrep', self.config.pgrep_path), None), 'ps': self._tool_record(self._resolve_tool('ps', self.config.ps_path), None), 'xz': self._tool_record(xz_path, ['--version'])}
        document = {'case_order': list(self.config.expected_cases), 'certifier': {**file_record(Path(__file__), self.workspace), 'path': str(Path(__file__).resolve())}, 'checker_contract': {'checker_argv_suffix': ['-i', '-w'], 'input_up_probe': 'empty binary proof stdin', 'required_status_line': 's VERIFIED', 'stream_argv': ['XZ', '--format=xz', '--decompress', '--stdout', '--', 'PROOF']}, 'child_environment': self.child_env, 'child_environment_sha256': self.child_env_sha, 'decision_manifest': file_record(self.config.manifest_path, self.workspace), 'decision_manifest_sha256': self.decision_sha, 'instances': self.bindings, 'lock_contract': {'driver': str(self.config.driver_lock_path.resolve()), 'host_heavy': str(self.config.host_lock_path.resolve()), 'inherit_driver_into_producer': True, 'inherit_host_into_verifier': True}, 'proof_format': 'kissat-default-binary-drat-inside-xz', 'rss_probe_argv': ['PS', '-axo', 'pid=,pgid=,rss='], 'schema': CERT_SCHEMA, 'tools': self.tools}
        payload = canonical_json(document)
        digest = sha256_bytes(payload)
        namespace = _directory(self.here / f'cert.v1.{digest}')
        self.proof_dir = _directory(namespace / 'proofs')
        self.log_dir = _directory(namespace / 'logs')
        self.start_dir = _directory(namespace / 'attempt-starts')
        self.recheck_dir = _directory(namespace / 'resume-checks')
        manifest_path = namespace / 'certificate-manifest.json'
        if manifest_path.is_symlink():
            raise Refused('certificate manifest must not be a symlink')
        if manifest_path.exists():
            if manifest_path.read_bytes() != payload:
                raise Refused('foreign certificate manifest in digest namespace')
        else:
            atomic_write(manifest_path, payload)
        _fsync(namespace)
        _fsync(namespace.parent)
        self.certificate_sha, self.namespace = (digest, namespace)
        self.ledger_path = namespace / 'ledger.json'

    def prepare(self) -> None:
        self._load_decision()
        self._build_certificate_manifest()
        self._verify_tools()
        self._load_ledger()
        self._scan_artifacts()

    def _artifact(self, path: Path) -> dict:
        _fsync(path)
        _fsync(path.parent)
        return file_record(path, self.workspace)

    def _artifact_path(self, record: dict) -> Path:
        if not isinstance(record, dict) or set(record) != {'bytes', 'path', 'sha256'}:
            raise Refused('malformed artifact record')
        if not isinstance(record['bytes'], int) or isinstance(record['bytes'], bool) or record['bytes'] < 0 or (not isinstance(record['path'], str)) or (not isinstance(record['sha256'], str)) or (not _SHA_RE.fullmatch(record['sha256'])):
            raise Refused('malformed artifact metadata')
        raw = Path(record['path'])
        if not raw.is_absolute():
            raw = self.workspace / raw
        if raw.is_symlink():
            raise Refused(f'artifact is a symlink: {raw}')
        path = _file(raw)
        assert self.namespace is not None
        try:
            path.relative_to(self.namespace)
        except ValueError as error:
            raise Refused(f'artifact escapes namespace: {path}') from error
        if path.stat().st_size != record['bytes'] or sha256_file(path) != record['sha256']:
            raise Refused(f'artifact bytes changed: {path}')
        return path

    def _new_ledger(self) -> dict:
        return {'certificate_manifest_sha256': self.certificate_sha, 'decision_manifest_sha256': self.decision_sha, 'rows': [], 'schema': LEDGER_SCHEMA}

    def _load_ledger(self) -> None:
        assert self.ledger_path is not None
        if self.ledger_path.is_symlink():
            raise Refused('ledger must not be a symlink')
        if not self.ledger_path.exists():
            self.ledger, self.ledger_bytes = (self._new_ledger(), None)
            return
        payload = self.ledger_path.read_bytes()
        try:
            document = json.loads(payload)
        except (UnicodeDecodeError, json.JSONDecodeError) as error:
            raise Refused(f'corrupt ledger: {error}') from error
        if not isinstance(document, dict) or canonical_json(document) != payload:
            raise Refused('ledger is not canonical JSON')
        if set(document) != {'certificate_manifest_sha256', 'decision_manifest_sha256', 'rows', 'schema'}:
            raise Refused('foreign ledger fields')
        if document['schema'] != LEDGER_SCHEMA or document['certificate_manifest_sha256'] != self.certificate_sha or document['decision_manifest_sha256'] != self.decision_sha or (not isinstance(document['rows'], list)):
            raise Refused('foreign ledger binding')
        attempts, terminals = (set(), set())
        for row in document['rows']:
            self._validate_row(row)
            if row['attempt_id'] in attempts:
                raise Refused(f"duplicate attempt: {row['attempt_id']}")
            attempts.add(row['attempt_id'])
            if row['terminal']:
                if row['case'] in terminals:
                    raise Refused(f"duplicate terminal row: {row['case']}")
                terminals.add(row['case'])
        self.ledger, self.ledger_bytes = (document, payload)

    @staticmethod
    def _time_ok(item: dict) -> bool:
        start = item.get('started_unix_ns')
        finish = item.get('finished_unix_ns')
        duration = item.get('duration_seconds')
        return isinstance(start, int) and (not isinstance(start, bool)) and isinstance(finish, int) and (not isinstance(finish, bool)) and (finish >= start) and (item.get('started_utc') == utc(start)) and (item.get('finished_utc') == utc(finish)) and isinstance(duration, (int, float)) and (not isinstance(duration, bool)) and math.isfinite(duration) and (duration >= 0)

    def _phase_logs(self, phase: object, pipeline: bool=False):
        if phase is None:
            return None
        common = {'duration_seconds', 'finished_unix_ns', 'finished_utc', 'launch_error', 'started_unix_ns', 'started_utc', 'stop_reason'}
        expected = common | ({'checker_argv', 'checker_exit_code', 'checker_log', 'xz_argv', 'xz_exit_code', 'xz_log'} if pipeline else {'argv', 'exit_code', 'log'})
        if not isinstance(phase, dict) or set(phase) != expected or (not self._time_ok(phase)):
            raise Refused('malformed phase')
        for key in ('launch_error', 'stop_reason'):
            if phase[key] is not None and (not isinstance(phase[key], str)):
                raise Refused(f'malformed phase {key}')
        argv_keys = ('checker_argv', 'xz_argv') if pipeline else ('argv',)
        exit_keys = ('checker_exit_code', 'xz_exit_code') if pipeline else ('exit_code',)
        for key in argv_keys:
            if not isinstance(phase[key], list) or not all((isinstance(value, str) for value in phase[key])):
                raise Refused(f'malformed phase {key}')
        for key in exit_keys:
            if phase[key] is not None and (not isinstance(phase[key], int) or isinstance(phase[key], bool)):
                raise Refused(f'malformed phase {key}')
        if pipeline:
            return (self._artifact_path(phase['xz_log']), self._artifact_path(phase['checker_log']))
        return (self._artifact_path(phase['log']),)

    def _paths(self, case: str, attempt: str) -> dict[str, Path]:
        assert self.proof_dir is not None and self.log_dir is not None
        prefix = f'{case}.{attempt}'
        return {'proof': self.proof_dir / f'{prefix}.drat.xz', 'solver': self.log_dir / f'{prefix}.solver.log', 'xz_integrity': self.log_dir / f'{prefix}.xz-integrity.log', 'input_up': self.log_dir / f'{prefix}.input-up.log', 'xz_stream': self.log_dir / f'{prefix}.xz-stream.log', 'checker': self.log_dir / f'{prefix}.checker.log'}

    def _start_path(self, case: str, attempt: str) -> Path:
        assert self.start_dir is not None
        return self.start_dir / f'{case}.{attempt}.json'
    def _artifact_identity(self, path: Path) -> tuple[str, str] | None:
        assert self.start_dir is not None and self.proof_dir is not None
        assert self.log_dir is not None
        parent = path.parent.resolve()
        if parent == self.start_dir.resolve():
            suffixes = ('.json',)
        elif parent == self.proof_dir.resolve():
            suffixes = ('.drat.xz',)
        elif parent == self.log_dir.resolve():
            suffixes = (
                '.solver.log', '.xz-integrity.log', '.input-up.log',
                '.xz-stream.log', '.checker.log',
            )
        else:
            return None
        for case in sorted(self.config.expected_cases, key=len, reverse=True):
            prefix = f'{case}.'
            if not path.name.startswith(prefix):
                continue
            for suffix in suffixes:
                if path.name.endswith(suffix):
                    attempt = path.name[len(prefix):-len(suffix)]
                    if _ATTEMPT_RE.fullmatch(attempt):
                        return case, attempt
        return None

    def _options(self) -> dict:
        c = self.config
        return {'checker_timeout_seconds': c.checker_timeout_seconds, 'disk_floor_bytes': c.disk_floor_bytes, 'disk_poll_seconds': c.disk_poll_seconds, 'kill_grace_seconds': c.kill_grace_seconds, 'memory_limit_bytes': c.memory_limit_bytes, 'process_poll_seconds': c.process_poll_seconds, 'resource_wait_seconds': c.resource_wait_seconds, 'solver_timeout_seconds': c.solver_timeout_seconds}

    @staticmethod
    def _options_ok(options: object) -> bool:
        keys = {'checker_timeout_seconds', 'disk_floor_bytes', 'disk_poll_seconds', 'kill_grace_seconds', 'memory_limit_bytes', 'process_poll_seconds', 'resource_wait_seconds', 'solver_timeout_seconds'}
        if not isinstance(options, dict) or set(options) != keys:
            return False
        if not isinstance(options['disk_floor_bytes'], int) or isinstance(options['disk_floor_bytes'], bool) or options['disk_floor_bytes'] < 0:
            return False
        memory = options['memory_limit_bytes']
        if memory is not None and (not isinstance(memory, int) or isinstance(memory, bool) or memory <= 0):
            return False
        for key in keys - {'disk_floor_bytes', 'memory_limit_bytes'}:
            value = options[key]
            if not isinstance(value, (int, float)) or isinstance(value, bool) or (not math.isfinite(value)) or (value < 0):
                return False
        return bool(options['disk_poll_seconds'] and options['process_poll_seconds'])

    def _classify(self, input_up: dict, proof_check: dict):
        input_data = self._artifact_path(input_up['log']).read_bytes()
        proof_data = self._artifact_path(proof_check['checker_log']).read_bytes()
        proof_verdict = checker_verdict(proof_data)
        if proof_check['xz_exit_code'] != 0:
            return (proof_verdict, 'XZ_STREAM_FAILED')
        if proof_verdict != 'VERIFIED':
            return (proof_verdict, 'CHECKER_REJECTED')
        input_verdict = checker_verdict(input_data)
        input_is_up = input_verdict == 'VERIFIED' and (input_up['exit_code'] == 0 or (input_up['exit_code'] == 1 and _has_line(input_data, b'c trivial UNSAT')))
        if input_is_up:
            if proof_check['checker_exit_code'] not in {0, 1}:
                return (proof_verdict, 'CHECKER_EXIT_FAILED')
            if proof_check['checker_exit_code'] == 1 and (not _has_line(proof_data, b'c trivial UNSAT')):
                return (proof_verdict, 'CHECKER_EXIT_FAILED')
            return (proof_verdict, 'VERIFIED_INPUT_UP_UNSAT')
        if input_verdict != 'NOT_VERIFIED' or input_up['exit_code'] != 1:
            return (proof_verdict, 'INPUT_UP_PROBE_FAILED')
        if proof_check['checker_exit_code'] != 0:
            return (proof_verdict, 'CHECKER_EXIT_FAILED')
        return (proof_verdict, 'CERTIFIED_UNSAT')

    def _validate_start(self, row: dict, path: Path) -> None:
        payload = path.read_bytes()
        try:
            start = json.loads(payload)
        except json.JSONDecodeError as error:
            raise Refused(f'corrupt attempt-start record: {error}') from error
        if canonical_json(start) != payload:
            raise Refused('attempt-start record is not canonical')
        case, attempt = row['case'], row['attempt_id']
        expected_paths = {
            key: str(value) for key, value in self._paths(case, attempt).items()
        }
        expected = {
            'artifact_targets': expected_paths,
            'attempt_id': attempt,
            'case': case,
            'certificate_manifest_sha256': self.certificate_sha,
            'child_environment': self.child_env,
            'child_environment_sha256': self.child_env_sha,
            'cnf': self.bindings[case],
            'decision_manifest_sha256': self.decision_sha,
            'options': row['options'],
            'schema': 'kobon-n12-unsat-certificate-attempt-start/1',
            'solver_argv': [
                self.tools['kissat']['path'],
                self.bindings[case]['resolved_path'],
                expected_paths['proof'],
            ],
            'started_unix_ns': row['started_unix_ns'],
            'started_utc': row['started_utc'],
            'tools': self.tools,
        }
        if start != expected:
            raise Refused('attempt-start record differs from ledger row')

    def _validate_row(self, row: object) -> None:
        keys = {
            'attempt_id', 'case', 'certificate_manifest_sha256',
            'checker_verdict', 'child_environment',
            'child_environment_sha256', 'cnf', 'decision_manifest_sha256',
            'duration_seconds', 'evidence_class', 'finished_unix_ns',
            'finished_utc', 'input_up', 'interruption_signal', 'options',
            'orphan_artifacts', 'outcome', 'proof', 'proof_check', 'schema',
            'solver', 'solver_verdict', 'start_record', 'started_unix_ns',
            'started_utc', 'terminal', 'tools', 'xz_integrity',
        }
        if not isinstance(row, dict) or set(row) != keys:
            raise Refused('malformed ledger row')
        case, attempt, outcome = row['case'], row['attempt_id'], row['outcome']
        if not isinstance(case, str) or case not in self.bindings:
            raise Refused(f'foreign case in ledger: {case!r}')
        if not isinstance(attempt, str) or not _ATTEMPT_RE.fullmatch(attempt):
            raise Refused('malformed attempt id')
        if not isinstance(outcome, str) or outcome not in OUTCOMES:
            raise Refused('unknown outcome')
        if (
            row['schema'] != ROW_SCHEMA
            or row['decision_manifest_sha256'] != self.decision_sha
            or row['certificate_manifest_sha256'] != self.certificate_sha
            or row['child_environment'] != self.child_env
            or row['child_environment_sha256'] != self.child_env_sha
            or row['cnf'] != self.bindings[case]
            or row['tools'] != self.tools
            or row['terminal'] != (outcome in SUCCESS)
            or row['evidence_class'] != EVIDENCE.get(outcome, 'NOT_CERTIFIED')
            or not self._time_ok(row)
            or not self._options_ok(row['options'])
            or not isinstance(row['orphan_artifacts'], list)
        ):
            raise Refused('foreign or inconsistent row binding')
        if outcome == 'ORPHANED_ATTEMPT':
            if any(
                row[key] is not None for key in (
                    'checker_verdict', 'input_up', 'interruption_signal',
                    'proof', 'proof_check', 'solver', 'solver_verdict',
                    'xz_integrity',
                )
            ) or not row['orphan_artifacts']:
                raise Refused('malformed orphaned-attempt row')
            paths = [self._artifact_path(item) for item in row['orphan_artifacts']]
            if len(paths) != len(set(paths)) or any(
                self._artifact_identity(path) != (case, attempt) for path in paths
            ):
                raise Refused('orphaned-attempt artifacts are not one attempt')
            if row['start_record'] is not None:
                start_path = self._artifact_path(row['start_record'])
                if start_path != self._start_path(case, attempt).resolve():
                    raise Refused('unbound orphan attempt-start record')
                if start_path not in paths:
                    raise Refused('orphan start record missing from artifacts')
                self._validate_start(row, start_path)
            return
        if row['orphan_artifacts']:
            raise Refused('ordinary row carries orphan artifacts')
        if row['solver_verdict'] not in {None, 'SAT', 'UNSAT', 'UNKNOWN'}:
            raise Refused('malformed solver verdict')
        if row['checker_verdict'] not in {
            None, 'VERIFIED', 'NOT_VERIFIED', 'DERIVATION'
        }:
            raise Refused('malformed checker verdict')
        expected = self._paths(case, attempt)
        solver_logs = self._phase_logs(row['solver'])
        solver_argv = [
            self.tools['kissat']['path'],
            self.bindings[case]['resolved_path'],
            str(expected['proof']),
        ]
        if (
            solver_logs is None
            or row['solver']['argv'] != solver_argv
            or solver_logs[0] != expected['solver'].resolve()
        ):
            raise Refused('unbound solver phase')
        if solver_verdict(solver_logs[0].read_bytes()) != row['solver_verdict']:
            raise Refused('solver verdict differs from log')
        start_path = self._artifact_path(row['start_record'])
        if start_path != self._start_path(case, attempt).resolve():
            raise Refused('unbound attempt-start record')
        self._validate_start(row, start_path)
        proof_path = (
            None if row['proof'] is None else self._artifact_path(row['proof'])
        )
        if proof_path is not None and proof_path != expected['proof'].resolve():
            raise Refused('unbound proof path')
        integrity_logs = self._phase_logs(row['xz_integrity'])
        input_logs = self._phase_logs(row['input_up'])
        check_logs = self._phase_logs(row['proof_check'], pipeline=True)
        xz = self.tools['xz']['path']
        checker_argv = [
            self.tools['checker']['path'],
            self.bindings[case]['resolved_path'], '-i', '-w',
        ]
        integrity_argv = [
            xz, '--format=xz', '--test', '--', str(expected['proof'])
        ]
        stream_argv = [
            xz, '--format=xz', '--decompress', '--stdout', '--',
            str(expected['proof']),
        ]
        if row['xz_integrity'] is not None and (
            row['xz_integrity']['argv'] != integrity_argv
            or integrity_logs[0] != expected['xz_integrity'].resolve()
        ):
            raise Refused('unbound XZ integrity phase')
        if row['input_up'] is not None and (
            row['input_up']['argv'] != checker_argv
            or input_logs[0] != expected['input_up'].resolve()
        ):
            raise Refused('unbound input-UP phase')
        if row['proof_check'] is not None and (
            row['proof_check']['xz_argv'] != stream_argv
            or row['proof_check']['checker_argv'] != checker_argv
            or check_logs[0] != expected['xz_stream'].resolve()
            or check_logs[1] != expected['checker'].resolve()
        ):
            raise Refused('unbound proof pipeline')
        for phase in (
            row['solver'], row['xz_integrity'], row['input_up'],
            row['proof_check'],
        ):
            if phase is not None and (
                phase['started_unix_ns'] < row['started_unix_ns']
                or phase['finished_unix_ns'] > row['finished_unix_ns']
            ):
                raise Refused('phase timing outside attempt')
        if check_logs is not None:
            if checker_verdict(check_logs[1].read_bytes()) != row['checker_verdict']:
                raise Refused('checker verdict differs from log')
        elif row['checker_verdict'] is not None:
            raise Refused('checker verdict has no log')
        if outcome == 'INTERRUPTED':
            if (
                not isinstance(row['interruption_signal'], int)
                or isinstance(row['interruption_signal'], bool)
                or row['interruption_signal'] <= 0
            ):
                raise Refused('interrupted row lacks signal')
        elif row['interruption_signal'] is not None:
            raise Refused('non-interrupted row carries signal')
        self._validate_outcome(row, proof_path)

    def _validate_outcome(self, row: dict, proof_path: Path | None) -> None:
        outcome, solver = (row['outcome'], row['solver'])
        if outcome in SUCCESS:
            if proof_path is None or solver['exit_code'] != 20 or row['solver_verdict'] != 'UNSAT' or (row['xz_integrity'] is None) or (row['xz_integrity']['exit_code'] != 0) or (row['input_up'] is None) or (row['proof_check'] is None):
                raise Refused('terminal row lacks certificate evidence')
            verdict, calculated = self._classify(row['input_up'], row['proof_check'])
            if verdict != row['checker_verdict'] or calculated != outcome:
                raise Refused('terminal outcome differs from checker evidence')
            return
        expected = {'DISK_FLOOR_ABORT': solver['stop_reason'] == 'DISK_FLOOR', 'MEMORY_LIMIT_ABORT': solver['stop_reason'] == 'MEMORY_LIMIT', 'SOLVER_TIMEOUT': solver['stop_reason'] == 'TIMEOUT', 'SOLVER_LAUNCH_FAILED': solver['launch_error'] is not None, 'SAT_NOT_CERTIFIED': solver['exit_code'] == 10 and row['solver_verdict'] == 'SAT', 'SOLVER_EVIDENCE_MISMATCH': solver['exit_code'] in {10, 20} and row['solver_verdict'] != {10: 'SAT', 20: 'UNSAT'}[solver['exit_code']], 'PROOF_MISSING_OR_EMPTY': solver['exit_code'] == 20 and row['solver_verdict'] == 'UNSAT' and (proof_path is None or row['proof']['bytes'] == 0), 'SOLVER_RESOURCE_LIMIT': solver['exit_code'] is not None and solver['exit_code'] < 0 or bool(_MEMORY_ERROR_RE.search(self._artifact_path(solver['log']).read_bytes())), 'SOLVER_FAILED': solver['exit_code'] not in {10, 20} and solver['launch_error'] is None and (solver['stop_reason'] is None), 'XZ_INTEGRITY_FAILED': row['xz_integrity'] is not None and row['xz_integrity']['exit_code'] != 0, 'RESOURCE_WAIT_TIMEOUT': solver['exit_code'] == 20 and row['proof'] is not None, 'INTERRUPTED': row['interruption_signal'] is not None, 'ROGUE_CAKEPB_ABORT': any((phase is not None and phase['stop_reason'] == 'ROGUE_CAKEPB' for phase in (row['xz_integrity'], row['input_up'], row['proof_check'])))}
        if outcome in expected and (not expected[outcome]):
            raise Refused(f'outcome does not match phases: {outcome}')
        checker_failures = {'CHECKER_EXIT_FAILED', 'CHECKER_REJECTED', 'INPUT_UP_PROBE_FAILED', 'XZ_STREAM_FAILED'}
        if outcome in checker_failures:
            if row['input_up'] is None or row['proof_check'] is None:
                raise Refused('checker failure lacks phases')
            verdict, calculated = self._classify(row['input_up'], row['proof_check'])
            if calculated != outcome or verdict != row['checker_verdict']:
                raise Refused('checker failure differs from evidence')

    def _referenced_paths(self) -> set[Path]:
        assert self.ledger is not None
        result = set()
        for row in self.ledger['rows']:
            if row['start_record'] is not None:
                result.add(self._artifact_path(row['start_record']))
            result.update(
                self._artifact_path(record)
                for record in row['orphan_artifacts']
            )
            if row['proof'] is not None:
                result.add(self._artifact_path(row['proof']))
            for phase, pipeline in (
                (row['solver'], False),
                (row['xz_integrity'], False),
                (row['input_up'], False),
                (row['proof_check'], True),
            ):
                logs = self._phase_logs(phase, pipeline)
                if logs:
                    result.update(logs)
        return result

    def _recheck_artifacts(self, manifest: Path) -> set[Path]:
        if not manifest.is_file() or manifest.is_symlink():
            raise Refused(f'foreign resume-check entry: {manifest}')
        payload = manifest.read_bytes()
        try:
            document = json.loads(payload)
        except json.JSONDecodeError as error:
            raise Refused(f'corrupt resume-check manifest: {error}') from error
        keys = {
            'artifacts', 'case', 'certificate_manifest_sha256',
            'decision_manifest_sha256', 'duration_seconds', 'error',
            'finished_unix_ns', 'finished_utc', 'original_attempt_id',
            'result_outcome', 'schema', 'started_unix_ns', 'started_utc',
        }
        case = document.get('case') if isinstance(document, dict) else None
        attempt = (
            document.get('original_attempt_id')
            if isinstance(document, dict) else None
        )
        terminal = [
            row for row in self.ledger['rows']
            if row['terminal']
            and row['case'] == case
            and row['attempt_id'] == attempt
        ]
        prefix = f'{case}.{attempt}.resume-'
        token = manifest.name[len(prefix):-5] if (
            manifest.name.startswith(prefix)
            and manifest.name.endswith('.json')
        ) else ''
        if (
            not isinstance(document, dict)
            or set(document) != keys
            or canonical_json(document) != payload
            or document['schema'] != RECHECK_SCHEMA
            or document['certificate_manifest_sha256'] != self.certificate_sha
            or document['decision_manifest_sha256'] != self.decision_sha
            or len(terminal) != 1
            or not re.fullmatch('[0-9a-f]{16}', token)
            or not self._time_ok(document)
            or not isinstance(document['artifacts'], dict)
            or document['error'] is not None
            and not isinstance(document['error'], str)
            or document['result_outcome'] is not None
            and document['result_outcome'] not in OUTCOMES
        ):
            raise Refused(f'foreign resume-check manifest: {manifest}')
        suffixes = {
            'xz_integrity': '.xz-integrity.log',
            'input_up': '.input-up.log',
            'xz_stream': '.xz-stream.log',
            'checker': '.checker.log',
        }
        if not set(document['artifacts']).issubset(suffixes):
            raise Refused(f'foreign resume-check artifacts: {manifest}')
        if document['error'] is None and (
            set(document['artifacts']) != set(suffixes)
            or document['result_outcome'] != terminal[0]['outcome']
        ):
            raise Refused(f'incomplete successful resume-check: {manifest}')
        result = {manifest.resolve()}
        stem = manifest.name[:-5]
        for name, record in document['artifacts'].items():
            path = self._artifact_path(record)
            if path != (self.log_dir / f'{stem}{suffixes[name]}').resolve():
                raise Refused(f'unbound resume-check artifact: {path}')
            result.add(path)
        return result
    def _scan_artifacts(self) -> None:
        """Validate known files and durably ledger every abandoned attempt."""
        assert all((
            self.start_dir, self.proof_dir, self.log_dir, self.recheck_dir
        ))
        referenced = self._referenced_paths()
        for manifest in self.recheck_dir.iterdir():
            referenced.update(self._recheck_artifacts(manifest))
        present = set()
        for directory in (
            self.start_dir, self.proof_dir, self.log_dir, self.recheck_dir
        ):
            for path in directory.iterdir():
                if not path.is_file() or path.is_symlink():
                    raise Refused(f'foreign artifact entry: {path}')
                present.add(path.resolve())
        groups: dict[tuple[str, str], list[Path]] = {}
        for path in present - referenced:
            identity = self._artifact_identity(path)
            if identity is None:
                raise Refused(f'unrecognized preserved artifact: {path}')
            groups.setdefault(identity, []).append(path)
        existing = {
            row['attempt_id'] for row in self.ledger['rows']
        }
        for (case, attempt), paths in sorted(groups.items()):
            if attempt in existing:
                raise Refused(
                    f'unledgered files collide with attempt {attempt}: '
                    + ', '.join(map(str, sorted(paths)))
                )
            records = [
                self._artifact(path) for path in sorted(paths)
            ]
            finished = time.time_ns()
            started = min(finished, min(path.stat().st_mtime_ns for path in paths))
            options = self._options()
            start_record = None
            start_path = self._start_path(case, attempt).resolve()
            candidate = next(
                (record for record, path in zip(records, sorted(paths))
                 if path == start_path),
                None,
            )
            if candidate is not None:
                try:
                    start_document = json.loads(start_path.read_bytes())
                    if (
                        canonical_json(start_document) == start_path.read_bytes()
                        and isinstance(start_document.get('started_unix_ns'), int)
                        and not isinstance(
                            start_document.get('started_unix_ns'), bool
                        )
                        and self._options_ok(start_document.get('options'))
                    ):
                        started = min(
                            finished, start_document['started_unix_ns']
                        )
                        options = start_document['options']
                        start_record = candidate
                except (UnicodeDecodeError, json.JSONDecodeError, AttributeError):
                    pass
            row = {
                'attempt_id': attempt,
                'case': case,
                'certificate_manifest_sha256': self.certificate_sha,
                'checker_verdict': None,
                'child_environment': self.child_env,
                'child_environment_sha256': self.child_env_sha,
                'cnf': self.bindings[case],
                'decision_manifest_sha256': self.decision_sha,
                'duration_seconds': max(0.0, (finished - started) / 1e9),
                'evidence_class': 'NOT_CERTIFIED',
                'finished_unix_ns': finished,
                'finished_utc': utc(finished),
                'input_up': None,
                'interruption_signal': None,
                'options': options,
                'orphan_artifacts': records,
                'outcome': 'ORPHANED_ATTEMPT',
                'proof': None,
                'proof_check': None,
                'schema': ROW_SCHEMA,
                'solver': None,
                'solver_verdict': None,
                'start_record': start_record,
                'started_unix_ns': started,
                'started_utc': utc(started),
                'terminal': False,
                'tools': self.tools,
                'xz_integrity': None,
            }
            if start_record is not None:
                try:
                    self._validate_start(row, start_path)
                except Refused:
                    row['start_record'] = None
                    row['options'] = self._options()
            self._commit(row)
            existing.add(attempt)

    def _commit(self, row: dict) -> None:
        assert self.ledger is not None and self.ledger_path is not None
        self._validate_row(row)
        if any(
            prior['attempt_id'] == row['attempt_id']
            for prior in self.ledger['rows']
        ):
            raise Refused(f"duplicate attempt: {row['attempt_id']}")
        if row['terminal'] and any(
            prior['terminal'] and prior['case'] == row['case']
            for prior in self.ledger['rows']
        ):
            raise Refused(f"duplicate terminal row: {row['case']}")
        current = (
            self.ledger_path.read_bytes() if self.ledger_path.exists() else None
        )
        if current != self.ledger_bytes:
            raise Refused('ledger changed concurrently')
        document = {**self.ledger, 'rows': [*self.ledger['rows'], row]}
        payload = canonical_json(document)
        atomic_write(self.ledger_path, payload)
        self.ledger, self.ledger_bytes = document, payload

    def _phase(self, argv: Sequence[str], log_path: Path, *, timeout: float, disk_guard: bool=False, cake_guard: bool=False, memory_guard: bool=False, stdin=None, preexec_fn=None) -> dict:
        started_ns, started_mono = (time.time_ns(), time.monotonic())
        stop_reason = launch_error = None
        process = None
        with log_path.open('xb') as log:
            try:
                process = self._spawn(argv, stdin=stdin if stdin is not None else subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT, env=self.child_env, preexec_fn=preexec_fn)
                deadline = None if timeout == 0 else started_mono + timeout
                next_disk = started_mono + self.config.disk_poll_seconds
                next_process = started_mono
                while process.poll() is None:
                    now = time.monotonic()
                    if self.stop_signal is not None:
                        stop_reason = f'SIGNAL_{self.stop_signal}'
                        self._terminate_all()
                        break
                    if now >= next_process and (cake_guard or memory_guard):
                        next_process = now + self.config.process_poll_seconds
                        if cake_guard:
                            pids = self._cake_processes()
                            if pids:
                                stop_reason = 'ROGUE_CAKEPB'
                                log.write(('ROGUE cake_pb: ' + ' '.join(map(str, pids)) + '\n').encode())
                                self._terminate_all()
                                break
                        if memory_guard and self._group_rss(process.pid) > self.config.memory_limit_bytes:
                            stop_reason = 'MEMORY_LIMIT'
                            log.write(b'PROCESS GROUP RSS EXCEEDED LIMIT\n')
                            self._terminate_all()
                            break
                    if disk_guard and now >= next_disk:
                        next_disk = now + self.config.disk_poll_seconds
                        assert self.proof_dir is not None
                        if self.disk_free(self.proof_dir) < self.config.disk_floor_bytes:
                            stop_reason = 'DISK_FLOOR'
                            self._terminate_all()
                            break
                    if deadline is not None and now >= deadline:
                        stop_reason = 'TIMEOUT'
                        self._terminate_all()
                        break
                    wake = now + 0.25
                    if cake_guard or memory_guard:
                        wake = min(wake, next_process)
                    if disk_guard:
                        wake = min(wake, next_disk)
                    if deadline is not None:
                        wake = min(wake, deadline)
                    time.sleep(max(0.001, wake - now))
            except RunInterrupted as error:
                stop_reason = f'SIGNAL_{error.signum}'
                self._terminate_all()
            except (OSError, subprocess.SubprocessError) as error:
                launch_error = f'{type(error).__name__}: {error}'
                log.write(f'LAUNCH ERROR: {error}\n'.encode('utf-8', 'replace'))
                self._terminate_all()
            exit_code = None
            if process is not None:
                if process.poll() is None:
                    self._terminate_all()
                exit_code = process.wait()
                self._complete(process)
            if self.stop_signal is not None and stop_reason is None:
                stop_reason = f'SIGNAL_{self.stop_signal}'
            log.flush()
            os.fsync(log.fileno())
        finished_ns = time.time_ns()
        return {'argv': list(argv), 'duration_seconds': max(0.0, time.monotonic() - started_mono), 'exit_code': exit_code, 'finished_unix_ns': finished_ns, 'finished_utc': utc(finished_ns), 'launch_error': launch_error, 'log': self._artifact(log_path), 'started_unix_ns': started_ns, 'started_utc': utc(started_ns), 'stop_reason': stop_reason}

    def _pipeline(self, xz_argv: Sequence[str], checker_argv: Sequence[str], xz_log_path: Path, checker_log_path: Path) -> dict:
        started_ns, started_mono = (time.time_ns(), time.monotonic())
        stop_reason = launch_error = None
        xz_process = checker_process = None
        with xz_log_path.open('xb') as xz_log, checker_log_path.open('xb') as checker_log:
            try:
                xz_process = self._spawn(xz_argv, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=xz_log, env=self.child_env)
                assert xz_process.stdout is not None
                checker_process = self._spawn(checker_argv, stdin=xz_process.stdout, stdout=checker_log, stderr=subprocess.STDOUT, env=self.child_env)
                xz_process.stdout.close()
                deadline = None if self.config.checker_timeout_seconds == 0 else started_mono + self.config.checker_timeout_seconds
                next_probe = started_mono
                while xz_process.poll() is None or checker_process.poll() is None:
                    now = time.monotonic()
                    if self.stop_signal is not None:
                        stop_reason = f'SIGNAL_{self.stop_signal}'
                        self._terminate_all()
                        break
                    if now >= next_probe:
                        next_probe = now + self.config.process_poll_seconds
                        pids = self._cake_processes()
                        if pids:
                            stop_reason = 'ROGUE_CAKEPB'
                            checker_log.write(('ROGUE cake_pb: ' + ' '.join(map(str, pids)) + '\n').encode())
                            self._terminate_all()
                            break
                    if deadline is not None and now >= deadline:
                        stop_reason = 'TIMEOUT'
                        self._terminate_all()
                        break
                    wake = min(now + 0.25, next_probe)
                    if deadline is not None:
                        wake = min(wake, deadline)
                    time.sleep(max(0.001, wake - now))
            except RunInterrupted as error:
                stop_reason = f'SIGNAL_{error.signum}'
                self._terminate_all()
            except (OSError, subprocess.SubprocessError) as error:
                launch_error = f'{type(error).__name__}: {error}'
                checker_log.write(f'LAUNCH ERROR: {error}\n'.encode('utf-8', 'replace'))
                self._terminate_all()
            exits = []
            for process in (xz_process, checker_process):
                if process is None:
                    exits.append(None)
                else:
                    if process.poll() is None:
                        self._terminate_all()
                    exits.append(process.wait())
                    self._complete(process)
            if self.stop_signal is not None and stop_reason is None:
                stop_reason = f'SIGNAL_{self.stop_signal}'
            for stream in (xz_log, checker_log):
                stream.flush()
                os.fsync(stream.fileno())
        finished_ns = time.time_ns()
        return {'checker_argv': list(checker_argv), 'checker_exit_code': exits[1], 'checker_log': self._artifact(checker_log_path), 'duration_seconds': max(0.0, time.monotonic() - started_mono), 'finished_unix_ns': finished_ns, 'finished_utc': utc(finished_ns), 'launch_error': launch_error, 'started_unix_ns': started_ns, 'started_utc': utc(started_ns), 'stop_reason': stop_reason, 'xz_argv': list(xz_argv), 'xz_exit_code': exits[0], 'xz_log': self._artifact(xz_log_path)}

    def _cake_processes(self) -> list[int]:
        self._verify_tool('pgrep')
        code, stdout, stderr = self._capture([self.tools['pgrep']['path'], '-x', 'cake_pb'], env=self.child_env)
        if code == 1 and (not stdout.strip()):
            return []
        if code != 0:
            detail = stderr.decode('utf-8', 'replace').strip()
            raise Refused(f'cake_pb probe failed ({code}): {detail}')
        fields = stdout.split()
        if not fields or any((not field.isdigit() for field in fields)):
            raise Refused('cake_pb probe returned malformed PIDs')
        return [int(field) for field in fields]

    def _group_rss(self, pgid: int) -> int:
        self._verify_tool('ps')
        code, stdout, stderr = self._capture([self.tools['ps']['path'], '-axo', 'pid=,pgid=,rss='], env=self.child_env)
        if code:
            detail = stderr.decode('utf-8', 'replace').strip()
            raise Refused(f'RSS probe failed ({code}): {detail}')
        total = 0
        for line in stdout.splitlines():
            fields = line.split()
            if len(fields) != 3 or any((not field.isdigit() for field in fields)):
                raise Refused('RSS probe returned malformed process data')
            if int(fields[1]) == pgid:
                total += int(fields[2]) * 1024
        return total

    def _wait_lock(self, lock: FileLock, deadline: float | None) -> None:
        while True:
            self._check_signal()
            try:
                lock.acquire()
                return
            except LockBusy:
                if deadline is not None and time.monotonic() >= deadline:
                    raise ResourceWaitTimeout('timed out waiting for host-heavy lock')
                time.sleep(self.config.process_poll_seconds)

    @contextlib.contextmanager
    def _host_exclusion(self) -> Iterator[None]:
        deadline = None if self.config.resource_wait_seconds == 0 else time.monotonic() + self.config.resource_wait_seconds
        lock = FileLock(self.config.host_lock_path)
        self._wait_lock(lock, deadline)
        self.host_lock = lock
        try:
            while self._cake_processes():
                if deadline is not None and time.monotonic() >= deadline:
                    raise ResourceWaitTimeout(
                        'timed out waiting for cake_pb under host lock')
                time.sleep(self.config.process_poll_seconds)
            yield
        except BaseException:
            # The inherited host-lock descriptor must remain locked until
            # every verifier child has terminated and its leader is reaped.
            self._terminate_all()
            raise
        finally:
            self.host_lock = None
            lock.release()

    def _verify_proof(self, case: str, proof: Path, paths: dict[str, Path]):
        self._verify_cnf(case)
        self._verify_tools()
        before = file_record(proof, self.workspace)
        cnf = self.bindings[case]['resolved_path']
        xz = self.tools['xz']['path']
        checker = self.tools['checker']['path']
        checker_argv = [checker, cnf, '-i', '-w']
        with self._host_exclusion():
            self._verify_cnf(case)
            self._verify_tools()
            integrity = self._phase([xz, '--format=xz', '--test', '--', str(proof)], paths['xz_integrity'], timeout=self.config.checker_timeout_seconds, cake_guard=True)
            if integrity['stop_reason'] == 'ROGUE_CAKEPB':
                return (integrity, None, None, None, 'ROGUE_CAKEPB_ABORT')
            if integrity['exit_code'] != 0:
                return (integrity, None, None, None, 'XZ_INTEGRITY_FAILED')
            input_up = self._phase(checker_argv, paths['input_up'], timeout=self.config.checker_timeout_seconds, cake_guard=True, stdin=subprocess.DEVNULL)
            if input_up['stop_reason'] == 'ROGUE_CAKEPB':
                return (integrity, input_up, None, None, 'ROGUE_CAKEPB_ABORT')
            proof_check = self._pipeline([xz, '--format=xz', '--decompress', '--stdout', '--', str(proof)], checker_argv, paths['xz_stream'], paths['checker'])
            verdict, outcome = self._classify(input_up, proof_check)
            if proof_check['stop_reason'] == 'ROGUE_CAKEPB' or self._cake_processes():
                proof_check['stop_reason'] = 'ROGUE_CAKEPB'
                outcome = 'ROGUE_CAKEPB_ABORT'
            self._verify_cnf(case)
            self._verify_tools()
            if file_record(proof, self.workspace) != before:
                raise Refused('proof changed during verification')
            return (integrity, input_up, proof_check, verdict, outcome)

    def _write_start(self, case: str, attempt: str, started_ns: int, solver_argv: list[str], paths: dict[str, Path], options: dict) -> dict:
        document = {'artifact_targets': {key: str(value) for key, value in paths.items()}, 'attempt_id': attempt, 'case': case, 'certificate_manifest_sha256': self.certificate_sha, 'child_environment': self.child_env, 'child_environment_sha256': self.child_env_sha, 'cnf': self.bindings[case], 'decision_manifest_sha256': self.decision_sha, 'options': options, 'schema': 'kobon-n12-unsat-certificate-attempt-start/1', 'solver_argv': solver_argv, 'started_unix_ns': started_ns, 'started_utc': utc(started_ns), 'tools': self.tools}
        path = self._start_path(case, attempt)
        atomic_write(path, canonical_json(document))
        return self._artifact(path)

    def _execute(self, case: str) -> dict:
        self._verify_cnf(case)
        self._verify_tools()
        assert self.proof_dir is not None
        if self.disk_free(self.proof_dir) < self.config.disk_floor_bytes:
            raise Refused(f'disk below floor before solver launch: {case}')
        attempt = time.strftime('%Y%m%dT%H%M%SZ', time.gmtime()) + f'.p{os.getpid()}.{secrets.token_hex(8)}'
        paths = self._paths(case, attempt)
        if any((path.exists() for path in paths.values())):
            raise Refused(f'attempt collision: {attempt}')
        started_ns, started_mono = (time.time_ns(), time.monotonic())
        options = self._options()
        solver_argv = [self.tools['kissat']['path'], self.bindings[case]['resolved_path'], str(paths['proof'])]
        start_record = self._write_start(case, attempt, started_ns, solver_argv, paths, options)
        solver = self._phase(solver_argv, paths['solver'], timeout=self.config.solver_timeout_seconds, disk_guard=True, memory_guard=self.config.memory_limit_bytes is not None, preexec_fn=_memory_limit(self.config.memory_limit_bytes))
        solver_data = self._artifact_path(solver['log']).read_bytes()
        solver_result = solver_verdict(solver_data)
        proof = self._artifact(paths['proof']) if paths['proof'].is_file() else None
        integrity = input_up = proof_check = None
        checker_result = None
        if solver['stop_reason'] == 'DISK_FLOOR':
            outcome = 'DISK_FLOOR_ABORT'
        elif solver['stop_reason'] == 'MEMORY_LIMIT':
            outcome = 'MEMORY_LIMIT_ABORT'
        elif solver['stop_reason'] == 'TIMEOUT':
            outcome = 'SOLVER_TIMEOUT'
        elif solver['stop_reason'] and solver['stop_reason'].startswith('SIGNAL_'):
            outcome = 'INTERRUPTED'
        elif solver['launch_error'] is not None:
            outcome = 'SOLVER_LAUNCH_FAILED'
        elif solver['exit_code'] == 10:
            outcome = 'SAT_NOT_CERTIFIED' if solver_result == 'SAT' else 'SOLVER_EVIDENCE_MISMATCH'
        elif solver['exit_code'] is not None and solver['exit_code'] < 0:
            outcome = 'SOLVER_RESOURCE_LIMIT'
        elif self.config.memory_limit_bytes is not None and _MEMORY_ERROR_RE.search(solver_data):
            outcome = 'SOLVER_RESOURCE_LIMIT'
        elif solver['exit_code'] != 20:
            outcome = 'SOLVER_FAILED'
        elif solver_result != 'UNSAT':
            outcome = 'SOLVER_EVIDENCE_MISMATCH'
        elif proof is None or proof['bytes'] == 0:
            outcome = 'PROOF_MISSING_OR_EMPTY'
        else:
            try:
                integrity, input_up, proof_check, checker_result, outcome = self._verify_proof(case, paths['proof'], paths)
            except ResourceWaitTimeout:
                outcome = 'RESOURCE_WAIT_TIMEOUT'
            except RunInterrupted:
                outcome = 'INTERRUPTED'
        if self.stop_signal is not None:
            outcome = 'INTERRUPTED'
        finished_ns = time.time_ns()
        row = {
            'attempt_id': attempt,
            'case': case,
            'certificate_manifest_sha256': self.certificate_sha,
            'checker_verdict': checker_result,
            'child_environment': self.child_env,
            'child_environment_sha256': self.child_env_sha,
            'cnf': self.bindings[case],
            'decision_manifest_sha256': self.decision_sha,
            'duration_seconds': max(0.0, time.monotonic() - started_mono),
            'evidence_class': EVIDENCE.get(outcome, 'NOT_CERTIFIED'),
            'finished_unix_ns': finished_ns,
            'finished_utc': utc(finished_ns),
            'input_up': input_up,
            'interruption_signal': (
                self.stop_signal if outcome == 'INTERRUPTED' else None
            ),
            'options': options,
            'orphan_artifacts': [],
            'outcome': outcome,
            'proof': proof,
            'proof_check': proof_check,
            'schema': ROW_SCHEMA,
            'solver': solver,
            'solver_verdict': solver_result,
            'start_record': start_record,
            'started_unix_ns': started_ns,
            'started_utc': utc(started_ns),
            'terminal': outcome in SUCCESS,
            'tools': self.tools,
            'xz_integrity': integrity,
        }
        self._commit(row)
        if self.stop_signal is not None:
            raise RunInterrupted(self.stop_signal)
        return row

    def _terminal_rows(self) -> dict[str, dict]:
        assert self.ledger is not None
        return {row['case']: row for row in self.ledger['rows'] if row['terminal']}

    def _attempted_cases(self) -> dict[str, list[dict]]:
        assert self.ledger is not None
        result = {case: [] for case in self.config.expected_cases}
        for row in self.ledger['rows']:
            result[row['case']].append(row)
        return result

    def _revalidate_terminal(self, row: dict) -> None:
        case, attempt = (row['case'], row['attempt_id'])
        assert self.log_dir is not None and self.recheck_dir is not None
        check_id = f'{case}.{attempt}.resume-{secrets.token_hex(8)}'
        paths = {'xz_integrity': self.log_dir / f'{check_id}.xz-integrity.log', 'input_up': self.log_dir / f'{check_id}.input-up.log', 'xz_stream': self.log_dir / f'{check_id}.xz-stream.log', 'checker': self.log_dir / f'{check_id}.checker.log'}
        started = time.time_ns()
        result = error = None
        try:
            result = self._verify_proof(case, self._artifact_path(row['proof']), paths)
            self._check_signal()
            integrity, input_up, proof_check, verdict, outcome = result
            if integrity is None or integrity['exit_code'] != 0 or input_up is None or (proof_check is None) or (verdict != 'VERIFIED') or (outcome != row['outcome']):
                raise Refused(f'resume recheck changed evidence: {case}')
        except BaseException as caught:
            error = f'{type(caught).__name__}: {caught}'
            raise
        finally:
            finished = time.time_ns()
            artifacts = {key: self._artifact(path) for key, path in paths.items() if path.is_file()}
            document = {'artifacts': artifacts, 'case': case, 'certificate_manifest_sha256': self.certificate_sha, 'decision_manifest_sha256': self.decision_sha, 'duration_seconds': max(0.0, (finished - started) / 1000000000.0), 'error': error, 'finished_unix_ns': finished, 'finished_utc': utc(finished), 'original_attempt_id': attempt, 'result_outcome': None if result is None else result[4], 'schema': RECHECK_SCHEMA, 'started_unix_ns': started, 'started_utc': utc(started)}
            atomic_write(self.recheck_dir / f'{check_id}.json', canonical_json(document))

    def _status(
        self, cases: Sequence[str], freshly_reverified: set[str]
    ) -> dict:
        attempted, terminals = self._attempted_cases(), self._terminal_rows()
        return {
            'cases': [
                {
                    'attempts': len(attempted[case]),
                    'case': case,
                    'status': (
                        'FRESHLY_REVERIFIED_TERMINAL'
                        if case in freshly_reverified
                        else (
                            'NONTERMINAL_RETRY_REQUIRED'
                            if attempted[case] else 'UNATTEMPTED'
                        )
                    ),
                }
                for case in cases
            ],
            'certificate_manifest_sha256': self.certificate_sha,
            'decision_manifest_sha256': self.decision_sha,
            'ledger': str(self.ledger_path),
            'schema': STATUS_SCHEMA,
            'validation': (
                'MANIFEST_LEDGER_TOOLS_ARTIFACTS_AND_SELECTED_TERMINAL_'
                'PROOFS_FRESHLY_REVALIDATED'
            ),
        }

    def run(self, cases: Sequence[str], *, all_cases: bool=False, status: bool=False, retry_nonterminal: bool=False) -> int:
        self.driver_lock.acquire()
        self._install_signals()
        try:
            self.prepare()
            requested = list(cases)
            if all_cases and requested:
                raise Refused('--all and named cases are mutually exclusive')
            if all_cases or (status and (not requested)):
                requested = list(self.config.expected_cases)
            if not requested:
                raise Refused('name manifest cases or use --all')
            if len(set(requested)) != len(requested):
                raise Refused('duplicate case names')
            if any((not isinstance(case, str) or case not in self.bindings for case in requested)):
                raise Refused('case not in manifest case_order')
            order = {case: index for index, case in enumerate(self.config.expected_cases)}
            requested.sort(key=order.__getitem__)
            if status:
                terminals = self._terminal_rows()
                reverified = set()
                for case in requested:
                    if case in terminals:
                        self._revalidate_terminal(terminals[case])
                        reverified.add(case)
                print(
                    canonical_json(
                        self._status(requested, reverified)
                    ).decode().rstrip(),
                    file=self.out,
                )
                return EXIT_OK
            terminals, attempted = (self._terminal_rows(), self._attempted_cases())
            incomplete = False
            for case in requested:
                if case in terminals:
                    self._revalidate_terminal(terminals[case])
                    event = {'case': case, 'event': 'RESUMED_EXACT_TERMINAL', 'outcome': terminals[case]['outcome']}
                elif attempted[case] and (not retry_nonterminal):
                    if not all_cases:
                        raise Refused(f'{case} has a preserved nonterminal attempt; pass --retry-nonterminal')
                    incomplete = True
                    event = {'case': case, 'event': 'NONTERMINAL_RETRY_REQUIRED', 'outcome': attempted[case][-1]['outcome']}
                    print(json.dumps(event, sort_keys=True), file=self.out)
                    if any(
                        row['outcome'] == 'DISK_FLOOR_ABORT'
                        for row in attempted[case]
                    ):
                        break
                    continue
                else:
                    row = self._execute(case)
                    incomplete |= not row['terminal']
                    event = {'attempt_id': row['attempt_id'], 'case': case, 'event': 'ATTEMPT_FINISHED', 'outcome': row['outcome']}
                    if row['outcome'] == 'DISK_FLOOR_ABORT':
                        print(json.dumps(event, sort_keys=True), file=self.out)
                        break
                print(json.dumps(event, sort_keys=True), file=self.out)
            return EXIT_INCOMPLETE if incomplete else EXIT_OK
        finally:
            self._terminate_all()
            self._restore_signals()
            self.driver_lock.release()

def _gib(value: float, name: str) -> int | None:
    if not math.isfinite(value) or value < 0:
        raise Refused(f'{name} must be finite and nonnegative')
    return None if value == 0 else int(value * (1 << 30))

def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument('cases', nargs='*', help='exact manifest case names')
    result.add_argument('--all', action='store_true', help='certify all cases')
    result.add_argument('--status', action='store_true')
    result.add_argument('--retry-nonterminal', action='store_true')
    result.add_argument('--disk-floor-gib', type=float, default=DEFAULT_DISK_FLOOR_GIB)
    result.add_argument('--memory-limit-gib', type=float, default=DEFAULT_MEMORY_LIMIT_GIB, help='solver RSS/address-space ceiling; 0 means unlimited')
    result.add_argument('--disk-poll-seconds', type=float, default=DEFAULT_DISK_POLL_SECONDS)
    result.add_argument('--process-poll-seconds', type=float, default=DEFAULT_PROCESS_POLL_SECONDS)
    result.add_argument('--solver-timeout-seconds', type=float, default=0.0)
    result.add_argument('--checker-timeout-seconds', type=float, default=0.0)
    result.add_argument('--resource-wait-seconds', type=float, default=0.0, help='0 waits indefinitely for the shared host lock')
    result.add_argument('--kill-grace-seconds', type=float, default=DEFAULT_KILL_GRACE_SECONDS)
    return result

def main(argv: Sequence[str] | None=None) -> int:
    args = parser().parse_args(argv)
    try:
        base = default_config()
        config = Config(**{**base.__dict__, 'disk_floor_bytes': _gib(args.disk_floor_gib, 'disk floor') or 0, 'memory_limit_bytes': _gib(args.memory_limit_gib, 'memory limit'), 'disk_poll_seconds': args.disk_poll_seconds, 'process_poll_seconds': args.process_poll_seconds, 'solver_timeout_seconds': args.solver_timeout_seconds, 'checker_timeout_seconds': args.checker_timeout_seconds, 'resource_wait_seconds': args.resource_wait_seconds, 'kill_grace_seconds': args.kill_grace_seconds})
        return Certifier(config).run(args.cases, all_cases=args.all, status=args.status, retry_nonterminal=args.retry_nonterminal)
    except RunInterrupted as error:
        print(f'INTERRUPTED: {error}', file=sys.stderr)
        return 128 + error.signum
    except (Refused, OSError, ValueError) as error:
        print(f'REFUSED: {error}', file=sys.stderr)
        return EXIT_REFUSED
if __name__ == '__main__':
    sys.exit(main())
