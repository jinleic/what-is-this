#!/usr/bin/env python3
"""Guarded collector tests and external, serial checkpoint-driven agent resumption."""

from __future__ import annotations

import argparse
import ast
import ctypes
import fcntl
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent
CHECKPOINT = HERE / "collector-integrity-checkpoint.json"
SOURCE_FILES = (
    "focused-deribit-gex/capture.py",
    "focused-liquidation-overlay/capture.py",
    "focused-deribit-gex/proofs_capture.py",
    "test_collector_integrity.py",
)
OUTPUT_LIMIT = 256 * 1024


def utc_now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def load_checkpoint():
    return json.loads(CHECKPOINT.read_text())


def save_checkpoint(checkpoint):
    checkpoint["updated_utc"] = utc_now()
    root = Path(checkpoint["constraints"]["artifact_root"])
    measured = allocated_bytes(checkpoint)
    checkpoint["constraints"]["artifact_allocated_bytes_last_measured"] = measured
    payload = json.dumps(checkpoint, indent=2, sort_keys=True) + "\n"
    staged_bytes = (len(payload) + 4095) // 4096 * 4096
    if measured + staged_bytes > checkpoint["constraints"]["artifact_limit_bytes"]:
        raise RuntimeError("insufficient artifact budget for an atomic checkpoint update")
    fd, staged = tempfile.mkstemp(prefix="collector-checkpoint-", suffix=".json", dir=root)
    with os.fdopen(fd, "w") as stream:
        stream.write(payload)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(staged, CHECKPOINT)
    directory_fd = os.open(CHECKPOINT.parent, os.O_RDONLY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)


def source_hashes():
    return {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest() for name in SOURCE_FILES}


def test_inventory():
    tree = ast.parse((HERE / "test_collector_integrity.py").read_text())
    return sorted(
        f"test_collector_integrity.{node.name}.{method.name}"
        for node in tree.body if isinstance(node, ast.ClassDef)
        for method in node.body
        if isinstance(method, ast.FunctionDef) and method.name.startswith("test_")
    )


def allocated_bytes(checkpoint):
    root = Path(checkpoint["constraints"]["artifact_root"])
    total = int(subprocess.check_output(["du", "-sk", str(root)], text=True).split()[0]) * 1024
    paths = {CHECKPOINT, Path(__file__).resolve()}
    paths.update(Path(path) for path in checkpoint["constraints"].get("additional_artifact_files", []))
    for path in paths:
        if path.exists() and not path.is_relative_to(root):
            total += path.stat().st_blocks * 512
    return total + checkpoint["constraints"]["retained_tempfile_probe_reserve_bytes"]


def resource_blocker(checkpoint):
    reserve = OUTPUT_LIMIT + CHECKPOINT.stat().st_size
    if allocated_bytes(checkpoint) >= checkpoint["constraints"]["artifact_limit_bytes"] - reserve:
        return "artifact budget reserve reached; deletion and budget reset are forbidden"
    if shutil.disk_usage(checkpoint["constraints"]["artifact_root"]).free < 100 * 1024**3:
        return "less than 100 GiB free disk"
    return None


class CpuSampler:
    def __init__(self):
        self.library = ctypes.CDLL("/usr/lib/libSystem.B.dylib")
        self.library.mach_host_self.restype = ctypes.c_uint
        self.library.host_statistics.argtypes = [
            ctypes.c_uint, ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_uint)
        ]
        self.host = self.library.mach_host_self()
        self.previous = self.ticks()

    def ticks(self):
        ticks = (ctypes.c_int * 4)()
        count = ctypes.c_uint(4)
        result = self.library.host_statistics(self.host, 3, ticks, ctypes.byref(count))
        if result != 0 or count.value != 4:
            raise RuntimeError("whole-machine CPU measurement failed")
        return tuple(ctypes.c_uint(value).value for value in ticks)

    def sample(self, interval):
        time.sleep(interval)
        current = self.ticks()
        delta = [(after - before) & 0xFFFFFFFF for before, after in zip(self.previous, current)]
        self.previous = current
        if not sum(delta):
            return None
        return 100 * (sum(delta) - delta[2]) / sum(delta)


def enforcement_mode(checkpoint):
    """Bind coverage to the runner and successor boundary, not a mutable label."""
    digest = hashlib.sha256(Path(__file__).read_bytes())
    digest.update((HERE / "collector_integrity_extension.ts").read_bytes())
    return "retained-offline-v1:" + digest.hexdigest()


def launch_gate():
    sampler = CpuSampler()
    samples = []
    consecutive = 0
    for _ in range(120):
        value = sampler.sample(1)
        samples.append(round(value, 2) if value is not None else None)
        consecutive = consecutive + 1 if value is not None and value < 30 else 0
        if consecutive == 3:
            return {"passed": True, "sample_interval_seconds": 1, "samples_taken": len(samples), "cpu_percent": samples[-12:]}
    return {"passed": False, "sample_interval_seconds": 1, "samples_taken": len(samples), "cpu_percent": samples[-12:]}


def signal_owned_group(process, sig):
    if process.poll() is None and os.getpgid(process.pid) == process.pid:
        os.killpg(process.pid, sig)


def run_process(command, checkpoint, wall_limit, worker_limits=False, launch_record=None, extra_env=None, cwd=HERE):
    environment = os.environ.copy()
    root = Path(checkpoint["constraints"]["artifact_root"])
    write_root = root
    if worker_limits:
        write_root = Path(tempfile.mkdtemp(prefix="collector-worker-", dir=root))
        write_root = write_root.resolve()
        (write_root / "sitecustomize.py").write_text(WORKER_POLICY)
        runtime_paths = [path for path in sys.path if path.startswith("/opt/homebrew/")]
        environment["PYTHONPATH"] = os.pathsep.join((str(write_root), str(HERE), *runtime_paths))
        environment["COLLECTOR_SCRATCH"] = str(write_root)
        command = ["/usr/bin/sandbox-exec", "-p", sandbox_profile(checkpoint, write_root), *command]
    environment.update({"TMPDIR": str(write_root), "PYTHONDONTWRITEBYTECODE": "1"})
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        environment[name] = "1"
    if extra_env:
        environment.update(extra_env)

    def configure_child():
        os.setpriority(os.PRIO_PROCESS, 0, 15)
        if worker_limits:
            resource.setrlimit(resource.RLIMIT_FSIZE, (1048576, 1048576))
            resource.setrlimit(resource.RLIMIT_CPU, (10, 10))

    sampler = CpuSampler()
    started = time.monotonic()
    process = subprocess.Popen(
        command, cwd=cwd, env=environment, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        stdin=subprocess.DEVNULL, start_new_session=True, preexec_fn=configure_child,
    )
    output = bytearray()
    peak = None
    paused = False
    pauses = 0
    unavailable_sample_pauses = 0
    stop_reason = None

    def drain():
        while True:
            try:
                chunk = os.read(process.stdout.fileno(), 65536)
            except BlockingIOError:
                break
            if not chunk:
                break
            output.extend(chunk)
            if len(output) > OUTPUT_LIMIT:
                del output[:-OUTPUT_LIMIT]

    try:
        os.set_blocking(process.stdout.fileno(), False)
        if launch_record is not None:
            launch_record.update({"worker_launched": True, "pid": process.pid, "outcome": "running"})
            save_checkpoint(checkpoint)
        while process.poll() is None:
            drain()
            value = sampler.sample(0.1)
            if value is not None:
                peak = value if peak is None else max(peak, value)
            if (value is None or value >= 45) and not paused and process.poll() is None:
                signal_owned_group(process, signal.SIGSTOP)
                paused = True
                pauses += 1
                unavailable_sample_pauses += value is None
            elif value is not None and value < 30 and paused:
                signal_owned_group(process, signal.SIGCONT)
                paused = False
            stop_reason = resource_blocker(checkpoint)
            if process.poll() is not None:
                break
            if time.monotonic() - started >= wall_limit:
                stop_reason = "owned worker wall deadline"
            if stop_reason:
                signal_owned_group(process, signal.SIGTERM)
                if paused:
                    signal_owned_group(process, signal.SIGCONT)
                process.wait(timeout=5)
                break
        drain()
    finally:
        if process.poll() is None:
            signal_owned_group(process, signal.SIGTERM)
            signal_owned_group(process, signal.SIGCONT)
            process.wait(timeout=5)
        process.stdout.close()
    return {
        "pid": process.pid, "exit_code": process.returncode, "nice": 15,
        "wall_seconds": round(time.monotonic() - started, 3),
        "sampled_peak_system_cpu_percent": round(peak, 2) if peak is not None else None,
        "cpu_guard_pauses": pauses,
        "unavailable_sample_pauses": unavailable_sample_pauses,
        "stop_reason": stop_reason, "artifact_bytes_total": allocated_bytes(checkpoint),
    }, output.decode("utf-8", errors="replace")


def sandbox_profile(checkpoint, write_root):
    """Kernel containment of worker reads, writes, network and process signals."""
    readable = ("/System", "/usr", "/opt/homebrew", "/private/etc", "/dev",
                str(HERE), str(Path(checkpoint["constraints"]["artifact_root"]).resolve()))
    executables = (Path(sys.executable).resolve(), (Path(sys.exec_prefix) / "Resources/Python.app/Contents/MacOS/Python").resolve())
    exec_paths = " ".join(f"(literal {json.dumps(str(path))})" for path in executables)
    paths = " ".join(f"(subpath {json.dumps(path)})" for path in readable)
    return "\n".join((
        "(version 1)", "(deny default)",
        '(import "dyld-support.sb")',
        "(allow file-read-metadata)",
        f"(allow file-read-data file-map-executable {paths})",
        f"(allow file-write* (subpath {json.dumps(str(write_root))}))",
        "(deny file-write-data file-write-unlink file-write-mode"
        f" (subpath {json.dumps(str(write_root / '.retained'))}))",
        "(allow process-fork sysctl-read)",
        f"(allow process-exec {exec_paths})",
        '(allow file-write-data (literal "/dev/null"))',
        "(allow signal (target same-sandbox))",
        '(allow network-bind network-inbound (local ip "localhost:*"))',
        '(allow network-outbound (remote ip "localhost:*"))',
    ))


WORKER_POLICY = r'''
import errno, os, stat, sys, tempfile
tempfile.tempdir = os.environ['COLLECTOR_SCRATCH']
scratch = os.path.realpath(tempfile.tempdir)
violations = []
current_test = None
retained_count = 0
root_fd = os.open(scratch, os.O_RDONLY | os.O_DIRECTORY)
root_stat = os.fstat(root_fd)
root_identity = (root_stat.st_dev, root_stat.st_ino)
retained = os.path.join(scratch, '.retained')
os.makedirs(retained, exist_ok=True)

def checked_parent(path, dir_fd=None):
    parent, name = os.path.split(os.fsdecode(path))
    if not name or name in ('.', '..'):
        raise PermissionError('cannot retire a directory root')
    fd = os.open(parent or '.', os.O_RDONLY | os.O_DIRECTORY, dir_fd=dir_fd)
    cursor = os.dup(fd)
    try:
        while True:
            metadata = os.fstat(cursor)
            identity = (metadata.st_dev, metadata.st_ino)
            if identity == root_identity:
                return fd, name
            ancestor = os.open('..', os.O_RDONLY | os.O_DIRECTORY, dir_fd=cursor)
            above = os.fstat(ancestor)
            os.close(cursor)
            cursor = ancestor
            if identity == (above.st_dev, above.st_ino):
                raise PermissionError('retirement outside worker scratch')
    except BaseException:
        os.close(fd)
        raise
    finally:
        os.close(cursor)

def retire(path, *, dir_fd=None, directory=False):
    global retained_count
    parent, name = checked_parent(path, dir_fd)
    try:
        metadata = os.stat(name, dir_fd=parent, follow_symlinks=False)
        is_directory = stat.S_ISDIR(metadata.st_mode)
        if directory:
            if not is_directory:
                raise NotADirectoryError(errno.ENOTDIR, 'not a directory', path)
            child = os.open(name, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=parent)
            try:
                if os.listdir(child):
                    raise OSError(errno.ENOTEMPTY, 'directory not empty', path)
            finally:
                os.close(child)
        elif is_directory:
            raise IsADirectoryError(errno.EISDIR, 'is a directory', path)
        slot = tempfile.mkdtemp(prefix='entry-', dir=retained)
        destination = os.open(slot, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.rename(name, 'value', src_dir_fd=parent, dst_dir_fd=destination)
            retained_count += 1
        finally:
            os.close(destination)
    finally:
        os.close(parent)

def retire_directory(path, *, dir_fd=None):
    return retire(path, dir_fd=dir_fd, directory=True)

os.unlink = os.remove = retire
os.rmdir = retire_directory
os.supports_dir_fd.update((retire, retire_directory))
loopback = {'localhost', '127.0.0.1', '::1'}
def guard(event, arguments):
    forbidden = event in {'os.remove', 'os.rmdir', 'os.system', 'os.posix_spawn', 'os.exec', 'ctypes.dlopen'}
    if event == 'subprocess.Popen':
        executable, argv, cwd, environment = arguments
        forbidden = (
            os.path.realpath(executable) != os.path.realpath(sys.executable)
            or any(flag in argv for flag in ('-S', '-I', '-E'))
            or (environment is not None and (
                environment.get('COLLECTOR_SCRATCH') != scratch
                or environment.get('PYTHONPATH') != os.environ['PYTHONPATH']
            ))
        )
    if event == 'socket.getaddrinfo':
        forbidden = arguments[0] not in loopback
    if event == 'socket.connect':
        address = arguments[1]
        forbidden = not isinstance(address, tuple) or address[0] not in loopback
    if forbidden:
        violations.append({'test': current_test, 'operation': event})
        raise PermissionError('retained offline worker policy: ' + event)
sys.addaudithook(guard)
'''


WORKER = r'''
import os, sys, json, unittest, sitecustomize as policy
assert policy.scratch == os.path.realpath(os.environ['TMPDIR'])
class Result(unittest.TextTestResult):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.passed_ids = []
    def startTest(self, test):
        policy.current_test = test.id()
        super().startTest(test)
    def addSuccess(self, test):
        self.passed_ids.append(test.id())
        super().addSuccess(test)
result = unittest.TextTestRunner(verbosity=1, resultclass=Result).run(
    unittest.defaultTestLoader.loadTestsFromNames(sys.argv[1:])
)
blocked = sorted({item['test'] for item in policy.violations if item['test'] is not None})
print('COLLECTOR_BATCH_RESULT ' + json.dumps({
    'ran': result.testsRun, 'passed': [name for name in result.passed_ids if name not in blocked],
    'failed': [test.id() for test, _ in result.failures],
    'errors': [test.id() for test, _ in result.errors],
    'skipped': [test.id() for test, _ in result.skipped],
    'policy_violations': policy.violations, 'blocked': blocked,
    'retained_cleanup_entries': policy.retained_count,
}))
raise SystemExit(not result.wasSuccessful() or bool(blocked))
'''


def update_coverage(checkpoint):
    inventory = test_inventory()
    hashes = source_hashes()
    latest = {}
    for attempt in checkpoint["attempts"]:
        if not attempt["counts_toward_full_suite"] or attempt.get("source_sha256") != hashes:
            continue
        if attempt.get("enforcement_mode") != enforcement_mode(checkpoint):
            continue
        for status in ("passed", "failed", "errors", "skipped", "blocked"):
            for test_id in attempt["test_result"].get(status, []):
                latest[test_id] = status
    passed = [name for name in inventory if latest.get(name) == "passed"]
    checkpoint["source_sha256"] = hashes
    checkpoint["verification"].update({
        "inventory_test_ids": inventory, "total_tests": len(inventory),
        "passed_test_ids": passed, "passed_count": len(passed),
        "pending_test_ids": [name for name in inventory if latest.get(name) != "passed"],
        "blocked_tests": sorted(name for name, status in latest.items() if status == "blocked"),
        "failed_test_ids": sorted(name for name, status in latest.items() if status in ("failed", "errors")),
        "full_suite": "passed" if len(passed) == len(inventory) else "incomplete",
    })
    return latest


def run_tests(names):
    checkpoint = load_checkpoint()
    root = Path(checkpoint["constraints"]["artifact_root"])
    with (root / "collector-execution.lock").open("a") as execution, (root / "collector-worker.lock").open("a") as lease:
        fcntl.flock(execution, fcntl.LOCK_EX | fcntl.LOCK_NB)
        fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
        latest = update_coverage(checkpoint)
        if not names:
            names = [name for name in checkpoint["verification"]["pending_test_ids"] if name not in latest][:checkpoint["constraints"]["batch_size"]]
        if not names:
            print(json.dumps({"outcome": "no_unattempted_tests", "verification": checkpoint["verification"]["full_suite"]}))
            return 0
        if any(name not in checkpoint["verification"]["inventory_test_ids"] for name in names):
            raise ValueError("requested test is not in the collector test inventory")
        attempt = {
            "attempt": len(checkpoint["attempts"]) + 1, "label": "guarded collector tests",
            "started_utc": utc_now(), "requested_test_ids": names,
            "runner_revision": "retained-offline-default-deny",
            "enforcement_mode": enforcement_mode(checkpoint),
            "source_sha256": checkpoint["source_sha256"], "worker_launched": False,
            "resources": None, "test_result": None, "exit_code": None,
            "counts_toward_full_suite": False,
        }
        attempt["outcome"] = "checking_launch_gate"
        checkpoint["attempts"].append(attempt)
        checkpoint["latest_attempt"] = attempt["attempt"]
        save_checkpoint(checkpoint)
        blocker = resource_blocker(checkpoint)
        attempt["gate"] = launch_gate() if not blocker else None
        if blocker or not attempt["gate"]["passed"]:
            attempt["outcome"] = "resource_blocked" if blocker else "deferred_without_launch"
            attempt["blocker"] = blocker or "three-consecutive-sample CPU launch gate did not pass"
        else:
            attempt.update({"worker_launched": None, "outcome": "launch_intent"})
            save_checkpoint(checkpoint)
            try:
                attempt["resources"], output = run_process(
                    [sys.executable, "-B", "-c", WORKER, *names], checkpoint, 25, True, attempt
                )
            except Exception as error:
                attempt["outcome"] = "interrupted_without_results"
                attempt["controller_error"] = f"{type(error).__name__}: {error}"
            else:
                attempt["exit_code"] = attempt["resources"]["exit_code"]
                for line in output.splitlines():
                    if line.startswith("COLLECTOR_BATCH_RESULT "):
                        attempt["test_result"] = json.loads(line.split(" ", 1)[1])
                attempt["counts_toward_full_suite"] = (
                    attempt["test_result"] is not None and attempt["resources"]["stop_reason"] is None
                )
                attempt["outcome"] = "passed" if attempt["exit_code"] == 0 else "tests_or_policy_failed"
                if attempt["test_result"] is None:
                    attempt["outcome"] = "interrupted_without_results"
                if attempt["exit_code"] != 0:
                    attempt["diagnostic_output"] = output[:8192]
        attempt["ended_utc"] = utc_now()
        update_coverage(checkpoint)
        checkpoint["status"] = "blocked" if (
            not attempt["worker_launched"] or attempt.get("exit_code") != 0
        ) else "running"
        checkpoint["next_action"]["prompt"] = (
            "Inspect the latest recorded collector result. Run a pending batch only through the guarded "
            "runner; investigate failed or policy-blocked cases without bypassing protections. Preserve "
            "all artifacts. If CPU or disk blocks runtime, inspect one scoped source risk, record an "
            "actionable next prompt, and return so the external controller can back off. "
            "Do not start captures or another supervisor. Do not claim unobserved results."
        )
        save_checkpoint(checkpoint)
        print(json.dumps({"attempt": attempt["attempt"], "outcome": attempt["outcome"],
                          "worker_launched": attempt["worker_launched"], "gate": attempt["gate"],
                          "resources": attempt["resources"], "passed": checkpoint["verification"]["passed_count"],
                          "pending": len(checkpoint["verification"]["pending_test_ids"])}))
        if attempt.get("diagnostic_output"):
            print(attempt["diagnostic_output"])
        if "controller_error" in attempt:
            print(attempt["controller_error"])
            return 70
        return attempt["exit_code"] if attempt["worker_launched"] else 75


def supervise(prove_resumption=False):
    checkpoint = load_checkpoint()
    root = Path(checkpoint["constraints"]["artifact_root"])
    with (root / "collector-supervisor.lock").open("a") as lease:
        fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
        print("COLLECTOR_SUPERVISOR_READY", flush=True)
        successful_handoffs = 0
        while True:
            checkpoint = load_checkpoint()
            supervision = checkpoint["supervisor"]
            if not supervision.get("enabled", False):
                return 0
            if not supervision.get("handoff_ready", False):
                time.sleep(30)
                continue
            if not supervision.get("operator_authorized_launch", False):
                raise PermissionError("successor launch is not authorized")
            if not supervision.get("tool_boundary_verified", False):
                supervision.update({"state": "boundary_blocked", "blocker": "bounded tool preflight has not passed"})
                save_checkpoint(checkpoint)
                return 78
            live_hashes = source_hashes()
            blocker = resource_blocker(checkpoint)
            if live_hashes != checkpoint["source_sha256"]:
                blocker = "source hashes changed outside a recorded controller edit"
            if not prove_resumption and not supervision.get("normal_exit_resumption_observed", False):
                blocker = "safe normal-exit resumption has not been observed"
            gate = launch_gate() if not blocker else None
            if blocker or not gate["passed"]:
                supervision.update({
                    "state": "backing_off", "latest_gate": gate,
                    "blocker": blocker or "three consecutive CPU samples below 30% unavailable",
                    "next_check_utc_epoch": time.time() + 300,
                })
                save_checkpoint(checkpoint)
                print("COLLECTOR_SUPERVISOR_BACKOFF", flush=True)
                time.sleep(300)
                continue
            omp = shutil.which("omp")
            if omp is None:
                raise RuntimeError("omp executable unavailable")
            agent_root = root / "collector-agent"
            agent_root.mkdir(exist_ok=True)
            config = root / "collector-agent-config.json"
            if not config.exists():
                config.write_text(json.dumps({"advisor": {"enabled": False}, "dev": {"autoqa": False}}) + "\n")
            token = os.urandom(16).hex()
            active_agent = {
                "token": token, "outcome": "launch_intent", "worker_launched": False,
                "started_utc": utc_now(), "started_epoch": time.time(),
                "source_sha256": live_hashes, "enforcement_mode": enforcement_mode(checkpoint),
                "input_prompt_sha256": hashlib.sha256(checkpoint["next_action"]["prompt"].encode()).hexdigest(),
                "proof_only": prove_resumption,
            }
            supervision.update({"state": "agent_running", "active_agent": active_agent, "latest_gate": gate})
            save_checkpoint(checkpoint)
            prompt = (
                "Use only collector_control. First call status. Execute the next_action from its result. "
                "Then call note with the precise next prompt before exiting. No captures, deletion, "
                "service restarts, new supervisors, or unrelated work. Only report observed evidence. "
            )
            if prove_resumption:
                prompt += (
                    "This is a read/note-only real successor handoff proof. Do not edit or run tests. "
                    "After status, call note with next_prompt equal to the status next_action followed by "
                    "' | observed successor " + token + "'. Exit immediately after note succeeds."
                )
            command = [
                omp, "-p", "--no-session", "--no-title", "--no-extensions", "--no-pty", "--no-lsp",
                "--no-skills", "--no-rules", "--no-tools", "--auto-approve",
                "-e", str(HERE / "collector_integrity_extension.ts"),
                "--max-time", "2m", "--thinking", "minimal", "--hide-thinking",
                "--model", "openai-codex/gpt-6-astra", "--config", str(config),
                "--system-prompt", "You are a bounded collector successor. Use the supplied tool and obey the task exactly.",
                prompt,
            ]
            with (root / "collector-execution.lock").open("a") as execution:
                fcntl.flock(execution, fcntl.LOCK_EX | fcntl.LOCK_NB)
                result, output = run_process(
                    command, checkpoint, 150, launch_record=active_agent,
                    extra_env={"COLLECTOR_SUCCESSOR_TOKEN": token}, cwd=agent_root,
                )
            checkpoint = load_checkpoint()
            supervision = checkpoint["supervisor"]
            receipt = supervision.get("last_status_receipt", {})
            note = supervision.get("last_note_receipt", {})
            handed_off = (
                result["exit_code"] == 0 and result["stop_reason"] is None
                and receipt.get("token") == token and note.get("token") == token
                and receipt.get("prompt_sha256") == active_agent["input_prompt_sha256"]
            )
            record = {
                **active_agent, "ended_utc": utc_now(), "ended_epoch": time.time(),
                "resources": result, "gate": gate, "handoff_observed": handed_off,
                "output_prompt_sha256": hashlib.sha256(checkpoint["next_action"]["prompt"].encode()).hexdigest(),
            }
            runs = supervision.setdefault("agent_runs", [])
            previous = runs[-1] if runs else None
            runs.append(record)
            resumed = (
                handed_off and previous is not None and previous.get("handoff_observed")
                and previous["ended_epoch"] <= record["started_epoch"]
                and previous["output_prompt_sha256"] == record["input_prompt_sha256"]
                and previous["resources"]["pid"] != result["pid"]
            )
            if resumed:
                supervision["normal_exit_resumption_observed"] = True
                supervision["resumption_evidence"] = {
                    "predecessor_pid": previous["resources"]["pid"], "successor_pid": result["pid"],
                    "predecessor_exit_epoch": previous["ended_epoch"], "successor_start_epoch": record["started_epoch"],
                    "handoff_prompt_sha256": record["input_prompt_sha256"], "verified_utc": utc_now(),
                }
            supervision["active_agent"] = None
            supervision["last_agent_output"] = output[-2048:]
            queued = checkpoint.get("pending_test_request")
            delay = 5 if prove_resumption else (0 if handed_off and queued else 300)
            supervision["state"] = "resuming" if delay < 300 else "backing_off"
            supervision["next_check_utc_epoch"] = time.time() + delay
            supervision["automatic_restart_verified"] = False
            supervision["restart_verification_note"] = "Only observed normal-exit handoffs are certified; no crash or service restart claim."
            save_checkpoint(checkpoint)
            print(json.dumps({"successor_pid": result["pid"], "exit": result["exit_code"],
                              "handoff_observed": handed_off, "resumption_observed": bool(resumed)}), flush=True)
            if handed_off and queued and not prove_resumption:
                if queued["source_sha256"] != source_hashes() or queued["enforcement_mode"] != enforcement_mode(checkpoint):
                    checkpoint["status"] = "blocked"
                    checkpoint["next_action"]["prompt"] = "Queued test provenance changed; inspect current sources before requesting a fresh batch."
                else:
                    checkpoint["pending_test_request"] = None
                    save_checkpoint(checkpoint)
                    run_tests(queued["tests"])
                    checkpoint = load_checkpoint()
                    checkpoint["supervisor"]["last_queued_execution"] = {
                        "successor_pid": result["pid"], "after_exit_epoch": record["ended_epoch"],
                        "attempt": checkpoint["latest_attempt"] if checkpoint["latest_attempt"] != queued.get("prior_attempt") else None,
                        "recorded_utc": utc_now(),
                    }
                checkpoint["supervisor"]["agent_runs"][-1]["output_prompt_sha256"] = hashlib.sha256(
                    checkpoint["next_action"]["prompt"].encode()
                ).hexdigest()
                save_checkpoint(checkpoint)
            successful_handoffs = successful_handoffs + 1 if handed_off else 0
            if prove_resumption and successful_handoffs == 2 and resumed:
                print("COLLECTOR_RESUMPTION_PROVED", flush=True)
                return 0
            if prove_resumption and not handed_off:
                return 1
            if delay:
                time.sleep(delay)



def scoped_path(name, writable=False):
    readable = (*SOURCE_FILES, CHECKPOINT.name, Path(__file__).name, "collector_integrity_extension.ts")
    if name not in (SOURCE_FILES if writable else readable):
        raise PermissionError("path is outside the approved collector scope")
    path = HERE / name
    if path.is_symlink() or path.resolve() != path or not path.is_file() or path.stat().st_nlink != 1:
        raise PermissionError("approved source must be a regular, unaliased file")
    return path


def tool_rpc():
    raw = sys.stdin.read(1024 * 1024 + 1)
    if len(raw) > 1024 * 1024:
        raise ValueError("tool request exceeds 1 MiB")
    request = json.loads(raw)
    operation = request["operation"]
    checkpoint = load_checkpoint()
    active = checkpoint["supervisor"].get("active_agent") or {}
    token = os.environ.get("COLLECTOR_SUCCESSOR_TOKEN")
    if token and token != active.get("token"):
        raise PermissionError("successor lease is no longer current")
    if active.get("proof_only") and operation not in ("status", "read", "search", "note"):
        raise PermissionError("resumption proof permits read and note only")
    if operation == "status":
        update_coverage(checkpoint)
        if token:
            checkpoint["supervisor"]["last_status_receipt"] = {
                "token": token, "recorded_utc": utc_now(),
                "prompt_sha256": hashlib.sha256(checkpoint["next_action"]["prompt"].encode()).hexdigest(),
            }
            save_checkpoint(checkpoint)
        print(json.dumps({
            "source_files": SOURCE_FILES, "source_sha256": checkpoint["source_sha256"],
            "verification": {key: checkpoint["verification"][key] for key in (
                "passed_count", "total_tests", "full_suite", "blocked_tests", "failed_test_ids"
            )},
            "latest_attempt": {key: checkpoint["attempts"][-1].get(key) for key in (
                "attempt", "label", "outcome", "resources", "counts_toward_full_suite"
            )},
            "enforcement_mode": enforcement_mode(checkpoint),
            "next_risk": checkpoint.get("next_risk"), "next_action": checkpoint["next_action"]["prompt"],
            "artifact_bytes": allocated_bytes(checkpoint),
        }))
        return 0
    if operation in ("read", "search"):
        path = scoped_path(request.get("path"))
        text = path.read_text()
        lines = text.splitlines()
        if operation == "read":
            start = max(1, int(request.get("start_line", 1)))
            count = min(200, max(1, int(request.get("line_count", 100))))
            selected = list(enumerate(lines[start - 1:start - 1 + count], start))
        else:
            needle = request.get("needle")
            if not isinstance(needle, str) or not 1 <= len(needle) <= 200:
                raise ValueError("search requires a nonempty literal needle of at most 200 characters")
            selected = [(number, line) for number, line in enumerate(lines, 1) if needle in line][:80]
        print(json.dumps({"path": request["path"], "sha256": hashlib.sha256(text.encode()).hexdigest(),
                          "lines": "\n".join(f"{number}:{line}" for number, line in selected)[:32768]}))
        return 0
    if operation in ("edit", "run"):
        if checkpoint["supervisor"].get("worker_boundary_mode") != enforcement_mode(checkpoint):
            raise PermissionError("current worker enforcement has not passed its boundary smoke")
    if operation == "run":
        names = request.get("tests", [])
        inventory = test_inventory()
        if not isinstance(names, list) or len(names) > checkpoint["constraints"]["batch_size"] or any(name not in inventory for name in names):
            raise ValueError("request requires a bounded list of collector test IDs")
        if token:
            checkpoint["pending_test_request"] = {
                "tests": names, "source_sha256": source_hashes(),
                "enforcement_mode": enforcement_mode(checkpoint), "requested_utc": utc_now(),
                "prior_attempt": checkpoint["latest_attempt"],
            }
            save_checkpoint(checkpoint)
            print(json.dumps({"operation": "run", "queued": True, "worker_launched": False,
                              "execution": "only after this successor exits and persists its next prompt"}))
            return 0
        return run_tests(names)
    if operation not in ("edit", "note"):
        raise ValueError("unknown bounded collector operation")
    root = Path(checkpoint["constraints"]["artifact_root"])
    with (root / "collector-worker.lock").open("a") as lease:
        fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
        checkpoint = load_checkpoint()
        if operation == "note":
            prompt = request.get("next_prompt")
            if not isinstance(prompt, str) or not 1 <= len(prompt.strip()) <= 4000:
                raise ValueError("next prompt must contain 1 to 4000 characters")
            checkpoint["next_action"]["prompt"] = prompt
            checkpoint["supervisor"]["last_note_receipt"] = {
                "token": token, "recorded_utc": utc_now(),
                "prompt_sha256": hashlib.sha256(prompt.encode()).hexdigest(),
            }
        else:
            path = scoped_path(request.get("path"), writable=True)
            original = path.read_text()
            expected = request.get("expected_sha256")
            if hashlib.sha256(original.encode()).hexdigest() != expected:
                raise ValueError("source hash changed; read the current file before editing")
            old, new = request.get("old_text"), request.get("new_text")
            if not isinstance(old, str) or not old or not isinstance(new, str) or original.count(old) != 1:
                raise ValueError("replacement requires one exact nonempty source match")
            revised = original.replace(old, new, 1)
            if not revised.strip():
                raise PermissionError("emptying or deleting source files is forbidden")
            ast.parse(revised, filename=str(path))
            payload = revised.encode()
            if allocated_bytes(checkpoint) + len(payload) + CHECKPOINT.stat().st_size + OUTPUT_LIMIT >= checkpoint["constraints"]["artifact_limit_bytes"]:
                raise RuntimeError("insufficient artifact budget for source replacement")
            fd, staged = tempfile.mkstemp(prefix="collector-source-", suffix=".py", dir=root)
            with os.fdopen(fd, "wb") as stream:
                os.fchmod(stream.fileno(), path.stat().st_mode & 0o777)
                stream.write(payload)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(staged, path)
            update_coverage(checkpoint)
            checkpoint.setdefault("source_edits", []).append({
                "path": request["path"], "before_sha256": expected,
                "after_sha256": checkpoint["source_sha256"][request["path"]], "recorded_utc": utc_now(),
            })
        save_checkpoint(checkpoint)
        print(json.dumps({"operation": operation, "recorded": True, "source_sha256": checkpoint["source_sha256"]}))
        return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("run", "supervise", "tool"))
    parser.add_argument("tests", nargs="*")
    parser.add_argument("--prove-resumption", action="store_true")
    args = parser.parse_args()
    os.setpriority(os.PRIO_PROCESS, 0, 15)
    if args.mode == "tool":
        return tool_rpc()
    return run_tests(args.tests) if args.mode == "run" else supervise(args.prove_resumption)


if __name__ == "__main__":
    raise SystemExit(main())
