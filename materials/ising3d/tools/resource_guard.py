#!/usr/bin/env python3
"""Serialize and pace one owned macOS research process group; never evict data."""
from __future__ import annotations

import argparse
import ctypes
import errno
import fcntl
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import sys
import time

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent.parent
RUNTIME = WORKSPACE / 'scratch' / 'ising3d-supervisor'
CPU_FRACTION = 0.35
RSS_CAP = 2 * 1024**3
FREE_FLOOR = 50 * 1024**3
FILE_CAP = 16 * 1024**2
WRITE_CAP = 64 * 1024**2
CPU_LIMIT = 1800
DENIED_GRACE_SECONDS = 30
CLEANUP_SECONDS = 5
THREAD_VARS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
               'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS', 'FLINT_NUM_THREADS',
               'BLIS_NUM_THREADS', 'GOTO_NUM_THREADS')


class Usage(ctypes.Structure):
    _fields_ = [('uuid', ctypes.c_ubyte * 16)] + [
        (name, ctypes.c_uint64) for name in (
            'user', 'system', 'wakeups', 'interrupts', 'pageins', 'wired',
            'resident', 'footprint', 'start', 'exit', 'child_user', 'child_system',
            'child_wakeups', 'child_interrupts', 'child_pageins', 'child_elapsed',
            'bytes_read', 'bytes_written',
        )
    ]


class Timebase(ctypes.Structure):
    _fields_ = [('numer', ctypes.c_uint32), ('denom', ctypes.c_uint32)]


lib = ctypes.CDLL('/usr/lib/libproc.dylib', use_errno=True)
lib.proc_pid_rusage.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
lib.proc_pid_rusage.restype = ctypes.c_int
lib.proc_listpgrppids.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_int]
lib.proc_listpgrppids.restype = ctypes.c_int
system = ctypes.CDLL('/usr/lib/libSystem.B.dylib')
system.mach_timebase_info.argtypes = [ctypes.POINTER(Timebase)]
system.mach_timebase_info.restype = ctypes.c_int
timebase = Timebase()
if system.mach_timebase_info(ctypes.byref(timebase)) or not timebase.denom:
    raise RuntimeError('mach_timebase_info failed')
TICK_SECONDS = timebase.numer / timebase.denom / 1e9


def sample_process(pid: int) -> dict:
    """Return current residency/footprint, distinct from process high-water RSS."""
    usage = Usage()
    if lib.proc_pid_rusage(pid, 2, ctypes.byref(usage)):
        raise OSError(ctypes.get_errno(), 'proc_pid_rusage', pid)
    return {
        'pid': pid, 'start_ticks': int(usage.start),
        'exit_ticks': int(usage.exit),
        'resident_bytes': int(usage.resident),
        'footprint_bytes': int(usage.footprint),
        'cpu_seconds': (usage.user + usage.system) * TICK_SECONDS,
        'children_cpu_seconds': (usage.child_user + usage.child_system) * TICK_SECONDS,
        'written_bytes': int(usage.bytes_written),
    }


def free_bytes(cwd: Path) -> int:
    return min(shutil.disk_usage(cwd).free, shutil.disk_usage(RUNTIME).free,
               shutil.disk_usage('/tmp').free)


def group_samples(pid: int, buffer) -> tuple[list[dict], list[int]]:
    """Split the owned group into measured members and unmeasurable live ones."""
    count = lib.proc_listpgrppids(pid, buffer, ctypes.sizeof(buffer))
    if count < 0:
        raise OSError(ctypes.get_errno(), 'proc_listpgrppids')
    if count >= len(buffer):
        raise RuntimeError('owned process group exceeds bounded sampling buffer')
    samples, denied = [], []
    for index in range(count):
        member = buffer[index]
        if member <= 0:
            continue
        try:
            samples.append(sample_process(member))
        except OSError as error:
            if error.errno == errno.EPERM:
                # A set-user-ID descendant (macOS /bin/ps is setuid root) can be
                # neither measured nor signalled by this unprivileged guard. Its
                # CPU still lands in the reaper's children_cpu_seconds, so the
                # blind spot is residency, footprint and written bytes only.
                # supervise() bounds how long one may live instead of trusting it.
                denied.append(member)
            elif error.errno != errno.ESRCH:
                raise
    return samples, denied


def exited_or_gone(pid: int) -> bool:
    """True when the pid provably stopped consuming resources under our watch."""
    try:
        return bool(sample_process(pid)['exit_ticks'])
    except ProcessLookupError:
        return True
    except PermissionError:
        return False


def signal_owned_group(pid: int, signum: int) -> None:
    try:
        os.killpg(pid, signum)
    except ProcessLookupError:
        pass
    except PermissionError:
        # macOS returns EPERM when no member accepted the signal: the unreaped
        # leader is already a zombie, it is reaped, or only an unsignallable
        # privileged member is left. A leader still running is fatal.
        if not exited_or_gone(pid):
            raise


def terminate_owned(pid: int):
    signal_owned_group(pid, signal.SIGTERM)
    signal_owned_group(pid, signal.SIGCONT)
    deadline = time.monotonic() + 2
    while time.monotonic() < deadline:
        ended, status, usage = os.wait4(pid, os.WNOHANG)
        if ended:
            return status, usage
        time.sleep(0.02)
    signal_owned_group(pid, signal.SIGKILL)
    _, status, usage = os.wait4(pid, 0)
    return status, usage


def stop_remaining_owned_group(pid: int, buffer) -> tuple[bool, list[int]]:
    """Kill what survived the leader; report whether the group actually went quiet.

    Measured survivors violate the finite-command contract and must not get an
    unpaced window. Unmeasurable members are signalled but never waited on: a
    genuinely privileged process cannot be killed by this uid, and spinning on
    one would hang the guard inside its own teardown and never free the slot.
    """
    samples, denied = group_samples(pid, buffer)
    if not denied and not any(not sample['exit_ticks'] for sample in samples):
        return True, []
    deadline = time.monotonic() + CLEANUP_SECONDS
    while True:
        signal_owned_group(pid, signal.SIGKILL)
        samples, denied = group_samples(pid, buffer)
        if not any(not sample['exit_ticks'] for sample in samples):
            return True, denied
        if time.monotonic() >= deadline:
            return False, denied
        time.sleep(0.01)


def supervise(command: list[str], cwd: Path) -> dict:
    RUNTIME.mkdir(parents=True, exist_ok=True)
    with (RUNTIME / 'worker.lock').open('a') as slot:
        try:
            fcntl.flock(slot.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('a healthy guarded worker holds the single worker slot')
        if free_bytes(cwd) < FREE_FLOOR:
            raise RuntimeError('free disk below 50 GiB reserve')
        started = time.monotonic()
        own_cpu = time.process_time()
        pid = os.fork()
        if pid == 0:
            try:
                os.setsid()
                os.chdir(cwd)
                os.nice(max(0, 19 - os.nice(0)))
                for name in THREAD_VARS:
                    os.environ[name] = '1'
                os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
                os.environ['PYTHONUNBUFFERED'] = '1'
                resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
                resource.setrlimit(resource.RLIMIT_FSIZE, (FILE_CAP, FILE_CAP))
                resource.setrlimit(resource.RLIMIT_CPU, (CPU_LIMIT, CPU_LIMIT + 1))
                os.kill(os.getpid(), signal.SIGSTOP)
                os.execvp(command[0], command)
            finally:
                os._exit(127)
        cancelled = False

        def cancel(signum, frame):
            nonlocal cancelled
            cancelled = True

        previous_handlers = {sig: signal.signal(sig, cancel)
                             for sig in (signal.SIGTERM, signal.SIGINT)}
        reason = 'completed'
        final_usage = None
        peak = 0
        cpu = 0.0
        written_by_process = {}
        current_resident = current_footprint = 0
        previous_report = previous_disk_check = started
        previous_cpu = 0.0
        max_window_percent = 0.0
        buffer = (ctypes.c_int * 1024)()
        cleanup_quiet = True
        denied_at_exit: list[int] = []
        denied_since: dict[int, float] = {}
        denied_seen: set[int] = set()
        try:
            _, status = os.waitpid(pid, os.WUNTRACED)
            if not os.WIFSTOPPED(status):
                raise RuntimeError('owned worker did not enter its initial stop')
            print(f'RESOURCE_GUARD READY pid={pid} target_cpu=35%_of_one_core '
                  f'free_gib={free_bytes(cwd) / 1024**3:.2f} reserve_gib=50 '
                  'rss_cap_gib=2 file_cap_mib=16 write_cap_mib=64', flush=True)
            while True:
                if cancelled:
                    reason = 'cancelled'
                    break
                cycle_started = time.monotonic()
                signal_owned_group(pid, signal.SIGCONT)
                time.sleep(0.04)
                signal_owned_group(pid, signal.SIGSTOP)
                # Read the unreaped leader's final counters even when it exits
                # inside the first quantum. wait4 destroys that native record.
                samples, denied = group_samples(pid, buffer)
                if pid not in {item['pid'] for item in samples}:
                    samples.append(sample_process(pid))
                denied_seen.update(denied)
                denied_since = {member: denied_since.get(member, cycle_started)
                                for member in denied}
                current_resident = sum(item['resident_bytes'] for item in samples)
                current_footprint = sum(item['footprint_bytes'] for item in samples)
                peak = max(peak, current_resident, current_footprint)
                # Reaped descendants transfer their accounting into the parent's
                # child totals. Keep a high-water total across that transition.
                cpu = max(cpu, sum(item['cpu_seconds'] + item['children_cpu_seconds']
                                   for item in samples))
                for item in samples:
                    key = (item['pid'], item['start_ticks'])
                    written_by_process[key] = max(written_by_process.get(key, 0),
                                                  item['written_bytes'])
                _, status, usage = os.wait4(pid, os.WUNTRACED)
                if os.WIFEXITED(status) or os.WIFSIGNALED(status):
                    final_usage = usage
                    remaining, remaining_denied = group_samples(pid, buffer)
                    # Witness a privileged member that only shows up at exit; it
                    # is not a contract breach, so it never fails a clean run.
                    denied_seen.update(remaining_denied)
                    if any(not item['exit_ticks'] for item in remaining):
                        reason = 'worker exited with background descendants'
                    break
                now = time.monotonic()
                elapsed = now - started
                if peak >= RSS_CAP:
                    reason = '2 GiB memory ceiling reached'
                    break
                if sum(written_by_process.values()) >= WRITE_CAP:
                    reason = '64 MiB write ceiling reached'
                    break
                if cpu >= CPU_LIMIT:
                    reason = '1800 CPU-second ceiling reached'
                    break
                if denied_since and now - min(denied_since.values()) >= DENIED_GRACE_SECONDS:
                    # Distinct wording: this is an accounting failure, not
                    # resource exhaustion, and the supervisor must not retire a
                    # legitimate trial over it.
                    reason = ('unmeasurable descendant outlived the '
                              f'{DENIED_GRACE_SECONDS:.0f} s accounting window')
                    break
                if now - previous_disk_check >= 1:
                    previous_disk_check = now
                    if free_bytes(cwd) < FREE_FLOOR:
                        reason = 'free disk below 50 GiB reserve'
                        break
                if now - previous_report >= 60:
                    percent = 100 * (cpu - previous_cpu) / (now - previous_report)
                    max_window_percent = max(max_window_percent, percent)
                    print('RESOURCE_GUARD SAMPLE ' + json.dumps({
                        'cpu_window_percent': round(percent, 2),
                        'cpu_total_seconds': round(cpu, 3),
                        'wall_seconds': round(elapsed, 2),
                        'current_resident_mib': round(current_resident / 1024**2, 2),
                        'current_footprint_mib': round(current_footprint / 1024**2, 2),
                        'peak_memory_mib': round(peak / 1024**2, 2),
                        'written_mib': round(sum(written_by_process.values()) / 1024**2, 2),
                        'free_gib': round(free_bytes(cwd) / 1024**3, 2),
                    }, sort_keys=True), flush=True)
                    previous_report, previous_cpu = now, cpu
                running_window = now - cycle_started
                debt = (cpu + time.process_time() - own_cpu) / CPU_FRACTION - elapsed
                time.sleep(max(running_window * (1 / CPU_FRACTION - 1), debt, 0.0))
        except BaseException as error:
            reason = f'guard failure: {type(error).__name__}: {error}'
        finally:
            if final_usage is None:
                status, final_usage = terminate_owned(pid)
            cleanup_quiet, denied_at_exit = stop_remaining_owned_group(pid, buffer)
            denied_seen.update(denied_at_exit)
            for sig, handler in previous_handlers.items():
                signal.signal(sig, handler)
        peak = max(peak, int(final_usage.ru_maxrss))  # macOS reports bytes.
        if reason == 'completed' and peak >= RSS_CAP:
            reason = '2 GiB memory ceiling reached at process exit'
        if reason == 'completed' and sum(written_by_process.values()) >= WRITE_CAP:
            reason = '64 MiB write ceiling reached at process exit'
        if reason == 'completed' and free_bytes(cwd) < FREE_FLOOR:
            reason = 'free disk below 50 GiB reserve at process exit'
        if reason == 'completed' and not cleanup_quiet:
            reason = 'owned group survived its bounded cleanup'
        total_cpu = max(cpu, final_usage.ru_utime + final_usage.ru_stime)
        total_cpu += time.process_time() - own_cpu
        command_wall = time.monotonic() - started
        # Hold the exclusive slot for the final off-time, including short jobs.
        # This is admission pacing, not a claim about instantaneous CPU usage.
        time.sleep(max(0.0, total_cpu / CPU_FRACTION - command_wall))
        wall = time.monotonic() - started
        child_exit = os.waitstatus_to_exitcode(status)
        guard_exit = (child_exit if child_exit >= 0 else 128 - child_exit) if reason == 'completed' else 1
        return {
            'guard_exit_code': guard_exit,
            'child_exit_code': child_exit,
            'reason': reason,
            'cpu_seconds': round(total_cpu, 6),
            'command_wall_seconds': round(command_wall, 6),
            'wall_seconds': round(wall, 6),
            'average_cpu_percent': round(100 * total_cpu / wall, 4),
            'max_60s_cpu_percent': round(max_window_percent, 4),
            'peak_memory_bytes': peak,
            'written_bytes': sum(written_by_process.values()),
            'unmeasurable_processes': len(denied_seen),
            'unmeasurable_at_exit': len(denied_at_exit),
            'free_disk_bytes': free_bytes(cwd),
        }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cwd', type=Path, default=ROOT.parent)
    parser.add_argument('--result', type=Path)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command:
        parser.error('a finite command is required after --')
    try:
        result = supervise(command, args.cwd.resolve())
    except Exception as error:
        result = {'guard_exit_code': 75, 'child_exit_code': None,
                  'reason': str(error), 'cpu_seconds': 0, 'wall_seconds': 0,
                  'average_cpu_percent': 0, 'peak_memory_bytes': 0,
                  'written_bytes': 0, 'free_disk_bytes': shutil.disk_usage(args.cwd).free}
    if args.result:
        args.result.parent.mkdir(parents=True, exist_ok=True)
        with args.result.open('x') as stream:
            json.dump(result, stream, sort_keys=True, separators=(',', ':'))
            stream.write('\n')
    print('RESOURCE_GUARD DONE ' + json.dumps(result, sort_keys=True), flush=True)
    return result['guard_exit_code']


if __name__ == '__main__':
    raise SystemExit(main())
